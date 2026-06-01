"""Distances analysis plugin.

Computes inter-atomic distances from MD trajectories for one or more
atom pairs, aggregates across replicates with SEM, and performs per-pair
cross-condition comparisons with t-tests, ANOVA, and rankings.

This plugin contains the **full distance calculation implementation**,
including trajectory loading, PBC-aware distance computation, KDE-based
distribution analysis, autocorrelation-corrected uncertainty, and
threshold-based contact analysis.

Unlike single-scalar analyses (RMSF, catalytic_triad), distances has **no
single primary metric** — each distance pair is compared independently since
averaging unrelated distances (e.g. H-bond distance + lid-opening distance)
is not semantically meaningful.  Therefore ``compare()`` is overridden
entirely and ``extract_metrics()`` is not used.

"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses._framework.cache_identity import settings_fingerprint
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    PlotContext,
    PluginContractError,
)
from polyzymd.analyses.distances._mda import (
    DistanceArtifactCollector,
    aggregate_distance_artifacts,
    build_distance_jobs,
)
from polyzymd.analyses.distances._plot_settings import DistancesPlotSettings
from polyzymd.analyses.distances._plotters import (
    _plot_distance_kde,
    _plot_distance_state_bars,
    _plot_distance_threshold_bars,
    build_distance_plot_data,
)
from polyzymd.analyses.mda import (
    ConditionArtifact,
    MDACollectorContext,
    MDAReplicateJobContext,
    ReplicateArtifact,
)

if TYPE_CHECKING:
    from polyzymd.analyses.distances._comparison_results import (
        DistanceComparisonResult,
    )

logger = logging.getLogger("polyzymd.analyses.distances")

NOT_TESTABLE_SINGLETON_NOTE = (
    "Inferential statistics require at least two replicates per condition."
)

# Default threshold from the existing settings module
DEFAULT_DISTANCE_THRESHOLD = 3.5


@dataclass(frozen=True)
class _DistancesTrajectoryWindow:
    """Distances trajectory window that carries loader-derived file metadata."""

    start: int
    stop: int
    step: int
    equilibration_start: int
    n_frames_total: int
    n_frames_selected: int
    timestep_ps: float
    equilibration_ps: float
    warning_message: str | None = None
    trajectory_files: tuple[Path, ...] = ()

    @classmethod
    def from_window(
        cls,
        window: Any,
        trajectory_files: Sequence[Path],
    ) -> _DistancesTrajectoryWindow:
        """Build a distances window wrapper from the shared trajectory window.

        Parameters
        ----------
        window : Any
            Shared trajectory window returned by the centralized resolver.
        trajectory_files : Sequence[Path]
            Trajectory files resolved by the existing loader instance.

        Returns
        -------
        _DistancesTrajectoryWindow
            Distances window wrapper that preserves run arguments and file metadata.
        """

        return cls(
            start=int(window.start),
            stop=int(window.stop),
            step=int(window.step),
            equilibration_start=int(window.equilibration_start),
            n_frames_total=int(window.n_frames_total),
            n_frames_selected=int(window.n_frames_selected),
            timestep_ps=float(window.timestep_ps),
            equilibration_ps=float(window.equilibration_ps),
            warning_message=getattr(window, "warning_message", None),
            trajectory_files=tuple(trajectory_files),
        )

    def run_kwargs(self) -> dict[str, int]:
        """Return keyword arguments for the runner ``run()`` call."""

        return {"start": self.start, "stop": self.stop, "step": self.step}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class DistancePairSettings(BaseModel):
    """Configuration for a single distance pair.

    Attributes
    ----------
    label : str
        Human-readable label (e.g. ``"Ser77-Substrate"``).
    selection_a : str
        MDAnalysis selection string for the first atom/group.
    selection_b : str
        MDAnalysis selection string for the second atom/group.
    threshold : float | None
        Per-pair distance threshold (Angstroms).  If ``None``, uses the
        global threshold from :class:`DistancesSettings`.
    below_label : str | None
        Display label for the "below threshold" state (e.g. ``"Bound"``).
    above_label : str | None
        Display label for the "above threshold" state (e.g. ``"Unbound"``).
    """

    label: str = Field(..., description="Human-readable label for this pair")
    selection_a: str = Field(..., description="First atom/point selection")
    selection_b: str = Field(..., description="Second atom/point selection")
    threshold: float | None = Field(
        default=None,
        description="Per-pair threshold (Angstroms). Falls back to global threshold.",
    )
    below_label: str | None = Field(
        default=None,
        description='Display label for "below threshold" state.',
    )
    above_label: str | None = Field(
        default=None,
        description='Display label for "above threshold" state.',
    )


class DistancesSettings(BaseModel):
    """Settings for distances analysis.

    Attributes
    ----------
    threshold : float | None
        Global distance threshold for contact analysis (Angstroms).
    pairs : list[DistancePairSettings]
        Atom pairs to measure distances between.
    use_pbc : bool
        Use PBC-aware minimum image distances (default ``True``).
    align_trajectory : bool
        Align trajectory before distance calculation (default ``True``).
    alignment_selection : str
        MDAnalysis selection for trajectory alignment.
    alignment_mode : str
        Reference mode for alignment: ``"centroid"``, ``"average"``, or
        ``"frame"``.
    alignment_frame : int | None
        Reference frame (1-indexed) when ``alignment_mode="frame"``.
    """

    threshold: float | None = Field(
        default=DEFAULT_DISTANCE_THRESHOLD,
        description="Global distance threshold for contact analysis (Angstroms)",
    )
    pairs: list[DistancePairSettings] = Field(
        default_factory=list,
        description="Distance pairs to monitor",
    )
    use_pbc: bool = Field(
        default=True,
        description="Use PBC-aware minimum image distances",
    )
    align_trajectory: bool = Field(
        default=True,
        description="Align trajectory before distance calculation",
    )
    alignment_selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for trajectory alignment",
    )
    alignment_mode: str = Field(
        default="centroid",
        description="Reference mode: centroid, average, or frame",
    )
    alignment_frame: int | None = Field(
        default=None,
        description="Reference frame (1-indexed) when alignment_mode='frame'",
    )

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs(cls, v: list[DistancePairSettings]) -> list[DistancePairSettings]:
        """Ensure at least one pair is defined."""
        if len(v) == 0:
            raise ValueError("At least one distance pair must be defined")
        return v

    @field_validator("alignment_mode", mode="after")
    @classmethod
    def validate_alignment_mode(cls, v: str) -> str:
        """Validate alignment mode."""
        valid = {"centroid", "average", "frame"}
        if v not in valid:
            raise ValueError(f"alignment_mode must be one of {valid}, got '{v}'")
        return v

    @model_validator(mode="after")
    def validate_alignment_frame_required(self) -> "DistancesSettings":
        """Ensure alignment_frame is set when alignment_mode is 'frame'."""
        if (
            self.align_trajectory
            and self.alignment_mode == "frame"
            and self.alignment_frame is None
        ):
            raise ValueError("alignment_frame is required when alignment_mode is 'frame'")
        return self

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of ``(selection_a, selection_b)`` tuples."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]

    def get_pair_thresholds(self) -> list[float | None]:
        """Get per-pair thresholds, falling back to the global threshold."""
        return [p.threshold if p.threshold is not None else self.threshold for p in self.pairs]

    def get_alignment_config(self) -> Any:
        """Build an ``AlignmentConfig`` from these settings.

        Returns
        -------
        AlignmentConfig
            Configuration for trajectory alignment.
        """
        from polyzymd.analyses.shared.alignment import AlignmentConfig

        return AlignmentConfig(
            enabled=self.align_trajectory,
            reference_mode=self.alignment_mode,
            reference_frame=self.alignment_frame,
            selection=self.alignment_selection,
        )


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class DistancesAnalysis(Analysis):
    """Distances analysis: inter-atomic distances from MD trajectories.

    This plugin performs full distance computation inline (no delegation
    to non-canonical calculator classes), aggregates across replicates, and
    performs per-pair cross-condition comparison with dual metrics (mean
    distance and fraction below threshold).

    The ``compare()`` method is **fully overridden** because:

    - Each distance pair is compared independently (no single scalar).
    - Two metrics per pair: mean distance (lower is better) and fraction
      below threshold (higher is better, optional).
    - Rankings are per-pair, not overall.

    Plots
    -----
    - ``distance_kde_*.png``: KDE distribution overlay per pair.
    - ``distance_threshold_bars.png``: Grouped bar chart of fraction below
      threshold.
    - ``distance_state_*.png``: Per-pair above/below state bar charts.
    """

    name: ClassVar[str] = "distances"
    Settings: ClassVar[type] = DistancesSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = DistancesPlotSettings
    AggregatedResultClass: ClassVar[type | None] = None
    ReplicateResultClass: ClassVar[type | None] = None
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 1

    # === Required methods ===

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel) -> str:
        """Build a short cache tag for settings.

        Parameters
        ----------
        settings : BaseModel
            Analysis settings model.

        Returns
        -------
        str
            First 8 hex characters from shared settings fingerprinting.
        """
        return settings_fingerprint(settings)

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[Any] | None:
        """Build the MDAnalysis-native pair-distance job for one replicate."""

        return build_distance_jobs(ctx, ctx.settings)

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the distances artifact collector."""

        del ctx
        return DistanceArtifactCollector()

    def _deserialize_result(self, path: Path) -> Any:
        """Load a canonical distances condition artifact."""

        if path.exists():
            try:
                loaded = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise PluginContractError(
                    f"Distances aggregate at {path} is not a valid canonical artifact. "
                    "Recompute the condition or clear stale caches before comparing."
                ) from exc
            if isinstance(loaded, dict) and loaded.get("artifact_type") == "condition":
                from polyzymd.analyses.mda import ArtifactStore

                return ArtifactStore(path.parent).read_condition_result(path.name)
        raise PluginContractError(
            f"Distances aggregate at {path} is not a canonical MDAnalysis condition artifact. "
            "Recompute the condition or clear stale non-canonical distances caches."
        )

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate distance replicate artifacts across one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : sequence of ReplicateArtifact
            Per-replicate MDAnalysis distance artifacts.

        Returns
        -------
        ConditionArtifact
            Aggregated condition artifact with canonical pair summaries.
        """

        if not results:
            raise ValueError(
                f"Distances aggregation for condition '{ctx.condition.label}' requires at least "
                "one replicate artifact. No replicate inputs were provided."
            )
        if not all(isinstance(result, ReplicateArtifact) for result in results):
            raise TypeError(
                "Distances aggregation expects MDAnalysis ReplicateArtifact inputs. Non-canonical "
                "distances replicate caches are incompatible with the MDAnalysis artifact "
                "lifecycle; recompute the condition or clear stale caches before aggregating."
            )
        return aggregate_distance_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            output_dir=ctx.output_dir,
            artifacts=results,
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
        )

    @staticmethod
    def _validate_replicate_pair_schema(
        ctx: AggregateContext,
        results: Sequence[Any],
        settings: DistancesSettings,
    ) -> None:
        """Require all replicate pair results to match the configured schema.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate distance results.
        settings : DistancesSettings
            Current distances settings that define the canonical pair schema.

        Raises
        ------
        ValueError
            Raised when any replicate result has a mismatched pair count, pair
            ordering, selections, labels, or effective threshold.
        """

        expected_thresholds = settings.get_pair_thresholds()
        schema_issues: list[str] = []

        for result in results:
            replicate = getattr(result, "replicate", None)
            pair_results = list(getattr(result, "pair_results", []))

            if len(pair_results) != len(settings.pairs):
                schema_issues.append(
                    f"replicate {replicate}: pair count {len(pair_results)} != {len(settings.pairs)}"
                )
                continue

            for pair_idx, (pair_result, pair_setting, threshold) in enumerate(
                zip(pair_results, settings.pairs, expected_thresholds, strict=True),
                start=1,
            ):
                mismatch_parts: list[str] = []
                expected_pair_label = pair_setting.label
                if pair_result.pair_label != expected_pair_label:
                    mismatch_parts.append(
                        f"label {pair_result.pair_label!r} != {expected_pair_label!r}"
                    )
                if pair_result.selection1 != pair_setting.selection_a:
                    mismatch_parts.append(
                        f"selection1 {pair_result.selection1!r} != {pair_setting.selection_a!r}"
                    )
                if pair_result.selection2 != pair_setting.selection_b:
                    mismatch_parts.append(
                        f"selection2 {pair_result.selection2!r} != {pair_setting.selection_b!r}"
                    )
                if pair_result.threshold != threshold:
                    mismatch_parts.append(f"threshold {pair_result.threshold!r} != {threshold!r}")

                if mismatch_parts:
                    schema_issues.append(
                        f"replicate {replicate} pair {pair_idx}: {', '.join(mismatch_parts)}"
                    )

        if schema_issues:
            issue_text = "; ".join(schema_issues)
            raise ValueError(
                f"Distances aggregation for condition '{ctx.condition.label}' requires every "
                f"replicate result to match the configured pair schema. Problems detected: "
                f"{issue_text}."
            )

    # === compare() — fully overridden ===

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare distance metrics across conditions.

        Each distance pair is compared independently — rankings and statistics
        are computed per-pair since averaging unrelated distances is not
        semantically meaningful.

        Primary ranking: mean distance (lower = closer = better).
        Secondary ranking: fraction below threshold (higher = more contact = better).

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided context.

        Returns
        -------
        DistanceComparisonResult | None
            Comparison result, or ``None`` if no conditions have data.
        """
        from polyzymd import __version__
        from polyzymd.analyses.distances._comparison_results import (
            DistanceComparisonResult,
            DistanceConditionSummary,
            DistancePairANOVA,
            DistancePairSummary,
            DistancePairwiseComparison,
        )
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            one_way_anova,
            percent_change,
        )

        settings = ctx.settings
        pair_labels = settings.get_pair_labels()

        logger.info(f"Starting distance comparison: {ctx.name}")
        logger.info(f"Conditions: {len(ctx.conditions)}, Pairs: {pair_labels}")

        # Load aggregated results for each condition
        summaries: list[DistanceConditionSummary] = []

        for cond in ctx.conditions:
            # Prefer in-memory results from orchestrator, fall back to disk
            agg_result = ctx.aggregated_results.get(cond.label)
            agg_source: Path | str = f"comparison context for {cond.label}"
            if agg_result is None:
                agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
                agg_source = self.aggregate_result_path(agg_dir)
                agg_result = self._load_aggregated_result(agg_dir)

            if agg_result is None:
                logger.warning(f"No aggregated result found for '{cond.label}' — skipping.")
                continue

            if not isinstance(agg_result, ConditionArtifact):
                raise TypeError(
                    f"Distances comparison for condition '{cond.label}' requires canonical "
                    "MDAnalysis condition artifacts. Non-artifact aggregate inputs are "
                    "incompatible; recompute the condition or clear stale caches before comparing."
                )
            agg_result = self.validate_aggregated_result(
                agg_result,
                condition=cond,
                settings=ctx.settings,
                equilibration=ctx.equilibration,
                source=agg_source,
                expected_replicates=cond.replicates,
                allow_replicate_subset=True,
            )

            pair_results = list(agg_result.payload.get("pair_results", []))
            pair_summaries = []
            for pair_index, pr in enumerate(pair_results):
                expected_label = pair_labels[pair_index] if pair_index < len(pair_labels) else None
                pair_label = str(pr.get("pair_label", f"Pair {pair_index}"))
                user_label = pair_label if pair_label == expected_label else expected_label
                if user_label is None:
                    user_label = pair_label
                pair_summaries.append(
                    DistancePairSummary(
                        label=user_label,
                        selection_a=str(pr.get("selection1", "")),
                        selection_b=str(pr.get("selection2", "")),
                        threshold=pr.get("threshold"),
                        mean_distance=float(pr.get("overall_mean", pr.get("mean_distance", 0.0))),
                        sem_distance=float(
                            pr.get("overall_sem", pr.get("sem_distance", 0.0)) or 0.0
                        ),
                        fraction_below_threshold=pr.get("overall_fraction_below"),
                        sem_fraction_below=pr.get("sem_fraction_below"),
                        per_replicate_means=[
                            float(value) for value in pr.get("per_replicate_means", [])
                        ],
                        per_replicate_fractions=(
                            [float(value) for value in pr.get("per_replicate_fractions_below", [])]
                            if pr.get("per_replicate_fractions_below") is not None
                            else None
                        ),
                    )
                )

            summaries.append(
                DistanceConditionSummary(
                    label=cond.label,
                    config_path=str(cond.config_path),
                    n_replicates=int(
                        agg_result.payload.get("n_replicates", len(agg_result.replicates))
                    ),
                    pair_summaries=pair_summaries,
                )
            )

        if not summaries:
            logger.warning("distances: no conditions have results — skipping comparison.")
            return None

        # Resolve effective control
        effective_control = ctx.effective_control

        # Per-pair rankings
        ranking_by_pair: dict[str, list[str]] = {}
        fraction_ranking_by_pair: dict[str, list[str]] = {}

        for pair_label in pair_labels:
            pair_data = []
            for summary in summaries:
                pair_summary = summary.get_pair(pair_label)
                pair_data.append((summary.label, pair_summary))

            # Rank by mean distance (ascending — lower is better)
            sorted_by_distance = sorted(pair_data, key=lambda x: x[1].mean_distance)
            ranking_by_pair[pair_label] = [label for label, _ in sorted_by_distance]

            # Rank by fraction below threshold (descending — higher is better)
            with_fraction = [
                (label, ps) for label, ps in pair_data if ps.fraction_below_threshold is not None
            ]
            if with_fraction:
                sorted_by_fraction = sorted(
                    with_fraction,
                    key=lambda x: x[1].fraction_below_threshold or 0,
                    reverse=True,
                )
                fraction_ranking_by_pair[pair_label] = [label for label, _ in sorted_by_fraction]

        # Pairwise comparisons (per-pair)
        comparisons: list[DistancePairwiseComparison] = []

        if len(summaries) >= 2:
            for pair_label in pair_labels:
                if effective_control:
                    control = next(s for s in summaries if s.label == effective_control)
                    control_pair = control.get_pair(pair_label)
                    for summary in summaries:
                        if summary.label == effective_control:
                            continue
                        treatment_pair = summary.get_pair(pair_label)
                        comp = self._compare_pair(
                            pair_label,
                            control.label,
                            summary.label,
                            control_pair,
                            treatment_pair,
                        )
                        comparisons.append(comp)
                else:
                    for i, summary_a in enumerate(summaries):
                        pair_a = summary_a.get_pair(pair_label)
                        for summary_b in summaries[i + 1 :]:
                            pair_b = summary_b.get_pair(pair_label)
                            comp = self._compare_pair(
                                pair_label,
                                summary_a.label,
                                summary_b.label,
                                pair_a,
                                pair_b,
                            )
                            comparisons.append(comp)

        # ANOVA (if 3+ conditions) — per-pair
        anova_by_pair: list[DistancePairANOVA] | None = None
        if len(summaries) >= 3:
            anova_by_pair = []
            for pair_label in pair_labels:
                distance_groups = []
                fraction_groups = []
                for summary in summaries:
                    pair_data = summary.get_pair(pair_label)
                    distance_groups.append(pair_data.per_replicate_means)
                    if pair_data.per_replicate_fractions:
                        fraction_groups.append(pair_data.per_replicate_fractions)

                anova_dist = one_way_anova(*distance_groups)
                distance_testable = all(len(group) >= 2 for group in distance_groups)

                fraction_f = None
                fraction_p = None
                fraction_sig = None
                fraction_testable = None
                fraction_note = None
                if len(fraction_groups) == len(summaries):
                    anova_frac = one_way_anova(*fraction_groups)
                    fraction_testable = all(len(group) >= 2 for group in fraction_groups)
                    fraction_f = anova_frac.f_statistic
                    fraction_p = anova_frac.p_value
                    fraction_sig = anova_frac.significant if fraction_testable else False
                    fraction_note = None if fraction_testable else NOT_TESTABLE_SINGLETON_NOTE

                anova_by_pair.append(
                    DistancePairANOVA(
                        pair_label=pair_label,
                        distance_f_statistic=anova_dist.f_statistic,
                        distance_p_value=anova_dist.p_value,
                        distance_significant=anova_dist.significant if distance_testable else False,
                        distance_testable=distance_testable,
                        distance_note=(None if distance_testable else NOT_TESTABLE_SINGLETON_NOTE),
                        fraction_f_statistic=fraction_f,
                        fraction_p_value=fraction_p,
                        fraction_significant=fraction_sig,
                        fraction_testable=fraction_testable,
                        fraction_note=fraction_note,
                    )
                )

        # Build result
        result = DistanceComparisonResult(
            metric="mean_distance",
            name=ctx.name,
            n_pairs=len(settings.pairs),
            pair_labels=pair_labels,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova_by_pair=anova_by_pair,
            ranking_by_pair=ranking_by_pair,
            fraction_ranking_by_pair=(
                fraction_ranking_by_pair if fraction_ranking_by_pair else None
            ),
            equilibration_time=ctx.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

        # Log rankings
        for pair_label in pair_labels:
            logger.info(f"Ranking for '{pair_label}': {ranking_by_pair[pair_label]}")

        return result

    # === plot() ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate distance comparison plots.

        Produces:

        - ``distance_kde_*.png``: KDE distribution overlay per pair.
        - ``distance_threshold_bars.png``: Grouped bar chart of fraction below threshold.
        - ``distance_state_*.png``: Per-pair above/below state bar charts.

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided context.

        Returns
        -------
        list[Path]
            Paths to generated figure files.
        """
        plots: list[Path] = []

        data, labels = self._build_plot_data(ctx, include_replicates=True)
        if not labels:
            return plots

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        plot_settings = ctx.plot_settings
        plot_data = build_distance_plot_data(data, labels)

        try:
            result = _plot_distance_kde(plot_data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except (ValueError, RuntimeError, OSError) as exc:
            logger.warning(f"Distance KDE plot failed: {exc}")

        try:
            result = _plot_distance_threshold_bars(plot_data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except (ValueError, RuntimeError, OSError) as exc:
            logger.warning(f"Distance threshold bars plot failed: {exc}")

        try:
            result = _plot_distance_state_bars(plot_data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except (ValueError, RuntimeError, OSError) as exc:
            logger.warning(f"Distance state bars plot failed: {exc}")

        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format distance comparison results without compatibility dispatch."""
        from polyzymd.analyses.distances._formatters import format_distances_result

        return format_distances_result(result, format=self._normalize_output_format(output_format))

    @staticmethod
    def _normalize_output_format(output_format: str) -> str:
        return "table" if output_format == "text" else output_format

    # === Private helpers ===

    @staticmethod
    def _compare_pair(
        pair_label: str,
        cond_a_label: str,
        cond_b_label: str,
        pair_a: Any,
        pair_b: Any,
    ) -> Any:
        """Compare two conditions statistically for a single distance pair.

        Parameters
        ----------
        pair_label : str
            Label of the distance pair being compared.
        cond_a_label : str
            Label of first condition (typically control).
        cond_b_label : str
            Label of second condition (typically treatment).
        pair_a : DistancePairSummary
            Pair data from condition A.
        pair_b : DistancePairSummary
            Pair data from condition B.

        Returns
        -------
        DistancePairwiseComparison
            Statistical comparison result.
        """
        from polyzymd.analyses.distances._comparison_results import DistancePairwiseComparison
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )
        from polyzymd.analyses.stats import interpret_direction

        # Distance metric comparison
        values_a = pair_a.per_replicate_means
        values_b = pair_b.per_replicate_means
        distance_testable = len(values_a) >= 2 and len(values_b) >= 2

        ttest_dist = independent_ttest(values_a, values_b)
        effect_dist = cohens_d(values_a, values_b)
        pct_dist = percent_change(pair_a.mean_distance, pair_b.mean_distance)

        # Direction: negative change = closer = improving
        direction_dist = interpret_direction(
            pct_dist,
            direction_labels=("closer", "unchanged", "farther"),
            threshold=1.0,
        )

        # Fraction metric comparison (optional)
        fraction_t = None
        fraction_p = None
        fraction_d = None
        fraction_interp = None
        fraction_dir = None
        fraction_sig = None
        fraction_pct = None
        fraction_testable = None
        fraction_note = None

        if pair_a.per_replicate_fractions and pair_b.per_replicate_fractions:
            frac_a = pair_a.per_replicate_fractions
            frac_b = pair_b.per_replicate_fractions
            fraction_testable = len(frac_a) >= 2 and len(frac_b) >= 2

            ttest_frac = independent_ttest(frac_a, frac_b)
            effect_frac = cohens_d(frac_a, frac_b)
            pct_frac = percent_change(
                pair_a.fraction_below_threshold or 0,
                pair_b.fraction_below_threshold or 0,
            )

            fraction_t = ttest_frac.t_statistic
            fraction_p = ttest_frac.p_value
            fraction_d = effect_frac.cohens_d
            fraction_interp = effect_frac.interpretation
            fraction_sig = ttest_frac.significant if fraction_testable else False
            fraction_note = None if fraction_testable else NOT_TESTABLE_SINGLETON_NOTE

            # Direction: positive change = more contact = improving
            fraction_dir = interpret_direction(
                pct_frac,
                direction_labels=("less_contact", "unchanged", "more_contact"),
                threshold=1.0,
            )

            fraction_pct = pct_frac

        return DistancePairwiseComparison(
            pair_label=pair_label,
            condition_a=cond_a_label,
            condition_b=cond_b_label,
            distance_t_statistic=ttest_dist.t_statistic,
            distance_p_value=ttest_dist.p_value,
            distance_cohens_d=effect_dist.cohens_d,
            distance_effect_interpretation=effect_dist.interpretation,
            distance_direction=direction_dist,
            distance_significant=ttest_dist.significant if distance_testable else False,
            distance_percent_change=pct_dist,
            distance_testable=distance_testable,
            distance_note=None if distance_testable else NOT_TESTABLE_SINGLETON_NOTE,
            fraction_t_statistic=fraction_t,
            fraction_p_value=fraction_p,
            fraction_cohens_d=fraction_d,
            fraction_effect_interpretation=fraction_interp,
            fraction_direction=fraction_dir,
            fraction_significant=fraction_sig,
            fraction_percent_change=fraction_pct,
            fraction_testable=fraction_testable,
            fraction_note=fraction_note,
        )
