"""Catalytic triad analysis plugin.

Analyses active-site geometry from MD trajectories by computing per-pair
distances and a simultaneous contact fraction (all pairs below threshold in the
same frame). Aggregation, comparison, and plotting remain in the plugin while
trajectory-native pair distances flow through the runner seam.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.catalytic_triad._plot_settings import TriadPlotSettings
from polyzymd.analyses.catalytic_triad._plotters import (
    plot_triad_kde_panel_from_data,
    plot_triad_threshold_bars_from_data,
)
from polyzymd.analyses.catalytic_triad._results import TriadAggregatedResult, TriadResult
from polyzymd.analyses.catalytic_triad._runner import CatalyticTriadReplicateRunner
from polyzymd.analyses.distances import _make_pair_label
from polyzymd.analyses.distances._results import (
    DistancePairAggregatedResult,
    DistancePairResult,
    DistanceResultMetadata,
)
from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
from polyzymd.analyses.shared.loader import TrajectoryLoader, parse_time_string

logger = logging.getLogger("polyzymd.analyses.catalytic_triad")


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class TriadPairSettings(BaseModel):
    """Configuration for a single distance pair in a catalytic triad.

    Attributes
    ----------
    label : str
        Human-readable pair label (e.g. ``"Asp133-His156"``).
    selection_a : str
        MDAnalysis/PolyzyMD selection string for the first atom/point.
    selection_b : str
        MDAnalysis/PolyzyMD selection string for the second atom/point.
    """

    label: str = Field(..., description="Human-readable label for this pair")
    selection_a: str = Field(..., description="First atom/point selection")
    selection_b: str = Field(..., description="Second atom/point selection")


class CatalyticTriadSettings(BaseModel):
    """Settings for catalytic triad analysis.

    Attributes
    ----------
    name : str
        Name of the triad/active site (e.g. ``"LipA_catalytic_triad"``).
    pairs : list[TriadPairSettings]
        Distance pairs to monitor.
    threshold : float
        Contact threshold in Angstroms (default 3.5).
    description : str | None
        Optional description.
    """

    name: str = Field(
        default="catalytic_triad",
        description="Name of the catalytic triad/active site",
    )
    pairs: list[TriadPairSettings] = Field(..., description="Distance pairs to monitor")
    threshold: float = Field(
        default=3.5,
        description="Contact threshold in Angstroms",
    )
    description: str | None = Field(
        default=None,
        description="Description of the active site",
    )

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs(cls, v: list[TriadPairSettings]) -> list[TriadPairSettings]:
        """Ensure at least one pair is defined."""
        if len(v) == 0:
            raise ValueError("At least one distance pair must be defined")
        return v

    @property
    def n_pairs(self) -> int:
        """Number of distance pairs."""
        return len(self.pairs)

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of ``(selection_a, selection_b)`` tuples."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class CatalyticTriadAnalysis(Analysis):
    """Catalytic triad analysis: active-site geometry from MD trajectories.

    Computes per-pair distances through the shared runner seam and derives a
    simultaneous contact fraction — the percentage of frames where ALL pairs are
    below the threshold at the same time.

    The ``compare()`` method is NOT overridden — it uses the default
    implementation which calls ``extract_metrics()`` to get
    ``simultaneous_contact_fraction`` as a single scalar metric with
    ``higher_is_better=True`` (more contact = better triad integrity).

    Plots
    -----
    - ``triad_kde_panel.png``: Multi-row KDE of per-pair distance distributions.
    - ``triad_threshold_bars.png``: Grouped bar chart of contact fractions.
    """

    name: ClassVar[str] = "catalytic_triad"
    Settings: ClassVar[type] = CatalyticTriadSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = TriadPlotSettings
    AggregatedResultClass: ClassVar[type] = TriadAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = TriadResult
    aliases: ClassVar[tuple[str, ...]] = ("triad",)
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 1

    # === Required methods ===

    def run_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Run triad analysis for a single replicate through the runner seam."""

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        result_file = ctx.output_dir / _make_result_filename(ctx.settings, eq_value, eq_unit)
        cached = self._check_cache(
            TriadResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=ctx.sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

        result = super().run_replicate(ctx, replicate)
        result.save(result_file)
        logger.info("Saved triad result to %s", result_file)
        return result

    @staticmethod
    def _settings_cache_tag(settings: BaseModel) -> str:
        """Return the catalytic-triad settings fingerprint."""

        return settings_fingerprint(settings)

    @classmethod
    def _validate_replicate_result_settings_identity(
        cls,
        ctx: AggregateContext,
        results: Sequence[TriadResult],
    ) -> None:
        """Validate settings fingerprints on per-replicate triad results.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[TriadResult]
            Per-replicate catalytic-triad results.

        Raises
        ------
        ValueError
            Raised when replicate results are missing settings fingerprints or
            were computed with different settings.
        """

        expected_fingerprint = cls._settings_cache_tag(ctx.settings)
        missing_fingerprint_replicates: list[int] = []
        mismatched_fingerprints: list[str] = []

        for result in results:
            replicate = getattr(result, "replicate", None)
            stored_fingerprint = getattr(result, "settings_fingerprint", None)

            if stored_fingerprint is None:
                if replicate is not None:
                    missing_fingerprint_replicates.append(replicate)
                continue

            if stored_fingerprint != expected_fingerprint:
                mismatched_fingerprints.append(
                    f"replicate {replicate}: stored={stored_fingerprint} current={expected_fingerprint}"
                )

        if missing_fingerprint_replicates:
            raise ValueError(
                f"Catalytic-triad aggregation for condition '{ctx.condition.label}' cannot use "
                "legacy cached replicate results missing settings fingerprints. Affected "
                f"replicates: {sorted(missing_fingerprint_replicates)}. Recompute the "
                "condition to refresh settings-sensitive caches before aggregating."
            )

        if mismatched_fingerprints:
            mismatch_text = "; ".join(mismatched_fingerprints)
            raise ValueError(
                f"Catalytic-triad aggregation for condition '{ctx.condition.label}' detected "
                f"settings fingerprint mismatches ({mismatch_text}). Recompute the condition "
                "or clear stale caches before aggregating."
            )

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return the triad loader class for the shared runner seam.

        Returns
        -------
        type[Any]
            Loader class patched by catalytic-triad unit tests.
        """

        return TrajectoryLoader

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any:
        """Build the runner-backed catalytic-triad execution object."""

        del replicate
        return CatalyticTriadReplicateRunner(
            universe=universe,
            pair_selections=ctx.settings.get_pair_selections(),
            threshold=ctx.settings.threshold,
            timestep_ps=window.timestep_ps,
            pair_label_func=_make_pair_label,
        )

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Serialize runner output into the legacy triad result schema."""

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

        settings = ctx.settings
        settings_tag = self._settings_cache_tag(settings)
        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        config_hash = compute_config_hash(ctx.sim_config)
        payload = runner.results
        metadata = DistanceResultMetadata(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
        )

        pair_results: list[DistancePairResult] = []
        pair_distance_arrays: list[np.ndarray] = []
        for pair_setting, pair_payload in zip(settings.pairs, payload.pair_payloads, strict=True):
            pair_results.append(
                DistancePairResult.from_runner_payload(
                    metadata,
                    pair_payload,
                    pair_label=pair_setting.label,
                )
            )
            pair_distance_arrays.append(np.asarray(pair_payload.distances, dtype=np.float64))

        all_below = np.ones(payload.n_frames_used, dtype=bool)
        for distance_array in pair_distance_arrays:
            all_below &= distance_array < settings.threshold

        simultaneous_fraction = float(all_below.mean())
        n_frames_simultaneous = int(all_below.sum())
        logger.info(
            "  Simultaneous contact: %.1f%% (%d/%d frames)",
            simultaneous_fraction * 100.0,
            n_frames_simultaneous,
            payload.n_frames_used,
        )

        contact_timeseries = all_below.astype(np.float64)
        sim_contact_sem: float | None = None
        sim_contact_tau: float | None = None
        sim_contact_tau_unit: str | None = None
        sim_contact_n_ind: int | None = None
        sim_contact_warning: str | None = None
        if payload.n_frames_used >= 20:
            try:
                tau_result = estimate_correlation_time(
                    contact_timeseries,
                    timestep=window.timestep_ps * window.step,
                    timestep_unit="ps",
                    method="integration",
                    n_frames=payload.n_frames_used,
                )
                sim_contact_tau = tau_result.tau
                sim_contact_tau_unit = tau_result.tau_unit
                sim_contact_n_ind = tau_result.n_independent
                sim_contact_warning = tau_result.warning
                if sim_contact_n_ind and sim_contact_n_ind > 0:
                    probability = simultaneous_fraction
                    sim_contact_sem = float(
                        np.sqrt(probability * (1.0 - probability) / float(sim_contact_n_ind))
                    )
            except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                logger.warning(
                    "Autocorrelation analysis for contact timeseries failed: %s",
                    exc,
                )
                probability = simultaneous_fraction
                sim_contact_sem = float(
                    np.sqrt(probability * (1.0 - probability) / float(payload.n_frames_used))
                )

        return TriadResult(
            config_hash=metadata.config_hash,
            polyzymd_version=metadata.polyzymd_version,
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string="; ".join(
                f"({pair.selection_a} : {pair.selection_b})" for pair in settings.pairs
            ),
            triad_name=settings.name,
            triad_description=settings.description,
            pair_results=pair_results,
            threshold=settings.threshold,
            simultaneous_contact_fraction=simultaneous_fraction,
            n_frames_simultaneous=n_frames_simultaneous,
            simultaneous_contact_timeseries=None,
            sim_contact_sem=sim_contact_sem,
            sim_contact_correlation_time=sim_contact_tau,
            sim_contact_correlation_time_unit=sim_contact_tau_unit,
            sim_contact_n_independent=sim_contact_n_ind,
            sim_contact_warning=sim_contact_warning,
            n_frames_total=payload.n_frames_total,
            n_frames_used=payload.n_frames_used,
            settings_fingerprint=settings_tag,
        )

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate triad results across replicates for one condition.

        Computes per-pair aggregated distance statistics and overall
        simultaneous contact fraction mean +/- SEM from the already-computed
        per-replicate results. Does NOT re-run the analyzer.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[TriadResult]
            Per-replicate triad results.

        Returns
        -------
        TriadAggregatedResult
            Aggregated result.
        """
        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.catalytic_triad._results import TriadAggregatedResult
        from polyzymd.analyses.shared.aggregation import aggregate_distance_pair_stats
        from polyzymd.analyses.shared.statistics import compute_sem

        settings = ctx.settings
        self._validate_replicate_result_settings_identity(ctx, results)
        self._validate_replicate_pair_schema(ctx, results, settings)
        first = results[0]
        settings_tag = self._settings_cache_tag(settings)
        metadata = DistanceResultMetadata(
            config_hash=first.config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=first.equilibration_time,
            equilibration_unit=first.equilibration_unit,
        )

        # Aggregate per-pair distance statistics
        n_pairs = len(first.pair_results)
        aggregated_pairs: list[DistancePairAggregatedResult] = []

        for pair_idx in range(n_pairs):
            stats = aggregate_distance_pair_stats(list(results), pair_idx)

            # Get pair config from the first result's pair_results
            pr = first.pair_results[pair_idx]

            agg_pair = DistancePairAggregatedResult.from_aggregated_stats(
                metadata,
                pr,
                stats,
                replicates=ctx.replicates,
                threshold=settings.threshold,
                pair_label=pr.pair_label,
            )
            aggregated_pairs.append(agg_pair)

        # Aggregate simultaneous contact fraction
        per_rep_simultaneous = [r.simultaneous_contact_fraction for r in results]
        sim_stats = compute_sem(per_rep_simultaneous)

        triad_name = settings.name
        triad_description = settings.description

        agg_result = TriadAggregatedResult(
            config_hash=first.config_hash,
            polyzymd_version=metadata.polyzymd_version,
            replicate=None,
            equilibration_time=first.equilibration_time,
            equilibration_unit=first.equilibration_unit,
            selection_string="; ".join(
                f"({pr.selection1} : {pr.selection2})" for pr in first.pair_results
            ),
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            triad_name=triad_name,
            triad_description=triad_description,
            pair_results=aggregated_pairs,
            threshold=settings.threshold,
            overall_simultaneous_contact=sim_stats.mean,
            sem_simultaneous_contact=sim_stats.sem,
            per_replicate_simultaneous=per_rep_simultaneous,
            source_result_files=[],  # Not tracked in plugin mode
            settings_fingerprint=settings_tag,
        )

        target_path = ctx.result_path
        if target_path is None:
            filename = self._make_aggregated_filename(ctx.replicates, first)
            target_path = ctx.output_dir / filename
        self.save_result(agg_result, target_path)
        logger.info(f"Saved aggregated triad result to {target_path}")

        return agg_result

    @staticmethod
    def _validate_replicate_pair_schema(
        ctx: AggregateContext,
        results: Sequence[TriadResult],
        settings: CatalyticTriadSettings,
    ) -> None:
        """Require all replicate triad pair results to match the configured schema.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[TriadResult]
            Per-replicate catalytic-triad results.
        settings : CatalyticTriadSettings
            Current triad settings that define the canonical ordered pair schema.

        Raises
        ------
        ValueError
            Raised when any replicate result has a mismatched threshold, pair
            count, ordering, labels, or selections.
        """

        schema_issues: list[str] = []

        for result in results:
            replicate = getattr(result, "replicate", None)
            pair_results = list(result.pair_results)

            if result.threshold != settings.threshold:
                schema_issues.append(
                    f"replicate {replicate}: threshold {result.threshold!r} != {settings.threshold!r}"
                )

            if len(pair_results) != len(settings.pairs):
                schema_issues.append(
                    f"replicate {replicate}: pair count {len(pair_results)} != {len(settings.pairs)}"
                )
                continue

            for pair_idx, (pair_result, pair_setting) in enumerate(
                zip(pair_results, settings.pairs, strict=True),
                start=1,
            ):
                mismatch_parts: list[str] = []
                if pair_result.pair_label != pair_setting.label:
                    mismatch_parts.append(
                        f"label {pair_result.pair_label!r} != {pair_setting.label!r}"
                    )
                if pair_result.selection1 != pair_setting.selection_a:
                    mismatch_parts.append(
                        f"selection1 {pair_result.selection1!r} != {pair_setting.selection_a!r}"
                    )
                if pair_result.selection2 != pair_setting.selection_b:
                    mismatch_parts.append(
                        f"selection2 {pair_result.selection2!r} != {pair_setting.selection_b!r}"
                    )
                if pair_result.threshold != settings.threshold:
                    mismatch_parts.append(
                        f"pair threshold {pair_result.threshold!r} != {settings.threshold!r}"
                    )

                if mismatch_parts:
                    schema_issues.append(
                        f"replicate {replicate} pair {pair_idx}: {', '.join(mismatch_parts)}"
                    )

        if schema_issues:
            issue_text = "; ".join(schema_issues)
            raise ValueError(
                f"Catalytic-triad aggregation for condition '{ctx.condition.label}' requires "
                f"every replicate result to match the configured pair schema. Problems "
                f"detected: {issue_text}."
            )

    # === Optional methods ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract simultaneous contact fraction for scalar comparison.

        Parameters
        ----------
        summary : TriadAggregatedResult
            Aggregated triad result.

        Returns
        -------
        dict[str, MetricValue]
            Single entry ``"simultaneous_contact_fraction"`` with
            ``higher_is_better=True`` (more contact = better integrity).
        """
        return {
            "simultaneous_contact_fraction": MetricValue(
                name="simultaneous_contact_fraction",
                mean=summary.overall_simultaneous_contact * 100,
                sem=summary.sem_simultaneous_contact * 100,
                replicate_values=[v * 100 for v in summary.per_replicate_simultaneous],
                higher_is_better=True,
                direction_labels=("worsening", "unchanged", "improving"),
            ),
        }

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format catalytic triad comparison result for CLI display.

        Parameters
        ----------
        result : ComparisonResult or BaseModel
            Comparison result to format.
        output_format : str
            ``"text"``, ``"markdown"``, or ``"json"``.

        Returns
        -------
        str
            Formatted output.
        """
        from polyzymd.analyses.base import ComparisonResult
        from polyzymd.analyses.stats import format_scalar_comparison

        if isinstance(result, ComparisonResult):
            return format_scalar_comparison(
                result,
                title="Catalytic Triad Comparison",
                metric_label="Simultaneous Contact",
                metric_unit="%",
                metric_key="simultaneous_contact_fraction",
                output_format=output_format,
                higher_is_better=True,
            )
        return super().format(result, output_format)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate triad comparison plots.

        Delegates to private module-level helpers that call the standalone
        plotting functions in ``_plotters``.

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

        # Resolve plot settings (guaranteed non-None by orchestrator)
        plot_settings = ctx.plot_settings

        # KDE panel plot
        result = plot_triad_kde_panel_from_data(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        # Threshold bars plot
        result = plot_triad_threshold_bars_from_data(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        return plots

    # === Private helpers ===

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Backward-compatible filename helper retained for tests."""
        eq_str = f"eq{first_result.equilibration_time:g}{first_result.equilibration_unit}"
        rep_str = Analysis._format_replicate_range(replicates)
        name_safe = first_result.triad_name.replace(" ", "_").replace("/", "-")
        return f"triad_{name_safe}_{rep_str}_{eq_str}.json"


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def _make_result_filename(settings: CatalyticTriadSettings, eq_value: float, eq_unit: str) -> str:
    """Generate filename for single-replicate result JSON."""
    eq_str = f"eq{eq_value:g}{eq_unit}"
    name_safe = settings.name.replace(" ", "_").replace("/", "-")
    settings_tag = settings_fingerprint(settings)
    return f"triad_{name_safe}_{eq_str}_s{settings_tag}.json"
