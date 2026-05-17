"""Catalytic triad analysis plugin.

Analyses active-site geometry from MD trajectories by computing per-pair
distances and a simultaneous contact fraction (all pairs below threshold in the
same frame). Trajectory-native pair distances run through the shared MDAnalysis
pair-distance primitive, while the plugin owns the triad-specific composite
reducer and artifact adapters.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ANOVAResult,
    BasePlotSettings,
    ComparisonResult,
    ConditionSummary,
    MetricValue,
    PairwiseResult,
    PlotContext,
)
from polyzymd.analyses.catalytic_triad._mda import (
    SIMULTANEOUS_CONTACT_METADATA,
    SIMULTANEOUS_CONTACT_METRIC,
    TriadArtifactCollector,
    aggregate_triad_artifacts,
    build_triad_jobs,
)
from polyzymd.analyses.catalytic_triad._plot_settings import TriadPlotSettings
from polyzymd.analyses.catalytic_triad._plotters import (
    plot_triad_kde_panel_from_data,
    plot_triad_threshold_bars_from_data,
)
from polyzymd.analyses.catalytic_triad._results import TriadAggregatedResult
from polyzymd.analyses.mda import (
    ArtifactStore,
    ComparisonArtifact,
    MDACollectorContext,
    MDAReplicateJobContext,
    ReplicateArtifact,
)
from polyzymd.analyses.shared.config_hash import settings_fingerprint

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

    Computes per-pair distances through the shared MDAnalysis pair-distance job
    and derives a simultaneous contact fraction — the percentage of frames where
    all pairs are below the threshold at the same time.

    The ``compare()`` method is NOT overridden. Canonical condition artifacts
    are compared through the MDAnalysis artifact lifecycle, while legacy
    in-memory summaries still expose ``simultaneous_contact_fraction`` as a
    percent-scaled scalar metric for the default comparison fallback.

    Plots
    -----
    - ``triad_kde_panel.png``: Multi-row KDE of per-pair distance distributions.
    - ``triad_threshold_bars.png``: Grouped bar chart of contact fractions.
    """

    name: ClassVar[str] = "catalytic_triad"
    Settings: ClassVar[type] = CatalyticTriadSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = TriadPlotSettings
    AggregatedResultClass: ClassVar[type] = TriadAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = None
    aliases: ClassVar[tuple[str, ...]] = ("triad",)
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 1

    # === Required methods ===

    @staticmethod
    def _settings_cache_tag(settings: BaseModel) -> str:
        """Return the catalytic-triad settings fingerprint."""

        return settings_fingerprint(settings)

    def aggregate_settings_fingerprint(self, settings: BaseModel | None) -> str | None:
        """Return the legacy catalytic-triad aggregate settings fingerprint.

        Parameters
        ----------
        settings : BaseModel or None
            Active catalytic-triad settings.

        Returns
        -------
        str or None
            Existing settings-only fingerprint used by saved aggregate results.
        """

        if settings is None:
            return None
        return self._settings_cache_tag(settings)

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[Any] | None:
        """Build MDAnalysis-native pair-distance jobs for one triad replicate."""

        return build_triad_jobs(ctx, ctx.settings)

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the catalytic-triad artifact collector."""

        del ctx
        return TriadArtifactCollector()

    def run_replicate(self, ctx: Any, replicate: int) -> Any:
        """Run one catalytic-triad replicate through the MDA artifact path.

        Parameters
        ----------
        ctx : Any
            Framework-provided replicate context.
        replicate : int
            One-indexed replicate number.

        Returns
        -------
        Any
            Canonical replicate artifact produced by the MDAnalysis lifecycle.
        """

        result = Analysis.run_replicate(self, ctx, replicate)
        if isinstance(result, ReplicateArtifact):
            target_path = ctx.result_path or ctx.output_dir / "result.json"
            ArtifactStore(target_path.parent).write_replicate_result(result, target_path.name)
        return result

    def _deserialize_result(self, path: Path) -> Any:
        """Load only canonical condition artifacts for catalytic-triad aggregates."""

        if path.exists():
            try:
                loaded = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(
                    f"Catalytic-triad aggregate at {path} is not a valid canonical artifact. "
                    "Recompute the condition or clear stale caches before comparing."
                ) from exc
            if isinstance(loaded, dict) and loaded.get("artifact_type") == "condition":
                return ArtifactStore(path.parent).read_condition_result(path.name)
        raise ValueError(
            f"Catalytic-triad aggregate at {path} is not a canonical MDAnalysis condition "
            "artifact. Recompute the condition or clear stale legacy triad caches."
        )

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate triad replicate artifacts across one condition.

        Computes per-pair aggregated distance statistics and overall
        simultaneous contact fraction mean +/- SEM from canonical replicate
        artifacts. Does not re-run trajectory analysis.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[ReplicateArtifact]
            Per-replicate MDAnalysis triad artifacts.

        Returns
        -------
        ConditionArtifact
            Aggregated condition artifact with legacy-compatible summaries.
        """
        if not results:
            raise ValueError(
                f"Catalytic-triad aggregation for condition '{ctx.condition.label}' requires at "
                "least one replicate artifact. No replicate inputs were provided."
            )
        if not all(isinstance(result, ReplicateArtifact) for result in results):
            raise TypeError(
                "Catalytic-triad aggregation expects MDAnalysis ReplicateArtifact inputs. Legacy "
                "triad replicate caches are incompatible with the MDAnalysis artifact lifecycle; "
                "recompute the condition or clear stale caches before aggregating."
            )
        target_path = ctx.result_path or ctx.output_dir / "result.json"
        aggregated = aggregate_triad_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            output_dir=ctx.output_dir,
            result_path=target_path,
            artifacts=results,
            settings_fingerprint=self._settings_cache_tag(ctx.settings),
        )
        logger.info("Saved aggregated catalytic-triad artifact to %s", target_path)
        return aggregated

    # === Optional methods ===

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
        from polyzymd.analyses.stats import format_scalar_comparison

        comparison_result = result
        if isinstance(result, ComparisonArtifact):
            if output_format == "json":
                return result.model_dump_json(indent=2)
            comparison_result = _comparison_artifact_to_result(result)

        if isinstance(comparison_result, ComparisonResult):
            return format_scalar_comparison(
                comparison_result,
                title="Catalytic Triad Comparison",
                metric_label="Simultaneous Contact",
                metric_unit="%",
                metric_key="simultaneous_contact_fraction",
                output_format=output_format,
                higher_is_better=True,
            )
        return super().format(result, output_format)

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract the percent-scaled simultaneous contact metric.

        Parameters
        ----------
        summary : Any
            Aggregated catalytic-triad summary with simultaneous contact fields.

        Returns
        -------
        dict[str, MetricValue]
            Mapping for the primary simultaneous contact metric.
        """

        return {
            SIMULTANEOUS_CONTACT_METRIC: MetricValue(
                name=SIMULTANEOUS_CONTACT_METRIC,
                mean=float(summary.overall_simultaneous_contact) * 100.0,
                sem=float(summary.sem_simultaneous_contact) * 100.0,
                replicate_values=[
                    float(value) * 100.0 for value in summary.per_replicate_simultaneous
                ],
                higher_is_better=bool(SIMULTANEOUS_CONTACT_METADATA["higher_is_better"]),
                direction_labels=tuple(SIMULTANEOUS_CONTACT_METADATA["direction_labels"]),
            )
        }

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


def _comparison_artifact_to_result(artifact: ComparisonArtifact) -> ComparisonResult:
    """Adapt a canonical comparison artifact for scalar CLI formatting.

    Parameters
    ----------
    artifact : ComparisonArtifact
        MDAnalysis comparison artifact generated from condition artifacts.

    Returns
    -------
    ComparisonResult
        Legacy scalar comparison model consumed by the shared formatter.
    """

    payload = artifact.payload
    statistical_parameters = payload.get("statistical_parameters", {})
    if not isinstance(statistical_parameters, dict):
        statistical_parameters = {}
    anova_payload = payload.get("anova") or []
    return ComparisonResult(
        analysis_type=artifact.analysis_name,
        name=str(artifact.metadata.get("project_name", artifact.analysis_name)),
        control_label=artifact.effective_control or artifact.control_label,
        fdr_alpha=statistical_parameters.get("fdr_alpha"),
        ttest_method=str(statistical_parameters.get("ttest_method", "student")),
        posthoc_method=str(statistical_parameters.get("posthoc_method", "ttest_bh")),
        conditions=[
            ConditionSummary.model_validate(summary)
            for summary in payload.get("condition_summaries", [])
        ],
        pairwise_comparisons=[
            PairwiseResult.model_validate(pairwise)
            for pairwise in payload.get("pairwise_comparisons", [])
        ],
        anova=[ANOVAResult.model_validate(entry) for entry in anova_payload] or None,
        ranking=list(payload.get("ranking", [])),
        rankings_by_metric=payload.get("rankings_by_metric"),
        equilibration_time=str(statistical_parameters.get("equilibration", "0ns")),
        created_at=str(artifact.metadata.get("created_at", "")),
        polyzymd_version=str(artifact.metadata.get("polyzymd_version", "")),
    )
