"""Secondary structure analysis plugin.

Computes per-residue DSSP categorical states through the MDAnalysis extension
layer, stores per-replicate state matrices in NPZ sidecars, aggregates state
occupancy across replicates, and compares conditions with the scalar artifact
comparison path.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses._framework.cache_identity import settings_fingerprint
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    MetricValue,
    PlotContext,
    SlurmResourceHint,
)
from polyzymd.analyses.mda import (
    ArtifactStore,
    ComparisonArtifact,
    ConditionArtifact,
    MDACollectorContext,
    ReplicateArtifact,
)
from polyzymd.analyses.secondary_structure._mda import (
    HELIX_FRACTION_METRIC,
    SecondaryStructureArtifactCollector,
    aggregate_secondary_structure_artifacts,
    build_secondary_structure_jobs,
)
from polyzymd.analyses.secondary_structure._plot_settings import SSPlotSettings
from polyzymd.analyses.secondary_structure._plotters import (
    _plot_ss_content_bars,
    _plot_ss_individual_bars,
    _plot_ss_persistence_diff_heatmap,
    _plot_ss_timeline_heatmap,
)

if TYPE_CHECKING:
    from polyzymd.analyses.mda import MDAAnalysisJob, MDAReplicateJobContext


def _legacy_secondary_structure_fingerprint(chain_id: str) -> str:
    """Return the legacy settings fingerprint for chain-only settings.

    Parameters
    ----------
    chain_id : str
        Effective protein chain ID.

    Returns
    -------
    str
        First eight hexadecimal characters of the SHA-256 digest.
    """

    serialized = json.dumps({"chain_id": chain_id}, sort_keys=True)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:8]


class SecondaryStructureSettings(BaseModel):
    """Settings for secondary structure analysis.

    Chain ``A`` is the default protein chain under the PolyzyMD chain
    convention. Set ``selection`` to an explicit MDAnalysis selection when
    topology files do not preserve chain IDs.
    """

    chain_id: str = Field(
        default="A",
        description="Chain letter for the protein chain",
    )
    selection: str | None = Field(
        default=None,
        description="Explicit MDAnalysis protein selection. Overrides chain_id when provided.",
    )

    @field_validator("chain_id")
    @classmethod
    def _validate_chain_id(cls, value: str) -> str:
        """Validate and normalize the protein chain ID.

        Parameters
        ----------
        value : str
            User-provided protein chain ID.

        Returns
        -------
        str
            Stripped chain ID.
        """

        stripped = value.strip()
        if not stripped:
            raise ValueError("Secondary-structure chain_id must not be empty")
        return stripped

    @field_validator("selection")
    @classmethod
    def _validate_selection(cls, value: str | None) -> str | None:
        """Validate and normalize an explicit protein selection.

        Parameters
        ----------
        value : str or None
            Optional MDAnalysis selection string.

        Returns
        -------
        str or None
            Stripped selection string, or ``None`` when not configured.
        """

        if value is None:
            return None
        stripped = value.strip()
        if not stripped:
            raise ValueError("Secondary-structure selection must not be empty when provided")
        return stripped


class SecondaryStructureAnalysis(Analysis):
    """DSSP secondary structure analysis using canonical MDA artifacts."""

    name: ClassVar[str] = "secondary_structure"
    Settings: ClassVar[type] = SecondaryStructureSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = SSPlotSettings
    AggregatedResultClass: ClassVar[type | None] = None
    ReplicateResultClass: ClassVar[type | None] = None
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 1
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = SlurmResourceHint(mem="16G")

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel) -> str:
        """Return the short settings fingerprint used by artifacts.

        Parameters
        ----------
        settings : BaseModel
            Analysis settings model.

        Returns
        -------
        str
            Stable settings fingerprint.
        """

        if isinstance(settings, SecondaryStructureSettings) and settings.selection is None:
            return _legacy_secondary_structure_fingerprint(settings.chain_id)
        return settings_fingerprint(settings)

    def aggregate_settings_fingerprint(self, settings: BaseModel | None) -> str | None:
        """Return the secondary-structure artifact settings fingerprint."""

        if settings is None:
            return None
        return self._make_settings_cache_tag(settings)

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[MDAAnalysisJob] | None:
        """Build the AnalysisBase-compatible DSSP job.

        Parameters
        ----------
        ctx : MDAReplicateJobContext
            Framework-provided MDAnalysis job context.

        Returns
        -------
        sequence of MDAAnalysisJob
            A single DSSP job.
        """

        return build_secondary_structure_jobs(ctx, ctx.settings)

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the secondary-structure artifact collector."""

        del ctx
        return SecondaryStructureArtifactCollector()

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate secondary-structure replicate artifacts across one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregate context.
        results : sequence of Any
            Per-replicate canonical artifacts.

        Returns
        -------
        ConditionArtifact
            Aggregated condition artifact.
        """

        if not results:
            raise ValueError(
                f"Secondary-structure aggregation for condition '{ctx.condition.label}' requires "
                "at least one replicate artifact. No replicate inputs were provided."
            )
        if not all(isinstance(result, ReplicateArtifact) for result in results):
            raise TypeError(
                "Secondary-structure aggregation expects MDAnalysis ReplicateArtifact inputs. "
                "Non-canonical secondary-structure replicate caches are incompatible with the "
                "MDAnalysis artifact lifecycle; recompute the condition or clear stale caches."
            )
        return aggregate_secondary_structure_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            output_dir=ctx.output_dir,
            artifacts=results,
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
        )

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract helix fraction directly from a condition artifact.

        Parameters
        ----------
        summary : Any
            Canonical condition artifact.

        Returns
        -------
        dict[str, MetricValue]
            Single ``helix_fraction`` metric with higher values treated as more
            structured.
        """

        if not isinstance(summary, ConditionArtifact):
            raise TypeError(
                "Secondary-structure metric extraction requires a canonical MDAnalysis "
                "condition artifact. Non-canonical aggregate summaries are incompatible; recompute "
                "the condition or clear stale caches before comparing."
            )
        metric = summary.payload.get("metrics", {}).get(HELIX_FRACTION_METRIC)
        if not isinstance(metric, dict):
            raise ValueError(
                "Secondary-structure condition artifact is missing the helix_fraction metric "
                "payload. Recompute the condition or clear stale caches before comparing."
            )
        return {
            HELIX_FRACTION_METRIC: MetricValue(
                name=HELIX_FRACTION_METRIC,
                mean=float(metric["mean"]),
                sem=float(metric["sem"]),
                replicate_values=[float(value) for value in metric["values"]],
                higher_is_better=bool(metric.get("higher_is_better", True)),
                direction_labels=tuple(
                    metric.get("direction_labels", ("destabilizing", "unchanged", "stabilizing"))
                ),
            ),
        }

    def compare(self, ctx: Any) -> Any:
        """Compare canonical condition artifacts and reject condition aggregates."""

        for label, summary in ctx.aggregated_results.items():
            if summary is not None and not isinstance(summary, ConditionArtifact):
                raise TypeError(
                    f"Secondary-structure comparison for condition '{label}' requires canonical "
                    "MDAnalysis condition artifacts. Non-canonical aggregate inputs are incompatible; "
                    "recompute the condition or clear stale caches before comparing."
                )
        return super().compare(ctx)

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format secondary-structure comparison output for CLI display."""

        from polyzymd.analyses.stats import format_scalar_comparison_artifact_payload

        if isinstance(result, ComparisonArtifact):
            if output_format == "json":
                return result.model_dump_json(indent=2)
            return format_scalar_comparison_artifact_payload(
                result.payload,
                title="Secondary Structure Comparison",
                metric_label="Helix Fraction",
                metric_unit="",
                metric_key=HELIX_FRACTION_METRIC,
                output_format=output_format,
                higher_is_better=True,
            )
        return super().format(result, output_format)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate secondary-structure plots from canonical artifacts only.

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided plot context.

        Returns
        -------
        list[Path]
            Paths to generated figure files.
        """

        data, labels = self._build_plot_data(ctx, include_replicates=True)
        if not labels:
            return []

        condition_by_label = {condition.label: condition for condition in ctx.conditions}
        for label in labels:
            cond_data = data.get(label)
            condition = condition_by_label.get(label)
            if cond_data is None or condition is None:
                continue
            agg_dir = Path(cond_data["aggregated_dir"])
            artifact = self._load_aggregated_result(agg_dir)
            if artifact is None:
                continue
            artifact = self.validate_aggregated_result(
                artifact,
                condition=condition,
                settings=ctx.settings,
                equilibration=ctx.equilibration,
                source=self.aggregate_result_path(agg_dir),
                expected_replicates=condition.replicates,
                allow_replicate_subset=True,
            )
            cond_data["condition_artifact"] = artifact

        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        plots: list[Path] = []
        plots.extend(_plot_ss_timeline_heatmap(data, labels, ctx.output_dir, ctx.plot_settings))
        plots.extend(_plot_ss_content_bars(data, labels, ctx.output_dir, ctx.plot_settings))
        plots.extend(_plot_ss_individual_bars(data, labels, ctx.output_dir, ctx.plot_settings))
        plots.extend(
            _plot_ss_persistence_diff_heatmap(data, labels, ctx.output_dir, ctx.plot_settings)
        )
        return plots

    def _deserialize_result(self, path: Path) -> Any:
        """Load only canonical secondary-structure condition artifacts."""

        if path.exists():
            try:
                loaded = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(
                    f"Secondary-structure aggregate at {path} is not a valid canonical artifact. "
                    "Recompute the condition or clear stale caches before comparing."
                ) from exc
            if isinstance(loaded, dict) and loaded.get("artifact_type") == "condition":
                return ArtifactStore(path.parent).read_condition_result(path.name)
        raise ValueError(
            f"Secondary-structure aggregate at {path} is not a canonical MDAnalysis condition "
            "artifact. Recompute the condition or clear stale non-canonical caches."
        )
