"""RMSF analysis plugin backed by MDAnalysis-compatible profile jobs."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    ComparisonResult,
    MetricValue,
    PlotContext,
)
from polyzymd.analyses.mda import (
    ArtifactStore,
    ComparisonArtifact,
    ConditionArtifact,
    ReplicateArtifact,
)
from polyzymd.analyses.rmsf._mda import (
    MEAN_RMSF_METRIC,
    RMSFArtifactCollector,
    aggregate_rmsf_artifacts,
    build_rmsf_jobs,
    condition_artifact_to_legacy_result,
    external_reference_file_identity,
)
from polyzymd.analyses.rmsf._plot_settings import RMSFPlotSettings
from polyzymd.analyses.rmsf._plotters import _plot_rmsf_comparison, _plot_rmsf_profile
from polyzymd.analyses.rmsf._results import RMSFAggregatedResult
from polyzymd.analyses.shared.config_hash import settings_fingerprint

if TYPE_CHECKING:
    from polyzymd.analyses.mda import MDACollectorContext, MDAReplicateJobContext

logger = logging.getLogger(__name__)


class RMSFSettings(BaseModel):
    """Settings for RMSF profile analysis.

    RMSF is computed as a per-residue fluctuation profile over the selected
    atoms after reference-based alignment.
    """

    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection string for RMSF calculation",
    )
    reference_mode: str = Field(
        default="centroid",
        description="Reference structure mode: centroid, average, frame, or external",
    )
    reference_frame: int | None = Field(
        default=None,
        description="Frame number if reference_mode is 'frame' (1-indexed)",
    )
    reference_file: str | None = Field(
        default=None,
        description="Path to external PDB file if reference_mode is 'external'",
    )
    alignment_selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for trajectory alignment",
    )
    centroid_selection: str = Field(
        default="protein",
        description="MDAnalysis selection for centroid finding",
    )

    @field_validator("reference_mode", mode="after")
    @classmethod
    def validate_reference_mode(cls, value: str) -> str:
        """Validate the reference mode setting."""

        valid = {"centroid", "average", "frame", "external"}
        if value not in valid:
            raise ValueError(f"reference_mode must be one of {valid}, got {value!r}")
        return value


class RMSFAnalysis(Analysis):
    """Per-residue RMSF analysis using canonical MDA artifacts."""

    name: ClassVar[str] = "rmsf"
    Settings: ClassVar[type] = RMSFSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = RMSFPlotSettings
    AggregatedResultClass: ClassVar[type] = RMSFAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = None
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 1

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel | Any) -> str:
        """Return a short settings fingerprint for RMSF artifacts.

        Parameters
        ----------
        settings : BaseModel or Any
            RMSF settings object.

        Returns
        -------
        str
            Stable short settings fingerprint.
        """

        if isinstance(settings, RMSFSettings):
            normalized = settings
        elif isinstance(settings, BaseModel):
            normalized = RMSFSettings.model_validate(settings.model_dump(mode="json"))
        else:
            normalized = RMSFSettings(
                selection=settings.selection,
                reference_mode=settings.reference_mode,
                reference_frame=settings.reference_frame,
                reference_file=settings.reference_file,
                alignment_selection=settings.alignment_selection,
                centroid_selection=settings.centroid_selection,
            )
        if normalized.reference_mode != "external":
            return settings_fingerprint(normalized)
        payload = normalized.model_dump(mode="json")
        payload["reference_file_identity"] = external_reference_file_identity(
            normalized.reference_file
        )
        serialized = json.dumps(payload, sort_keys=True)
        return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:8]

    def aggregate_settings_fingerprint(self, settings: BaseModel | None) -> str | None:
        """Return the RMSF artifact settings fingerprint."""

        if settings is None:
            return None
        return self._make_settings_cache_tag(settings)

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[Any] | None:
        """Build the RMSF MDAnalysis-compatible profile job."""

        return build_rmsf_jobs(ctx, ctx.settings)

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the RMSF artifact collector."""

        del ctx
        return RMSFArtifactCollector()

    def run_replicate(self, ctx: Any, replicate: int) -> Any:
        """Run one RMSF replicate through the MDA artifact lifecycle.

        Parameters
        ----------
        ctx : Any
            Framework-provided replicate context.
        replicate : int
            One-indexed replicate number.

        Returns
        -------
        Any
            Canonical replicate artifact.
        """

        result = Analysis.run_replicate(self, ctx, replicate)
        if isinstance(result, ReplicateArtifact):
            target_path = ctx.result_path or ctx.output_dir / "result.json"
            ArtifactStore(target_path.parent).write_replicate_result(result, target_path.name)
        return result

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate RMSF replicate artifacts across one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : sequence of Any
            Per-replicate canonical RMSF artifacts.

        Returns
        -------
        ConditionArtifact
            Aggregated condition artifact.
        """

        if not results:
            raise ValueError(
                f"RMSF aggregation for condition '{ctx.condition.label}' requires at least one "
                "replicate artifact. No replicate inputs were provided."
            )
        if not all(isinstance(result, ReplicateArtifact) for result in results):
            raise TypeError(
                "RMSF aggregation expects MDAnalysis ReplicateArtifact inputs. Legacy RMSF "
                "replicate caches are incompatible with the MDAnalysis artifact lifecycle; "
                "recompute the condition or clear stale caches before aggregating."
            )
        target_path = ctx.result_path or ctx.output_dir / "result.json"
        aggregated = aggregate_rmsf_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            output_dir=ctx.output_dir,
            result_path=target_path,
            artifacts=results,
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
        )
        logger.info("Saved aggregated RMSF artifact to %s", target_path)
        return aggregated

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract mean RMSF from a legacy aggregate adapter.

        Parameters
        ----------
        summary : Any
            Aggregated RMSF result or condition artifact.

        Returns
        -------
        dict[str, MetricValue]
            Mean RMSF metric for the scalar comparison path.
        """

        if isinstance(summary, ConditionArtifact):
            summary = condition_artifact_to_legacy_result(summary)
        return {
            MEAN_RMSF_METRIC: MetricValue(
                name=MEAN_RMSF_METRIC,
                mean=summary.overall_mean_rmsf,
                sem=summary.overall_sem_rmsf,
                replicate_values=summary.per_replicate_mean_rmsf,
                higher_is_better=False,
                direction_labels=("stabilizing", "unchanged", "destabilizing"),
            )
        }

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare RMSF condition artifacts with legacy-cache rejection."""

        for label, summary in ctx.aggregated_results.items():
            if summary is not None and not isinstance(summary, ConditionArtifact):
                raise TypeError(
                    f"RMSF comparison for condition '{label}' requires canonical MDAnalysis "
                    "condition artifacts. Legacy RMSF aggregate inputs are incompatible; "
                    "recompute the condition or clear stale caches before comparing."
                )
        return super().compare(ctx)

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format RMSF comparison output for CLI display."""

        from polyzymd.analyses.stats import format_scalar_comparison

        comparison_result = result
        if isinstance(result, ComparisonArtifact):
            if output_format == "json":
                return result.model_dump_json(indent=2)
            comparison_result = _comparison_artifact_to_result(result)

        if isinstance(comparison_result, ComparisonResult):
            return format_scalar_comparison(
                comparison_result,
                title="RMSF Comparison",
                metric_label="Mean RMSF",
                metric_unit="A",
                metric_key=MEAN_RMSF_METRIC,
                output_format=output_format,
                higher_is_better=False,
            )
        return super().format(result, output_format)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate RMSF plots from cached artifacts and sidecars only."""

        data, labels = self._build_plot_data(ctx)
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
            cond_data["aggregated_result"] = condition_artifact_to_legacy_result(artifact)
            cond_data["condition_artifact"] = artifact

        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        plots: list[Path] = []
        plots.extend(_plot_rmsf_comparison(data, labels, ctx.output_dir, ctx.plot_settings))
        plots.extend(_plot_rmsf_profile(data, labels, ctx.output_dir, ctx.plot_settings))
        return plots

    def _deserialize_result(self, path: Path) -> Any:
        """Load only canonical RMSF condition artifacts for aggregate results."""

        if path.exists():
            try:
                loaded = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(
                    f"RMSF aggregate at {path} is not a valid canonical artifact. Recompute "
                    "the condition or clear stale caches before comparing."
                ) from exc
            if isinstance(loaded, dict) and loaded.get("artifact_type") == "condition":
                return ArtifactStore(path.parent).read_condition_result(path.name)
        raise ValueError(
            f"RMSF aggregate at {path} is not a canonical MDAnalysis condition artifact. "
            "Recompute the condition or clear stale legacy RMSF caches."
        )


def _comparison_artifact_to_result(artifact: ComparisonArtifact) -> ComparisonResult:
    """Convert a comparison artifact to the generic scalar comparison model."""

    payload = artifact.payload
    statistical_parameters = payload.get("statistical_parameters", {})
    return ComparisonResult(
        analysis_type=artifact.analysis_name,
        name=str(statistical_parameters.get("project_name", artifact.analysis_name)),
        control_label=artifact.effective_control,
        fdr_alpha=statistical_parameters.get("fdr_alpha"),
        ttest_method=statistical_parameters.get("ttest_method"),
        posthoc_method=statistical_parameters.get("posthoc_method"),
        conditions=payload.get("condition_summaries", []),
        pairwise_comparisons=payload.get("pairwise_comparisons", []),
        anova=payload.get("anova") or None,
        ranking=payload.get("ranking", []),
        rankings_by_metric=payload.get("rankings_by_metric"),
        equilibration_time=statistical_parameters.get("equilibration", "0ns"),
        created_at=str(artifact.metadata.get("created_at", "")),
        polyzymd_version=str(artifact.metadata.get("polyzymd_version", "")),
    )
