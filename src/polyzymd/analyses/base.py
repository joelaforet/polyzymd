"""Public facade for the PolyzyMD analysis plugin system.

Every analysis in PolyzyMD inherits from :class:`Analysis`. The framework
discovers subclasses automatically and owns replicate iteration, caching,
dependency ordering, comparison, plotting, and CLI wiring.

This module remains the stable public import surface. Implementation details
live in private framework modules so plugins and tests can keep importing all
public symbols from ``polyzymd.analyses.base``.
"""

from __future__ import annotations

import json
from abc import ABC
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel

from polyzymd.analyses._aggregate_validation import AggregateValidationError
from polyzymd.analyses._aggregate_validation import (
    aggregate_settings_fingerprint as _aggregate_settings_fingerprint_impl,
)
from polyzymd.analyses._aggregate_validation import (
    validate_aggregated_result as _validate_aggregated_result_impl,
)
from polyzymd.analyses._analysis_compare import default_compare as _default_compare
from polyzymd.analyses._analysis_contract import validate_analysis_subclass
from polyzymd.analyses._analysis_io import (
    aggregate_result_path as _aggregate_result_path,
)
from polyzymd.analyses._analysis_io import (
    build_plot_data as _build_plot_data_impl,
)
from polyzymd.analyses._analysis_io import (
    check_cache as _check_cache_impl,
)
from polyzymd.analyses._analysis_io import (
    comparison_result_path as _comparison_result_path,
)
from polyzymd.analyses._analysis_io import (
    deserialize_replicate_result,
    format_replicate_range,
    load_aggregated_result,
    load_replicate_result,
)
from polyzymd.analyses._analysis_io import (
    deserialize_result as _deserialize_result_impl,
)
from polyzymd.analyses._analysis_io import (
    figures_output_dir as _figures_output_dir,
)
from polyzymd.analyses._analysis_io import (
    replicate_result_path as _replicate_result_path,
)
from polyzymd.analyses._analysis_io import (
    resolve_output_dir as _resolve_output_dir,
)
from polyzymd.analyses._analysis_io import (
    save_result as _save_result_impl,
)
from polyzymd.analyses._analysis_runner import (
    run_replicate_default,
    run_replicate_via_mda_jobs,
    run_replicate_via_runner,
)
from polyzymd.analyses._comparison_models import (
    ANOVAResult,
    BaseComparisonResult,
    BaseConditionSummary,
    BasePlotSettings,
    ComparisonResult,
    ConditionSummary,
    MetricValue,
    PairwiseResult,
    SlurmResourceHint,
    TConditionSummary,
    TPairwiseResult,
)
from polyzymd.analyses._contexts import (
    AggregateContext,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
    _default_plot_settings,
)
from polyzymd.analyses._measurement import CacheIdentity, Measurement, MetricSpec, ScalarMeasurement
from polyzymd.analyses.exceptions import PluginContractError

if TYPE_CHECKING:
    from polyzymd.analyses.mda import (
        MDAAnalysisJob,
        MDAArtifactCollector,
        MDACollectorContext,
        MDAReplicateJobContext,
    )
    from polyzymd.config.schema import SimulationConfig

__all__ = [
    "ANOVAResult",
    "AggregateContext",
    "AggregateValidationError",
    "Analysis",
    "BaseComparisonResult",
    "BaseConditionSummary",
    "BasePlotSettings",
    "CacheIdentity",
    "ComparisonContext",
    "ComparisonResult",
    "Condition",
    "ConditionSummary",
    "Measurement",
    "MetricValue",
    "MetricSpec",
    "PairwiseResult",
    "PlotContext",
    "PluginContractError",
    "ReplicateContext",
    "ScalarMeasurement",
    "ScalarMeasurementAnalysis",
    "SlurmResourceHint",
    "TConditionSummary",
    "TPairwiseResult",
]

for _public_class in (
    BasePlotSettings,
    SlurmResourceHint,
    Condition,
    ReplicateContext,
    AggregateContext,
    ComparisonContext,
    PlotContext,
    MetricValue,
    ConditionSummary,
    PairwiseResult,
    ANOVAResult,
    ComparisonResult,
    BaseConditionSummary,
    BaseComparisonResult,
    MetricSpec,
    CacheIdentity,
    Measurement,
    ScalarMeasurement,
):
    _public_class.__module__ = __name__
del _public_class


class Analysis(ABC):
    """Base class for all PolyzyMD analyses.

    Subclasses represent a complete analysis lifecycle: per-replicate
    computation, aggregation across replicates, cross-condition comparison,
    plotting, and CLI formatting.
    """

    name: ClassVar[str]
    Settings: ClassVar[type]
    PlotSettingsModel: ClassVar[type[BasePlotSettings] | None] = None
    AggregatedResultClass: ClassVar[type | None] = None
    ReplicateResultClass: ClassVar[type | None] = None
    execution_cost_hint: ClassVar[str] = "medium"
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2
    has_compute_stage: ClassVar[bool] = True
    has_aggregate_stage: ClassVar[bool] = True
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = None
    settings_path_fields: ClassVar[tuple[str, ...]] = ()

    def run_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Run analysis for a single replicate.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context with paths, config, and settings.
        replicate : int
            One-indexed replicate number.

        Returns
        -------
        Any
            Per-replicate result, or ``None`` when compute is disabled.
        """
        if not type(self).has_compute_stage:
            return None
        return self._run_replicate_default(ctx, replicate)

    def _run_replicate_default(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Run the default per-replicate dispatch path.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context with paths, config, and settings.
        replicate : int
            One-indexed replicate number.

        Returns
        -------
        Any
            Per-replicate result, or ``None`` when compute is disabled.
        """
        return run_replicate_default(self, ctx, replicate, base_cls=Analysis)

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate results across replicates for one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate results.

        Returns
        -------
        Any
            Aggregated result, or ``None`` when aggregation is disabled.
        """
        if not type(self).has_aggregate_stage:
            return None
        raise NotImplementedError(
            f"{type(self).__name__} must implement aggregate() or set has_aggregate_stage = False."
        )

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any | None:
        """Build a trajectory-native runner for one replicate.

        The base implementation returns ``None`` to indicate that the plugin is
        not runner-backed.
        """
        del ctx, replicate, universe, window
        return None

    def build_mda_jobs(
        self,
        ctx: MDAReplicateJobContext,
    ) -> Sequence[MDAAnalysisJob] | None:
        """Build MDAnalysis-compatible jobs for one replicate.

        Parameters
        ----------
        ctx : MDAReplicateJobContext
            Framework-provided MDAnalysis job context with a loaded universe,
            frame selection, universe policy, and artifact store.

        Returns
        -------
        sequence of MDAAnalysisJob or None
            Jobs to execute for the replicate, or ``None`` to decline the MDA
            path and allow legacy fallback.
        """

        del ctx
        return None

    def build_mda_collector(self, ctx: MDACollectorContext) -> MDAArtifactCollector:
        """Build the artifact collector for completed MDAnalysis jobs.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context for one replicate.

        Returns
        -------
        MDAArtifactCollector
            Collector that maps completed job results to a replicate artifact.
        """

        del ctx
        from polyzymd.analyses.mda import StrictJSONMDAResultCollector

        return StrictJSONMDAResultCollector()

    def _mda_universe_provider_factory(self) -> type[Any]:
        """Return the universe-provider class used by the MDA job seam.

        Returns
        -------
        type[Any]
            Universe-provider class for the current simulation configuration.
        """

        from polyzymd.analyses.mda import UniverseProvider

        return UniverseProvider

    def _mda_artifact_store_factory(self) -> type[Any]:
        """Return the artifact-store class used by the MDA job seam.

        Returns
        -------
        type[Any]
            Artifact-store class rooted at the replicate output directory.
        """

        from polyzymd.analyses.mda import ArtifactStore

        return ArtifactStore

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return the trajectory-loader class used by the runner seam.

        Returns
        -------
        type[Any]
            Loader class for the current simulation configuration.
        """
        from polyzymd.analyses.shared.loader import TrajectoryLoader

        return TrajectoryLoader

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the frame window for a runner-based replicate analysis.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        loader : Any
            Trajectory loader used for the replicate.
        universe : Any
            Loaded trajectory universe.

        Returns
        -------
        Any
            Resolved trajectory window object.
        """
        from polyzymd.analyses.shared.window import resolve_replicate_trajectory_window

        return resolve_replicate_trajectory_window(
            loader=loader,
            replicate=replicate,
            equilibration=ctx.equilibration,
            n_frames_total=len(universe.trajectory),
        )

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Convert runner output into a PolyzyMD replicate result."""
        del ctx, replicate, runner, window
        raise PluginContractError(
            f"{type(self).__name__} must implement summarize_replicate() when using build_runner()."
        )

    def filter_conditions(
        self,
        conditions: list[Condition],
        settings: BaseModel | None = None,
    ) -> list[Condition]:
        """Filter conditions before comparison.

        Parameters
        ----------
        conditions : list[Condition]
            All conditions from the comparison config.
        settings : BaseModel or None
            Resolved plugin settings.

        Returns
        -------
        list[Condition]
            Conditions to include in analysis.
        """
        del settings
        return list(conditions)

    def compare(self, ctx: ComparisonContext) -> BaseModel | None:
        """Compare results across conditions.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided comparison context.

        Returns
        -------
        BaseModel | None
            Comparison result, or ``None`` if comparison is not supported.
        """
        return _default_compare(self, ctx)

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract scalar metrics from an aggregated result for comparison.

        Parameters
        ----------
        summary : Any
            Aggregated result.

        Returns
        -------
        dict[str, MetricValue]
            Mapping from metric name to metric value.
        """
        del summary
        return {}

    def aggregate_settings_fingerprint(self, settings: BaseModel | None) -> str | None:
        """Return the settings fingerprint expected on aggregate results.

        Parameters
        ----------
        settings : BaseModel or None
            Active analysis settings.

        Returns
        -------
        str or None
            Fingerprint used to validate aggregate caches, or ``None`` to skip
            settings identity checks.
        """
        return _aggregate_settings_fingerprint_impl(settings)

    def validate_aggregated_result(
        self,
        result: Any,
        *,
        condition: Condition | None,
        settings: BaseModel | None,
        equilibration: str,
        source: str | Path | None = None,
        expected_replicates: Sequence[int] | None = None,
        allow_replicate_subset: bool = False,
    ) -> Any:
        """Validate an aggregate result against the active framework context.

        Parameters
        ----------
        result : Any
            Loaded or newly computed aggregate result.
        condition : Condition or None
            Condition providing configuration context.
        settings : BaseModel or None
            Active analysis settings.
        equilibration : str
            Requested equilibration window.
        source : str or Path or None, optional
            Cache path or description used in diagnostics.
        expected_replicates : sequence of int or None, optional
            Replicate IDs expected in the aggregate.
        allow_replicate_subset : bool, optional
            Whether a successful subset of requested replicates is acceptable.

        Returns
        -------
        Any
            Validated aggregate result, potentially coerced through the plugin's
            ``AggregatedResultClass``.
        """
        return _validate_aggregated_result_impl(
            self,
            result,
            condition=condition,
            settings=settings,
            equilibration=equilibration,
            source=source,
            expected_replicates=expected_replicates,
            allow_replicate_subset=allow_replicate_subset,
        )

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate comparison figures.

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided plot context.

        Returns
        -------
        list[Path]
            Paths to generated figure files.
        """
        del ctx
        return []

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format a comparison result for CLI display.

        Parameters
        ----------
        result : Any
            Comparison result to format.
        output_format : str, optional
            Output format.

        Returns
        -------
        str
            Formatted string ready for CLI display.
        """
        if output_format == "json":
            if hasattr(result, "model_dump_json"):
                return result.model_dump_json(indent=2)
            return json.dumps(result, indent=2, default=str)
        return str(result)

    def _load_aggregated_result(self, aggregated_dir: Path) -> Any:
        """Load the aggregated result from disk.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.

        Returns
        -------
        Any
            Loaded result, or ``None`` if no file exists.
        """
        return load_aggregated_result(self, aggregated_dir)

    def _run_replicate_via_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Execute the opt-in runner-based replicate path."""
        return run_replicate_via_runner(self, ctx, replicate, base_cls=Analysis)

    def _run_replicate_via_mda_jobs(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Execute the opt-in MDAnalysis job-based replicate path."""
        return run_replicate_via_mda_jobs(self, ctx, replicate, base_cls=Analysis)

    def _deserialize_result(self, path: Path) -> Any:
        """Load an aggregated result from a JSON file."""
        return _deserialize_result_impl(self, path)

    def _deserialize_replicate_result(self, path: Path) -> Any:
        """Load a single replicate result from disk."""
        return deserialize_replicate_result(self, path)

    def _load_replicate_result(self, run_dir: Path) -> Any | None:
        """Load a replicate result from a run directory."""
        return load_replicate_result(self, run_dir)

    def _check_cache(
        self,
        result_cls: type,
        cache_path: Path,
        *,
        recompute: bool,
        sim_config: SimulationConfig | None = None,
        settings: BaseModel | None = None,
    ) -> Any | None:
        """Load a cached result from disk if valid."""
        return _check_cache_impl(
            self,
            result_cls,
            cache_path,
            recompute=recompute,
            sim_config=sim_config,
            settings=settings,
        )

    @staticmethod
    def replicate_result_path(output_dir: Path) -> Path:
        """Return the canonical per-replicate cache path."""
        return _replicate_result_path(output_dir)

    @staticmethod
    def aggregate_result_path(output_dir: Path) -> Path:
        """Return the canonical aggregated cache path."""
        return _aggregate_result_path(output_dir)

    @staticmethod
    def _format_replicate_range(replicates: Sequence[int]) -> str:
        """Format replicate numbers as a compact string."""
        return format_replicate_range(replicates)

    @staticmethod
    def _build_plot_data(
        ctx: PlotContext,
        *,
        include_replicates: bool = False,
    ) -> tuple[dict[str, Any], list[str]]:
        """Build the data and labels consumed by plotter functions."""
        return _build_plot_data_impl(ctx, include_replicates=include_replicates)

    def comparison_result_path(self, results_dir: Path) -> Path:
        """Return the canonical comparison cache path."""
        return _comparison_result_path(results_dir)

    def figures_output_dir(self, figures_root: Path) -> Path:
        """Return the analysis-specific figure directory."""
        return _figures_output_dir(self, figures_root)

    def save_result(self, result: Any, path: Path) -> Path:
        """Save a result object to disk using a common contract."""
        return _save_result_impl(result, path)

    def resolve_output_dir(
        self,
        analysis_root: Path,
        condition_label: str,
    ) -> Path:
        """Build the analysis output directory for a condition."""
        return _resolve_output_dir(self, analysis_root, condition_label)

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """Validate that subclasses satisfy the analysis contract."""
        super().__init_subclass__(**kwargs)
        validate_analysis_subclass(cls, base_cls=Analysis, kwargs=kwargs)

    def __repr__(self) -> str:
        """Return a concise representation for debugging."""
        return f"<{type(self).__name__}(name={self.name!r})>"


Analysis.__module__ = __name__
from polyzymd.analyses._scalar_measurement_analysis import ScalarMeasurementAnalysis  # noqa: E402

ScalarMeasurementAnalysis.__module__ = __name__
