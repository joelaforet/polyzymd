"""RMSF analysis plugin.

Computes per-residue Root Mean Square Fluctuation from MD trajectories,
aggregates across replicates, compares conditions via the default scalar
comparison pipeline, and produces comparison/profile plots.

All heavy computation is self-contained within this plugin package.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.rmsf._plot_settings import RMSFPlotSettings
from polyzymd.analyses.rmsf._plotters import _plot_rmsf_comparison, _plot_rmsf_profile
from polyzymd.analyses.rmsf._results import RMSFAggregatedResult, RMSFResult
from polyzymd.analyses.rmsf._runner import (
    RMSFReplicateRunner,
    aggregate_per_residue,
    compute_rmsd_timeseries,
    compute_rmsf,
)
from polyzymd.analyses.shared.alignment import align_trajectory
from polyzymd.analyses.shared.centroid import ReferenceMode
from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
from polyzymd.analyses.shared.diagnostics import get_selection_diagnostics
from polyzymd.analyses.shared.loader import TrajectoryLoader, parse_time_string
from polyzymd.analyses.shared.statistics import aggregate_per_residue_stats, compute_sem

if TYPE_CHECKING:
    from numpy.typing import NDArray

logger = logging.getLogger("polyzymd.analyses.rmsf")


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class RMSFSettings(BaseModel):
    """Settings for RMSF analysis.

    Attributes
    ----------
    selection : str
        MDAnalysis selection string for RMSF calculation.
    reference_mode : str
        Reference structure mode: centroid, average, frame, or external.
    reference_frame : int | None
        Frame number if reference_mode is 'frame' (1-indexed).
    reference_file : str | None
        Path to external PDB file if reference_mode is 'external'.
    alignment_selection : str
        MDAnalysis selection for trajectory alignment.
    centroid_selection : str
        MDAnalysis selection for centroid finding (K-Means clustering).
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
        description="MDAnalysis selection for centroid finding (K-Means clustering)",
    )

    @field_validator("reference_mode", mode="after")
    @classmethod
    def validate_reference_mode(cls, v: str) -> str:
        valid = {"centroid", "average", "frame", "external"}
        if v not in valid:
            raise ValueError(f"reference_mode must be one of {valid}, got {v!r}")
        return v


@dataclass(frozen=True)
class _RMSFTrajectoryWindow:
    """RMSF trajectory window that carries loader-derived file metadata."""

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
    ) -> _RMSFTrajectoryWindow:
        """Build an RMSF window wrapper from the shared trajectory window.

        Parameters
        ----------
        window : Any
            Shared trajectory window returned by the centralized resolver.
        trajectory_files : Sequence[Path]
            Trajectory files resolved by the existing loader instance.

        Returns
        -------
        _RMSFTrajectoryWindow
            RMSF window wrapper that preserves run arguments and file metadata.
        """

        return cls(
            start=window.start,
            stop=window.stop,
            step=window.step,
            equilibration_start=window.equilibration_start,
            n_frames_total=window.n_frames_total,
            n_frames_selected=window.n_frames_selected,
            timestep_ps=window.timestep_ps,
            equilibration_ps=window.equilibration_ps,
            warning_message=window.warning_message,
            trajectory_files=tuple(trajectory_files),
        )

    def run_kwargs(self) -> dict[str, int]:
        """Return keyword arguments for the runner ``run()`` call.

        Returns
        -------
        dict[str, int]
            ``start``, ``stop``, and ``step`` values for ``run()``.
        """

        return {"start": self.start, "stop": self.stop, "step": self.step}


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class RMSFAnalysis(Analysis):
    """RMSF analysis: per-residue flexibility from MD trajectories.

    Performs RMSF through the runner seam while preserving the historical
    cache layout:

    1. Load trajectories from config
    2. Apply equilibration offset
    3. Find reference frame and align trajectory
    4. Compute autocorrelation and select independent frames
    5. Calculate per-residue RMSF
    6. Aggregate across replicates with SEM

    The ``compare()`` method adds RMSF-specific completeness guards before
    delegating to the default scalar comparison path so missing or incomplete
    condition results fail loudly.

    Plots
    -----
    - ``rmsf_comparison.png``: Horizontal bar chart of mean RMSF per condition.
    - ``rmsf_profile.png``: Per-residue RMSF overlay with optional SS annotation.
    """

    name: ClassVar[str] = "rmsf"
    Settings: ClassVar[type] = RMSFSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = RMSFPlotSettings
    AggregatedResultClass: ClassVar[type] = RMSFAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = RMSFResult
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel | Any) -> str:
        """Build a short cache tag from RMSF settings.

        Parameters
        ----------
        settings : BaseModel or Any
            RMSF settings model or a legacy settings-like object.

        Returns
        -------
        str
            First 8 hex characters from shared settings fingerprinting.
        """
        if isinstance(settings, BaseModel):
            return settings_fingerprint(settings)

        normalized = RMSFSettings(
            selection=settings.selection,
            reference_mode=settings.reference_mode,
            reference_frame=settings.reference_frame,
            reference_file=settings.reference_file,
            alignment_selection=settings.alignment_selection,
            centroid_selection=settings.centroid_selection,
        )
        return settings_fingerprint(normalized)

    @staticmethod
    def _validate_aggregate_input_completeness(
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> None:
        """Validate that aggregation inputs cover the configured replicates.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate RMSF results.

        Raises
        ------
        ValueError
            Raised when configured replicates are missing, duplicated, or
            unexpected in the aggregation inputs.
        """
        if not results:
            raise ValueError(
                f"RMSF aggregation for condition '{ctx.condition.label}' requires at least one "
                "replicate result. No replicate inputs were provided."
            )

        expected_replicates = sorted(ctx.replicates)
        observed_replicates = [
            result.replicate for result in results if getattr(result, "replicate", None) is not None
        ]
        observed_counter = Counter(observed_replicates)
        duplicate_replicates = sorted(
            replicate for replicate, count in observed_counter.items() if count > 1
        )
        missing_replicates = sorted(set(expected_replicates) - set(observed_replicates))
        unexpected_replicates = sorted(set(observed_replicates) - set(expected_replicates))
        missing_metadata_count = len(results) - len(observed_replicates)
        if (
            missing_replicates
            or unexpected_replicates
            or duplicate_replicates
            or missing_metadata_count
        ):
            details: list[str] = []
            if missing_replicates:
                details.append(f"missing replicates {missing_replicates}")
            if unexpected_replicates:
                details.append(f"unexpected replicates {unexpected_replicates}")
            if duplicate_replicates:
                details.append(f"duplicate replicates {duplicate_replicates}")
            if missing_metadata_count:
                details.append(f"results without replicate identifiers {missing_metadata_count}")
            detail_text = "; ".join(details)
            raise ValueError(
                f"RMSF aggregation for condition '{ctx.condition.label}' is incomplete. "
                f"Expected replicate results for {expected_replicates}; observed "
                f"{sorted(observed_replicates)} ({detail_text}). Recompute missing replicates "
                "or clear stale caches before aggregating."
            )

    @staticmethod
    def _order_aggregate_results_by_replicate(
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> list[Any]:
        """Return aggregate inputs in declared replicate order.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate RMSF results that already passed completeness checks.

        Returns
        -------
        list[Any]
            Replicate results ordered to match ``ctx.replicates`` exactly.
        """
        replicate_to_result = {result.replicate: result for result in results}
        return [replicate_to_result[replicate] for replicate in ctx.replicates]

    @classmethod
    def _validate_replicate_result_settings_identity(
        cls,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> None:
        """Validate settings fingerprints on per-replicate RMSF results.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate RMSF results.

        Raises
        ------
        ValueError
            Raised when replicate results are missing settings fingerprints or
            were computed with different settings.
        """
        expected_fingerprint = cls._make_settings_cache_tag(ctx.settings)
        missing_fingerprint_replicates: list[int] = []
        mismatched_fingerprints: list[str] = []

        for result in results:
            replicate = getattr(result, "replicate", None)
            stored_fingerprint = getattr(result, "settings_fingerprint", None)
            if stored_fingerprint is None:
                stored_fingerprint = getattr(result, "settings_fp", None)

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
                f"RMSF aggregation for condition '{ctx.condition.label}' cannot use legacy cached "
                "replicate results missing settings fingerprints. Affected replicates: "
                f"{sorted(missing_fingerprint_replicates)}. Recompute the condition to refresh "
                "settings-sensitive caches before aggregating."
            )

        if mismatched_fingerprints:
            mismatch_text = "; ".join(mismatched_fingerprints)
            raise ValueError(
                f"RMSF aggregation for condition '{ctx.condition.label}' detected settings "
                f"fingerprint mismatches ({mismatch_text}). Recompute the condition or clear "
                "stale caches before aggregating."
            )

    @classmethod
    def _validate_aggregated_result_completeness(
        cls,
        condition: Any,
        agg_result: Any,
        settings: BaseModel | Any,
    ) -> None:
        """Validate that an aggregated RMSF result is complete for comparison.

        Parameters
        ----------
        condition : Any
            Condition associated with the aggregated result.
        agg_result : Any
            Aggregated RMSF result to validate.
        settings : BaseModel or Any
            Current RMSF settings used for comparison.

        Raises
        ------
        ValueError
            Raised when the aggregated result is missing settings identity or
            replicate coverage required for comparison.
        """
        expected_replicates = sorted(condition.replicates)
        observed_replicates = sorted(getattr(agg_result, "replicates", []))
        n_replicates = getattr(agg_result, "n_replicates", None)
        if n_replicates != len(expected_replicates) or observed_replicates != expected_replicates:
            raise ValueError(
                f"Aggregated RMSF result for condition '{condition.label}' has incomplete "
                f"replicate coverage. Expected replicates {expected_replicates}, found "
                f"{observed_replicates} with n_replicates={n_replicates}. Recompute the "
                "condition or clear stale caches before comparing."
            )

        per_replicate_values = list(getattr(agg_result, "per_replicate_mean_rmsf", []))
        if len(per_replicate_values) != len(expected_replicates):
            raise ValueError(
                f"Aggregated RMSF result for condition '{condition.label}' has incomplete "
                f"replicate values: per_replicate_mean_rmsf has {len(per_replicate_values)} "
                f"entries, expected {len(expected_replicates)}. Recompute the condition or "
                "clear stale caches before comparing."
            )

        stored_fingerprint = getattr(agg_result, "settings_fingerprint", None)
        if stored_fingerprint is None:
            stored_fingerprint = getattr(agg_result, "settings_fp", None)
        current_fingerprint = cls._make_settings_cache_tag(settings)
        if stored_fingerprint is None:
            raise ValueError(
                f"Aggregated RMSF result for condition '{condition.label}' is missing a settings "
                "fingerprint. Legacy RMSF aggregated caches are not compatible with "
                "settings-sensitive comparison. Recompute the condition before comparing."
            )
        if stored_fingerprint != current_fingerprint:
            raise ValueError(
                f"Aggregated RMSF result for condition '{condition.label}' was computed with "
                f"settings fingerprint {stored_fingerprint}, but current settings require "
                f"{current_fingerprint}. Recompute the condition or clear stale caches before "
                "comparing."
            )

    def run_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Run RMSF for a single replicate.

        Performs settings validation and delegates trajectory-native execution
        through the runner seam.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        RMSFResult
            Per-replicate RMSF result.
        """
        settings = ctx.settings

        selection = settings.selection
        reference_mode: ReferenceMode = settings.reference_mode
        reference_frame = settings.reference_frame
        reference_file = settings.reference_file
        alignment_selection = settings.alignment_selection
        centroid_selection = settings.centroid_selection

        # Validate reference_mode + reference_frame/file combinations
        if reference_mode == "frame" and reference_frame is None:
            raise ValueError("reference_frame is required when reference_mode='frame'")
        if reference_mode == "external" and reference_file is None:
            raise ValueError(
                "reference_file is required when reference_mode='external'. "
                "Provide a path to the external PDB reference structure."
            )
        if reference_mode == "external" and reference_file is not None:
            ref_path = Path(reference_file)
            if not ref_path.exists():
                raise ValueError(
                    f"reference_file does not exist: {ref_path}. "
                    "Provide a valid path to the external PDB reference structure."
                )

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        output_dir = ctx.output_dir
        eq_str = f"eq{eq_value:g}{eq_unit}"
        settings_tag = self._make_settings_cache_tag(ctx.settings)
        result_filename = f"rmsf_{eq_str}_{settings_tag}.json"
        result_file = output_dir / result_filename
        legacy_result_file = output_dir / f"rmsf_{eq_str}.json"

        if not ctx.recompute and legacy_result_file.exists() and not result_file.exists():
            logger.info(
                "Ignoring legacy RMSF cache without settings fingerprint at %s and recomputing "
                "with settings-sensitive cache identity",
                legacy_result_file,
            )

        cached = self._check_cache(
            RMSFResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=ctx.sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

        result = super().run_replicate(ctx, replicate)
        result.save(result_file)
        logger.info(f"Saved result to {result_file}")
        return result

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return the RMSF loader class for the shared runner seam.

        Returns
        -------
        type[Any]
            Loader class patched by RMSF unit tests.
        """

        return TrajectoryLoader

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any:
        """Build the runner-backed RMSF execution object.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        universe : Any
            Loaded universe for the replicate.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        Any
            Runner object compatible with the trajectory seam.
        """

        del replicate
        return RMSFReplicateRunner(
            universe=universe,
            settings=ctx.settings,
            timestep_ps=window.timestep_ps,
            align_trajectory_func=align_trajectory,
            get_selection_diagnostics_func=get_selection_diagnostics,
            compute_rmsd_timeseries_func=_compute_rmsd_timeseries,
            compute_rmsf_func=_compute_rmsf,
            aggregate_per_residue_func=_aggregate_per_residue,
        )

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the RMSF window and retain trajectory file metadata.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        loader : Any
            Trajectory loader already constructed for this replicate.
        universe : Any
            Loaded universe for the replicate.

        Returns
        -------
        Any
            Shared trajectory window augmented with trajectory file metadata.
        """

        window = super().get_trajectory_window(ctx, replicate, loader, universe)
        traj_info = loader.get_trajectory_info(replicate)
        return _RMSFTrajectoryWindow.from_window(window, traj_info.trajectory_files)

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Serialize runner output into the legacy RMSF result schema.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        runner : Any
            Executed RMSF runner.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        RMSFResult
            Cache-compatible per-replicate RMSF result.
        """

        from polyzymd.analyses._results_base import get_polyzymd_version

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        settings_tag = self._make_settings_cache_tag(ctx.settings)
        trajectory_files = getattr(window, "trajectory_files", None)
        if trajectory_files is None:
            trajectory_files = getattr(runner, "trajectory_files", ())
        payload = runner.results

        return RMSFResult(
            config_hash=compute_config_hash(ctx.sim_config),
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string=payload.selection,
            correlation_time=payload.correlation_time,
            correlation_time_unit=payload.correlation_time_unit,
            n_independent_frames=payload.n_independent_frames,
            residue_ids=payload.residue_ids,
            residue_names=payload.residue_names,
            rmsf_values=payload.rmsf_values.tolist(),
            mean_rmsf=payload.mean_rmsf,
            std_rmsf=payload.std_rmsf,
            min_rmsf=payload.min_rmsf,
            max_rmsf=payload.max_rmsf,
            reference_mode=payload.reference_mode,
            reference_frame=payload.reference_frame,
            alignment_selection=payload.alignment_selection,
            reference_file=payload.reference_file,
            n_frames_total=payload.n_frames_total,
            n_frames_used=payload.n_frames_used,
            settings_fingerprint=settings_tag,
            trajectory_files=[str(path) for path in trajectory_files],
        )

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate RMSF across replicates for one condition.

        Computes per-residue mean +/- SEM and overall statistics from
        the already-computed per-replicate results.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[RMSFResult]
            Per-replicate RMSF results.

        Returns
        -------
        RMSFAggregatedResult
            Aggregated result with per-residue and overall statistics.
        """
        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rmsf._results import RMSFAggregatedResult

        settings = ctx.settings
        self._validate_aggregate_input_completeness(ctx, results)
        self._validate_replicate_result_settings_identity(ctx, results)
        ordered_results = self._order_aggregate_results_by_replicate(ctx, results)
        settings_tag = self._make_settings_cache_tag(settings)

        # Collect per-residue RMSF arrays from each replicate
        per_replicate_rmsf = [np.array(result.rmsf_values) for result in ordered_results]

        # Aggregate per-residue statistics
        per_residue_stats = aggregate_per_residue_stats(
            per_replicate_rmsf,
            residue_ids=np.array(ordered_results[0].residue_ids),
        )

        # Aggregate whole-protein statistics
        per_replicate_means = [result.mean_rmsf for result in ordered_results]
        overall_stats = compute_sem(per_replicate_means)

        config_hash = ordered_results[0].config_hash

        agg_result = RMSFAggregatedResult(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=ordered_results[0].equilibration_time,
            equilibration_unit=ordered_results[0].equilibration_unit,
            selection_string=settings.selection,
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            residue_ids=ordered_results[0].residue_ids,
            residue_names=ordered_results[0].residue_names,
            mean_rmsf_per_residue=per_residue_stats.means.tolist(),
            sem_rmsf_per_residue=per_residue_stats.sems.tolist(),
            per_replicate_mean_rmsf=per_replicate_means,
            overall_mean_rmsf=overall_stats.mean,
            overall_sem_rmsf=overall_stats.sem,
            overall_min_rmsf=float(np.min(per_residue_stats.means)),
            overall_max_rmsf=float(np.max(per_residue_stats.means)),
            settings_fingerprint=settings_tag,
            source_result_files=[],
        )

        target_path = ctx.result_path
        if target_path is None:
            filename = self._make_aggregated_filename(ctx.replicates, ordered_results[0])
            target_path = ctx.output_dir / filename
        self.save_result(agg_result, target_path)
        logger.info(f"Saved aggregated RMSF to {target_path}")

        return agg_result

    # === Optional methods ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract the mean RMSF metric for scalar comparison.

        Parameters
        ----------
        summary : RMSFAggregatedResult
            Aggregated RMSF result.

        Returns
        -------
        dict[str, MetricValue]
            Single entry ``"mean_rmsf"`` with ``higher_is_better=False``
            (lower RMSF = more stable).
        """
        return {
            "mean_rmsf": MetricValue(
                name="mean_rmsf",
                mean=summary.overall_mean_rmsf,
                sem=summary.overall_sem_rmsf,
                replicate_values=summary.per_replicate_mean_rmsf,
                higher_is_better=False,
                direction_labels=("stabilizing", "unchanged", "destabilizing"),
            ),
        }

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare RMSF across conditions with fail-loud completeness checks.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided comparison context.

        Returns
        -------
        Any
            Default scalar comparison result, or ``None`` when no conditions are
            available.
        """
        from polyzymd.analyses.stats import default_scalar_comparison

        metrics_by_condition: dict[str, dict[str, MetricValue]] = {}
        for condition in ctx.conditions:
            summary = ctx.aggregated_results.get(condition.label)
            source: Path | str = f"comparison context for {condition.label}"
            if summary is None:
                analysis_dir = ctx.analysis_dirs.get(condition.label)
                if analysis_dir is None:
                    raise ValueError(
                        f"RMSF comparison requires an aggregated result for condition "
                        f"'{condition.label}', but no analysis directory was found. Recompute "
                        "the condition or clear stale caches before comparing."
                    )

                agg_dir = analysis_dir / "aggregated"
                source = self.aggregate_result_path(agg_dir)
                summary = self._load_aggregated_result(agg_dir)

            if summary is None:
                raise ValueError(
                    f"RMSF comparison requires an aggregated result for condition "
                    f"'{condition.label}'. Recompute the condition or clear stale caches before "
                    "comparing."
                )

            summary = self.validate_aggregated_result(
                summary,
                condition=condition,
                settings=ctx.settings,
                equilibration=ctx.equilibration,
                source=source,
                expected_replicates=condition.replicates,
                allow_replicate_subset=True,
            )
            self._validate_aggregated_result_completeness(condition, summary, ctx.settings)

            extracted = self.extract_metrics(summary)
            if not isinstance(extracted, dict):
                raise PluginContractError(
                    f"RMSF compare expected extract_metrics() to return dict[str, MetricValue] "
                    f"for condition '{condition.label}', got {type(extracted).__name__}."
                )
            if not extracted:
                raise PluginContractError(
                    f"RMSF compare expected at least one metric for condition '{condition.label}'."
                )
            for metric_key, metric_value in extracted.items():
                if not isinstance(metric_value, MetricValue):
                    raise PluginContractError(
                        f"RMSF compare expected MetricValue for key '{metric_key}' in condition "
                        f"'{condition.label}', got {type(metric_value).__name__}."
                    )
            metrics_by_condition[condition.label] = extracted

        if not metrics_by_condition:
            logger.warning("%s: no conditions have metrics — skipping comparison.", self.name)
            return None

        return default_scalar_comparison(
            analysis_name=self.name,
            project_name=ctx.name,
            metrics_by_condition=metrics_by_condition,
            control_label=ctx.effective_control,
            equilibration=ctx.equilibration,
            fdr_alpha=ctx.fdr_alpha,
            ttest_method=ctx.ttest_method,
            posthoc_method=ctx.posthoc_method,
        )

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format RMSF comparison result for CLI display.

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
                title="RMSF Comparison",
                metric_label="Mean RMSF",
                metric_unit="A",
                metric_key="mean_rmsf",
                output_format=output_format,
                higher_is_better=False,
            )
        return super().format(result, output_format)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate RMSF comparison and profile plots.

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

        data, labels = self._build_plot_data(ctx)
        if not labels:
            return plots

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        plot_settings = ctx.plot_settings

        result = _plot_rmsf_comparison(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        result = _plot_rmsf_profile(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        return plots

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Backward-compatible filename helper retained for tests."""
        eq_str = f"eq{first_result.equilibration_time:g}{first_result.equilibration_unit}"
        rep_str = Analysis._format_replicate_range(replicates)
        return f"rmsf_{rep_str}_{eq_str}.json"


# ---------------------------------------------------------------------------
# Private helper functions (extracted from legacy RMSFCalculator)
# ---------------------------------------------------------------------------


def _compute_rmsd_timeseries(
    u: Any,
    atoms: Any,
    start_frame: int,
    stop_frame: int | None = None,
    step: int = 1,
) -> NDArray[np.float64]:
    """Compute RMSD timeseries for autocorrelation analysis.

    This compatibility wrapper preserves the historical helper name while the
    implementation lives in ``rmsf._runner``.
    """

    return compute_rmsd_timeseries(
        u,
        atoms,
        start_frame=start_frame,
        stop_frame=stop_frame,
        step=step,
    )


def _compute_rmsf(
    u: Any,
    atoms: Any,
    frame_indices: NDArray[np.int64],
    reference_positions: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Compute RMSF using selected frames.

    This compatibility wrapper preserves the historical helper name while the
    implementation lives in ``rmsf._runner``.
    """

    return compute_rmsf(u, atoms, frame_indices, reference_positions)


def _aggregate_per_residue(
    atoms: Any,
    atom_rmsf: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Aggregate per-atom RMSF to per-residue (mean within residue).

    This compatibility wrapper preserves the historical helper name while the
    implementation lives in ``rmsf._runner``.
    """

    return aggregate_per_residue(atoms, atom_rmsf)
