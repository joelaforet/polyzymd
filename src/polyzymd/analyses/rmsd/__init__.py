"""RMSD analysis plugin.

Computes per-frame Root Mean Square Deviation for one or more named selections,
aggregates across replicates, performs per-run cross-condition comparisons,
and exposes plotting/formatting hooks.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.rmsd._plot_settings import RMSDPlotSettings
from polyzymd.analyses.rmsd._results import RMSDAggregatedResult, RMSDResult
from polyzymd.analyses.rmsd._runner import RMSDReplicateRunner, compute_rmsd_run
from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
from polyzymd.analyses.shared.loader import TrajectoryLoader, parse_time_string
from polyzymd.analyses.shared.multi_run_comparison import (
    apply_fdr_correction,
    build_condition_pairs,
)
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.rmsd._comparison_results import RMSDComparisonResult

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _RMSDTrajectoryWindow:
    """RMSD trajectory window that carries loader-derived file metadata.

    This keeps summarize-time metadata lookup on the same loader seam used by
    the base runner orchestration.
    """

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
    ) -> _RMSDTrajectoryWindow:
        """Build an RMSD window wrapper from the shared trajectory window.

        Parameters
        ----------
        window : Any
            Shared trajectory window returned by the centralized resolver.
        trajectory_files : Sequence[Path]
            Trajectory files resolved by the existing loader instance.

        Returns
        -------
        _RMSDTrajectoryWindow
            RMSD window wrapper that preserves run arguments and file metadata.
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


class RMSDRunSettings(BaseModel):
    """Settings for a single RMSD run.

    Attributes
    ----------
    label : str
        Human-readable run label.
    selection : str
        MDAnalysis selection used for RMSD calculation.
    alignment_selection : str
        MDAnalysis selection used for trajectory alignment.
    reference_mode : str
        Reference mode: ``"centroid"``, ``"average"``, ``"frame"``, or ``"external"``.
    reference_frame : int
        0-indexed frame index used when ``reference_mode="frame"``.
    reference_file : str | None
        Path to external PDB file when ``reference_mode="external"``.
    centroid_selection : str | None
        Selection used for centroid finding when ``reference_mode="centroid"``.
        If ``None``, ``alignment_selection`` is used.
    """

    label: str = Field(..., description="Human-readable run label")
    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for RMSD calculation",
    )
    alignment_selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for trajectory alignment",
    )
    reference_mode: str = Field(
        default="centroid",
        description="Reference mode: centroid, average, frame, or external",
    )
    reference_frame: int = Field(
        default=0,
        description="0-indexed reference frame for reference_mode='frame'",
    )
    reference_file: str | None = Field(
        default=None,
        description="Path to external PDB file for reference_mode='external'",
    )
    centroid_selection: str | None = Field(
        default=None,
        description="Selection for centroid mode; defaults to alignment_selection",
    )
    convergence_window_size_ns: float = Field(
        default=15.0,
        gt=0,
        description="Sliding window size for convergence detection in ns",
    )
    convergence_step_size_ns: float = Field(
        default=5.0,
        gt=0,
        description="Sliding window step size for convergence detection in ns",
    )
    convergence_slope_threshold: float = Field(
        default=0.0005,
        gt=0,
        description="Absolute slope threshold for convergence detection",
    )
    convergence_sustained_for_ns: float = Field(
        default=15.0,
        gt=0,
        description="Required sustained converged duration in ns",
    )

    @field_validator("reference_mode", mode="after")
    @classmethod
    def validate_reference_mode(cls, value: str) -> str:
        """Validate reference mode value."""
        valid_modes = {"centroid", "average", "frame", "external"}
        if value not in valid_modes:
            raise ValueError(f"reference_mode must be one of {valid_modes}, got {value!r}")
        return value

    @model_validator(mode="after")
    def validate_external_reference(self) -> "RMSDRunSettings":
        """Validate external reference requirements."""
        if self.reference_mode != "external":
            return self

        if self.reference_file is None:
            raise ValueError(
                "reference_file is required when reference_mode='external'. "
                "Provide a path to the external PDB reference structure."
            )

        ref_path = Path(self.reference_file)
        if not ref_path.exists():
            raise ValueError(
                f"reference_file does not exist: {ref_path}. "
                "Provide a valid path to the external PDB reference structure."
            )

        return self

    @model_validator(mode="after")
    def validate_convergence_window_step(self) -> "RMSDRunSettings":
        """Validate convergence window and step compatibility."""
        if self.convergence_step_size_ns > self.convergence_window_size_ns:
            raise ValueError(
                "convergence_step_size_ns must be less than or equal to convergence_window_size_ns"
            )
        return self


class RMSDSettings(BaseModel):
    """Top-level RMSD settings.

    Attributes
    ----------
    runs : list[RMSDRunSettings]
        Named RMSD runs to compute.
    """

    runs: list[RMSDRunSettings] = Field(
        default_factory=list,
        description="RMSD runs to compute",
    )

    @field_validator("runs", mode="after")
    @classmethod
    def validate_runs_nonempty(cls, value: list[RMSDRunSettings]) -> list[RMSDRunSettings]:
        """Ensure at least one run exists."""
        if not value:
            raise ValueError("At least one RMSD run must be defined")
        return value

    @field_validator("runs", mode="after")
    @classmethod
    def validate_unique_labels(cls, value: list[RMSDRunSettings]) -> list[RMSDRunSettings]:
        """Ensure run labels are unique."""
        labels = [run.label for run in value]
        if len(labels) != len(set(labels)):
            raise ValueError("RMSD run labels must be unique")
        return value


class RMSDAnalysis(Analysis):
    """RMSD analysis plugin using a multi-run comparison model."""

    name: ClassVar[str] = "rmsd"
    min_replicates: ClassVar[int] = 1
    Settings: ClassVar[type] = RMSDSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = RMSDPlotSettings
    AggregatedResultClass: ClassVar[type] = RMSDAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = RMSDResult
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel) -> str:
        """Build a short cache tag for settings and equilibration.

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

    @classmethod
    def _coerce_and_validate_aggregated_result(
        cls,
        result: Any,
        settings: RMSDSettings,
        *,
        condition_label: str | None = None,
        source: Path | None = None,
    ) -> RMSDAggregatedResult:
        """Coerce an aggregated result and validate its settings identity.

        Parameters
        ----------
        result : Any
            Aggregated result loaded from disk or supplied in memory.
        settings : RMSDSettings
            Current RMSD settings for comparison or plotting.
        condition_label : str | None, optional
            Condition label for error reporting.
        source : Path | None, optional
            Source file path for diagnostics.

        Returns
        -------
        RMSDAggregatedResult
            Validated aggregated result.

        Raises
        ------
        ValueError
            Raised when the aggregated result is missing a settings
            fingerprint or was computed with different settings.
        """
        if isinstance(result, dict):
            result = RMSDAggregatedResult.model_validate(result)

        if not isinstance(result, RMSDAggregatedResult):
            raise TypeError(
                f"RMSD aggregated result loader expected RMSDAggregatedResult, got "
                f"{type(result).__name__}"
            )

        stored_fingerprint = getattr(result, "settings_fingerprint", None)
        if stored_fingerprint is None:
            stored_fingerprint = getattr(result, "settings_fp", None)

        current_fingerprint = cls._make_settings_cache_tag(settings)
        condition_text = (
            f" for condition '{condition_label}'" if condition_label is not None else ""
        )
        source_text = f" at {source}" if source is not None else ""
        if stored_fingerprint is None:
            raise ValueError(
                "Aggregated RMSD result"
                f"{condition_text} is missing a settings fingerprint{source_text}. "
                "Legacy RMSD aggregated caches are not compatible with "
                "settings-sensitive compare/plot loading. Recompute the condition before "
                "comparing or plotting."
            )
        if stored_fingerprint != current_fingerprint:
            raise ValueError(
                "Aggregated RMSD result"
                f"{condition_text} was computed with settings fingerprint "
                f"{stored_fingerprint}, but current settings require {current_fingerprint}"
                f"{source_text}. Recompute the condition or clear stale caches before "
                "comparing or plotting."
            )
        return result

    def _resolve_aggregated_result_path(self, aggregated_dir: Path) -> Path | None:
        """Resolve the aggregated RMSD result path.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.

        Returns
        -------
        Path | None
            Path to the selected JSON result, or ``None`` when no result file
            exists.
        """
        if not aggregated_dir.exists():
            return None
        canonical = self.aggregate_result_path(aggregated_dir)
        if canonical.exists():
            return canonical

        json_files = sorted(aggregated_dir.glob("*.json"), key=lambda p: p.stat().st_mtime)
        if not json_files:
            return None

        chosen = json_files[-1]
        logger.warning(
            "%s: canonical result.json not found in %s — falling back to %s "
            "(%d JSON file(s) present)",
            self.name,
            aggregated_dir,
            chosen.name,
            len(json_files),
        )
        return chosen

    def _load_aggregated_result(
        self,
        aggregated_dir: Path,
        *,
        settings: RMSDSettings | None = None,
        condition_label: str | None = None,
    ) -> Any:
        """Load and optionally validate an aggregated RMSD result.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.
        settings : RMSDSettings | None, optional
            Current settings used to validate settings-sensitive aggregated
            caches. When omitted, the result is loaded without settings
            identity validation.
        condition_label : str | None, optional
            Condition label for validation diagnostics.

        Returns
        -------
        Any
            Loaded aggregated result, or ``None`` when no result file exists.
        """
        result_path = self._resolve_aggregated_result_path(aggregated_dir)
        if result_path is None:
            return None

        result = self._deserialize_result(result_path)
        if settings is None:
            return result

        return self._coerce_and_validate_aggregated_result(
            result,
            settings,
            condition_label=condition_label,
            source=result_path,
        )

    @staticmethod
    def _validate_aggregate_input_completeness(
        ctx: AggregateContext,
        results: Sequence[Any],
        configured_run_labels: Sequence[str],
    ) -> None:
        """Validate that aggregation inputs cover all configured replicates and runs.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate RMSD results.
        configured_run_labels : Sequence[str]
            Run labels defined in the RMSD settings.

        Raises
        ------
        ValueError
            Raised when configured replicates or runs are missing from the
            aggregation inputs.
        """
        expected_replicates = sorted(ctx.replicates)
        observed_replicates = sorted(
            result.replicate for result in results if getattr(result, "replicate", None) is not None
        )
        if observed_replicates != expected_replicates:
            raise ValueError(
                f"RMSD aggregation for condition '{ctx.condition.label}' is incomplete. Expected "
                f"replicate results for {expected_replicates}, found {observed_replicates}. "
                "Recompute missing replicates or clear stale caches before aggregating."
            )

        for run_label in configured_run_labels:
            missing_replicates: list[int] = []
            duplicate_replicates: list[int] = []
            for result in results:
                replicate = getattr(result, "replicate", None)
                matches = [
                    run_result
                    for run_result in result.run_results
                    if run_result.run_label == run_label
                ]
                if not matches:
                    if replicate is not None:
                        missing_replicates.append(replicate)
                    continue
                if len(matches) > 1 and replicate is not None:
                    duplicate_replicates.append(replicate)

            if missing_replicates:
                raise ValueError(
                    f"Configured RMSD run '{run_label}' is missing replicate entries in condition "
                    f"'{ctx.condition.label}'. Missing replicates: {sorted(missing_replicates)}. "
                    "Recompute missing replicates or clear stale caches before aggregating."
                )

            if duplicate_replicates:
                raise ValueError(
                    f"Configured RMSD run '{run_label}' has duplicate replicate entries in "
                    f"condition '{ctx.condition.label}' for replicates "
                    f"{sorted(duplicate_replicates)}. Clear stale caches and recompute before "
                    "aggregating."
                )

    @staticmethod
    def _order_results_by_replicate(
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> list[Any]:
        """Return aggregate inputs ordered to match ``ctx.replicates``.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate RMSD results that already passed completeness checks.

        Returns
        -------
        list[Any]
            Replicate results in the declared replicate order.

        Raises
        ------
        ValueError
            Raised when a replicate identifier is missing or duplicated while
            normalizing the aggregate inputs.
        """
        ordered_results: dict[int, Any] = {}
        for result in results:
            replicate = getattr(result, "replicate", None)
            if replicate is None:
                raise ValueError(
                    f"RMSD aggregation for condition '{ctx.condition.label}' encountered a "
                    "replicate result without a replicate identifier while normalizing "
                    "aggregate inputs."
                )
            if replicate in ordered_results:
                raise ValueError(
                    f"RMSD aggregation for condition '{ctx.condition.label}' encountered "
                    f"duplicate replicate {replicate} while normalizing aggregate inputs."
                )
            ordered_results[replicate] = result

        missing_replicates = [
            replicate for replicate in ctx.replicates if replicate not in ordered_results
        ]
        if missing_replicates:
            raise ValueError(
                f"RMSD aggregation for condition '{ctx.condition.label}' cannot order aggregate "
                f"inputs because replicates {missing_replicates} are missing."
            )

        return [ordered_results[replicate] for replicate in ctx.replicates]

    @staticmethod
    def _validate_aggregated_result_completeness(
        condition: Any,
        agg_result: RMSDAggregatedResult,
        configured_run_labels: Sequence[str],
    ) -> None:
        """Validate that an aggregated RMSD result is complete for comparison.

        Parameters
        ----------
        condition : Any
            Condition associated with the aggregated result.
        agg_result : RMSDAggregatedResult
            Aggregated RMSD result to validate.
        configured_run_labels : Sequence[str]
            Run labels defined in the RMSD settings.

        Raises
        ------
        ValueError
            Raised when the aggregated result omits configured runs or contains
            incomplete per-run replicate data.
        """
        expected_run_labels = set(configured_run_labels)
        observed_run_labels = {run_result.run_label for run_result in agg_result.run_results}
        missing_runs = sorted(expected_run_labels - observed_run_labels)
        unexpected_runs = sorted(observed_run_labels - expected_run_labels)
        if missing_runs or unexpected_runs:
            details: list[str] = []
            if missing_runs:
                details.append(f"missing runs {missing_runs}")
            if unexpected_runs:
                details.append(f"unexpected runs {unexpected_runs}")
            detail_text = "; ".join(details)
            raise ValueError(
                f"Aggregated RMSD result for condition '{condition.label}' is incomplete: "
                f"{detail_text}. Recompute the condition or clear stale caches before "
                "comparing."
            )

        expected_replicates = sorted(condition.replicates)
        observed_replicates = sorted(agg_result.replicates)
        if (
            agg_result.n_replicates != len(expected_replicates)
            or observed_replicates != expected_replicates
        ):
            raise ValueError(
                f"Aggregated RMSD result for condition '{condition.label}' has incomplete "
                f"replicate coverage. Expected replicates {expected_replicates}, found "
                f"{observed_replicates} with n_replicates={agg_result.n_replicates}. Recompute "
                "the condition or clear stale caches before comparing."
            )

        for run_result in agg_result.run_results:
            run_replicates = sorted(run_result.replicates)
            counts = {
                "per_replicate_means": len(run_result.per_replicate_means),
                "per_replicate_stds": len(run_result.per_replicate_stds),
                "per_replicate_medians": len(run_result.per_replicate_medians),
                "per_replicate_convergence_times_ns": len(
                    run_result.per_replicate_convergence_times_ns
                ),
                "per_replicate_convergence_assessable": len(
                    run_result.per_replicate_convergence_assessable
                ),
            }
            mismatched_fields = {
                name: count for name, count in counts.items() if count != len(expected_replicates)
            }
            if (
                run_result.n_replicates != len(expected_replicates)
                or run_replicates != expected_replicates
            ):
                raise ValueError(
                    f"Aggregated RMSD run '{run_result.run_label}' for condition "
                    f"'{condition.label}' has incomplete replicate metadata. Expected "
                    f"replicates {expected_replicates}, found {run_replicates} with "
                    f"n_replicates={run_result.n_replicates}. Recompute the condition or "
                    "clear stale caches before comparing."
                )

            if mismatched_fields:
                raise ValueError(
                    f"Aggregated RMSD run '{run_result.run_label}' for condition "
                    f"'{condition.label}' has incomplete replicate values: {mismatched_fields}. "
                    f"Expected {len(expected_replicates)} entries per field. Recompute the "
                    "condition or clear stale caches before comparing."
                )

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Run RMSD for all configured runs for a single replicate.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        RMSDResult
            Per-replicate RMSD result containing all run outputs.
        """

        settings = ctx.settings
        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        eq_str = f"eq{eq_value:g}{eq_unit}"
        settings_tag = self._make_settings_cache_tag(settings)
        result_file = ctx.output_dir / f"rmsd_{eq_str}_{settings_tag}.json"

        cached = self._check_cache(
            RMSDResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=ctx.sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

        result = super().run_replicate(ctx, replicate)
        result.save(result_file)
        logger.info("Saved RMSD result to %s", result_file)
        return result

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return the RMSD loader class for the shared runner seam.

        Returns
        -------
        type[Any]
            Loader class patched by RMSD unit tests.
        """

        return TrajectoryLoader

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the RMSD window and retain trajectory file metadata.

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
        return _RMSDTrajectoryWindow.from_window(window, traj_info.trajectory_files)

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any:
        """Build the runner-backed RMSD execution object.

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

        return RMSDReplicateRunner(
            sim_config=ctx.sim_config,
            replicate=replicate,
            runs=list(ctx.settings.runs),
            loader_factory=self._trajectory_loader_factory(),
            n_frames_total=len(universe.trajectory),
            timestep_ps=window.timestep_ps,
        )

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Serialize runner output into the legacy RMSD result schema.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        runner : Any
            Executed RMSD runner.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        RMSDResult
            Cache-compatible per-replicate RMSD result.
        """
        import numpy as np

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rmsd._results import RMSDResult, RMSDRunResult

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        eq_str = f"eq{eq_value:g}{eq_unit}"
        settings_tag = self._make_settings_cache_tag(ctx.settings)
        config_hash = compute_config_hash(ctx.sim_config)
        trajectory_files = getattr(window, "trajectory_files", ())

        run_results: list[RMSDRunResult] = []
        for payload in runner.results.run_payloads:
            npz_filename = f"rmsd_{payload.run_label}_{eq_str}_{settings_tag}_timeseries.npz"
            npz_path = ctx.output_dir / npz_filename
            np.savez_compressed(
                npz_path,
                rmsd_values=payload.rmsd_values,
                time_ns=payload.time_ns,
                frames=payload.frames,
                convergence_window_start_ns=np.asarray(
                    payload.convergence_result.window_start_times_ns,
                    dtype=np.float64,
                ),
                convergence_window_mean_rmsd=np.asarray(
                    payload.convergence_result.window_mean_values,
                    dtype=np.float64,
                ),
                convergence_slope_time_ns=np.asarray(
                    payload.convergence_result.slope_times_ns,
                    dtype=np.float64,
                ),
                convergence_slope=np.asarray(payload.convergence_result.slopes, dtype=np.float64),
                convergence_converged=np.asarray(
                    payload.convergence_result.converged,
                    dtype=np.bool_,
                ),
                convergence_time_ns=np.asarray(
                    (
                        np.nan
                        if payload.convergence_result.convergence_time_ns is None
                        else payload.convergence_result.convergence_time_ns
                    ),
                    dtype=np.float64,
                ),
                raw_timestep_ps=np.asarray(payload.raw_timestep_ps, dtype=np.float64),
                frame_stride=np.asarray(payload.frame_stride, dtype=np.int64),
                effective_timestep_ps=np.asarray(payload.effective_timestep_ps, dtype=np.float64),
            )
            run_results.append(
                RMSDRunResult(
                    config_hash=config_hash,
                    polyzymd_version=get_polyzymd_version(),
                    replicate=replicate,
                    equilibration_time=eq_value,
                    equilibration_unit=eq_unit,
                    selection_string=payload.selection,
                    correlation_time=payload.correlation_time,
                    n_independent_frames=payload.n_independent_frames,
                    run_label=payload.run_label,
                    selection=payload.selection,
                    alignment_selection=payload.alignment_selection,
                    reference_mode=payload.reference_mode,
                    reference_frame=payload.reference_frame,
                    reference_file=payload.reference_file,
                    mean_rmsd=payload.mean_rmsd,
                    std_rmsd=payload.std_rmsd,
                    median_rmsd=payload.median_rmsd,
                    min_rmsd=payload.min_rmsd,
                    max_rmsd=payload.max_rmsd,
                    final_rmsd=payload.final_rmsd,
                    sem_rmsd=payload.sem_rmsd,
                    correlation_time_unit=payload.correlation_time_unit,
                    statistical_inefficiency=payload.statistical_inefficiency,
                    autocorrelation_warning=payload.autocorrelation_warning,
                    converged=payload.convergence_result.converged,
                    convergence_assessable=payload.convergence_result.assessable,
                    convergence_time_ns=payload.convergence_result.convergence_time_ns,
                    convergence_message=payload.convergence_result.message,
                    n_frames_total=runner.results.n_frames_total,
                    n_frames_used=window.n_frames_selected,
                    npz_path=str(npz_path),
                    time_unit="ns",
                    timestep_ps=payload.effective_timestep_ps,
                    raw_timestep_ps=payload.raw_timestep_ps,
                    frame_stride=payload.frame_stride,
                )
            )

        return RMSDResult(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string="; ".join(run.selection for run in ctx.settings.runs),
            run_results=run_results,
            settings_fingerprint=settings_tag,
            n_frames_total=runner.results.n_frames_total,
            n_frames_used=window.n_frames_selected,
            trajectory_files=[str(path) for path in trajectory_files],
        )

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate RMSD results across replicates for one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[RMSDResult]
            Per-replicate RMSD results.

        Returns
        -------
        RMSDAggregatedResult
            Aggregated RMSD result for all configured runs.
        """
        import numpy as np

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rmsd._results import RMSDAggregatedResult, RMSDRunAggregatedResult

        if not results:
            raise ValueError(
                f"RMSD aggregation for condition '{ctx.condition.label}' requires at least one "
                "replicate result. No replicate inputs were provided."
            )

        run_labels = [run.label for run in ctx.settings.runs]
        self._validate_aggregate_input_completeness(ctx, results, run_labels)
        ordered_results = self._order_results_by_replicate(ctx, results)
        first = ordered_results[0]
        replicate_order = list(ctx.replicates)

        if len(ctx.replicates) == 1:
            logger.warning(
                "Only one replicate available for RMSD aggregation in condition '%s'; "
                "replicate-level SEM is reported as 0.0",
                ctx.condition.label,
            )

        aggregated_runs: list[RMSDRunAggregatedResult] = []
        for run_label in run_labels:
            run_entries = []
            for result in ordered_results:
                matches = [
                    run_result
                    for run_result in result.run_results
                    if run_result.run_label == run_label
                ]
                if len(matches) != 1:
                    raise ValueError(
                        f"Configured RMSD run '{run_label}' has invalid aggregate inputs in "
                        f"condition '{ctx.condition.label}'. Expected one entry per replicate, "
                        f"found {len(matches)} for replicate {result.replicate}."
                    )
                run_entries.append(matches[0])

            per_means = [entry.mean_rmsd for entry in run_entries]
            per_stds = [entry.std_rmsd for entry in run_entries]
            per_medians = [entry.median_rmsd for entry in run_entries]
            per_convergence_times = [entry.convergence_time_ns for entry in run_entries]
            per_assessable = [entry.convergence_assessable for entry in run_entries]

            n_converged = sum(time is not None for time in per_convergence_times)
            n_assessable = sum(per_assessable)
            convergence_fraction = (
                float(n_converged) / float(n_assessable) if n_assessable > 0 else 0.0
            )
            all_converged = n_assessable > 0 and n_converged == n_assessable
            finite_convergence_times = [time for time in per_convergence_times if time is not None]
            mean_convergence_time_ns = (
                float(np.mean(np.asarray(finite_convergence_times, dtype=np.float64)))
                if finite_convergence_times
                else None
            )
            median_convergence_time_ns = (
                float(np.median(np.asarray(finite_convergence_times, dtype=np.float64)))
                if finite_convergence_times
                else None
            )

            mean_stats = compute_sem(per_means)
            overall_median = float(np.median(np.asarray(per_medians, dtype=np.float64)))

            template = run_entries[0]
            aggregated_runs.append(
                RMSDRunAggregatedResult(
                    config_hash=first.config_hash,
                    polyzymd_version=get_polyzymd_version(),
                    replicate=None,
                    equilibration_time=first.equilibration_time,
                    equilibration_unit=first.equilibration_unit,
                    selection_string=template.selection,
                    replicates=replicate_order,
                    n_replicates=len(replicate_order),
                    run_label=run_label,
                    selection=template.selection,
                    alignment_selection=template.alignment_selection,
                    overall_mean=mean_stats.mean,
                    overall_sem=mean_stats.sem,
                    overall_median=overall_median,
                    per_replicate_means=per_means,
                    per_replicate_stds=per_stds,
                    per_replicate_medians=per_medians,
                    per_replicate_convergence_times_ns=per_convergence_times,
                    per_replicate_convergence_assessable=per_assessable,
                    n_converged_replicates=n_converged,
                    n_assessable_replicates=n_assessable,
                    convergence_fraction=convergence_fraction,
                    all_converged=all_converged,
                    mean_convergence_time_ns=mean_convergence_time_ns,
                    median_convergence_time_ns=median_convergence_time_ns,
                )
            )

        agg_result = RMSDAggregatedResult(
            config_hash=first.config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=first.equilibration_time,
            equilibration_unit=first.equilibration_unit,
            selection_string=first.selection_string,
            replicates=replicate_order,
            n_replicates=len(replicate_order),
            run_results=aggregated_runs,
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
            source_result_files=[],
        )

        target_path = ctx.result_path
        if target_path is None:
            settings_tag = self._make_settings_cache_tag(ctx.settings)
            target_path = ctx.output_dir / self._make_aggregated_filename(
                ctx.replicates,
                first,
                settings_tag,
            )
        self.save_result(agg_result, target_path)
        logger.info("Saved aggregated RMSD result to %s", target_path)

        return agg_result

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare RMSD runs across conditions.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided comparison context.

        Returns
        -------
        RMSDComparisonResult | None
            Comparison result, or ``None`` if no conditions are available to compare.
        """
        from polyzymd import __version__
        from polyzymd.analyses.rmsd._comparison_results import (
            RMSDComparisonResult,
            RMSDConditionSummary,
            RMSDRunANOVA,
            RMSDRunPairwiseComparison,
            RMSDRunSummary,
        )
        from polyzymd.analyses.shared.inferential_statistics import one_way_anova

        run_labels = [run.label for run in ctx.settings.runs]

        summaries: list[RMSDConditionSummary] = []
        for condition in ctx.conditions:
            agg_result = ctx.aggregated_results.get(condition.label)
            if agg_result is not None:
                agg_result = self._coerce_and_validate_aggregated_result(
                    agg_result,
                    ctx.settings,
                    condition_label=condition.label,
                    source=self._resolve_aggregated_result_path(
                        ctx.analysis_dirs[condition.label] / "aggregated"
                    ),
                )
            else:
                agg_dir = ctx.analysis_dirs[condition.label] / "aggregated"
                agg_result = self._load_aggregated_result(
                    agg_dir,
                    settings=ctx.settings,
                    condition_label=condition.label,
                )

            if agg_result is None:
                raise ValueError(
                    f"Missing aggregated RMSD result for condition '{condition.label}'. "
                    "Comparison requires aggregated results for every configured "
                    "condition. Recompute the condition or clear stale caches before "
                    "comparing."
                )

            self._validate_aggregated_result_completeness(condition, agg_result, run_labels)

            run_summaries = [
                RMSDRunSummary(
                    label=run_result.run_label,
                    selection=run_result.selection,
                    mean_rmsd=run_result.overall_mean,
                    sem_rmsd=run_result.overall_sem,
                    per_replicate_means=run_result.per_replicate_means,
                    n_converged_replicates=run_result.n_converged_replicates,
                    n_assessable_replicates=run_result.n_assessable_replicates,
                    convergence_fraction=run_result.convergence_fraction,
                    all_converged=run_result.all_converged,
                    mean_convergence_time_ns=run_result.mean_convergence_time_ns,
                    median_convergence_time_ns=run_result.median_convergence_time_ns,
                )
                for run_result in agg_result.run_results
            ]

            summaries.append(
                RMSDConditionSummary(
                    label=condition.label,
                    config_path=str(condition.config_path),
                    n_replicates=agg_result.n_replicates,
                    run_summaries=run_summaries,
                )
            )

        if not summaries:
            logger.warning("RMSD comparison skipped because no conditions have results")
            return None

        effective_control = ctx.effective_control
        summaries_by_label = {summary.label: summary for summary in summaries}

        ranking_by_run: dict[str, list[str]] = {}
        for run_label in run_labels:
            ranked_labels = sorted(
                summaries_by_label,
                key=lambda label: summaries_by_label[label].get_run(run_label).mean_rmsd,
            )
            ranking_by_run[run_label] = ranked_labels

        pairwise_comparisons: list[RMSDRunPairwiseComparison] = []
        if len(summaries) >= 2:
            for run_label in run_labels:
                condition_pairs = build_condition_pairs(
                    list(summaries_by_label.keys()),
                    effective_control,
                    on_control_missing="skip",
                    logger=logger,
                )

                for condition_a, condition_b in condition_pairs:
                    pairwise_comparisons.append(
                        self._compare_run(
                            run_label=run_label,
                            condition_a=condition_a,
                            condition_b=condition_b,
                            run_a=summaries_by_label[condition_a].get_run(run_label),
                            run_b=summaries_by_label[condition_b].get_run(run_label),
                        )
                    )

        anova_by_run: list[RMSDRunANOVA] | None = None
        if len(summaries) >= 3:
            anova_by_run = []
            for run_label in run_labels:
                groups = [
                    summary.get_run(run_label).per_replicate_means
                    for summary in summaries_by_label.values()
                ]

                if len(groups) < 3 or any(len(group) < 2 for group in groups):
                    anova_by_run.append(
                        RMSDRunANOVA(
                            run_label=run_label,
                            significant=False,
                            testable=False,
                            note="Insufficient replicates (n < 2) for inferential statistics",
                        )
                    )
                    continue

                anova_result = one_way_anova(*groups)
                anova_by_run.append(
                    RMSDRunANOVA(
                        run_label=run_label,
                        f_statistic=anova_result.f_statistic,
                        p_value=anova_result.p_value,
                        significant=anova_result.significant,
                    )
                )

        fdr_alpha = getattr(ctx, "fdr_alpha", 0.05)
        self._apply_fdr_correction(pairwise_comparisons, anova_by_run, fdr_alpha)

        return RMSDComparisonResult(
            metric="mean_rmsd",
            name=ctx.name,
            n_runs=len(run_labels),
            run_labels=run_labels,
            control_label=effective_control,
            fdr_alpha=fdr_alpha,
            conditions=summaries,
            pairwise_comparisons=pairwise_comparisons,
            anova_by_run=anova_by_run,
            ranking_by_run=ranking_by_run,
            equilibration_time=ctx.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate RMSD comparison plots via delegated plotters."""
        comparison_path = ctx.comparison_path or self.comparison_result_path(ctx.results_dir)
        if not comparison_path.exists():
            logger.warning(
                "RMSD comparison result not found at %s; skipping plots", comparison_path
            )
            return []

        comparison_result = self._deserialize_comparison(comparison_path)
        if comparison_result is None:
            return []

        data, labels = self._build_plot_data(ctx)
        for label in labels:
            aggregated_dir = data[label]["aggregated_dir"]
            self._load_aggregated_result(
                aggregated_dir,
                settings=ctx.settings,
                condition_label=label,
            )

        try:
            from polyzymd.analyses.rmsd._plotters import (
                plot_rmsd_comparison_bars,
                plot_rmsd_convergence_diagnostics,
                plot_rmsd_timeseries,
            )
        except ImportError as exc:
            logger.warning("RMSD plotter module unavailable: %s", exc)
            return []

        plots: list[Path] = []
        plots.extend(plot_rmsd_timeseries(ctx, comparison_result))
        plots.extend(plot_rmsd_comparison_bars(ctx, comparison_result))
        plots.extend(plot_rmsd_convergence_diagnostics(ctx, comparison_result))
        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format RMSD comparison output via delegated formatter."""
        try:
            from polyzymd.analyses.rmsd._formatters import format_rmsd_comparison
        except ImportError as exc:
            logger.warning("RMSD formatter module unavailable: %s", exc)
            return super().format(result, output_format)

        return format_rmsd_comparison(result, output_format)

    def _compute_single_run(
        self,
        *,
        ctx: ReplicateContext,
        replicate: int,
        run: RMSDRunSettings,
        loader: TrajectoryLoader,
        config_hash: str,
        eq_value: float,
        eq_unit: str,
        eq_str: str,
        settings_tag: str,
        start_frame: int,
        n_frames_total: int,
        n_frames_used: int,
        timestep_ps: float,
    ) -> Any:
        """Compatibility shim for one RMSD run.

        This helper remains for focused unit tests while delegating the actual
        trajectory-native work to ``rmsd._runner``.
        """

        del ctx, config_hash, eq_value, eq_unit, eq_str, settings_tag, n_frames_used
        universe = loader.load_universe(replicate, cache=False)
        try:
            return compute_rmsd_run(
                universe=universe,
                run=run,
                start=start_frame,
                stop=n_frames_total,
                step=1,
                timestep_ps=timestep_ps,
            )
        except ValueError as exc:
            if "selection matched no atoms" not in str(exc):
                raise
            logger.warning("%s", exc)
            return None

    @staticmethod
    def _compare_run(
        *,
        run_label: str,
        condition_a: str,
        condition_b: str,
        run_a: Any,
        run_b: Any,
    ) -> Any:
        """Compare a single RMSD run between two conditions."""
        from polyzymd.analyses.rmsd._comparison_results import RMSDRunPairwiseComparison
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )
        from polyzymd.analyses.stats import interpret_direction

        run_a_values = run_a.per_replicate_means
        run_b_values = run_b.per_replicate_means
        pct_change = percent_change(run_a.mean_rmsd, run_b.mean_rmsd)

        direction = interpret_direction(
            pct_change,
            direction_labels=("stabilizing", "unchanged", "destabilizing"),
            threshold=1.0,
        )

        if len(run_a_values) < 2 or len(run_b_values) < 2:
            return RMSDRunPairwiseComparison(
                run_label=run_label,
                condition_a=condition_a,
                condition_b=condition_b,
                effect_interpretation="not_testable",
                direction=direction,
                significant=False,
                percent_change=pct_change,
                testable=False,
                note="Insufficient replicates (n < 2) for inferential statistics",
            )

        t_result = independent_ttest(run_a_values, run_b_values)
        d_result = cohens_d(run_a_values, run_b_values)

        return RMSDRunPairwiseComparison(
            run_label=run_label,
            condition_a=condition_a,
            condition_b=condition_b,
            t_statistic=t_result.t_statistic,
            p_value=t_result.p_value,
            cohens_d=d_result.cohens_d,
            effect_interpretation=d_result.interpretation,
            direction=direction,
            significant=t_result.significant,
            percent_change=pct_change,
        )

    @staticmethod
    def _apply_fdr_correction(
        pairwise: list[Any],
        anova_by_run: list[Any] | None,
        fdr_alpha: float,
    ) -> None:
        """Apply Benjamini-Hochberg FDR correction to pairwise and ANOVA p-values.

        Treats all pairwise comparisons as one family and ANOVA tests as
        a separate family.

        Parameters
        ----------
        pairwise : list
            Pairwise comparison results (mutated in place).
        anova_by_run : list or None
            ANOVA results (mutated in place).
        fdr_alpha : float
            FDR significance threshold.
        """

        def _set_corrected(result: Any, bh_result: Any) -> None:
            if hasattr(result, "p_value_adjusted"):
                result.p_value_adjusted = bh_result.adjusted_p_value
            result.significant = bh_result.significant

        apply_fdr_correction(
            pairwise,
            anova_by_run,
            fdr_alpha,
            get_p_value=lambda result: result.p_value if result.testable else None,
            set_corrected=lambda result, bh: _set_corrected(result, bh),
        )

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
        settings_tag: str,
    ) -> str:
        """Generate an aggregated RMSD filename."""
        eq_str = f"eq{first_result.equilibration_time:g}{first_result.equilibration_unit}"
        rep_str = Analysis._format_replicate_range(replicates)
        return f"rmsd_{rep_str}_{eq_str}_{settings_tag}.json"

    @staticmethod
    def _deserialize_comparison(path: Path) -> RMSDComparisonResult | None:
        """Load RMSD comparison result from disk."""
        try:
            from polyzymd.analyses.rmsd._comparison_results import RMSDComparisonResult
        except ImportError as exc:
            logger.warning("Cannot import RMSD comparison model: %s", exc)
            return None

        if not path.exists():
            return None
        return RMSDComparisonResult.load(path)
