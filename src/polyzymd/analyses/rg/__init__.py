"""Rg analysis plugin.

Computes per-frame Radius of Gyration for one or more named selections,
aggregates across replicates, performs per-run cross-condition comparisons,
and exposes plotting/formatting hooks.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Literal, Sequence

from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.rg._plot_settings import RgPlotSettings
from polyzymd.analyses.rg._results import RgAggregatedResult, RgResult
from polyzymd.analyses.rg._runner import RgReplicateRunner, compute_rg_run
from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
from polyzymd.analyses.shared.loader import TrajectoryLoader, parse_time_string
from polyzymd.analyses.shared.multi_run_comparison import (
    apply_fdr_correction,
    build_condition_pairs,
)
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.rg._comparison_results import RgComparisonResult

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _RgTrajectoryWindow:
    """Rg trajectory window that carries loader-derived file metadata.

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
    ) -> _RgTrajectoryWindow:
        """Build an Rg window wrapper from the shared trajectory window.

        Parameters
        ----------
        window : Any
            Shared trajectory window returned by the centralized resolver.
        trajectory_files : Sequence[Path]
            Trajectory files resolved by the existing loader instance.

        Returns
        -------
        _RgTrajectoryWindow
            Rg window wrapper that preserves run arguments and file metadata.
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


class RgRunSettings(BaseModel):
    """Settings for a single Rg run.

    Attributes
    ----------
    label : str
        Human-readable run label.
    selection : str
        MDAnalysis selection used for Rg calculation.
    calculation_mode : Literal["selection", "fragments"]
        Rg calculation mode for either whole-selection or per-fragment reduction.
    fragment_weighting : Literal["equal", "mass"]
        Fragment reduction weighting scheme when fragment mode is enabled.
    save_fragment_distribution : bool
        Whether to save per-fragment Rg distribution data in NPZ sidecar output.
    histogram_bins : int
        Number of histogram bins used for fragment distribution summaries.
    """

    label: str = Field(..., description="Human-readable run label")
    selection: str = Field(..., description="MDAnalysis selection for Rg calculation")
    calculation_mode: Literal["selection", "fragments"] = Field(
        default="selection",
        description='Rg calculation mode: "selection" (whole group) or "fragments" (per-fragment with reduction)',
    )
    fragment_weighting: Literal["equal", "mass"] = Field(
        default="equal",
        description='Fragment reduction weighting: "equal" (arithmetic mean) or "mass" (mass-weighted mean). Only used when calculation_mode="fragments".',
    )
    save_fragment_distribution: bool = Field(
        default=True,
        description="Save per-fragment Rg distribution data in NPZ sidecar. Only used when calculation_mode='fragments'.",
    )
    histogram_bins: int = Field(
        default=50,
        description="Number of histogram bins for fragment Rg distribution summaries.",
    )

    @field_validator("histogram_bins")
    @classmethod
    def validate_histogram_bins(cls, value: int) -> int:
        """Ensure histogram bin count is at least 2."""
        if value < 2:
            raise ValueError("histogram_bins must be >= 2")
        return value

    @model_validator(mode="after")
    def validate_fragment_options(self) -> RgRunSettings:
        """Validate fragment options against selected calculation mode."""
        if self.calculation_mode == "selection" and self.fragment_weighting != "equal":
            raise ValueError(
                "fragment_weighting is only meaningful when calculation_mode='fragments'"
            )
        return self


class RgSettings(BaseModel):
    """Top-level Rg settings.

    Attributes
    ----------
    runs : list[RgRunSettings]
        Named Rg runs to compute.
    """

    runs: list[RgRunSettings] = Field(
        default_factory=list,
        description="Rg runs to compute",
    )

    @field_validator("runs", mode="after")
    @classmethod
    def validate_runs_nonempty(cls, value: list[RgRunSettings]) -> list[RgRunSettings]:
        """Ensure at least one run exists."""
        if not value:
            raise ValueError("At least one Rg run must be defined")
        return value

    @field_validator("runs", mode="after")
    @classmethod
    def validate_unique_labels(cls, value: list[RgRunSettings]) -> list[RgRunSettings]:
        """Ensure run labels are unique."""
        labels = [run.label for run in value]
        if len(labels) != len(set(labels)):
            raise ValueError("Rg run labels must be unique")
        return value


class RgAnalysis(Analysis):
    """Rg analysis plugin using a multi-run comparison model."""

    name: ClassVar[str] = "rg"
    min_replicates: ClassVar[int] = 1
    Settings: ClassVar[type] = RgSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = RgPlotSettings
    AggregatedResultClass: ClassVar[type] = RgAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = RgResult
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel) -> str:
        """Build a short cache tag from analysis settings.

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
        settings: RgSettings,
        *,
        condition_label: str | None = None,
        source: Path | None = None,
    ) -> RgAggregatedResult:
        """Coerce an aggregated result and validate its settings identity.

        Parameters
        ----------
        result : Any
            Aggregated result loaded from disk or supplied in memory.
        settings : RgSettings
            Current Rg settings for comparison or plotting.
        condition_label : str | None, optional
            Condition label for error reporting.
        source : Path | None, optional
            Source file path for diagnostics.

        Returns
        -------
        RgAggregatedResult
            Validated aggregated result.

        Raises
        ------
        ValueError
            Raised when the aggregated result is missing a settings
            fingerprint or was computed with different settings.
        """
        if isinstance(result, dict):
            result = RgAggregatedResult.model_validate(result)

        if not isinstance(result, RgAggregatedResult):
            raise TypeError(
                f"Rg aggregated result loader expected RgAggregatedResult, got "
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
                "Aggregated Rg result"
                f"{condition_text} is missing a settings fingerprint{source_text}. "
                "Legacy Rg aggregated caches are not compatible with "
                "settings-sensitive compare/plot loading. Recompute the condition before "
                "comparing or plotting."
            )
        if stored_fingerprint != current_fingerprint:
            raise ValueError(
                "Aggregated Rg result"
                f"{condition_text} was computed with settings fingerprint "
                f"{stored_fingerprint}, but current settings require {current_fingerprint}"
                f"{source_text}. Recompute the condition or clear stale caches before "
                "comparing or plotting."
            )
        return result

    def _resolve_aggregated_result_path(self, aggregated_dir: Path) -> Path | None:
        """Resolve the aggregated Rg result path.

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
        settings: RgSettings | None = None,
        condition_label: str | None = None,
    ) -> Any:
        """Load and optionally validate an aggregated Rg result.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.
        settings : RgSettings | None, optional
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
            Per-replicate Rg results.
        configured_run_labels : Sequence[str]
            Run labels defined in the Rg settings.

        Raises
        ------
        ValueError
            Raised when configured replicates or runs are missing from the
            aggregation inputs.
        """
        if not results:
            raise ValueError(
                f"Rg aggregation for condition '{ctx.condition.label}' requires at least one "
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
                f"Rg aggregation for condition '{ctx.condition.label}' is incomplete. Expected "
                f"replicate results for {expected_replicates}; observed {sorted(observed_replicates)} "
                f"({detail_text}). Recompute missing replicates or clear stale caches before "
                "aggregating."
            )

        for run_label in configured_run_labels:
            missing_run_replicates: list[int] = []
            duplicate_run_replicates: list[int] = []
            for result in results:
                replicate = getattr(result, "replicate", None)
                matches = [
                    run_result
                    for run_result in result.run_results
                    if run_result.run_label == run_label
                ]
                if not matches:
                    if replicate is not None:
                        missing_run_replicates.append(replicate)
                    continue
                if len(matches) > 1 and replicate is not None:
                    duplicate_run_replicates.append(replicate)

            if missing_run_replicates:
                raise ValueError(
                    f"Configured Rg run '{run_label}' is missing replicate entries in condition "
                    f"'{ctx.condition.label}'. Missing replicates: "
                    f"{sorted(missing_run_replicates)}. Recompute missing replicates or clear "
                    "stale caches before aggregating."
                )

            if duplicate_run_replicates:
                raise ValueError(
                    f"Configured Rg run '{run_label}' has duplicate replicate entries in "
                    f"condition '{ctx.condition.label}' for replicates "
                    f"{sorted(duplicate_run_replicates)}. Clear stale caches and recompute "
                    "before aggregating."
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
            Per-replicate Rg results that already passed completeness checks.

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
                    f"Rg aggregation for condition '{ctx.condition.label}' encountered a "
                    "replicate result without a replicate identifier while normalizing "
                    "aggregate inputs."
                )
            if replicate in ordered_results:
                raise ValueError(
                    f"Rg aggregation for condition '{ctx.condition.label}' encountered "
                    f"duplicate replicate {replicate} while normalizing aggregate inputs."
                )
            ordered_results[replicate] = result

        missing_replicates = [
            replicate for replicate in ctx.replicates if replicate not in ordered_results
        ]
        if missing_replicates:
            raise ValueError(
                f"Rg aggregation for condition '{ctx.condition.label}' cannot order aggregate "
                f"inputs because replicates {missing_replicates} are missing."
            )

        return [ordered_results[replicate] for replicate in ctx.replicates]

    @staticmethod
    def _validate_aggregated_result_completeness(
        condition: Any,
        agg_result: RgAggregatedResult,
        configured_run_labels: Sequence[str],
    ) -> None:
        """Validate that an aggregated Rg result is complete for comparison.

        Parameters
        ----------
        condition : Any
            Condition associated with the aggregated result.
        agg_result : RgAggregatedResult
            Aggregated Rg result to validate.
        configured_run_labels : Sequence[str]
            Run labels defined in the Rg settings.

        Raises
        ------
        ValueError
            Raised when the aggregated result omits configured runs or contains
            incomplete replicate data.
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
                f"Aggregated Rg result for condition '{condition.label}' is incomplete: "
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
                f"Aggregated Rg result for condition '{condition.label}' has incomplete "
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
            }
            if run_result.per_replicate_mean_fragments_per_frame is not None:
                counts["per_replicate_mean_fragments_per_frame"] = len(
                    run_result.per_replicate_mean_fragments_per_frame
                )
            mismatched_fields = {
                name: count for name, count in counts.items() if count != len(expected_replicates)
            }
            if (
                run_result.n_replicates != len(expected_replicates)
                or run_replicates != expected_replicates
            ):
                raise ValueError(
                    f"Aggregated Rg run '{run_result.run_label}' for condition "
                    f"'{condition.label}' has incomplete replicate metadata. Expected "
                    f"replicates {expected_replicates}, found {run_replicates} with "
                    f"n_replicates={run_result.n_replicates}. Recompute the condition or clear "
                    "stale caches before comparing."
                )

            if mismatched_fields:
                raise ValueError(
                    f"Aggregated Rg run '{run_result.run_label}' for condition "
                    f"'{condition.label}' has incomplete replicate values: {mismatched_fields}. "
                    f"Expected {len(expected_replicates)} entries per field. Recompute the "
                    "condition or clear stale caches before comparing."
                )

    @staticmethod
    def _describe_sidecar_aggregation_contract(run: RgRunSettings) -> str:
        """Describe sidecar-backed aggregate outputs required for a configured run.

        Parameters
        ----------
        run : RgRunSettings
            Configured Rg run definition.

        Returns
        -------
        str
            Human-readable description of the aggregate outputs that require
            NPZ sidecars for every replicate.
        """
        required_outputs = ["reduced-series histograms"]
        if run.calculation_mode == "fragments" and run.save_fragment_distribution:
            required_outputs.append("fragment distributions")
        return " and ".join(required_outputs)

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Run Rg for all configured runs for a single replicate.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        RgResult
            Per-replicate Rg result containing all run outputs.
        """
        settings = ctx.settings

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        eq_str = f"eq{eq_value:g}{eq_unit}"
        settings_tag = self._make_settings_cache_tag(settings)
        result_file = ctx.output_dir / f"rg_{eq_str}_{settings_tag}.json"

        cached = self._check_cache(
            RgResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=ctx.sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

        result = super().run_replicate(ctx, replicate)
        result.save(result_file)
        logger.info("Saved Rg result to %s", result_file)
        return result

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return the Rg loader class for the shared runner seam.

        Returns
        -------
        type[Any]
            Loader class patched by Rg unit tests.
        """

        return TrajectoryLoader

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the Rg window and retain trajectory file metadata.

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
        return _RgTrajectoryWindow.from_window(window, traj_info.trajectory_files)

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any:
        """Build the runner-backed Rg execution object.

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

        return RgReplicateRunner(
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
        """Serialize runner output into the legacy Rg result schema.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        runner : Any
            Executed Rg runner.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        RgResult
            Cache-compatible per-replicate Rg result.
        """
        import numpy as np

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rg._results import RgResult, RgRunResult

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        eq_str = f"eq{eq_value:g}{eq_unit}"
        settings_tag = self._make_settings_cache_tag(ctx.settings)
        config_hash = compute_config_hash(ctx.sim_config)
        trajectory_files = getattr(window, "trajectory_files", ())

        run_results: list[RgRunResult] = []
        for payload in runner.results.run_payloads:
            npz_filename = f"rg_{payload.run_label}_{eq_str}_{settings_tag}_timeseries.npz"
            npz_path = ctx.output_dir / npz_filename
            npz_data: dict[str, np.ndarray] = {
                "rg_values": payload.rg_values,
                "time_ns": payload.time_ns,
                "frames": payload.frames,
            }
            if payload.fragment_rg_values is not None:
                npz_data["fragment_rg_values"] = payload.fragment_rg_values
            if payload.fragment_counts_per_frame is not None:
                npz_data["fragment_counts_per_frame"] = payload.fragment_counts_per_frame
            if payload.fragment_masses is not None:
                npz_data["fragment_masses"] = payload.fragment_masses
            np.savez_compressed(npz_path, **npz_data)

            run_results.append(
                RgRunResult(
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
                    calculation_mode=payload.calculation_mode,
                    fragment_weighting=payload.fragment_weighting,
                    mean_rg=payload.mean_rg,
                    std_rg=payload.std_rg,
                    median_rg=payload.median_rg,
                    min_rg=payload.min_rg,
                    max_rg=payload.max_rg,
                    final_rg=payload.final_rg,
                    sem_rg=payload.sem_rg,
                    correlation_time_unit=payload.correlation_time_unit,
                    statistical_inefficiency=payload.statistical_inefficiency,
                    autocorrelation_warning=payload.autocorrelation_warning,
                    n_frames_total=runner.results.n_frames_total,
                    n_frames_used=window.n_frames_selected,
                    npz_path=str(npz_path),
                    time_unit="ns",
                    timestep_ps=window.timestep_ps,
                    **payload.frag_metadata,
                )
            )

        return RgResult(
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
        """Aggregate Rg results across replicates for one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[RgResult]
            Per-replicate Rg results.

        Returns
        -------
        RgAggregatedResult
            Aggregated Rg result for all configured runs.
        """
        import numpy as np

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rg._results import RgAggregatedResult, RgRunAggregatedResult

        run_labels = [run.label for run in ctx.settings.runs]
        self._validate_aggregate_input_completeness(ctx, results, run_labels)
        ordered_results = self._order_results_by_replicate(ctx, results)
        replicate_order = list(ctx.replicates)

        first = ordered_results[0]

        if len(ctx.replicates) == 1:
            logger.warning(
                "Only one replicate available for Rg aggregation in condition '%s'; "
                "replicate-level SEM is reported as 0.0",
                ctx.condition.label,
            )

        aggregated_runs: list[RgRunAggregatedResult] = []
        for run in ctx.settings.runs:
            run_label = run.label
            run_entries = []
            for result in ordered_results:
                matches = [
                    run_result
                    for run_result in result.run_results
                    if run_result.run_label == run_label
                ]
                if len(matches) != 1:
                    raise ValueError(
                        f"Configured Rg run '{run_label}' has invalid aggregate inputs in "
                        f"condition '{ctx.condition.label}'. Expected one entry per replicate, "
                        f"found {len(matches)} for replicate {result.replicate}."
                    )
                run_entries.append(matches[0])

            per_means = [entry.mean_rg for entry in run_entries]
            per_stds = [entry.std_rg for entry in run_entries]
            per_medians = [entry.median_rg for entry in run_entries]

            mean_stats = compute_sem(per_means)
            overall_median = float(np.median(np.asarray(per_medians, dtype=np.float64)))

            template = run_entries[0]

            per_replicate_mean_fragments_per_frame: list[float] | None = None
            overall_mean_fragments_per_frame: float | None = None
            fragment_histogram_edges: list[float] | None = None
            fragment_histogram_density_mean: list[float] | None = None
            fragment_histogram_density_sem: list[float] | None = None
            reduced_histogram_edges: list[float] | None = None
            reduced_histogram_density_mean: list[float] | None = None
            reduced_histogram_density_sem: list[float] | None = None

            all_reduced_rg_per_rep: list[np.ndarray] = []
            all_fragment_rg_per_rep: list[np.ndarray] = []
            requires_fragment_distribution = (
                run.calculation_mode == "fragments" and run.save_fragment_distribution
            )
            required_sidecar_outputs = self._describe_sidecar_aggregation_contract(run)
            missing_sidecar_replicates = [
                replicate
                for replicate, entry in zip(replicate_order, run_entries)
                if entry.npz_path is None
            ]
            if missing_sidecar_replicates:
                raise ValueError(
                    f"Rg aggregation for run '{run_label}' in condition "
                    f"'{ctx.condition.label}' requires NPZ sidecar metadata for "
                    f"{required_sidecar_outputs} in replicates "
                    f"{missing_sidecar_replicates}. Recompute those replicates or clear "
                    "stale caches before aggregating."
                )

            for replicate, entry in zip(replicate_order, run_entries):
                npz_path_value = entry.npz_path
                if npz_path_value is None:
                    raise ValueError(
                        f"Rg aggregation for run '{run_label}' in condition "
                        f"'{ctx.condition.label}' is missing NPZ sidecar metadata for "
                        f"replicate {replicate}."
                    )
                npz_path = Path(npz_path_value)
                if not npz_path.exists():
                    raise ValueError(
                        f"Rg aggregation for run '{run_label}' in condition "
                        f"'{ctx.condition.label}' expected NPZ sidecar {npz_path} for "
                        f"replicate {replicate}, but the file is missing. Recompute that "
                        "replicate or clear stale caches before aggregating."
                    )

                with np.load(npz_path) as npz_data:
                    if "rg_values" not in npz_data:
                        raise ValueError(
                            f"Rg aggregation for run '{run_label}' in condition "
                            f"'{ctx.condition.label}' expected 'rg_values' in NPZ sidecar "
                            f"{npz_path} for replicate {replicate}."
                        )

                    reduced_values = np.asarray(npz_data["rg_values"], dtype=np.float64)
                    if reduced_values.size == 0:
                        raise ValueError(
                            f"Rg aggregation for run '{run_label}' in condition "
                            f"'{ctx.condition.label}' found empty reduced-series data in "
                            f"NPZ sidecar {npz_path} for replicate {replicate}."
                        )
                    all_reduced_rg_per_rep.append(reduced_values)

                    if requires_fragment_distribution:
                        if "fragment_rg_values" not in npz_data:
                            raise ValueError(
                                f"Rg aggregation for run '{run_label}' in condition "
                                f"'{ctx.condition.label}' expected 'fragment_rg_values' in "
                                f"NPZ sidecar {npz_path} for replicate {replicate}."
                            )

                        fragment_values = np.asarray(
                            npz_data["fragment_rg_values"],
                            dtype=np.float64,
                        )
                        if fragment_values.size == 0:
                            raise ValueError(
                                f"Rg aggregation for run '{run_label}' in condition "
                                f"'{ctx.condition.label}' found empty fragment distribution "
                                f"data in NPZ sidecar {npz_path} for replicate {replicate}."
                            )
                        all_fragment_rg_per_rep.append(fragment_values)

            if template.calculation_mode == "fragments":
                missing_fragment_metrics = [
                    replicate
                    for replicate, entry in zip(replicate_order, run_entries)
                    if entry.mean_fragments_per_frame is None
                ]
                if missing_fragment_metrics:
                    raise ValueError(
                        f"Rg aggregation for run '{run_label}' in condition "
                        f"'{ctx.condition.label}' is missing mean_fragments_per_frame for "
                        f"replicates {missing_fragment_metrics}. Recompute those replicates or "
                        "clear stale caches before aggregating."
                    )

                per_replicate_mean_fragments_per_frame = [
                    float(entry.mean_fragments_per_frame) for entry in run_entries
                ]
                overall_mean_fragments_per_frame = float(
                    np.mean(
                        np.asarray(
                            per_replicate_mean_fragments_per_frame,
                            dtype=np.float64,
                        )
                    )
                )

                if all_fragment_rg_per_rep:
                    pooled_fragment_rg = np.concatenate(all_fragment_rg_per_rep)
                    fragment_min = float(np.min(pooled_fragment_rg))
                    fragment_max = float(np.max(pooled_fragment_rg))
                    if fragment_min == fragment_max:
                        fragment_min -= 1.0e-6
                        fragment_max += 1.0e-6

                    fragment_edges = np.linspace(
                        fragment_min,
                        fragment_max,
                        run.histogram_bins + 1,
                        dtype=np.float64,
                    )
                    fragment_densities = np.asarray(
                        [
                            np.histogram(rep_data, bins=fragment_edges, density=True)[0]
                            for rep_data in all_fragment_rg_per_rep
                        ],
                        dtype=np.float64,
                    )
                    fragment_histogram_edges = fragment_edges.tolist()
                    fragment_histogram_density_mean = np.mean(fragment_densities, axis=0).tolist()
                    if len(all_fragment_rg_per_rep) > 1:
                        fragment_histogram_density_sem = (
                            np.std(fragment_densities, axis=0, ddof=1)
                            / np.sqrt(len(all_fragment_rg_per_rep))
                        ).tolist()
                    else:
                        fragment_histogram_density_sem = [
                            0.0 for _ in range(len(fragment_histogram_density_mean))
                        ]

            # --- Reduced-series histogram (both selection and fragment modes) ---
            if all_reduced_rg_per_rep:
                pooled_reduced_rg = np.concatenate(all_reduced_rg_per_rep)
                reduced_min = float(np.min(pooled_reduced_rg))
                reduced_max = float(np.max(pooled_reduced_rg))
                if reduced_min == reduced_max:
                    reduced_min -= 1.0e-6
                    reduced_max += 1.0e-6

                reduced_edges = np.linspace(
                    reduced_min,
                    reduced_max,
                    run.histogram_bins + 1,
                    dtype=np.float64,
                )
                reduced_densities = np.asarray(
                    [
                        np.histogram(rep_data, bins=reduced_edges, density=True)[0]
                        for rep_data in all_reduced_rg_per_rep
                    ],
                    dtype=np.float64,
                )
                reduced_histogram_edges = reduced_edges.tolist()
                reduced_histogram_density_mean = np.mean(reduced_densities, axis=0).tolist()
                if len(all_reduced_rg_per_rep) > 1:
                    reduced_histogram_density_sem = (
                        np.std(reduced_densities, axis=0, ddof=1)
                        / np.sqrt(len(all_reduced_rg_per_rep))
                    ).tolist()
                else:
                    reduced_histogram_density_sem = [
                        0.0 for _ in range(len(reduced_histogram_density_mean))
                    ]

            aggregated_runs.append(
                RgRunAggregatedResult(
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
                    overall_mean=mean_stats.mean,
                    overall_sem=mean_stats.sem,
                    overall_median=overall_median,
                    per_replicate_means=per_means,
                    per_replicate_stds=per_stds,
                    per_replicate_medians=per_medians,
                    calculation_mode=template.calculation_mode,
                    fragment_weighting=template.fragment_weighting,
                    overall_mean_fragments_per_frame=overall_mean_fragments_per_frame,
                    per_replicate_mean_fragments_per_frame=per_replicate_mean_fragments_per_frame,
                    fragment_histogram_edges=fragment_histogram_edges,
                    fragment_histogram_density_mean=fragment_histogram_density_mean,
                    fragment_histogram_density_sem=fragment_histogram_density_sem,
                    reduced_histogram_edges=reduced_histogram_edges,
                    reduced_histogram_density_mean=reduced_histogram_density_mean,
                    reduced_histogram_density_sem=reduced_histogram_density_sem,
                )
            )

        agg_result = RgAggregatedResult(
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
        logger.info("Saved aggregated Rg result to %s", target_path)

        return agg_result

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare Rg runs across conditions.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided comparison context.

        Returns
        -------
        RgComparisonResult | None
            Comparison result, or ``None`` if no conditions have data.
        """
        from polyzymd import __version__
        from polyzymd.analyses.rg._comparison_results import (
            RgComparisonResult,
            RgConditionSummary,
            RgRunANOVA,
            RgRunPairwiseComparison,
            RgRunSummary,
        )
        from polyzymd.analyses.shared.inferential_statistics import one_way_anova

        configured_run_labels = [run.label for run in ctx.settings.runs]
        summaries: list[RgConditionSummary] = []
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
                    f"Rg comparison requires an aggregated result for condition "
                    f"'{condition.label}'. Recompute the condition or clear stale caches "
                    "before comparing."
                )

            self._validate_aggregated_result_completeness(
                condition,
                agg_result,
                configured_run_labels,
            )

            run_summaries = [
                RgRunSummary(
                    label=run_result.run_label,
                    selection=run_result.selection,
                    mean_rg=run_result.overall_mean,
                    sem_rg=run_result.overall_sem,
                    per_replicate_means=run_result.per_replicate_means,
                    calculation_mode=run_result.calculation_mode,
                    fragment_weighting=run_result.fragment_weighting,
                    mean_fragments_per_frame=run_result.overall_mean_fragments_per_frame,
                )
                for run_result in agg_result.run_results
            ]

            summaries.append(
                RgConditionSummary(
                    label=condition.label,
                    config_path=str(condition.config_path),
                    n_replicates=agg_result.n_replicates,
                    run_summaries=run_summaries,
                )
            )

        if not summaries:
            logger.warning("Rg comparison skipped because no conditions have results")
            return None

        summaries_by_label = {summary.label: summary for summary in summaries}

        effective_control = ctx.effective_control

        run_labels: list[str] = []
        ranking_by_run: dict[str, list[str]] = {}
        pairwise_comparisons: list[RgRunPairwiseComparison] = []
        anova_by_run: list[RgRunANOVA] | None = []

        for run_label in configured_run_labels:
            run_labels.append(run_label)
            ranked = sorted(
                summaries_by_label,
                key=lambda label: summaries_by_label[label].get_run(run_label).mean_rg,
            )
            ranking_by_run[run_label] = ranked

            if len(summaries_by_label) >= 2:
                condition_pairs = build_condition_pairs(
                    list(summaries_by_label.keys()),
                    effective_control,
                    on_control_missing="all_pairs",
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

            if len(summaries_by_label) >= 3:
                groups = [
                    summary.get_run(run_label).per_replicate_means
                    for summary in summaries_by_label.values()
                ]
                if any(len(group) < 2 for group in groups):
                    anova_by_run.append(
                        RgRunANOVA(
                            run_label=run_label,
                            significant=False,
                            testable=False,
                            note="Insufficient replicates (n < 2) for inferential statistics",
                        )
                    )
                    continue

                anova_result = one_way_anova(*groups)
                anova_by_run.append(
                    RgRunANOVA(
                        run_label=run_label,
                        f_statistic=anova_result.f_statistic,
                        p_value=anova_result.p_value,
                        significant=anova_result.significant,
                    )
                )

        if not anova_by_run:
            anova_by_run = None

        fdr_alpha = getattr(ctx, "fdr_alpha", 0.05)
        self._apply_fdr_correction(pairwise_comparisons, anova_by_run, fdr_alpha)

        return RgComparisonResult(
            metric="mean_rg",
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
        """Generate Rg comparison plots via delegated plotters."""
        comparison_path = ctx.comparison_path or self.comparison_result_path(ctx.results_dir)
        if not comparison_path.exists():
            logger.warning("Rg comparison result not found at %s; skipping plots", comparison_path)
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
            from polyzymd.analyses.rg._plotters import (
                plot_rg_comparison_bars,
                plot_rg_distributions,
                plot_rg_timeseries,
            )
        except ImportError as exc:
            logger.warning("Rg plotter module unavailable: %s", exc)
            return []

        plots: list[Path] = []
        plots.extend(plot_rg_timeseries(ctx, comparison_result))
        plots.extend(plot_rg_comparison_bars(ctx, comparison_result))
        plots.extend(plot_rg_distributions(ctx, comparison_result))
        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format Rg comparison output via delegated formatter."""
        try:
            from polyzymd.analyses.rg._formatters import format_rg_comparison
        except ImportError as exc:
            logger.warning("Rg formatter module unavailable: %s", exc)
            return super().format(result, output_format)

        return format_rg_comparison(result, output_format)

    def _compute_single_run(
        self,
        *,
        ctx: ReplicateContext,
        replicate: int,
        run: RgRunSettings,
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
        """Compatibility shim for one Rg run.

        This helper remains for focused unit tests while delegating the actual
        trajectory-native work to ``rg._runner``. Empty selections now fail
        immediately so the replicate exits before any partial result files are
        written.
        """

        del ctx, config_hash, eq_value, eq_unit, eq_str, settings_tag, n_frames_used
        universe = loader.load_universe(replicate, cache=False)
        return compute_rg_run(
            universe=universe,
            run=run,
            replicate=replicate,
            start=start_frame,
            stop=n_frames_total,
            step=1,
            timestep_ps=timestep_ps,
        )

    @staticmethod
    def _has_run_summary(summary: Any, run_label: str) -> bool:
        """Check whether a condition summary includes a given run

        Parameters
        ----------
        summary : Any
            Condition summary object
        run_label : str
            Run label to query

        Returns
        -------
        bool
            ``True`` when the run exists in the condition summary
        """
        try:
            summary.get_run(run_label)
            return True
        except KeyError:
            return False

    @staticmethod
    def _compare_run(
        *,
        run_label: str,
        condition_a: str,
        condition_b: str,
        run_a: Any,
        run_b: Any,
    ) -> Any:
        """Compare a single Rg run between two conditions.

        Parameters
        ----------
        run_label : str
            Label of the run being compared.
        condition_a : str
            First condition label.
        condition_b : str
            Second condition label.
        run_a : Any
            Condition A run summary object.
        run_b : Any
            Condition B run summary object.

        Returns
        -------
        RgRunPairwiseComparison
            Pairwise statistics and direction for this run.
        """
        from polyzymd.analyses.rg._comparison_results import RgRunPairwiseComparison
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )
        from polyzymd.analyses.stats import interpret_direction

        run_a_values = run_a.per_replicate_means
        run_b_values = run_b.per_replicate_means
        pct_change = percent_change(run_a.mean_rg, run_b.mean_rg)

        direction = interpret_direction(
            pct_change,
            direction_labels=("compaction", "unchanged", "expansion"),
            threshold=1.0,
        )

        if len(run_a_values) < 2 or len(run_b_values) < 2:
            return RgRunPairwiseComparison(
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

        return RgRunPairwiseComparison(
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
        """Generate an aggregated Rg filename."""
        eq_str = f"eq{first_result.equilibration_time:g}{first_result.equilibration_unit}"
        rep_str = Analysis._format_replicate_range(replicates)
        return f"rg_{rep_str}_{eq_str}_{settings_tag}.json"

    @staticmethod
    def _deserialize_comparison(path: Path) -> RgComparisonResult | None:
        """Load Rg comparison result from disk."""
        try:
            from polyzymd.analyses.rg._comparison_results import RgComparisonResult
        except ImportError as exc:
            logger.warning("Cannot import Rg comparison model: %s", exc)
            return None

        if not path.exists():
            return None
        return RgComparisonResult.load(path)
