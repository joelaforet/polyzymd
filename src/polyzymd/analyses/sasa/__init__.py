"""SASA analysis plugin."""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence, cast

from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    PlotContext,
    ReplicateContext,
    SlurmResourceHint,
)
from polyzymd.analyses.sasa._plot_settings import SASAPlotSettings
from polyzymd.analyses.sasa._results import SASAAggregatedResult, SASAResult
from polyzymd.analyses.sasa._runner import SASAReplicateRunner
from polyzymd.analyses.shared.config_hash import compute_config_hash
from polyzymd.analyses.shared.loader import parse_time_string
from polyzymd.analyses.shared.multi_run_comparison import (
    apply_fdr_correction,
    build_condition_pairs,
    filter_summaries_with_run,
)
from polyzymd.analyses.shared.sasa import compute_sasa, save_sasa_artifacts
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.sasa._comparison_results import SASAComparisonResult

LOGGER = logging.getLogger(__name__)
NOT_TESTABLE_SINGLETON_NOTE = "Inferential statistics require at least two replicates per condition."


@dataclass(frozen=True)
class _SASATrajectoryWindow:
    """SASA trajectory window that carries loader-derived file metadata."""

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
    ) -> _SASATrajectoryWindow:
        """Build an SASA window wrapper from the shared trajectory window."""

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


class SASARunSettings(BaseModel):
    """Settings for a single SASA run."""

    label: str = Field(..., description="Human-readable run label")
    target_selection: str = Field(..., description="Selection of atoms whose SASA is reported")
    context_selection: str | None = Field(
        default=None,
        description="Selection of atoms considered during SASA computation",
    )
    stride: int = Field(default=1, description="Frame stride (1 = every frame)")

    @field_validator("label", mode="after")
    @classmethod
    def validate_label(cls, value: str) -> str:
        """Validate run label content."""
        stripped = value.strip()
        if not stripped:
            raise ValueError("SASA run label must not be blank")
        if "/" in stripped or "\\" in stripped:
            raise ValueError("SASA run label must not contain '/' or '\\'")
        return stripped

    @field_validator("target_selection", mode="after")
    @classmethod
    def validate_target_selection(cls, value: str) -> str:
        """Validate target selection content."""
        stripped = value.strip()
        if not stripped:
            raise ValueError("SASA target_selection must not be blank")
        return stripped

    @field_validator("context_selection", mode="after")
    @classmethod
    def normalize_context_selection(cls, value: str | None) -> str | None:
        """Normalize context selection string."""
        if value is None:
            return None
        stripped = value.strip()
        return stripped if stripped else None

    @model_validator(mode="after")
    def apply_context_default(self) -> SASARunSettings:
        """Default context selection to target selection when omitted."""
        if self.context_selection is None:
            self.context_selection = self.target_selection
        return self

    @field_validator("stride", mode="after")
    @classmethod
    def validate_stride(cls, value: int) -> int:
        """Validate stride value."""
        if value <= 0:
            raise ValueError("stride must be >= 1")
        return value


class SASASettings(BaseModel):
    """Top-level SASA settings."""

    runs: list[SASARunSettings] = Field(default_factory=list, description="SASA runs to compute")
    probe_radius_nm: float = Field(default=0.14, description="Shrake-Rupley probe radius in nm")
    n_sphere_points: int = Field(default=960, description="Shrake-Rupley sphere point count")
    chunk_size: int = Field(
        default=100,
        description="Frames per chunk for memory-efficient SASA computation",
    )

    @field_validator("runs", mode="after")
    @classmethod
    def validate_runs_nonempty(cls, value: list[SASARunSettings]) -> list[SASARunSettings]:
        """Ensure at least one SASA run exists."""
        if not value:
            raise ValueError("At least one SASA run must be defined")
        labels = [run.label for run in value]
        if len(labels) != len(set(labels)):
            raise ValueError("SASA run labels must be unique")
        return value

    @field_validator("probe_radius_nm", mode="after")
    @classmethod
    def validate_probe_radius(cls, value: float) -> float:
        """Validate probe radius."""
        if value <= 0.0:
            raise ValueError("probe_radius_nm must be > 0")
        return value

    @field_validator("n_sphere_points", mode="after")
    @classmethod
    def validate_sphere_points(cls, value: int) -> int:
        """Validate sphere points."""
        if value < 100:
            raise ValueError("n_sphere_points must be >= 100")
        return value

    @field_validator("chunk_size", mode="after")
    @classmethod
    def validate_chunk_size(cls, value: int) -> int:
        """Validate chunk size."""
        if value <= 0:
            raise ValueError("chunk_size must be >= 1")
        return value


class SASAAnalysis(Analysis):
    """SASA analysis plugin using a multi-run comparison model."""

    name: ClassVar[str] = "sasa"
    min_replicates: ClassVar[int] = 1
    Settings: ClassVar[type] = SASASettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = SASAPlotSettings
    AggregatedResultClass: ClassVar[type | None] = SASAAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = SASAResult
    execution_cost_hint: ClassVar[str] = "high"
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = SlurmResourceHint(
        mem="8G", time="02:00:00"
    )
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    # SASA is a mean-based observable (all frames contribute, SEM corrected via N_eff)

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Run SASA for all configured runs for a single replicate."""

        settings = cast(SASASettings, ctx.settings)
        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        eq_str = f"eq{eq_value:g}{eq_unit}"
        settings_token = self._settings_cache_token(settings)
        result_file = ctx.output_dir / f"sasa_{eq_str}_{settings_token}.json"

        cached = self._check_cache(
            SASAResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=ctx.sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

        result = super().run_replicate(ctx, replicate)
        result.save(result_file)
        return result

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the SASA window and retain trajectory file metadata."""

        window = super().get_trajectory_window(ctx, replicate, loader, universe)
        traj_info = loader.get_trajectory_info(replicate)
        return _SASATrajectoryWindow.from_window(window, traj_info.trajectory_files)

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any:
        """Build the runner-backed SASA execution object."""

        from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

        del replicate
        settings = cast(SASASettings, ctx.settings)
        return SASAReplicateRunner(
            universe=universe,
            runs=settings.runs,
            probe_radius_nm=settings.probe_radius_nm,
            n_sphere_points=settings.n_sphere_points,
            chunk_size=settings.chunk_size,
            timestep_ps=window.timestep_ps,
            output_dir=ctx.output_dir,
            equilibration=ctx.equilibration,
            trajectory_files=getattr(window, "trajectory_files", ()),
            compute_sasa_func=compute_sasa,
            save_sasa_artifacts_func=save_sasa_artifacts,
            estimate_correlation_time_func=estimate_correlation_time,
            run_cache_token_func=self._run_cache_token,
        )

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Serialize runner output into the legacy SASA result schema."""

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.sasa._results import SASARunResult

        settings = cast(SASASettings, ctx.settings)
        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        config_hash = compute_config_hash(ctx.sim_config)
        payload = runner.results

        run_results = [
            SASARunResult(
                config_hash=config_hash,
                polyzymd_version=get_polyzymd_version(),
                replicate=replicate,
                equilibration_time=eq_value,
                equilibration_unit=eq_unit,
                selection_string=run_payload.target_selection,
                correlation_time=run_payload.correlation_time,
                correlation_time_unit=run_payload.correlation_time_unit,
                n_independent_frames=run_payload.n_independent_frames,
                statistical_inefficiency=run_payload.statistical_inefficiency,
                autocorrelation_warning=run_payload.autocorrelation_warning,
                run_label=run_payload.run_label,
                target_selection=run_payload.target_selection,
                context_selection=run_payload.context_selection,
                mean_sasa=run_payload.mean_sasa,
                std_sasa=run_payload.std_sasa,
                median_sasa=run_payload.median_sasa,
                min_sasa=run_payload.min_sasa,
                max_sasa=run_payload.max_sasa,
                final_sasa=run_payload.final_sasa,
                sem_sasa=run_payload.sem_sasa,
                n_frames_total=run_payload.n_frames_total,
                n_frames_used=run_payload.n_frames_used,
                n_target_atoms=run_payload.n_target_atoms,
                n_context_atoms=run_payload.n_context_atoms,
                n_target_residues=run_payload.n_target_residues,
                zero_atom_selection=run_payload.zero_atom_selection,
                raw_npz_path=run_payload.raw_npz_path,
                raw_metadata_path=run_payload.raw_metadata_path,
                npz_path=run_payload.npz_path,
                metadata_path=run_payload.metadata_path,
                time_unit=run_payload.time_unit,
                timestep_ps=run_payload.timestep_ps,
                raw_timestep_ps=run_payload.raw_timestep_ps,
                frame_stride=run_payload.frame_stride,
            )
            for run_payload in payload.run_payloads
        ]

        del window
        return SASAResult(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string="; ".join(
                f"{run.target_selection}|{run.context_selection or run.target_selection}"
                for run in settings.runs
            ),
            run_results=run_results,
            n_frames_total=payload.n_frames_total,
            n_frames_used=payload.n_frames_used,
            trajectory_files=[str(path) for path in payload.trajectory_files],
        )

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate SASA results across replicates for one condition."""
        import numpy as np

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.sasa._results import SASAAggregatedResult, SASARunAggregatedResult

        first = results[0]
        aggregated_runs: list[SASARunAggregatedResult] = []
        settings = cast(SASASettings, ctx.settings)
        settings_token = self._settings_cache_token(settings)
        if len(ctx.replicates) == 1:
            LOGGER.warning(
                "Only one replicate available for SASA aggregation in condition '%s'; "
                "replicate-level SEM is reported as 0.0",
                ctx.condition.label,
            )

        self._validate_replicate_run_coverage(ctx, results, settings)
        for run in settings.runs:
            entries = [
                next(entry for entry in result.run_results if entry.run_label == run.label)
                for result in results
            ]
            self._validate_structural_metadata_consistency(ctx, run.label, entries)

            per_means = np.asarray([entry.mean_sasa for entry in entries], dtype=np.float64)
            per_stds = [float(entry.std_sasa) for entry in entries]
            per_medians = [float(entry.median_sasa) for entry in entries]
            per_mins = [float(entry.min_sasa) for entry in entries]
            per_maxs = [float(entry.max_sasa) for entry in entries]
            per_finals = [float(entry.final_sasa) for entry in entries]

            finite_means = per_means[np.isfinite(per_means)]
            if finite_means.size:
                mean_stats = compute_sem(finite_means)
                overall_mean = mean_stats.mean
                overall_sem = mean_stats.sem
            else:
                overall_mean = float("nan")
                overall_sem = float("nan")

            overall_median = float(np.nanmean(np.asarray(per_medians, dtype=np.float64)))
            overall_min = float(np.nanmin(np.asarray(per_mins, dtype=np.float64)))
            overall_max = float(np.nanmax(np.asarray(per_maxs, dtype=np.float64)))
            overall_final = float(np.nanmean(np.asarray(per_finals, dtype=np.float64)))

            (
                residue_keys,
                residue_chainids,
                residue_resids,
                residue_resnames,
                residue_matrix,
            ) = self._load_per_residue_contributions(ctx, run.label, entries)

            if residue_matrix:
                stacked = np.stack(residue_matrix, axis=0)
                per_residue_mean = np.nanmean(stacked, axis=0).astype(np.float64)
                if stacked.shape[0] > 1:
                    per_residue_sem = np.nanstd(stacked, axis=0, ddof=1) / np.sqrt(
                        float(stacked.shape[0])
                    )
                else:
                    per_residue_sem = np.zeros(stacked.shape[1], dtype=np.float64)
            else:
                per_residue_mean = np.asarray([], dtype=np.float64)
                per_residue_sem = np.asarray([], dtype=np.float64)

            template = entries[0]
            context_atom_counts = [int(entry.n_context_atoms) for entry in entries]
            context_count_is_variable = len(set(context_atom_counts)) > 1
            aggregated_runs.append(
                SASARunAggregatedResult(
                    config_hash=first.config_hash,
                    polyzymd_version=get_polyzymd_version(),
                    replicate=None,
                    equilibration_time=first.equilibration_time,
                    equilibration_unit=first.equilibration_unit,
                    selection_string=template.target_selection,
                    replicates=list(ctx.replicates),
                    n_replicates=len(ctx.replicates),
                    run_label=run.label,
                    target_selection=template.target_selection,
                    context_selection=template.context_selection,
                    overall_mean=overall_mean,
                    overall_sem=overall_sem,
                    overall_median=overall_median,
                    overall_min=overall_min,
                    overall_max=overall_max,
                    overall_final=overall_final,
                    per_replicate_means=[float(v) for v in per_means.tolist()],
                    per_replicate_stds=per_stds,
                    per_replicate_medians=per_medians,
                    per_replicate_mins=per_mins,
                    per_replicate_maxs=per_maxs,
                    per_replicate_finals=per_finals,
                    n_target_atoms=template.n_target_atoms,
                    n_context_atoms=(None if context_count_is_variable else context_atom_counts[0]),
                    per_replicate_context_atom_counts=context_atom_counts,
                    n_context_atoms_variable=context_count_is_variable,
                    n_target_residues=template.n_target_residues,
                    zero_atom_selection=all(entry.zero_atom_selection for entry in entries),
                    residue_keys=residue_keys,
                    residue_chainids=residue_chainids,
                    residue_resids=residue_resids,
                    residue_resnames=residue_resnames,
                    per_residue_mean_sasa=per_residue_mean.tolist(),
                    per_residue_sem_sasa=per_residue_sem.tolist(),
                )
            )

        agg_result = SASAAggregatedResult(
            config_hash=first.config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=first.equilibration_time,
            equilibration_unit=first.equilibration_unit,
            selection_string=first.selection_string,
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            run_results=aggregated_runs,
            settings_fingerprint=settings_token,
            source_result_files=[],
        )

        target_path = ctx.result_path
        if target_path is None:
            target_path = ctx.output_dir / self._make_aggregated_filename(
                ctx.replicates,
                first,
                settings_token,
            )
        self.save_result(agg_result, target_path)
        return agg_result

    @staticmethod
    def _load_per_residue_contributions(
        ctx: AggregateContext,
        run_label: str,
        entries: Sequence[Any],
    ) -> tuple[list[str], list[str], list[int], list[str], list[Any]]:
        """Load per-residue SASA contributions with strict sidecar validation.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        run_label : str
            SASA run label being aggregated.
        entries : Sequence[Any]
            Replicate-level run results for the current run.

        Returns
        -------
        tuple[list[str], list[str], list[int], list[str], list[Any]]
            Residue metadata and per-replicate mean residue SASA arrays.

        Raises
        ------
        ValueError
            Raised when an expected residue-level sidecar is missing or when
            residue metadata differs between contributing replicates.
        """

        import numpy as np

        from polyzymd.analyses.shared.sasa import load_sasa_artifacts

        residue_keys: list[str] = []
        residue_chainids: list[str] = []
        residue_resids: list[int] = []
        residue_resnames: list[str] = []
        residue_matrix: list[Any] = []

        for entry in entries:
            if entry.zero_atom_selection or entry.n_target_residues == 0:
                continue

            if entry.raw_npz_path is None or entry.raw_metadata_path is None:
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"requires residue-level sidecars for replicate {entry.replicate}, "
                    "but the run result does not record them."
                )

            npz_file = Path(entry.raw_npz_path)
            metadata_file = Path(entry.raw_metadata_path)
            if not npz_file.exists() or not metadata_file.exists():
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"requires residue-level sidecars for replicate {entry.replicate}, "
                    f"but at least one sidecar is missing: npz={npz_file}, metadata={metadata_file}."
                )

            raw, _metadata = load_sasa_artifacts(npz_file, metadata_file)
            if raw.residue_sasa_a2.size == 0:
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"requires residue-level SASA data for replicate {entry.replicate}, "
                    "but the stored sidecar payload is empty."
                )

            current_keys = list(raw.residue_keys)
            current_chainids = list(raw.residue_chainids)
            current_resids = list(raw.residue_resids)
            current_resnames = list(raw.residue_resnames)
            residue_mean = np.nanmean(raw.residue_sasa_a2, axis=0).astype(np.float64)

            if residue_mean.shape[0] != len(current_keys):
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"found inconsistent residue SASA array width for replicate {entry.replicate}."
                )

            if not residue_keys:
                residue_keys = current_keys
                residue_chainids = current_chainids
                residue_resids = current_resids
                residue_resnames = current_resnames
            elif (
                current_keys != residue_keys
                or current_chainids != residue_chainids
                or current_resids != residue_resids
                or current_resnames != residue_resnames
            ):
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"found residue metadata mismatch in replicate {entry.replicate}."
                )

            residue_matrix.append(residue_mean)

        return residue_keys, residue_chainids, residue_resids, residue_resnames, residue_matrix

    @staticmethod
    def _validate_replicate_run_coverage(
        ctx: AggregateContext,
        results: Sequence[Any],
        settings: SASASettings,
    ) -> None:
        """Require every configured SASA run in every replicate result.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate SASA results.
        settings : SASASettings
            Current SASA settings with configured run labels.

        Raises
        ------
        ValueError
            Raised when any replicate result is missing a configured run or
            contains duplicate configured run labels.
        """

        expected_labels = [run.label for run in settings.runs]
        expected_set = set(expected_labels)
        coverage_issues: list[str] = []

        for result in results:
            run_labels = [entry.run_label for entry in result.run_results]
            missing = [label for label in expected_labels if label not in run_labels]
            duplicates = sorted(
                {
                    label
                    for label in run_labels
                    if label in expected_set and run_labels.count(label) > 1
                }
            )
            if not missing and not duplicates:
                continue

            issue_parts: list[str] = []
            if missing:
                issue_parts.append(f"missing runs {missing}")
            if duplicates:
                issue_parts.append(f"duplicate runs {duplicates}")
            coverage_issues.append(
                f"replicate {getattr(result, 'replicate', None)}: {', '.join(issue_parts)}"
            )

        if coverage_issues:
            issue_text = "; ".join(coverage_issues)
            raise ValueError(
                f"SASA aggregation for condition '{ctx.condition.label}' requires complete "
                f"configured run coverage for every replicate. Problems detected: {issue_text}."
            )

    @staticmethod
    def _validate_structural_metadata_consistency(
        ctx: AggregateContext,
        run_label: str,
        entries: Sequence[Any],
    ) -> None:
        """Require target SASA counts to match across contributing replicates.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        run_label : str
            SASA run label being aggregated.
        entries : Sequence[Any]
            Replicate-level run results contributing to the aggregated run.

        Raises
        ------
        ValueError
            Raised when target atom or target residue counts do not match
            across contributing replicates.
        """

        first = entries[0]
        expected_counts = (
            int(first.n_target_atoms),
            int(first.n_target_residues),
        )
        mismatch_details: list[str] = []

        for entry in entries[1:]:
            current_counts = (
                int(entry.n_target_atoms),
                int(entry.n_target_residues),
            )
            if current_counts == expected_counts:
                continue

            mismatch_details.append(
                "replicate "
                f"{entry.replicate}: target counts {current_counts} != {expected_counts} "
                f"(replicate {first.replicate})"
            )

        if mismatch_details:
            issue_text = "; ".join(mismatch_details)
            raise ValueError(
                f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                f"found target metadata mismatch across replicates. Problems detected: "
                f"{issue_text}."
            )

    @classmethod
    def _coerce_and_validate_aggregated_result(
        cls,
        result: Any,
        settings: SASASettings,
        condition: Any,
        *,
        source: Path | None = None,
    ) -> SASAAggregatedResult:
        """Coerce an aggregate and validate it against the current config.

        Parameters
        ----------
        result : Any
            Aggregated result loaded from disk or supplied in memory.
        settings : SASASettings
            Current SASA settings from the comparison configuration.
        condition : Any
            Condition whose replicate coverage must match the aggregate.
        source : Path | None, optional
            Source path for diagnostics.

        Returns
        -------
        SASAAggregatedResult
            Validated aggregate result.

        Raises
        ------
        ValueError
            Raised when cache identity, replicate coverage, or run selections
            do not match the current configuration.
        """

        if isinstance(result, dict):
            result = SASAAggregatedResult.model_validate(result)
        if not isinstance(result, SASAAggregatedResult):
            raise TypeError(
                "SASA aggregated result loader expected SASAAggregatedResult, "
                f"got {type(result).__name__}"
            )

        cls._validate_aggregated_result_identity(
            result,
            settings,
            condition,
            source=source,
        )
        return result

    @classmethod
    def _validate_aggregated_result_identity(
        cls,
        result: SASAAggregatedResult,
        settings: SASASettings,
        condition: Any,
        *,
        source: Path | None = None,
    ) -> None:
        """Validate an SASA aggregate cache against active comparison inputs.

        Parameters
        ----------
        result : SASAAggregatedResult
            Aggregated SASA result to validate.
        settings : SASASettings
            Current SASA settings.
        condition : Any
            Current comparison condition.
        source : Path | None, optional
            Source path for diagnostics.

        Raises
        ------
        ValueError
            Raised when the aggregate was produced for different replicates,
            settings, or run selections.
        """

        source_text = f" at {source}" if source is not None else ""
        condition_text = f" for condition '{condition.label}'"
        expected_replicates = sorted(condition.replicates)
        observed_replicates = sorted(result.replicates)
        if (
            result.n_replicates != len(expected_replicates)
            or observed_replicates != expected_replicates
        ):
            raise ValueError(
                "Aggregated SASA result"
                f"{condition_text} has stale replicate coverage{source_text}. Expected "
                f"replicates {expected_replicates}, found {observed_replicates} with "
                f"n_replicates={result.n_replicates}. Recompute the condition or clear "
                "stale caches before comparing."
            )

        stored_fingerprint = result.settings_fingerprint
        current_fingerprint = cls._settings_cache_token(settings)
        if stored_fingerprint is None:
            raise ValueError(
                "Aggregated SASA result"
                f"{condition_text} is missing a settings fingerprint{source_text}. Legacy "
                "SASA aggregated caches are not compatible with settings-sensitive "
                "comparison loading. Recompute the condition before comparing."
            )
        if stored_fingerprint != current_fingerprint:
            raise ValueError(
                "Aggregated SASA result"
                f"{condition_text} was computed with settings fingerprint "
                f"{stored_fingerprint}, but current settings require {current_fingerprint}"
                f"{source_text}. Recompute the condition or clear stale caches before comparing."
            )

        expected_runs = {run.label: run for run in settings.runs}
        observed_labels = [run_result.run_label for run_result in result.run_results]
        missing_runs = [label for label in expected_runs if label not in observed_labels]
        unexpected_runs = sorted(set(observed_labels) - set(expected_runs))
        duplicate_runs = sorted(
            {label for label in observed_labels if observed_labels.count(label) > 1}
        )
        if missing_runs or unexpected_runs or duplicate_runs:
            details: list[str] = []
            if missing_runs:
                details.append(f"missing runs {missing_runs}")
            if unexpected_runs:
                details.append(f"unexpected runs {unexpected_runs}")
            if duplicate_runs:
                details.append(f"duplicate runs {duplicate_runs}")
            raise ValueError(
                "Aggregated SASA result"
                f"{condition_text} has stale run coverage{source_text}: "
                f"{'; '.join(details)}. Recompute the condition before comparing."
            )

        for run_result in result.run_results:
            expected_run = expected_runs[run_result.run_label]
            if (
                run_result.target_selection != expected_run.target_selection
                or run_result.context_selection != expected_run.context_selection
            ):
                raise ValueError(
                    "Aggregated SASA run "
                    f"'{run_result.run_label}'{condition_text} has stale selection context"
                    f"{source_text}. Expected target/context "
                    f"{expected_run.target_selection!r}/{expected_run.context_selection!r}, "
                    f"found {run_result.target_selection!r}/{run_result.context_selection!r}. "
                    "Recompute the condition before comparing."
                )

            run_replicates = sorted(run_result.replicates)
            value_count = len(run_result.per_replicate_means)
            if (
                run_result.n_replicates != len(expected_replicates)
                or run_replicates != expected_replicates
                or value_count != len(expected_replicates)
            ):
                raise ValueError(
                    "Aggregated SASA run "
                    f"'{run_result.run_label}'{condition_text} has stale replicate metadata"
                    f"{source_text}. Expected replicates {expected_replicates}, found "
                    f"{run_replicates} with n_replicates={run_result.n_replicates} and "
                    f"{value_count} per-replicate values. Recompute the condition before "
                    "comparing."
                )

    def _resolve_aggregated_result_path(self, aggregated_dir: Path) -> Path | None:
        """Resolve the aggregated SASA result path.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated SASA result files.

        Returns
        -------
        Path | None
            Path to the selected JSON result, or ``None`` when absent.
        """

        if not aggregated_dir.exists():
            return None
        canonical = self.aggregate_result_path(aggregated_dir)
        if canonical.exists():
            return canonical

        json_files = sorted(aggregated_dir.glob("*.json"), key=lambda path: path.stat().st_mtime)
        if not json_files:
            return None

        chosen = json_files[-1]
        LOGGER.warning(
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
        settings: SASASettings | None = None,
        condition: Any | None = None,
    ) -> Any:
        """Load and optionally validate an aggregated SASA result.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated SASA result files.
        settings : SASASettings | None, optional
            Current SASA settings for cache identity validation.
        condition : Any | None, optional
            Current comparison condition for replicate coverage validation.

        Returns
        -------
        Any
            Loaded aggregate, or ``None`` when no result file exists.
        """

        result_path = self._resolve_aggregated_result_path(aggregated_dir)
        if result_path is None:
            return None

        result = self._deserialize_result(result_path)
        if settings is None or condition is None:
            return result

        return self._coerce_and_validate_aggregated_result(
            result,
            settings,
            condition,
            source=result_path,
        )

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare SASA runs across conditions."""
        from polyzymd import __version__
        from polyzymd.analyses.sasa._comparison_results import (
            SASAComparisonResult,
            SASAConditionSummary,
            SASARunANOVA,
            SASARunPairwiseComparison,
            SASARunSummary,
        )
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            one_way_anova,
            percent_change,
        )

        settings = cast(SASASettings, ctx.settings)
        run_labels_configured = [run.label for run in settings.runs]
        summaries: list[SASAConditionSummary] = []

        for condition in ctx.conditions:
            agg_result = ctx.aggregated_results.get(condition.label)
            if agg_result is not None:
                agg_result = self._coerce_and_validate_aggregated_result(
                    agg_result,
                    settings,
                    condition,
                    source=self._resolve_aggregated_result_path(
                        ctx.analysis_dirs[condition.label] / "aggregated"
                    ),
                )
            else:
                agg_result = self._load_aggregated_result(
                    ctx.analysis_dirs[condition.label] / "aggregated",
                    settings=settings,
                    condition=condition,
                )
            if agg_result is None:
                continue

            run_summaries = [
                SASARunSummary(
                    label=run_result.run_label,
                    target_selection=run_result.target_selection,
                    context_selection=run_result.context_selection,
                    mean_sasa=run_result.overall_mean,
                    sem_sasa=run_result.overall_sem,
                    per_replicate_means=run_result.per_replicate_means,
                    zero_atom_selection=run_result.zero_atom_selection,
                )
                for run_result in agg_result.run_results
            ]

            summaries.append(
                SASAConditionSummary(
                    label=condition.label,
                    config_path=str(condition.config_path),
                    n_replicates=agg_result.n_replicates,
                    run_summaries=run_summaries,
                )
            )

        if not summaries:
            LOGGER.warning("SASA comparison skipped because no conditions have results")
            return None

        summaries_by_label = {summary.label: summary for summary in summaries}

        run_labels: list[str] = []
        ranking_by_run: dict[str, list[str]] = {}
        pairwise: list[SASARunPairwiseComparison] = []
        anova_by_run: list[SASARunANOVA] | None = []
        effective_control = ctx.effective_control

        for run_label in run_labels_configured:
            available = filter_summaries_with_run(
                summaries_by_label,
                run_label,
                lambda summary, label: summary.get_run(label),
                logger=LOGGER,
            )
            if not available:
                continue

            comparable = []
            for summary in available.values():
                run_summary = summary.get_run(run_label)
                finite_values = [
                    value for value in run_summary.per_replicate_means if self._is_finite(value)
                ]
                if len(finite_values) >= 2 and self._is_finite(run_summary.mean_sasa):
                    comparable.append(summary)

            if not comparable:
                LOGGER.warning(
                    "Run '%s' has no finite replicate data across conditions; skipping",
                    run_label,
                )
                continue

            run_labels.append(run_label)
            ranked = sorted(
                comparable,
                key=lambda summary: (
                    not self._is_finite(summary.get_run(run_label).mean_sasa),
                    summary.get_run(run_label).mean_sasa,
                ),
            )
            ranking_by_run[run_label] = [summary.label for summary in ranked]

            if len(comparable) >= 2:
                comparable_by_label = {summary.label: summary for summary in comparable}
                condition_pairs = build_condition_pairs(
                    list(comparable_by_label.keys()),
                    effective_control,
                    on_control_missing="all_pairs",
                    logger=LOGGER,
                )
                for condition_a, condition_b in condition_pairs:
                    candidate = self._compare_run(
                        run_label=run_label,
                        condition_a=condition_a,
                        condition_b=condition_b,
                        run_a=comparable_by_label[condition_a].get_run(run_label),
                        run_b=comparable_by_label[condition_b].get_run(run_label),
                        independent_ttest=independent_ttest,
                        cohens_d=cohens_d,
                        percent_change=percent_change,
                    )
                    if candidate is not None:
                        pairwise.append(candidate)

            if len(comparable) >= 3:
                groups = []
                for summary in comparable:
                    values = [
                        v
                        for v in summary.get_run(run_label).per_replicate_means
                        if self._is_finite(v)
                    ]
                    if values:
                        groups.append(values)
                if len(groups) >= 3:
                    anova = one_way_anova(*groups)
                    testable = all(len(group) >= 2 for group in groups)
                    anova_by_run.append(
                        SASARunANOVA(
                            run_label=run_label,
                            f_statistic=anova.f_statistic if testable else None,
                            p_value=anova.p_value if testable else None,
                            significant=anova.significant if testable else False,
                            testable=testable,
                            note=None if testable else NOT_TESTABLE_SINGLETON_NOTE,
                        )
                    )

        if not run_labels:
            LOGGER.warning("SASA comparison skipped because no runs had comparable data")
            return None

        if not anova_by_run:
            anova_by_run = None

        fdr_alpha = getattr(ctx, "fdr_alpha", 0.05)
        self._apply_fdr_correction(pairwise, anova_by_run, fdr_alpha)

        return SASAComparisonResult(
            metric="mean_sasa",
            name=ctx.name,
            n_runs=len(run_labels),
            run_labels=run_labels,
            control_label=effective_control,
            fdr_alpha=fdr_alpha,
            conditions=summaries,
            pairwise_comparisons=pairwise,
            anova_by_run=anova_by_run,
            ranking_by_run=ranking_by_run,
            equilibration_time=ctx.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate SASA comparison plots via delegated plotters."""
        comparison_path = ctx.comparison_path or self.comparison_result_path(ctx.results_dir)
        if not comparison_path.exists():
            LOGGER.warning(
                "SASA comparison result not found at %s; skipping plots", comparison_path
            )
            return []

        comparison_result = self._deserialize_comparison(comparison_path)
        if comparison_result is None:
            return []

        from polyzymd.analyses.sasa._plotters import (
            plot_sasa_comparison_bars,
            plot_sasa_normalized_control_bars,
            plot_sasa_residue_profiles,
            plot_sasa_timeseries,
        )

        paths: list[Path] = []
        paths.extend(plot_sasa_comparison_bars(ctx, comparison_result))
        paths.extend(plot_sasa_timeseries(ctx, comparison_result))
        paths.extend(plot_sasa_residue_profiles(ctx, comparison_result))
        paths.extend(plot_sasa_normalized_control_bars(ctx, comparison_result))
        return paths

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format SASA comparison output via delegated formatter."""
        from polyzymd.analyses.sasa._formatters import format_sasa_comparison

        return format_sasa_comparison(result, output_format)

    @staticmethod
    def _run_cache_token(
        *,
        label: str,
        target_selection: str,
        context_selection: str,
        probe_radius_nm: float,
        n_sphere_points: int,
        stride: int,
        equilibration: str,
    ) -> str:
        """Build a stable cache token for raw SASA artifacts."""
        payload = (
            f"{label}|{target_selection}|{context_selection}|{probe_radius_nm:.6f}|"
            f"{n_sphere_points}|{stride}|{equilibration}"
        )
        digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]
        safe_label = label.replace(" ", "_").replace("-", "_").replace("/", "_").lower()
        return f"{safe_label}_{digest}"

    @staticmethod
    def _settings_cache_token(settings: SASASettings) -> str:
        """Build a stable token for SASA settings-based cache invalidation."""
        payload = json.dumps(settings.model_dump(mode="json"), sort_keys=True)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]

    @staticmethod
    def _compare_run(
        *,
        run_label: str,
        condition_a: str,
        condition_b: str,
        run_a: Any,
        run_b: Any,
        independent_ttest: Any,
        cohens_d: Any,
        percent_change: Any,
    ) -> Any:
        """Compare one SASA run between two conditions."""
        from polyzymd.analyses.sasa._comparison_results import SASARunPairwiseComparison
        from polyzymd.analyses.stats import interpret_direction

        values_a = [value for value in run_a.per_replicate_means if SASAAnalysis._is_finite(value)]
        values_b = [value for value in run_b.per_replicate_means if SASAAnalysis._is_finite(value)]
        testable = len(values_a) >= 2 and len(values_b) >= 2

        t_result = independent_ttest(values_a, values_b)
        d_result = cohens_d(values_a, values_b)
        pct_change = percent_change(run_a.mean_sasa, run_b.mean_sasa)
        direction = interpret_direction(
            pct_change,
            direction_labels=("shielding", "unchanged", "exposure"),
            threshold=1.0,
        )

        return SASARunPairwiseComparison(
            run_label=run_label,
            condition_a=condition_a,
            condition_b=condition_b,
            t_statistic=t_result.t_statistic if testable else None,
            p_value=t_result.p_value if testable else None,
            cohens_d=d_result.cohens_d if testable else None,
            effect_interpretation=d_result.interpretation,
            direction=direction,
            significant=t_result.significant if testable else False,
            percent_change=pct_change,
            testable=testable,
            note=None if testable else NOT_TESTABLE_SINGLETON_NOTE,
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
        apply_fdr_correction(pairwise, anova_by_run, fdr_alpha)

    @staticmethod
    def _has_run_summary(summary: Any, run_label: str) -> bool:
        """Check whether a condition summary includes a given run."""
        try:
            summary.get_run(run_label)
            return True
        except KeyError:
            return False

    @staticmethod
    def _is_finite(value: float) -> bool:
        """Return True when value is finite."""
        return value == value and value not in (float("inf"), float("-inf"))

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
        settings_token: str,
    ) -> str:
        """Generate an aggregated SASA filename."""
        eq_str = f"eq{first_result.equilibration_time:g}{first_result.equilibration_unit}"
        rep_str = Analysis._format_replicate_range(replicates)
        return f"sasa_{rep_str}_{eq_str}_{settings_token}.json"

    @staticmethod
    def _deserialize_comparison(path: Path) -> SASAComparisonResult | None:
        """Load SASA comparison result from disk."""
        from polyzymd.analyses.sasa._comparison_results import SASAComparisonResult

        if not path.exists():
            return None
        return SASAComparisonResult.load(path)
