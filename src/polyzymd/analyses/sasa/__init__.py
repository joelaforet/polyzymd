"""SASA analysis plugin."""

from __future__ import annotations

import hashlib
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence, cast

from pydantic import BaseModel, Field, field_validator, model_validator

import polyzymd.analyses.sasa._plot_settings as _plot_settings  # noqa: F401
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.sasa._results import SASAAggregatedResult, SASAResult
from polyzymd.analyses.shared.config_hash import compute_config_hash
from polyzymd.analyses.shared.loader import (
    TrajectoryLoader,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analyses.shared.sasa import compute_sasa, save_sasa_artifacts
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.sasa._comparison_results import SASAComparisonResult

LOGGER = logging.getLogger(__name__)


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
    AggregatedResultClass: ClassVar[type | None] = SASAAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = SASAResult
    execution_cost_hint: ClassVar[str] = "high"
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    # SASA is a mean-based observable (all frames contribute, SEM corrected via N_eff)

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Compute SASA for all configured runs for a single replicate."""
        import numpy as np

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.sasa._results import SASAResult, SASARunResult
        from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

        settings = cast(SASASettings, ctx.settings)
        sim_config = ctx.sim_config

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        eq_str = f"eq{eq_value:.2f}{eq_unit}"
        settings_token = self._settings_cache_token(settings)
        result_file = ctx.output_dir / f"sasa_{eq_str}_{settings_token}.json"

        cached = self._check_cache(
            SASAResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=sim_config,
        )
        if cached is not None:
            return cached

        loader = TrajectoryLoader(sim_config)
        config_hash = compute_config_hash(sim_config)

        universe = loader.load_universe(replicate, cache=False)
        traj_info = loader.get_trajectory_info(replicate)
        timestep_ps = loader.get_timestep(replicate, unit="ps")

        n_frames_total = len(universe.trajectory)
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")
        start_frame = time_to_frame(eq_time_ps, "ps", timestep_ps, "ps")
        n_frames_used = n_frames_total - start_frame
        if n_frames_used <= 0:
            raise ValueError(
                "Equilibration removed all frames for SASA analysis. "
                f"Got start_frame={start_frame}, n_frames_total={n_frames_total}."
            )

        run_results: list[SASARunResult] = []
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        for run in settings.runs:
            run_token = self._run_cache_token(
                label=run.label,
                target_selection=run.target_selection,
                context_selection=run.context_selection or run.target_selection,
                probe_radius_nm=settings.probe_radius_nm,
                n_sphere_points=settings.n_sphere_points,
                equilibration=ctx.equilibration,
            )
            npz_path = ctx.output_dir / f"sasa_{run_token}.npz"
            metadata_path = ctx.output_dir / f"sasa_{run_token}.json"

            raw = compute_sasa(
                universe,
                run_label=run.label,
                target_selection=run.target_selection,
                context_selection=run.context_selection or run.target_selection,
                probe_radius_nm=settings.probe_radius_nm,
                n_sphere_points=settings.n_sphere_points,
                start_frame=start_frame,
                stop_frame=n_frames_total,
                timestep_ps=timestep_ps,
                chunk_size=settings.chunk_size,
                stride=run.stride,
            )

            save_sasa_artifacts(
                npz_path=npz_path,
                metadata_path=metadata_path,
                result=raw,
                run_label=run.label,
                target_selection=run.target_selection,
                context_selection=run.context_selection or run.target_selection,
                probe_radius_nm=settings.probe_radius_nm,
                n_sphere_points=settings.n_sphere_points,
                equilibration=ctx.equilibration,
            )

            total = raw.total_sasa_a2
            finite_total = total[np.isfinite(total)]
            zero_atom = raw.target_atom_indices.size == 0 or raw.context_atom_indices.size == 0
            if zero_atom:
                LOGGER.warning(
                    "Run '%s' selection matched zero atoms in replicate %d; "
                    "recording NaN SASA metrics",
                    run.label,
                    replicate,
                )

            if finite_total.size:
                mean_sasa = float(np.mean(finite_total))
                std_sasa = float(np.std(finite_total, ddof=0))
                median_sasa = float(np.median(finite_total))
                min_sasa = float(np.min(finite_total))
                max_sasa = float(np.max(finite_total))
                final_sasa = float(total[-1]) if np.isfinite(total[-1]) else float("nan")
            else:
                mean_sasa = float("nan")
                std_sasa = float("nan")
                median_sasa = float("nan")
                min_sasa = float("nan")
                max_sasa = float("nan")
                final_sasa = float("nan")

            sem_sasa: float | None = None
            correlation_time: float | None = None
            correlation_time_unit: str | None = None
            n_independent_frames: int | None = None
            statistical_inefficiency: float | None = None
            autocorrelation_warning: str | None = None
            if finite_total.size >= 20:
                tau = estimate_correlation_time(
                    finite_total,
                    timestep=timestep_ps,
                    timestep_unit="ps",
                    method="integration",
                    n_frames=len(finite_total),
                )
                correlation_time = tau.tau
                correlation_time_unit = tau.tau_unit
                n_independent_frames = tau.n_independent
                statistical_inefficiency = tau.statistical_inefficiency
                autocorrelation_warning = tau.warning
                if n_independent_frames > 0 and np.isfinite(std_sasa):
                    sem_sasa = float(std_sasa / np.sqrt(float(n_independent_frames)))

            run_results.append(
                SASARunResult(
                    config_hash=config_hash,
                    polyzymd_version=get_polyzymd_version(),
                    replicate=replicate,
                    equilibration_time=eq_value,
                    equilibration_unit=eq_unit,
                    selection_string=run.target_selection,
                    correlation_time=correlation_time,
                    correlation_time_unit=correlation_time_unit,
                    n_independent_frames=n_independent_frames,
                    statistical_inefficiency=statistical_inefficiency,
                    autocorrelation_warning=autocorrelation_warning,
                    run_label=run.label,
                    target_selection=run.target_selection,
                    context_selection=run.context_selection or run.target_selection,
                    mean_sasa=mean_sasa,
                    std_sasa=std_sasa,
                    median_sasa=median_sasa,
                    min_sasa=min_sasa,
                    max_sasa=max_sasa,
                    final_sasa=final_sasa,
                    sem_sasa=sem_sasa,
                    n_frames_total=n_frames_total,
                    n_frames_used=n_frames_used,
                    n_target_atoms=int(raw.target_atom_indices.size),
                    n_context_atoms=int(raw.context_atom_indices.size),
                    n_target_residues=len(raw.residue_keys),
                    zero_atom_selection=zero_atom,
                    raw_npz_path=str(npz_path),
                    raw_metadata_path=str(metadata_path),
                    npz_path=str(npz_path),
                    metadata_path=str(metadata_path),
                    time_unit="ns",
                    timestep_ps=timestep_ps,
                )
            )

        result = SASAResult(
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
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
            trajectory_files=[str(path) for path in traj_info.trajectory_files],
        )
        result.save(result_file)
        return result

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate SASA results across replicates for one condition."""
        import numpy as np

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.sasa._results import SASAAggregatedResult, SASARunAggregatedResult
        from polyzymd.analyses.shared.sasa import load_sasa_artifacts

        first = results[0]
        aggregated_runs: list[SASARunAggregatedResult] = []
        if len(ctx.replicates) == 1:
            LOGGER.warning(
                "Only one replicate available for SASA aggregation in condition '%s'; "
                "replicate-level SEM is reported as 0.0",
                ctx.condition.label,
            )

        settings = cast(SASASettings, ctx.settings)
        for run in settings.runs:
            entries = []
            for result in results:
                match = next(
                    (entry for entry in result.run_results if entry.run_label == run.label), None
                )
                if match is not None:
                    entries.append(match)
            if not entries:
                LOGGER.warning("No SASA entries found for run '%s'; skipping", run.label)
                continue

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

            residue_matrix: list[np.ndarray] = []
            residue_keys: list[str] = []
            residue_chainids: list[str] = []
            residue_resids: list[int] = []
            residue_resnames: list[str] = []
            for entry in entries:
                if (
                    entry.zero_atom_selection
                    or entry.raw_npz_path is None
                    or entry.raw_metadata_path is None
                ):
                    continue
                npz_file = Path(entry.raw_npz_path)
                metadata_file = Path(entry.raw_metadata_path)
                if not npz_file.exists() or not metadata_file.exists():
                    continue
                raw, _metadata = load_sasa_artifacts(npz_file, metadata_file)
                if raw.residue_sasa_a2.size == 0:
                    continue
                if not residue_keys:
                    residue_keys = raw.residue_keys
                    residue_chainids = raw.residue_chainids
                    residue_resids = raw.residue_resids
                    residue_resnames = raw.residue_resnames
                residue_matrix.append(np.nanmean(raw.residue_sasa_a2, axis=0))

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
                    n_context_atoms=template.n_context_atoms,
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
            source_result_files=[],
        )

        target_path = ctx.result_path
        if target_path is None:
            target_path = ctx.output_dir / self._make_aggregated_filename(ctx.replicates, first)
        self.save_result(agg_result, target_path)
        return agg_result

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
        from polyzymd.compare.statistics import (
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
            if agg_result is None:
                agg_result = self._load_aggregated_result(
                    ctx.analysis_dirs[condition.label] / "aggregated"
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

        if len(summaries) < 2:
            LOGGER.warning("SASA comparison skipped because fewer than 2 conditions have results")
            return None

        run_labels: list[str] = []
        ranking_by_run: dict[str, list[str]] = {}
        pairwise: list[SASARunPairwiseComparison] = []
        anova_by_run: list[SASARunANOVA] | None = []
        effective_control = ctx.effective_control

        for run_label in run_labels_configured:
            available = [
                summary for summary in summaries if self._has_run_summary(summary, run_label)
            ]
            if len(available) < 2:
                continue

            comparable = []
            for summary in available:
                run_summary = summary.get_run(run_label)
                finite_values = [
                    value for value in run_summary.per_replicate_means if self._is_finite(value)
                ]
                if len(finite_values) >= 2 and self._is_finite(run_summary.mean_sasa):
                    comparable.append(summary)

            if len(comparable) < 2:
                LOGGER.warning(
                    "Run '%s' has insufficient finite replicate data across conditions; skipping",
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

            if effective_control:
                control_summary = next(
                    (summary for summary in comparable if summary.label == effective_control),
                    None,
                )
                if control_summary is not None:
                    control_run = control_summary.get_run(run_label)
                    for summary in comparable:
                        if summary.label == effective_control:
                            continue
                        candidate = self._compare_run(
                            run_label=run_label,
                            condition_a=control_summary.label,
                            condition_b=summary.label,
                            run_a=control_run,
                            run_b=summary.get_run(run_label),
                            independent_ttest=independent_ttest,
                            cohens_d=cohens_d,
                            percent_change=percent_change,
                        )
                        if candidate is not None:
                            pairwise.append(candidate)
                else:
                    for i, summary_a in enumerate(comparable):
                        run_a = summary_a.get_run(run_label)
                        for summary_b in comparable[i + 1 :]:
                            candidate = self._compare_run(
                                run_label=run_label,
                                condition_a=summary_a.label,
                                condition_b=summary_b.label,
                                run_a=run_a,
                                run_b=summary_b.get_run(run_label),
                                independent_ttest=independent_ttest,
                                cohens_d=cohens_d,
                                percent_change=percent_change,
                            )
                            if candidate is not None:
                                pairwise.append(candidate)
            else:
                for i, summary_a in enumerate(comparable):
                    run_a = summary_a.get_run(run_label)
                    for summary_b in comparable[i + 1 :]:
                        candidate = self._compare_run(
                            run_label=run_label,
                            condition_a=summary_a.label,
                            condition_b=summary_b.label,
                            run_a=run_a,
                            run_b=summary_b.get_run(run_label),
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
                    anova_by_run.append(
                        SASARunANOVA(
                            run_label=run_label,
                            f_statistic=anova.f_statistic,
                            p_value=anova.p_value,
                            significant=anova.significant,
                        )
                    )

        if not run_labels:
            LOGGER.warning("SASA comparison skipped because no runs had comparable data")
            return None

        if not anova_by_run:
            anova_by_run = None

        return SASAComparisonResult(
            metric="mean_sasa",
            name=ctx.name,
            n_runs=len(run_labels),
            run_labels=run_labels,
            control_label=effective_control,
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
            plot_sasa_residue_profiles,
            plot_sasa_timeseries,
        )

        paths: list[Path] = []
        paths.extend(plot_sasa_comparison_bars(ctx, comparison_result))
        paths.extend(plot_sasa_timeseries(ctx, comparison_result))
        paths.extend(plot_sasa_residue_profiles(ctx, comparison_result))
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
        equilibration: str,
    ) -> str:
        """Build a stable cache token for raw SASA artifacts."""
        payload = (
            f"{label}|{target_selection}|{context_selection}|{probe_radius_nm:.6f}|"
            f"{n_sphere_points}|{equilibration}"
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

        values_a = [value for value in run_a.per_replicate_means if SASAAnalysis._is_finite(value)]
        values_b = [value for value in run_b.per_replicate_means if SASAAnalysis._is_finite(value)]
        if len(values_a) < 2 or len(values_b) < 2:
            return None

        t_result = independent_ttest(values_a, values_b)
        d_result = cohens_d(values_a, values_b)
        pct_change = percent_change(run_a.mean_sasa, run_b.mean_sasa)
        if pct_change < -1.0:
            direction = "shielding"
        elif pct_change > 1.0:
            direction = "exposure"
        else:
            direction = "unchanged"

        return SASARunPairwiseComparison(
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
    ) -> str:
        """Generate an aggregated SASA filename."""
        eq_str = f"eq{first_result.equilibration_time:.2f}{first_result.equilibration_unit}"
        rep_str = Analysis._format_replicate_range(replicates)
        return f"sasa_{rep_str}_{eq_str}.json"

    @staticmethod
    def _deserialize_comparison(path: Path) -> SASAComparisonResult | None:
        """Load SASA comparison result from disk."""
        from polyzymd.analyses.sasa._comparison_results import SASAComparisonResult

        if not path.exists():
            return None
        return SASAComparisonResult.load(path)
