"""RMSD analysis plugin.

Computes per-frame Root Mean Square Deviation for one or more named selections,
aggregates across replicates, performs per-run cross-condition comparisons,
and exposes plotting/formatting hooks.
"""

from __future__ import annotations

import hashlib
import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator, model_validator

import polyzymd.analyses.rmsd._plot_settings as _plot_settings  # noqa: F401
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.rmsd._results import RMSDAggregatedResult
from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory
from polyzymd.analyses.shared.config_hash import compute_config_hash
from polyzymd.analyses.shared.loader import (
    TrajectoryLoader,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from polyzymd.analyses.rmsd._comparison_results import RMSDComparisonResult
    from polyzymd.analyses.rmsd._results import RMSDResult

logger = logging.getLogger(__name__)


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
    AggregatedResultClass: ClassVar[type] = RMSDAggregatedResult
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel, equilibration: str) -> str:
        """Build a short cache tag for settings and equilibration.

        Parameters
        ----------
        settings : BaseModel
            Analysis settings model.
        equilibration : str
            Equilibration string from analysis context.

        Returns
        -------
        str
            First 8 hex characters of the MD5 fingerprint.
        """
        payload = f"{settings.model_dump_json()}|{equilibration}".encode("utf-8")
        return hashlib.md5(payload).hexdigest()[:8]

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Compute RMSD for all configured runs for a single replicate.

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
        import numpy as np

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rmsd._results import RMSDResult

        settings = ctx.settings
        sim_config = ctx.sim_config

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        eq_str = f"eq{eq_value:.2f}{eq_unit}"
        settings_tag = self._make_settings_cache_tag(settings, ctx.equilibration)
        result_file = ctx.output_dir / f"rmsd_{eq_str}_{settings_tag}.json"

        cached = self._check_cache(
            RMSDResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=sim_config,
        )
        if cached is not None:
            return cached

        loader = TrajectoryLoader(sim_config)
        config_hash = compute_config_hash(sim_config)

        u0 = loader.load_universe(replicate, cache=False)
        traj_info = loader.get_trajectory_info(replicate)
        timestep_ps = loader.get_timestep(replicate, unit="ps")

        n_frames_total = len(u0.trajectory)
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")
        start_frame = time_to_frame(eq_time_ps, "ps", timestep_ps, "ps")
        n_frames_used = n_frames_total - start_frame

        if n_frames_used <= 0:
            raise ValueError(
                "Equilibration removed all frames for RMSD analysis. "
                f"Got start_frame={start_frame}, n_frames_total={n_frames_total}."
            )

        run_results = []
        for run in settings.runs:
            run_result = self._compute_single_run(
                ctx=ctx,
                replicate=replicate,
                run=run,
                loader=loader,
                config_hash=config_hash,
                eq_value=eq_value,
                eq_unit=eq_unit,
                eq_str=eq_str,
                settings_tag=settings_tag,
                start_frame=start_frame,
                n_frames_total=n_frames_total,
                n_frames_used=n_frames_used,
                timestep_ps=timestep_ps,
            )
            run_results.append(run_result)

        result = RMSDResult(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string="; ".join(run.selection for run in settings.runs),
            run_results=run_results,
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
            trajectory_files=[str(path) for path in traj_info.trajectory_files],
        )

        result.save(result_file)
        logger.info("Saved RMSD result to %s", result_file)

        return result

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

        first = results[0]
        run_labels = [run.label for run in ctx.settings.runs]

        if len(ctx.replicates) == 1:
            logger.warning(
                "Only one replicate available for RMSD aggregation in condition '%s'; "
                "replicate-level SEM is reported as 0.0",
                ctx.condition.label,
            )

        aggregated_runs: list[RMSDRunAggregatedResult] = []
        for run_label in run_labels:
            run_entries = []
            for result in results:
                for run_result in result.run_results:
                    if run_result.run_label == run_label:
                        run_entries.append(run_result)
                        break

            if not run_entries:
                logger.warning(
                    "No RMSD entries found for run '%s'; skipping in aggregate", run_label
                )
                continue

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
                    replicates=list(ctx.replicates),
                    n_replicates=len(ctx.replicates),
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
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            run_results=aggregated_runs,
            source_result_files=[],
        )

        target_path = ctx.result_path
        if target_path is None:
            target_path = ctx.output_dir / self._make_aggregated_filename(ctx.replicates, first)
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
            Comparison result, or ``None`` if no conditions have data.
        """
        from polyzymd import __version__
        from polyzymd.analyses.rmsd._comparison_results import (
            RMSDComparisonResult,
            RMSDConditionSummary,
            RMSDRunANOVA,
            RMSDRunPairwiseComparison,
            RMSDRunSummary,
        )
        from polyzymd.compare.statistics import (
            cohens_d,
            independent_ttest,
            one_way_anova,
            percent_change,
        )

        run_labels = [run.label for run in ctx.settings.runs]

        summaries: list[RMSDConditionSummary] = []
        for condition in ctx.conditions:
            agg_result = ctx.aggregated_results.get(condition.label)
            if agg_result is None:
                agg_dir = ctx.analysis_dirs[condition.label] / "aggregated"
                agg_result = self._load_aggregated_result(agg_dir)

            if agg_result is None:
                logger.warning(
                    "No aggregated RMSD result for condition '%s'; skipping", condition.label
                )
                continue

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

        ranking_by_run: dict[str, list[str]] = {}
        for run_label in run_labels:
            summaries_with_run = []
            for summary in summaries:
                try:
                    run_summary = summary.get_run(run_label)
                except KeyError:
                    logger.warning(
                        "Run '%s' missing for condition '%s'; excluding from ranking",
                        run_label,
                        summary.label,
                    )
                    continue
                summaries_with_run.append((summary, run_summary))

            ranked = sorted(summaries_with_run, key=lambda entry: entry[1].mean_rmsd)
            ranking_by_run[run_label] = [summary.label for summary, _ in ranked]

        pairwise_comparisons: list[RMSDRunPairwiseComparison] = []
        if len(summaries) >= 2:
            for run_label in run_labels:
                summaries_with_run = []
                for summary in summaries:
                    try:
                        run_summary = summary.get_run(run_label)
                    except KeyError:
                        logger.warning(
                            "Run '%s' missing for condition '%s'; excluding from pairwise comparison",
                            run_label,
                            summary.label,
                        )
                        continue
                    summaries_with_run.append((summary, run_summary))

                if len(summaries_with_run) < 2:
                    logger.warning(
                        "Run '%s' has fewer than two conditions with data; skipping pairwise comparison",
                        run_label,
                    )
                    continue

                if effective_control:
                    control_entry = next(
                        (
                            (summary, run_summary)
                            for summary, run_summary in summaries_with_run
                            if summary.label == effective_control
                        ),
                        None,
                    )
                    if control_entry is None:
                        logger.warning(
                            "Run '%s' missing in control condition '%s'; skipping run comparisons",
                            run_label,
                            effective_control,
                        )
                        continue

                    control_summary, control_run = control_entry
                    for summary, run_summary in summaries_with_run:
                        if summary.label == control_summary.label:
                            continue
                        pairwise_comparisons.append(
                            self._compare_run(
                                run_label=run_label,
                                condition_a=control_summary.label,
                                condition_b=summary.label,
                                run_a=control_run,
                                run_b=run_summary,
                            )
                        )
                else:
                    for i, (summary_a, run_a) in enumerate(summaries_with_run):
                        for summary_b, run_b in summaries_with_run[i + 1 :]:
                            pairwise_comparisons.append(
                                self._compare_run(
                                    run_label=run_label,
                                    condition_a=summary_a.label,
                                    condition_b=summary_b.label,
                                    run_a=run_a,
                                    run_b=run_b,
                                )
                            )

        anova_by_run: list[RMSDRunANOVA] | None = None
        if len(summaries) >= 3:
            anova_by_run = []
            for run_label in run_labels:
                groups = []
                for summary in summaries:
                    try:
                        run_summary = summary.get_run(run_label)
                    except KeyError:
                        logger.warning(
                            "Run '%s' missing for condition '%s'; excluding from ANOVA",
                            run_label,
                            summary.label,
                        )
                        continue
                    groups.append(run_summary.per_replicate_means)

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
        """Compute one RMSD run for a single replicate."""
        import numpy as np
        from MDAnalysis.analysis.rms import RMSD

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rmsd._results import RMSDRunResult
        from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time
        from polyzymd.analyses.shared.convergence import find_convergence_time

        u = loader.load_universe(replicate, cache=False)

        centroid_selection = run.centroid_selection
        if run.reference_mode == "centroid" and centroid_selection is None:
            centroid_selection = run.alignment_selection
            logger.info(
                "Run '%s': centroid_selection not set, using alignment_selection='%s'",
                run.label,
                centroid_selection,
            )

        reference_frame_1indexed: int | None
        if run.reference_mode == "frame":
            reference_frame_1indexed = run.reference_frame + 1
        else:
            reference_frame_1indexed = None

        alignment_config = AlignmentConfig(
            enabled=True,
            reference_mode=run.reference_mode,
            reference_frame=reference_frame_1indexed,
            selection=run.alignment_selection,
            centroid_selection=centroid_selection or run.alignment_selection,
            reference_file=(Path(run.reference_file) if run.reference_file is not None else None),
        )
        ref_frame_idx = align_trajectory(
            u,
            alignment_config,
            start_frame=start_frame,
            stop_frame=n_frames_total,
        )

        # Build RMSD reference according to requested mode
        atom_group = u.select_atoms(run.selection)
        if len(atom_group) == 0:
            raise ValueError(f"Run '{run.label}' selection matched no atoms: {run.selection!r}")

        reference_universe, reference_atom_group = self._build_reference_structure(
            universe=u,
            atom_group=atom_group,
            run=run,
            start_frame=start_frame,
            stop_frame=n_frames_total,
            ref_frame_idx=ref_frame_idx,
        )

        rmsd_analysis = RMSD(
            atom_group,
            reference=reference_atom_group,
            select="all",
            ref_frame=0,
        ).run(start=start_frame, stop=n_frames_total)

        rmsd_values = rmsd_analysis.results.rmsd[:, 2].astype(np.float64)
        frames = np.arange(start_frame, n_frames_total, dtype=np.int64)
        time_ns = (frames.astype(np.float64) * timestep_ps) / 1000.0

        mean_rmsd = float(np.mean(rmsd_values))
        std_rmsd = float(np.std(rmsd_values, ddof=0))
        median_rmsd = float(np.median(rmsd_values))
        min_rmsd = float(np.min(rmsd_values))
        max_rmsd = float(np.max(rmsd_values))
        final_rmsd = float(rmsd_values[-1])

        sem_rmsd: float | None = None
        correlation_time: float | None = None
        correlation_time_unit: str | None = None
        n_independent_frames: int | None = None
        statistical_inefficiency: float | None = None
        autocorrelation_warning: str | None = None

        if len(rmsd_values) >= 20:
            tau_result = estimate_correlation_time(
                rmsd_values,
                timestep=timestep_ps,
                timestep_unit="ps",
                method="integration",
                n_frames=len(rmsd_values),
            )
            correlation_time = tau_result.tau
            correlation_time_unit = tau_result.tau_unit
            n_independent_frames = tau_result.n_independent
            statistical_inefficiency = tau_result.statistical_inefficiency
            autocorrelation_warning = tau_result.warning
            if n_independent_frames > 0:
                sem_rmsd = float(std_rmsd / np.sqrt(float(n_independent_frames)))

        convergence_result = find_convergence_time(
            time_ns,
            rmsd_values,
            window_size_ns=run.convergence_window_size_ns,
            step_size_ns=run.convergence_step_size_ns,
            slope_threshold=run.convergence_slope_threshold,
            sustained_for_ns=run.convergence_sustained_for_ns,
        )

        npz_filename = f"rmsd_{run.label}_{eq_str}_{settings_tag}_timeseries.npz"
        npz_path = ctx.output_dir / npz_filename
        np.savez_compressed(
            npz_path,
            rmsd_values=rmsd_values,
            time_ns=time_ns,
            frames=frames,
            convergence_window_start_ns=np.asarray(
                convergence_result.window_start_times_ns,
                dtype=np.float64,
            ),
            convergence_window_mean_rmsd=np.asarray(
                convergence_result.window_mean_values,
                dtype=np.float64,
            ),
            convergence_slope_time_ns=np.asarray(
                convergence_result.slope_times_ns,
                dtype=np.float64,
            ),
            convergence_slope=np.asarray(convergence_result.slopes, dtype=np.float64),
            convergence_converged=np.asarray(convergence_result.converged, dtype=np.bool_),
            convergence_time_ns=np.asarray(
                (
                    np.nan
                    if convergence_result.convergence_time_ns is None
                    else convergence_result.convergence_time_ns
                ),
                dtype=np.float64,
            ),
        )

        return RMSDRunResult(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string=run.selection,
            correlation_time=correlation_time,
            n_independent_frames=n_independent_frames,
            run_label=run.label,
            selection=run.selection,
            alignment_selection=run.alignment_selection,
            reference_mode=run.reference_mode,
            reference_frame=(ref_frame_idx + 1 if ref_frame_idx is not None else None),
            reference_file=run.reference_file,
            mean_rmsd=mean_rmsd,
            std_rmsd=std_rmsd,
            median_rmsd=median_rmsd,
            min_rmsd=min_rmsd,
            max_rmsd=max_rmsd,
            final_rmsd=final_rmsd,
            sem_rmsd=sem_rmsd,
            correlation_time_unit=correlation_time_unit,
            statistical_inefficiency=statistical_inefficiency,
            autocorrelation_warning=autocorrelation_warning,
            converged=convergence_result.converged,
            convergence_assessable=convergence_result.assessable,
            convergence_time_ns=convergence_result.convergence_time_ns,
            convergence_message=convergence_result.message,
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
            npz_path=str(npz_path),
            time_unit="ns",
            timestep_ps=timestep_ps,
        )

    def _build_reference_structure(
        self,
        *,
        universe: Any,
        atom_group: Any,
        run: RMSDRunSettings,
        start_frame: int,
        stop_frame: int,
        ref_frame_idx: int | None,
    ) -> tuple[Any, Any]:
        """Build reference universe and atom group for RMSD calculations."""
        import MDAnalysis as mda
        import numpy as np
        from MDAnalysis.coordinates.memory import MemoryReader

        if run.reference_mode in {"centroid", "frame"}:
            if ref_frame_idx is None:
                raise ValueError(
                    f"Run '{run.label}' expected a reference frame for mode '{run.reference_mode}'"
                )
            universe.trajectory[ref_frame_idx]
            ref_positions = atom_group.positions.copy().astype(np.float64)
        elif run.reference_mode == "average":
            positions = []
            for frame_idx in range(start_frame, stop_frame):
                universe.trajectory[frame_idx]
                positions.append(atom_group.positions.copy().astype(np.float64))
            ref_positions = np.mean(np.stack(positions, axis=0), axis=0)
        elif run.reference_mode == "external":
            if run.reference_file is None:
                raise ValueError(
                    f"Run '{run.label}' requires reference_file when reference_mode='external'"
                )

            ref_path = Path(run.reference_file)
            logger.info("Using external reference from: %s", ref_path)

            ref_universe = mda.Universe(str(ref_path))
            ref_atoms = ref_universe.select_atoms(run.selection)

            if len(ref_atoms) == 0:
                raise ValueError(
                    f"Run '{run.label}' external PDB '{ref_path.name}' has no atoms matching "
                    f"selection {run.selection!r}."
                )

            if len(ref_atoms) != len(atom_group):
                logger.warning(
                    "Run '%s' external reference atom count mismatch for selection %r "
                    "(trajectory=%d, external=%d)",
                    run.label,
                    run.selection,
                    len(atom_group),
                    len(ref_atoms),
                )
                raise ValueError(
                    f"Run '{run.label}' atom count mismatch between trajectory "
                    f"({len(atom_group)}) and external PDB ({len(ref_atoms)}) for "
                    f"selection {run.selection!r}."
                )

            ref_positions = ref_atoms.positions.copy().astype(np.float64)
        else:
            raise ValueError(f"Unsupported RMSD reference_mode: {run.reference_mode!r}")

        reference_universe = mda.Merge(atom_group)
        reference_universe.load_new(ref_positions[np.newaxis, :, :], format=MemoryReader)
        reference_atom_group = reference_universe.atoms
        return reference_universe, reference_atom_group

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
        from polyzymd.analyses.stats import interpret_direction
        from polyzymd.compare.statistics import cohens_d, independent_ttest, percent_change

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
        from polyzymd.compare.statistics import benjamini_hochberg

        if pairwise:
            raw_p = [comp.p_value if comp.testable else None for comp in pairwise]
            bh_results = benjamini_hochberg(raw_p, alpha=fdr_alpha)
            for comp, bh in zip(pairwise, bh_results, strict=False):
                comp.p_value_adjusted = bh.adjusted_p_value
                comp.significant = bh.significant

        if anova_by_run:
            raw_p = [a.p_value if a.testable else None for a in anova_by_run]
            bh_results = benjamini_hochberg(raw_p, alpha=fdr_alpha)
            for a, bh in zip(anova_by_run, bh_results, strict=False):
                a.significant = bh.significant

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Generate an aggregated RMSD filename."""
        eq_str = f"eq{first_result.equilibration_time:.2f}{first_result.equilibration_unit}"
        rep_str = Analysis._format_replicate_range(replicates)
        return f"rmsd_{rep_str}_{eq_str}.json"

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
