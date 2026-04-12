"""Rg analysis plugin.

Computes per-frame Radius of Gyration for one or more named selections,
aggregates across replicates, performs per-run cross-condition comparisons,
and exposes plotting/formatting hooks.
"""

from __future__ import annotations

import logging
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
from polyzymd.analyses.rg._results import RgAggregatedResult
from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
from polyzymd.analyses.shared.loader import (
    TrajectoryLoader,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analyses.shared.multi_run_comparison import (
    apply_fdr_correction,
    build_condition_pairs,
    filter_summaries_with_run,
)
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.rg._comparison_results import RgComparisonResult
    from polyzymd.analyses.rg._results import RgResult, RgRunResult

logger = logging.getLogger(__name__)


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

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Compute Rg for all configured runs for a single replicate.

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
        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rg._results import RgResult

        settings = ctx.settings
        sim_config = ctx.sim_config

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        eq_str = f"eq{eq_value:g}{eq_unit}"
        settings_tag = self._make_settings_cache_tag(settings)
        result_file = ctx.output_dir / f"rg_{eq_str}_{settings_tag}.json"

        cached = self._check_cache(
            RgResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=sim_config,
            settings=ctx.settings,
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
                "Equilibration removed all frames for Rg analysis. "
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
            if run_result is not None:
                run_results.append(run_result)

        result = RgResult(
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
        logger.info("Saved Rg result to %s", result_file)

        return result

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

        first = results[0]
        run_labels = [run.label for run in ctx.settings.runs]

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
            for result in results:
                for run_result in result.run_results:
                    if run_result.run_label == run_label:
                        run_entries.append(run_result)
                        break

            if not run_entries:
                logger.warning(
                    "No Rg entries found for run '%s' in condition '%s' — selection may "
                    "match no atoms in this condition; skipping in aggregate",
                    run_label,
                    ctx.condition.label,
                )
                continue

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

            # --- Load NPZ sidecars for histogram aggregation (both modes) ---
            all_reduced_rg_per_rep: list[np.ndarray] = []
            all_fragment_rg_per_rep: list[np.ndarray] = []
            for entry in run_entries:
                if entry.npz_path is None:
                    continue
                npz_path = Path(entry.npz_path)
                if not npz_path.exists():
                    continue

                with np.load(npz_path) as npz_data:
                    if "rg_values" in npz_data:
                        reduced_values = np.asarray(npz_data["rg_values"], dtype=np.float64)
                        if reduced_values.size > 0:
                            all_reduced_rg_per_rep.append(reduced_values)

                    if "fragment_rg_values" in npz_data:
                        fragment_values = np.asarray(
                            npz_data["fragment_rg_values"], dtype=np.float64
                        )
                        if fragment_values.size > 0:
                            all_fragment_rg_per_rep.append(fragment_values)

            # --- Fragment-specific aggregation (counts + fragment histogram) ---
            if template.calculation_mode == "fragments":
                per_replicate_mean_fragments_per_frame = [
                    entry.mean_fragments_per_frame
                    for entry in run_entries
                    if entry.mean_fragments_per_frame is not None
                ]
                if per_replicate_mean_fragments_per_frame:
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
                    replicates=list(ctx.replicates),
                    n_replicates=len(ctx.replicates),
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
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            run_results=aggregated_runs,
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
            if agg_result is None:
                agg_dir = ctx.analysis_dirs[condition.label] / "aggregated"
                agg_result = self._load_aggregated_result(agg_dir)

            if agg_result is None:
                logger.warning(
                    "No aggregated Rg result for condition '%s'; skipping", condition.label
                )
                continue

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
            summaries_with_run = filter_summaries_with_run(
                summaries_by_label,
                run_label,
                lambda summary, label: summary.get_run(label),
                logger=logger,
            )
            if not summaries_with_run:
                logger.warning(
                    "Run '%s' has no condition data; skipping run-level comparison",
                    run_label,
                )
                continue

            run_labels.append(run_label)
            ranked = sorted(
                summaries_with_run,
                key=lambda label: summaries_with_run[label].get_run(run_label).mean_rg,
            )
            ranking_by_run[run_label] = ranked

            if len(summaries_with_run) >= 2:
                condition_pairs = build_condition_pairs(
                    list(summaries_with_run.keys()),
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
                            run_a=summaries_with_run[condition_a].get_run(run_label),
                            run_b=summaries_with_run[condition_b].get_run(run_label),
                        )
                    )

            if len(summaries_with_run) >= 3:
                groups = [
                    summary.get_run(run_label).per_replicate_means
                    for summary in summaries_with_run.values()
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

        if not run_labels:
            logger.warning("Rg comparison skipped because no runs had data in any condition")
            return None

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
    ) -> RgRunResult | None:
        """Compute one Rg run for a single replicate.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            1-indexed replicate number.
        run : RgRunSettings
            Run definition containing label and selection.
        loader : TrajectoryLoader
            Trajectory loader for this simulation condition.
        config_hash : str
            Hash of simulation configuration.
        eq_value : float
            Equilibration value.
        eq_unit : str
            Equilibration unit.
        start_frame : int
            First frame to include (0-indexed).
        n_frames_total : int
            Total number of frames in trajectory.
        n_frames_used : int
            Number of frames analyzed after equilibration.
        timestep_ps : float
            Trajectory timestep in picoseconds.

        Returns
        -------
        RgRunResult | None
            Per-run Rg result for this replicate, or ``None`` when the
            selection matches no atoms.
        """
        import numpy as np

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rg._results import RgRunResult
        from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

        u = loader.load_universe(replicate, cache=False)
        atom_group = u.select_atoms(run.selection)
        if len(atom_group) == 0:
            logger.warning(
                "Run '%s' selection matched no atoms for replicate %d: %r — skipping run",
                run.label,
                replicate,
                run.selection,
            )
            return None

        rg_values = np.empty(n_frames_used, dtype=np.float64)
        all_fragment_rg: list[float] = []
        fragment_counts: list[int] = []
        frag_masses: np.ndarray | None = None
        frag_metadata: dict[str, float | int] = {}

        if run.calculation_mode == "fragments":
            fragments = list(atom_group.fragments) if atom_group.fragments else [atom_group]
            if len(fragments) == 1:
                logger.warning(
                    "Run '%s' selection produced only 1 fragment; fragment mode will behave "
                    "identically to selection mode",
                    run.label,
                )

            if run.fragment_weighting == "mass":
                frag_masses = np.asarray(
                    [frag.total_mass() for frag in fragments], dtype=np.float64
                )
                if np.any(frag_masses <= 0) or np.any(np.isnan(frag_masses)):
                    raise ValueError(
                        f"Run '{run.label}': fragment masses contain zero or NaN values. "
                        "This suggests a problem with the MDAnalysis universe topology. "
                        f"Fragment masses: {frag_masses.tolist()}"
                    )

            for idx, frame_idx in enumerate(range(start_frame, n_frames_total)):
                u.trajectory[frame_idx]
                frag_rg = np.asarray(
                    [frag.radius_of_gyration() for frag in fragments], dtype=np.float64
                )
                fragment_counts.append(len(fragments))
                if run.save_fragment_distribution:
                    all_fragment_rg.extend(frag_rg.tolist())
                if frag_masses is not None:
                    rg_values[idx] = float(np.average(frag_rg, weights=frag_masses))
                else:
                    rg_values[idx] = float(np.mean(frag_rg))

            if all_fragment_rg:
                frag_arr = np.asarray(all_fragment_rg, dtype=np.float64)
                frag_metadata = {
                    "fragment_mean_rg": float(np.mean(frag_arr)),
                    "fragment_std_rg": float(np.std(frag_arr, ddof=0)),
                    "fragment_median_rg": float(np.median(frag_arr)),
                    "fragment_min_rg": float(np.min(frag_arr)),
                    "fragment_max_rg": float(np.max(frag_arr)),
                    "fragment_rg_p10": float(np.percentile(frag_arr, 10)),
                    "fragment_rg_p25": float(np.percentile(frag_arr, 25)),
                    "fragment_rg_p50": float(np.percentile(frag_arr, 50)),
                    "fragment_rg_p75": float(np.percentile(frag_arr, 75)),
                    "fragment_rg_p90": float(np.percentile(frag_arr, 90)),
                }

            frag_metadata["mean_fragments_per_frame"] = float(np.mean(fragment_counts))
            frag_metadata["min_fragments_per_frame"] = int(np.min(fragment_counts))
            frag_metadata["max_fragments_per_frame"] = int(np.max(fragment_counts))
        else:
            for idx, frame_idx in enumerate(range(start_frame, n_frames_total)):
                u.trajectory[frame_idx]
                rg_values[idx] = float(atom_group.radius_of_gyration())

        frames = np.arange(start_frame, n_frames_total, dtype=np.int64)
        time_ns = (frames.astype(np.float64) * timestep_ps) / 1000.0

        mean_rg = float(np.mean(rg_values))
        std_rg = float(np.std(rg_values, ddof=0))
        median_rg = float(np.median(rg_values))
        min_rg = float(np.min(rg_values))
        max_rg = float(np.max(rg_values))
        final_rg = float(rg_values[-1])

        sem_rg: float | None = None
        correlation_time: float | None = None
        correlation_time_unit: str | None = None
        n_independent_frames: int | None = None
        statistical_inefficiency: float | None = None
        autocorrelation_warning: str | None = None

        if len(rg_values) >= 20:
            tau_result = estimate_correlation_time(
                rg_values,
                timestep=timestep_ps,
                timestep_unit="ps",
                method="integration",
                n_frames=len(rg_values),
            )
            correlation_time = tau_result.tau
            correlation_time_unit = tau_result.tau_unit
            n_independent_frames = tau_result.n_independent
            statistical_inefficiency = tau_result.statistical_inefficiency
            autocorrelation_warning = tau_result.warning
            if n_independent_frames > 0:
                sem_rg = float(std_rg / np.sqrt(float(n_independent_frames)))

        npz_filename = f"rg_{run.label}_{eq_str}_{settings_tag}_timeseries.npz"
        npz_path = ctx.output_dir / npz_filename
        npz_data: dict[str, np.ndarray] = {
            "rg_values": rg_values,
            "time_ns": time_ns,
            "frames": frames,
        }
        if run.calculation_mode == "fragments" and run.save_fragment_distribution:
            npz_data["fragment_rg_values"] = np.asarray(all_fragment_rg, dtype=np.float64)
            npz_data["fragment_counts_per_frame"] = np.asarray(fragment_counts, dtype=np.int64)
            if frag_masses is not None:
                npz_data["fragment_masses"] = frag_masses

        np.savez_compressed(npz_path, **npz_data)

        return RgRunResult(
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
            calculation_mode=run.calculation_mode,
            fragment_weighting=(
                run.fragment_weighting if run.calculation_mode == "fragments" else None
            ),
            mean_rg=mean_rg,
            std_rg=std_rg,
            median_rg=median_rg,
            min_rg=min_rg,
            max_rg=max_rg,
            final_rg=final_rg,
            sem_rg=sem_rg,
            correlation_time_unit=correlation_time_unit,
            statistical_inefficiency=statistical_inefficiency,
            autocorrelation_warning=autocorrelation_warning,
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
            npz_path=str(npz_path),
            time_unit="ns",
            timestep_ps=timestep_ps,
            **frag_metadata,
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
