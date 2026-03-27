"""Distances analysis plugin.

Computes inter-atomic distances from MD trajectories for one or more
atom pairs, aggregates across replicates with SEM, and performs per-pair
cross-condition comparisons with t-tests, ANOVA, and rankings.

This plugin delegates heavy computation to
:class:`polyzymd.analysis.distances.calculator.DistanceCalculator`
and wraps the existing plotter classes for figure generation.

Unlike single-scalar analyses (RMSF, catalytic_triad), distances has **no
single primary metric** — each distance pair is compared independently since
averaging unrelated distances (e.g. H-bond distance + lid-opening distance)
is not semantically meaningful.  Therefore ``compare()`` is overridden
entirely and ``extract_metrics()`` is not used.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    PlotContext,
    ReplicateContext,
)

if TYPE_CHECKING:
    from polyzymd.analysis.results.distances import (
        DistanceAggregatedResult,
        DistanceResult,
    )
    from polyzymd.compare.results.distances import (
        DistanceComparisonResult,
    )

logger = logging.getLogger("polyzymd.analyses.distances")

# Default threshold from the existing settings module
DEFAULT_DISTANCE_THRESHOLD = 3.5


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class DistancePairSettings(BaseModel):
    """Configuration for a single distance pair.

    Attributes
    ----------
    label : str
        Human-readable label (e.g. ``"Ser77-Substrate"``).
    selection_a : str
        MDAnalysis selection string for the first atom/group.
    selection_b : str
        MDAnalysis selection string for the second atom/group.
    threshold : float | None
        Per-pair distance threshold (Angstroms).  If ``None``, uses the
        global threshold from :class:`DistancesSettings`.
    below_label : str | None
        Display label for the "below threshold" state (e.g. ``"Bound"``).
    above_label : str | None
        Display label for the "above threshold" state (e.g. ``"Unbound"``).
    """

    label: str = Field(..., description="Human-readable label for this pair")
    selection_a: str = Field(..., description="First atom/point selection")
    selection_b: str = Field(..., description="Second atom/point selection")
    threshold: float | None = Field(
        default=None,
        description="Per-pair threshold (Angstroms). Falls back to global threshold.",
    )
    below_label: str | None = Field(
        default=None,
        description='Display label for "below threshold" state.',
    )
    above_label: str | None = Field(
        default=None,
        description='Display label for "above threshold" state.',
    )


class DistancesSettings(BaseModel):
    """Settings for distances analysis.

    Attributes
    ----------
    threshold : float | None
        Global distance threshold for contact analysis (Angstroms).
    pairs : list[DistancePairSettings]
        Atom pairs to measure distances between.
    use_pbc : bool
        Use PBC-aware minimum image distances (default ``True``).
    align_trajectory : bool
        Align trajectory before distance calculation (default ``True``).
    alignment_selection : str
        MDAnalysis selection for trajectory alignment.
    alignment_mode : str
        Reference mode for alignment: ``"centroid"``, ``"average"``, or
        ``"frame"``.
    alignment_frame : int | None
        Reference frame (1-indexed) when ``alignment_mode="frame"``.
    """

    threshold: float | None = Field(
        default=DEFAULT_DISTANCE_THRESHOLD,
        description="Global distance threshold for contact analysis (Angstroms)",
    )
    pairs: list[DistancePairSettings] = Field(
        default_factory=list,
        description="Distance pairs to monitor",
    )
    use_pbc: bool = Field(
        default=True,
        description="Use PBC-aware minimum image distances",
    )
    align_trajectory: bool = Field(
        default=True,
        description="Align trajectory before distance calculation",
    )
    alignment_selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for trajectory alignment",
    )
    alignment_mode: str = Field(
        default="centroid",
        description="Reference mode: centroid, average, or frame",
    )
    alignment_frame: int | None = Field(
        default=None,
        description="Reference frame (1-indexed) when alignment_mode='frame'",
    )

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs(cls, v: list[DistancePairSettings]) -> list[DistancePairSettings]:
        """Ensure at least one pair is defined."""
        if len(v) == 0:
            raise ValueError("At least one distance pair must be defined")
        return v

    @field_validator("alignment_mode", mode="after")
    @classmethod
    def validate_alignment_mode(cls, v: str) -> str:
        """Validate alignment mode."""
        valid = {"centroid", "average", "frame"}
        if v not in valid:
            raise ValueError(f"alignment_mode must be one of {valid}, got '{v}'")
        return v

    @model_validator(mode="after")
    def validate_alignment_frame_required(self) -> "DistancesSettings":
        """Ensure alignment_frame is set when alignment_mode is 'frame'."""
        if (
            self.align_trajectory
            and self.alignment_mode == "frame"
            and self.alignment_frame is None
        ):
            raise ValueError("alignment_frame is required when alignment_mode is 'frame'")
        return self

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of ``(selection_a, selection_b)`` tuples."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]

    def get_pair_thresholds(self) -> list[float | None]:
        """Get per-pair thresholds, falling back to the global threshold."""
        return [p.threshold if p.threshold is not None else self.threshold for p in self.pairs]

    def get_alignment_config(self) -> Any:
        """Build an ``AlignmentConfig`` from these settings.

        Returns
        -------
        AlignmentConfig
            Configuration for trajectory alignment.
        """
        from polyzymd.analysis.core.alignment import AlignmentConfig

        return AlignmentConfig(
            enabled=self.align_trajectory,
            reference_mode=self.alignment_mode,
            reference_frame=self.alignment_frame,
            selection=self.alignment_selection,
        )


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class DistancesAnalysis(Analysis):
    """Distances analysis: inter-atomic distances from MD trajectories.

    This plugin wraps :class:`~polyzymd.analysis.distances.calculator.DistanceCalculator`
    for per-replicate computation, aggregates across replicates, and performs
    per-pair cross-condition comparison with dual metrics (mean distance and
    fraction below threshold).

    The ``compare()`` method is **fully overridden** because:

    - Each distance pair is compared independently (no single scalar).
    - Two metrics per pair: mean distance (lower is better) and fraction
      below threshold (higher is better, optional).
    - Rankings are per-pair, not overall.

    Plots
    -----
    - ``distance_kde_*.png``: KDE distribution overlay per pair.
    - ``distance_threshold_bars.png``: Grouped bar chart of fraction below
      threshold.
    - ``distance_state_*.png``: Per-pair above/below state bar charts.
    """

    name: ClassVar[str] = "distances"
    Settings: ClassVar[type] = DistancesSettings
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute distances for a single replicate.

        Delegates to :class:`~polyzymd.analysis.distances.calculator.DistanceCalculator`.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        DistanceResult
            Per-replicate distance result.
        """
        from polyzymd.analysis.distances.calculator import DistanceCalculator

        settings = ctx.settings

        pairs = settings.get_pair_selections()
        thresholds = settings.get_pair_thresholds()
        alignment = settings.get_alignment_config()

        calc = DistanceCalculator(
            config=ctx.sim_config,
            pairs=pairs,
            equilibration=ctx.equilibration,
            thresholds=thresholds,
            use_pbc=getattr(settings, "use_pbc", True),
            alignment=alignment,
        )

        result = calc.compute(
            replicate=replicate,
            save=True,
            output_dir=ctx.output_dir,
            recompute=ctx.recompute,
        )

        return result

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate distance results across replicates for one condition.

        Computes per-pair aggregated distance statistics from the
        already-computed per-replicate results. Does NOT re-run the
        calculator.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[DistanceResult]
            Per-replicate distance results.

        Returns
        -------
        DistanceAggregatedResult
            Aggregated result with per-pair statistics and SEM.
        """
        import numpy as np

        from polyzymd.analysis.core.aggregation import aggregate_distance_pair_stats
        from polyzymd.analysis.results.base import get_polyzymd_version
        from polyzymd.analysis.results.distances import (
            DistanceAggregatedResult,
            DistancePairAggregatedResult,
        )

        settings = ctx.settings
        first = results[0]

        n_pairs = len(first.pair_results)
        aggregated_pairs: list[DistancePairAggregatedResult] = []

        for pair_idx in range(n_pairs):
            stats = aggregate_distance_pair_stats(list(results), pair_idx)

            pr = first.pair_results[pair_idx]
            thresholds = settings.get_pair_thresholds()
            threshold = thresholds[pair_idx] if pair_idx < len(thresholds) else None

            agg_pair = DistancePairAggregatedResult(
                config_hash=first.config_hash,
                polyzymd_version=get_polyzymd_version(),
                replicate=None,
                equilibration_time=first.equilibration_time,
                equilibration_unit=first.equilibration_unit,
                selection_string=f"{pr.selection1} : {pr.selection2}",
                replicates=list(ctx.replicates),
                n_replicates=len(ctx.replicates),
                pair_label=pr.pair_label,
                selection1=pr.selection1,
                selection2=pr.selection2,
                overall_mean=stats.mean_stats.mean,
                overall_sem=stats.mean_stats.sem,
                overall_median=stats.median_stats.mean,
                per_replicate_means=stats.per_rep_means,
                per_replicate_stds=stats.per_rep_stds,
                per_replicate_medians=stats.per_rep_medians,
                threshold=threshold,
                overall_fraction_below=(
                    stats.fraction_stats.mean if stats.fraction_stats else None
                ),
                sem_fraction_below=(stats.fraction_stats.sem if stats.fraction_stats else None),
                per_replicate_fractions_below=(
                    stats.per_rep_fractions if stats.per_rep_fractions else None
                ),
                overall_kde_peak=(stats.kde_peak_stats.mean if stats.kde_peak_stats else None),
                sem_kde_peak=(stats.kde_peak_stats.sem if stats.kde_peak_stats else None),
                per_replicate_kde_peaks=(
                    stats.per_rep_kde_peaks if stats.per_rep_kde_peaks else None
                ),
            )
            aggregated_pairs.append(agg_pair)

        # Create aggregated result
        selection_strs = [f"({pr.selection1} : {pr.selection2})" for pr in first.pair_results]
        combined_selection = "; ".join(selection_strs)

        agg_result = DistanceAggregatedResult(
            config_hash=first.config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=first.equilibration_time,
            equilibration_unit=first.equilibration_unit,
            selection_string=combined_selection,
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            pair_results=aggregated_pairs,
            source_result_files=[],  # Not tracked in plugin mode
        )

        # Save
        filename = self._make_aggregated_filename(ctx.replicates, first)
        result_file = ctx.output_dir / filename
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        agg_result.save(result_file)
        logger.info(f"Saved aggregated distances to {result_file}")

        return agg_result

    # === compare() — fully overridden ===

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare distance metrics across conditions.

        Each distance pair is compared independently — rankings and statistics
        are computed per-pair since averaging unrelated distances is not
        semantically meaningful.

        Primary ranking: mean distance (lower = closer = better).
        Secondary ranking: fraction below threshold (higher = more contact = better).

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided context.

        Returns
        -------
        DistanceComparisonResult | None
            Comparison result, or ``None`` if fewer than 2 conditions.
        """
        from polyzymd import __version__
        from polyzymd.compare.results.distances import (
            DistanceComparisonResult,
            DistanceConditionSummary,
            DistancePairANOVA,
            DistancePairSummary,
            DistancePairwiseComparison,
        )
        from polyzymd.compare.statistics import (
            cohens_d,
            independent_ttest,
            one_way_anova,
            percent_change,
        )

        settings = ctx.settings
        pair_labels = settings.get_pair_labels()

        logger.info(f"Starting distance comparison: {ctx.name}")
        logger.info(f"Conditions: {len(ctx.conditions)}, Pairs: {pair_labels}")

        # Load aggregated results for each condition
        summaries: list[DistanceConditionSummary] = []

        for cond in ctx.conditions:
            agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
            agg_result = self._load_aggregated_result(agg_dir)

            if agg_result is None:
                logger.warning(f"No aggregated result found for '{cond.label}' — skipping.")
                continue

            # Map auto-generated pair labels to user-defined labels from settings
            selection_to_label: dict[tuple[str, str], str] = {
                (p.selection_a, p.selection_b): p.label for p in settings.pairs
            }

            pair_summaries = []
            for pr in agg_result.pair_results:
                user_label = selection_to_label.get((pr.selection1, pr.selection2), pr.pair_label)
                pair_summaries.append(
                    DistancePairSummary(
                        label=user_label,
                        selection_a=pr.selection1,
                        selection_b=pr.selection2,
                        threshold=pr.threshold,
                        mean_distance=pr.overall_mean,
                        sem_distance=pr.overall_sem,
                        fraction_below_threshold=pr.overall_fraction_below,
                        sem_fraction_below=pr.sem_fraction_below,
                        per_replicate_means=pr.per_replicate_means,
                        per_replicate_fractions=pr.per_replicate_fractions_below,
                    )
                )

            summaries.append(
                DistanceConditionSummary(
                    label=cond.label,
                    config_path=str(cond.config_path),
                    n_replicates=agg_result.n_replicates,
                    pair_summaries=pair_summaries,
                )
            )

        if len(summaries) < 2:
            logger.warning("distances: fewer than 2 conditions have results — skipping comparison.")
            return None

        # Resolve effective control
        effective_control = ctx.effective_control

        # Per-pair rankings
        ranking_by_pair: dict[str, list[str]] = {}
        fraction_ranking_by_pair: dict[str, list[str]] = {}

        for pair_label in pair_labels:
            pair_data = []
            for summary in summaries:
                pair_summary = summary.get_pair(pair_label)
                pair_data.append((summary.label, pair_summary))

            # Rank by mean distance (ascending — lower is better)
            sorted_by_distance = sorted(pair_data, key=lambda x: x[1].mean_distance)
            ranking_by_pair[pair_label] = [label for label, _ in sorted_by_distance]

            # Rank by fraction below threshold (descending — higher is better)
            with_fraction = [
                (label, ps) for label, ps in pair_data if ps.fraction_below_threshold is not None
            ]
            if with_fraction:
                sorted_by_fraction = sorted(
                    with_fraction,
                    key=lambda x: x[1].fraction_below_threshold or 0,
                    reverse=True,
                )
                fraction_ranking_by_pair[pair_label] = [label for label, _ in sorted_by_fraction]

        # Pairwise comparisons (per-pair)
        comparisons: list[DistancePairwiseComparison] = []

        for pair_label in pair_labels:
            if effective_control:
                control = next(s for s in summaries if s.label == effective_control)
                control_pair = control.get_pair(pair_label)
                for summary in summaries:
                    if summary.label == effective_control:
                        continue
                    treatment_pair = summary.get_pair(pair_label)
                    comp = self._compare_pair(
                        pair_label,
                        control.label,
                        summary.label,
                        control_pair,
                        treatment_pair,
                    )
                    comparisons.append(comp)
            else:
                for i, summary_a in enumerate(summaries):
                    pair_a = summary_a.get_pair(pair_label)
                    for summary_b in summaries[i + 1 :]:
                        pair_b = summary_b.get_pair(pair_label)
                        comp = self._compare_pair(
                            pair_label,
                            summary_a.label,
                            summary_b.label,
                            pair_a,
                            pair_b,
                        )
                        comparisons.append(comp)

        # ANOVA (if 3+ conditions) — per-pair
        anova_by_pair: list[DistancePairANOVA] | None = None
        if len(summaries) >= 3:
            anova_by_pair = []
            for pair_label in pair_labels:
                distance_groups = []
                fraction_groups = []
                for summary in summaries:
                    pair_data = summary.get_pair(pair_label)
                    distance_groups.append(pair_data.per_replicate_means)
                    if pair_data.per_replicate_fractions:
                        fraction_groups.append(pair_data.per_replicate_fractions)

                anova_dist = one_way_anova(*distance_groups)

                fraction_f = None
                fraction_p = None
                fraction_sig = None
                if len(fraction_groups) == len(summaries):
                    anova_frac = one_way_anova(*fraction_groups)
                    fraction_f = anova_frac.f_statistic
                    fraction_p = anova_frac.p_value
                    fraction_sig = anova_frac.significant

                anova_by_pair.append(
                    DistancePairANOVA(
                        pair_label=pair_label,
                        distance_f_statistic=anova_dist.f_statistic,
                        distance_p_value=anova_dist.p_value,
                        distance_significant=anova_dist.significant,
                        fraction_f_statistic=fraction_f,
                        fraction_p_value=fraction_p,
                        fraction_significant=fraction_sig,
                    )
                )

        # Build result
        result = DistanceComparisonResult(
            metric="mean_distance",
            name=ctx.name,
            n_pairs=len(settings.pairs),
            pair_labels=pair_labels,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova_by_pair=anova_by_pair,
            ranking_by_pair=ranking_by_pair,
            fraction_ranking_by_pair=(
                fraction_ranking_by_pair if fraction_ranking_by_pair else None
            ),
            equilibration_time=ctx.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

        # Log rankings
        for pair_label in pair_labels:
            logger.info(f"Ranking for '{pair_label}': {ranking_by_pair[pair_label]}")

        return result

    # === plot() ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate distance comparison plots.

        Delegates to the existing plotter classes:
        - :class:`DistanceKDEPlotter`
        - :class:`DistanceThresholdBarsPlotter`
        - :class:`DistanceStateBarsPlotter`

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

        # Build the data dict expected by the old plotter API
        data: dict[str, Any] = {}
        labels: list[str] = []

        for cond in ctx.conditions:
            label = cond.label
            labels.append(label)
            analysis_dir = ctx.analysis_dirs.get(label)
            if analysis_dir is not None:
                data[label] = {
                    "analysis_dir": analysis_dir,
                    "aggregated_dir": analysis_dir / "aggregated",
                    "replicates": list(cond.replicates),
                }

        if not labels:
            return plots

        # Add __meta__ with results_dir for comparison result lookup
        data["__meta__"] = {"results_dir": ctx.results_dir}

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve plot settings
        plot_settings = ctx.plot_settings
        if plot_settings is None:
            from polyzymd.compare.config import PlotSettings

            plot_settings = PlotSettings()

        # KDE plotter
        try:
            from polyzymd.compare.plotters.distances import DistanceKDEPlotter

            plotter = DistanceKDEPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Distance KDE plot failed: {exc}")

        # Threshold bars plotter
        try:
            from polyzymd.compare.plotters.distances import DistanceThresholdBarsPlotter

            plotter = DistanceThresholdBarsPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Distance threshold bars plot failed: {exc}")

        # State bars plotter
        try:
            from polyzymd.compare.plotters.distances import DistanceStateBarsPlotter

            plotter = DistanceStateBarsPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Distance state bars plot failed: {exc}")

        return plots

    # === Framework hooks ===

    def _deserialize_result(self, path: Path) -> Any:
        """Load an aggregated distance result from JSON.

        Parameters
        ----------
        path : Path
            Path to the JSON result file.

        Returns
        -------
        DistanceAggregatedResult
            Loaded result.
        """
        from polyzymd.analysis.results.distances import DistanceAggregatedResult

        return DistanceAggregatedResult.load(path)

    # === Private helpers ===

    @staticmethod
    def _compare_pair(
        pair_label: str,
        cond_a_label: str,
        cond_b_label: str,
        pair_a: Any,
        pair_b: Any,
    ) -> Any:
        """Compare two conditions statistically for a single distance pair.

        Parameters
        ----------
        pair_label : str
            Label of the distance pair being compared.
        cond_a_label : str
            Label of first condition (typically control).
        cond_b_label : str
            Label of second condition (typically treatment).
        pair_a : DistancePairSummary
            Pair data from condition A.
        pair_b : DistancePairSummary
            Pair data from condition B.

        Returns
        -------
        DistancePairwiseComparison
            Statistical comparison result.
        """
        from polyzymd.compare.results.distances import DistancePairwiseComparison
        from polyzymd.compare.statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )

        # Distance metric comparison
        values_a = pair_a.per_replicate_means
        values_b = pair_b.per_replicate_means

        ttest_dist = independent_ttest(values_a, values_b)
        effect_dist = cohens_d(values_a, values_b)
        pct_dist = percent_change(pair_a.mean_distance, pair_b.mean_distance)

        # Direction: negative change = closer = improving
        if pct_dist < -1:
            direction_dist = "closer"
        elif pct_dist > 1:
            direction_dist = "farther"
        else:
            direction_dist = "unchanged"

        # Fraction metric comparison (optional)
        fraction_t = None
        fraction_p = None
        fraction_d = None
        fraction_interp = None
        fraction_dir = None
        fraction_sig = None
        fraction_pct = None

        if pair_a.per_replicate_fractions and pair_b.per_replicate_fractions:
            frac_a = pair_a.per_replicate_fractions
            frac_b = pair_b.per_replicate_fractions

            ttest_frac = independent_ttest(frac_a, frac_b)
            effect_frac = cohens_d(frac_a, frac_b)
            pct_frac = percent_change(
                pair_a.fraction_below_threshold or 0,
                pair_b.fraction_below_threshold or 0,
            )

            fraction_t = ttest_frac.t_statistic
            fraction_p = ttest_frac.p_value
            fraction_d = effect_frac.cohens_d
            fraction_interp = effect_frac.interpretation
            fraction_sig = ttest_frac.significant

            # Direction: positive change = more contact = improving
            if pct_frac > 1:
                fraction_dir = "more_contact"
            elif pct_frac < -1:
                fraction_dir = "less_contact"
            else:
                fraction_dir = "unchanged"

            fraction_pct = pct_frac

        return DistancePairwiseComparison(
            pair_label=pair_label,
            condition_a=cond_a_label,
            condition_b=cond_b_label,
            distance_t_statistic=ttest_dist.t_statistic,
            distance_p_value=ttest_dist.p_value,
            distance_cohens_d=effect_dist.cohens_d,
            distance_effect_interpretation=effect_dist.interpretation,
            distance_direction=direction_dist,
            distance_significant=ttest_dist.significant,
            distance_percent_change=pct_dist,
            fraction_t_statistic=fraction_t,
            fraction_p_value=fraction_p,
            fraction_cohens_d=fraction_d,
            fraction_effect_interpretation=fraction_interp,
            fraction_direction=fraction_dir,
            fraction_significant=fraction_sig,
            fraction_percent_change=fraction_pct,
        )

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Generate filename for aggregated result JSON.

        Matches the naming convention of ``DistanceCalculator`` for backward
        compatibility.

        Parameters
        ----------
        replicates : tuple[int, ...] | Sequence[int]
            Replicate numbers included.
        first_result : DistanceResult
            First per-replicate result (for metadata).

        Returns
        -------
        str
            Filename like ``distances_reps1-5_eq100ns_pbc_align-centroid.json``.
        """
        eq_str = f"eq{first_result.equilibration_time:.0f}{first_result.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
        return f"distances_{rep_str}_{eq_str}.json"
