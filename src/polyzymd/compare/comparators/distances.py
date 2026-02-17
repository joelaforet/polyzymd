"""Distances comparator for comparing distance metrics across conditions.

This module provides the DistancesComparator class that orchestrates
distance analysis and statistical comparison across multiple conditions.

The primary ranking metric is mean distance (lower = closer interactions).
Secondary metric is fraction below threshold (if threshold specified).

The comparator inherits from BaseComparator and implements the Template Method
pattern for DRY comparison logic.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from polyzymd import __version__
from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.compare.core.registry import ComparatorRegistry
from polyzymd.compare.results.distances import (
    DistanceANOVASummary,
    DistanceComparisonResult,
    DistanceConditionSummary,
    DistancePairSummary,
    DistancePairwiseComparison,
)
from polyzymd.compare.settings import DistancesAnalysisSettings
from polyzymd.compare.statistics import (
    cohens_d,
    independent_ttest,
    one_way_anova,
    percent_change,
)

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")


# Type alias for condition data (dict returned by _load_or_compute)
DistanceConditionData = dict[str, Any]


@ComparatorRegistry.register("distances")
class DistancesComparator:
    """Compare distance metrics across multiple simulation conditions.

    This class loads distance analysis results for each condition (computing them
    if necessary), then performs statistical comparisons including t-tests,
    ANOVA, and effect size calculations on both mean distance and fraction
    below threshold.

    The primary ranking metric is mean distance (lower = closer = better).

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions.
    analysis_settings : DistancesAnalysisSettings
        Distance analysis settings (from config.analysis_settings.get("distances")).
    equilibration : str, optional
        Equilibration time override (e.g., "10ns"). If None, uses
        config.defaults.equilibration_time.

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> dist_settings = config.analysis_settings.get("distances")
    >>> comparator = DistancesComparator(config, dist_settings, equilibration="10ns")
    >>> result = comparator.compare()
    >>> print(result.ranking)  # Sorted by mean distance (ascending)
    ["100% SBMA", "No Polymer", "50/50 Mix", "100% EGMA"]

    Notes
    -----
    Lower mean distance is better (closer interactions).
    Higher fraction below threshold is better (more time in contact).
    """

    comparison_type: ClassVar[str] = "distances"

    def __init__(
        self,
        config: "ComparisonConfig",
        analysis_settings: DistancesAnalysisSettings | Any,
        equilibration: str | None = None,
    ):
        self.config = config
        # Cast to concrete type if needed
        if not isinstance(analysis_settings, DistancesAnalysisSettings):
            self.analysis_settings = DistancesAnalysisSettings.model_validate(
                analysis_settings.model_dump()
            )
        else:
            self.analysis_settings = analysis_settings
        self.equilibration = equilibration or config.defaults.equilibration_time

    @classmethod
    def comparison_type_name(cls) -> str:
        """Return the comparison type identifier."""
        return "distances"

    @property
    def metric_type(self) -> MetricType:
        """Distance analysis is a mean-based metric.

        The mean distance is an average over frames. The mean converges
        regardless of autocorrelation, but we need to correct the uncertainty
        using N_eff (effective sample size = N/g where g is the statistical
        inefficiency).

        Returns
        -------
        MetricType
            MetricType.MEAN_BASED
        """
        return MetricType.MEAN_BASED

    def compare(self, recompute: bool = False) -> DistanceComparisonResult:
        """Run the comparison across all conditions.

        Parameters
        ----------
        recompute : bool
            Force recompute even if cached results exist.

        Returns
        -------
        DistanceComparisonResult
            Complete comparison result.
        """
        logger.info(f"Starting distance comparison: {self.config.name}")
        logger.info(f"Conditions: {len(self.config.conditions)}")
        logger.info(f"Equilibration: {self.equilibration}")
        logger.info(f"Pairs: {self.analysis_settings.get_pair_labels()}")

        # Load/compute data for each condition
        condition_data: list[tuple["ConditionConfig", DistanceConditionData]] = []
        for cond in self.config.conditions:
            try:
                data = self._load_or_compute(cond, recompute)
                condition_data.append((cond, data))
            except Exception as e:
                logger.warning(f"Skipping condition '{cond.label}': {e}")
                continue

        if len(condition_data) < 2:
            raise ValueError(
                f"Need at least 2 conditions for comparison, got {len(condition_data)}"
            )

        # Build condition summaries
        summaries = [self._build_condition_summary(cond, data) for cond, data in condition_data]

        # Determine control
        effective_control = self._get_effective_control(summaries)

        # Pairwise comparisons
        comparisons = self._compute_pairwise_comparisons(summaries, effective_control)

        # ANOVA (if 3+ conditions)
        anova = None
        if len(summaries) >= 3:
            anova = self._compute_anova(summaries)

        # Rank by mean distance (ascending = lowest first)
        ranking = [s.label for s in sorted(summaries, key=lambda s: s.overall_mean_distance)]

        # Rank by fraction (descending = highest first) if threshold specified
        ranking_by_fraction = None
        if self.analysis_settings.threshold is not None:
            # Filter to conditions with fraction data
            with_fraction = [s for s in summaries if s.overall_fraction_below is not None]
            if with_fraction:
                ranking_by_fraction = [
                    s.label
                    for s in sorted(
                        with_fraction, key=lambda s: s.overall_fraction_below or 0, reverse=True
                    )
                ]

        # Build result
        result = DistanceComparisonResult(
            metric="mean_distance",
            name=self.config.name,
            n_pairs=len(self.analysis_settings.pairs),
            pair_labels=self.analysis_settings.get_pair_labels(),
            threshold=self.analysis_settings.threshold,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova,
            ranking=ranking,
            ranking_by_fraction=ranking_by_fraction,
            equilibration_time=self.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

        logger.info(f"Comparison complete. Ranking by distance: {ranking}")
        if ranking_by_fraction:
            logger.info(f"Ranking by contact fraction: {ranking_by_fraction}")

        return result

    def _load_or_compute(
        self,
        cond: "ConditionConfig",
        recompute: bool,
    ) -> DistanceConditionData:
        """Load existing distance results or compute them.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze.
        recompute : bool
            Force recompute even if cached.

        Returns
        -------
        dict
            Dictionary with pair summaries, overall stats, and replicate values.
        """
        from polyzymd.analysis.distances.calculator import DistanceCalculator
        from polyzymd.analysis.results.distances import DistanceAggregatedResult
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Processing condition: {cond.label}")

        # Load simulation config
        sim_config = SimulationConfig.from_yaml(cond.config)

        # Try to find existing aggregated result
        result_path = self._find_aggregated_result(sim_config, cond.replicates)

        if result_path and result_path.exists() and not recompute:
            logger.info(f"  Loading cached result: {result_path}")
            agg_result = DistanceAggregatedResult.load(result_path)
        else:
            # Compute distance analysis
            logger.info(f"  Computing distance analysis for replicates {cond.replicates}...")

            # Get pair selections from settings
            pairs = self.analysis_settings.get_pair_selections()

            calculator = DistanceCalculator(
                config=sim_config,
                pairs=pairs,
                equilibration=self.equilibration,
                threshold=self.analysis_settings.threshold,
            )
            agg_result = calculator.compute_aggregated(
                replicates=cond.replicates,
                save=True,
                recompute=recompute,
            )

        # Build pair summaries from aggregated result
        pair_summaries = []
        per_replicate_overall_distances = []  # Average across pairs per replicate
        per_replicate_overall_fractions = []  # Average fraction per replicate

        # Initialize per-replicate accumulators
        n_reps = agg_result.n_replicates
        rep_distance_sums = [0.0] * n_reps
        rep_fraction_sums = [0.0] * n_reps
        n_pairs_with_fraction = 0

        for pr in agg_result.pair_results:
            pair_summary = DistancePairSummary(
                label=pr.pair_label,
                selection_a=pr.selection1,
                selection_b=pr.selection2,
                mean_distance=pr.overall_mean,
                sem_distance=pr.overall_sem,
                fraction_below_threshold=pr.overall_fraction_below,
                sem_fraction_below=pr.sem_fraction_below,
                per_replicate_means=pr.per_replicate_means,
            )
            pair_summaries.append(pair_summary)

            # Accumulate per-replicate values
            for i, mean in enumerate(pr.per_replicate_means):
                rep_distance_sums[i] += mean

            if pr.per_replicate_fractions_below:
                n_pairs_with_fraction += 1
                for i, frac in enumerate(pr.per_replicate_fractions_below):
                    rep_fraction_sums[i] += frac

        # Compute per-replicate overall averages
        n_pairs = len(agg_result.pair_results)
        per_replicate_overall_distances = [s / n_pairs for s in rep_distance_sums]

        if n_pairs_with_fraction > 0:
            per_replicate_overall_fractions = [s / n_pairs_with_fraction for s in rep_fraction_sums]
        else:
            per_replicate_overall_fractions = None

        # Compute overall statistics
        overall_mean_distance = float(np.mean(per_replicate_overall_distances))
        overall_sem_distance = float(
            np.std(per_replicate_overall_distances, ddof=1) / np.sqrt(n_reps)
        )

        overall_fraction_below = None
        overall_sem_fraction_below = None
        if per_replicate_overall_fractions:
            overall_fraction_below = float(np.mean(per_replicate_overall_fractions))
            overall_sem_fraction_below = float(
                np.std(per_replicate_overall_fractions, ddof=1) / np.sqrt(n_reps)
            )

        return {
            "pair_summaries": pair_summaries,
            "n_replicates": agg_result.n_replicates,
            "overall_mean_distance": overall_mean_distance,
            "overall_sem_distance": overall_sem_distance,
            "overall_fraction_below": overall_fraction_below,
            "overall_sem_fraction_below": overall_sem_fraction_below,
            "per_replicate_distances": per_replicate_overall_distances,
            "per_replicate_fractions": per_replicate_overall_fractions,
        }

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: DistanceConditionData,
    ) -> DistanceConditionSummary:
        """Build a distance condition summary from raw data.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : dict
            Raw analysis data from _load_or_compute.

        Returns
        -------
        DistanceConditionSummary
            Structured condition summary.
        """
        return DistanceConditionSummary(
            label=cond.label,
            config_path=str(cond.config),
            n_replicates=data["n_replicates"],
            pair_summaries=data["pair_summaries"],
            overall_mean_distance=data["overall_mean_distance"],
            overall_sem_distance=data["overall_sem_distance"],
            overall_fraction_below=data["overall_fraction_below"],
            overall_sem_fraction_below=data["overall_sem_fraction_below"],
            replicate_values=data["per_replicate_distances"],
            replicate_fractions=data["per_replicate_fractions"],
        )

    def _get_effective_control(self, summaries: list[DistanceConditionSummary]) -> str | None:
        """Determine the effective control condition.

        Parameters
        ----------
        summaries : list[DistanceConditionSummary]
            Condition summaries.

        Returns
        -------
        str or None
            Control label, or None if no control specified.
        """
        if self.config.control:
            # Verify control exists
            labels = [s.label for s in summaries]
            if self.config.control in labels:
                return self.config.control
            else:
                logger.warning(
                    f"Control '{self.config.control}' not found in conditions. Available: {labels}"
                )
        return None

    def _compute_pairwise_comparisons(
        self,
        summaries: list[DistanceConditionSummary],
        control_label: str | None,
    ) -> list[DistancePairwiseComparison]:
        """Compute pairwise statistical comparisons.

        If control is specified, compare all vs control.
        Otherwise, compare all pairs.

        Parameters
        ----------
        summaries : list[DistanceConditionSummary]
            Condition summaries.
        control_label : str or None
            Control condition label.

        Returns
        -------
        list[DistancePairwiseComparison]
            Pairwise comparison results.
        """
        comparisons = []

        if control_label:
            # Compare all vs control
            control = next(s for s in summaries if s.label == control_label)
            for summary in summaries:
                if summary.label == control_label:
                    continue
                comp = self._compare_pair(control, summary)
                comparisons.append(comp)
        else:
            # Compare all pairs
            for i, summary_a in enumerate(summaries):
                for summary_b in summaries[i + 1 :]:
                    comp = self._compare_pair(summary_a, summary_b)
                    comparisons.append(comp)

        return comparisons

    def _compare_pair(
        self,
        cond_a: DistanceConditionSummary,
        cond_b: DistanceConditionSummary,
    ) -> DistancePairwiseComparison:
        """Compare two conditions statistically.

        Parameters
        ----------
        cond_a : DistanceConditionSummary
            First condition (typically control).
        cond_b : DistanceConditionSummary
            Second condition (typically treatment).

        Returns
        -------
        DistancePairwiseComparison
            Statistical comparison result.
        """
        # Distance metric comparison
        values_a = cond_a.replicate_values
        values_b = cond_b.replicate_values

        ttest_dist = independent_ttest(values_a, values_b)
        effect_dist = cohens_d(values_a, values_b)
        pct_dist = percent_change(cond_a.overall_mean_distance, cond_b.overall_mean_distance)

        # Direction for distance: negative change = closer = improving
        if pct_dist < -1:  # 1% threshold for "closer"
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

        if cond_a.replicate_fractions and cond_b.replicate_fractions:
            frac_a = cond_a.replicate_fractions
            frac_b = cond_b.replicate_fractions

            ttest_frac = independent_ttest(frac_a, frac_b)
            effect_frac = cohens_d(frac_a, frac_b)
            pct_frac = percent_change(
                cond_a.overall_fraction_below or 0, cond_b.overall_fraction_below or 0
            )

            fraction_t = ttest_frac.t_statistic
            fraction_p = ttest_frac.p_value
            fraction_d = effect_frac.cohens_d
            fraction_interp = effect_frac.interpretation
            fraction_sig = ttest_frac.significant

            # Direction for fraction: positive change = more contact = improving
            if pct_frac > 1:
                fraction_dir = "more_contact"
            elif pct_frac < -1:
                fraction_dir = "less_contact"
            else:
                fraction_dir = "unchanged"

            fraction_pct = pct_frac

        return DistancePairwiseComparison(
            condition_a=cond_a.label,
            condition_b=cond_b.label,
            # Distance metric
            distance_t_statistic=ttest_dist.t_statistic,
            distance_p_value=ttest_dist.p_value,
            distance_cohens_d=effect_dist.cohens_d,
            distance_effect_interpretation=effect_dist.interpretation,
            distance_direction=direction_dist,
            distance_significant=ttest_dist.significant,
            distance_percent_change=pct_dist,
            # Fraction metric
            fraction_t_statistic=fraction_t,
            fraction_p_value=fraction_p,
            fraction_cohens_d=fraction_d,
            fraction_effect_interpretation=fraction_interp,
            fraction_direction=fraction_dir,
            fraction_significant=fraction_sig,
            fraction_percent_change=fraction_pct,
        )

    def _compute_anova(self, summaries: list[DistanceConditionSummary]) -> DistanceANOVASummary:
        """Compute ANOVA across all conditions.

        Parameters
        ----------
        summaries : list[DistanceConditionSummary]
            Condition summaries.

        Returns
        -------
        DistanceANOVASummary
            ANOVA results for both metrics.
        """
        # Distance ANOVA
        distance_groups = [s.replicate_values for s in summaries]
        anova_dist = one_way_anova(*distance_groups)

        # Fraction ANOVA (if available)
        fraction_f = None
        fraction_p = None
        fraction_sig = None

        fraction_groups = [s.replicate_fractions for s in summaries if s.replicate_fractions]
        if len(fraction_groups) == len(summaries):
            anova_frac = one_way_anova(*fraction_groups)
            fraction_f = anova_frac.f_statistic
            fraction_p = anova_frac.p_value
            fraction_sig = anova_frac.significant

        return DistanceANOVASummary(
            distance_f_statistic=anova_dist.f_statistic,
            distance_p_value=anova_dist.p_value,
            distance_significant=anova_dist.significant,
            fraction_f_statistic=fraction_f,
            fraction_p_value=fraction_p,
            fraction_significant=fraction_sig,
        )

    def _find_aggregated_result(
        self,
        sim_config: Any,
        replicates: list[int],
    ) -> Path | None:
        """Find path to existing aggregated distance result.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        replicates : list[int]
            Replicate numbers.

        Returns
        -------
        Path or None
            Path to result file if it might exist.
        """
        # Parse equilibration time
        eq_str = self.equilibration.lower()
        if eq_str.endswith("ns"):
            eq_value = float(eq_str[:-2])
            eq_unit = "ns"
        elif eq_str.endswith("ps"):
            eq_value = float(eq_str[:-2])
            eq_unit = "ps"
        else:
            eq_value = float(eq_str)
            eq_unit = "ns"

        # Build expected filename pattern
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))

        # The DistanceCalculator uses a pattern like:
        # distances_reps1-3_eq100ns.json
        filename = f"distances_{rep_str}_eq{eq_value:.0f}{eq_unit}.json"

        result_path = (
            sim_config.output.projects_directory
            / "analysis"
            / "distances"
            / "aggregated"
            / filename
        )

        return result_path
