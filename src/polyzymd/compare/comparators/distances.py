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
from polyzymd.compare.core.base import BaseComparator
from polyzymd.compare.core.registry import ComparatorRegistry
from polyzymd.compare.results.distances import (
    DistanceComparisonResult,
    DistanceConditionSummary,
    DistancePairANOVA,
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
    from polyzymd.analysis.results.distances import (
        DistanceAggregatedResult,
        DistancePairAggregatedResult,
        DistanceResult,
    )
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")


# Type alias for condition data (dict returned by _load_or_compute)
DistanceConditionData = dict[str, Any]


@ComparatorRegistry.register("distances")
class DistancesComparator(
    BaseComparator[
        DistancesAnalysisSettings,
        DistanceConditionData,
        DistanceConditionSummary,
        DistanceComparisonResult,
    ]
):
    """Compare distance metrics across multiple simulation conditions.

    This class loads distance analysis results for each condition (computing them
    if necessary), then performs statistical comparisons including t-tests,
    ANOVA, and effect size calculations on both mean distance and fraction
    below threshold.

    Each distance pair is compared independently - there is no cross-pair
    averaging since different pairs measure fundamentally different physical
    quantities (e.g., H-bond distances vs lid-opening distances).

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
    >>> print(result.ranking_by_pair["Catalytic H-bond"])  # Per-pair ranking
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
        # Cast to concrete type if needed (before super().__init__)
        if not isinstance(analysis_settings, DistancesAnalysisSettings):
            analysis_settings = DistancesAnalysisSettings.model_validate(
                analysis_settings.model_dump()
            )
        super().__init__(config, analysis_settings, equilibration)

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

    @property
    def _direction_labels(self) -> tuple[str, str, str]:
        """Direction labels for distance comparisons.

        Returns
        -------
        tuple[str, str, str]
            (improving, unchanged, worsening) labels.

        Notes
        -----
        Distances uses inline direction logic in ``_compare_pair_data``
        rather than the base ``_interpret_direction`` method, because it
        has two independent metrics (distance and fraction) with opposite
        directions. This property satisfies the BaseComparator contract
        for interface consistency.
        """
        return ("closer", "unchanged", "farther")

    def compare(self, recompute: bool = False) -> DistanceComparisonResult:
        """Run the comparison across all conditions.

        Each distance pair is compared independently - rankings and statistics
        are computed per-pair since averaging unrelated distances (e.g., H-bond
        + lid-opening) is not semantically meaningful.

        Parameters
        ----------
        recompute : bool
            Force recompute even if cached results exist.

        Returns
        -------
        DistanceComparisonResult
            Complete comparison result with per-pair rankings.
        """
        pair_labels = self.analysis_settings.get_pair_labels()
        logger.info(f"Starting distance comparison: {self.config.name}")
        logger.info(f"Conditions: {len(self.config.conditions)}")
        logger.info(f"Equilibration: {self.equilibration}")
        logger.info(f"Pairs: {pair_labels}")

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
        effective_control = self._resolve_effective_control(summaries)

        # Compute per-pair rankings
        # For each pair, rank conditions by mean distance (ascending = lowest first)
        ranking_by_pair: dict[str, list[str]] = {}
        fraction_ranking_by_pair: dict[str, list[str]] = {}

        for pair_label in pair_labels:
            # Get pair data from each condition
            pair_data = []
            for summary in summaries:
                pair_summary = summary.get_pair(pair_label)
                pair_data.append((summary.label, pair_summary))

            # Rank by mean distance (ascending)
            sorted_by_distance = sorted(pair_data, key=lambda x: x[1].mean_distance)
            ranking_by_pair[pair_label] = [label for label, _ in sorted_by_distance]

            # Rank by fraction below threshold (descending) if threshold specified
            with_fraction = [
                (label, ps) for label, ps in pair_data if ps.fraction_below_threshold is not None
            ]
            if with_fraction:
                sorted_by_fraction = sorted(
                    with_fraction, key=lambda x: x[1].fraction_below_threshold or 0, reverse=True
                )
                fraction_ranking_by_pair[pair_label] = [label for label, _ in sorted_by_fraction]

        # Pairwise comparisons (now per-pair)
        comparisons = self._compute_distance_pairwise_comparisons(
            summaries, effective_control, pair_labels
        )

        # ANOVA (if 3+ conditions) - now per-pair
        anova_by_pair: list[DistancePairANOVA] | None = None
        if len(summaries) >= 3:
            anova_by_pair = self._compute_distance_anova(summaries, pair_labels)

        # Build result
        result = DistanceComparisonResult(
            metric="mean_distance",
            name=self.config.name,
            n_pairs=len(self.analysis_settings.pairs),
            pair_labels=pair_labels,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova_by_pair=anova_by_pair,
            ranking_by_pair=ranking_by_pair,
            fraction_ranking_by_pair=fraction_ranking_by_pair if fraction_ranking_by_pair else None,
            equilibration_time=self.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

        # Log rankings per pair
        for pair_label in pair_labels:
            logger.info(f"Ranking for '{pair_label}': {ranking_by_pair[pair_label]}")
            if pair_label in fraction_ranking_by_pair:
                logger.info(f"  Contact fraction ranking: {fraction_ranking_by_pair[pair_label]}")

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
            Dictionary with pair_summaries and n_replicates.
        """
        from polyzymd.analysis.distances.calculator import DistanceCalculator
        from polyzymd.analysis.results.distances import DistanceAggregatedResult
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Processing condition: {cond.label}")

        # Load simulation config
        sim_config = SimulationConfig.from_yaml(cond.config)

        # Resolve condition-specific output directory (None in standalone mode)
        condition_output_dir = self._resolve_condition_output_dir(cond.label, "distances")

        # Get per-pair thresholds
        pair_thresholds = self.analysis_settings.get_pair_thresholds()

        # Try to find existing aggregated result
        result_path = self._find_aggregated_result(
            sim_config, cond.replicates, condition_output_dir=condition_output_dir
        )

        agg_result: DistanceAggregatedResult | None = None

        if result_path and result_path.exists() and not recompute:
            logger.info(f"  Loading cached result: {result_path}")
            cached_result = DistanceAggregatedResult.load(result_path)

            # Validate and update thresholds if needed
            agg_result = self._update_aggregated_thresholds_if_needed(
                cached_result, sim_config, cond.replicates, pair_thresholds
            )
            if agg_result is None:
                # Threshold update failed, need full recompute
                logger.info("  Threshold update failed, forcing full recompute...")

        if agg_result is None:
            # Compute distance analysis
            logger.info(f"  Computing distance analysis for replicates {cond.replicates}...")

            # Get pair selections from settings
            pairs = self.analysis_settings.get_pair_selections()

            calculator = DistanceCalculator(
                config=sim_config,
                pairs=pairs,
                equilibration=self.equilibration,
                thresholds=pair_thresholds,
                use_pbc=self.analysis_settings.use_pbc,
                alignment=self.analysis_settings.get_alignment_config(),
            )
            agg_output_dir = condition_output_dir / "aggregated" if condition_output_dir else None
            agg_result = calculator.compute_aggregated(
                replicates=cond.replicates,
                save=True,
                output_dir=agg_output_dir,
                recompute=recompute,
            )

        # Build pair summaries from aggregated result
        # Note: We do NOT compute cross-pair averages here. Each pair is compared
        # independently since averaging unrelated distances (e.g., H-bond distance
        # + lid-opening distance) is not semantically meaningful.

        # Map auto-generated pair labels to user-defined labels from settings.
        # The DistanceCalculator auto-generates labels from selection strings
        # (e.g., "resid77_OG-RBY"), but comparison.yaml defines human-readable
        # labels (e.g., "Ser77-Substrate"). We match by selection strings.
        selection_to_label: dict[tuple[str, str], str] = {
            (p.selection_a, p.selection_b): p.label for p in self.analysis_settings.pairs
        }

        pair_summaries = []

        for pr in agg_result.pair_results:
            # Use the user-defined label if selections match, else keep auto-generated
            user_label = selection_to_label.get((pr.selection1, pr.selection2), pr.pair_label)
            pair_summary = DistancePairSummary(
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
            pair_summaries.append(pair_summary)

        return {
            "pair_summaries": pair_summaries,
            "n_replicates": agg_result.n_replicates,
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
        )

    def _resolve_effective_control(self, summaries: list[DistanceConditionSummary]) -> str | None:
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

    def _compute_distance_pairwise_comparisons(
        self,
        summaries: list[DistanceConditionSummary],
        control_label: str | None,
        pair_labels: list[str],
    ) -> list[DistancePairwiseComparison]:
        """Compute pairwise statistical comparisons for each distance pair.

        For each distance pair, compare conditions either all-vs-control or
        all pairwise combinations.

        Parameters
        ----------
        summaries : list[DistanceConditionSummary]
            Condition summaries.
        control_label : str or None
            Control condition label.
        pair_labels : list[str]
            Labels of distance pairs to compare.

        Returns
        -------
        list[DistancePairwiseComparison]
            Pairwise comparison results (one per pair per condition comparison).
        """
        comparisons = []

        for pair_label in pair_labels:
            if control_label:
                # Compare all vs control for this pair
                control = next(s for s in summaries if s.label == control_label)
                control_pair = control.get_pair(pair_label)
                for summary in summaries:
                    if summary.label == control_label:
                        continue
                    treatment_pair = summary.get_pair(pair_label)
                    comp = self._compare_pair_data(
                        pair_label=pair_label,
                        cond_a_label=control.label,
                        cond_b_label=summary.label,
                        pair_a=control_pair,
                        pair_b=treatment_pair,
                    )
                    comparisons.append(comp)
            else:
                # Compare all pairs of conditions for this distance pair
                for i, summary_a in enumerate(summaries):
                    pair_a = summary_a.get_pair(pair_label)
                    for summary_b in summaries[i + 1 :]:
                        pair_b = summary_b.get_pair(pair_label)
                        comp = self._compare_pair_data(
                            pair_label=pair_label,
                            cond_a_label=summary_a.label,
                            cond_b_label=summary_b.label,
                            pair_a=pair_a,
                            pair_b=pair_b,
                        )
                        comparisons.append(comp)

        return comparisons

    def _compare_pair_data(
        self,
        pair_label: str,
        cond_a_label: str,
        cond_b_label: str,
        pair_a: DistancePairSummary,
        pair_b: DistancePairSummary,
    ) -> DistancePairwiseComparison:
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
        # Distance metric comparison using per-replicate values
        values_a = pair_a.per_replicate_means
        values_b = pair_b.per_replicate_means

        ttest_dist = independent_ttest(values_a, values_b)
        effect_dist = cohens_d(values_a, values_b)
        pct_dist = percent_change(pair_a.mean_distance, pair_b.mean_distance)

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

        if pair_a.per_replicate_fractions and pair_b.per_replicate_fractions:
            frac_a = pair_a.per_replicate_fractions
            frac_b = pair_b.per_replicate_fractions

            ttest_frac = independent_ttest(frac_a, frac_b)
            effect_frac = cohens_d(frac_a, frac_b)
            pct_frac = percent_change(
                pair_a.fraction_below_threshold or 0, pair_b.fraction_below_threshold or 0
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
            pair_label=pair_label,
            condition_a=cond_a_label,
            condition_b=cond_b_label,
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

    def _compute_distance_anova(
        self,
        summaries: list[DistanceConditionSummary],
        pair_labels: list[str],
    ) -> list[DistancePairANOVA]:
        """Compute ANOVA across all conditions for each distance pair.

        Parameters
        ----------
        summaries : list[DistanceConditionSummary]
            Condition summaries.
        pair_labels : list[str]
            Labels of distance pairs.

        Returns
        -------
        list[DistancePairANOVA]
            ANOVA results for each pair.
        """
        anova_results = []

        for pair_label in pair_labels:
            # Get per-replicate values for this pair from each condition
            distance_groups = []
            fraction_groups = []

            for summary in summaries:
                pair_data = summary.get_pair(pair_label)
                distance_groups.append(pair_data.per_replicate_means)
                if pair_data.per_replicate_fractions:
                    fraction_groups.append(pair_data.per_replicate_fractions)

            # Distance ANOVA
            anova_dist = one_way_anova(*distance_groups)

            # Fraction ANOVA (if available for all conditions)
            fraction_f = None
            fraction_p = None
            fraction_sig = None

            if len(fraction_groups) == len(summaries):
                anova_frac = one_way_anova(*fraction_groups)
                fraction_f = anova_frac.f_statistic
                fraction_p = anova_frac.p_value
                fraction_sig = anova_frac.significant

            anova_results.append(
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

        return anova_results

    def _update_aggregated_thresholds_if_needed(
        self,
        agg_result: "DistanceAggregatedResult",
        sim_config: Any,
        replicates: list[int],
        expected_thresholds: list[float | None],
    ) -> "DistanceAggregatedResult | None":
        """Update contact fractions in aggregated result if thresholds changed.

        If the cached aggregated result used different thresholds than currently
        requested, attempts to reload individual replicate results and recompute
        the contact fractions from the stored distances. This avoids expensive
        full trajectory reprocessing when only threshold parameters change.

        Parameters
        ----------
        agg_result : DistanceAggregatedResult
            Cached aggregated result to potentially update.
        sim_config : SimulationConfig
            Simulation configuration for locating replicate results.
        replicates : list[int]
            Replicate numbers included in the aggregation.
        expected_thresholds : list[float | None]
            Expected thresholds for each pair (from analysis settings).

        Returns
        -------
        DistanceAggregatedResult or None
            Updated aggregated result with recomputed contact fractions,
            or None if the update failed and full recomputation is needed.
        """
        from polyzymd.analysis.results.distances import (
            DistanceAggregatedResult,
            DistancePairAggregatedResult,
            DistanceResult,
        )

        # Check if any thresholds mismatch
        needs_update = False
        for idx, pr in enumerate(agg_result.pair_results):
            expected = expected_thresholds[idx] if idx < len(expected_thresholds) else None
            cached = pr.threshold
            if expected != cached:
                needs_update = True
                logger.info(
                    f"Threshold mismatch for {pr.pair_label}: cached={cached}, expected={expected}"
                )
                break

        if not needs_update:
            return agg_result  # No update needed

        logger.info("Attempting to recompute contact fractions from cached replicate results...")

        # Try to load individual replicate results
        individual_results: list[DistanceResult] = []
        for rep in replicates:
            result_path = self._find_replicate_result(sim_config, rep)
            if result_path is None or not result_path.exists():
                logger.warning(
                    f"Cannot find replicate {rep} result file for threshold update. "
                    f"Full recomputation required."
                )
                return None

            try:
                result = DistanceResult.load(result_path)
                individual_results.append(result)
            except Exception as e:
                logger.warning(f"Failed to load replicate {rep} result: {e}")
                return None

        # Check that all replicate results have stored distances
        for result in individual_results:
            for pr in result.pair_results:
                if pr.distances is None or len(pr.distances) == 0:
                    logger.warning(
                        f"Replicate {result.replicate} pair {pr.pair_label} has no stored "
                        f"distances. Full recomputation required."
                    )
                    return None

        # Recompute aggregated pair results with new thresholds
        updated_pair_results: list[DistancePairAggregatedResult] = []

        for pair_idx, agg_pr in enumerate(agg_result.pair_results):
            new_threshold = (
                expected_thresholds[pair_idx] if pair_idx < len(expected_thresholds) else None
            )

            # Recompute per-replicate fractions from stored distances
            per_rep_fractions: list[float] = []
            for result in individual_results:
                pr = result.pair_results[pair_idx]
                if new_threshold is not None and pr.distances:
                    distances_arr = np.array(pr.distances)
                    fraction = float(np.mean(distances_arr < new_threshold))
                    per_rep_fractions.append(fraction)

            # Compute aggregated fraction statistics
            overall_fraction = None
            sem_fraction = None
            per_rep_fractions_out = None

            if per_rep_fractions and new_threshold is not None:
                overall_fraction = float(np.mean(per_rep_fractions))
                if len(per_rep_fractions) > 1:
                    sem_fraction = float(
                        np.std(per_rep_fractions, ddof=1) / np.sqrt(len(per_rep_fractions))
                    )
                else:
                    sem_fraction = 0.0
                per_rep_fractions_out = per_rep_fractions

            # Create updated pair result
            updated_pr = agg_pr.model_copy(
                update={
                    "threshold": new_threshold,
                    "overall_fraction_below": overall_fraction,
                    "sem_fraction_below": sem_fraction,
                    "per_replicate_fractions_below": per_rep_fractions_out,
                }
            )
            updated_pair_results.append(updated_pr)

        # Create updated aggregated result
        updated_agg = agg_result.model_copy(update={"pair_results": updated_pair_results})
        logger.info("Successfully recomputed contact fractions from cached replicate results.")

        return updated_agg

    def _find_replicate_result(
        self,
        sim_config: Any,
        replicate: int,
    ) -> Path | None:
        """Find path to existing single replicate distance result.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        replicate : int
            Replicate number.

        Returns
        -------
        Path or None
            Path to result file if it might exist.
        """
        # Parse equilibration time
        from polyzymd.compare.comparators._utils import parse_equilibration_time

        eq_value, eq_unit = parse_equilibration_time(self.equilibration)

        # Build expected filename pattern (matches _make_result_filename in calculator)
        pairs = self.analysis_settings.get_pair_selections()
        if pairs:
            # Create short label from first pair
            sel1, sel2 = pairs[0]
            # Simplified label extraction (matches calculator logic)
            import re

            def _sel_to_label(sel: str) -> str:
                label = sel.lower()
                label = re.sub(r"\b(and|or|not|protein)\b", "", label)
                resid_match = re.search(r"resid\s*(\d+)", label)
                name_match = re.search(r"name\s+(\w+)", label)
                parts = []
                if resid_match:
                    parts.append(f"resid{resid_match.group(1)}")
                if "midpoint" in sel.lower():
                    parts.append("mid")
                elif "com" in sel.lower():
                    parts.append("com")
                elif name_match:
                    parts.append(name_match.group(1).upper())
                if parts:
                    return "_".join(parts)
                label = re.sub(r"[^a-z0-9]+", "_", label)
                return label.strip("_")

            l1 = _sel_to_label(sel1)
            l2 = _sel_to_label(sel2)
            pair_label = f"{l1}-{l2}"
            if len(pairs) > 1:
                pair_label += f"_and{len(pairs) - 1}more"
        else:
            pair_label = "nopairs"

        filename = f"distances_{pair_label}_eq{eq_value:.0f}{eq_unit}.json"

        result_path = (
            sim_config.output.projects_directory
            / "analysis"
            / "distances"
            / f"run_{replicate}"
            / filename
        )

        return result_path

    def _find_aggregated_result(
        self,
        sim_config: Any,
        replicates: list[int],
        condition_output_dir: Path | None = None,
    ) -> Path | None:
        """Find path to existing aggregated distance result.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        replicates : list[int]
            Replicate numbers.
        condition_output_dir : Path, optional
            Condition-specific output directory (from comparison mode).
            Checked first before falling back to ``projects_directory``.

        Returns
        -------
        Path or None
            Path to result file if it might exist.
        """
        # Parse equilibration time
        from polyzymd.compare.comparators._utils import (
            format_replicate_range,
            parse_equilibration_time,
        )

        eq_value, eq_unit = parse_equilibration_time(self.equilibration)

        # Build expected filename pattern
        rep_str = format_replicate_range(replicates)

        # Build settings suffix to match DistanceCalculator._make_aggregated_filename()
        settings_parts = []

        # PBC setting
        pbc_str = "pbc" if self.analysis_settings.use_pbc else "nopbc"
        settings_parts.append(pbc_str)

        # Alignment setting
        alignment_config = self.analysis_settings.get_alignment_config()
        if alignment_config.enabled:
            align_str = f"align-{alignment_config.reference_mode}"
        else:
            align_str = "noalign"
        settings_parts.append(align_str)

        settings_suffix = "_".join(settings_parts)

        # The DistanceCalculator uses a pattern like:
        # distances_reps1-3_eq100ns_pbc_align-centroid.json
        filename = f"distances_{rep_str}_eq{eq_value:.0f}{eq_unit}_{settings_suffix}.json"

        # Check condition-specific path first (comparison mode)
        if condition_output_dir is not None:
            cond_path = condition_output_dir / "aggregated" / filename
            if cond_path.exists():
                return cond_path
            # In comparison mode, do NOT fall back to the shared
            # projects_directory — all conditions share the same path and
            # the cached file would belong to whichever condition wrote it
            # first.  Return None to trigger recomputation into the
            # condition-specific directory.
            return None

        # Fallback to projects_directory (standalone mode only)
        result_path = (
            sim_config.output.projects_directory
            / "analysis"
            / "distances"
            / "aggregated"
            / filename
        )

        return result_path
