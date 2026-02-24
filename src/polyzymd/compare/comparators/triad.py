"""Catalytic triad comparator for comparing active site geometry across conditions.

This module provides the TriadComparator class that orchestrates
catalytic triad analysis and statistical comparison across multiple conditions.

The key metric is "simultaneous contact fraction" - the percentage of frames
where ALL pairs in the triad are below the contact threshold simultaneously.
Higher values indicate better triad integrity and potentially better catalytic
competence.

The comparator inherits from BaseComparator and implements the Template Method
pattern for DRY comparison logic.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

from polyzymd import __version__
from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.compare.core.base import ANOVASummary, BaseComparator, PairwiseComparison
from polyzymd.compare.core.registry import ComparatorRegistry
from polyzymd.compare.results.triad import (
    TriadANOVASummary,
    TriadComparisonResult,
    TriadConditionSummary,
    TriadPairSummary,
    TriadPairwiseComparison,
)
from polyzymd.compare.settings import CatalyticTriadAnalysisSettings
from polyzymd.compare.statistics import (
    cohens_d,
    independent_ttest,
    percent_change,
)

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")


# Type alias for condition data (dict returned by _load_or_compute)
TriadConditionData = dict[str, Any]


@ComparatorRegistry.register("triad")
class TriadComparator(
    BaseComparator[
        CatalyticTriadAnalysisSettings,
        TriadConditionData,
        TriadConditionSummary,
        TriadComparisonResult,
    ]
):
    """Compare catalytic triad geometry across multiple simulation conditions.

    This class loads triad analysis results for each condition (computing them
    if necessary), then performs statistical comparisons including t-tests,
    ANOVA, and effect size calculations on the simultaneous contact fraction.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions.
    analysis_settings : CatalyticTriadAnalysisSettings
        Catalytic triad analysis settings (from config.analysis_settings.get("catalytic_triad")).
    equilibration : str, optional
        Equilibration time override (e.g., "10ns"). If None, uses
        config.defaults.equilibration_time.

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> triad_settings = config.analysis_settings.get("catalytic_triad")
    >>> comparator = TriadComparator(config, triad_settings, equilibration="10ns")
    >>> result = comparator.compare()
    >>> print(result.ranking)
    ["100% SBMA", "100% EGMA", "No Polymer", "50/50 Mix"]

    Notes
    -----
    Higher simultaneous contact fraction is better (triad is more intact).
    """

    comparison_type: ClassVar[str] = "triad"

    def __init__(
        self,
        config: "ComparisonConfig",
        analysis_settings: CatalyticTriadAnalysisSettings,
        equilibration: str | None = None,
    ):
        super().__init__(config, analysis_settings, equilibration)

    @classmethod
    def comparison_type_name(cls) -> str:
        """Return the comparison type identifier."""
        return "triad"

    @property
    def metric_type(self) -> MetricType:
        """Catalytic triad contact fraction is a mean-based metric.

        The simultaneous contact fraction is an average over frames
        (fraction of frames where all pairs are in contact). The mean
        converges regardless of autocorrelation, but we need to correct
        the uncertainty using N_eff (effective sample size = N/g where
        g is the statistical inefficiency).

        Returns
        -------
        MetricType
            MetricType.MEAN_BASED
        """
        return MetricType.MEAN_BASED

    # ========================================================================
    # Abstract Method Implementations
    # ========================================================================

    def _load_or_compute(
        self,
        cond: "ConditionConfig",
        recompute: bool,
    ) -> TriadConditionData:
        """Load existing triad results or compute them.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze.
        recompute : bool
            Force recompute even if cached.

        Returns
        -------
        dict
            Dictionary with mean/sem simultaneous contact, n_replicates,
            replicate_values, and pair_summaries.
        """
        from polyzymd.analysis.results.triad import TriadAggregatedResult
        from polyzymd.analysis.triad import CatalyticTriadAnalyzer
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Processing condition: {cond.label}")

        # Load simulation config
        sim_config = SimulationConfig.from_yaml(cond.config)

        # Try to find existing aggregated result
        result_path = self._find_aggregated_result(sim_config, cond.replicates)

        if result_path and result_path.exists() and not recompute:
            logger.info(f"  Loading cached result: {result_path}")
            agg_result = TriadAggregatedResult.load(result_path)
        else:
            # Compute triad analysis
            logger.info(f"  Computing triad analysis for replicates {cond.replicates}...")
            analyzer = CatalyticTriadAnalyzer(
                config=sim_config,
                triad_config=self.analysis_settings,
                equilibration=self.equilibration,
            )
            agg_result = analyzer.compute_aggregated(
                replicates=cond.replicates,
                save=True,
                recompute=recompute,
            )

        # Build pair summaries
        pair_summaries = []
        for pr in agg_result.pair_results:
            pair_summary = TriadPairSummary(
                label=pr.pair_label,
                mean_distance=pr.overall_mean,
                sem_distance=pr.overall_sem,
                mean_fraction_below=pr.overall_fraction_below,
                sem_fraction_below=pr.sem_fraction_below,
            )
            pair_summaries.append(pair_summary)

        return {
            "mean_simultaneous_contact": agg_result.overall_simultaneous_contact,
            "sem_simultaneous_contact": agg_result.sem_simultaneous_contact,
            "n_replicates": agg_result.n_replicates,
            "replicate_values": agg_result.per_replicate_simultaneous,
            "pair_summaries": pair_summaries,
        }

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: TriadConditionData,
    ) -> TriadConditionSummary:
        """Build a triad condition summary from raw data.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : dict
            Raw analysis data from _load_or_compute.

        Returns
        -------
        TriadConditionSummary
            Structured condition summary.
        """
        return TriadConditionSummary(
            label=cond.label,
            config_path=str(cond.config),
            n_replicates=data["n_replicates"],
            mean_simultaneous_contact=data["mean_simultaneous_contact"],
            sem_simultaneous_contact=data["sem_simultaneous_contact"],
            replicate_values=data["replicate_values"],
            pair_summaries=data["pair_summaries"],
        )

    def _build_result(
        self,
        summaries: list[TriadConditionSummary],
        comparisons: list[Any],
        anova: ANOVASummary | None,
        ranking: list[str],
        effective_control: str | None,
        excluded_conditions: list["ConditionConfig"],
    ) -> TriadComparisonResult:
        """Build the final triad comparison result.

        Parameters
        ----------
        summaries : list[TriadConditionSummary]
            Condition summaries.
        comparisons : list
            Pairwise comparison results (TriadPairwiseComparison).
        anova : ANOVASummary or None
            ANOVA result (will be converted to TriadANOVASummary).
        ranking : list[str]
            Ranked condition labels.
        effective_control : str or None
            Effective control label.
        excluded_conditions : list[ConditionConfig]
            Conditions that were excluded.

        Returns
        -------
        TriadComparisonResult
            Complete comparison result.
        """
        # Convert base ANOVA to triad-specific format
        triad_anova = None
        if anova:
            triad_anova = TriadANOVASummary(
                f_statistic=anova.f_statistic,
                p_value=anova.p_value,
                significant=anova.significant,
            )

        return TriadComparisonResult(
            metric="simultaneous_contact_fraction",
            name=self.config.name,
            triad_name=self.analysis_settings.name,
            triad_description=self.analysis_settings.description,
            threshold=self.analysis_settings.threshold,
            n_pairs=self.analysis_settings.n_pairs,
            pair_labels=self.analysis_settings.get_pair_labels(),
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=triad_anova,
            ranking=ranking,
            equilibration_time=self.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def _get_replicate_values(self, summary: TriadConditionSummary) -> list[float]:
        """Extract per-replicate simultaneous contact values for statistical tests."""
        return summary.replicate_values

    def _get_mean_value(self, summary: TriadConditionSummary) -> float:
        """Get the mean simultaneous contact value."""
        return summary.mean_simultaneous_contact

    @property
    def _direction_labels(self) -> tuple[str, str, str]:
        """Positive triad contact change = improving (more contact)."""
        return ("worsening", "unchanged", "improving")

    def _rank_summaries(
        self, summaries: list[TriadConditionSummary]
    ) -> list[TriadConditionSummary]:
        """Sort summaries by simultaneous contact (highest first = best integrity)."""
        return sorted(summaries, key=lambda s: s.mean_simultaneous_contact, reverse=True)

    # ========================================================================
    # Override _compare_pair to return TriadPairwiseComparison
    # ========================================================================

    def _compare_pair(
        self,
        cond_a: TriadConditionSummary,
        cond_b: TriadConditionSummary,
    ) -> TriadPairwiseComparison:
        """Compare two conditions statistically.

        Overrides base to return TriadPairwiseComparison instead of PairwiseComparison.

        Parameters
        ----------
        cond_a : TriadConditionSummary
            First condition (typically control).
        cond_b : TriadConditionSummary
            Second condition (typically treatment).

        Returns
        -------
        TriadPairwiseComparison
            Statistical comparison result.
        """
        values_a = self._get_replicate_values(cond_a)
        values_b = self._get_replicate_values(cond_b)
        mean_a = self._get_mean_value(cond_a)
        mean_b = self._get_mean_value(cond_b)

        # T-test
        ttest = independent_ttest(values_a, values_b)

        # Effect size
        effect = cohens_d(values_a, values_b)

        # Percent change
        from polyzymd.compare.statistics import percent_change as pct_change

        pct = pct_change(mean_a, mean_b)

        # Direction interpretation
        direction = self._interpret_direction(pct)

        return TriadPairwiseComparison(
            condition_a=cond_a.label,
            condition_b=cond_b.label,
            t_statistic=ttest.t_statistic,
            p_value=ttest.p_value,
            cohens_d=effect.cohens_d,
            effect_size_interpretation=effect.interpretation,
            direction=direction,
            significant=ttest.significant,
            percent_change=pct,
        )

    # ========================================================================
    # Helper Methods
    # ========================================================================

    def _find_aggregated_result(
        self,
        sim_config: Any,
        replicates: list[int],
    ) -> Path | None:
        """Find path to existing aggregated triad result.

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
        from polyzymd.compare.comparators._utils import (
            format_replicate_range,
            parse_equilibration_time,
        )

        eq_value, eq_unit = parse_equilibration_time(self.equilibration)

        # Build expected filename
        rep_str = format_replicate_range(replicates)

        name_safe = self.analysis_settings.name.replace(" ", "_").replace("/", "-")
        filename = f"triad_{name_safe}_{rep_str}_eq{eq_value:.0f}{eq_unit}.json"

        result_path = (
            sim_config.output.projects_directory
            / "analysis"
            / "catalytic_triad"
            / "aggregated"
            / filename
        )

        return result_path
