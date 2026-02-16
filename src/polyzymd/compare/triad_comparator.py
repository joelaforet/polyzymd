"""Catalytic triad comparator for comparing active site geometry across conditions.

This module provides the TriadComparator class that orchestrates
catalytic triad analysis and statistical comparison across multiple conditions.

The key metric is "simultaneous contact fraction" - the percentage of frames
where ALL pairs in the triad are below the contact threshold simultaneously.
Higher values indicate better triad integrity and potentially better catalytic
competence.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

from polyzymd import __version__
from polyzymd.compare.config import ComparisonConfig, ConditionConfig
from polyzymd.compare.results import (
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
    one_way_anova,
    percent_change,
)

LOGGER = logging.getLogger("polyzymd.compare")


class TriadComparator:
    """Compare catalytic triad geometry across multiple simulation conditions.

    This class loads triad analysis results for each condition (computing them
    if necessary), then performs statistical comparisons including t-tests,
    ANOVA, and effect size calculations on the simultaneous contact fraction.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions.
    triad_settings : CatalyticTriadAnalysisSettings
        Catalytic triad analysis settings.
    equilibration : str, optional
        Equilibration time override (e.g., "10ns")

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

    def __init__(
        self,
        config: ComparisonConfig,
        triad_settings: CatalyticTriadAnalysisSettings,
        equilibration: str | None = None,
    ):
        self.config = config
        self.triad_settings = triad_settings
        self.equilibration = equilibration or config.defaults.equilibration_time

    def compare(self, recompute: bool = False) -> TriadComparisonResult:
        """Run comparison across all conditions.

        Parameters
        ----------
        recompute : bool, optional
            If True, force recompute triad analysis even if cached results exist.
            Default is False.

        Returns
        -------
        TriadComparisonResult
            Complete comparison results with statistics and rankings
        """
        LOGGER.info(f"Starting triad comparison: {self.config.name}")
        LOGGER.info(f"Triad: {self.triad_settings.name}")
        LOGGER.info(f"Conditions: {len(self.config.conditions)}")
        LOGGER.info(f"Equilibration: {self.equilibration}")

        # 1. Load or compute triad analysis for each condition
        condition_data: list[tuple[ConditionConfig, dict]] = []
        for cond in self.config.conditions:
            data = self._load_or_compute_triad(cond, recompute)
            condition_data.append((cond, data))

        # 2. Build condition summaries
        summaries = []
        for cond, data in condition_data:
            summary = TriadConditionSummary(
                label=cond.label,
                config_path=str(cond.config),
                n_replicates=data["n_replicates"],
                mean_simultaneous_contact=data["mean_simultaneous_contact"],
                sem_simultaneous_contact=data["sem_simultaneous_contact"],
                replicate_values=data["replicate_values"],
                pair_summaries=data["pair_summaries"],
            )
            summaries.append(summary)

        # 3. Compute pairwise comparisons
        comparisons = self._compute_pairwise_comparisons(summaries)

        # 4. ANOVA if 3+ conditions
        anova = None
        if len(summaries) >= 3:
            anova = self._compute_anova(summaries)

        # 5. Rank by simultaneous contact (highest first = best triad integrity)
        ranked = sorted(summaries, key=lambda s: s.mean_simultaneous_contact, reverse=True)
        ranking = [s.label for s in ranked]

        return TriadComparisonResult(
            metric="simultaneous_contact_fraction",
            name=self.config.name,
            triad_name=self.triad_settings.name,
            triad_description=self.triad_settings.description,
            threshold=self.triad_settings.threshold,
            n_pairs=self.triad_settings.n_pairs,
            pair_labels=self.triad_settings.get_pair_labels(),
            control_label=self.config.control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova,
            ranking=ranking,
            equilibration_time=self.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def _load_or_compute_triad(
        self,
        cond: ConditionConfig,
        recompute: bool,
    ) -> dict:
        """Load existing triad results or compute them.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze
        recompute : bool
            Force recompute even if cached

        Returns
        -------
        dict
            Dictionary with mean/sem simultaneous contact, n_replicates,
            replicate_values, and pair_summaries
        """
        from polyzymd.analysis.results.triad import TriadAggregatedResult
        from polyzymd.analysis.triad import CatalyticTriadAnalyzer
        from polyzymd.config.schema import SimulationConfig

        LOGGER.info(f"Processing condition: {cond.label}")

        # Load simulation config
        sim_config = SimulationConfig.from_yaml(cond.config)

        # Try to find existing aggregated result
        result_path = self._find_aggregated_result(sim_config, cond.replicates)

        if result_path and result_path.exists() and not recompute:
            LOGGER.info(f"  Loading cached result: {result_path}")
            agg_result = TriadAggregatedResult.load(result_path)
        else:
            # Compute triad analysis
            LOGGER.info(f"  Computing triad analysis for replicates {cond.replicates}...")
            analyzer = CatalyticTriadAnalyzer(
                config=sim_config,
                triad_config=self.triad_settings,
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

    def _find_aggregated_result(
        self,
        sim_config: "SimulationConfig",
        replicates: list[int],
    ) -> Optional[Path]:
        """Find path to existing aggregated triad result.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration
        replicates : list[int]
            Replicate numbers

        Returns
        -------
        Path or None
            Path to result file if it might exist
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

        # Build expected filename
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))

        name_safe = self.triad_settings.name.replace(" ", "_").replace("/", "-")
        filename = f"triad_{name_safe}_{rep_str}_eq{eq_value:.0f}{eq_unit}.json"

        result_path = (
            sim_config.output.projects_directory / "analysis" / "triad" / "aggregated" / filename
        )

        return result_path

    def _compute_pairwise_comparisons(
        self,
        summaries: list[TriadConditionSummary],
    ) -> list[TriadPairwiseComparison]:
        """Compute pairwise statistical comparisons.

        If a control is specified, compares all conditions vs control.
        Otherwise, compares all pairs.

        Parameters
        ----------
        summaries : list[TriadConditionSummary]
            Condition summaries

        Returns
        -------
        list[TriadPairwiseComparison]
            Pairwise comparison results
        """
        comparisons = []

        if self.config.control:
            # Compare all vs control
            control = next(s for s in summaries if s.label == self.config.control)
            treatments = [s for s in summaries if s.label != self.config.control]

            for treatment in treatments:
                comp = self._compare_pair(control, treatment)
                comparisons.append(comp)
        else:
            # Compare all pairs (first condition as reference for each)
            for i, cond_a in enumerate(summaries):
                for cond_b in summaries[i + 1 :]:
                    comp = self._compare_pair(cond_a, cond_b)
                    comparisons.append(comp)

        return comparisons

    def _compare_pair(
        self,
        cond_a: TriadConditionSummary,
        cond_b: TriadConditionSummary,
    ) -> TriadPairwiseComparison:
        """Compare two conditions statistically.

        Parameters
        ----------
        cond_a : TriadConditionSummary
            First condition (typically control)
        cond_b : TriadConditionSummary
            Second condition (typically treatment)

        Returns
        -------
        TriadPairwiseComparison
            Statistical comparison result
        """
        # T-test
        ttest = independent_ttest(cond_a.replicate_values, cond_b.replicate_values)

        # Effect size
        effect = cohens_d(cond_a.replicate_values, cond_b.replicate_values)

        # Percent change (note: for contact fraction, higher is better)
        pct = percent_change(cond_a.mean_simultaneous_contact, cond_b.mean_simultaneous_contact)

        # For triad contact: positive change = improving (more contact)
        # Negative change = worsening (less contact)
        if pct > 0:
            direction = "improving"
        elif pct < 0:
            direction = "worsening"
        else:
            direction = "unchanged"

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

    def _compute_anova(
        self,
        summaries: list[TriadConditionSummary],
    ) -> TriadANOVASummary:
        """Compute one-way ANOVA across all conditions.

        Parameters
        ----------
        summaries : list[TriadConditionSummary]
            Condition summaries

        Returns
        -------
        TriadANOVASummary
            ANOVA result
        """
        groups = [s.replicate_values for s in summaries]
        result = one_way_anova(*groups)

        return TriadANOVASummary(
            f_statistic=result.f_statistic,
            p_value=result.p_value,
            significant=result.significant,
        )
