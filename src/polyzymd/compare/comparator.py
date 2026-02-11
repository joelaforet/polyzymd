"""RMSF comparator for comparing flexibility across conditions.

This module provides the RMSFComparator class that orchestrates
RMSF analysis and statistical comparison across multiple conditions.
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
    ANOVASummary,
    ComparisonResult,
    ConditionSummary,
    PairwiseComparison,
)
from polyzymd.compare.statistics import (
    cohens_d,
    independent_ttest,
    one_way_anova,
    percent_change,
)

LOGGER = logging.getLogger("polyzymd.compare")


class RMSFComparator:
    """Compare RMSF across multiple simulation conditions.

    This class loads RMSF results for each condition (computing them
    if necessary), then performs statistical comparisons including
    t-tests, ANOVA, and effect size calculations.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions
    equilibration : str, optional
        Equilibration time override (e.g., "10ns")
    selection : str, optional
        Atom selection override

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> comparator = RMSFComparator(config, equilibration="10ns")
    >>> result = comparator.compare()
    >>> print(result.ranking)
    ["100% SBMA", "100% EGMA", "No Polymer", "50/50 Mix"]
    """

    def __init__(
        self,
        config: ComparisonConfig,
        equilibration: Optional[str] = None,
        selection: Optional[str] = None,
    ):
        self.config = config
        self.equilibration = equilibration or config.defaults.equilibration_time
        self.selection = selection or config.defaults.selection

    def compare(self, recompute: bool = False) -> ComparisonResult:
        """Run comparison across all conditions.

        Parameters
        ----------
        recompute : bool, optional
            If True, force recompute RMSF even if cached results exist.
            Default is False.

        Returns
        -------
        ComparisonResult
            Complete comparison results with statistics and rankings
        """
        LOGGER.info(f"Starting RMSF comparison: {self.config.name}")
        LOGGER.info(f"Conditions: {len(self.config.conditions)}")
        LOGGER.info(f"Equilibration: {self.equilibration}")

        # 1. Load or compute RMSF for each condition
        condition_data: list[tuple[ConditionConfig, dict]] = []
        for cond in self.config.conditions:
            data = self._load_or_compute_rmsf(cond, recompute)
            condition_data.append((cond, data))

        # 2. Build condition summaries
        summaries = []
        for cond, data in condition_data:
            summary = ConditionSummary(
                label=cond.label,
                config_path=str(cond.config),
                n_replicates=data["n_replicates"],
                mean_rmsf=data["mean_rmsf"],
                sem_rmsf=data["sem_rmsf"],
                replicate_values=data["replicate_values"],
            )
            summaries.append(summary)

        # 3. Compute pairwise comparisons
        comparisons = self._compute_pairwise_comparisons(summaries)

        # 4. ANOVA if 3+ conditions
        anova = None
        if len(summaries) >= 3:
            anova = self._compute_anova(summaries)

        # 5. Rank by RMSF (lowest first = most stable)
        ranked = sorted(summaries, key=lambda s: s.mean_rmsf)
        ranking = [s.label for s in ranked]

        return ComparisonResult(
            metric="rmsf",
            name=self.config.name,
            control_label=self.config.control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova,
            ranking=ranking,
            equilibration_time=self.equilibration,
            selection=self.selection,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def _load_or_compute_rmsf(
        self,
        cond: ConditionConfig,
        recompute: bool,
    ) -> dict:
        """Load existing RMSF results or compute them.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze
        recompute : bool
            Force recompute even if cached

        Returns
        -------
        dict
            Dictionary with mean_rmsf, sem_rmsf, n_replicates, replicate_values
        """
        from polyzymd.analysis.rmsf import RMSFCalculator
        from polyzymd.analysis.results import RMSFAggregatedResult
        from polyzymd.config.schema import SimulationConfig

        LOGGER.info(f"Processing condition: {cond.label}")

        # Load simulation config
        sim_config = SimulationConfig.from_yaml(cond.config)

        # Try to find existing aggregated result
        result_path = self._find_aggregated_result(sim_config, cond.replicates)

        if result_path and result_path.exists() and not recompute:
            LOGGER.info(f"  Loading cached result: {result_path}")
            agg_result = RMSFAggregatedResult.load(result_path)
        else:
            # Compute RMSF
            LOGGER.info(f"  Computing RMSF for replicates {cond.replicates}...")
            calc = RMSFCalculator(
                config=sim_config,
                selection=self.selection,
                equilibration=self.equilibration,
            )
            agg_result = calc.compute_aggregated(
                replicates=cond.replicates,
                save=True,
                recompute=recompute,
            )

        return {
            "mean_rmsf": agg_result.overall_mean_rmsf,
            "sem_rmsf": agg_result.overall_sem_rmsf,
            "n_replicates": agg_result.n_replicates,
            "replicate_values": agg_result.per_replicate_mean_rmsf,
        }

    def _find_aggregated_result(
        self,
        sim_config: "SimulationConfig",
        replicates: list[int],
    ) -> Optional[Path]:
        """Find path to existing aggregated RMSF result.

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
        elif eq_str.endswith("ps"):
            eq_value = float(eq_str[:-2]) / 1000
        else:
            eq_value = float(eq_str)

        # Build expected filename
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))

        filename = f"rmsf_{rep_str}_eq{eq_value:.0f}ns.json"

        result_path = (
            sim_config.output.projects_directory / "analysis" / "rmsf" / "aggregated" / filename
        )

        return result_path

    def _compute_pairwise_comparisons(
        self,
        summaries: list[ConditionSummary],
    ) -> list[PairwiseComparison]:
        """Compute pairwise statistical comparisons.

        If a control is specified, compares all conditions vs control.
        Otherwise, compares all pairs.

        Parameters
        ----------
        summaries : list[ConditionSummary]
            Condition summaries

        Returns
        -------
        list[PairwiseComparison]
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
        cond_a: ConditionSummary,
        cond_b: ConditionSummary,
    ) -> PairwiseComparison:
        """Compare two conditions statistically.

        Parameters
        ----------
        cond_a : ConditionSummary
            First condition (typically control)
        cond_b : ConditionSummary
            Second condition (typically treatment)

        Returns
        -------
        PairwiseComparison
            Statistical comparison result
        """
        # T-test
        ttest = independent_ttest(cond_a.replicate_values, cond_b.replicate_values)

        # Effect size
        effect = cohens_d(cond_a.replicate_values, cond_b.replicate_values)

        # Percent change
        pct = percent_change(cond_a.mean_rmsf, cond_b.mean_rmsf)

        return PairwiseComparison(
            condition_a=cond_a.label,
            condition_b=cond_b.label,
            t_statistic=ttest.t_statistic,
            p_value=ttest.p_value,
            cohens_d=effect.cohens_d,
            effect_size_interpretation=effect.interpretation,
            direction=effect.direction,
            significant=ttest.significant,
            percent_change=pct,
        )

    def _compute_anova(
        self,
        summaries: list[ConditionSummary],
    ) -> ANOVASummary:
        """Compute one-way ANOVA across all conditions.

        Parameters
        ----------
        summaries : list[ConditionSummary]
            Condition summaries

        Returns
        -------
        ANOVASummary
            ANOVA result
        """
        groups = [s.replicate_values for s in summaries]
        result = one_way_anova(*groups)

        return ANOVASummary(
            f_statistic=result.f_statistic,
            p_value=result.p_value,
            significant=result.significant,
        )
