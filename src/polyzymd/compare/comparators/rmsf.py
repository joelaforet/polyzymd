"""RMSF comparator for comparing flexibility across conditions.

This module provides the RMSFComparator class that orchestrates
RMSF analysis and statistical comparison across multiple conditions.

The comparator inherits from BaseComparator and implements the
Template Method pattern for DRY comparison logic.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

from polyzymd import __version__
from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.compare.core.base import ANOVASummary, BaseComparator
from polyzymd.compare.core.registry import ComparatorRegistry
from polyzymd.compare.results.rmsf import RMSFComparisonResult, RMSFConditionSummary
from polyzymd.compare.settings import RMSFAnalysisSettings

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")


# Type alias for condition data (dict returned by _load_or_compute_rmsf)
RMSFConditionData = dict[str, Any]


@ComparatorRegistry.register("rmsf")
class RMSFComparator(
    BaseComparator[
        RMSFAnalysisSettings, RMSFConditionData, RMSFConditionSummary, RMSFComparisonResult
    ]
):
    """Compare RMSF across multiple simulation conditions.

    This class loads RMSF results for each condition (computing them
    if necessary), then performs statistical comparisons including
    t-tests, ANOVA, and effect size calculations.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions.
    analysis_settings : RMSFAnalysisSettings
        RMSF analysis settings (from config.analysis_settings.get("rmsf")).
    equilibration : str, optional
        Equilibration time override (e.g., "10ns"). If None, uses
        config.defaults.equilibration_time.
    selection_override : str, optional
        Override for atom selection (requires --override flag on CLI).
    reference_mode_override : str, optional
        Override for reference mode (requires --override flag on CLI).
    reference_frame_override : int, optional
        Override for reference frame (requires --override flag on CLI).

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> rmsf_settings = config.analysis_settings.get("rmsf")
    >>> comparator = RMSFComparator(config, rmsf_settings, equilibration="10ns")
    >>> result = comparator.compare()
    >>> print(result.ranking)
    ["100% SBMA", "100% EGMA", "No Polymer", "50/50 Mix"]
    """

    comparison_type: ClassVar[str] = "rmsf"

    def __init__(
        self,
        config: "ComparisonConfig",
        analysis_settings: RMSFAnalysisSettings,
        equilibration: str | None = None,
        selection_override: str | None = None,
        reference_mode_override: str | None = None,
        reference_frame_override: int | None = None,
    ):
        super().__init__(config, analysis_settings, equilibration)

        # Apply overrides (CLI --override flag)
        self.selection = selection_override or analysis_settings.selection
        self.reference_mode = reference_mode_override or analysis_settings.reference_mode
        self.reference_frame = reference_frame_override or analysis_settings.reference_frame

    @classmethod
    def comparison_type_name(cls) -> str:
        """Return the comparison type identifier."""
        return "rmsf"

    @property
    def metric_type(self) -> MetricType:
        """RMSF is a variance-based metric.

        RMSF measures root-mean-square fluctuations, which are inherently
        variance-based. Correlated frames lead to biased variance estimates,
        so independent subsampling (2τ separation) is required for accurate
        uncertainty quantification.

        Returns
        -------
        MetricType
            MetricType.VARIANCE_BASED
        """
        return MetricType.VARIANCE_BASED

    # ========================================================================
    # Abstract Method Implementations
    # ========================================================================

    def _load_or_compute(
        self,
        cond: "ConditionConfig",
        recompute: bool,
    ) -> RMSFConditionData:
        """Load existing RMSF results or compute them.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze.
        recompute : bool
            Force recompute even if cached.

        Returns
        -------
        dict
            Dictionary with mean_rmsf, sem_rmsf, n_replicates, replicate_values.
        """
        from polyzymd.analysis.results import RMSFAggregatedResult
        from polyzymd.analysis.rmsf import RMSFCalculator
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Processing condition: {cond.label}")

        # Load simulation config
        sim_config = SimulationConfig.from_yaml(cond.config)

        # Try to find existing aggregated result
        result_path = self._find_aggregated_result(sim_config, cond.replicates)

        if result_path and result_path.exists() and not recompute:
            logger.info(f"  Loading cached result: {result_path}")
            agg_result = RMSFAggregatedResult.load(result_path)
        else:
            # Compute RMSF with full settings
            logger.info(f"  Computing RMSF for replicates {cond.replicates}...")
            calc = RMSFCalculator(
                config=sim_config,
                selection=self.selection,
                equilibration=self.equilibration,
                reference_mode=self.reference_mode,
                reference_frame=self.reference_frame,
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

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: RMSFConditionData,
    ) -> RMSFConditionSummary:
        """Build an RMSF condition summary from raw data.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : dict
            Raw analysis data from _load_or_compute.

        Returns
        -------
        RMSFConditionSummary
            Structured condition summary.
        """
        return RMSFConditionSummary(
            label=cond.label,
            config_path=str(cond.config),
            n_replicates=data["n_replicates"],
            mean_rmsf=data["mean_rmsf"],
            sem_rmsf=data["sem_rmsf"],
            replicate_values=data["replicate_values"],
        )

    def _build_result(
        self,
        summaries: list[RMSFConditionSummary],
        comparisons: list[Any],
        anova: ANOVASummary | None,
        ranking: list[str],
        effective_control: str | None,
        excluded_conditions: list["ConditionConfig"],
    ) -> RMSFComparisonResult:
        """Build the final RMSF comparison result.

        Parameters
        ----------
        summaries : list[RMSFConditionSummary]
            Condition summaries.
        comparisons : list
            Pairwise comparison results.
        anova : ANOVASummary or None
            ANOVA result.
        ranking : list[str]
            Ranked condition labels.
        effective_control : str or None
            Effective control label.
        excluded_conditions : list[ConditionConfig]
            Conditions that were excluded.

        Returns
        -------
        RMSFComparisonResult
            Complete comparison result.
        """
        return RMSFComparisonResult(
            metric="rmsf",
            name=self.config.name,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova,
            ranking=ranking,
            equilibration_time=self.equilibration,
            selection=self.selection,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def _get_replicate_values(self, summary: RMSFConditionSummary) -> list[float]:
        """Extract per-replicate RMSF values for statistical tests."""
        return summary.replicate_values

    def _get_mean_value(self, summary: RMSFConditionSummary) -> float:
        """Get the mean RMSF value."""
        return summary.mean_rmsf

    @property
    def _direction_labels(self) -> tuple[str, str, str]:
        """Negative RMSF change = lower flexibility = stabilizing."""
        return ("stabilizing", "unchanged", "destabilizing")

    def _rank_summaries(self, summaries: list[RMSFConditionSummary]) -> list[RMSFConditionSummary]:
        """Sort summaries by RMSF (lowest first = most stable)."""
        return sorted(summaries, key=lambda s: s.mean_rmsf)

    def _use_rmsf_mode_for_cohens_d(self) -> bool:
        """Use RMSF-specific Cohen's d interpretation."""
        return True

    # ========================================================================
    # Helper Methods
    # ========================================================================

    def _find_aggregated_result(
        self,
        sim_config: Any,
        replicates: list[int],
    ) -> Path | None:
        """Find path to existing aggregated RMSF result.

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
        # RMSF filenames always use ns — convert ps if needed
        if eq_unit == "ps":
            eq_value = eq_value / 1000

        # Build expected filename
        rep_str = format_replicate_range(replicates)

        filename = f"rmsf_{rep_str}_eq{eq_value:.0f}ns.json"

        result_path = (
            sim_config.output.projects_directory / "analysis" / "rmsf" / "aggregated" / filename
        )

        return result_path
