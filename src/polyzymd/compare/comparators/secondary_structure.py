"""Secondary structure comparator for comparing SS content across conditions.

This module provides the SecondaryStructureComparator class that orchestrates
DSSP analysis and statistical comparison of helix/strand/coil content across
multiple simulation conditions.

The comparator inherits from BaseComparator and implements the Template Method
pattern.  The **primary metric** is helix fraction (helix loss being the most
common signature of thermal unfolding in globular enzymes).
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
from polyzymd.compare.results.secondary_structure import (
    SSComparisonResult,
    SSConditionSummary,
)

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig
    from polyzymd.compare.settings import SecondaryStructureAnalysisSettings

logger = logging.getLogger("polyzymd.compare")

# Type alias for condition data dict returned by _load_or_compute
SSConditionData = dict[str, Any]


@ComparatorRegistry.register("secondary_structure")
class SecondaryStructureComparator(
    BaseComparator[
        "SecondaryStructureAnalysisSettings",
        SSConditionData,
        SSConditionSummary,
        SSComparisonResult,
    ]
):
    """Compare secondary structure content across multiple simulation conditions.

    Loads or computes DSSP secondary structure results for each condition,
    then performs statistical comparisons on helix fraction (primary metric).

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions.
    analysis_settings : SecondaryStructureAnalysisSettings
        SS analysis settings (from ``config.analysis_settings.get("secondary_structure")``).
    equilibration : str, optional
        Equilibration time override (e.g., ``"200ns"``).  If ``None``, uses
        ``config.defaults.equilibration_time``.

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> ss_settings = config.analysis_settings.get("secondary_structure")
    >>> comp = SecondaryStructureComparator(config, ss_settings)
    >>> result = comp.compare()
    >>> print(result.ranking)
    ["100% SBMA", "100% EGMA", "No Polymer"]
    """

    comparison_type: ClassVar[str] = "secondary_structure"

    def __init__(
        self,
        config: "ComparisonConfig",
        analysis_settings: "SecondaryStructureAnalysisSettings",
        equilibration: str | None = None,
    ):
        super().__init__(config, analysis_settings, equilibration)
        self.chain_id = analysis_settings.chain_id

    @property
    def metric_type(self) -> MetricType:
        """SS content is a mean-based metric.

        Helix/strand/coil fractions are frame averages.  The mean converges
        regardless of autocorrelation, but naive SEM underestimates
        uncertainty.  Use ``N_eff = N/g`` for corrected SEM.

        Returns
        -------
        MetricType
            ``MetricType.MEAN_BASED``
        """
        return MetricType.MEAN_BASED

    # ====================================================================
    # Abstract method implementations
    # ====================================================================

    def _load_or_compute(
        self,
        cond: "ConditionConfig",
        recompute: bool,
    ) -> SSConditionData:
        """Load existing SS results or compute them.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyse.
        recompute : bool
            Force recompute even if cached.

        Returns
        -------
        dict
            Keys: ``mean_helix``, ``sem_helix``, ``mean_strand``,
            ``sem_strand``, ``mean_coil``, ``sem_coil``, ``n_replicates``,
            ``per_replicate_helix``, ``per_replicate_strand``,
            ``per_replicate_coil``.
        """
        from polyzymd.analysis.secondary_structure import SecondaryStructureCalculator
        from polyzymd.analysis.secondary_structure.results import (
            SecondaryStructureAggregatedResult,
        )
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Processing condition: {cond.label}")

        sim_config = SimulationConfig.from_yaml(cond.config)

        # Resolve condition-specific output directory
        condition_output_dir = self._resolve_condition_output_dir(cond.label, "secondary_structure")

        # Try to find existing aggregated result
        result_path = self._find_aggregated_result(
            sim_config, cond.replicates, condition_output_dir=condition_output_dir
        )

        if result_path and result_path.exists() and not recompute:
            logger.info(f"  Loading cached result: {result_path}")
            agg_result = SecondaryStructureAggregatedResult.load(result_path)
        else:
            logger.info(f"  Computing SS for replicates {cond.replicates}...")
            calc = SecondaryStructureCalculator(
                config=sim_config,
                chain_id=self.chain_id,
                equilibration=self.equilibration,
            )
            agg_output_dir = condition_output_dir / "aggregated" if condition_output_dir else None
            agg_result = calc.compute_aggregated(
                replicates=cond.replicates,
                save=True,
                output_dir=agg_output_dir,
                recompute=recompute,
            )

        return {
            "mean_helix": agg_result.mean_overall_helix,
            "sem_helix": agg_result.sem_overall_helix,
            "mean_strand": agg_result.mean_overall_strand,
            "sem_strand": agg_result.sem_overall_strand,
            "mean_coil": agg_result.mean_overall_coil,
            "sem_coil": agg_result.sem_overall_coil,
            "n_replicates": agg_result.n_replicates,
            "per_replicate_helix": agg_result.per_replicate_helix,
            "per_replicate_strand": agg_result.per_replicate_strand,
            "per_replicate_coil": agg_result.per_replicate_coil,
        }

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: SSConditionData,
    ) -> SSConditionSummary:
        """Build an SS condition summary from raw data.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : dict
            Raw analysis data from ``_load_or_compute``.

        Returns
        -------
        SSConditionSummary
        """
        return SSConditionSummary(
            label=cond.label,
            config_path=str(cond.config),
            n_replicates=data["n_replicates"],
            # Primary metric for BaseConditionSummary: helix values
            replicate_values=data["per_replicate_helix"],
            # SS-specific fields
            mean_helix=data["mean_helix"],
            sem_helix=data["sem_helix"],
            mean_strand=data["mean_strand"],
            sem_strand=data["sem_strand"],
            mean_coil=data["mean_coil"],
            sem_coil=data["sem_coil"],
            per_replicate_helix=data["per_replicate_helix"],
            per_replicate_strand=data["per_replicate_strand"],
            per_replicate_coil=data["per_replicate_coil"],
        )

    def _build_result(
        self,
        summaries: list[SSConditionSummary],
        comparisons: list[Any],
        anova: ANOVASummary | None,
        ranking: list[str],
        effective_control: str | None,
        excluded_conditions: list["ConditionConfig"],
    ) -> SSComparisonResult:
        """Build the final SS comparison result.

        Parameters
        ----------
        summaries : list[SSConditionSummary]
            Condition summaries.
        comparisons : list
            Pairwise comparison results (on helix fraction).
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
        SSComparisonResult
        """
        return SSComparisonResult(
            metric="helix_fraction",
            name=self.config.name,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova,
            ranking=ranking,
            equilibration_time=self.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def _get_replicate_values(self, summary: SSConditionSummary) -> list[float]:
        """Extract per-replicate helix fractions for statistical tests."""
        return summary.per_replicate_helix

    def _get_mean_value(self, summary: SSConditionSummary) -> float:
        """Get the mean helix fraction."""
        return summary.mean_helix

    @property
    def _direction_labels(self) -> tuple[str, str, str]:
        """Positive helix change = more structure = stabilising."""
        return ("destabilizing", "unchanged", "stabilizing")

    def _rank_summaries(self, summaries: list[SSConditionSummary]) -> list[SSConditionSummary]:
        """Sort summaries by helix fraction (highest first = most structured)."""
        return sorted(summaries, key=lambda s: s.mean_helix, reverse=True)

    # ====================================================================
    # Helper methods
    # ====================================================================

    def _find_aggregated_result(
        self,
        sim_config: Any,
        replicates: list[int],
        condition_output_dir: Path | None = None,
    ) -> Path | None:
        """Find path to existing aggregated SS result.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        replicates : list[int]
            Replicate numbers.
        condition_output_dir : Path, optional
            Condition-specific output directory (from comparison mode).

        Returns
        -------
        Path or None
            Path to result file if it might exist.
        """
        from polyzymd.compare.comparators._utils import (
            format_replicate_range,
            parse_equilibration_time,
        )

        eq_value, eq_unit = parse_equilibration_time(self.equilibration)
        if eq_unit == "ps":
            eq_value = eq_value / 1000

        rep_str = format_replicate_range(replicates)
        filename = f"secondary_structure_{rep_str}_eq{eq_value:.0f}ns.json"

        # Check condition-specific path first (comparison mode)
        if condition_output_dir is not None:
            cond_path = condition_output_dir / "aggregated" / filename
            if cond_path.exists():
                return cond_path
            return None

        # Fallback to projects_directory (standalone mode only)
        result_path = (
            sim_config.output.projects_directory
            / "analysis"
            / "secondary_structure"
            / "aggregated"
            / filename
        )
        return result_path
