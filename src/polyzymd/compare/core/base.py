"""Base classes for comparison analysis.

This module provides abstract base classes that consolidate common patterns
across all comparator types, following the Template Method design pattern.

Classes
-------
BaseConditionSummary
    Abstract base for condition-level summary statistics.
BaseComparisonResult
    Abstract base for complete comparison results with save/load.
PairwiseComparison
    Shared model for statistical comparison between two conditions.
ANOVASummary
    Shared model for ANOVA results.
BaseComparator
    Abstract base implementing the Template Method pattern for comparisons.

Design Principles
-----------------
1. Open-Closed Principle: New comparators extend base classes without modifying them.
2. Template Method: `compare()` defines the algorithm skeleton; subclasses fill in specifics.
3. DRY: Statistical tests, pairwise logic, and serialization are implemented once.
"""

from __future__ import annotations

import json
import logging
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Generic, Self, TypeVar

from pydantic import BaseModel, Field

from polyzymd import __version__
from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.compare.statistics import (
    cohens_d,
    independent_ttest,
    one_way_anova,
    percent_change,
)

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")


# ============================================================================
# Shared Result Models (DRY - used by all comparators)
# ============================================================================


class PairwiseComparison(BaseModel):
    """Statistical comparison between two conditions.

    This is the standard pairwise comparison result used across all
    comparator types. For comparators that need additional fields
    (e.g., multiple metrics), subclass this model.

    Attributes
    ----------
    condition_a : str
        Label of first condition (typically control/reference).
    condition_b : str
        Label of second condition (typically treatment).
    metric : str
        Name of the metric being compared.
    t_statistic : float
        T-test statistic.
    p_value : float
        Two-tailed p-value.
    cohens_d : float
        Effect size (Cohen's d).
    effect_size_interpretation : str
        "negligible", "small", "medium", or "large".
    direction : str
        Interpretation of change (e.g., "stabilizing", "improving").
    significant : bool
        Whether p < 0.05.
    percent_change : float
        Percent change from condition_a to condition_b.
    """

    condition_a: str
    condition_b: str
    metric: str = "default"
    t_statistic: float
    p_value: float
    cohens_d: float
    effect_size_interpretation: str
    direction: str
    significant: bool
    percent_change: float


class ANOVASummary(BaseModel):
    """One-way ANOVA result summary.

    Attributes
    ----------
    metric : str
        Name of the metric tested (e.g., "rmsf", "coverage").
    f_statistic : float
        F-statistic from ANOVA.
    p_value : float
        P-value for the test.
    significant : bool
        Whether p < 0.05.
    """

    metric: str = "default"
    f_statistic: float
    p_value: float
    significant: bool


# ============================================================================
# Abstract Base Classes for Results
# ============================================================================


class BaseConditionSummary(BaseModel, ABC):
    """Abstract base class for condition-level summary statistics.

    All condition summaries share these common fields. Subclasses add
    analysis-specific fields (e.g., mean_rmsf, coverage_mean).

    Attributes
    ----------
    label : str
        Display name for this condition.
    config_path : str
        Path to the simulation config file.
    n_replicates : int
        Number of replicates included.
    replicate_values : list[float]
        Per-replicate values of the primary metric (for statistical tests).
    """

    label: str
    config_path: str
    n_replicates: int
    replicate_values: list[float]

    @property
    @abstractmethod
    def primary_metric_value(self) -> float:
        """Return the primary metric value for ranking/comparison.

        This is used by BaseComparator for sorting and statistical tests.
        """
        ...

    @property
    @abstractmethod
    def primary_metric_sem(self) -> float:
        """Return the SEM of the primary metric."""
        ...


# Type variable for condition summary subtypes
TConditionSummary = TypeVar("TConditionSummary", bound=BaseConditionSummary)
TPairwiseComparison = TypeVar("TPairwiseComparison", bound=PairwiseComparison)


class BaseComparisonResult(BaseModel, ABC, Generic[TConditionSummary, TPairwiseComparison]):
    """Abstract base class for comparison results.

    Provides common serialization (save/load) and accessor methods.
    Subclasses define analysis-specific fields.

    Attributes
    ----------
    metric : str
        The primary metric being compared (e.g., "rmsf", "simultaneous_contact_fraction").
    name : str
        Name of the comparison project.
    control_label : str, optional
        Label of the control condition.
    conditions : list[TConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[TPairwiseComparison]
        Statistical comparisons (all vs control, or all pairs).
    anova : ANOVASummary, optional
        ANOVA result if 3+ conditions.
    ranking : list[str]
        Labels sorted by primary metric.
    equilibration_time : str
        Equilibration time used.
    created_at : datetime
        When the analysis was run.
    polyzymd_version : str
        Version of polyzymd used.
    """

    # Class variable - subclasses should override
    comparison_type: ClassVar[str] = "base"

    metric: str
    name: str
    control_label: str | None = None
    conditions: list[Any]  # Will be overridden in subclasses with specific type
    pairwise_comparisons: list[Any]  # Will be overridden in subclasses
    anova: ANOVASummary | list[ANOVASummary] | None = None
    ranking: list[str]
    equilibration_time: str
    created_at: datetime
    polyzymd_version: str = __version__

    def save(self, path: Path | str) -> Path:
        """Save result to JSON file.

        Parameters
        ----------
        path : Path or str
            Output path.

        Returns
        -------
        Path
            Path to saved file.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path | str) -> Self:
        """Load result from JSON file.

        Parameters
        ----------
        path : Path or str
            Path to JSON file.

        Returns
        -------
        Self
            Loaded result.
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())

    def get_condition(self, label: str) -> Any:
        """Get a condition by label.

        Parameters
        ----------
        label : str
            Condition label.

        Returns
        -------
        BaseConditionSummary
            The matching condition.

        Raises
        ------
        KeyError
            If condition not found.
        """
        for cond in self.conditions:
            if cond.label == label:
                return cond
        raise KeyError(f"Condition '{label}' not found")

    def get_comparison(self, label: str) -> Any | None:
        """Get pairwise comparison for a condition vs control.

        Parameters
        ----------
        label : str
            Treatment condition label.

        Returns
        -------
        PairwiseComparison or None
            The comparison, or None if not found.
        """
        for comp in self.pairwise_comparisons:
            if comp.condition_b == label:
                return comp
        return None


# ============================================================================
# Base Comparator (Template Method Pattern)
# ============================================================================


# Type variables for generic comparator
TAnalysisSettings = TypeVar("TAnalysisSettings")
TComparisonSettings = TypeVar("TComparisonSettings")
TConditionData = TypeVar("TConditionData")
TResult = TypeVar("TResult", bound=BaseComparisonResult)


class BaseComparator(ABC, Generic[TAnalysisSettings, TConditionData, TConditionSummary, TResult]):
    """Abstract base class for all comparators using Template Method pattern.

    The `compare()` method defines the comparison algorithm skeleton:
    1. Load/compute analysis for each condition
    2. Build condition summaries
    3. Compute pairwise statistical comparisons
    4. Compute ANOVA (if 3+ conditions)
    5. Rank conditions
    6. Build and return result

    Subclasses implement the abstract methods to customize each step.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions.
    analysis_settings : TAnalysisSettings
        Analysis-specific settings.
    equilibration : str, optional
        Equilibration time override.

    Type Parameters
    ---------------
    TAnalysisSettings
        Type of analysis settings (e.g., RMSFAnalysisSettings).
    TConditionData
        Type of raw data loaded for each condition.
    TConditionSummary
        Type of condition summary (e.g., RMSFConditionSummary).
    TResult
        Type of comparison result (e.g., RMSFComparisonResult).
    """

    # Class variable - subclasses should override
    comparison_type: ClassVar[str] = "base"

    def __init__(
        self,
        config: "ComparisonConfig",
        analysis_settings: TAnalysisSettings,
        equilibration: str | None = None,
    ):
        self.config = config
        self.analysis_settings = analysis_settings
        self.equilibration = equilibration or config.defaults.equilibration_time

    @classmethod
    @abstractmethod
    def comparison_type_name(cls) -> str:
        """Return the comparison type identifier (e.g., "rmsf", "contacts").

        Returns
        -------
        str
            Type identifier used in registry and CLI.
        """
        ...

    @property
    @abstractmethod
    def metric_type(self) -> MetricType:
        """Declare whether this comparator's metric is mean or variance-based.

        This determines how autocorrelation is handled in the underlying analysis:

        - **MEAN_BASED**: Use all frames for computation, correct uncertainty
          using N_eff (effective sample size). Examples: average distance,
          contact fraction, catalytic triad proximity.

        - **VARIANCE_BASED**: Subsample to independent frames separated by 2τ
          (correlation time) to avoid bias in variance estimates. Examples:
          RMSF, fluctuation metrics.

        Contributors implementing new comparators MUST declare the appropriate
        metric type to ensure correct statistical treatment per LiveCoMS
        best practices (Grossfield et al., 2018).

        Returns
        -------
        MetricType
            The metric type for this comparator.

        References
        ----------
        - Grossfield et al. (2018) LiveCoMS 1:5067 (Best Practices for Uncertainty)
        - GitHub: dmzuckerman/Sampling-Uncertainty
        """
        ...

    def compare(self, recompute: bool = False) -> TResult:
        """Run comparison across all conditions (Template Method).

        This method defines the algorithm skeleton. Subclasses customize
        behavior by implementing the abstract hook methods.

        Parameters
        ----------
        recompute : bool, optional
            If True, force recompute even if cached results exist.

        Returns
        -------
        TResult
            Complete comparison results with statistics and rankings.
        """
        logger.info(f"Starting {self.comparison_type_name()} comparison: {self.config.name}")
        logger.info(f"Conditions: {len(self.config.conditions)}")
        logger.info(f"Equilibration: {self.equilibration}")

        # Step 1: Filter conditions (optional hook - default returns all)
        valid_conditions, excluded_conditions = self._filter_conditions()

        if excluded_conditions:
            logger.warning(
                f"Excluding {len(excluded_conditions)} condition(s): "
                f"{[c.label for c in excluded_conditions]}"
            )

        # Step 2: Load or compute analysis for each condition
        condition_data: list[tuple["ConditionConfig", TConditionData]] = []
        for cond in valid_conditions:
            data = self._load_or_compute(cond, recompute)
            condition_data.append((cond, data))

        # Step 3: Build condition summaries
        summaries: list[TConditionSummary] = []
        for cond, data in condition_data:
            summary = self._build_condition_summary(cond, data)
            summaries.append(summary)

        # Step 4: Determine effective control
        effective_control = self._get_effective_control(excluded_conditions)

        # Step 5: Compute pairwise comparisons
        comparisons = self._compute_pairwise_comparisons(summaries, effective_control)

        # Step 6: ANOVA if 3+ conditions
        anova = None
        if len(summaries) >= 3:
            anova = self._compute_anova(summaries)

        # Step 7: Rank conditions
        ranking = self._compute_ranking(summaries)

        # Step 8: Build result
        return self._build_result(
            summaries=summaries,
            comparisons=comparisons,
            anova=anova,
            ranking=ranking,
            effective_control=effective_control,
            excluded_conditions=excluded_conditions,
        )

    # ========================================================================
    # Abstract Methods (must be implemented by subclasses)
    # ========================================================================

    @abstractmethod
    def _load_or_compute(
        self,
        cond: "ConditionConfig",
        recompute: bool,
    ) -> TConditionData:
        """Load existing results or compute analysis for a condition.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze.
        recompute : bool
            Force recompute even if cached.

        Returns
        -------
        TConditionData
            Raw analysis data for this condition.
        """
        ...

    @abstractmethod
    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: TConditionData,
    ) -> TConditionSummary:
        """Build a condition summary from raw data.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : TConditionData
            Raw analysis data.

        Returns
        -------
        TConditionSummary
            Structured condition summary.
        """
        ...

    @abstractmethod
    def _build_result(
        self,
        summaries: list[TConditionSummary],
        comparisons: list[Any],
        anova: ANOVASummary | list[ANOVASummary] | None,
        ranking: list[str],
        effective_control: str | None,
        excluded_conditions: list["ConditionConfig"],
    ) -> TResult:
        """Build the final comparison result.

        Parameters
        ----------
        summaries : list[TConditionSummary]
            Condition summaries.
        comparisons : list
            Pairwise comparison results.
        anova : ANOVASummary or list or None
            ANOVA result(s).
        ranking : list[str]
            Ranked condition labels.
        effective_control : str or None
            Effective control label.
        excluded_conditions : list[ConditionConfig]
            Conditions that were excluded.

        Returns
        -------
        TResult
            Complete comparison result.
        """
        ...

    @abstractmethod
    def _get_replicate_values(self, summary: TConditionSummary) -> list[float]:
        """Extract per-replicate values for statistical tests.

        Parameters
        ----------
        summary : TConditionSummary
            Condition summary.

        Returns
        -------
        list[float]
            Per-replicate values of the primary metric.
        """
        ...

    @abstractmethod
    def _get_mean_value(self, summary: TConditionSummary) -> float:
        """Get the mean value of the primary metric.

        Parameters
        ----------
        summary : TConditionSummary
            Condition summary.

        Returns
        -------
        float
            Mean value.
        """
        ...

    @abstractmethod
    def _interpret_direction(self, percent_change: float) -> str:
        """Interpret the direction of change for this metric.

        Parameters
        ----------
        percent_change : float
            Percent change from control to treatment.

        Returns
        -------
        str
            Direction interpretation (e.g., "stabilizing", "improving").
        """
        ...

    @abstractmethod
    def _rank_summaries(self, summaries: list[TConditionSummary]) -> list[TConditionSummary]:
        """Sort summaries by the primary metric.

        Parameters
        ----------
        summaries : list[TConditionSummary]
            Condition summaries to rank.

        Returns
        -------
        list[TConditionSummary]
            Sorted summaries (best first).
        """
        ...

    # ========================================================================
    # Hook Methods (can be overridden by subclasses)
    # ========================================================================

    def _filter_conditions(
        self,
    ) -> tuple[list["ConditionConfig"], list["ConditionConfig"]]:
        """Filter conditions before analysis.

        Override this to exclude certain conditions (e.g., no-polymer conditions
        for contacts analysis).

        Returns
        -------
        tuple[list[ConditionConfig], list[ConditionConfig]]
            (valid_conditions, excluded_conditions)
        """
        return self.config.conditions, []

    def _get_effective_control(
        self,
        excluded_conditions: list["ConditionConfig"],
    ) -> str | None:
        """Determine the effective control label.

        If the configured control was excluded, returns None.

        Parameters
        ----------
        excluded_conditions : list[ConditionConfig]
            Conditions that were excluded.

        Returns
        -------
        str or None
            Effective control label.
        """
        if not self.config.control:
            return None

        excluded_labels = {c.label for c in excluded_conditions}
        if self.config.control in excluded_labels:
            logger.warning(
                f"Control '{self.config.control}' was excluded. Comparisons will be pairwise."
            )
            return None

        return self.config.control

    def _use_rmsf_mode_for_cohens_d(self) -> bool:
        """Whether to use RMSF-specific Cohen's d interpretation.

        Override in RMSF comparator to return True.

        Returns
        -------
        bool
            True if negative d should be "stabilizing".
        """
        return False

    # ========================================================================
    # Shared Implementation Methods (DRY)
    # ========================================================================

    def _compute_pairwise_comparisons(
        self,
        summaries: list[TConditionSummary],
        effective_control: str | None,
    ) -> list[PairwiseComparison]:
        """Compute pairwise statistical comparisons.

        If a control is specified, compares all conditions vs control.
        Otherwise, compares all pairs.

        Parameters
        ----------
        summaries : list[TConditionSummary]
            Condition summaries.
        effective_control : str or None
            Control condition label.

        Returns
        -------
        list[PairwiseComparison]
            Pairwise comparison results.
        """
        comparisons = []

        if effective_control:
            # Compare all vs control
            control = next(s for s in summaries if s.label == effective_control)
            treatments = [s for s in summaries if s.label != effective_control]

            for treatment in treatments:
                comp = self._compare_pair(control, treatment)
                comparisons.append(comp)
        else:
            # Compare all pairs
            for i, cond_a in enumerate(summaries):
                for cond_b in summaries[i + 1 :]:
                    comp = self._compare_pair(cond_a, cond_b)
                    comparisons.append(comp)

        return comparisons

    def _compare_pair(
        self,
        cond_a: TConditionSummary,
        cond_b: TConditionSummary,
    ) -> PairwiseComparison:
        """Compare two conditions statistically.

        Override this method for comparators that need custom comparison
        logic (e.g., multiple metrics like contacts).

        Parameters
        ----------
        cond_a : TConditionSummary
            First condition (typically control).
        cond_b : TConditionSummary
            Second condition (typically treatment).

        Returns
        -------
        PairwiseComparison
            Statistical comparison result.
        """
        values_a = self._get_replicate_values(cond_a)
        values_b = self._get_replicate_values(cond_b)
        mean_a = self._get_mean_value(cond_a)
        mean_b = self._get_mean_value(cond_b)

        # T-test
        ttest = independent_ttest(values_a, values_b)

        # Effect size
        effect = cohens_d(values_a, values_b, rmsf_mode=self._use_rmsf_mode_for_cohens_d())

        # Percent change
        pct = percent_change(mean_a, mean_b)

        # Direction interpretation
        direction = self._interpret_direction(pct)

        return PairwiseComparison(
            condition_a=cond_a.label,
            condition_b=cond_b.label,
            metric=self.comparison_type_name(),
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
        summaries: list[TConditionSummary],
    ) -> ANOVASummary:
        """Compute one-way ANOVA across all conditions.

        Override this for comparators that test multiple metrics.

        Parameters
        ----------
        summaries : list[TConditionSummary]
            Condition summaries.

        Returns
        -------
        ANOVASummary
            ANOVA result.
        """
        groups = [self._get_replicate_values(s) for s in summaries]
        result = one_way_anova(*groups)

        return ANOVASummary(
            metric=self.comparison_type_name(),
            f_statistic=result.f_statistic,
            p_value=result.p_value,
            significant=result.significant,
        )

    def _compute_ranking(self, summaries: list[TConditionSummary]) -> list[str]:
        """Compute ranking of conditions.

        Parameters
        ----------
        summaries : list[TConditionSummary]
            Condition summaries.

        Returns
        -------
        list[str]
            Labels in ranked order (best first).
        """
        ranked = self._rank_summaries(summaries)
        return [s.label for s in ranked]
