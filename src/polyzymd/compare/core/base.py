"""Base classes for comparison results.

This module provides abstract base classes for comparison result models
used by the comparison analysis system.

Classes
-------
PairwiseComparison
    Shared model for statistical comparison between two conditions.
ANOVASummary
    Shared model for ANOVA results.
BaseConditionSummary
    Abstract base for condition-level summary statistics.
BaseComparisonResult
    Abstract base for complete comparison results with save/load.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Any, ClassVar, Generic, Self, TypeVar

from pydantic import BaseModel, ConfigDict

from polyzymd import __version__

# ============================================================================
# Shared Result Models (DRY - used by all analysis plugins)
# ============================================================================


class PairwiseComparison(BaseModel):
    """Statistical comparison between two conditions.

    This is the standard pairwise comparison result used across all
    analysis plugins. For plugins that need additional fields
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

    model_config = ConfigDict(ser_json_inf_nan="strings")

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
        """Return the primary metric value for ranking/comparison."""
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
        The primary metric being compared.
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

    model_config = ConfigDict(ser_json_inf_nan="strings")

    # Class variable - subclasses should override
    comparison_type: ClassVar[str] = "base"

    metric: str
    name: str
    control_label: str | None = None
    # These fields are typed as list[Any] because Pydantic v2 does not
    # resolve Generic type variables for field validation.  Subclasses
    # override these with concrete types (e.g. list[RMSFConditionSummary]).
    conditions: list[Any]
    pairwise_comparisons: list[Any]
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
