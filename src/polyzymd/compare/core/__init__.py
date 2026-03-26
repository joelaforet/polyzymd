"""Core infrastructure for comparison analysis.

This module provides base classes for comparison results.

Classes
-------
ANOVASummary
    Shared model for ANOVA results.
PairwiseComparison
    Shared model for pairwise statistical comparisons.
BaseConditionSummary
    Abstract base class for condition summaries.
BaseComparisonResult
    Abstract base class for comparison results.
"""

from polyzymd.compare.core.base import (
    ANOVASummary,
    BaseComparisonResult,
    BaseConditionSummary,
    PairwiseComparison,
)

__all__ = [
    "ANOVASummary",
    "BaseComparisonResult",
    "BaseConditionSummary",
    "PairwiseComparison",
]
