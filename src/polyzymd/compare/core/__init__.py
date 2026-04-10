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
    ANOVAResult,
    ANOVASummary,
    BaseComparisonResult,
    BaseConditionSummary,
    PairwiseComparison,
    PairwiseResult,
)

__all__ = [
    "ANOVASummary",
    "ANOVAResult",
    "BaseComparisonResult",
    "BaseConditionSummary",
    "PairwiseComparison",
    "PairwiseResult",
]
