"""Backward compatibility shim for comparison base classes.

Canonical comparison base classes now live in ``polyzymd.analyses.base``.
This module re-exports legacy names for older imports.
"""

from polyzymd.analyses.base import (
    ANOVAResult,
    BaseComparisonResult,
    BaseConditionSummary,
    PairwiseResult,
)

# Legacy aliases
PairwiseComparison = PairwiseResult
ANOVASummary = ANOVAResult

__all__ = [
    "ANOVAResult",
    "ANOVASummary",
    "BaseComparisonResult",
    "BaseConditionSummary",
    "PairwiseComparison",
    "PairwiseResult",
]
