"""Result models for distance comparison analysis.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    distance comparison result models is now
    :mod:`polyzymd.analyses.distances._comparison_results`.

All symbols are re-exported so that existing ``from polyzymd.compare.results.distances import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.distances._comparison_results import (  # noqa: F401
    DistanceComparisonResult,
    DistanceConditionSummary,
    DistancePairANOVA,
    DistancePairSummary,
    DistancePairwiseComparison,
)

__all__ = [
    "DistanceComparisonResult",
    "DistanceConditionSummary",
    "DistancePairANOVA",
    "DistancePairSummary",
    "DistancePairwiseComparison",
]
