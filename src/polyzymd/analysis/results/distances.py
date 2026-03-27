"""Backward-compatibility shim — canonical location is ``polyzymd.analyses._results_distances``.

All public names are re-exported so that existing ``from polyzymd.analysis.results.distances
import …`` statements continue to work during the migration period.
"""

from polyzymd.analyses._results_distances import (  # noqa: F401
    DistanceAggregatedResult,
    DistancePairAggregatedResult,
    DistancePairResult,
    DistanceResult,
)

__all__ = [
    "DistancePairResult",
    "DistanceResult",
    "DistancePairAggregatedResult",
    "DistanceAggregatedResult",
]
