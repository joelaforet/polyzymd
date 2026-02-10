"""Analysis result models.

This module provides Pydantic models for serializing and storing
analysis results with proper metadata for reproducibility.

Classes
-------
BaseAnalysisResult
    Base class for all result types
RMSFResult
    Single replicate RMSF results
RMSFAggregatedResult
    Multi-replicate aggregated RMSF results
DistanceResult
    Single replicate distance analysis
DistanceAggregatedResult
    Multi-replicate aggregated distance results
"""

from polyzymd.analysis.results.base import (
    AggregatedResultMixin,
    BaseAnalysisResult,
    get_polyzymd_version,
)
from polyzymd.analysis.results.distances import (
    DistanceAggregatedResult,
    DistancePairAggregatedResult,
    DistancePairResult,
    DistanceResult,
)
from polyzymd.analysis.results.rmsf import (
    RMSFAggregatedResult,
    RMSFResult,
)

__all__ = [
    # Base
    "BaseAnalysisResult",
    "AggregatedResultMixin",
    "get_polyzymd_version",
    # RMSF
    "RMSFResult",
    "RMSFAggregatedResult",
    # Distances
    "DistancePairResult",
    "DistanceResult",
    "DistancePairAggregatedResult",
    "DistanceAggregatedResult",
]
