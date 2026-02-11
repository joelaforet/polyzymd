"""Core analysis infrastructure.

This module provides foundational utilities for trajectory analysis:
- Config hashing for cache validation
- Statistical functions for replicate aggregation
- Autocorrelation analysis for independent sampling
- Trajectory loading from PolyzyMD outputs
- Centroid/representative frame finding for trajectory alignment
"""

from polyzymd.analysis.core.config_hash import (
    compute_config_hash,
    validate_config_hash,
)
from polyzymd.analysis.core.statistics import (
    StatResult,
    PerResidueStats,
    compute_sem,
    aggregate_per_residue_stats,
    aggregate_region_stats,
    weighted_mean_with_sem,
)
from polyzymd.analysis.core.autocorrelation import (
    ACFResult,
    CorrelationTimeResult,
    MIN_RECOMMENDED_N_INDEPENDENT,
    compute_acf,
    estimate_correlation_time,
    get_independent_indices,
)
from polyzymd.analysis.core.loader import (
    TrajectoryInfo,
    TrajectoryLoader,
    parse_time_string,
    convert_time,
    time_to_frame,
)
from polyzymd.analysis.core.centroid import (
    ReferenceMode,
    find_centroid_frame,
    find_reference_frame,
    get_reference_mode_description,
)

__all__ = [
    # Config hashing
    "compute_config_hash",
    "validate_config_hash",
    # Statistics
    "StatResult",
    "PerResidueStats",
    "compute_sem",
    "aggregate_per_residue_stats",
    "aggregate_region_stats",
    "weighted_mean_with_sem",
    # Autocorrelation
    "ACFResult",
    "CorrelationTimeResult",
    "MIN_RECOMMENDED_N_INDEPENDENT",
    "compute_acf",
    "estimate_correlation_time",
    "get_independent_indices",
    # Trajectory loading
    "TrajectoryInfo",
    "TrajectoryLoader",
    "parse_time_string",
    "convert_time",
    "time_to_frame",
    # Centroid/reference frame finding
    "ReferenceMode",
    "find_centroid_frame",
    "find_reference_frame",
    "get_reference_mode_description",
]
