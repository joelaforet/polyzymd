"""Core analysis infrastructure.

This module provides foundational utilities for trajectory analysis:
- Config hashing for cache validation
- Statistical functions for replicate aggregation
- Autocorrelation analysis for independent sampling
- Trajectory loading from PolyzyMD outputs
- Centroid/representative frame finding for trajectory alignment
- Trajectory alignment utilities (shared across all analysis modules)
- PBC-aware distance calculations
- Registry pattern for extensible analysis types
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
    statistical_inefficiency,
    statistical_inefficiency_multiple,
    n_effective,
    check_statistical_reliability,
)
from polyzymd.analysis.core.loader import (
    TrajectoryInfo,
    TrajectoryLoader,
    parse_time_string,
    convert_time,
    time_to_frame,
)
from polyzymd.analysis.core.centroid import (
    find_centroid_frame,
    find_reference_frame,
    get_reference_mode_description,
)
from polyzymd.analysis.core.selections import (
    SelectionMode,
    ParsedSelection,
    translate_selection,
    parse_selection_string,
    select_atoms,
    get_position,
    get_position_from_selection,
    validate_selection,
    format_selection_for_label,
)
from polyzymd.analysis.core.diagnostics import (
    get_selection_diagnostics,
    get_residue_info,
    get_protein_residue_range,
    format_diagnostic_message,
    validate_equilibration_time,
)
from polyzymd.analysis.core.registry import (
    BaseAnalysisSettings,
    BaseComparisonSettings,
    BaseAnalyzer,
    AnalysisSettingsRegistry,
    ComparisonSettingsRegistry,
    AnalyzerRegistry,
)
from polyzymd.analysis.core.metric_type import (
    MetricType,
    AutocorrelationStrategy,
    get_autocorrelation_strategy,
)
from polyzymd.analysis.core.pbc import (
    minimum_image_distance,
    pairwise_distances_pbc,
)
from polyzymd.analysis.core.alignment import (
    AlignmentConfig,
    ReferenceMode,
    align_trajectory,
    get_alignment_description,
)
from polyzymd.analysis.core.constants import (
    DEFAULT_CONTACT_CUTOFF,
    DEFAULT_DISTANCE_THRESHOLD,
    DEFAULT_SURFACE_EXPOSURE_THRESHOLD,
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
    "statistical_inefficiency",
    "statistical_inefficiency_multiple",
    "n_effective",
    "check_statistical_reliability",
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
    # Selection parsing
    "SelectionMode",
    "ParsedSelection",
    "translate_selection",
    "parse_selection_string",
    "select_atoms",
    "get_position",
    "get_position_from_selection",
    "validate_selection",
    "format_selection_for_label",
    # Diagnostics
    "get_selection_diagnostics",
    "get_residue_info",
    "get_protein_residue_range",
    "format_diagnostic_message",
    "validate_equilibration_time",
    # Registry pattern
    "BaseAnalysisSettings",
    "BaseComparisonSettings",
    "BaseAnalyzer",
    "AnalysisSettingsRegistry",
    "ComparisonSettingsRegistry",
    "AnalyzerRegistry",
    # Metric type for autocorrelation handling
    "MetricType",
    "AutocorrelationStrategy",
    "get_autocorrelation_strategy",
    # PBC utilities
    "minimum_image_distance",
    "pairwise_distances_pbc",
    # Trajectory alignment
    "AlignmentConfig",
    "align_trajectory",
    "get_alignment_description",
    # Shared constants
    "DEFAULT_CONTACT_CUTOFF",
    "DEFAULT_DISTANCE_THRESHOLD",
    "DEFAULT_SURFACE_EXPOSURE_THRESHOLD",
]
