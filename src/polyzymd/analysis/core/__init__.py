"""Core analysis infrastructure."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_MODULE_EXPORTS: dict[str, list[str]] = {
    "polyzymd.analysis.core.config_hash": ["compute_config_hash", "validate_config_hash"],
    "polyzymd.analysis.core.statistics": [
        "StatResult",
        "PerResidueStats",
        "compute_sem",
        "aggregate_per_residue_stats",
        "aggregate_region_stats",
        "weighted_mean_with_sem",
    ],
    "polyzymd.analysis.core.autocorrelation": [
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
    ],
    "polyzymd.analysis.core.loader": [
        "TrajectoryInfo",
        "TrajectoryLoader",
        "parse_time_string",
        "convert_time",
        "time_to_frame",
    ],
    "polyzymd.analysis.core.centroid": [
        "find_centroid_frame",
        "find_reference_frame",
        "get_reference_mode_description",
    ],
    "polyzymd.analysis.core.selections": [
        "SelectionMode",
        "ParsedSelection",
        "translate_selection",
        "parse_selection_string",
        "select_atoms",
        "get_position",
        "get_position_from_selection",
        "validate_selection",
        "format_selection_for_label",
    ],
    "polyzymd.analysis.core.diagnostics": [
        "get_selection_diagnostics",
        "get_residue_info",
        "get_protein_residue_range",
        "format_diagnostic_message",
        "validate_equilibration_time",
    ],
    "polyzymd.analysis.core.registry": [
        "BaseAnalysisSettings",
        "BaseComparisonSettings",
        "BaseAnalyzer",
        "BasePlotSettings",
        "AnalysisSettingsRegistry",
        "ComparisonSettingsRegistry",
        "AnalyzerRegistry",
        "PlotSettingsRegistry",
    ],
    "polyzymd.analysis.core.metric_type": [
        "MetricType",
        "AutocorrelationStrategy",
        "get_autocorrelation_strategy",
    ],
    "polyzymd.analysis.core.pbc": ["minimum_image_distance", "pairwise_distances_pbc"],
    "polyzymd.analysis.core.alignment": [
        "AlignmentConfig",
        "ReferenceMode",
        "align_trajectory",
        "get_alignment_description",
    ],
    "polyzymd.analysis.core.constants": [
        "DEFAULT_CONTACT_CUTOFF",
        "DEFAULT_DISTANCE_THRESHOLD",
        "DEFAULT_SURFACE_EXPOSURE_THRESHOLD",
    ],
}

_EXPORTS = {
    name: (module_name, name) for module_name, names in _MODULE_EXPORTS.items() for name in names
}

__all__ = [
    "compute_config_hash",
    "validate_config_hash",
    "StatResult",
    "PerResidueStats",
    "compute_sem",
    "aggregate_per_residue_stats",
    "aggregate_region_stats",
    "weighted_mean_with_sem",
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
    "TrajectoryInfo",
    "TrajectoryLoader",
    "parse_time_string",
    "convert_time",
    "time_to_frame",
    "ReferenceMode",
    "find_centroid_frame",
    "find_reference_frame",
    "get_reference_mode_description",
    "SelectionMode",
    "ParsedSelection",
    "translate_selection",
    "parse_selection_string",
    "select_atoms",
    "get_position",
    "get_position_from_selection",
    "validate_selection",
    "format_selection_for_label",
    "get_selection_diagnostics",
    "get_residue_info",
    "get_protein_residue_range",
    "format_diagnostic_message",
    "validate_equilibration_time",
    "BaseAnalysisSettings",
    "BaseComparisonSettings",
    "BaseAnalyzer",
    "BasePlotSettings",
    "AnalysisSettingsRegistry",
    "ComparisonSettingsRegistry",
    "AnalyzerRegistry",
    "PlotSettingsRegistry",
    "MetricType",
    "AutocorrelationStrategy",
    "get_autocorrelation_strategy",
    "minimum_image_distance",
    "pairwise_distances_pbc",
    "AlignmentConfig",
    "align_trajectory",
    "get_alignment_description",
    "DEFAULT_CONTACT_CUTOFF",
    "DEFAULT_DISTANCE_THRESHOLD",
    "DEFAULT_SURFACE_EXPOSURE_THRESHOLD",
]


def __getattr__(name: str) -> Any:
    """Lazily import analysis-core exports to avoid heavy optional deps."""
    try:
        module_name, attr_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc

    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__() -> list[str]:
    """Return module attributes for tab completion and introspection."""
    return sorted(set(globals()) | set(__all__))
