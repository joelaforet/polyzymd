"""PolyzyMD Comparison Module.

This module provides tools for comparing analysis results across multiple
simulation conditions (e.g., different polymer compositions, temperatures).

Main components:
- ComparisonConfig: Schema for comparison.yaml configuration files
- RMSFComparator: Compare RMSF across conditions with statistical analysis
- ComparisonResult: Structured results with statistics and rankings
- Formatters: Output formatting for console, markdown, and JSON

Usage:
    # Initialize a comparison project
    polyzymd compare init my_polymer_study

    # Edit comparison.yaml to define conditions

    # Run RMSF comparison
    polyzymd compare rmsf --eq-time 10ns
"""

from polyzymd.compare.config import ComparisonConfig, ConditionConfig
from polyzymd.compare.comparator import RMSFComparator
from polyzymd.compare.results import ComparisonResult
from polyzymd.compare.formatters import (
    format_console_table,
    format_markdown,
    format_result,
    to_json,
)
from polyzymd.compare.plotting import (
    plot_rmsf_comparison,
    plot_percent_change,
    plot_effect_sizes,
    plot_summary_panel,
)

__all__ = [
    "ComparisonConfig",
    "ConditionConfig",
    "RMSFComparator",
    "ComparisonResult",
    "format_console_table",
    "format_markdown",
    "format_result",
    "to_json",
    # Plotting functions
    "plot_rmsf_comparison",
    "plot_percent_change",
    "plot_effect_sizes",
    "plot_summary_panel",
]
