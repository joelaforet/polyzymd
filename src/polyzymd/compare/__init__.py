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

__all__ = [
    "ComparisonConfig",
    "ConditionConfig",
    "RMSFComparator",
    "ComparisonResult",
    "format_console_table",
    "format_markdown",
    "format_result",
    "to_json",
]
