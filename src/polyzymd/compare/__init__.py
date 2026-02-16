"""PolyzyMD Comparison Module.

This module provides tools for comparing analysis results across multiple
simulation conditions (e.g., different polymer compositions, temperatures).

Main components:
- ComparisonConfig: Schema for comparison.yaml configuration files
- RMSFComparator: Compare RMSF across conditions with statistical analysis
- TriadComparator: Compare catalytic triad geometry across conditions
- ContactsComparator: Compare polymer-protein contacts across conditions
- ComparisonResult: Structured results with statistics and rankings
- TriadComparisonResult: Structured results for triad comparison
- ContactsComparisonResult: Structured results for contacts comparison
- Formatters: Output formatting for console, markdown, and JSON

Usage:
    # Initialize a comparison project
    polyzymd compare init my_polymer_study

    # Edit comparison.yaml to define conditions

    # Run RMSF comparison
    polyzymd compare rmsf --eq-time 10ns

    # Run triad comparison
    polyzymd compare triad --eq-time 10ns

    # Run contacts comparison
    polyzymd compare contacts --eq-time 10ns
"""

from polyzymd.compare.comparator import RMSFComparator
from polyzymd.compare.config import (
    ComparisonConfig,
    ConditionConfig,
)
from polyzymd.compare.contacts_comparator import ContactsComparator
from polyzymd.compare.contacts_formatters import (
    contacts_to_json,
    format_contacts_console_table,
    format_contacts_markdown,
    format_contacts_result,
)
from polyzymd.compare.formatters import (
    format_console_table,
    format_markdown,
    format_result,
    to_json,
)
from polyzymd.compare.plotting import (
    plot_effect_sizes,
    plot_percent_change,
    plot_rmsf_comparison,
    plot_summary_panel,
)
from polyzymd.compare.results import (
    AggregateComparisonResult,
    ComparisonResult,
    ContactsComparisonResult,
    ContactsConditionSummary,
    ContactsPairwiseComparison,
    TriadComparisonResult,
    TriadConditionSummary,
)
from polyzymd.compare.settings import (
    CatalyticTriadAnalysisSettings,
    CatalyticTriadComparisonSettings,
    ContactsAnalysisSettings,
    ContactsComparisonSettings,
    RMSFAnalysisSettings,
    RMSFComparisonSettings,
)
from polyzymd.compare.triad_comparator import TriadComparator
from polyzymd.compare.triad_formatters import (
    format_triad_console_table,
    format_triad_markdown,
    format_triad_result,
    triad_to_json,
)

__all__ = [
    # Config
    "ComparisonConfig",
    "ConditionConfig",
    # Analysis Settings
    "RMSFAnalysisSettings",
    "RMSFComparisonSettings",
    "CatalyticTriadAnalysisSettings",
    "CatalyticTriadComparisonSettings",
    "ContactsAnalysisSettings",
    "ContactsComparisonSettings",
    # RMSF comparison
    "RMSFComparator",
    "ComparisonResult",
    "format_console_table",
    "format_markdown",
    "format_result",
    "to_json",
    # Triad comparison
    "TriadComparator",
    "TriadComparisonResult",
    "TriadConditionSummary",
    "format_triad_console_table",
    "format_triad_markdown",
    "format_triad_result",
    "triad_to_json",
    # Contacts comparison
    "ContactsComparator",
    "ContactsComparisonResult",
    "ContactsConditionSummary",
    "ContactsPairwiseComparison",
    "AggregateComparisonResult",
    "format_contacts_console_table",
    "format_contacts_markdown",
    "format_contacts_result",
    "contacts_to_json",
    # Plotting functions
    "plot_rmsf_comparison",
    "plot_percent_change",
    "plot_effect_sizes",
    "plot_summary_panel",
]
