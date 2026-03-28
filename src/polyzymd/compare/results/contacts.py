"""Result models for polymer-protein contacts comparison analysis.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    contacts comparison result models is now
    :mod:`polyzymd.analyses.contacts._comparison_results`.

All symbols are re-exported so that existing ``from polyzymd.compare.results.contacts import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.contacts._comparison_results import (  # noqa: F401
    AggregateComparisonResult,
    BindingPreferenceComparisonEntry,
    BindingPreferenceComparisonSummary,
    ContactsANOVASummary,
    ContactsComparisonResult,
    ContactsConditionSummary,
    ContactsPairwiseComparison,
)

__all__ = [
    "AggregateComparisonResult",
    "BindingPreferenceComparisonEntry",
    "BindingPreferenceComparisonSummary",
    "ContactsANOVASummary",
    "ContactsComparisonResult",
    "ContactsConditionSummary",
    "ContactsPairwiseComparison",
]
