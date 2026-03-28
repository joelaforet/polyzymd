"""Result models for catalytic triad comparison analysis.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    catalytic triad comparison result models is now
    :mod:`polyzymd.analyses.catalytic_triad._comparison_results`.

All symbols are re-exported so that existing ``from polyzymd.compare.results.triad import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.catalytic_triad._comparison_results import (  # noqa: F401
    TriadANOVASummary,
    TriadComparisonResult,
    TriadConditionSummary,
    TriadPairSummary,
    TriadPairwiseComparison,
)

__all__ = [
    "TriadANOVASummary",
    "TriadComparisonResult",
    "TriadConditionSummary",
    "TriadPairSummary",
    "TriadPairwiseComparison",
]
