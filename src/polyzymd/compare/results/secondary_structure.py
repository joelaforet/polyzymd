"""Secondary structure condition summary and comparison result models.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    secondary structure comparison result models is now
    :mod:`polyzymd.analyses.secondary_structure._comparison_results`.

All symbols are re-exported so that existing
``from polyzymd.compare.results.secondary_structure import …`` statements
continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.secondary_structure._comparison_results import (  # noqa: F401
    SSComparisonResult,
    SSConditionSummary,
)

__all__ = [
    "SSComparisonResult",
    "SSConditionSummary",
]
