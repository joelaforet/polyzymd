"""RMSF comparison result models (legacy).

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    legacy RMSF comparison result models is now
    :mod:`polyzymd.analyses.rmsf._comparison_results_legacy`.

All symbols are re-exported so that existing
``from polyzymd.compare.results.rmsf_legacy import …`` statements continue
to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.rmsf._comparison_results_legacy import (  # noqa: F401
    ANOVASummary,
    ComparisonResult,
    ConditionSummary,
    PairwiseComparison,
)

__all__ = [
    "ANOVASummary",
    "ComparisonResult",
    "ConditionSummary",
    "PairwiseComparison",
]
