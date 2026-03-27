"""Backward-compatibility shim — canonical location is ``polyzymd.analyses._results_rmsf``.

All public names are re-exported so that existing ``from polyzymd.analysis.results.rmsf
import …`` statements continue to work during the migration period.
"""

from polyzymd.analyses._results_rmsf import (  # noqa: F401
    RMSFAggregatedResult,
    RMSFResult,
)

__all__ = [
    "RMSFResult",
    "RMSFAggregatedResult",
]
