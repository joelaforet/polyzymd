"""Backward-compatibility shim — canonical location is ``polyzymd.analyses._results_triad``.

All public names are re-exported so that existing ``from polyzymd.analysis.results.triad
import …`` statements continue to work during the migration period.
"""

from polyzymd.analyses._results_triad import (  # noqa: F401
    TriadAggregatedResult,
    TriadResult,
)

__all__ = [
    "TriadResult",
    "TriadAggregatedResult",
]
