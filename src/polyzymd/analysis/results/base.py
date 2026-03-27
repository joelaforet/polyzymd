"""Backward-compatibility shim — canonical location is ``polyzymd.analyses._results_base``.

All public names are re-exported so that existing ``from polyzymd.analysis.results.base
import …`` statements continue to work during the migration period.
"""

from polyzymd.analyses._results_base import (  # noqa: F401
    AggregatedResultMixin,
    BaseAnalysisResult,
    get_polyzymd_version,
)

__all__ = [
    "BaseAnalysisResult",
    "AggregatedResultMixin",
    "get_polyzymd_version",
]
