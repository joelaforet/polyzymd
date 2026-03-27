"""Backward-compatibility shim — canonical location is ``polyzymd.analyses._results_secondary_structure``.

All public names are re-exported so that existing
``from polyzymd.analysis.secondary_structure.results import …`` statements
continue to work during the migration period.
"""

from polyzymd.analyses._results_secondary_structure import (  # noqa: F401
    SS_CHAR_TO_INT,
    SS_INT_TO_CHAR,
    SS_LABELS,
    SecondaryStructureAggregatedResult,
    SecondaryStructureResult,
)

__all__ = [
    "SecondaryStructureResult",
    "SecondaryStructureAggregatedResult",
    "SS_CHAR_TO_INT",
    "SS_INT_TO_CHAR",
    "SS_LABELS",
]
