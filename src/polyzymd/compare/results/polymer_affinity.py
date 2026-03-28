"""Result models for polymer affinity score comparison analysis.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    polymer affinity score comparison result models is now
    :mod:`polyzymd.analyses.polymer_affinity._comparison_results`.

All symbols are re-exported so that existing ``from polyzymd.compare.results.polymer_affinity import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.polymer_affinity._comparison_results import (  # noqa: F401
    AffinityScoreConditionSummary,
    AffinityScoreEntry,
    AffinityScorePairwiseEntry,
    PolymerAffinityScoreResult,
    PolymerTypeScore,
)

__all__ = [
    "AffinityScoreConditionSummary",
    "AffinityScoreEntry",
    "AffinityScorePairwiseEntry",
    "PolymerAffinityScoreResult",
    "PolymerTypeScore",
]
