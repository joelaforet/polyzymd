"""Result models for binding free energy comparison analysis.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    binding free energy comparison result models is now
    :mod:`polyzymd.analyses.binding_free_energy._comparison_results`.

All symbols are re-exported so that existing ``from polyzymd.compare.results.binding_free_energy import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.binding_free_energy._comparison_results import (  # noqa: F401
    BindingFreeEnergyResult,
    FreeEnergyConditionSummary,
    FreeEnergyEntry,
    FreeEnergyPairwiseEntry,
)

__all__ = [
    "BindingFreeEnergyResult",
    "FreeEnergyConditionSummary",
    "FreeEnergyEntry",
    "FreeEnergyPairwiseEntry",
]
