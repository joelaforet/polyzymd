"""Backward-compatibility re-export shim.

Canonical location: ``polyzymd.analyses._exposure_dynamics``

This module re-exports exposure dynamics classes and functions from their new
home so that existing imports continue to work during the migration period.
Will be removed in Phase 7.
"""

from polyzymd.analyses._exposure_dynamics import (
    ExposureDynamicsResult,
    ResidueExposureSummary,
    analyze_exposure_dynamics,
)

__all__ = [
    "ResidueExposureSummary",
    "ExposureDynamicsResult",
    "analyze_exposure_dynamics",
]
