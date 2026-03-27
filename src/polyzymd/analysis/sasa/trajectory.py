"""Backward-compatibility re-export shim.

Canonical location: ``polyzymd.analyses._sasa_trajectory``

This module re-exports ``SASATrajectoryResult`` and ``compute_trajectory_sasa``
from their new home so that existing imports continue to work during the
migration period.  Will be removed in Phase 7.
"""

from polyzymd.analyses._sasa_trajectory import (
    SASATrajectoryResult,
    compute_trajectory_sasa,
)

__all__ = [
    "SASATrajectoryResult",
    "compute_trajectory_sasa",
]
