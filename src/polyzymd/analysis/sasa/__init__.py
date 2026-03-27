"""Backward-compatibility re-export shim for SASA analysis module.

Canonical locations are now under ``polyzymd.analyses._sasa_*``.
This package re-exports all public names so that existing
``from polyzymd.analysis.sasa import ...`` statements continue to work.
Will be removed in Phase 7.
"""

from polyzymd.analyses._sasa_config import SASAConfig
from polyzymd.analyses._sasa_trajectory import SASATrajectoryResult, compute_trajectory_sasa

__all__ = [
    "SASAConfig",
    "SASATrajectoryResult",
    "compute_trajectory_sasa",
]
