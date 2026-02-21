"""SASA trajectory analysis module.

Provides per-frame, per-residue solvent-accessible surface area computation
using MDTraj's shrake_rupley algorithm. SASA is computed on protein-only
trajectories to measure intrinsic exposure independent of polymer position.

Key classes
-----------
SASATrajectoryResult
    Per-frame, per-residue SASA with classification helpers.
SASAConfig
    Pydantic configuration model.

Key functions
-------------
compute_trajectory_sasa
    Compute (or load cached) trajectory SASA.
"""

from polyzymd.analysis.sasa.config import SASAConfig
from polyzymd.analysis.sasa.trajectory import SASATrajectoryResult, compute_trajectory_sasa

__all__ = [
    "SASAConfig",
    "SASATrajectoryResult",
    "compute_trajectory_sasa",
]
