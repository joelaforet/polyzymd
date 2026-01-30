"""Trajectory analysis and post-processing tools."""

from polyzymd.analysis.make_whole import (
    find_production_trajectories,
    make_whole_trajectory,
)
from polyzymd.analysis.unwrap import (
    ClusterAroundProtein,
    box_dimensions_to_matrix,
    cluster_around_point,
    compute_residue_coms,
    make_molecules_whole,
    minimum_image_shift,
)

__all__ = [
    # make_whole.py - high-level API
    "find_production_trajectories",
    "make_whole_trajectory",
    # unwrap.py - transformation class
    "ClusterAroundProtein",
    # unwrap.py - utility functions (for advanced users)
    "box_dimensions_to_matrix",
    "minimum_image_shift",
    "make_molecules_whole",
    "compute_residue_coms",
    "cluster_around_point",
]
