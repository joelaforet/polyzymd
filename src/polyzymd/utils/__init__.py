"""
Utility functions for PolyzyMD.

This module provides internal utilities for:
- Force group assignment for energy decomposition
- OpenFF topology utilities (SDF loading, molecule extraction)
- Molecular charging with ML models (NAGL, Espaloma)
- Box vector calculations (bounding boxes, padding, volumes)
- Unit conversion between OpenFF and OpenMM unit systems
"""

from polyzymd.utils.boxvectors import (
    get_box_volume,
    get_topology_bbox,
    get_topology_bbox_bounds,
    pad_box_vectors_uniform,
)
from polyzymd.utils.charging import (
    AM1BCCCharger,
    EspalomaCharger,
    MoleculeCharger,
    NAGLCharger,
    get_charger,
)
from polyzymd.utils.forcegroups import impose_unique_force_groups
from polyzymd.utils.topology import get_largest_offmol, topology_from_sdf
from polyzymd.utils.units import openff_to_openmm, openmm_to_openff

__all__ = [
    # Force groups
    "impose_unique_force_groups",
    # Topology
    "get_largest_offmol",
    "topology_from_sdf",
    # Charging
    "MoleculeCharger",
    "NAGLCharger",
    "EspalomaCharger",
    "AM1BCCCharger",
    "get_charger",
    # Box vectors
    "get_topology_bbox",
    "get_topology_bbox_bounds",
    "pad_box_vectors_uniform",
    "get_box_volume",
    # Units
    "openff_to_openmm",
    "openmm_to_openff",
]
