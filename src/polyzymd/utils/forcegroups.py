"""
Force group utilities for OpenMM systems.

This module provides functions for assigning and managing force groups
in OpenMM System objects, enabling per-force energy decomposition.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import openmm


def impose_unique_force_groups(system: "openmm.System") -> None:
    """Assign unique force group IDs to each force in the system.

    This enables per-force-type energy decomposition during simulation
    by assigning each force to its own group (0, 1, 2, ...).

    After calling this function, you can get energy contributions from
    individual forces using:

        state = context.getState(getEnergy=True, groups={i})
        energy = state.getPotentialEnergy()

    Args:
        system: OpenMM System object to modify in-place.

    Example:
        >>> from openmm import System
        >>> from polyzymd.utils.forcegroups import impose_unique_force_groups
        >>> system = System()
        >>> # ... add forces to system ...
        >>> impose_unique_force_groups(system)
        >>> # Now each force has a unique group ID
    """
    for i in range(system.getNumForces()):
        system.getForce(i).setForceGroup(i)
