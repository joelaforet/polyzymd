"""
Positional restraints for equilibration protocols.

This module provides classes for applying harmonic positional restraints
to atoms during equilibration, restraining them to their initial coordinates.
This is distinct from distance restraints (in restraints.py) which restrain
the distance between two atoms.

Typical usage in MD equilibration protocols:
- Restrain protein heavy atoms during heating
- Restrain protein backbone while allowing sidechains to relax
- Gradually release restraints on different system components
"""

from __future__ import annotations

import logging
from typing import List, Optional

from openmm import CustomExternalForce, System
from openmm.unit import Quantity, nanometer, kilojoule_per_mole

logger = logging.getLogger(__name__)

# Unit conversion constant
# 1.0 kcal/mol/A^2 = 4.184 kJ/mol/A^2 = 418.4 kJ/mol/nm^2
# But we typically express as kJ/mol/nm^2, so:
# 1.0 kcal/mol/A^2 = 4184.0 kJ/mol/nm^2
KCAL_MOL_ANGSTROM2_TO_KJ_MOL_NM2 = 4184.0


class PositionalRestraintForce:
    """Creates and manages harmonic positional restraints.

    Restrains atoms to their reference positions using a harmonic potential:
        U = 0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)

    where (x0, y0, z0) is the reference position and k is the force constant.

    This uses OpenMM's CustomExternalForce, which applies a force to particles
    based on their absolute coordinates (not relative to other particles).

    The force constant k is stored as a per-particle parameter, allowing
    different atoms to have different restraint strengths. This also avoids
    OpenMM errors when multiple restraint forces are added to the same system
    (which would occur with global parameters).

    Attributes:
        default_force_constant: Default force constant in kJ/mol/nm^2
        particle_count: Number of particles added to the restraint

    Example:
        >>> restraint = PositionalRestraintForce(4184.0)  # 1 kcal/mol/A^2
        >>> restraint.add_particles_from_positions([0, 1, 2], positions)
        >>> force_idx = system.addForce(restraint.force)
    """

    def __init__(self, default_force_constant: float = 4184.0):
        """Initialize positional restraint force.

        Args:
            default_force_constant: Default force constant in kJ/mol/nm^2.
                                   Common value: 4184.0 = 1.0 kcal/mol/A^2
                                   Can be overridden per-particle when adding atoms.
        """
        self._default_force_constant = default_force_constant
        self._particle_count = 0

        # Create CustomExternalForce with harmonic potential
        # All parameters are per-particle to avoid global parameter conflicts
        expression = "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
        self._force = CustomExternalForce(expression)

        # Add per-particle parameters (k is per-particle, not global)
        self._force.addPerParticleParameter("k")
        self._force.addPerParticleParameter("x0")
        self._force.addPerParticleParameter("y0")
        self._force.addPerParticleParameter("z0")

    @property
    def force(self) -> CustomExternalForce:
        """Get the underlying OpenMM CustomExternalForce object."""
        return self._force

    @property
    def force_constant(self) -> float:
        """Get the default force constant in kJ/mol/nm^2."""
        return self._default_force_constant

    @property
    def particle_count(self) -> int:
        """Get the number of particles in this restraint."""
        return self._particle_count

    def add_particle(
        self,
        atom_index: int,
        x0: float,
        y0: float,
        z0: float,
        force_constant: Optional[float] = None,
    ) -> int:
        """Add a particle to restrain at a specific position.

        Args:
            atom_index: OpenMM atom index (0-indexed)
            x0: Reference x coordinate in nanometers
            y0: Reference y coordinate in nanometers
            z0: Reference z coordinate in nanometers
            force_constant: Force constant in kJ/mol/nm^2 (uses default if None)

        Returns:
            Index of the particle in this force
        """
        k = force_constant if force_constant is not None else self._default_force_constant
        particle_idx = self._force.addParticle(atom_index, [k, x0, y0, z0])
        self._particle_count += 1
        return particle_idx

    def add_particles_from_positions(
        self,
        atom_indices: List[int],
        positions: Quantity,
        force_constant: Optional[float] = None,
    ) -> int:
        """Add multiple particles, restraining each to its current position.

        Args:
            atom_indices: List of OpenMM atom indices (0-indexed)
            positions: OpenMM Quantity with positions for all atoms in the system.
                      Must be indexable by atom index.
            force_constant: Force constant in kJ/mol/nm^2 (uses default if None)

        Returns:
            Number of particles added
        """
        k = force_constant if force_constant is not None else self._default_force_constant
        count = 0
        for idx in atom_indices:
            pos = positions[idx]
            # Convert to nanometers
            x0 = pos[0].value_in_unit(nanometer)
            y0 = pos[1].value_in_unit(nanometer)
            z0 = pos[2].value_in_unit(nanometer)
            self.add_particle(idx, x0, y0, z0, force_constant=k)
            count += 1

        logger.debug(f"Added {count} particles to positional restraint (k={k:.1f} kJ/mol/nm^2)")
        return count


def create_position_restraints(
    atom_indices: List[int],
    positions: Quantity,
    force_constant: float,
) -> CustomExternalForce:
    """Convenience function to create a positional restraint force.

    Args:
        atom_indices: List of atom indices to restrain
        positions: Current positions of all atoms
        force_constant: Force constant in kJ/mol/nm^2

    Returns:
        OpenMM CustomExternalForce ready to add to a System
    """
    restraint = PositionalRestraintForce(force_constant)
    restraint.add_particles_from_positions(atom_indices, positions)
    return restraint.force


def add_position_restraints_to_system(
    system: System,
    atom_indices: List[int],
    positions: Quantity,
    force_constant: float,
    force_group: int = 0,
) -> int:
    """Add positional restraints to an OpenMM System.

    Args:
        system: OpenMM System to modify
        atom_indices: List of atom indices to restrain
        positions: Current positions of all atoms
        force_constant: Force constant in kJ/mol/nm^2
        force_group: Force group for the restraint (default 0)

    Returns:
        Index of the added force in the System
    """
    if not atom_indices:
        logger.warning("No atom indices provided for position restraints")
        return -1

    force = create_position_restraints(atom_indices, positions, force_constant)
    force.setForceGroup(force_group)
    force_idx = system.addForce(force)

    logger.info(
        f"Added positional restraints to {len(atom_indices)} atoms "
        f"(k={force_constant:.1f} kJ/mol/nm^2, force_idx={force_idx})"
    )

    return force_idx


def remove_position_restraints_from_system(
    system: System,
    force_indices: List[int],
) -> int:
    """Remove positional restraint forces from an OpenMM System.

    Forces must be removed in reverse order to maintain valid indices.

    Args:
        system: OpenMM System to modify
        force_indices: List of force indices to remove

    Returns:
        Number of forces removed
    """
    # Sort in reverse order so we remove from the end first
    # (removing forces changes indices of subsequent forces)
    sorted_indices = sorted(force_indices, reverse=True)

    removed = 0
    for idx in sorted_indices:
        if 0 <= idx < system.getNumForces():
            system.removeForce(idx)
            removed += 1
            logger.debug(f"Removed positional restraint force at index {idx}")
        else:
            logger.warning(f"Invalid force index {idx}, skipping removal")

    if removed > 0:
        logger.info(f"Removed {removed} positional restraint forces from system")

    return removed
