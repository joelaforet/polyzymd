"""
Restraint definitions and application for OpenMM systems.

This module provides classes for defining and applying various types
of restraints (flat-bottom, harmonic, etc.) to OpenMM simulations.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union

from openmm import CustomBondForce, HarmonicBondForce, System
from openmm.app import Topology as OpenMMTopology
from openmm.unit import Quantity, angstrom, kilojoule_per_mole, nanometer

logger = logging.getLogger(__name__)


class RestraintType(str, Enum):
    """Types of restraints that can be applied."""

    FLAT_BOTTOM = "flat_bottom"
    HARMONIC = "harmonic"
    UPPER_WALL = "upper_wall"
    LOWER_WALL = "lower_wall"


@dataclass
class AtomSelection:
    """Represents a selection of atoms using MDAnalysis-style syntax.

    This class provides a flexible way to specify atoms for restraints
    using selection strings that are compatible with MDAnalysis.

    Attributes:
        selection: MDAnalysis-compatible selection string
        description: Human-readable description of what this selects

    Example:
        >>> sel = AtomSelection("resid 77 and name OG", "Catalytic serine oxygen")
        >>> indices = sel.resolve(topology)
    """

    selection: str
    description: Optional[str] = None

    def resolve(self, topology: OpenMMTopology) -> List[int]:
        """Resolve the selection to atom indices.

        For OpenMM topologies, we parse the selection string and
        find matching atoms. Supports basic MDAnalysis-style selections.

        Args:
            topology: OpenMM Topology object

        Returns:
            List of atom indices matching the selection

        Raises:
            ValueError: If selection syntax is invalid or no atoms match
        """
        return _parse_selection(self.selection, topology)


def _parse_selection(selection: str, topology: OpenMMTopology) -> List[int]:
    """Parse an MDAnalysis-style selection string for OpenMM topology.

    Supports a subset of MDAnalysis selection syntax:
    - resid N: Select atoms in residue with ID N (1-indexed, like PDB)
    - resname XXX: Select atoms in residue with name XXX
    - name XXX: Select atoms with name XXX
    - index N: Select atom with OpenMM index N (0-indexed)
    - pdbindex N: Select atom with PDB serial number N (1-indexed, auto-converts)
    - and: Intersection of selections
    - or: Union of selections

    Index conventions:
    - `resid` uses 1-indexed residue numbers (matches PDB/PyMOL display)
    - `index` uses 0-indexed atom indices (matches OpenMM internal indexing)
    - `pdbindex` uses 1-indexed atom serial (matches PDB ATOM column, PyMOL display)

    Example: If PyMOL shows atom serial 2740 and residue 77:
    - Use "pdbindex 2740" or "index 2739" to select that atom
    - Use "resid 77" to select all atoms in that residue

    Args:
        selection: Selection string
        topology: OpenMM Topology

    Returns:
        List of matching atom indices
    """
    # Tokenize the selection
    tokens = selection.lower().replace("(", " ( ").replace(")", " ) ").split()

    # Build list of all atoms with their properties
    atoms_data = []
    for atom in topology.atoms():
        atoms_data.append(
            {
                "index": atom.index,
                "name": atom.name.lower() if atom.name else "",
                "resname": atom.residue.name.lower() if atom.residue else "",
                "resid": atom.residue.index,  # 0-indexed internally
                "chain": atom.residue.chain.id if atom.residue and atom.residue.chain else "",
            }
        )

    def evaluate_simple(tokens: List[str], start: int) -> Tuple[set, int]:
        """Evaluate a simple selection (keyword value)."""
        if start >= len(tokens):
            return set(), start

        keyword = tokens[start]

        if keyword == "(":
            # Recurse into parentheses
            result, end = evaluate_or(tokens, start + 1)
            if end < len(tokens) and tokens[end] == ")":
                return result, end + 1
            return result, end

        if start + 1 >= len(tokens):
            raise ValueError(f"Missing value after keyword '{keyword}'")

        value = tokens[start + 1]

        matching = set()

        if keyword == "resid":
            # resid is 1-indexed in selection (matches PDB), 0-indexed internally
            target_resid = int(value) - 1
            for atom in atoms_data:
                if atom["resid"] == target_resid:
                    matching.add(atom["index"])

        elif keyword == "resname":
            for atom in atoms_data:
                if atom["resname"] == value.lower():
                    matching.add(atom["index"])

        elif keyword == "name":
            for atom in atoms_data:
                if atom["name"] == value.lower():
                    matching.add(atom["index"])

        elif keyword == "index":
            # index is 0-indexed (matches OpenMM internal indexing)
            target_idx = int(value)
            if 0 <= target_idx < len(atoms_data):
                matching.add(target_idx)

        elif keyword == "pdbindex":
            # pdbindex is 1-indexed (matches PDB ATOM serial number / PyMOL display)
            # Converts to 0-indexed for internal use
            target_idx = int(value) - 1
            if 0 <= target_idx < len(atoms_data):
                matching.add(target_idx)

        elif keyword == "chainid" or keyword == "chain":
            for atom in atoms_data:
                if atom["chain"].lower() == value.lower():
                    matching.add(atom["index"])
        else:
            raise ValueError(f"Unknown selection keyword: '{keyword}'")

        return matching, start + 2

    def evaluate_and(tokens: List[str], start: int) -> Tuple[set, int]:
        """Evaluate AND expressions."""
        result, pos = evaluate_simple(tokens, start)

        while pos < len(tokens) and tokens[pos] == "and":
            right, pos = evaluate_simple(tokens, pos + 1)
            result = result & right

        return result, pos

    def evaluate_or(tokens: List[str], start: int) -> Tuple[set, int]:
        """Evaluate OR expressions."""
        result, pos = evaluate_and(tokens, start)

        while pos < len(tokens) and tokens[pos] == "or":
            right, pos = evaluate_and(tokens, pos + 1)
            result = result | right

        return result, pos

    if not tokens:
        raise ValueError("Empty selection string")

    result, _ = evaluate_or(tokens, 0)

    if not result:
        raise ValueError(f"No atoms match selection: '{selection}'")

    return sorted(list(result))


@dataclass
class RestraintDefinition:
    """Definition of a single restraint to be applied to a system.

    Attributes:
        restraint_type: Type of restraint (flat_bottom, harmonic, etc.)
        name: Human-readable identifier
        atom1: First atom selection
        atom2: Second atom selection
        distance: Target or threshold distance
        force_constant: Force constant for the restraint
        enabled: Whether this restraint should be applied

    Example:
        >>> restraint = RestraintDefinition(
        ...     restraint_type=RestraintType.FLAT_BOTTOM,
        ...     name="catalytic_serine",
        ...     atom1=AtomSelection("resid 77 and name OG"),
        ...     atom2=AtomSelection("resname LIG and name C12"),
        ...     distance=3.3 * angstrom,
        ...     force_constant=10000 * kilojoule_per_mole / nanometer**2
        ... )
    """

    restraint_type: RestraintType
    name: str
    atom1: AtomSelection
    atom2: AtomSelection
    distance: Quantity = field(default_factory=lambda: 3.3 * angstrom)
    force_constant: Quantity = field(
        default_factory=lambda: 10000 * kilojoule_per_mole / nanometer**2
    )
    enabled: bool = True

    def apply(self, topology: OpenMMTopology, system: System) -> Optional[int]:
        """Apply this restraint to an OpenMM system.

        Args:
            topology: OpenMM Topology for resolving atom selections
            system: OpenMM System to add the force to

        Returns:
            Index of the added force, or None if restraint is disabled
        """
        if not self.enabled:
            logger.info(f"Skipping disabled restraint: {self.name}")
            return None

        # Resolve atom selections
        atom1_indices = self.atom1.resolve(topology)
        atom2_indices = self.atom2.resolve(topology)

        if len(atom1_indices) != 1 or len(atom2_indices) != 1:
            raise ValueError(
                f"Restraint '{self.name}' requires exactly one atom per selection. "
                f"Got {len(atom1_indices)} for atom1, {len(atom2_indices)} for atom2"
            )

        atom1_idx = atom1_indices[0]
        atom2_idx = atom2_indices[0]

        # Create the appropriate force
        if self.restraint_type == RestraintType.FLAT_BOTTOM:
            force = self._create_flat_bottom_force(atom1_idx, atom2_idx)
        elif self.restraint_type == RestraintType.HARMONIC:
            force = self._create_harmonic_force(atom1_idx, atom2_idx)
        elif self.restraint_type == RestraintType.UPPER_WALL:
            force = self._create_upper_wall_force(atom1_idx, atom2_idx)
        elif self.restraint_type == RestraintType.LOWER_WALL:
            force = self._create_lower_wall_force(atom1_idx, atom2_idx)
        else:
            raise ValueError(f"Unknown restraint type: {self.restraint_type}")

        force_idx = system.addForce(force)

        logger.info(
            f"Applied {self.restraint_type.value} restraint '{self.name}' "
            f"between atoms {atom1_idx} and {atom2_idx} "
            f"(r0={self.distance}, k={self.force_constant})"
        )

        return force_idx

    def _create_flat_bottom_force(self, atom1_idx: int, atom2_idx: int) -> CustomBondForce:
        """Create a flat-bottom potential force.

        U(r) = 0 if r < r0
               0.5 * k * (r - r0)^2 if r >= r0
        """
        expression = "step(r - r0) * 0.5 * k * (r - r0)^2"
        force = CustomBondForce(expression)
        force.addGlobalParameter("k", self.force_constant)
        force.addGlobalParameter("r0", self.distance)
        force.addBond(atom1_idx, atom2_idx, [])
        return force

    def _create_harmonic_force(self, atom1_idx: int, atom2_idx: int) -> HarmonicBondForce:
        """Create a harmonic bond force.

        U(r) = 0.5 * k * (r - r0)^2
        """
        force = HarmonicBondForce()
        # Convert distance to nanometers for OpenMM
        r0_nm = self.distance.value_in_unit(nanometer)
        # Convert force constant to kJ/mol/nm^2
        k_value = self.force_constant.value_in_unit(kilojoule_per_mole / nanometer**2)
        force.addBond(atom1_idx, atom2_idx, r0_nm, k_value)
        return force

    def _create_upper_wall_force(self, atom1_idx: int, atom2_idx: int) -> CustomBondForce:
        """Create an upper wall potential (prevent distance exceeding r0).

        U(r) = 0 if r < r0
               0.5 * k * (r - r0)^2 if r >= r0

        (Same as flat bottom)
        """
        return self._create_flat_bottom_force(atom1_idx, atom2_idx)

    def _create_lower_wall_force(self, atom1_idx: int, atom2_idx: int) -> CustomBondForce:
        """Create a lower wall potential (prevent distance below r0).

        U(r) = 0.5 * k * (r0 - r)^2 if r < r0
               0 if r >= r0
        """
        expression = "step(r0 - r) * 0.5 * k * (r0 - r)^2"
        force = CustomBondForce(expression)
        force.addGlobalParameter("k", self.force_constant)
        force.addGlobalParameter("r0", self.distance)
        force.addBond(atom1_idx, atom2_idx, [])
        return force


class RestraintFactory:
    """Factory for creating restraints from configuration.

    This class bridges the configuration schema with the restraint
    implementation, creating RestraintDefinition objects from config.
    """

    @staticmethod
    def from_config(config: Dict[str, Any]) -> RestraintDefinition:
        """Create a RestraintDefinition from a configuration dictionary.

        Args:
            config: Dictionary with restraint configuration

        Returns:
            RestraintDefinition instance
        """
        # Parse restraint type
        type_str = config.get("type", "flat_bottom")
        try:
            restraint_type = RestraintType(type_str)
        except ValueError:
            raise ValueError(f"Unknown restraint type: {type_str}")

        # Parse atom selections
        atom1_config = config.get("atom1", {})
        atom2_config = config.get("atom2", {})

        atom1 = AtomSelection(
            selection=atom1_config.get("selection", ""), description=atom1_config.get("description")
        )
        atom2 = AtomSelection(
            selection=atom2_config.get("selection", ""), description=atom2_config.get("description")
        )

        # Parse distance (default unit: angstrom)
        distance_value = config.get("distance", 3.3)
        distance = distance_value * angstrom

        # Parse force constant (default unit: kJ/mol/nm^2)
        k_value = config.get("force_constant", 10000.0)
        force_constant = k_value * kilojoule_per_mole / nanometer**2

        return RestraintDefinition(
            restraint_type=restraint_type,
            name=config.get("name", "unnamed_restraint"),
            atom1=atom1,
            atom2=atom2,
            distance=distance,
            force_constant=force_constant,
            enabled=config.get("enabled", True),
        )

    @staticmethod
    def create_flat_bottom(
        name: str,
        atom1_selection: str,
        atom2_selection: str,
        distance: float,
        force_constant: float = 10000.0,
    ) -> RestraintDefinition:
        """Convenience method to create a flat-bottom restraint.

        Args:
            name: Restraint identifier
            atom1_selection: Selection string for first atom
            atom2_selection: Selection string for second atom
            distance: Threshold distance in Angstroms
            force_constant: Force constant in kJ/mol/nm^2

        Returns:
            RestraintDefinition for flat-bottom potential
        """
        return RestraintDefinition(
            restraint_type=RestraintType.FLAT_BOTTOM,
            name=name,
            atom1=AtomSelection(atom1_selection),
            atom2=AtomSelection(atom2_selection),
            distance=distance * angstrom,
            force_constant=force_constant * kilojoule_per_mole / nanometer**2,
        )

    @staticmethod
    def create_harmonic(
        name: str,
        atom1_selection: str,
        atom2_selection: str,
        distance: float,
        force_constant: float = 10000.0,
    ) -> RestraintDefinition:
        """Convenience method to create a harmonic restraint.

        Args:
            name: Restraint identifier
            atom1_selection: Selection string for first atom
            atom2_selection: Selection string for second atom
            distance: Equilibrium distance in Angstroms
            force_constant: Force constant in kJ/mol/nm^2

        Returns:
            RestraintDefinition for harmonic potential
        """
        return RestraintDefinition(
            restraint_type=RestraintType.HARMONIC,
            name=name,
            atom1=AtomSelection(atom1_selection),
            atom2=AtomSelection(atom2_selection),
            distance=distance * angstrom,
            force_constant=force_constant * kilojoule_per_mole / nanometer**2,
        )


def apply_restraints(
    restraints: List[RestraintDefinition], topology: OpenMMTopology, system: System
) -> List[int]:
    """Apply multiple restraints to a system.

    Args:
        restraints: List of restraint definitions
        topology: OpenMM Topology for resolving selections
        system: OpenMM System to modify

    Returns:
        List of force indices for the added restraints
    """
    force_indices = []
    for restraint in restraints:
        idx = restraint.apply(topology, system)
        if idx is not None:
            force_indices.append(idx)
    return force_indices
