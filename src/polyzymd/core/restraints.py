"""
Restraint definitions and application for OpenMM systems.

This module provides classes for defining and applying various types
of restraints (flat-bottom, harmonic, etc.) to OpenMM simulations.
"""

from __future__ import annotations

import ast
import logging
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, Any, Dict, List, Optional

if TYPE_CHECKING:
    from openmm import CustomBondForce, HarmonicBondForce, System
    from openmm.app import Topology as OpenMMTopology
    from openmm.unit import Quantity

logger = logging.getLogger(__name__)


def _distance_in_angstroms(value: float) -> Quantity:
    """Create an OpenMM distance quantity in Angstroms.

    Parameters
    ----------
    value : float
        Distance magnitude in Angstroms.

    Returns
    -------
    Quantity
        OpenMM quantity with Angstrom units.
    """
    from openmm.unit import angstrom

    return value * angstrom


def _force_constant_in_kj_per_mol_nm2(value: float) -> Quantity:
    """Create an OpenMM force constant quantity.

    Parameters
    ----------
    value : float
        Force constant magnitude in kJ/mol/nm^2.

    Returns
    -------
    Quantity
        OpenMM quantity with kJ/mol/nm^2 units.
    """
    from openmm.unit import kilojoule_per_mole, nanometer

    return value * kilojoule_per_mole / nanometer**2


def _quantity_in_nanometers(quantity: Quantity) -> float:
    """Convert an OpenMM distance quantity to nanometers.

    Parameters
    ----------
    quantity : Quantity
        Distance quantity to convert.

    Returns
    -------
    float
        Distance magnitude in nanometers.
    """
    from openmm.unit import nanometer

    return quantity.value_in_unit(nanometer)


def _force_constant_value_in_openmm_units(quantity: Quantity) -> float:
    """Convert an OpenMM force constant to kJ/mol/nm^2.

    Parameters
    ----------
    quantity : Quantity
        Force constant quantity to convert.

    Returns
    -------
    float
        Force constant magnitude in kJ/mol/nm^2.
    """
    from openmm.unit import kilojoule_per_mole, nanometer

    return quantity.value_in_unit(kilojoule_per_mole / nanometer**2)


class RestraintType(str, Enum):
    """Types of restraints that can be applied."""

    FLAT_BOTTOM = "flat_bottom"
    HARMONIC = "harmonic"
    UPPER_WALL = "upper_wall"
    LOWER_WALL = "lower_wall"


@dataclass
class AtomSelection:
    """An atom selection for a restraint, resolved with MDTraj on the built topology.

    ``resid``, ``chain`` and ``pdbindex`` keep their MDAnalysis meaning; see
    :func:`_parse_selection`.

    Attributes:
        selection: Selection string, such as ``"protein and resid 77 and name OG"``
        description: Human-readable description of what this selects

    Example:
        >>> sel = AtomSelection("resid 77 and name OG", "Catalytic serine oxygen")
        >>> indices = sel.resolve(topology)
    """

    selection: str
    description: Optional[str] = None

    def resolve(self, topology: OpenMMTopology) -> List[int]:
        """Resolve the selection to atom indices.

        Args:
            topology: OpenMM Topology object

        Returns:
            List of atom indices matching the selection

        Raises:
            ValueError: If selection syntax is invalid or no atoms match
        """
        return _parse_selection(self.selection, topology)


def _parse_selection(selection: str, topology: OpenMMTopology) -> List[int]:
    """Return the sorted 0-based indices of the atoms ``selection`` picks in ``topology``.

    The selection is an MDTraj selection with three MDAnalysis spellings
    translated first: ``resid N`` is the residue number of the built PDB
    (MDTraj ``resSeq N``), ``chain X`` or ``chainid X`` is the chain letter,
    and ``pdbindex N`` is the PDB atom serial, counted from 1 (``index
    N-1``). ``index N`` counts from 0, as OpenMM does. ``protein``, ``not``,
    ``element``, ranges (``resid 70 to 80``) and parentheses work as in MDTraj.
    The words ``and``, ``or``, ``not`` and ``to`` may be written in any case.

    Raises:
        ValueError: If the selection cannot be parsed, gives a word to a
            numeric keyword such as ``index`` or ``resid``, or matches no atom.
        ImportError: If MDTraj is not installed.
    """
    try:
        import mdtraj
        from mdtraj.core.selection import parse_selection
    except ImportError as error:
        raise ImportError(
            "Restraint selections need MDTraj. Use the pixi 'build' environment "
            "or install the 'analysis' extra: pip install 'polyzymd[analysis]'."
        ) from error

    md_topology = mdtraj.Topology.from_openmm(topology)

    def chain(match: re.Match) -> str:
        found = [str(c.index) for c in md_topology.chains if c.chain_id == match.group(1)]
        return "(" + " or ".join(f"chainid {i}" for i in found) + ")" if found else "none"

    translated = re.sub(r"\bpdbindex\s+(\d+)", lambda m: f"index {int(m.group(1)) - 1}", selection)
    translated = re.sub(r"\bchain(?:id)?\s+([^\s()]+)", chain, translated)
    translated = re.sub(r"\bresid\b", "resSeq", translated)
    translated = re.sub(
        r"\b(and|or|not|to)\b", lambda m: m.group(1).lower(), translated, flags=re.IGNORECASE
    )
    try:
        astnode = parse_selection(translated).astnode
        indices = md_topology.select(translated)
    except Exception as error:  # MDTraj raises several types for a bad selection
        raise ValueError(f"Cannot parse selection {selection!r}: {error}") from error
    # MDTraj reads any word after a keyword as one more value, so a misspelled
    # operator such as "index 4 x" would otherwise still select atom 4.
    for node in ast.walk(astnode):
        if isinstance(node, ast.Compare) and any(
            isinstance(n, ast.Attribute) and n.attr in ("index", "resSeq") for n in ast.walk(node)
        ):
            words = [
                n.value
                for n in ast.walk(node)
                if isinstance(n, ast.Constant) and isinstance(n.value, str)
            ]
            if words:
                raise ValueError(
                    f"Cannot parse selection {selection!r}: index, resid, residue, "
                    f"pdbindex and chainid take numbers, not {', '.join(map(repr, words))}"
                )
    if len(indices) == 0:
        raise ValueError(f"No atoms match selection: '{selection}'")
    return sorted(int(i) for i in indices)


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
    distance: Quantity = field(default_factory=lambda: _distance_in_angstroms(3.3))
    force_constant: Quantity = field(
        default_factory=lambda: _force_constant_in_kj_per_mol_nm2(10000)
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
        from openmm import CustomBondForce

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
        from openmm import HarmonicBondForce

        force = HarmonicBondForce()
        # Convert distance to nanometers for OpenMM
        r0_nm = _quantity_in_nanometers(self.distance)
        # Convert force constant to kJ/mol/nm^2
        k_value = _force_constant_value_in_openmm_units(self.force_constant)
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
        from openmm import CustomBondForce

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
        distance = _distance_in_angstroms(distance_value)

        # Parse force constant (default unit: kJ/mol/nm^2)
        k_value = config.get("force_constant", 10000.0)
        force_constant = _force_constant_in_kj_per_mol_nm2(k_value)

        return RestraintDefinition(
            restraint_type=restraint_type,
            name=config.get("name", "unnamed_restraint"),
            atom1=atom1,
            atom2=atom2,
            distance=distance,
            force_constant=force_constant,
            enabled=config.get("enabled", True),
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
