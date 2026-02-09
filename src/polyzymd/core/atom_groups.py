"""
Predefined atom group resolution for equilibration restraints.

This module provides functionality to resolve predefined atom group names
(like 'protein_heavy', 'polymer_heavy', etc.) to atom indices for use
with positional restraints during equilibration.

The atom groups are based on the chain assignment convention used in
SystemBuilder:
- Chain A: Protein/Enzyme
- Chain B: Substrate/Ligand
- Chain C: Polymers
- Chain D+: Solvent (water, ions, co-solvents)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, FrozenSet, List, Optional, Set

from openmm.app import Topology as OpenMMTopology

logger = logging.getLogger(__name__)

# Predefined atom group names
PREDEFINED_GROUPS: FrozenSet[str] = frozenset(
    {
        # Protein/Enzyme groups
        "protein_heavy",
        "protein_backbone",
        "protein_calpha",
        # Ligand/Substrate groups
        "ligand_heavy",
        # Polymer groups
        "polymer_heavy",
        # Solvent groups
        "solvent",
        "water_only",
        "ions_only",
        "cosolvents_only",
    }
)

# Backbone atom names (standard amino acids)
BACKBONE_ATOM_NAMES: FrozenSet[str] = frozenset({"N", "CA", "C", "O"})

# Common water residue names
WATER_RESIDUE_NAMES: FrozenSet[str] = frozenset(
    {
        "HOH",
        "WAT",
        "TIP3",
        "TIP4",
        "TIP5",
        "SPC",
        "SPCE",
        "OPC",
    }
)

# Common ion residue names
ION_RESIDUE_NAMES: FrozenSet[str] = frozenset(
    {
        "NA",
        "CL",
        "K",
        "MG",
        "CA",
        "ZN",
        "FE",
        "CU",
        "NA+",
        "CL-",
        "K+",
        "MG2+",
        "CA2+",
        "ZN2+",
        "FE2+",
        "FE3+",
        "CU2+",
        # OpenMM/AMBER naming
        "SOD",
        "CLA",
        "POT",
        "MG2",
        "CAL",
    }
)

# Hydrogen element indicators
HYDROGEN_ELEMENTS: FrozenSet[str] = frozenset({"H"})


@dataclass
class SystemComponentInfo:
    """Information about system composition for atom group resolution.

    This dataclass stores the atom counts and chain assignments for each
    system component, enabling efficient resolution of predefined atom groups.

    Attributes:
        n_protein_atoms: Total atoms in protein/enzyme
        n_substrate_atoms: Total atoms in substrate/ligand
        n_polymer_atoms: Total atoms in all polymer chains
        protein_chain_id: Chain ID for protein (default "A")
        substrate_chain_id: Chain ID for substrate (default "B")
        polymer_chain_id: Chain ID for polymers (default "C")
        solvent_start_chain_id: First chain ID for solvent (default "D")
    """

    n_protein_atoms: int = 0
    n_substrate_atoms: int = 0
    n_polymer_atoms: int = 0
    protein_chain_id: str = "A"
    substrate_chain_id: str = "B"
    polymer_chain_id: str = "C"
    solvent_start_chain_id: str = "D"

    @property
    def has_protein(self) -> bool:
        """Check if system has protein atoms."""
        return self.n_protein_atoms > 0

    @property
    def has_substrate(self) -> bool:
        """Check if system has substrate atoms."""
        return self.n_substrate_atoms > 0

    @property
    def has_polymers(self) -> bool:
        """Check if system has polymer atoms."""
        return self.n_polymer_atoms > 0

    @classmethod
    def from_topology(cls, topology: OpenMMTopology) -> "SystemComponentInfo":
        """Reconstruct SystemComponentInfo from PDB topology chain assignments.

        Uses the PolyzyMD chain ID convention:
        - Chain A: Protein/Enzyme
        - Chain B: Substrate/Ligand
        - Chain C: Polymers
        - Chain D+: Solvent (water, ions, co-solvents)

        This method enables position restraints to work with pre-built systems
        loaded from PDB files (--skip-build mode).

        Args:
            topology: OpenMM Topology object (e.g., from PDBFile)

        Returns:
            SystemComponentInfo with atom counts derived from chain IDs

        Raises:
            ValueError: If topology has no chain information or is missing
                chain A (protein), indicating it was not built by PolyzyMD
        """
        n_protein = 0
        n_substrate = 0
        n_polymer = 0
        found_chains: set = set()

        for atom in topology.atoms():
            chain_id = atom.residue.chain.id if atom.residue.chain else ""
            found_chains.add(chain_id)

            if chain_id == "A":
                n_protein += 1
            elif chain_id == "B":
                n_substrate += 1
            elif chain_id == "C":
                n_polymer += 1
            # D+ are solvent - we don't need to track counts

        # Validate that we found expected chain structure
        if not found_chains:
            raise ValueError(
                "Topology has no chain information. Ensure the PDB was built by PolyzyMD "
                "with proper chain assignments (A=protein, B=substrate, C=polymers, D+=solvent)."
            )

        if "A" not in found_chains:
            raise ValueError(
                f"Topology missing chain A (protein). Found chains: {sorted(found_chains)}. "
                "Ensure the PDB was built by PolyzyMD with proper chain assignments."
            )

        logger.info(
            f"Parsed topology: {n_protein} protein atoms (chain A), "
            f"{n_substrate} substrate atoms (chain B), "
            f"{n_polymer} polymer atoms (chain C)"
        )

        return cls(
            n_protein_atoms=n_protein,
            n_substrate_atoms=n_substrate,
            n_polymer_atoms=n_polymer,
        )


class AtomGroupResolver:
    """Resolves predefined atom group names to atom indices.

    This class takes an OpenMM topology and system component information
    to resolve predefined group names (like 'protein_heavy') to lists of
    atom indices that can be used with positional restraints.

    The resolver uses chain IDs to identify system components:
    - Protein: Chain A (configurable)
    - Substrate: Chain B (configurable)
    - Polymers: Chain C (configurable)
    - Solvent: Chain D and beyond

    Example:
        >>> resolver = AtomGroupResolver(topology, component_info)
        >>> protein_heavy_atoms = resolver.resolve("protein_heavy")
        >>> backbone_atoms = resolver.resolve("protein_backbone")
    """

    def __init__(
        self,
        topology: OpenMMTopology,
        component_info: SystemComponentInfo,
    ):
        """Initialize the atom group resolver.

        Args:
            topology: OpenMM Topology object
            component_info: Information about system composition
        """
        self._topology = topology
        self._info = component_info

        # Build atom index lookup tables (lazy-initialized)
        self._atoms_by_chain: Optional[Dict[str, List[int]]] = None
        self._atom_data: Optional[List[Dict]] = None

    def _build_lookup_tables(self) -> None:
        """Build internal lookup tables for atom properties."""
        if self._atom_data is not None:
            return  # Already built

        self._atoms_by_chain = {}
        self._atom_data = []

        for atom in self._topology.atoms():
            chain_id = atom.residue.chain.id if atom.residue.chain else ""
            residue_name = atom.residue.name if atom.residue else ""
            atom_name = atom.name if atom.name else ""

            # Get element symbol
            element_symbol = ""
            if atom.element is not None:
                element_symbol = atom.element.symbol

            atom_info = {
                "index": atom.index,
                "name": atom_name.upper(),
                "element": element_symbol.upper(),
                "residue_name": residue_name.upper(),
                "chain_id": chain_id,
            }

            self._atom_data.append(atom_info)

            # Group by chain
            if chain_id not in self._atoms_by_chain:
                self._atoms_by_chain[chain_id] = []
            self._atoms_by_chain[chain_id].append(atom.index)

    def resolve(self, group_name: str) -> List[int]:
        """Resolve a predefined group name to atom indices.

        Args:
            group_name: One of the predefined group names

        Returns:
            List of atom indices (0-indexed). Empty list if group has no atoms.

        Raises:
            ValueError: If group_name is not a recognized predefined group
        """
        if group_name not in PREDEFINED_GROUPS:
            raise ValueError(
                f"Unknown atom group: '{group_name}'. Valid groups: {sorted(PREDEFINED_GROUPS)}"
            )

        # Build lookup tables on first use
        self._build_lookup_tables()

        # Dispatch to specific resolver method
        if group_name == "protein_heavy":
            return self._get_protein_heavy()
        elif group_name == "protein_backbone":
            return self._get_protein_backbone()
        elif group_name == "protein_calpha":
            return self._get_protein_calpha()
        elif group_name == "ligand_heavy":
            return self._get_ligand_heavy()
        elif group_name == "polymer_heavy":
            return self._get_polymer_heavy()
        elif group_name == "solvent":
            return self._get_solvent()
        elif group_name == "water_only":
            return self._get_water_only()
        elif group_name == "ions_only":
            return self._get_ions_only()
        elif group_name == "cosolvents_only":
            return self._get_cosolvents_only()
        else:
            # Should not reach here due to validation above
            raise ValueError(f"Unhandled group name: {group_name}")

    def _is_heavy_atom(self, atom_info: Dict) -> bool:
        """Check if an atom is a heavy atom (not hydrogen)."""
        return atom_info["element"] not in HYDROGEN_ELEMENTS

    def _get_chain_atoms(self, chain_id: str) -> List[int]:
        """Get all atom indices for a specific chain."""
        return self._atoms_by_chain.get(chain_id, [])

    def _get_protein_heavy(self) -> List[int]:
        """Get all non-hydrogen protein atoms."""
        if not self._info.has_protein:
            logger.warning("No protein atoms in system - 'protein_heavy' group is empty")
            return []

        chain_id = self._info.protein_chain_id
        heavy_atoms = []

        for idx in self._get_chain_atoms(chain_id):
            if self._is_heavy_atom(self._atom_data[idx]):
                heavy_atoms.append(idx)

        logger.debug(f"Resolved 'protein_heavy': {len(heavy_atoms)} atoms")
        return heavy_atoms

    def _get_protein_backbone(self) -> List[int]:
        """Get protein backbone atoms (N, CA, C, O)."""
        if not self._info.has_protein:
            logger.warning("No protein atoms in system - 'protein_backbone' group is empty")
            return []

        chain_id = self._info.protein_chain_id
        backbone_atoms = []

        for idx in self._get_chain_atoms(chain_id):
            atom_name = self._atom_data[idx]["name"]
            if atom_name in BACKBONE_ATOM_NAMES:
                backbone_atoms.append(idx)

        logger.debug(f"Resolved 'protein_backbone': {len(backbone_atoms)} atoms")
        return backbone_atoms

    def _get_protein_calpha(self) -> List[int]:
        """Get protein alpha carbon atoms only."""
        if not self._info.has_protein:
            logger.warning("No protein atoms in system - 'protein_calpha' group is empty")
            return []

        chain_id = self._info.protein_chain_id
        calpha_atoms = []

        for idx in self._get_chain_atoms(chain_id):
            if self._atom_data[idx]["name"] == "CA":
                calpha_atoms.append(idx)

        logger.debug(f"Resolved 'protein_calpha': {len(calpha_atoms)} atoms")
        return calpha_atoms

    def _get_ligand_heavy(self) -> List[int]:
        """Get all non-hydrogen ligand/substrate atoms."""
        if not self._info.has_substrate:
            logger.warning("No substrate atoms in system - 'ligand_heavy' group is empty")
            return []

        chain_id = self._info.substrate_chain_id
        heavy_atoms = []

        for idx in self._get_chain_atoms(chain_id):
            if self._is_heavy_atom(self._atom_data[idx]):
                heavy_atoms.append(idx)

        logger.debug(f"Resolved 'ligand_heavy': {len(heavy_atoms)} atoms")
        return heavy_atoms

    def _get_polymer_heavy(self) -> List[int]:
        """Get all non-hydrogen polymer atoms."""
        if not self._info.has_polymers:
            logger.warning("No polymer atoms in system - 'polymer_heavy' group is empty")
            return []

        chain_id = self._info.polymer_chain_id
        heavy_atoms = []

        for idx in self._get_chain_atoms(chain_id):
            if self._is_heavy_atom(self._atom_data[idx]):
                heavy_atoms.append(idx)

        logger.debug(f"Resolved 'polymer_heavy': {len(heavy_atoms)} atoms")
        return heavy_atoms

    def _get_solvent_chain_ids(self) -> Set[str]:
        """Get all chain IDs that belong to solvent."""
        # Solvent starts at the configured chain and goes through all remaining
        start_ord = ord(self._info.solvent_start_chain_id)
        solvent_chains = set()

        for chain_id in self._atoms_by_chain.keys():
            if len(chain_id) == 1 and ord(chain_id) >= start_ord:
                solvent_chains.add(chain_id)

        return solvent_chains

    def _get_solvent(self) -> List[int]:
        """Get all solvent atoms (water + ions + co-solvents)."""
        solvent_chains = self._get_solvent_chain_ids()
        solvent_atoms = []

        for chain_id in sorted(solvent_chains):
            solvent_atoms.extend(self._get_chain_atoms(chain_id))

        logger.debug(f"Resolved 'solvent': {len(solvent_atoms)} atoms")
        return solvent_atoms

    def _get_water_only(self) -> List[int]:
        """Get water molecule atoms only."""
        solvent_chains = self._get_solvent_chain_ids()
        water_atoms = []

        for chain_id in solvent_chains:
            for idx in self._get_chain_atoms(chain_id):
                residue_name = self._atom_data[idx]["residue_name"]
                if residue_name in WATER_RESIDUE_NAMES:
                    water_atoms.append(idx)

        logger.debug(f"Resolved 'water_only': {len(water_atoms)} atoms")
        return water_atoms

    def _get_ions_only(self) -> List[int]:
        """Get ion atoms only."""
        solvent_chains = self._get_solvent_chain_ids()
        ion_atoms = []

        for chain_id in solvent_chains:
            for idx in self._get_chain_atoms(chain_id):
                residue_name = self._atom_data[idx]["residue_name"]
                if residue_name in ION_RESIDUE_NAMES:
                    ion_atoms.append(idx)

        logger.debug(f"Resolved 'ions_only': {len(ion_atoms)} atoms")
        return ion_atoms

    def _get_cosolvents_only(self) -> List[int]:
        """Get co-solvent atoms only (not water, not ions)."""
        solvent_chains = self._get_solvent_chain_ids()
        cosolvent_atoms = []

        for chain_id in solvent_chains:
            for idx in self._get_chain_atoms(chain_id):
                residue_name = self._atom_data[idx]["residue_name"]
                if (
                    residue_name not in WATER_RESIDUE_NAMES
                    and residue_name not in ION_RESIDUE_NAMES
                ):
                    cosolvent_atoms.append(idx)

        logger.debug(f"Resolved 'cosolvents_only': {len(cosolvent_atoms)} atoms")
        return cosolvent_atoms

    def get_group_summary(self) -> Dict[str, int]:
        """Get a summary of atom counts for all predefined groups.

        Returns:
            Dictionary mapping group names to atom counts
        """
        self._build_lookup_tables()

        summary = {}
        for group_name in sorted(PREDEFINED_GROUPS):
            try:
                atoms = self.resolve(group_name)
                summary[group_name] = len(atoms)
            except Exception as e:
                logger.warning(f"Error resolving group '{group_name}': {e}")
                summary[group_name] = 0

        return summary
