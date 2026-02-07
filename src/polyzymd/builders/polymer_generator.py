"""
Polymer generator for dynamic polymer generation.

This module handles the building of complete polymer chains from MonomerGroup
fragments, including 3D structure generation, ring-piercing validation,
and partial charge assignment.

The workflow mirrors the notebook process:
1. Set terminal orientations based on sequence head/tail
2. Call build_linear_polymer() with middle sequence
3. Validate no ring-piercing (retry up to max_retries times)
4. Create OpenFF Topology and partition
5. Assign partial charges (NAGL/Espaloma/AM1BCC)
6. Cache as charged SDF

Made by PolyzyMD, by Joseph R. Laforet Jr.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, TYPE_CHECKING

from rdkit import Chem

from polymerist.polymers.monomers import MonomerGroup
from polymerist.polymers.building import build_linear_polymer, mbmol_to_openmm_pdb
from polymerist.polymers.building import mbmol_to_rdmol
from polymerist.rdutils.rdcoords.piercing import summarize_ring_piercing
from polymerist.mdtools.openfftools.topology import topology_from_sdf, topology_to_sdf
from polymerist.mdtools.openfftools.topology import get_largest_offmol
from polymerist.mdtools.openfftools.partition import partition
from polymerist.genutils.fileutils.pathutils import assemble_path

if TYPE_CHECKING:
    from openff.toolkit import Topology as OFFTopology
    from openff.toolkit import Molecule as OFFMolecule

logger = logging.getLogger(__name__)


class PolymerGenerationError(Exception):
    """Raised when polymer generation fails after all retries."""

    pass


class PolymerGenerator:
    """Generates complete polymer chains from MonomerGroup fragments.

    This class handles the full polymer building workflow:
    1. Building 3D structures using Polymerist's build_linear_polymer
    2. Validating no ring-piercing events occur
    3. Creating OpenFF topologies
    4. Assigning partial charges
    5. Caching results as SDF files

    Attributes:
        monomer_group: MonomerGroup containing polymer fragments
        cache_directory: Directory for caching generated polymers
        max_retries: Maximum attempts for building a valid polymer
        charger_type: Method for partial charge assignment
    """

    def __init__(
        self,
        monomer_group: MonomerGroup,
        cache_directory: Path,
        max_retries: int = 10,
        charger_type: str = "nagl",
    ):
        """Initialize the polymer generator.

        Args:
            monomer_group: MonomerGroup containing all named fragments
            cache_directory: Directory for caching polymer SDF files
            max_retries: Maximum attempts for building (ring-piercing failures)
            charger_type: Charge method ("nagl", "espaloma", "am1bcc")
        """
        self.monomer_group = monomer_group
        self.cache_directory = Path(cache_directory)
        self.max_retries = max_retries
        self.charger_type = charger_type.lower()

        # Ensure cache directory exists
        self.cache_directory.mkdir(parents=True, exist_ok=True)

        # Lazy-loaded charger
        self._charger = None

    @property
    def charger(self):
        """Get the molecular charger (lazy-loaded)."""
        if self._charger is None:
            self._charger = self._create_charger()
        return self._charger

    def _create_charger(self):
        """Create the appropriate charger based on charger_type."""
        if self.charger_type == "nagl":
            from polymerist.mdtools.openfftools.partialcharge.molchargers import (
                NAGLCharger,
            )

            return NAGLCharger()
        elif self.charger_type == "espaloma":
            from polymerist.mdtools.openfftools.partialcharge.molchargers import (
                EspalomaCharger,
            )

            return EspalomaCharger()
        elif self.charger_type == "am1bcc":
            from polymerist.mdtools.openfftools.partialcharge.molchargers import (
                AM1BCCCharger,
            )

            return AM1BCCCharger()
        else:
            raise ValueError(f"Unknown charger type: {self.charger_type}")

    def _get_terminal_fragment_name(self, monomer_name: str) -> str:
        """Get the 1-site fragment name for a monomer (terminal position).

        Args:
            monomer_name: Base monomer name (e.g., "SBMA")

        Returns:
            Fragment name for 1-site (e.g., "SBMA-1a")
        """
        # Look for matching 1-site fragment in monomer_group
        for frag_name in self.monomer_group.monomers.keys():
            if frag_name.startswith(monomer_name) and "-1" in frag_name:
                return frag_name

        raise ValueError(f"No 1-site terminal fragment found for monomer: {monomer_name}")

    def _get_middle_fragment_name(self, monomer_name: str) -> str:
        """Get the 2-site fragment name for a monomer (middle position).

        Args:
            monomer_name: Base monomer name (e.g., "SBMA")

        Returns:
            Fragment name for 2-site (e.g., "SBMA-2a")
        """
        # Look for matching 2-site fragment in monomer_group
        for frag_name in self.monomer_group.monomers.keys():
            if frag_name.startswith(monomer_name) and "-2" in frag_name:
                return frag_name

        raise ValueError(f"No 2-site middle fragment found for monomer: {monomer_name}")

    def _build_polymer_structure(
        self,
        sequence: str,
        monomer_names: Dict[str, str],
        residue_names: Optional[Dict[str, str]] = None,
    ) -> Tuple[any, Path]:
        """Build a polymer 3D structure from sequence.

        Args:
            sequence: Polymer sequence string (e.g., "ABCAB")
            monomer_names: Mapping of sequence labels to monomer names
            residue_names: Optional mapping of monomer names to 3-char residue names

        Returns:
            Tuple of (mbuild compound, path to PDB file)

        Raises:
            PolymerGenerationError: If building fails after max_retries
        """
        # Parse sequence
        head_label = sequence[0]
        tail_label = sequence[-1]
        middle_labels = sequence[1:-1] if len(sequence) > 2 else ""

        head_monomer = monomer_names[head_label]
        tail_monomer = monomer_names[tail_label]

        # Set terminal orientations
        monogrp_local = MonomerGroup(monomers=self.monomer_group.monomers)
        monogrp_local.term_orient = {
            "head": self._get_terminal_fragment_name(head_monomer),
            "tail": self._get_terminal_fragment_name(tail_monomer),
        }

        # Build middle sequence string
        middle_sequence = (
            "".join(self._get_middle_fragment_name(monomer_names[label]) for label in middle_labels)
            if middle_labels
            else ""
        )

        # Attempt building with retries for ring-piercing
        for attempt in range(self.max_retries):
            logger.debug(f"Building polymer attempt {attempt + 1}/{self.max_retries}")

            chain = build_linear_polymer(
                monomers=monogrp_local,
                n_monomers=len(sequence),
                sequence=middle_sequence,
                energy_minimize=True,
                allow_partial_sequences=True,
            )

            # Check for ring-piercing
            poly_mol = mbmol_to_rdmol(chain)
            piercing_summary = summarize_ring_piercing(poly_mol)

            if not piercing_summary:
                logger.info(f"Polymer built successfully on attempt {attempt + 1}")
                break
            else:
                logger.warning(
                    f"Ring-piercing detected on attempt {attempt + 1}: {piercing_summary}"
                )
        else:
            raise PolymerGenerationError(
                f"Failed to build polymer after {self.max_retries} attempts due to ring-piercing"
            )

        # Save PDB
        pdb_filename = self._make_polymer_filename(sequence, monomer_names, charged=False)
        pdb_path = self.cache_directory / f"{pdb_filename}.pdb"

        resname_map = self._build_resname_map(monomer_names, residue_names)
        mbmol_to_openmm_pdb(pdb_path, chain, resname_map=resname_map)
        logger.info(f"Saved polymer PDB: {pdb_path}")

        return chain, pdb_path

    def _build_resname_map(
        self,
        monomer_names: Dict[str, str],
        residue_names: Optional[Dict[str, str]] = None,
    ) -> Dict[str, str]:
        """Build residue name mapping for PDB output.

        Args:
            monomer_names: Mapping of labels to monomer names
            residue_names: Optional custom residue names

        Returns:
            Mapping of fragment names to 3-char residue names
        """
        resname_map = {}
        for label, monomer_name in monomer_names.items():
            # Get 3-char residue name
            if residue_names and monomer_name in residue_names:
                base_resname = residue_names[monomer_name]
            else:
                base_resname = monomer_name[:3].upper()

            # Map both 1-site and 2-site fragments
            for frag_name in self.monomer_group.monomers.keys():
                if frag_name.startswith(monomer_name):
                    if "-1" in frag_name:
                        resname_map[frag_name] = f"{base_resname[:2]}1"
                    elif "-2" in frag_name:
                        resname_map[frag_name] = f"{base_resname[:2]}2"

        return resname_map

    def _make_polymer_filename(
        self,
        sequence: str,
        monomer_names: Dict[str, str],
        charged: bool = True,
    ) -> str:
        """Create a descriptive filename for a polymer.

        Args:
            sequence: Polymer sequence string
            monomer_names: Mapping of labels to monomer names
            charged: Whether this is a charged structure

        Returns:
            Filename without extension
        """
        # Get unique monomers in sequence order
        unique_labels = []
        for label in sequence:
            if label not in unique_labels:
                unique_labels.append(label)

        # Build monomer prefix
        monomers_used = [monomer_names[label] for label in unique_labels]
        monomer_prefix = "-".join(monomers_used) if monomers_used else "NO_MONOMERS"

        length = len(sequence)
        filename = f"{monomer_prefix}_seq={sequence}_{length}-mer"
        if charged:
            filename += "_charged"

        return filename

    def generate_polymer(
        self,
        sequence: str,
        monomer_names: Dict[str, str],
        residue_names: Optional[Dict[str, str]] = None,
        force_regenerate: bool = False,
    ) -> "OFFMolecule":
        """Generate a complete, charged polymer molecule.

        This method handles the full workflow:
        1. Check cache for existing charged SDF
        2. Build 3D structure if needed
        3. Create OpenFF topology
        4. Assign partial charges
        5. Cache the result

        Args:
            sequence: Polymer sequence string (e.g., "ABCAB")
            monomer_names: Mapping of sequence labels to monomer names
            residue_names: Optional mapping of monomer names to 3-char residue names
            force_regenerate: If True, regenerate even if cached

        Returns:
            OpenFF Molecule with partial charges assigned

        Raises:
            PolymerGenerationError: If generation fails
        """
        from openff.toolkit import Topology as OFFTopology

        # Check cache
        charged_filename = self._make_polymer_filename(sequence, monomer_names, charged=True)
        charged_sdf_path = self.cache_directory / f"{charged_filename}.sdf"

        if charged_sdf_path.exists() and not force_regenerate:
            logger.info(f"Loading cached polymer from {charged_sdf_path}")
            off_top = topology_from_sdf(charged_sdf_path)
            return get_largest_offmol(off_top)

        # Build structure
        chain, pdb_path = self._build_polymer_structure(sequence, monomer_names, residue_names)

        # Create topology
        uncharged_filename = self._make_polymer_filename(sequence, monomer_names, charged=False)
        uncharged_sdf_path = self.cache_directory / f"{uncharged_filename}.sdf"

        logger.info("Creating OpenFF topology...")
        off_top = OFFTopology.from_pdb(
            str(pdb_path), _custom_substructures=self.monomer_group.monomers
        )
        was_partitioned = partition(off_top)
        if not was_partitioned:
            raise PolymerGenerationError("Failed to partition polymer topology")

        # Fix residue names (truncate to 3 chars)
        for mol in off_top.molecules:
            for atom in mol.atoms:
                if "residue_name" in atom.metadata:
                    atom.metadata["extended_name"] = atom.metadata["residue_name"]
                    atom.metadata["residue_name"] = atom.metadata["residue_name"][:3]

        # Save uncharged SDF
        topology_to_sdf(uncharged_sdf_path, off_top)
        logger.info(f"Saved uncharged SDF: {uncharged_sdf_path}")

        # Assign partial charges
        logger.info(f"Assigning partial charges using {self.charger_type}...")
        off_mol = get_largest_offmol(off_top)
        charged_mol = self.charger.charge_molecule(off_mol)

        # Save charged SDF
        charged_top = charged_mol.to_topology()
        topology_to_sdf(charged_sdf_path, charged_top)
        logger.info(f"Saved charged SDF: {charged_sdf_path}")

        return charged_mol

    def generate_polymers_batch(
        self,
        sequences: List[str],
        monomer_names: Dict[str, str],
        residue_names: Optional[Dict[str, str]] = None,
    ) -> Dict[str, "OFFMolecule"]:
        """Generate multiple polymers (sequential processing).

        Args:
            sequences: List of polymer sequences to generate
            monomer_names: Mapping of sequence labels to monomer names
            residue_names: Optional mapping of monomer names to 3-char residue names

        Returns:
            Dictionary mapping sequence -> OpenFF Molecule
        """
        results = {}
        for i, sequence in enumerate(sequences):
            logger.info(f"Generating polymer {i + 1}/{len(sequences)}: {sequence}")
            try:
                mol = self.generate_polymer(sequence, monomer_names, residue_names)
                results[sequence] = mol
            except Exception as e:
                logger.error(f"Failed to generate polymer {sequence}: {e}")
                raise

        return results

    def get_cached_polymer(
        self,
        sequence: str,
        monomer_names: Dict[str, str],
    ) -> Optional["OFFMolecule"]:
        """Load a polymer from cache if it exists.

        Args:
            sequence: Polymer sequence string
            monomer_names: Mapping of labels to monomer names

        Returns:
            OpenFF Molecule if cached, None otherwise
        """
        charged_filename = self._make_polymer_filename(sequence, monomer_names, charged=True)
        charged_sdf_path = self.cache_directory / f"{charged_filename}.sdf"

        if charged_sdf_path.exists():
            off_top = topology_from_sdf(charged_sdf_path)
            return get_largest_offmol(off_top)

        return None
