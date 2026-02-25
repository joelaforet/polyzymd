"""
System builder for orchestrating the complete simulation system construction.

This module coordinates all builders to create a complete, solvated
molecular system and generate the OpenFF Interchange for simulation.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, Topology

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig

from polyzymd.builders.enzyme import EnzymeBuilder
from polyzymd.builders.polymer import PolymerBuilder
from polyzymd.builders.solvent import SolventBuilder, SolventComposition
from polyzymd.builders.substrate import SubstrateBuilder

LOGGER = logging.getLogger(__name__)


class SystemBuilder:
    """Orchestrator for building complete simulation systems.

    This class coordinates the 5-stage build pipeline:
    1. Enzyme loading and partitioning
    2. Substrate loading and charging
    3. Polymer generation and packing
    4. Solvation with water, ions, and co-solvents
    5. Interchange creation with optimized combining

    The Interchange.combine() optimization is used when polymers are present,
    which significantly speeds up parameterization for systems with many
    unique molecules.

    Example:
        >>> builder = SystemBuilder.from_config(config)
        >>> interchange = builder.build()
        >>> # Or step by step:
        >>> builder = SystemBuilder()
        >>> builder.build_enzyme("enzyme.pdb")
        >>> builder.build_substrate("docked.sdf")
        >>> builder.build_polymers(...)
        >>> builder.solvate()
        >>> interchange = builder.create_interchange()
    """

    def __init__(
        self,
        protein_forcefield: str = "ff14sb_off_impropers_0.0.4.offxml",
        small_molecule_forcefield: str = "openff-2.0.0.offxml",
    ) -> None:
        """Initialize the SystemBuilder.

        Args:
            protein_forcefield: Force field for proteins.
            small_molecule_forcefield: Force field for small molecules.
        """
        self._protein_ff = protein_forcefield
        self._sm_ff = small_molecule_forcefield

        # Sub-builders
        self._enzyme_builder = EnzymeBuilder()
        self._substrate_builder = SubstrateBuilder()
        self._polymer_builder: Optional[PolymerBuilder] = None
        self._solvent_builder = SolventBuilder()

        # State
        self._enzyme_topology: Optional[Topology] = None
        self._substrate_molecule: Optional[Molecule] = None
        self._polymer_molecules: List[Molecule] = []
        self._polymer_counts: List[int] = []
        self._combined_topology: Optional[Topology] = None
        self._solvated_topology: Optional[Topology] = None
        self._interchange: Optional[Interchange] = None
        self._working_dir: Optional[Path] = None

        # Molecule count tracking for PDB chain/residue assignment
        # These are set during build and used by _assign_pdb_identifiers()
        self._n_enzyme_molecules: int = 0
        self._n_substrate_molecules: int = 0
        self._n_polymer_chains: int = 0

    @property
    def interchange(self) -> Optional[Interchange]:
        """Get the created Interchange object."""
        return self._interchange

    @property
    def solvated_topology(self) -> Optional[Topology]:
        """Get the solvated topology."""
        return self._solvated_topology

    @classmethod
    def from_config(cls, config: "SimulationConfig") -> "SystemBuilder":
        """Create a SystemBuilder from a configuration object.

        Args:
            config: SimulationConfig with all system settings.

        Returns:
            Configured SystemBuilder instance.
        """
        builder = cls(
            protein_forcefield=config.force_field.protein,
            small_molecule_forcefield=config.force_field.small_molecule,
        )
        builder._config = config
        return builder

    def build_enzyme(self, pdb_path: Union[str, Path]) -> Topology:
        """Build the enzyme component.

        Args:
            pdb_path: Path to enzyme PDB file.

        Returns:
            Enzyme topology.
        """
        self._enzyme_topology = self._enzyme_builder.build(pdb_path)
        self._n_enzyme_molecules = 1
        return self._enzyme_topology

    def build_substrate(
        self,
        sdf_path: Union[str, Path],
        conformer_index: int = 0,
        charge_method: str = "nagl",
        residue_name: str = "LIG",
    ) -> Molecule:
        """Build the substrate component.

        Args:
            sdf_path: Path to docked conformers SDF.
            conformer_index: Index of conformer to use.
            charge_method: Charge assignment method.
            residue_name: 3-letter residue name.

        Returns:
            Charged substrate molecule.
        """
        self._substrate_molecule = self._substrate_builder.build(
            sdf_path=sdf_path,
            conformer_index=conformer_index,
            charge_method=charge_method,
            residue_name=residue_name,
        )
        self._n_substrate_molecules = 1
        return self._substrate_molecule

    def build_polymers(
        self,
        characters: List[str],
        probabilities: List[float],
        length: int,
        count: int,
        type_prefix: str,
        sdf_directory: Optional[Union[str, Path]] = None,
        seed: Optional[int] = None,
        # New parameters for dynamic generation
        generation_mode: str = "cached",
        monomer_smiles: Optional[Dict[str, str]] = None,
        monomer_names: Optional[Dict[str, str]] = None,
        residue_names: Optional[Dict[str, str]] = None,
        reactions: Optional[Any] = None,
        charger_type: str = "nagl",
        max_retries: int = 10,
        cache_directory: Optional[Union[str, Path]] = None,
    ) -> Tuple[List[Molecule], List[int]]:
        """Build polymer components.

        Supports two generation modes:
        - "cached": Load pre-built polymer SDFs from sdf_directory
        - "dynamic": Generate polymers on-the-fly from monomer SMILES

        Args:
            characters: Monomer labels.
            probabilities: Monomer selection probabilities.
            length: Monomers per chain.
            count: Number of polymer chains.
            type_prefix: Filename prefix.
            sdf_directory: Directory with pre-built polymer SDFs (cached mode).
            seed: Random seed for reproducibility.
            generation_mode: "cached" or "dynamic".
            monomer_smiles: Dict of monomer name -> SMILES (dynamic mode).
            monomer_names: Dict of label -> monomer name (dynamic mode).
            residue_names: Dict of monomer name -> 3-char residue name.
            reactions: ReactionConfig with ATRP reaction paths (dynamic mode).
            charger_type: Charge method ("nagl", "espaloma", "am1bcc").
            max_retries: Max retries for polymer generation.
            cache_directory: Directory for caching generated polymers.

        Returns:
            Tuple of (unique polymer molecules, counts).
        """
        self._polymer_builder = PolymerBuilder(
            characters=characters,
            probabilities=probabilities,
            length=length,
            type_prefix=type_prefix,
            sdf_directory=sdf_directory,
            cache_directory=cache_directory,
            generation_mode=generation_mode,
            monomer_smiles=monomer_smiles,
            monomer_names=monomer_names,
            residue_names=residue_names,
            reactions=reactions,
            charger_type=charger_type,
            max_retries=max_retries,
        )

        molecules, counts = self._polymer_builder.build(count=count, seed=seed)
        self._polymer_molecules = molecules
        self._polymer_counts = counts
        self._n_polymer_chains = sum(counts)

        return molecules, counts

    def combine_solutes(self) -> Topology:
        """Combine enzyme, substrate, and polymers into a single topology.

        Returns:
            Combined topology ready for solvation.

        Raises:
            RuntimeError: If enzyme has not been built.
        """
        if self._enzyme_topology is None:
            raise RuntimeError("Enzyme must be built before combining solutes")

        LOGGER.info("Combining solute components")

        # Start with enzyme
        molecules = [self._enzyme_topology.molecule(0)]

        # Add substrate if present
        if self._substrate_molecule is not None:
            molecules.append(self._substrate_molecule)

        # Create combined topology
        self._combined_topology = Topology.from_molecules(molecules)

        # Re-number chains
        self._renumber_chains(self._combined_topology)

        LOGGER.info(
            f"Combined topology: {self._combined_topology.n_molecules} molecules, "
            f"{self._combined_topology.n_atoms} atoms"
        )

        return self._combined_topology

    def pack_polymers(
        self,
        padding: float = 2.0,
        tolerance: float = 2.0,
        movebadrandom: bool = False,
        working_directory: Optional[Union[str, Path]] = None,
    ) -> Topology:
        """Pack polymers around the combined solute topology.

        Args:
            padding: Box padding in nm. Larger values give polymers more room
                and can significantly speed up PACKMOL convergence.
            tolerance: PACKMOL tolerance in Angstrom.
            movebadrandom: When True, pass the ``movebadrandom`` keyword to
                PACKMOL. Improves convergence for dense or heterogeneous
                polymer systems (many unique chain types) by placing
                badly-packed molecules at random positions in the box.
            working_directory: Directory for PACKMOL files.

        Returns:
            Topology with polymers packed.

        Raises:
            RuntimeError: If solutes not combined or no polymers built.
        """
        if self._combined_topology is None:
            self.combine_solutes()

        if not self._polymer_molecules:
            LOGGER.info("No polymers to pack, returning combined topology")
            return self._combined_topology

        from openff.units import Quantity

        from polyzymd.utils import boxvectors
        from polyzymd.utils.packmol import pack_polymers as _pack_polymers

        LOGGER.info(
            f"Packing {sum(self._polymer_counts)} polymer chains "
            f"({len(self._polymer_molecules)} unique types), "
            f"padding={padding} nm, tolerance={tolerance} A"
            + (" [movebadrandom]" if movebadrandom else "")
        )

        # Calculate box vectors
        padding_qty = Quantity(padding, "nanometer")
        min_box_vecs = boxvectors.get_topology_bbox(self._combined_topology)
        box_vecs = boxvectors.pad_box_vectors_uniform(min_box_vecs, padding_qty)

        # Pack polymers using our custom PACKMOL runner (supports movebadrandom)
        packed_top = _pack_polymers(
            molecules=self._polymer_molecules,
            number_of_copies=self._polymer_counts,
            solute=self._combined_topology,
            box_vectors=box_vecs,
            tolerance_angstrom=tolerance,
            movebadrandom=movebadrandom,
            working_directory=str(working_directory) if working_directory else None,
            retain_working_files=True,
        )

        # Re-number chains
        self._renumber_chains(packed_top)

        self._combined_topology = packed_top

        LOGGER.info(f"Packed topology: {packed_top.n_atoms} atoms")

        return packed_top

    def solvate(
        self,
        composition: Optional[SolventComposition] = None,
        padding: float = 0.9,
        box_shape: str = "rhombic_dodecahedron",
    ) -> Topology:
        """Solvate the system with water and ions.

        Args:
            composition: Solvent composition specification.
            padding: Box padding in nm.
            box_shape: Box geometry.

        Returns:
            Solvated topology.

        Raises:
            RuntimeError: If solutes not combined.
        """
        if self._combined_topology is None:
            self.combine_solutes()

        LOGGER.info("Solvating system")

        self._solvated_topology = self._solvent_builder.solvate(
            topology=self._combined_topology,
            composition=composition,
            padding=padding,
            box_shape=box_shape,
        )

        return self._solvated_topology

    def create_interchange(
        self,
        use_optimized_combining: bool = True,
    ) -> Interchange:
        """Create the OpenFF Interchange for simulation.

        When polymers are present and use_optimized_combining is True,
        uses the Interchange.combine() optimization which is significantly
        faster for systems with many unique molecules.

        This method uses pre-computed charges from PolyzyMD's solvent cache
        to ensure all copies of water and co-solvent molecules have identical
        parameters, avoiding per-molecule AM1BCC calculations.

        Args:
            use_optimized_combining: Use optimized combining for polymers.

        Returns:
            OpenFF Interchange object.

        Raises:
            RuntimeError: If system not solvated.
        """
        if self._solvated_topology is None:
            raise RuntimeError("System must be solvated before creating Interchange")

        from polyzymd.data.solvent_molecules import get_solvent_molecule

        LOGGER.info("Creating Interchange")

        ff = ForceField(self._protein_ff, self._sm_ff)

        # Get water molecule with pre-computed charges for charge_from_molecules
        # Use PolyzyMD's cached water instead of Polymerist's for consistency
        water_model = "tip3p"
        if self._solvent_builder._composition:
            water_model = self._solvent_builder._composition.water_model
        water_mol = get_solvent_molecule(water_model)

        # Decide whether to use optimized combining
        use_combining = (
            use_optimized_combining and self._polymer_molecules and len(self._polymer_molecules) > 0
        )

        if use_combining:
            self._interchange = self._create_interchange_optimized(ff, water_mol)
        else:
            self._interchange = self._create_interchange_simple(ff, water_mol)

        LOGGER.info("Interchange created successfully")

        return self._interchange

    def _build_molecule_name_mapping(self) -> Dict[str, str]:
        """Build SMILES -> friendly name mapping for logging.

        Creates a mapping from molecule SMILES to human-readable names
        for use in log messages during Interchange creation.

        Returns:
            Dictionary mapping SMILES strings to display names.

        Example output:
            {
                "[large-enzyme-smiles]": "LipA",
                "[substrate-smiles]": "ResorufinButyrate",
                "[polymer-smiles]": "EGPMA-SBMA_AAABA",
                "CS(=O)C": "dmso",
            }
        """
        name_map: Dict[str, str] = {}

        # Add enzyme name from config or default
        if self._enzyme_topology:
            enzyme_smiles = self._enzyme_topology.molecule(0).to_smiles()
            config = getattr(self, "_config", None)
            if config and hasattr(config, "enzyme") and config.enzyme:
                name_map[enzyme_smiles] = config.enzyme.name
            else:
                name_map[enzyme_smiles] = "Enzyme"

        # Add substrate name from config or default
        if self._substrate_molecule:
            substrate_smiles = self._substrate_molecule.to_smiles()
            config = getattr(self, "_config", None)
            if config and hasattr(config, "substrate") and config.substrate:
                name_map[substrate_smiles] = config.substrate.name
            else:
                name_map[substrate_smiles] = "Substrate"

        # Add polymer names: type_prefix + sequence (e.g., "EGPMA-SBMA_AAABA")
        if self._polymer_builder:
            for sequence, mol in self._polymer_builder.loaded_molecules.items():
                polymer_smiles = mol.to_smiles()
                name_map[polymer_smiles] = f"{self._polymer_builder._type_prefix}_{sequence}"

        # Add co-solvent names
        if self._solvent_builder._composition:
            for cosolvent in self._solvent_builder._composition.co_solvents:
                if cosolvent.molecule:
                    cosolvent_smiles = cosolvent.molecule.to_smiles()
                    # Use name if available, otherwise use SMILES
                    name_map[cosolvent_smiles] = cosolvent.name

        return name_map

    def _build_charge_template_mapping(self, water_mol: Molecule) -> Dict[str, Molecule]:
        """Build SMILES -> charge template molecule mapping.

        Creates a mapping from molecule SMILES to template molecules with
        pre-computed charges for use with charge_from_molecules parameter.

        Args:
            water_mol: Water molecule with pre-computed charges.

        Returns:
            Dictionary mapping SMILES strings to charged template molecules.
        """
        smiles_to_template: Dict[str, Molecule] = {}

        # Add water
        smiles_to_template[water_mol.to_smiles()] = water_mol

        # Add polymers (already have charges from SDF)
        for mol in self._polymer_molecules:
            smiles_to_template[mol.to_smiles()] = mol

        # Add substrate (already has charges from NAGL/AM1BCC)
        if self._substrate_molecule:
            smiles_to_template[self._substrate_molecule.to_smiles()] = self._substrate_molecule

        # Add co-solvents with pre-computed charges
        if self._solvent_builder._composition:
            for cosolvent in self._solvent_builder._composition.co_solvents:
                if cosolvent.molecule is not None:
                    smiles_to_template[cosolvent.molecule.to_smiles()] = cosolvent.molecule

        return smiles_to_template

    def _create_interchange_simple(self, ff: ForceField, water_mol: Molecule) -> Interchange:
        """Create Interchange using batched molecule processing.

        Groups molecules by type (SMILES) and creates one Interchange per
        unique molecule type, then combines them. This is much more efficient
        than creating one Interchange per molecule instance.

        Args:
            ff: OpenFF ForceField.
            water_mol: Water molecule with pre-computed charges.

        Returns:
            Interchange object.
        """
        # Use the same optimized batching logic
        return self._create_interchange_batched(ff, water_mol)

    def _create_interchange_optimized(self, ff: ForceField, water_mol: Molecule) -> Interchange:
        """Create Interchange using optimized batched combining.

        Groups molecules by type (SMILES) and creates one Interchange per
        unique molecule type, then combines them. This dramatically reduces
        the number of Interchange.combine() calls needed.

        For a system with 1133 DMSO molecules, this creates ~7 Interchanges
        instead of ~1140, providing ~160x fewer combine operations.

        Uses pre-computed charges for water and co-solvents to avoid running
        AM1BCC charge calculation for every solvent molecule.

        Args:
            ff: OpenFF ForceField.
            water_mol: Water molecule with pre-computed charges.

        Returns:
            Combined Interchange object.
        """
        return self._create_interchange_batched(ff, water_mol)

    def _create_interchange_batched(self, ff: ForceField, water_mol: Molecule) -> Interchange:
        """Create Interchange using batched molecule processing with preserved order.

        This implementation batches molecules by type for efficiency while preserving
        the exact molecule order from self._solvated_topology. This is critical for
        DCD trajectory compatibility - the atom order in the Interchange must match
        the atom order in solvated_system.pdb.

        Batching strategy:
        1. Build an ``id(molecule) -> SMILES`` cache in ONE pass (avoids calling
           ``to_smiles()`` three separate times per molecule).
        2. Parameterize each unique molecule type ONCE (cache parameters).
        3. Create individual Interchanges in original topology order, grouping
           consecutive same-type molecules into single batches.
        4. Combine Interchanges using an **adjacent-pair tree reduction** instead
           of a serial left-fold.  This keeps the atom order identical to the
           original topology while reducing the asymptotic cost of
           ``Interchange.combine()`` from O(N²) to O(N log N).

        Note:
            When combining Interchanges from molecules parameterized by different
            force fields (e.g., ff14SB for proteins + OpenFF 2.0 for small molecules),
            you may see "Key collision with different parameters" warnings. This is
            expected behavior - the same SMIRKS pattern (e.g., [#6X4:1]-[#1:2] for
            C-H bonds) can have different parameters in different force fields.
            OpenFF handles this by appending "_DUPLICATE" to the key, allowing both
            parameter sets to coexist. See OpenFF's "Sharp Edges" documentation:
            https://docs.openforcefield.org/projects/interchange/en/stable/using/edges.html

        Args:
            ff: OpenFF ForceField.
            water_mol: Water molecule with pre-computed charges.

        Returns:
            Combined Interchange object.
        """
        LOGGER.info("Using batched Interchange creation (preserving molecule order)")

        # Build helper mappings
        smiles_to_name = self._build_molecule_name_mapping()
        smiles_to_template = self._build_charge_template_mapping(water_mol)

        # ── Step 0: Cache to_smiles() for every molecule in ONE pass ──
        # molecule.to_smiles() is expensive; previous code called it 3×.
        t0 = time.perf_counter()
        mol_id_to_smiles: Dict[int, str] = {}
        smiles_counts: Dict[str, int] = {}
        first_mol_for_smiles: Dict[str, Molecule] = {}

        for molecule in self._solvated_topology.molecules:
            mid = id(molecule)
            mol_smiles = molecule.to_smiles()
            mol_id_to_smiles[mid] = mol_smiles
            smiles_counts[mol_smiles] = smiles_counts.get(mol_smiles, 0) + 1
            if mol_smiles not in first_mol_for_smiles:
                first_mol_for_smiles[mol_smiles] = molecule

        LOGGER.info(
            f"SMILES cache built in {time.perf_counter() - t0:.1f}s "
            f"({len(mol_id_to_smiles)} molecules, {len(smiles_counts)} unique types)"
        )

        # ── Step 1: Prepare per-type metadata (no parameterisation yet) ──
        smiles_info: Dict[str, dict] = {}
        for mol_smiles, count in smiles_counts.items():
            display_name = smiles_to_name.get(mol_smiles, mol_smiles[:50] + "...")
            LOGGER.info(f"  {display_name}: {count} instance(s)")

            charge_from = []
            if mol_smiles in smiles_to_template:
                charge_from = [smiles_to_template[mol_smiles]]

            smiles_info[mol_smiles] = {
                "template": first_mol_for_smiles[mol_smiles],
                "charge_from": charge_from,
                "display_name": display_name,
            }

        # ── Step 2: Create Interchanges in ORIGINAL TOPOLOGY ORDER ──
        # Group consecutive molecules of the same type into single batches.
        t1 = time.perf_counter()
        all_interchanges: List[Interchange] = []
        interchange_names: List[str] = []

        current_smiles: str | None = None
        current_batch: List[Molecule] = []

        def flush_batch():
            """Create interchange for current batch and add to list."""
            nonlocal current_batch, current_smiles
            if not current_batch:
                return

            info = smiles_info[current_smiles]
            display_name = info["display_name"]
            charge_from = info["charge_from"]

            LOGGER.debug(f"Creating Interchange for {len(current_batch)} {display_name}")

            inc = ff.create_interchange(
                topology=current_batch,
                charge_from_molecules=charge_from,
            )
            all_interchanges.append(inc)
            interchange_names.append(f"{len(current_batch)}x {display_name}")

            current_batch = []

        for molecule in self._solvated_topology.molecules:
            mol_smiles = mol_id_to_smiles[id(molecule)]

            if mol_smiles != current_smiles:
                flush_batch()
                current_smiles = mol_smiles

            current_batch.append(molecule)

        flush_batch()

        t2 = time.perf_counter()
        LOGGER.info(f"Created {len(all_interchanges)} Interchange batch(es) in {t2 - t1:.1f}s")
        if len(interchange_names) <= 10:
            LOGGER.info(f"  Batches: {', '.join(interchange_names)}")
        else:
            LOGGER.info(f"  First 5: {', '.join(interchange_names[:5])}")
            LOGGER.info(f"  Last 5: {', '.join(interchange_names[-5:])}")

        if not all_interchanges:
            raise RuntimeError("No molecules found in solvated topology")

        # ── Step 3: Serial combine (preserving order) ──
        # OpenFF Interchange.combine() creates _DUPLICATE parameter keys when
        # SMIRKS patterns resolve to different force constants in different
        # molecular contexts.  Combining two objects that *both* contain
        # _DUPLICATE keys raises UnsupportedCombinationError, which rules out
        # tree-based (pairwise) reduction.  We therefore keep the serial
        # left-fold: combined = ((A + B) + C) + D.
        #
        # The O(N²) cost is mitigated by having very few batches (one per
        # contiguous same-SMILES run), typically 5-10 for a solvated system.
        combined = all_interchanges[0]
        for idx, inc in enumerate(all_interchanges[1:], 1):
            t_step = time.perf_counter()
            combined = combined.combine(inc)
            elapsed = time.perf_counter() - t_step
            LOGGER.info(
                f"  Combined batch {idx}/{len(all_interchanges) - 1} "
                f"({interchange_names[idx]}) in {elapsed:.1f}s"
            )

        t3 = time.perf_counter()
        LOGGER.info(f"Interchange combination complete ({t3 - t2:.1f}s total)")

        # Reset box vectors
        combined.box = self._solvated_topology.box_vectors

        return combined

    def _renumber_chains(self, topology: Topology) -> None:
        """Re-number chain IDs in the topology.

        Note: This is the legacy method. For PDB-compliant output, use
        _assign_pdb_identifiers() which is called by save_topology().

        Args:
            topology: Topology to modify.
        """
        for i, mol in enumerate(topology.molecules):
            for atom in mol.atoms:
                atom.metadata["chain_id"] = str(i + 1)

    def _assign_pdb_identifiers(self) -> None:
        """Assign PDB-compliant chain IDs and residue numbers to the solvated topology.

        This method assigns chain IDs and residue numbers based on molecule types
        tracked during the build process. The config YAML serves as the single
        source of truth for what each molecule represents.

        Chain assignment:
        - Protein: First chain (A), preserves original residue numbers from input PDB
        - Substrate: Next chain (B if present), residue 1
        - Polymers: Next chain (C if present), preserves per-monomer residue numbers
        - Solvent: Remaining chains with overflow at 9999 residues per chain

        This ensures every atom can be uniquely identified by the tuple
        (chain_id, residue_number, residue_name, atom_name) for downstream
        analysis tools like MDAnalysis and PyMOL.

        Raises:
            RuntimeError: If no solvated topology exists.
        """
        if self._solvated_topology is None:
            raise RuntimeError("No solvated topology. Call solvate() first.")

        CHAIN_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        chain_idx = 0
        mol_idx = 0

        # 1. Protein: Assign chain letter, preserve original residue numbers
        if self._n_enzyme_molecules > 0:
            chain_id = CHAIN_LETTERS[chain_idx]
            LOGGER.debug(f"Assigning chain {chain_id} to protein")

            for _ in range(self._n_enzyme_molecules):
                mol = self._solvated_topology.molecule(mol_idx)
                for atom in mol.atoms:
                    atom.metadata["chain_id"] = chain_id
                    # Ensure residue_number is a string (PDB loader may store as int)
                    # OpenMM's addResidue(id=...) expects a string
                    if "residue_number" in atom.metadata:
                        atom.metadata["residue_number"] = str(atom.metadata["residue_number"])
                mol_idx += 1
            chain_idx += 1

        # 2. Substrate: Assign next chain letter, residue 1
        if self._n_substrate_molecules > 0:
            chain_id = CHAIN_LETTERS[chain_idx]
            LOGGER.debug(f"Assigning chain {chain_id} to substrate")

            for _ in range(self._n_substrate_molecules):
                mol = self._solvated_topology.molecule(mol_idx)
                for atom in mol.atoms:
                    atom.metadata["chain_id"] = chain_id
                    atom.metadata["residue_number"] = "1"
                mol_idx += 1
            chain_idx += 1

        # 3. Polymers: Assign next chain letter, continue residue numbering across chains
        if self._n_polymer_chains > 0:
            chain_id = CHAIN_LETTERS[chain_idx]
            LOGGER.debug(f"Assigning chain {chain_id} to {self._n_polymer_chains} polymer chain(s)")

            # Track residue number across all polymer chains (continue, don't restart)
            polymer_residue_num = 1

            for _ in range(self._n_polymer_chains):
                mol = self._solvated_topology.molecule(mol_idx)

                # Group atoms by their current residue_number to identify monomers
                # Polymer molecules have per-monomer residue metadata from SDF
                current_monomer_residue = None

                for atom in mol.atoms:
                    atom.metadata["chain_id"] = chain_id

                    # Check if this atom belongs to a new monomer
                    atom_residue = atom.metadata.get("residue_number", "0")
                    if atom_residue != current_monomer_residue:
                        # New monomer - increment our counter
                        if current_monomer_residue is not None:
                            polymer_residue_num += 1
                        current_monomer_residue = atom_residue

                    # Assign the sequential residue number
                    atom.metadata["residue_number"] = str(polymer_residue_num)

                # After processing this polymer chain, increment for next chain's first monomer
                polymer_residue_num += 1
                mol_idx += 1

            chain_idx += 1

        # 4. Solvent: Assign remaining chains with overflow at 9999
        self._assign_solvent_identifiers(
            start_mol_idx=mol_idx,
            start_chain_idx=chain_idx,
            chain_letters=CHAIN_LETTERS,
        )

        LOGGER.info(
            f"PDB identifiers assigned: protein={self._n_enzyme_molecules}, "
            f"substrate={self._n_substrate_molecules}, polymers={self._n_polymer_chains}, "
            f"solvent molecules start at chain {CHAIN_LETTERS[chain_idx]}"
        )

    def _assign_solvent_identifiers(
        self,
        start_mol_idx: int,
        start_chain_idx: int,
        chain_letters: str = "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
        max_residue: int = 9999,
    ) -> None:
        """Assign chain IDs and residue numbers to solvent molecules.

        Solvent molecules (water, ions, co-solvents) are numbered sequentially.
        When residue number exceeds max_residue (9999), overflow to the next
        chain letter and restart residue numbering at 1.

        Args:
            start_mol_idx: Index of first solvent molecule in topology.
            start_chain_idx: Index into chain_letters for first solvent chain.
            chain_letters: Available chain ID characters.
            max_residue: Maximum residue number before overflow (default 9999).
        """
        chain_idx = start_chain_idx
        residue_num = 1

        n_solvent = self._solvated_topology.n_molecules - start_mol_idx
        LOGGER.debug(f"Assigning identifiers to {n_solvent} solvent molecules")

        for mol_idx in range(start_mol_idx, self._solvated_topology.n_molecules):
            # Check for overflow before processing this molecule
            if residue_num > max_residue:
                chain_idx += 1
                residue_num = 1

                if chain_idx >= len(chain_letters):
                    LOGGER.warning(
                        f"Exceeded {len(chain_letters)} chain letters - cycling. "
                        "Consider using a topology format with larger chain ID capacity."
                    )
                    chain_idx = chain_idx % len(chain_letters)

            chain_id = chain_letters[chain_idx]
            mol = self._solvated_topology.molecule(mol_idx)

            for atom in mol.atoms:
                atom.metadata["chain_id"] = chain_id
                atom.metadata["residue_number"] = str(residue_num)

            residue_num += 1

        # Log summary
        n_chains_used = chain_idx - start_chain_idx + 1
        if n_chains_used > 1:
            LOGGER.info(
                f"Solvent required {n_chains_used} chains "
                f"({chain_letters[start_chain_idx]}-{chain_letters[chain_idx]})"
            )

    def save_topology(
        self,
        path: Union[str, Path],
        topology: Optional[Topology] = None,
    ) -> None:
        """Save a topology to PDB file.

        Assigns PDB-compliant chain IDs and residue numbers before writing
        to ensure downstream analysis tools can uniquely identify all atoms.

        Args:
            path: Output path.
            topology: Topology to save (defaults to solvated).
        """
        if topology is None:
            topology = self._solvated_topology

        if topology is None:
            raise RuntimeError("No topology to save")

        # Assign PDB-compliant identifiers before writing
        # This ensures unique (chain_id, residue_number, atom_name) tuples
        if topology is self._solvated_topology:
            self._assign_pdb_identifiers()

        path = Path(path)
        topology.to_file(str(path), keep_ids=True)
        LOGGER.info(f"Saved topology to {path}")

    def build_from_config(
        self,
        config: "SimulationConfig",
        working_dir: Optional[Union[str, Path]] = None,
        polymer_seed: Optional[int] = None,
    ) -> Interchange:
        """Build the complete system from a configuration.

        This is the main entry point for config-driven builds.

        Args:
            config: SimulationConfig with all settings.
            working_dir: Working directory for output files.
            polymer_seed: Random seed for polymer generation. This is used as a
                fallback if config.polymers.random_seed is not set.

        Returns:
            OpenFF Interchange ready for simulation.
        """
        self._working_dir = Path(working_dir) if working_dir else None

        # 1. Build enzyme
        LOGGER.info(f"Building enzyme: {config.enzyme.name}")
        self.build_enzyme(config.enzyme.pdb_path)

        # 2. Build substrate (if configured)
        if config.substrate:
            LOGGER.info(f"Building substrate: {config.substrate.name}")
            self._substrate_builder.build_from_config(config.substrate)
            self._substrate_molecule = self._substrate_builder.molecule
            self._n_substrate_molecules = 1

        # 3. Combine enzyme + substrate
        self.combine_solutes()

        # 4. Build and pack polymers (if configured)
        if config.polymers and config.polymers.enabled:
            LOGGER.info(f"Building polymers: {config.polymers.type_prefix}")

            characters = [m.label for m in config.polymers.monomers]
            probabilities = [m.probability for m in config.polymers.monomers]

            # Determine effective seed: config.random_seed takes precedence over polymer_seed
            effective_seed = config.polymers.random_seed
            if effective_seed is None:
                effective_seed = polymer_seed
            LOGGER.info(f"Using polymer random seed: {effective_seed}")

            # Extract dynamic mode parameters if applicable
            generation_mode = config.polymers.generation_mode.value
            monomer_smiles = None
            monomer_names = None
            residue_names = None
            reactions = None

            if generation_mode == "dynamic":
                # Build monomer_smiles: name -> SMILES
                monomer_smiles = {
                    m.name: m.smiles for m in config.polymers.monomers if m.smiles is not None
                }
                # Build monomer_names: label -> name
                monomer_names = {m.label: m.name for m in config.polymers.monomers}
                # Build residue_names: name -> 3-char residue name
                residue_names = {
                    m.name: m.residue_name
                    for m in config.polymers.monomers
                    if m.residue_name is not None
                }
                # Get reaction config
                reactions = config.polymers.reactions

                LOGGER.info(
                    f"Dynamic mode: monomers={list(monomer_smiles.keys())}, "
                    f"charger={config.polymers.charger.value}"
                )

            self.build_polymers(
                characters=characters,
                probabilities=probabilities,
                length=config.polymers.length,
                count=config.polymers.count,
                type_prefix=config.polymers.type_prefix,
                sdf_directory=config.polymers.sdf_directory,
                seed=effective_seed,
                generation_mode=generation_mode,
                monomer_smiles=monomer_smiles,
                monomer_names=monomer_names,
                residue_names=residue_names,
                reactions=reactions,
                charger_type=config.polymers.charger.value,
                max_retries=config.polymers.max_retries,
                cache_directory=config.polymers.cache_directory,
            )

            # Get packing config (uses defaults if not specified)
            packing = config.polymers.packing
            self.pack_polymers(
                padding=packing.padding,
                tolerance=packing.tolerance,
                movebadrandom=packing.movebadrandom,
                working_directory=self._working_dir,
            )

        # 5. Solvate
        LOGGER.info("Solvating system")
        self._solvent_builder.solvate_from_config(self._combined_topology, config.solvent)
        self._solvated_topology = self._solvent_builder.solvated_topology

        # Save solvated PDB if working dir specified
        if self._working_dir:
            self._working_dir.mkdir(parents=True, exist_ok=True)
            pdb_path = self._working_dir / "solvated_system.pdb"
            self.save_topology(pdb_path)

        # 6. Create Interchange
        use_optimized = config.polymers is not None and config.polymers.enabled
        self.create_interchange(use_optimized_combining=use_optimized)

        return self._interchange

    def get_openmm_components(self) -> Tuple[Any, Any, Any]:
        """Extract OpenMM components from the Interchange.

        Returns:
            Tuple of (topology, system, positions).

        Raises:
            RuntimeError: If Interchange not created.
        """
        if self._interchange is None:
            raise RuntimeError("Interchange not created. Call create_interchange() first.")

        from openff.interchange.interop.openmm._positions import to_openmm_positions

        omm_topology = self._interchange.to_openmm_topology()
        omm_system = self._interchange.to_openmm(combine_nonbonded_forces=False)
        omm_positions = to_openmm_positions(self._interchange, include_virtual_sites=True)

        return omm_topology, omm_system, omm_positions

    def export_to_gromacs(
        self,
        output_dir: Union[str, Path],
        prefix: Optional[str] = None,
        gmx_command: str = "gmx",
        generate_mdps: bool = True,
    ) -> Dict[str, Any]:
        """Export the system to GROMACS format with full simulation setup.

        Generates a complete GROMACS simulation setup including:
        - .gro (coordinates) and .top (topology) files
        - MDP files for energy minimization, equilibration, and production
        - Position restraint files for equilibration stages
        - Run script for executing the full workflow

        The topology is split into separate .itp files for each molecule type
        (monolithic=False), which is cleaner for multi-component systems.

        MDP files are generated from config.yaml parameters to match OpenMM
        simulation settings (temperature, pressure, duration, etc.). OpenFF
        defaults are used for force field parameters (rcoulomb=0.9, rvdw=0.9,
        PME, etc.) to ensure 1:1 parity with OpenMM.

        Args:
            output_dir: Directory to write GROMACS files. Will be created if it
                doesn't exist.
            prefix: Filename prefix for output files. If None, generates a
                descriptive name from the config (e.g., "LipA_EGPMA-SBMA").
            gmx_command: GROMACS command/path for the run script (default "gmx").
            generate_mdps: If True (default), generate MDP files from config.
                If False, only export coordinates and topology.

        Returns:
            Dictionary with paths to generated files:
            - "gro": Path to coordinate file
            - "top": Path to topology file
            - "em_mdp": Path to energy minimization MDP (if generate_mdps=True)
            - "eq_mdps": List of equilibration MDP paths (if generate_mdps=True)
            - "prod_mdp": Path to production MDP (if generate_mdps=True)
            - "posres": Dict of position restraint files (if applicable)
            - "run_script": Path to run script (if generate_mdps=True)

        Raises:
            RuntimeError: If Interchange has not been created.

        Example:
            >>> builder = SystemBuilder.from_config(config)
            >>> builder.build_from_config(config, working_dir)
            >>> result = builder.export_to_gromacs(
            ...     output_dir="gromacs/",
            ...     prefix="my_system"
            ... )
            >>> print(f"Run: cd {result['gro'].parent} && ./{result['run_script'].name}")
        """
        if self._interchange is None:
            raise RuntimeError(
                "Interchange not created. Call create_interchange() or build_from_config() first."
            )

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate prefix from config if not provided
        if prefix is None:
            prefix = self._generate_gromacs_prefix()

        LOGGER.info(f"Exporting to GROMACS format: {output_dir / prefix}.*")

        # Check if we have config for MDP generation
        config = getattr(self, "_config", None)

        if generate_mdps and config is not None:
            # Use GromacsExporter for full export with MDP generation
            from polyzymd.exporters.gromacs import GromacsExporter

            # Get component info for position restraints
            try:
                component_info = self.get_component_info()
            except RuntimeError:
                component_info = None
                LOGGER.warning(
                    "Could not get component info for position restraints. "
                    "Position restraint files will not be generated."
                )

            exporter = GromacsExporter(
                interchange=self._interchange,
                config=config,
                component_info=component_info,
            )

            result = exporter.export(
                output_dir=output_dir,
                prefix=prefix,
                gmx_command=gmx_command,
            )

            return result
        else:
            # Legacy mode: just export coordinates and topology
            return self._export_gromacs_minimal(output_dir, prefix)

    def _export_gromacs_minimal(
        self,
        output_dir: Path,
        prefix: str,
    ) -> Dict[str, Any]:
        """Export minimal GROMACS files (coordinates and topology only).

        This is the legacy export mode used when no config is available
        or when generate_mdps=False.

        Args:
            output_dir: Output directory.
            prefix: Filename prefix.

        Returns:
            Dictionary with gro and top paths.
        """
        # Fix 0-indexed residues (GROMACS requires 1-indexed)
        self._fix_zero_indexed_residues()

        # Export using OpenFF Interchange
        output_prefix = str(output_dir / prefix)
        self._interchange.to_gromacs(
            prefix=output_prefix,
            monolithic=False,
            _merge_atom_types=True,
        )

        gro_path = output_dir / f"{prefix}.gro"
        top_path = output_dir / f"{prefix}.top"
        mdp_path = output_dir / f"{prefix}.mdp"

        LOGGER.info("Generated GROMACS files:")
        LOGGER.info(f"  Coordinates: {gro_path}")
        LOGGER.info(f"  Topology:    {top_path}")
        LOGGER.info(f"  MDP stub:    {mdp_path}")
        LOGGER.info("")
        LOGGER.info("Note: The .mdp file is a stub for single-point energy calculation.")
        LOGGER.info("For MD simulations, use generate_mdps=True with a config file,")
        LOGGER.info("or create custom MDP files.")

        return {
            "gro": gro_path,
            "top": top_path,
            "mdp_stub": mdp_path,
        }

    def _generate_gromacs_prefix(self) -> str:
        """Generate a descriptive filename prefix from the config.

        Creates a prefix in the format: {enzyme_name}_{polymer_prefix}
        Falls back to "system" if config information is not available.

        Returns:
            Filename prefix string (without extension).
        """
        parts = []

        # Get enzyme name from config
        config = getattr(self, "_config", None)
        if config and hasattr(config, "enzyme") and config.enzyme:
            parts.append(config.enzyme.name)

        # Get polymer type prefix from config
        if config and hasattr(config, "polymers") and config.polymers:
            if config.polymers.enabled and config.polymers.type_prefix:
                parts.append(config.polymers.type_prefix)

        if parts:
            return "_".join(parts)
        else:
            return "system"

    def _fix_zero_indexed_residues(self) -> None:
        """Fix 0-indexed residues in the topology for GROMACS compatibility.

        GROMACS requires residue numbers to be 1-indexed. This method checks
        all atoms in the topology and increments any 0-indexed residue numbers.
        """
        if self._interchange is None:
            return

        found_zero_indexed = False

        # First pass: check for 0-indexed residues
        for molecule in self._interchange.topology.molecules:
            for atom in molecule.atoms:
                residue_num = atom.metadata.get("residue_number")
                if residue_num is not None:
                    # Handle both string and int representations
                    if isinstance(residue_num, str):
                        try:
                            if int(residue_num) == 0:
                                found_zero_indexed = True
                                break
                        except ValueError:
                            pass
                    elif residue_num == 0:
                        found_zero_indexed = True
                        break
            if found_zero_indexed:
                break

        if not found_zero_indexed:
            LOGGER.debug("No 0-indexed residues found")
            return

        # Second pass: fix 0-indexed residues
        LOGGER.info("Fixing 0-indexed residues for GROMACS compatibility")
        for molecule in self._interchange.topology.molecules:
            for atom in molecule.atoms:
                residue_num = atom.metadata.get("residue_number")
                if residue_num is not None:
                    if isinstance(residue_num, str):
                        try:
                            atom.metadata["residue_number"] = str(int(residue_num) + 1)
                        except ValueError:
                            pass
                    else:
                        atom.metadata["residue_number"] = residue_num + 1

        LOGGER.info("Fixed all 0-indexed residues")

    def get_component_info(self) -> "SystemComponentInfo":
        """Get system component information for atom group resolution.

        This method returns a SystemComponentInfo dataclass containing the
        atom counts and chain assignments for each system component (protein,
        substrate, polymers, solvent). This information is needed by the
        AtomGroupResolver to resolve predefined atom group names to indices.

        Returns:
            SystemComponentInfo with atom counts and chain assignments

        Raises:
            RuntimeError: If solvated topology not created.
        """
        from polyzymd.core.atom_groups import SystemComponentInfo

        if self._solvated_topology is None:
            raise RuntimeError("No solvated topology. Call solvate() first.")

        # Count atoms per component from topology
        n_protein_atoms = 0
        n_substrate_atoms = 0
        n_polymer_atoms = 0

        mol_idx = 0

        # Protein atoms
        for _ in range(self._n_enzyme_molecules):
            mol = self._solvated_topology.molecule(mol_idx)
            n_protein_atoms += mol.n_atoms
            mol_idx += 1

        # Substrate atoms
        for _ in range(self._n_substrate_molecules):
            mol = self._solvated_topology.molecule(mol_idx)
            n_substrate_atoms += mol.n_atoms
            mol_idx += 1

        # Polymer atoms
        for _ in range(self._n_polymer_chains):
            mol = self._solvated_topology.molecule(mol_idx)
            n_polymer_atoms += mol.n_atoms
            mol_idx += 1

        LOGGER.debug(
            f"Component info: protein={n_protein_atoms} atoms, "
            f"substrate={n_substrate_atoms} atoms, polymer={n_polymer_atoms} atoms"
        )

        return SystemComponentInfo(
            n_protein_atoms=n_protein_atoms,
            n_substrate_atoms=n_substrate_atoms,
            n_polymer_atoms=n_polymer_atoms,
            protein_chain_id="A",
            substrate_chain_id="B",
            polymer_chain_id="C",
            solvent_start_chain_id="D",
        )
