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

if TYPE_CHECKING:
    from openff.interchange import Interchange
    from openff.toolkit import ForceField, Molecule, Topology

    from polyzymd.config.schema import SimulationConfig

from polyzymd.builders.enzyme import EnzymeBuilder
from polyzymd.builders.polymer import PolymerBuilder
from polyzymd.builders.solvent import SolventBuilder, SolventComposition
from polyzymd.builders.substrate import SubstrateBuilder

LOGGER = logging.getLogger(__name__)


class SystemBuilder:
    """Orchestrator for building complete simulation systems.

    This class coordinates enzyme loading, substrate charging, polymer
    generation, solvation, and Interchange creation through the canonical
    single-call OpenFF pathway.

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
        self._preserve_enzyme_chain_ids: bool = False

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
        characters: List[str] | None = None,
        probabilities: List[float] | None = None,
        length: int | None = None,
        count: int | None = None,
        type_prefix: str | None = None,
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
        fragments: Optional[Dict[str, Any]] = None,
        provided_molecules: Optional[List[Any]] = None,
        polymer_random_seed: Optional[int] = None,
    ) -> Tuple[List[Molecule], List[int]]:
        """Build polymer components.

        Supports three sequence-derived generation modes and one provided-only mode:

        - "cached": Load legacy sequence-derived polymer SDFs from sdf_directory
        - "dynamic": Generate polymers on-the-fly from monomer SMILES
        - "fragments": Assemble explicit terminal/middle fragments natively
        - "provided": Pack only explicit charged SDF molecule pools

        ``provided_molecules`` adds opaque user-provided charged SDF molecules
        to the merged molecule/count list before packing in generated modes.

        Args:
            characters: Monomer labels.
            probabilities: Monomer selection probabilities.
            length: Monomers per chain.
            count: Number of polymer chains.
            type_prefix: Filename prefix.
            sdf_directory: Deprecated directory with sequence-derived polymer SDFs.
            seed: Random seed for reproducibility.
            generation_mode: "cached", "dynamic", "fragments", or "provided".
            monomer_smiles: Dict of monomer name -> SMILES (dynamic mode).
            monomer_names: Dict of label -> monomer name (dynamic mode).
            residue_names: Dict of monomer name -> 3-char residue name.
            reactions: ReactionConfig with literal "default" markers for native dynamic mode;
                custom .rxn files are unsupported.
            charger_type: Charge method ("nagl", "espaloma", "am1bcc").
            max_retries: Max retries for polymer generation.
            cache_directory: Directory for caching generated polymers.
            fragments: Explicit terminal/middle fragment specs for fragments mode.
            provided_molecules: Additive opaque user-provided charged SDF molecule pools.
            polymer_random_seed: Polymer-level seed for provided molecule precedence.

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
            fragments=fragments,
            provided_molecules=provided_molecules,
            polymer_random_seed=polymer_random_seed,
        )

        molecules, counts = self._polymer_builder.build(count=count or 0, seed=seed)
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

        from openff.toolkit import Topology

        # Start with enzyme
        molecules = [self._enzyme_topology.molecule(0)]

        # Add substrate if present
        if self._substrate_molecule is not None:
            molecules.append(self._substrate_molecule)

        # Create combined topology
        self._combined_topology = Topology.from_molecules(molecules)

        # Re-number legacy solutes. Pablo-ingested conjugates already carry
        # protein/polymer chain IDs that must survive packing and solvation.
        if not self._preserve_enzyme_chain_ids:
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
        box_vectors_nm: Optional[List[float]] = None,
    ) -> Topology:
        """Pack polymers around the combined solute topology.

        Args:
            padding: Box padding in nm. Larger values give polymers more room
                and can significantly speed up PACKMOL convergence.  Ignored
                when *box_vectors_nm* is provided.
            tolerance: PACKMOL tolerance in Angstrom.
            movebadrandom: When True, pass the ``movebadrandom`` keyword to
                PACKMOL. Improves convergence for dense or heterogeneous
                polymer systems (many unique chain types) by placing
                badly-packed molecules at random positions in the box.
            working_directory: Directory for PACKMOL files.
            box_vectors_nm: Optional explicit box dimensions ``[Lx, Ly, Lz]``
                in nanometers.  When provided, overrides the auto-computed
                bounding box + *padding*.  The protein is centered at the
                midpoint of this box.

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

        import numpy as np
        from openff.units import Quantity

        from polyzymd.utils import boxvectors
        from polyzymd.utils.packmol import pack_polymers as _pack_polymers

        # Calculate box vectors — explicit override or auto from bbox + padding
        if box_vectors_nm is not None:
            LOGGER.info(
                f"Using explicit box vectors: "
                f"[{box_vectors_nm[0]:.2f}, {box_vectors_nm[1]:.2f}, "
                f"{box_vectors_nm[2]:.2f}] nm"
            )
            box_vecs = Quantity(np.diag(box_vectors_nm), "nanometer")
        else:
            padding_qty = Quantity(padding, "nanometer")
            min_box_vecs = boxvectors.get_topology_bbox(self._combined_topology)
            box_vecs = boxvectors.pad_box_vectors_uniform(min_box_vecs, padding_qty)

        LOGGER.info(
            f"Packing {sum(self._polymer_counts)} polymer chains "
            f"({len(self._polymer_molecules)} unique types), "
            f"tolerance={tolerance} A" + (" [movebadrandom]" if movebadrandom else "")
        )

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

        # Re-number legacy systems only. Conjugate workflows carry mixed
        # protein/polymer chain IDs inside the first molecule and must preserve
        # that metadata through packing.
        if not self._preserve_enzyme_chain_ids:
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

    def create_interchange(self) -> Interchange:
        """Create the OpenFF Interchange for simulation.

        By default, passes the entire solvated topology to a single
        ``ForceField.create_interchange()`` call.  OpenFF's internal
        ``identical_molecule_groups`` mechanism deduplicates SMIRKS
        matching and charge assignment so that work scales with the
        number of *unique* molecule types, not total molecules.

        Pre-computed charges are supplied via ``charge_from_molecules``
        for water, polymers, substrate, and co-solvents so that no
        AM1BCC calculations are triggered at runtime.

        Returns:
            OpenFF Interchange object.

        Raises:
            RuntimeError: If system not solvated.
        """
        if self._solvated_topology is None:
            raise RuntimeError("System must be solvated before creating Interchange")

        from openff.toolkit import ForceField

        from polyzymd.data.solvent_molecules import get_solvent_molecule

        LOGGER.info("Creating Interchange")

        ff = ForceField(self._protein_ff, self._sm_ff)

        # Get water molecule with pre-computed charges for charge_from_molecules
        # Use PolyzyMD's cached water for consistent residue metadata
        water_model = "tip3p"
        if self._solvent_builder._composition:
            water_model = self._solvent_builder._composition.water_model
        water_mol = get_solvent_molecule(water_model)

        LOGGER.info("Using single-call Interchange creation")
        self._interchange = self._create_interchange_single_call(ff, water_mol)

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

    def collect_standard_charge_templates(self) -> tuple[Any, ...]:
        """Collect charged templates for non-conjugate molecules.

        Returns
        -------
        tuple[Any, ...]
            Charged templates for water, free polymers, substrate, and
            co-solvents in first-seen order.
        """
        from polyzymd.data.solvent_molecules import get_solvent_molecule

        water_model = "tip3p"
        if self._solvent_builder._composition:
            water_model = self._solvent_builder._composition.water_model
        water_mol = get_solvent_molecule(water_model)

        charge_from: list[Any] = [water_mol]
        seen_smiles: set[str] = set()
        for mol in self._polymer_molecules:
            smi = mol.to_smiles()
            if smi in seen_smiles:
                continue
            charge_from.append(mol)
            seen_smiles.add(smi)

        if self._substrate_molecule:
            charge_from.append(self._substrate_molecule)

        if self._solvent_builder._composition:
            for cosolvent in self._solvent_builder._composition.co_solvents:
                if cosolvent.molecule is not None:
                    charge_from.append(cosolvent.molecule)

        return tuple(charge_from)

    def _create_interchange_single_call(self, ff: ForceField, water_mol: Molecule) -> Interchange:
        """Create Interchange with a single ``ForceField.create_interchange()`` call.

        Passes the entire solvated topology at once.  OpenFF internally
        groups identical molecules (``Topology.identical_molecule_groups``)
        so SMIRKS matching and charge assignment scale with the number of
        *unique* molecule types, not total molecules.

        This avoids ``Interchange.combine()`` entirely — previously the
        dominant bottleneck (84 min for a 308 K-atom system).

        Args:
            ff: OpenFF ForceField.
            water_mol: Water molecule with pre-computed charges.

        Returns:
            Interchange object with box vectors set.
        """
        # Build charge templates: one per unique species
        t0 = time.perf_counter()
        charge_from = list(self.collect_standard_charge_templates())
        polymer_template_count = max(len(charge_from) - 1, 0)

        LOGGER.info(
            f"Charge templates: {len(charge_from)} "
            f"(water + {polymer_template_count} standard non-water template(s))"
        )
        LOGGER.info(f"  Template prep took {time.perf_counter() - t0:.1f}s")

        # Single parameterization call
        t1 = time.perf_counter()
        LOGGER.info(
            f"Calling ff.create_interchange() on full topology "
            f"({self._solvated_topology.n_molecules} molecules, "
            f"{self._solvated_topology.n_atoms} atoms) ..."
        )

        # Suppress per-atom INFO logging from OpenFF's LibraryCharges handler.
        # It emits one line per atom during charge assignment (~200K lines for
        # a typical solvated system), which floods logs without adding value.
        _nonbonded_logger = logging.getLogger("openff.interchange.smirnoff._nonbonded")
        _prev_level = _nonbonded_logger.level
        _nonbonded_logger.setLevel(logging.WARNING)
        try:
            interchange = ff.create_interchange(
                self._solvated_topology,
                charge_from_molecules=charge_from,
            )
        finally:
            _nonbonded_logger.setLevel(_prev_level)

        t2 = time.perf_counter()
        LOGGER.info(f"  create_interchange() completed in {t2 - t1:.1f}s")

        # Set box vectors
        interchange.box = self._solvated_topology.box_vectors

        return interchange

    def _renumber_chains(self, topology: Topology) -> None:
        """Normalize chain IDs before canonical PDB identifier assignment.

        This lightweight normalization keeps intermediate OpenFF metadata
        deterministic before ``_assign_pdb_identifiers()`` applies the
        canonical PolyzyMD chain convention during topology export.

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

        Chain assignment uses FIXED letters regardless of component presence:
        - Chain A: Protein (preserves original residue numbers from input PDB)
        - Chain B: Substrate (residue 1; letter reserved even if no substrate)
        - Chain C: Polymers (preserves per-monomer residue numbers)
        - Chain D+: Solvent (overflow at 9999 residues per chain)

        Using fixed letters ensures consistency with SystemComponentInfo and
        AtomGroupResolver, which hardcode A=protein, B=substrate, C=polymer.
        Without this, absent components (e.g., no substrate) would shift later
        chain letters, causing restraints to target wrong atoms.

        This ensures every atom can be uniquely identified by the tuple
        (chain_id, residue_number, residue_name, atom_name) for downstream
        analysis tools like MDAnalysis and PyMOL.

        Raises:
            RuntimeError: If no solvated topology exists.
        """
        if self._solvated_topology is None:
            raise RuntimeError("No solvated topology. Call solvate() first.")

        CHAIN_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        mol_idx = 0

        # Fixed chain letter assignments per component type.
        # A=protein, B=substrate, C=polymer, D+=solvent — regardless of whether
        # a component is present. This ensures downstream code (SystemComponentInfo,
        # AtomGroupResolver, from_topology()) always sees the expected chain IDs.
        PROTEIN_CHAIN = "A"
        SUBSTRATE_CHAIN = "B"
        POLYMER_CHAIN = "C"
        SOLVENT_START_IDX = 3  # index of 'D' in CHAIN_LETTERS

        # 1. Protein: Always chain A, preserve original residue numbers
        if self._n_enzyme_molecules > 0:
            LOGGER.debug(f"Assigning chain {PROTEIN_CHAIN} to protein")

            for _ in range(self._n_enzyme_molecules):
                mol = self._solvated_topology.molecule(mol_idx)
                for atom in mol.atoms:
                    if self._preserve_enzyme_chain_ids:
                        atom.metadata.setdefault("chain_id", PROTEIN_CHAIN)
                    else:
                        atom.metadata["chain_id"] = PROTEIN_CHAIN
                    # Ensure residue_number is a string (PDB loader may store as int)
                    # OpenMM's addResidue(id=...) expects a string
                    if "residue_number" in atom.metadata:
                        atom.metadata["residue_number"] = str(atom.metadata["residue_number"])
                mol_idx += 1

        # 2. Substrate: Always chain B, residue 1
        if self._n_substrate_molecules > 0:
            LOGGER.debug(f"Assigning chain {SUBSTRATE_CHAIN} to substrate")

            for _ in range(self._n_substrate_molecules):
                mol = self._solvated_topology.molecule(mol_idx)
                for atom in mol.atoms:
                    atom.metadata["chain_id"] = SUBSTRATE_CHAIN
                    atom.metadata["residue_number"] = "1"
                mol_idx += 1

        # 3. Polymers: Always chain C, continue residue numbering across chains
        if self._n_polymer_chains > 0:
            LOGGER.debug(
                f"Assigning chain {POLYMER_CHAIN} to {self._n_polymer_chains} polymer chain(s)"
            )

            # Track residue number across all polymer chains (continue, don't restart)
            polymer_residue_num = self._next_residue_number_for_chain(
                POLYMER_CHAIN,
                end_mol_idx=mol_idx,
            )

            for _ in range(self._n_polymer_chains):
                mol = self._solvated_topology.molecule(mol_idx)

                # Group atoms by their current residue_number to identify monomers
                # Polymer molecules have per-monomer residue metadata from SDF
                current_monomer_residue = None

                for atom in mol.atoms:
                    atom.metadata["chain_id"] = POLYMER_CHAIN

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

        # 4. Solvent: Always starts at chain D (index 3)
        self._assign_solvent_identifiers(
            start_mol_idx=mol_idx,
            start_chain_idx=SOLVENT_START_IDX,
            chain_letters=CHAIN_LETTERS,
        )

        LOGGER.info(
            f"PDB identifiers assigned: protein={self._n_enzyme_molecules}, "
            f"substrate={self._n_substrate_molecules}, polymers={self._n_polymer_chains}, "
            f"solvent molecules start at chain {CHAIN_LETTERS[SOLVENT_START_IDX]}"
        )

    def _next_residue_number_for_chain(self, chain_id: str, *, end_mol_idx: int) -> int:
        """Return the next residue number after preserved solute residues on a chain."""
        if self._solvated_topology is None or not self._preserve_enzyme_chain_ids:
            return 1

        max_residue_number = 0
        for mol_idx in range(end_mol_idx):
            mol = self._solvated_topology.molecule(mol_idx)
            for atom in mol.atoms:
                if str(atom.metadata.get("chain_id", "")).upper() != chain_id:
                    continue
                residue_number = _metadata_residue_number(atom.metadata)
                if residue_number is not None:
                    max_residue_number = max(max_residue_number, residue_number)
        return max_residue_number + 1

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
                fragments=config.polymers.fragments,
                provided_molecules=config.polymers.provided_molecules,
                polymer_random_seed=config.polymers.random_seed,
            )

            # Get packing config (uses defaults if not specified)
            packing = config.polymers.packing
            self.pack_polymers(
                padding=packing.padding,
                tolerance=packing.tolerance,
                movebadrandom=packing.movebadrandom,
                working_directory=self._working_dir,
                box_vectors_nm=packing.box_vectors,
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
        self.create_interchange()

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

        t0 = time.perf_counter()
        omm_topology = self._interchange.to_openmm_topology()
        t_topo = time.perf_counter() - t0
        LOGGER.info("  to_openmm_topology: %.1fs", t_topo)

        t0 = time.perf_counter()
        omm_system = self._interchange.to_openmm(combine_nonbonded_forces=True)
        t_sys = time.perf_counter() - t0
        LOGGER.info("  to_openmm (system): %.1fs", t_sys)

        t0 = time.perf_counter()
        omm_positions = to_openmm_positions(self._interchange, include_virtual_sites=True)
        t_pos = time.perf_counter() - t0
        LOGGER.info("  to_openmm_positions: %.1fs", t_pos)

        LOGGER.info(
            "  OpenMM extraction total: %.1fs (topo=%.1f, sys=%.1f, pos=%.1f)",
            t_topo + t_sys + t_pos,
            t_topo,
            t_sys,
            t_pos,
        )

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

        if generate_mdps and config is None:
            raise RuntimeError(
                "GROMACS MDP generation requires a validated SimulationConfig. "
                "Pass generate_mdps=False to export coordinates and topology only."
            )

        if generate_mdps:
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
        # Build-only export without MDP generation
        return self._export_gromacs_minimal(output_dir, prefix)

    def _export_gromacs_minimal(
        self,
        output_dir: Path,
        prefix: str,
    ) -> Dict[str, Any]:
        """Export GROMACS coordinates and topology only.

        This build-only mode lets users export ``.gro``, ``.top``, and
        ``.itp`` files and then manage downstream GROMACS inputs themselves.

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


def _metadata_residue_number(metadata: dict[str, Any]) -> int | None:
    """Return an integer residue number from OpenFF atom metadata when present."""
    for key in ("residue_number", "res_seq"):
        value = metadata.get(key)
        if value is None:
            continue
        try:
            return int(str(value).strip())
        except ValueError:
            continue
    return None
