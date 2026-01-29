"""
System builder for orchestrating the complete simulation system construction.

This module coordinates all builders to create a complete, solvated
molecular system and generate the OpenFF Interchange for simulation.
"""

from __future__ import annotations

import logging
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
    ) -> Tuple[List[Molecule], List[int]]:
        """Build polymer components.

        Args:
            characters: Monomer labels.
            probabilities: Monomer selection probabilities.
            length: Monomers per chain.
            count: Number of polymer chains.
            type_prefix: Filename prefix.
            sdf_directory: Directory with pre-built polymer SDFs.
            seed: Random seed for reproducibility.

        Returns:
            Tuple of (unique polymer molecules, counts).
        """
        self._polymer_builder = PolymerBuilder(
            characters=characters,
            probabilities=probabilities,
            length=length,
            type_prefix=type_prefix,
            sdf_directory=sdf_directory,
        )

        molecules, counts = self._polymer_builder.build(count=count, seed=seed)
        self._polymer_molecules = molecules
        self._polymer_counts = counts

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
        padding: float = 1.5,
        working_directory: Optional[Union[str, Path]] = None,
    ) -> Topology:
        """Pack polymers around the combined solute topology.

        Args:
            padding: Box padding in nm.
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

        from openff.interchange.components import _packmol as packmol
        from openff.units import Quantity
        from polymerist.mdtools.openfftools import boxvectors

        LOGGER.info(
            f"Packing {sum(self._polymer_counts)} polymer chains "
            f"({len(self._polymer_molecules)} unique types)"
        )

        # Calculate box vectors
        padding_qty = Quantity(padding, "nanometer")
        min_box_vecs = boxvectors.get_topology_bbox(self._combined_topology)
        box_vecs = boxvectors.pad_box_vectors_uniform(min_box_vecs, padding_qty)

        # Pack polymers
        packed_top = packmol.pack_box(
            self._polymer_molecules,
            self._polymer_counts,
            self._combined_topology,
            box_vectors=box_vecs,
            box_shape=packmol.UNIT_CUBE,
            center_solute="BRICK",
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

    def _create_interchange_simple(self, ff: ForceField, water_mol: Molecule) -> Interchange:
        """Create Interchange using the simple method.

        Args:
            ff: OpenFF ForceField.
            water_mol: Water molecule with pre-computed charges.

        Returns:
            Interchange object.
        """
        charge_from = []
        if self._substrate_molecule:
            charge_from.append(self._substrate_molecule)
        charge_from.append(water_mol)

        # Add co-solvents with pre-computed charges
        if self._solvent_builder._composition:
            for cosolvent in self._solvent_builder._composition.co_solvents:
                if cosolvent.molecule is not None:
                    charge_from.append(cosolvent.molecule)

        return ff.create_interchange(
            self._solvated_topology,
            charge_from_molecules=charge_from,
        )

    def _create_interchange_optimized(self, ff: ForceField, water_mol: Molecule) -> Interchange:
        """Create Interchange using optimized combining.

        This method creates separate Interchanges for each unique molecule
        type and then combines them, which is much faster for systems with
        many polymers.

        Uses pre-computed charges for water and co-solvents to avoid running
        AM1BCC charge calculation for every solvent molecule.

        Args:
            ff: OpenFF ForceField.
            water_mol: Water molecule with pre-computed charges.

        Returns:
            Combined Interchange object.
        """
        LOGGER.info("Using optimized Interchange combining")

        # Build SMILES to molecule mapping for charge templates
        smiles_to_mol: Dict[str, Molecule] = {}

        # Add polymers
        for mol in self._polymer_molecules:
            smiles_to_mol[mol.to_smiles()] = mol

        # Add substrate
        if self._substrate_molecule:
            smiles_to_mol[self._substrate_molecule.to_smiles()] = self._substrate_molecule

        # Add co-solvents with pre-computed charges
        # This ensures DMSO, ethanol, etc. use cached charges instead of AM1BCC
        if self._solvent_builder._composition:
            for cosolvent in self._solvent_builder._composition.co_solvents:
                if cosolvent.molecule is not None:
                    smiles_to_mol[cosolvent.molecule.to_smiles()] = cosolvent.molecule
                    LOGGER.debug(f"Added co-solvent {cosolvent.name} to charge templates")

        # Create interchanges for non-solvent molecules
        all_interchanges: List[Interchange] = []

        for molecule in self._solvated_topology.molecules:
            mol_smiles = molecule.to_smiles()

            # Skip water and ions (handle separately in batch)
            if mol_smiles in ("[H][O][H]", "[Na+]", "[Cl-]"):
                continue

            # Get charge template if available
            charge_from = []
            if mol_smiles in smiles_to_mol:
                charge_from = [smiles_to_mol[mol_smiles]]

            inc = ff.create_interchange(
                molecule.to_topology(),
                charge_from_molecules=charge_from,
            )
            all_interchanges.append(inc)

        # Process water + ions together
        LOGGER.info("Processing water + ions")
        water_ion_mols = []
        for mol in self._solvated_topology.molecules:
            if mol.to_smiles() in ("[H][O][H]", "[Na+]", "[Cl-]"):
                water_ion_mols.append(mol)

        if water_ion_mols:
            water_ion_inc = ff.create_interchange(
                topology=water_ion_mols,
                charge_from_molecules=[water_mol],
            )
            all_interchanges.append(water_ion_inc)

        # Combine all interchanges
        LOGGER.info(f"Combining {len(all_interchanges)} component Interchanges")

        combined = all_interchanges[0]
        for inc in all_interchanges[1:]:
            combined = combined.combine(inc)

        # Reset box vectors
        combined.box = self._solvated_topology.box_vectors

        return combined

    def _renumber_chains(self, topology: Topology) -> None:
        """Re-number chain IDs in the topology.

        Args:
            topology: Topology to modify.
        """
        for i, mol in enumerate(topology.molecules):
            for atom in mol.atoms:
                atom.metadata["chain_id"] = str(i + 1)

    def save_topology(
        self,
        path: Union[str, Path],
        topology: Optional[Topology] = None,
    ) -> None:
        """Save a topology to PDB file.

        Args:
            path: Output path.
            topology: Topology to save (defaults to solvated).
        """
        if topology is None:
            topology = self._solvated_topology

        if topology is None:
            raise RuntimeError("No topology to save")

        path = Path(path)
        topology.to_file(str(path))
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
            polymer_seed: Random seed for polymer generation.

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

        # 3. Combine enzyme + substrate
        self.combine_solutes()

        # 4. Build and pack polymers (if configured)
        if config.polymers and config.polymers.enabled:
            LOGGER.info(f"Building polymers: {config.polymers.type_prefix}")

            characters = [m.label for m in config.polymers.monomers]
            probabilities = [m.probability for m in config.polymers.monomers]

            self.build_polymers(
                characters=characters,
                probabilities=probabilities,
                length=config.polymers.length,
                count=config.polymers.count,
                type_prefix=config.polymers.type_prefix,
                sdf_directory=config.polymers.sdf_directory,
                seed=polymer_seed,
            )

            self.pack_polymers(
                padding=config.solvent.box.padding,
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
