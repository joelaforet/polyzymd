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

    from polyzymd.config.schema import ProteinSolute, SimulationConfig
    from polyzymd.core.components import ChargeSource, ComponentRecord

from polyzymd.builders.enzyme import EnzymeBuilder
from polyzymd.builders.modified_protein import ModifiedProteinBuild
from polyzymd.builders.polymer import PolymerBuilder
from polyzymd.builders.solvent import SolventBuilder, SolventComposition
from polyzymd.builders.substrate import SubstrateBuilder

LOGGER = logging.getLogger(__name__)

_ORIGINAL_RESIDUE_TOKEN_KEY = "_polyzymd_original_residue_token"


def _default_nagl_model() -> str:
    from polyzymd.utils.charging import DEFAULT_NAGL_MODEL

    return DEFAULT_NAGL_MODEL


def _cosolvent_charges(cosolvent: Any, nagl_model: Optional[str]) -> "ChargeSource":
    """Return where a co-solvent's charges come from: its bundled SDF, or its charge method."""
    from polyzymd.core.components import ChargeSource
    from polyzymd.data.solvent_molecules import _SOLVENTS_DIR, is_bundled_solvent

    library = _SOLVENTS_DIR / f"{cosolvent.name.lower().strip()}.sdf"
    method = getattr(cosolvent, "charge_method", None)
    bundled = is_bundled_solvent(cosolvent.name) and library.exists()
    if method in (None, "preset") and bundled and cosolvent.molecule is not None:
        from openff.toolkit import Molecule

        # The solvent loader uses the bundled charges whenever the molecule is the
        # library's, whether or not the config spelled out its SMILES.
        if Molecule.from_file(str(library)).is_isomorphic_with(cosolvent.molecule):
            return ChargeSource("preset", source_file=str(library))
    method = method if method not in (None, "preset") else "nagl"
    model = (nagl_model or _default_nagl_model()) if method == "nagl" else None
    return ChargeSource(method, model=model)


def _optional_str(value: Any) -> Optional[str]:
    return None if value is None else str(value)


def _cosolvent_smiles(cosolvent: Any) -> str:
    """The SMILES of the co-solvent molecule the build packed, without counter-ions."""
    from polyzymd.data.solvent_molecules import split_counter_ions

    if cosolvent.smiles:
        return split_counter_ions(cosolvent.smiles, cosolvent.name)[0]
    return cosolvent.molecule.to_smiles()


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
        nagl_model: Optional[str] = None,
    ) -> None:
        """Initialize the SystemBuilder.

        Args:
            protein_forcefield: Force field for proteins.
            small_molecule_forcefield: Force field for small molecules.
            nagl_model: NAGL model for every molecule charged with nagl; the
                default model when None.
        """
        self._protein_ff = protein_forcefield
        self._sm_ff = small_molecule_forcefield
        self._nagl_model = nagl_model

        # Sub-builders
        self._enzyme_builder = EnzymeBuilder()
        self._substrate_builder = SubstrateBuilder()
        self._polymer_builder: Optional[PolymerBuilder] = None
        self._solvent_builder = SolventBuilder()

        # State
        self._enzyme_topology: Optional[Topology] = None
        self._substrate_molecule: Optional[Molecule] = None
        self._substrate_charge_method: Optional[str] = None
        self._substrate_conformer_index: Optional[int] = None
        #: (name, segid) of each solute type, from the config's solutes.
        self._solute_handles: Dict[str, Tuple[str, str]] = {}
        self._polymer_molecules: List[Molecule] = []
        self._polymer_counts: List[int] = []
        self._combined_topology: Optional[Topology] = None
        self._solvated_topology: Optional[Topology] = None
        self._interchange: Optional[Interchange] = None
        self._working_dir: Optional[Path] = None
        self._build_provenance: Dict[str, Any] = {}
        #: The protein after its mutations and modifications, when it has any.
        self._modified_protein: Optional[ModifiedProteinBuild] = None

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
            nagl_model=config.force_field.nagl_model,
        )
        builder._config = config
        return builder

    def build_enzyme(
        self,
        pdb_path: Union[str, Path],
        custom_substructures_path: Union[str, Path, None] = None,
    ) -> Topology:
        """Build the enzyme component.

        Parameters
        ----------
        pdb_path : str or Path
            Path to enzyme PDB file.
        custom_substructures_path : str or Path, optional
            JSON file of residue templates for residues OpenFF does not know.

        Returns
        -------
        Topology
            Enzyme topology.

        Raises
        ------
        RuntimeError
            If OpenFF loads the enzyme PDB without any molecules.
        """
        self._enzyme_topology = self._enzyme_builder.build(pdb_path, custom_substructures_path)
        self._n_enzyme_molecules = self._enzyme_topology.n_molecules
        if self._n_enzyme_molecules <= 0:
            raise RuntimeError("OpenFF enzyme topology contains no molecules")
        if self._n_enzyme_molecules > 1:
            LOGGER.info(
                "Enzyme topology contains %d protein molecules; retaining all on chain A",
                self._n_enzyme_molecules,
            )
        return self._enzyme_topology

    def build_modified_enzyme(
        self, solute: "ProteinSolute", base_dir: Union[str, Path, None] = None
    ) -> Topology:
        """Build a protein with mutations or modifications, parameterized site by site.

        Every change is applied (and so checked against the structure) before
        anything is packed. The protein is written to the run folder as
        ``modified_protein.pdb``.
        """
        from openff.toolkit import Topology

        self._modified_protein = ModifiedProteinBuild.build(
            solute,
            constraints_from=self._sm_ff,
            working_dir=self._working_dir,
            base_dir=base_dir,
        )
        molecules = self._modified_protein.molecules()
        self._enzyme_topology = Topology.from_molecules(molecules)
        self._n_enzyme_molecules = len(molecules)
        self._build_provenance["protein_modifications"] = self._modified_protein.provenance()
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
            nagl_model=self._nagl_model,
        )
        self._substrate_charge_method = charge_method
        self._substrate_conformer_index = conformer_index
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
        end_groups: Optional[Dict[str, str]] = None,
        force_field: Optional[List[str]] = None,
        charger_type: str = "nagl",
        max_retries: int = 10,
        cache_directory: Optional[Union[str, Path]] = None,
    ) -> Tuple[List[Molecule], List[int]]:
        """Build polymer components.

        Supports two generation modes:
        - "cached": Load pre-built polymer SDFs from sdf_directory
        - "dynamic": Build chains with mBuild from star-sited monomer SMILES

        Args:
            characters: Monomer labels.
            probabilities: Monomer selection probabilities.
            length: Monomers per chain.
            count: Number of polymer chains.
            type_prefix: Filename prefix.
            sdf_directory: Directory with pre-built polymer SDFs (cached mode).
            seed: Random seed for reproducibility.
            generation_mode: "cached" or "dynamic".
            monomer_smiles: Dict of monomer name -> star-sited repeat-unit SMILES (dynamic mode).
            monomer_names: Dict of label -> monomer name (dynamic mode).
            residue_names: Dict of monomer name -> 3-char residue name.
            end_groups: {"head": smiles, "tail": smiles}, the star-sited caps (dynamic mode).
            force_field: Force field files the chains are relaxed with (dynamic mode).
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
            end_groups=end_groups,
            force_field=force_field,
            charger_type=charger_type,
            max_retries=max_retries,
            nagl_model=self._nagl_model,
        )

        molecules, counts = self._polymer_builder.build(count=count, seed=seed)
        self._polymer_molecules = molecules
        self._polymer_counts = counts
        self._n_polymer_chains = sum(counts)

        return molecules, counts

    def _polymer_chain_builders(self) -> List[Dict[str, Any]]:
        """Return the builder records of the polymer chains (mBuild), if any built them."""
        builder = self._polymer_builder
        return builder.chain_builders() if builder is not None else []

    def combine_solutes(self) -> Topology:
        """Combine enzyme, substrate, and polymers into a single topology.

        All OpenFF enzyme topology molecules are retained in their original
        order before substrate and polymer components are added. During export,
        every enzyme molecule is assigned to PolyzyMD protein chain ``A``.

        Returns
        -------
        Topology
            Combined topology ready for solvation.

        Raises
        ------
        RuntimeError
            If enzyme has not been built or contains no molecules.
        """
        if self._enzyme_topology is None:
            raise RuntimeError("Enzyme must be built before combining solutes")

        LOGGER.info("Combining solute components")

        from openff.toolkit import Topology

        enzyme_molecule_count = self._enzyme_topology.n_molecules
        if enzyme_molecule_count <= 0:
            raise RuntimeError("OpenFF enzyme topology contains no molecules")
        self._n_enzyme_molecules = enzyme_molecule_count

        # Retain all enzyme molecules before substrate and polymers
        molecules = [self._enzyme_topology.molecule(i) for i in range(self._n_enzyme_molecules)]

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
        box_vectors: Any,
        padding: float = 2.0,
        tolerance: float = 2.0,
        movebadrandom: bool = False,
        working_directory: Optional[Union[str, Path]] = None,
        seed: Optional[int] = None,
        confine_to_sphere: bool = True,
        nloop: int = 200,
    ) -> Topology:
        """Pack polymers around the combined solute topology.

        Args:
            box_vectors: Final periodic box vectors (3x3 Quantity) of the
                simulation cell.  Chains are packed inside the rectangular
                brick of this cell so that no atom can overlap its own
                periodic image.
            padding: Padding in nm of the confinement sphere around the solute.
            tolerance: PACKMOL tolerance in Angstrom.
            movebadrandom: When True, pass the ``movebadrandom`` keyword to
                PACKMOL. Improves convergence for dense or heterogeneous
                polymer systems (many unique chain types) by placing
                badly-packed molecules at random positions in the box.
            working_directory: Directory for PACKMOL files.
            seed: Packmol random seed (typically the replicate index).
                ``None`` leaves Packmol on its fixed built-in default.
            confine_to_sphere: Confine chains to a sphere around the solute
                (radius = solute bounding-box circumradius + *padding*) while
                packing inside the final periodic brick. Default ``True``.
            nloop: Maximum PACKMOL GENCAN loops per molecule type.

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

        from polyzymd.utils.packmol import pack_polymers as _pack_polymers
        from polyzymd.utils.packmol import solute_sphere_constraint

        LOGGER.info(
            "Packing polymers inside the final periodic brick "
            f"(box vectors {np.diagonal(np.asarray(box_vectors.m_as('nanometer'))).round(3)} nm "
            "on the diagonal)"
        )
        LOGGER.info(
            f"Packing {sum(self._polymer_counts)} polymer chains "
            f"({len(self._polymer_molecules)} unique types), "
            f"tolerance={tolerance} A" + (" [movebadrandom]" if movebadrandom else "")
        )

        # Pack polymers using our custom PACKMOL runner (supports movebadrandom)
        solute_topology = self._combined_topology
        packed_top = _pack_polymers(
            molecules=self._polymer_molecules,
            number_of_copies=self._polymer_counts,
            solute=self._combined_topology,
            box_vectors=box_vectors,
            tolerance_angstrom=tolerance,
            movebadrandom=movebadrandom,
            nloop=nloop,
            seed=seed,
            confine_to_sphere=confine_to_sphere,
            sphere_padding_angstrom=padding * 10.0,
            working_directory=str(working_directory) if working_directory else None,
            retain_working_files=True,
        )
        self._build_provenance["polymer_packmol_seed"] = seed
        if confine_to_sphere:
            sphere = solute_sphere_constraint(solute_topology, padding_angstrom=padding * 10.0)
            self._build_provenance["polymer_sphere_radius_nm"] = float(sphere[3]) / 10.0

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

        if self._nagl_model is not None:
            import dataclasses

            composition = dataclasses.replace(
                composition or SolventComposition(), nagl_model=self._nagl_model
            )
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

        # Canonical metadata must be present before OpenFF serializes topology
        # identifiers into downstream coordinate and topology exports
        self._assign_pdb_identifiers()

        # Loaded last, the small-molecule force field gives every non-protein
        # molecule its parameters, and its generic SMIRKS also override every
        # ff14SB bonded and Lennard-Jones term of the protein; only the ff14SB
        # library charges survive. _apply_protein_force_field restores the
        # protein terms. Loading ff14SB last instead is no fix: its cap and
        # amide SMIRKS then match acetone, acetates and amide polymers.
        ff = ForceField(self._protein_ff, self._sm_ff)

        # Get water molecule with pre-computed charges for charge_from_molecules
        # Use PolyzyMD's cached water instead of Polymerist's for consistency
        water_model = "tip3p"
        if self._solvent_builder._composition:
            water_model = self._solvent_builder._composition.water_model
        water_mol = get_solvent_molecule(water_model)

        LOGGER.info("Using single-call Interchange creation")
        self._interchange = self._create_interchange_single_call(ff, water_mol)
        self._apply_protein_force_field(self._interchange)

        LOGGER.info("Interchange created successfully")

        return self._interchange

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
        charge_from: List[Molecule] = [water_mol]
        if self._modified_protein is not None:
            # A modified protein has no ff14SB library charges for its sites;
            # its charges are the ones parameterization assigned, site by site.
            from polyzymd.modifying_proteins.parameters import charge_templates

            charge_from.extend(charge_templates(list(self._enzyme_topology.molecules)))

        seen_smiles: set = set()
        for mol in self._polymer_molecules:
            smi = mol.to_smiles()
            if smi not in seen_smiles:
                charge_from.append(mol)
                seen_smiles.add(smi)

        parts = ["water"]
        if seen_smiles:
            parts.append(f"{len(seen_smiles)} polymer type(s)")

        if self._substrate_molecule:
            charge_from.append(self._substrate_molecule)
            parts.append("substrate")

        if self._solvent_builder._composition:
            for cosolvent in self._solvent_builder._composition.co_solvents:
                if cosolvent.molecule is not None:
                    charge_from.append(cosolvent.molecule)
                    parts.append(cosolvent.name)

        LOGGER.info(f"Charge templates: {len(charge_from)} ({' + '.join(parts)})")
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

    def _apply_protein_force_field(self, interchange: Interchange) -> None:
        """Replace the protein's terms with those of the protein force field alone.

        The protein molecules come first in the solvated topology. They are
        parameterized on their own with the protein force field, and every
        bonded, Lennard-Jones and constraint term wholly inside them is swapped
        for its counterpart, so the protein no longer depends on load order.
        Terms of other molecules and all charges stay as they are; the
        protein's charges are already the ff14SB library charges. Constraints
        follow the small-molecule force field's bonds-to-hydrogen policy, with
        the protein force field's bond lengths.

        Impropers are not checked for coverage, because the two force fields
        define them differently; every other handler must cover the same
        protein atoms as the combined force field did.

        Args:
            interchange: Interchange of the whole solvated topology, changed
                in place.

        Raises:
            RuntimeError: If the protein force field does not cover a protein
                term the combined force field covered, the global nonbonded
                settings of the two force fields differ, or the first
                molecules of the solvated topology are not the enzyme.
        """
        from openff.toolkit import ForceField, Topology

        if self._n_enzyme_molecules == 0:
            LOGGER.warning("No enzyme molecules recorded; the protein force field is not applied")
            return
        protein_atoms = sum(
            self._solvated_topology.molecule(index).n_atoms
            for index in range(self._n_enzyme_molecules)
        )
        if self._enzyme_topology is not None and protein_atoms != self._enzyme_topology.n_atoms:
            raise RuntimeError(
                f"The first {self._n_enzyme_molecules} molecule(s) of the solvated topology hold "
                f"{protein_atoms} atoms, the enzyme {self._enzyme_topology.n_atoms}"
            )

        t0 = time.perf_counter()
        if self._modified_protein is not None:
            # ff14SB outside the sites, each site from its parameter set, as
            # parameterization made them; atoms are in the same order.
            reference = self._modified_protein.interchange
        else:
            protein_ff = ForceField(self._protein_ff)
            small_molecule_ff = ForceField(self._sm_ff)
            # The Amber port has no Constraints handler.
            if "Constraints" not in protein_ff.registered_parameter_handlers:
                protein_ff.register_parameter_handler(small_molecule_ff["Constraints"])
            protein_topology = Topology.from_molecules(
                [
                    self._solvated_topology.molecule(index)
                    for index in range(self._n_enzyme_molecules)
                ]
            )
            reference = protein_ff.create_interchange(protein_topology)

        for name in ("vdW", "Electrostatics"):
            for attribute in ("scale_12", "scale_13", "scale_14", "scale_15"):
                if (
                    abs(getattr(interchange[name], attribute) - getattr(reference[name], attribute))
                    > 1e-6
                ):
                    raise RuntimeError(
                        f"{self._protein_ff} and {self._sm_ff} use different {name} {attribute}"
                    )
        if interchange["vdW"].mixing_rule != reference["vdW"].mixing_rule:
            raise RuntimeError(f"{self._protein_ff} and {self._sm_ff} use different mixing rules")

        for name in ("Bonds", "Angles", "ProperTorsions", "ImproperTorsions", "vdW", "Constraints"):
            target = interchange[name]
            source = reference[name]
            removed = [key for key in target.key_map if max(key.atom_indices) < protein_atoms]
            for key in removed:
                del target.key_map[key]
            for key, potential_key in source.key_map.items():
                # A potential is keyed by its SMIRKS. The prefix keeps a SMIRKS
                # the two force fields share from replacing the potential that
                # non-protein molecules use.
                copied = potential_key.model_copy(
                    update={"id": f"{self._protein_ff}::{potential_key.id}"}
                )
                target.key_map[key] = copied
                target.potentials[copied] = source.potentials[potential_key]
            # Parameterization checked a modified protein's coverage site by
            # site (GLYCAM06j impropers are written as proper torsions).
            if name != "ImproperTorsions" and self._modified_protein is None:
                before = {tuple(sorted(key.atom_indices)) for key in removed}
                after = {tuple(sorted(key.atom_indices)) for key in source.key_map}
                if before != after:
                    example = sorted(before - after or after - before)[0]
                    raise RuntimeError(
                        f"{self._protein_ff} does not cover the protein {name} terms "
                        f"({len(before - after)} missing, {len(after - before)} extra, "
                        f"for example atoms {example})"
                    )
            used = set(target.key_map.values())
            target.potentials = {
                key: value for key, value in target.potentials.items() if key in used
            }
            source_name = (
                f"{self._sm_ff} policy at {self._protein_ff} bond lengths"
                if name == "Constraints"
                else self._protein_ff
            )
            LOGGER.info(f"  {name}: {len(source.key_map)} protein terms from {source_name}")

        LOGGER.info(f"  Protein force field applied in {time.perf_counter() - t0:.1f}s")

    def _renumber_chains(self, topology: Topology) -> None:
        """Normalize chain IDs before canonical PDB identifier assignment.

        This lightweight normalization keeps intermediate OpenFF metadata
        deterministic before ``_assign_pdb_identifiers()`` applies the
        canonical PolyzyMD chain convention during topology export.

        Args:
            topology: Topology to modify.
        """
        for i, mol in enumerate(topology.molecules):
            for atom_index, atom in enumerate(mol.atoms):
                # Record the input residue, chain included, before the chain is replaced.
                self._get_original_residue_token(atom, atom_index)
                atom.metadata["chain_id"] = str(i + 1)

    def _assign_pdb_identifiers(self) -> None:
        """Assign PDB-compliant chain IDs and residue numbers to the solvated topology.

        This method assigns chain IDs and residue numbers based on molecule types
        tracked during the build process. The config YAML serves as the single
        source of truth for what each molecule represents.

        Chain assignment uses FIXED letters regardless of component presence:

        - Chain A: Protein/enzyme molecules (continuous residues across molecules)
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

        # Fixed chain letter assignments per component type
        # A=protein, B=substrate, C=polymer, D+=solvent — regardless of whether
        # a component is present. This ensures downstream code (SystemComponentInfo,
        # AtomGroupResolver, from_topology()) always sees the expected chain IDs.
        PROTEIN_CHAIN = "A"
        SUBSTRATE_CHAIN = "B"
        POLYMER_CHAIN = "C"
        SOLVENT_START_IDX = 3  # index of 'D' in CHAIN_LETTERS

        # 1. Protein molecules always use chain A with continuous unique residues
        if self._n_enzyme_molecules > 0:
            LOGGER.debug(
                f"Assigning chain {PROTEIN_CHAIN} to {self._n_enzyme_molecules} protein molecule(s)"
            )

            protein_residue_num = 1

            # A PDB can hold several protein molecules (a homodimer such as
            # https://www.rcsb.org/structure/1HXW), or several chains joined into one
            # molecule by disulfides (insulin). All of them go on chain A, numbered one
            # after the other: two 100-residue chains become residues 1-100 and
            # 101-200. Adjust your downstream analysis accordingly!!!

            # A modified protein's attachments (glycans, grafted polymers, labels)
            # are segments of their own (G001, P001, M001, ...): their residues are
            # numbered from 1 within the segment, so "segid G001 and resid 1" is a
            # glycan's first residue, and the protein's own residues keep the numbers
            # an unmodified build gives them.
            segments = self._protein_segments()
            segment_numbers: dict[str, int] = {}
            protein_atom = 0
            for _ in range(self._n_enzyme_molecules):
                mol = self._solvated_topology.molecule(mol_idx)
                residue_map: dict[tuple[str, str], int] = {}
                for atom_index, atom in enumerate(mol.atoms):
                    residue_token = self._get_original_residue_token(atom, atom_index)
                    segment = segments[protein_atom] if segments else "PROT"
                    protein_atom += 1
                    atom.metadata["chain_id"] = PROTEIN_CHAIN
                    key = (segment, residue_token)
                    if key not in residue_map:
                        if segment == "PROT":
                            residue_map[key] = protein_residue_num
                            protein_residue_num += 1
                        else:
                            segment_numbers[segment] = segment_numbers.get(segment, 0) + 1
                            residue_map[key] = segment_numbers[segment]
                    atom.metadata["residue_number"] = str(residue_map[key])
                mol_idx += 1

        # 2. Substrate always uses chain B and residue 1
        if self._n_substrate_molecules > 0:
            LOGGER.debug(f"Assigning chain {SUBSTRATE_CHAIN} to substrate")

            for _ in range(self._n_substrate_molecules):
                mol = self._solvated_topology.molecule(mol_idx)
                for atom in mol.atoms:
                    atom.metadata["chain_id"] = SUBSTRATE_CHAIN
                    atom.metadata["residue_number"] = "1"
                mol_idx += 1

        # 3. Polymers always use chain C and continue residue numbering across chains
        if self._n_polymer_chains > 0:
            LOGGER.debug(
                f"Assigning chain {POLYMER_CHAIN} to {self._n_polymer_chains} polymer chain(s)"
            )

            # Track residue number across all polymer chains (continue, don't restart)
            polymer_residue_num = 1

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

        # 4. Solvent always starts at chain D
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

    @staticmethod
    def _get_original_residue_token(atom: Any, atom_index: int) -> str:
        """Return a stable residue grouping token for a protein atom.

        The first call stores source residue metadata before PolyzyMD writes
        canonical output residue numbers. Later calls reuse that stored token so
        protein residue numbering is idempotent and does not cascade from
        already-renumbered output metadata.

        Parameters
        ----------
        atom : Any
            OpenFF-like atom with mutable ``metadata``.
        atom_index : int
            Atom index within the molecule, used only for deterministic fallback
            metadata when residue fields are absent.

        Returns
        -------
        str
            Stable token identifying the atom's original residue group.
        """
        metadata = atom.metadata
        residue_token = metadata.get(_ORIGINAL_RESIDUE_TOKEN_KEY)
        if residue_token is not None:
            return str(residue_token)

        chain_id = metadata.get("chain_id") or ""
        residue_number = metadata.get("residue_number")
        residue_name = metadata.get("residue_name") or metadata.get("residue_name_3") or ""
        insertion_code = metadata.get("insertion_code") or ""
        if residue_number is None:
            residue_token = f"missing-residue-metadata:{atom_index}"
        else:
            # Chains joined by a disulfide are one molecule; A7 and B7 are two residues.
            residue_token = f"{chain_id!s}|{residue_number!s}|{residue_name!s}|{insertion_code!s}"

        metadata[_ORIGINAL_RESIDUE_TOKEN_KEY] = residue_token
        return residue_token

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
                    # Cycling back would give solvent the chain IDs that mark the
                    # protein, substrate and polymer, which analyses select on.
                    raise ValueError(
                        f"The solvent needs more than {len(chain_letters) - start_chain_idx} "
                        f"chains of {max_residue} residues each ({n_solvent} solvent "
                        "molecules), and chain IDs would repeat those of the solute."
                    )

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
        if topology is self._solvated_topology:
            from polyzymd.core.components import atom_segids, write_pdb_segids

            components = self.components()
            attachments = {segid for c in components for segid in c.segments}
            write_pdb_segids(
                path, atom_segids(components, topology.n_atoms), separate=attachments
            )
        LOGGER.info(f"Saved topology to {path}")

    def build_from_config(
        self,
        config: "SimulationConfig",
        working_dir: Optional[Union[str, Path]] = None,
        polymer_seed: Optional[int] = None,
        publish_topology: bool = True,
    ) -> Interchange:
        """Build the complete system from a configuration.

        This is the main entry point for config-driven builds.

        Args:
            config: SimulationConfig with all settings.
            working_dir: Working directory for output files.
            polymer_seed: Random seed for polymer generation. This is used as a
                fallback if config.polymers.random_seed is not set. Callers pass
                the replicate index. The same value seeds Packmol for polymer
                packing and solvation so that replicates start from independent
                coordinates (Packmol is otherwise deterministic for identical
                inputs).

        Returns:
            OpenFF Interchange ready for simulation.

        Raises:
            ValueError: If the config names something the build cannot make
                correctly (see ``SimulationConfig.require_buildable``).
        """
        config.require_buildable()
        self._working_dir = Path(working_dir) if working_dir else None
        self._solute_handles = {
            solute.type: (solute.name, solute.segid) for solute in getattr(config, "solutes", [])
        }
        from polyzymd.utils.charging import nagl_model_provenance

        self._build_provenance = {"packmol_seed": polymer_seed}
        # Record a NAGL model only when a component is charged with one.
        nagl_model = getattr(config, "nagl_model", None)
        if nagl_model is not None:
            self._build_provenance.update(nagl_model_provenance(nagl_model))
        if config.polymers and config.polymers.enabled:
            self._build_provenance["polymer_seed"] = polymer_seed
        if polymer_seed is None:
            LOGGER.warning(
                "No replicate seed supplied: Packmol will use its fixed default seed and "
                "replicates will share identical starting coordinates."
            )
        else:
            LOGGER.info(f"Packmol seed (from replicate): {polymer_seed}")

        # 1. Build enzyme
        LOGGER.info(f"Building enzyme: {config.enzyme.name}")
        if getattr(config.enzyme, "is_modified", False):
            self.build_modified_enzyme(config.enzyme, getattr(config, "_config_dir", None))
        else:
            self.build_enzyme(config.enzyme.pdb_path, config.enzyme.custom_substructures_path)

        # 2. Build substrate (if configured)
        if config.substrate:
            LOGGER.info(f"Building substrate: {config.substrate.name}")
            self._substrate_builder.build_from_config(config.substrate, nagl_model=self._nagl_model)
            self._substrate_molecule = self._substrate_builder.molecule
            self._substrate_charge_method = config.substrate.charge_method.value
            self._substrate_conformer_index = config.substrate.conformer_index
            self._n_substrate_molecules = 1

        # 3. Combine enzyme + substrate
        self.combine_solutes()

        # 3b. Compute the periodic cell FIRST, from the protein + substrate
        # alone.  The cell must not depend on where Packmol happens to put the
        # polymer chains, otherwise replicates of one condition end up with
        # different box volumes and therefore different water and ion counts.
        # Polymers get their room by reserving polymers.packing.padding on top
        # of solvent.box.padding.
        polymers_enabled = bool(config.polymers and config.polymers.enabled)
        box_vectors = None
        if polymers_enabled:
            box_vectors = self._solvent_builder.compute_box_vectors_from_config(
                self._combined_topology,
                config.solvent,
                extra_padding_nm=config.polymers.packing.padding,
            )
            LOGGER.info(
                "Deterministic periodic cell from the protein + substrate diameter "
                f"+ 2 x ({config.polymers.packing.padding} + {config.solvent.box.padding}) nm: "
                f"{box_vectors}"
            )

        # 4. Build and pack polymers (if configured)
        if polymers_enabled:
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
            end_groups = None
            force_field = None

            # A monomer without a name is named by its label. The SMILES and caps
            # build dynamic chains; in cached mode they are checked against the
            # build metadata of cached chains that dynamic generation made.
            monomers = config.polymers.monomers
            if any(m.smiles for m in monomers) or config.polymers.end_groups is not None:
                monomer_smiles = {(m.name or m.label): m.smiles for m in monomers if m.smiles}
                monomer_names = {m.label: (m.name or m.label) for m in monomers}
                residue_names = {
                    (m.name or m.label): m.residue_name for m in monomers if m.residue_name
                }
            if config.polymers.end_groups is not None:
                end_groups = {
                    end: [
                        option.model_dump() for option in getattr(config.polymers.end_groups, end)
                    ]
                    for end in ("head", "tail")
                }
            if generation_mode == "dynamic":
                force_field = [str(name) for name in config.polymers.force_field]
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
                end_groups=end_groups,
                force_field=force_field,
                charger_type=config.polymers.charger.value,
                max_retries=config.polymers.max_retries,
                cache_directory=config.polymers.cache_directory,
            )
            self._build_provenance["polymer_sequence_seed"] = effective_seed
            builders = self._polymer_chain_builders()
            if builders:
                from polyzymd.builders.polymer_generator import builder_summary

                self._build_provenance["polymer_builder"] = "; ".join(
                    builder_summary(record) for record in builders
                )

            # Get packing config (uses defaults if not specified)
            packing = config.polymers.packing
            self.pack_polymers(
                box_vectors=box_vectors,
                padding=packing.padding,
                tolerance=packing.tolerance,
                movebadrandom=packing.movebadrandom,
                working_directory=self._working_dir,
                seed=polymer_seed,
                confine_to_sphere=packing.confine_to_sphere,
                nloop=packing.nloop,
            )

        # 5. Solvate.  With a precomputed cell the packed topology is already
        # framed in the brick, so the solvent builder must not re-centre it.
        LOGGER.info("Solvating system")
        self._solvent_builder.solvate_from_config(
            self._combined_topology,
            config.solvent,
            seed=polymer_seed,
            box_vectors=box_vectors,
            nagl_model=self._nagl_model,
        )
        self._solvated_topology = self._solvent_builder.solvated_topology
        self._build_provenance["solvent_packmol_seed"] = self._solvent_builder.packmol_seed
        self._record_box_provenance(polymers_enabled)

        # Save solvated PDB if working dir specified
        self._assign_pdb_identifiers()
        if self._working_dir and publish_topology:
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
        if self.glycam_sites:
            from polyzymd.modifying_proteins.glycam import apply_to_openmm

            # The protein is the first molecule: its indices are the system's.
            changed = apply_to_openmm(omm_system, self.glycam_sites, offset=0)
            LOGGER.info("  GLYCAM06j 1-4 pairs and masses written: %d pairs", changed)
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

    def _protein_segments(self) -> List[str]:
        """The segment of each protein atom (``PROT`` or an attachment's), in atom order;
        empty for a protein without modifications."""
        if self._modified_protein is None:
            return []
        return list(self._modified_protein.record.segments)

    @property
    def protein_molecule_ranges(self) -> List[List[int]]:
        """``[start, stop)`` of each molecule of a modified protein (the system's first)."""
        if self._modified_protein is None:
            return []
        return self._modified_protein.record.molecules or [
            [0, len(self._modified_protein.record.segments)]
        ]

    @property
    def glycam_sites(self) -> List[Any]:
        """The protein's GLYCAM06j sites, whose 1-4 pairs and masses an Interchange
        cannot hold; empty without GLYCAM glycans. The protein is the first molecule."""
        if self._modified_protein is None:
            return []
        return self._modified_protein.glycam_sites

    def _record_box_provenance(self, deterministic_box: bool) -> None:
        """Record the final cell in the build provenance.

        ``box_vectors_nm`` and ``brick_nm`` let a user confirm that replicates
        of one condition really do share a box (and therefore the same water
        and ion counts) by diffing their ``build_provenance.json`` files.
        """
        import numpy as np

        self._build_provenance["deterministic_box"] = bool(deterministic_box)
        box_vecs = self._solvent_builder.box_vectors
        if box_vecs is None:
            return
        arr = np.asarray(box_vecs.m_as("nanometer"), dtype=float)
        self._build_provenance["box_vectors_nm"] = [[float(v) for v in row] for row in arr]
        self._build_provenance["brick_nm"] = [float(v) for v in np.diagonal(arr)]

    @property
    def build_provenance(self) -> Dict[str, Any]:
        """Seeds and settings that determine the build's starting coordinates.

        Populated by :meth:`build_from_config`; recorded in ``build_provenance.json``
        by :func:`polyzymd.simulation.artifact_integrity.publish_build_bundle`.
        """
        return dict(self._build_provenance)

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

    def components(self) -> List["ComponentRecord"]:
        """Return the components of the solvated system, in atom order.

        Each record names the atoms the component owns, where their parameters
        and charges come from (as the build was asked to make them), and the
        inputs it was made from. The physics audit checks the built system
        against these records.

        Raises:
            RuntimeError: If the system has not been solvated.
        """
        from polyzymd.core.components import (
            DEFAULT_SEGIDS,
            ChargeSource,
            ComponentRecord,
            Domain,
            check_layout,
            file_record,
            molecule_ranges,
        )

        if self._solvated_topology is None:
            raise RuntimeError("No solvated topology. Call solvate() first.")

        def handle(type_: str) -> Tuple[str, str]:
            # The config names its solutes; a builder driven from Python uses the type.
            return self._solute_handles.get(type_, (type_, DEFAULT_SEGIDS[type_]))

        sizes = [molecule.n_atoms for molecule in self._solvated_topology.molecules]
        records: List[ComponentRecord] = []
        cursor = 0  # index of the next molecule

        def take(n_molecules: int):
            nonlocal cursor
            start = sum(sizes[:cursor])
            ranges = molecule_ranges(sizes[cursor : cursor + n_molecules], start)
            cursor += n_molecules
            return ranges

        def charges(method: Optional[str], source_file=None) -> ChargeSource:
            model = (self._nagl_model or _default_nagl_model()) if method == "nagl" else None
            return ChargeSource(method=method, model=model, source_file=source_file)

        if self._n_enzyme_molecules and self._modified_protein is not None:
            take(self._n_enzyme_molecules)
            name, segid = handle("protein")
            records.append(self._modified_protein.component(name=name, segid=segid, start=0))
        elif self._n_enzyme_molecules:
            atoms = take(self._n_enzyme_molecules)
            records.append(
                ComponentRecord(
                    name=handle("protein")[0],
                    type="protein",
                    segid=handle("protein")[1],
                    placement="fixed",
                    molecules=atoms,
                    domains=[
                        Domain(
                            "protein",
                            tuple(atoms),
                            (self._protein_ff,),
                            ChargeSource("library"),
                            reference={
                                "pdb": str(self._enzyme_builder.pdb_path),
                                "custom_substructures": _optional_str(
                                    getattr(self._enzyme_builder, "custom_substructures_path", None)
                                ),
                            },
                        )
                    ],
                    provenance={"pdb": file_record(self._enzyme_builder.pdb_path)},
                )
            )
        if self._n_substrate_molecules:
            atoms = take(self._n_substrate_molecules)
            method = self._substrate_charge_method or "nagl"
            records.append(
                ComponentRecord(
                    name=handle("small_molecule")[0],
                    type="small_molecule",
                    segid=handle("small_molecule")[1],
                    placement="fixed",
                    molecules=atoms,
                    domains=[
                        Domain(
                            "substrate",
                            tuple(atoms),
                            (self._sm_ff,),
                            charges(method),
                            reference={
                                "sdf": str(self._substrate_builder.sdf_path),
                                "conformer_index": self._substrate_conformer_index or 0,
                            },
                        )
                    ],
                    provenance={
                        "sdf": file_record(self._substrate_builder.sdf_path),
                        "conformer_index": self._substrate_conformer_index,
                    },
                )
            )
        if self._n_polymer_chains:
            molecules, domains, sources = [], [], {}
            from polyzymd.utils.packmol import packing_order

            sequences = list(self._polymer_builder.sequence_counts) if self._polymer_builder else []
            for index in packing_order(self._polymer_counts):
                count = self._polymer_counts[index]
                atoms = take(count)
                molecules += atoms
                sequence = sequences[index] if index < len(sequences) else f"type{index}"
                source = (
                    self._polymer_builder.molecule_source(sequence) if self._polymer_builder else {}
                )
                sdf = source.get("sdf_path")
                domains.append(
                    Domain(
                        f"polymer:{sequence}",
                        tuple(atoms),
                        (self._sm_ff,),
                        charges(source.get("charge_method", "preset"), str(sdf) if sdf else None),
                        reference={"sdf": str(sdf)} if sdf else None,
                    )
                )
                sources[sequence] = {"chains": count, "sdf": file_record(sdf)}
            provenance = {"sequences": sources}
            builders = self._polymer_chain_builders()
            if builders:
                provenance["builder"] = builders[0] if len(builders) == 1 else builders
            chemistry = self._polymer_builder.chain_chemistry() if self._polymer_builder else {}
            if chemistry:
                provenance.update(chemistry)
            records.append(
                ComponentRecord(
                    name=handle("polymer")[0],
                    type="polymer",
                    segid=handle("polymer")[1],
                    placement="packed",
                    molecules=molecules,
                    domains=domains,
                    provenance=provenance,
                )
            )
        counts = self._solvent_builder._solvation_counts
        composition = self._solvent_builder._composition
        if counts is not None:
            water_model = composition.water_model if composition else "tip3p"
            water = take(counts.water)
            na, cl = take(counts.na), take(counts.cl)
            solvent_molecules = list(water)
            domains = [
                Domain(
                    "water",
                    tuple(water),
                    (self._sm_ff,),
                    ChargeSource("water_model", water_model),
                    reference={"smiles": "O"},
                )
            ]
            for (name, count), cosolvent in zip(
                counts.co_solvents, composition.co_solvents if composition else []
            ):
                atoms = take(count)
                solvent_molecules += atoms
                source = _cosolvent_charges(cosolvent, self._nagl_model)
                domains.append(
                    Domain(
                        f"co_solvent:{name}",
                        tuple(atoms),
                        (self._sm_ff,),
                        source,
                        reference=(
                            {"sdf": source.source_file}
                            if source.method == "preset"
                            else {"smiles": _cosolvent_smiles(cosolvent)}
                        ),
                    )
                )
            records.append(
                ComponentRecord(
                    name="solvent",
                    type="solvent",
                    segid=DEFAULT_SEGIDS["solvent"],
                    placement="packed",
                    molecules=solvent_molecules,
                    domains=domains,
                    provenance={
                        "water_model": water_model,
                        "waters": counts.water,
                        "co_solvents": dict(counts.co_solvents),
                    },
                )
            )
            if na or cl:
                records.append(
                    ComponentRecord(
                        name="ions",
                        type="ions",
                        segid=DEFAULT_SEGIDS["ions"],
                        placement="packed",
                        molecules=na + cl,
                        domains=[
                            Domain(
                                ion,
                                tuple(atoms),
                                (self._sm_ff,),
                                ChargeSource("library"),
                                reference={"smiles": smiles},
                            )
                            for ion, atoms, smiles in (("Na+", na, "[Na+]"), ("Cl-", cl, "[Cl-]"))
                            if atoms
                        ],
                        provenance={"Na+": counts.na, "Cl-": counts.cl},
                    )
                )
        if cursor != len(sizes):
            raise RuntimeError(
                f"{len(sizes) - cursor} molecules of the solvated system belong to no component"
            )
        check_layout(records, sum(sizes))
        return records
