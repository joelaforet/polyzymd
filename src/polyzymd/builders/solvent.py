"""
Builder for solvent components.

This module handles solvation with water, ions, and optional co-solvents
using PACKMOL for molecular packing.

Solvent Parameterization
------------------------
This module uses pre-computed partial charges for solvent molecules to ensure:

1. **Consistency**: All copies of a solvent molecule have identical parameters
2. **Speed**: Charges are computed once, not N times for N molecules
3. **Reproducibility**: Same charges across different simulation runs

For built-in solvents (water models, DMSO, ethanol, etc.), charges are loaded
from pre-computed SDF files. For custom solvents, AM1BCC charges are computed
once and cached in ~/.polyzymd/solvent_cache/ for future use.

See Also
--------
polyzymd.data.solvent_molecules : Module providing pre-parameterized solvents
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from openff.toolkit import Molecule, Topology
    from openff.units import Quantity

    from polyzymd.config.schema import SolventConfig

LOGGER = logging.getLogger(__name__)
AVOGADRO_CONSTANT = 6.02214076e23

# Water model type
WaterModelType = Literal["tip3p", "spce", "tip4p", "tip4pew", "opc"]

# Box shape type
BoxShapeType = Literal["cube", "rhombic_dodecahedron", "truncated_octahedron"]


@dataclass
class CoSolvent:
    """Specification for a co-solvent component.

    Supports two specification methods:

    - ``mole_fraction``: Specify as molecule fraction (0-1)
    - ``concentration``: Specify as molarity (mol/L)

    Note: The molecule is NOT created in __post_init__ to avoid running
    charge calculations prematurely. Instead, molecules are loaded via
    get_solvent_molecule() in SolventBuilder.solvate() which uses cached
    charges to ensure all copies have identical parameters.

    Attributes:
        name: Identifier for the co-solvent.
        smiles: SMILES string for the molecule.
        mole_fraction: Mole fraction (0-1), mutually exclusive with concentration.
        concentration: Molar concentration (mol/L), mutually exclusive with mole_fraction.
        density: Optional density in g/mL for metadata.
        residue_name: 3-letter residue name.
        molecule: OpenFF Molecule (assigned later with cached charges).
    """

    name: str
    smiles: str
    mole_fraction: float | None = None
    concentration: float | None = None
    density: float | None = None  # g/mL
    residue_name: str = "COS"
    molecule: Molecule | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        """Validate co-solvent specification.

        Note: We intentionally do NOT create the molecule here. This is done
        later in SolventBuilder.solvate() using get_solvent_molecule() which
        provides proper charge caching.
        """
        # Set residue name if not already 3 chars
        if len(self.residue_name) > 3:
            self.residue_name = self.residue_name[:3].upper()

        # Validate that exactly one composition method is used
        if self.mole_fraction is None and self.concentration is None:
            raise ValueError(
                f"CoSolvent '{self.name}': Must specify either mole_fraction or concentration"
            )
        if self.mole_fraction is not None and self.concentration is not None:
            raise ValueError(
                f"CoSolvent '{self.name}': Cannot specify both mole_fraction and concentration"
            )
        if self.mole_fraction is not None and not 0.0 < self.mole_fraction < 1.0:
            raise ValueError(f"CoSolvent '{self.name}': mole_fraction must be between 0 and 1")
        if self.concentration is not None and self.concentration <= 0.0:
            raise ValueError(f"CoSolvent '{self.name}': concentration must be greater than 0")


@dataclass
class SolvationCounts:
    """Molecule counts from solvation for PDB chain/residue assignment.

    This dataclass stores the number of each molecule type added during
    solvation. It is used by SystemBuilder to assign PDB-compliant chain
    IDs and residue numbers.

    The molecule order in the solvated topology is:
    [solute molecules] + [water] + [Na+] + [Cl-] + [co-solvent 1] + [co-solvent 2] + ...

    Attributes:
        water: Number of water molecules added.
        na: Number of Na+ ions added.
        cl: Number of Cl- ions added.
        co_solvents: List of (name, count) tuples for each co-solvent type.
    """

    water: int
    na: int
    cl: int
    co_solvents: List[Tuple[str, int]] = field(default_factory=list)

    @property
    def total_solvent_molecules(self) -> int:
        """Total number of solvent molecules (water + ions + co-solvents)."""
        cosolvent_total = sum(count for _, count in self.co_solvents)
        return self.water + self.na + self.cl + cosolvent_total


@dataclass
class SolventComposition:
    """Complete specification of solvent composition.

    Attributes:
        water_model: Water model to use.
        co_solvents: List of co-solvent specifications.
        nacl_concentration: NaCl concentration in mol/L.
        kcl_concentration: KCl concentration in mol/L.
        mgcl2_concentration: MgCl2 concentration in mol/L.
        neutralize: Whether to neutralize system charge.
    """

    water_model: WaterModelType = "tip3p"
    co_solvents: List[CoSolvent] = field(default_factory=list)
    nacl_concentration: float = 0.0  # mol/L
    kcl_concentration: float = 0.0
    mgcl2_concentration: float = 0.0
    neutralize: bool = True

    def __post_init__(self) -> None:
        """Validate co-solvent mole fractions."""
        total_cosolvent = sum(
            cs.mole_fraction for cs in self.co_solvents if cs.mole_fraction is not None
        )
        if total_cosolvent >= 1.0:
            raise ValueError(f"Total co-solvent mole fraction must be < 1.0, got {total_cosolvent}")

    @property
    def water_mole_fraction(self) -> float:
        """Get the water mole fraction in the neutral solvent mixture.

        Note: Only counts co-solvents specified by mole_fraction.
        Co-solvents specified by concentration don't reduce water count.
        """
        total_cosolvent = sum(
            cs.mole_fraction for cs in self.co_solvents if cs.mole_fraction is not None
        )
        return 1.0 - total_cosolvent


class SolventBuilder:
    """Builder for solvating molecular systems.

    This class handles:
    - Water solvation with different water models
    - Ion addition for neutralization and ionic strength
    - Co-solvent addition with specified mole fractions or concentrations
    - Box shape and padding configuration

    Example:
        >>> builder = SolventBuilder()
        >>> composition = SolventComposition(
        ...     water_model="tip3p",
        ...     nacl_concentration=0.15,
        ...     co_solvents=[CoSolvent("dmso", "CS(=O)C", mole_fraction=0.1)]
        ... )
        >>> solvated = builder.solvate(
        ...     topology=solute_topology,
        ...     composition=composition,
        ...     padding=1.2,  # nm
        ... )
    """

    def __init__(self) -> None:
        """Initialize the SolventBuilder."""
        self._solvated_topology: Optional[Topology] = None
        self._box_vectors: Optional[NDArray] = None
        self._composition: Optional[SolventComposition] = None
        self._solvation_counts: Optional[SolvationCounts] = None

    @property
    def solvated_topology(self) -> Optional[Topology]:
        """Get the solvated topology."""
        return self._solvated_topology

    @property
    def box_vectors(self) -> Optional[NDArray]:
        """Get the box vectors of the solvated system."""
        return self._box_vectors

    @property
    def solvation_counts(self) -> Optional[SolvationCounts]:
        """Get the molecule counts from solvation.

        Returns:
            SolvationCounts with water, ion, and co-solvent molecule counts,
            or None if solvate() has not been called.
        """
        return self._solvation_counts

    def solvate(
        self,
        topology: Topology,
        composition: Optional[SolventComposition] = None,
        padding: float = 1.2,
        box_shape: BoxShapeType = "rhombic_dodecahedron",
        target_density: float = 1.0,
        tolerance: float = 2.0,
    ) -> Topology:
        """Solvate a topology with water, ions, and optional co-solvents.

        Args:
            topology: OpenFF Topology to solvate.
            composition: Solvent composition specification.
            padding: Distance from solute to box edge in nm.
            box_shape: Box geometry.
            target_density: Target density in g/mL.
            tolerance: Minimum molecular spacing for PACKMOL in Angstrom.

        Returns:
            Solvated OpenFF Topology.
        """
        from openff.toolkit import Molecule
        from openff.units import Quantity

        from polyzymd.data.solvent_molecules import get_solvent_molecule
        from polyzymd.utils import boxvectors
        from polyzymd.utils.packmol import solvate_with_packmol

        if composition is None:
            composition = SolventComposition()

        self._composition = composition

        LOGGER.info(
            f"Solvating system with {composition.water_model} water, "
            f"padding={padding} nm, shape={box_shape}"
        )

        # Get box shape matrix
        box_shape_matrix = self._get_box_shape_matrix(box_shape)

        # Calculate box vectors from solute bounding box + padding
        padding_qty = Quantity(padding, "nanometer")
        min_box_vecs = boxvectors.get_topology_bbox(topology)
        box_vecs = boxvectors.pad_box_vectors_uniform(min_box_vecs, padding_qty)

        # Apply box shape transformation
        box_vecs = box_shape_matrix @ box_vecs

        LOGGER.info(f"Computed box vectors: {box_vecs}")

        # Center topology in box
        self._center_topology_in_box(topology, box_vecs)

        # Calculate box volume and target masses
        box_vol = boxvectors.get_box_volume(box_vecs, units_as_openmm=False)
        target_density_qty = Quantity(target_density, "gram / milliliter")
        target_mass = box_vol * target_density_qty
        solute_mass = self._calculate_topology_mass(topology)
        solvent_mass = target_mass - solute_mass

        LOGGER.info(f"Target solvent mass: {solvent_mass}")

        # Get water molecule with charges
        water = self._get_water_molecule(composition.water_model)

        # Calculate numbers of each component
        water_mass = sum(atom.mass for atom in water.atoms)
        molarity_pure_water = Quantity(55.5, "mole / liter")

        # Calculate ions
        na = Molecule.from_smiles("[Na+]")
        cl = Molecule.from_smiles("[Cl-]")
        nacl_mass = sum(atom.mass for atom in na.atoms) + sum(atom.mass for atom in cl.atoms)

        nacl_conc = Quantity(composition.nacl_concentration, "mole / liter")
        nacl_mass_fraction = (nacl_conc * nacl_mass) / (molarity_pure_water * water_mass)
        nacl_mass_to_add = solvent_mass * nacl_mass_fraction
        nacl_to_add = (nacl_mass_to_add / nacl_mass).m_as("dimensionless").round()

        # Neutralize system
        solute_charge = sum(mol.total_charge for mol in topology.molecules)
        if composition.neutralize:
            na_to_add = int(np.ceil(nacl_to_add - solute_charge.m / 2.0))
            cl_to_add = int(np.floor(nacl_to_add + solute_charge.m / 2.0))
        else:
            na_to_add = int(nacl_to_add)
            cl_to_add = int(nacl_to_add)

        # Load co-solvents before count calculations so molar masses are known
        cosolvent_masses: list[tuple[str, float, Any]] = []
        for cosolvent in composition.co_solvents:
            if cosolvent.molecule is None:
                # Load molecule with pre-computed charges to ensure all copies
                # of the same co-solvent have identical parameters
                cosolvent.molecule = get_solvent_molecule(
                    name=cosolvent.name,
                    smiles=cosolvent.smiles,
                    residue_name=cosolvent.residue_name,
                )

            cosolvent_molar_mass = sum(atom.mass for atom in cosolvent.molecule.atoms)
            if cosolvent.mole_fraction is not None:
                cosolvent_masses.append(
                    (cosolvent.name, cosolvent.mole_fraction, cosolvent_molar_mass)
                )
            elif cosolvent.concentration is None:
                # Should not reach here due to validation in CoSolvent.__post_init__
                raise ValueError(
                    f"CoSolvent '{cosolvent.name}' has neither mole_fraction nor concentration"
                )

        neutral_solvent_mass = solvent_mass - nacl_mass_to_add
        if cosolvent_masses:
            water_to_add, mole_fraction_counts = self._calculate_mole_fraction_counts(
                neutral_solvent_mass=neutral_solvent_mass,
                water_mass=water_mass,
                cosolvent_mole_fractions=cosolvent_masses,
            )
        else:
            water_to_add = self._round_dimensionless_to_int(neutral_solvent_mass / water_mass)
            mole_fraction_counts = []

        LOGGER.info(f"Adding {int(water_to_add)} water, {na_to_add} Na+, {cl_to_add} Cl-")

        # Build molecule and count lists
        solvent_molecules = [water, na, cl]
        solvent_counts = [int(water_to_add), na_to_add, cl_to_add]

        # Track co-solvent counts for SolvationCounts
        cosolvent_counts_list: List[Tuple[str, int]] = []

        # Add co-solvents (using cached charges for consistency)
        mole_fraction_count_index = 0
        for cosolvent in composition.co_solvents:
            if cosolvent.mole_fraction is not None:
                _, n_cosolvent = mole_fraction_counts[mole_fraction_count_index]
                mole_fraction_count_index += 1

                LOGGER.info(
                    f"Adding {n_cosolvent} {cosolvent.name} molecules "
                    f"({cosolvent.mole_fraction * 100:.1f} mol%)"
                )

            elif cosolvent.concentration is not None:
                # =============================================================
                # CONCENTRATION METHOD
                # =============================================================
                # Formula: n = C × V × N_A
                # Where:
                #   C   = concentration (mol/L)
                #   V   = box volume (L)
                #   N_A = Avogadro's number
                #
                # Convert volume to liters explicitly before applying N_A
                # =============================================================
                box_volume_liters = self._box_volume_liters(box_vol)
                n_cosolvent = self._count_concentration_molecules(
                    concentration_molar=cosolvent.concentration,
                    box_volume_liters=box_volume_liters,
                )

                LOGGER.info(
                    f"Adding {n_cosolvent} {cosolvent.name} molecules ({cosolvent.concentration} M)"
                )

            else:
                # Should not reach here due to validation in CoSolvent.__post_init__
                raise ValueError(
                    f"CoSolvent '{cosolvent.name}' has neither mole_fraction nor concentration"
                )

            solvent_molecules.append(cosolvent.molecule)
            solvent_counts.append(n_cosolvent)
            cosolvent_counts_list.append((cosolvent.name, n_cosolvent))

        # Pack the box using PACKMOL (with CONECT-record overflow protection)
        solvated_top = solvate_with_packmol(
            molecules=solvent_molecules,
            number_of_copies=solvent_counts,
            solute=topology,
            box_vectors=box_vecs,
            tolerance_angstrom=tolerance,
        )

        # Set residue names
        self._set_solvent_residue_names(solvated_top, composition)

        self._solvated_topology = solvated_top
        self._box_vectors = box_vecs

        # Store solvation counts for PDB chain/residue assignment
        self._solvation_counts = SolvationCounts(
            water=int(water_to_add),
            na=na_to_add,
            cl=cl_to_add,
            co_solvents=cosolvent_counts_list,
        )

        LOGGER.info(
            f"Solvation complete: {solvated_top.n_molecules} molecules, "
            f"{solvated_top.n_atoms} atoms"
        )

        return solvated_top

    def solvate_from_config(
        self,
        topology: Topology,
        config: "SolventConfig",
    ) -> Topology:
        """Solvate using configuration object.

        Args:
            topology: OpenFF Topology to solvate.
            config: SolventConfig with solvent settings.

        Returns:
            Solvated OpenFF Topology.
        """
        # Build composition from config
        co_solvents = [
            CoSolvent(
                name=cs.name,
                smiles=cs.smiles,  # type: ignore[arg-type]  # Validated in schema
                mole_fraction=cs.mole_fraction,
                concentration=cs.concentration,
                density=cs.density,
                residue_name=cs.residue_name or cs.name[:3].upper(),
            )
            for cs in config.co_solvents
        ]

        composition = SolventComposition(
            water_model=config.primary.model.value,
            co_solvents=co_solvents,
            nacl_concentration=config.ions.nacl_concentration,
            kcl_concentration=config.ions.kcl_concentration,
            mgcl2_concentration=config.ions.mgcl2_concentration,
            neutralize=config.ions.neutralize,
        )

        return self.solvate(
            topology=topology,
            composition=composition,
            padding=config.box.padding,
            box_shape=config.box.shape.value,
            target_density=config.box.target_density,
            tolerance=config.box.tolerance,
        )

    @staticmethod
    def _round_dimensionless_to_int(value: Any) -> int:
        """Round a dimensionless value to an integer molecule count.

        Parameters
        ----------
        value : Any
            Dimensionless scalar or quantity-like object with ``m_as``.

        Returns
        -------
        int
            Rounded integer count.
        """
        if hasattr(value, "m_as"):
            value = value.m_as("dimensionless")
        return int(round(float(value)))

    @staticmethod
    def _box_volume_liters(box_volume: Any) -> float:
        """Convert a box volume to liters.

        Parameters
        ----------
        box_volume : Any
            Volume as an OpenFF quantity-like object or a numeric value already
            expressed in liters.

        Returns
        -------
        float
            Box volume in liters.
        """
        if hasattr(box_volume, "m_as"):
            return float(box_volume.m_as("liter"))
        return float(box_volume)

    @staticmethod
    def _count_concentration_molecules(concentration_molar: float, box_volume_liters: float) -> int:
        """Count additive co-solvent molecules from concentration.

        Parameters
        ----------
        concentration_molar : float
            Co-solvent concentration in mol/L.
        box_volume_liters : float
            Simulation box volume in liters.

        Returns
        -------
        int
            Number of co-solvent molecules to add.
        """
        return int(round(concentration_molar * box_volume_liters * AVOGADRO_CONSTANT))

    @classmethod
    def _calculate_mole_fraction_counts(
        cls,
        neutral_solvent_mass: Any,
        water_mass: Any,
        cosolvent_mole_fractions: Sequence[tuple[str, float, Any]],
    ) -> tuple[int, list[tuple[str, int]]]:
        """Count water and mole-fraction co-solvents from a mass budget.

        Parameters
        ----------
        neutral_solvent_mass : Any
            Mass allocated to water plus mole-fraction co-solvents.
        water_mass : Any
            Molecular mass of one water molecule.
        cosolvent_mole_fractions : Sequence[tuple[str, float, Any]]
            Tuples of co-solvent name, mole fraction, and molecular mass.

        Returns
        -------
        water_count : int
            Number of water molecules in the neutral solvent mixture.
        cosolvent_counts : list[tuple[str, int]]
            Number of molecules for each mole-fraction co-solvent.
        """
        total_mole_fraction = sum(mole_fraction for _, mole_fraction, _ in cosolvent_mole_fractions)
        if total_mole_fraction >= 1.0:
            raise ValueError(
                f"Total co-solvent mole fraction must be < 1.0, got {total_mole_fraction}"
            )

        water_mole_fraction = 1.0 - total_mole_fraction
        average_mass = water_mass * water_mole_fraction
        for _, mole_fraction, cosolvent_mass in cosolvent_mole_fractions:
            average_mass += cosolvent_mass * mole_fraction

        total_neutral_molecules = cls._round_dimensionless_to_int(
            neutral_solvent_mass / average_mass
        )
        cosolvent_counts = [
            (name, int(round(mole_fraction * total_neutral_molecules)))
            for name, mole_fraction, _ in cosolvent_mole_fractions
        ]
        water_count = total_neutral_molecules - sum(count for _, count in cosolvent_counts)
        return water_count, cosolvent_counts

    def _get_box_shape_matrix(self, shape: BoxShapeType) -> NDArray:
        """Get the transformation matrix for the box shape.

        Args:
            shape: Box shape identifier.

        Returns:
            3x3 transformation matrix.
        """
        import openff.packmol as packmol

        if shape == "cube":
            return packmol.UNIT_CUBE
        elif shape == "rhombic_dodecahedron":
            return packmol.RHOMBIC_DODECAHEDRON
        elif shape == "truncated_octahedron":
            return packmol.TRUNCATED_OCTAHEDRON
        else:
            LOGGER.warning(f"Unknown box shape '{shape}', using cube")
            return packmol.UNIT_CUBE

    def _get_water_molecule(self, model: WaterModelType) -> Molecule:
        """Get a water molecule with the appropriate charges.

        Uses pre-computed charges from the solvent cache for consistency.

        Args:
            model: Water model identifier.

        Returns:
            OpenFF Molecule for water with correct partial charges.
        """
        from polyzymd.data.solvent_molecules import get_solvent_molecule

        # Map water model name to canonical form
        # get_solvent_molecule handles: tip3p, spce, tip4pew, opc, etc.
        return get_solvent_molecule(model)

    def _center_topology_in_box(self, topology: Topology, box_vecs: NDArray) -> None:
        """Center the topology at the center of the box.

        Args:
            topology: Topology to center.
            box_vecs: Box vectors (3x3 array with units).
        """
        from openff.toolkit import unit

        # Calculate center of mass
        com = self._calculate_center_of_mass(topology)

        # Calculate box center
        box_center = np.sum(box_vecs, axis=0) / 2

        # Translation vector
        trans_vec = box_center - com

        # Apply translation
        old_positions = topology.get_positions()
        new_positions = old_positions + trans_vec
        topology.set_positions(new_positions)

        LOGGER.debug(f"Centered topology, translation: {trans_vec}")

    def _calculate_center_of_mass(self, topology: Topology) -> NDArray:
        """Calculate the center of mass of a topology.

        Args:
            topology: Topology to analyze.

        Returns:
            Center of mass as array with units.
        """
        from openff.toolkit import unit

        all_coords = []
        all_masses = []

        for molecule in topology.molecules:
            if molecule.conformers:
                coords = molecule.conformers[0].m_as(unit.angstrom)
                all_coords.append(coords)

                for atom in molecule.atoms:
                    all_masses.append(atom.mass.magnitude)

        if not all_coords:
            return np.array([0.0, 0.0, 0.0]) * unit.angstrom

        positions = np.concatenate(all_coords, axis=0)
        masses = np.array(all_masses)

        com = np.average(positions, axis=0, weights=masses) * unit.angstrom
        return com

    def _calculate_topology_mass(self, topology: Topology) -> Quantity:
        """Calculate the total mass of a topology.

        Args:
            topology: Topology to analyze.

        Returns:
            Total mass with units.
        """
        total_mass = sum(sum(atom.mass for atom in mol.atoms) for mol in topology.molecules)
        return total_mass

    def _set_solvent_residue_names(
        self, topology: Topology, composition: SolventComposition
    ) -> None:
        """Set residue names for solvent molecules.

        Args:
            topology: Solvated topology.
            composition: Solvent composition info.
        """
        # Build SMILES to residue name mapping
        smiles_to_residue = {
            "[H][O][H]": "HOH",
            "[Na+]": "Na+",
            "[Cl-]": "Cl-",
        }

        # Add co-solvents
        for cs in composition.co_solvents:
            if cs.molecule:
                smiles_to_residue[cs.molecule.to_smiles()] = cs.residue_name

        # Apply to topology
        for mol in topology.molecules:
            mol_smiles = mol.to_smiles()
            if mol_smiles in smiles_to_residue:
                residue_name = smiles_to_residue[mol_smiles]
                for atom in mol.atoms:
                    atom.metadata["residue_name"] = residue_name

    def validate(self) -> bool:
        """Validate the solvated topology.

        Returns:
            True if validation passes.

        Raises:
            RuntimeError: If no topology has been solvated.
            ValueError: If validation fails.
        """
        if self._solvated_topology is None:
            raise RuntimeError("No solvated topology. Call solvate() first.")

        if self._solvated_topology.n_atoms == 0:
            raise ValueError("Solvated topology contains no atoms")

        if self._box_vectors is None:
            raise ValueError("No box vectors defined")

        LOGGER.info("Solvated topology validation passed")
        return True
