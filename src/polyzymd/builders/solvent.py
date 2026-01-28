"""
Builder for solvent components.

This module handles solvation with water, ions, and optional co-solvents
using PACKMOL for molecular packing.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
from numpy.typing import NDArray
from openff.toolkit import Molecule, Topology
from openff.units import Quantity

if TYPE_CHECKING:
    from polyzymd.config.schema import SolventConfig

LOGGER = logging.getLogger(__name__)

# Water model type
WaterModelType = Literal["tip3p", "spce", "tip4p", "tip4pew", "opc"]

# Box shape type
BoxShapeType = Literal["cube", "rhombic_dodecahedron", "truncated_octahedron"]


@dataclass
class CoSolvent:
    """Specification for a co-solvent component.

    Supports two specification methods:
    - volume_fraction: Specify as fraction (0-1), e.g., 0.30 for 30% v/v
    - concentration: Specify as molarity (mol/L)

    Attributes:
        name: Identifier for the co-solvent.
        smiles: SMILES string for the molecule.
        volume_fraction: Volume fraction (0-1), mutually exclusive with concentration.
        concentration: Molar concentration (mol/L), mutually exclusive with volume_fraction.
        density: Density in g/mL (required for volume_fraction calculation).
        residue_name: 3-letter residue name.
        molecule: OpenFF Molecule (created from SMILES).
    """

    name: str
    smiles: str
    volume_fraction: Optional[float] = None
    concentration: Optional[float] = None
    density: Optional[float] = None  # g/mL
    residue_name: str = "COS"
    molecule: Optional[Molecule] = field(default=None, repr=False)

    def __post_init__(self) -> None:
        """Create the molecule from SMILES if not provided."""
        if self.molecule is None:
            self.molecule = Molecule.from_smiles(self.smiles)

        # Set residue name if not already 3 chars
        if len(self.residue_name) > 3:
            self.residue_name = self.residue_name[:3].upper()

        # Validate that exactly one specification method is used
        if self.volume_fraction is None and self.concentration is None:
            raise ValueError(
                f"CoSolvent '{self.name}': Must specify either volume_fraction or concentration"
            )
        if self.volume_fraction is not None and self.concentration is not None:
            raise ValueError(
                f"CoSolvent '{self.name}': Cannot specify both volume_fraction and concentration"
            )

        # Density is required for volume_fraction
        if self.volume_fraction is not None and self.density is None:
            raise ValueError(
                f"CoSolvent '{self.name}': density (g/mL) is required for volume_fraction calculation"
            )


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
    nacl_concentration: float = 0.1  # mol/L
    kcl_concentration: float = 0.0
    mgcl2_concentration: float = 0.0
    neutralize: bool = True

    @property
    def water_volume_fraction(self) -> float:
        """Get the water volume fraction (1 - sum of co-solvent fractions).

        Note: Only counts co-solvents specified by volume_fraction.
        Co-solvents specified by concentration don't reduce water volume.
        """
        total_cosolvent = sum(
            cs.volume_fraction for cs in self.co_solvents if cs.volume_fraction is not None
        )
        return 1.0 - total_cosolvent


class SolventBuilder:
    """Builder for solvating molecular systems.

    This class handles:
    - Water solvation with different water models
    - Ion addition for neutralization and ionic strength
    - Co-solvent addition with specified volume fractions
    - Box shape and padding configuration

    Example:
        >>> builder = SolventBuilder()
        >>> composition = SolventComposition(
        ...     water_model="tip3p",
        ...     nacl_concentration=0.15,
        ...     co_solvents=[CoSolvent("dmso", "CS(=O)C", 0.1)]
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

    @property
    def solvated_topology(self) -> Optional[Topology]:
        """Get the solvated topology."""
        return self._solvated_topology

    @property
    def box_vectors(self) -> Optional[NDArray]:
        """Get the box vectors of the solvated system."""
        return self._box_vectors

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
        from openff.interchange.components import _packmol as packmol
        from polymerist.mdtools.openfftools import boxvectors
        from polymerist.mdtools.openfftools.solvation import solvents

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
        box_vol = boxvectors.get_box_volume(box_vecs, units_as_openm=False)
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

        # Calculate water to add
        water_mass_to_add = solvent_mass - nacl_mass_to_add
        water_to_add = (water_mass_to_add / water_mass).m_as("dimensionless").round()

        # Adjust for co-solvents (reduce water proportionally)
        if composition.co_solvents:
            water_fraction = composition.water_volume_fraction
            water_to_add = int(water_to_add * water_fraction)

        # Neutralize system
        solute_charge = sum(mol.total_charge for mol in topology.molecules)
        if composition.neutralize:
            na_to_add = int(np.ceil(nacl_to_add - solute_charge.m / 2.0))
            cl_to_add = int(np.floor(nacl_to_add + solute_charge.m / 2.0))
        else:
            na_to_add = int(nacl_to_add)
            cl_to_add = int(nacl_to_add)

        LOGGER.info(f"Adding {int(water_to_add)} water, {na_to_add} Na+, {cl_to_add} Cl-")

        # Build molecule and count lists
        solvent_molecules = [water, na, cl]
        solvent_counts = [int(water_to_add), na_to_add, cl_to_add]

        # Add co-solvents
        for cosolvent in composition.co_solvents:
            if cosolvent.molecule is None:
                cosolvent.molecule = Molecule.from_smiles(cosolvent.smiles)

            # Get molar mass of co-solvent
            cosolvent_molar_mass = sum(atom.mass for atom in cosolvent.molecule.atoms)

            if cosolvent.volume_fraction is not None:
                # =============================================================
                # VOLUME FRACTION METHOD
                # =============================================================
                # Formula: n = (V_box × φ × ρ) / M
                # Where:
                #   V_box = simulation box volume
                #   φ     = volume fraction (e.g., 0.30 for 30%)
                #   ρ     = co-solvent density (g/mL)
                #   M     = molar mass (g/mol)
                #
                # This assumes ideal mixing (volumes are additive).
                # =============================================================
                cosolvent_density = Quantity(cosolvent.density, "gram / milliliter")
                # Volume of co-solvent = box_volume × volume_fraction
                cosolvent_volume = box_vol * cosolvent.volume_fraction
                # Mass of co-solvent = volume × density
                cosolvent_mass_to_add = cosolvent_volume * cosolvent_density
                # Number of molecules = mass / molar_mass
                n_cosolvent = int(
                    (cosolvent_mass_to_add / cosolvent_molar_mass).m_as("dimensionless").round()
                )

                LOGGER.info(
                    f"Adding {n_cosolvent} {cosolvent.name} molecules "
                    f"({cosolvent.volume_fraction * 100:.1f}% v/v, ρ={cosolvent.density} g/mL)"
                )

            elif cosolvent.concentration is not None:
                # =============================================================
                # CONCENTRATION METHOD
                # =============================================================
                # Formula: n = C × V × N_A
                # Where:
                #   C   = concentration (mol/L)
                #   V   = box volume (L)
                #   N_A = Avogadro's number (implicit in unit conversion)
                #
                # This directly gives the number of moles, then molecules.
                # =============================================================
                conc = Quantity(cosolvent.concentration, "mole / liter")
                n_moles = conc * box_vol
                n_cosolvent = int(n_moles.m_as("dimensionless").round())

                LOGGER.info(
                    f"Adding {n_cosolvent} {cosolvent.name} molecules ({cosolvent.concentration} M)"
                )

            else:
                # Should not reach here due to validation in CoSolvent.__post_init__
                raise ValueError(
                    f"CoSolvent '{cosolvent.name}' has neither volume_fraction nor concentration"
                )

            solvent_molecules.append(cosolvent.molecule)
            solvent_counts.append(n_cosolvent)

        # Pack the box using PACKMOL
        tolerance_qty = Quantity(tolerance, "angstrom")
        solvated_top = packmol.pack_box(
            solvent_molecules,
            solvent_counts,
            solute=topology,
            tolerance=tolerance_qty,
            box_vectors=box_vecs,
        )

        # Set residue names
        self._set_solvent_residue_names(solvated_top, composition)

        self._solvated_topology = solvated_top
        self._box_vectors = box_vecs

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
                volume_fraction=cs.volume_fraction,
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

    def _get_box_shape_matrix(self, shape: BoxShapeType) -> NDArray:
        """Get the transformation matrix for the box shape.

        Args:
            shape: Box shape identifier.

        Returns:
            3x3 transformation matrix.
        """
        from openff.interchange.components import _packmol as packmol

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

        Args:
            model: Water model identifier.

        Returns:
            OpenFF Molecule for water.
        """
        from polymerist.mdtools.openfftools.solvation import solvents

        if model == "tip3p":
            return solvents.water_TIP3P
        else:
            # For other models, just use basic water
            # The force field will apply the correct parameters
            return Molecule.from_smiles("O")

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
