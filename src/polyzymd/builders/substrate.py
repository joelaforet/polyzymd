"""
Builder for substrate/ligand components.

This module handles loading docked conformers from SDF files and
assigning partial charges using various methods (NAGL, Espaloma, AM1BCC).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, List, Literal, Optional, Union

from openff.toolkit import Molecule

if TYPE_CHECKING:
    from polyzymd.config.schema import ChargeMethod, SubstrateConfig

LOGGER = logging.getLogger(__name__)

# Type alias for charge methods
ChargeMethodType = Literal["nagl", "espaloma", "am1bcc"]


class SubstrateBuilder:
    """Builder for loading and preparing substrate/ligand structures.

    This class handles:
    - Loading docked conformers from SDF files
    - Selecting specific conformers
    - Assigning partial charges using NAGL, Espaloma, or AM1BCC
    - Setting residue metadata

    Example:
        >>> builder = SubstrateBuilder()
        >>> molecule = builder.build(
        ...     sdf_path="docked.sdf",
        ...     conformer_index=0,
        ...     charge_method="nagl"
        ... )
        >>> print(f"Loaded substrate with {molecule.n_atoms} atoms")
    """

    def __init__(self) -> None:
        """Initialize the SubstrateBuilder."""
        self._molecule: Optional[Molecule] = None
        self._all_conformers: Optional[List[Molecule]] = None
        self._sdf_path: Optional[Path] = None
        self._conformer_index: Optional[int] = None

    @property
    def molecule(self) -> Optional[Molecule]:
        """Get the loaded and charged substrate molecule."""
        return self._molecule

    @property
    def all_conformers(self) -> Optional[List[Molecule]]:
        """Get all conformers loaded from the SDF file."""
        return self._all_conformers

    @property
    def sdf_path(self) -> Optional[Path]:
        """Get the path to the loaded SDF file."""
        return self._sdf_path

    @property
    def conformer_index(self) -> Optional[int]:
        """Get the index of the selected conformer."""
        return self._conformer_index

    def build(
        self,
        sdf_path: Union[str, Path],
        conformer_index: int = 0,
        charge_method: ChargeMethodType = "nagl",
        residue_name: str = "LIG",
    ) -> Molecule:
        """Load a substrate conformer and assign partial charges.

        Args:
            sdf_path: Path to the SDF file containing docked conformers.
            conformer_index: Index of the conformer to use (0-indexed).
            charge_method: Method for assigning partial charges.
            residue_name: 3-letter residue name for the substrate.

        Returns:
            OpenFF Molecule with assigned partial charges.

        Raises:
            FileNotFoundError: If the SDF file does not exist.
            IndexError: If the conformer index is out of range.
            ValueError: If the charge method is not supported.
        """
        sdf_path = Path(sdf_path)

        if not sdf_path.exists():
            raise FileNotFoundError(f"Substrate SDF file not found: {sdf_path}")

        LOGGER.info(f"Loading substrate conformers from {sdf_path}")

        # Load all conformers from the SDF file
        all_mols = Molecule.from_file(str(sdf_path))

        # Handle both single molecule and list of molecules
        if isinstance(all_mols, Molecule):
            all_mols = [all_mols]

        self._all_conformers = all_mols
        self._sdf_path = sdf_path

        LOGGER.info(f"Loaded {len(all_mols)} conformer(s) from SDF")

        # Select the requested conformer
        if conformer_index >= len(all_mols):
            raise IndexError(
                f"Conformer index {conformer_index} out of range. "
                f"SDF contains {len(all_mols)} conformer(s)."
            )

        selected_mol = all_mols[conformer_index]
        self._conformer_index = conformer_index

        LOGGER.info(f"Selected conformer {conformer_index}")

        # Assign partial charges
        charged_mol = self._assign_charges(selected_mol, charge_method)

        # Set residue metadata
        self._set_residue_metadata(charged_mol, residue_name)

        self._molecule = charged_mol

        LOGGER.info(
            f"Successfully prepared substrate: {charged_mol.n_atoms} atoms, "
            f"total charge = {charged_mol.total_charge}"
        )

        return charged_mol

    def build_from_config(self, config: "SubstrateConfig") -> Molecule:
        """Load substrate from a configuration object.

        Args:
            config: SubstrateConfig with SDF path and options.

        Returns:
            OpenFF Molecule with assigned partial charges.
        """
        LOGGER.info(f"Building substrate: {config.name}")
        return self.build(
            sdf_path=config.sdf_path,
            conformer_index=config.conformer_index,
            charge_method=config.charge_method.value,
            residue_name=config.residue_name,
        )

    def _assign_charges(self, molecule: Molecule, charge_method: ChargeMethodType) -> Molecule:
        """Assign partial charges to the molecule.

        Args:
            molecule: OpenFF Molecule to charge.
            charge_method: Method for assigning charges.

        Returns:
            Molecule with assigned partial charges.

        Raises:
            ValueError: If the charge method is not supported.
        """
        LOGGER.info(f"Assigning partial charges using {charge_method}")

        from polyzymd.utils import get_charger

        try:
            charger = get_charger(charge_method)
            return charger.charge_molecule(molecule)
        except ValueError as e:
            raise ValueError(
                f"Unsupported charge method: {charge_method}. "
                f"Supported methods: nagl, espaloma, am1bcc"
            ) from e

    def _set_residue_metadata(self, molecule: Molecule, residue_name: str) -> None:
        """Set residue metadata on all atoms of the molecule.

        Args:
            molecule: OpenFF Molecule to modify.
            residue_name: 3-letter residue name to assign.
        """
        for atom in molecule.atoms:
            atom.metadata["residue_name"] = residue_name

        LOGGER.debug(f"Set residue name to '{residue_name}' for all atoms")

    def get_n_conformers(self) -> int:
        """Get the number of conformers available in the loaded SDF.

        Returns:
            Number of conformers.

        Raises:
            RuntimeError: If no SDF has been loaded.
        """
        if self._all_conformers is None:
            raise RuntimeError("No SDF file loaded. Call build() first.")
        return len(self._all_conformers)

    def validate(self) -> bool:
        """Validate the loaded substrate molecule.

        Returns:
            True if validation passes.

        Raises:
            RuntimeError: If no molecule has been loaded.
            ValueError: If validation fails.
        """
        if self._molecule is None:
            raise RuntimeError("No substrate molecule loaded. Call build() first.")

        # Check that the molecule has atoms
        if self._molecule.n_atoms == 0:
            raise ValueError("Substrate molecule contains no atoms")

        # Check that partial charges have been assigned
        if self._molecule.partial_charges is None:
            raise ValueError("Substrate molecule has no partial charges assigned")

        LOGGER.info("Substrate molecule validation passed")
        return True
