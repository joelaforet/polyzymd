"""
Builder for enzyme/protein components.

This module handles loading PDB structures and partitioning them
for use with OpenFF force fields.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

from openff.toolkit import Topology

if TYPE_CHECKING:
    from polyzymd.config.schema import EnzymeConfig

LOGGER = logging.getLogger(__name__)


class EnzymeBuilder:
    """Builder for loading and preparing enzyme structures.

    This class handles:
    - Loading PDB structures into OpenFF Topology
    - Partitioning the protein for OpenFF compatibility
    - Basic validation of the loaded structure

    Example:
        >>> builder = EnzymeBuilder()
        >>> topology = builder.build("path/to/enzyme.pdb")
        >>> print(f"Loaded enzyme with {topology.n_atoms} atoms")
    """

    def __init__(self) -> None:
        """Initialize the EnzymeBuilder."""
        self._topology: Optional[Topology] = None
        self._pdb_path: Optional[Path] = None

    @property
    def topology(self) -> Optional[Topology]:
        """Get the loaded enzyme topology."""
        return self._topology

    @property
    def pdb_path(self) -> Optional[Path]:
        """Get the path to the loaded PDB file."""
        return self._pdb_path

    def build(self, pdb_path: Union[str, Path]) -> Topology:
        """Load and partition an enzyme structure from a PDB file.

        Args:
            pdb_path: Path to the enzyme PDB file.

        Returns:
            OpenFF Topology with the partitioned enzyme.

        Raises:
            FileNotFoundError: If the PDB file does not exist.
            ValueError: If partitioning fails.
        """
        pdb_path = Path(pdb_path)

        if not pdb_path.exists():
            raise FileNotFoundError(f"Enzyme PDB file not found: {pdb_path}")

        LOGGER.info(f"Loading enzyme from {pdb_path}")

        # Load the PDB into an OpenFF Topology
        topology = Topology.from_pdb(str(pdb_path))

        # Partition the protein for OpenFF compatibility
        # This assigns residue templates and ensures proper atom typing
        from polymerist.mdtools.openfftools.partition import partition

        was_partitioned = partition(topology)
        if not was_partitioned:
            raise ValueError(
                f"Failed to partition enzyme topology from {pdb_path}. "
                "Ensure the PDB contains a valid protein structure."
            )

        LOGGER.info(
            f"Successfully loaded enzyme: {topology.n_molecules} molecule(s), "
            f"{topology.n_atoms} atoms"
        )

        self._topology = topology
        self._pdb_path = pdb_path

        return topology

    def build_from_config(self, config: "EnzymeConfig") -> Topology:
        """Load enzyme from a configuration object.

        Args:
            config: EnzymeConfig with PDB path.

        Returns:
            OpenFF Topology with the partitioned enzyme.
        """
        LOGGER.info(f"Building enzyme: {config.name}")
        return self.build(config.pdb_path)

    def get_molecule(self) -> "Topology":
        """Get the first (and typically only) molecule from the topology.

        Returns:
            The enzyme molecule.

        Raises:
            RuntimeError: If no topology has been loaded.
        """
        if self._topology is None:
            raise RuntimeError("No enzyme topology loaded. Call build() first.")

        return self._topology.molecule(0)

    def validate(self) -> bool:
        """Validate the loaded enzyme topology.

        Returns:
            True if validation passes.

        Raises:
            RuntimeError: If no topology has been loaded.
            ValueError: If validation fails.
        """
        if self._topology is None:
            raise RuntimeError("No enzyme topology loaded. Call build() first.")

        # Check that we have at least one molecule
        if self._topology.n_molecules == 0:
            raise ValueError("Enzyme topology contains no molecules")

        # Check that the molecule has atoms
        if self._topology.n_atoms == 0:
            raise ValueError("Enzyme topology contains no atoms")

        LOGGER.info("Enzyme topology validation passed")
        return True
