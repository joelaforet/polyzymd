"""
Builder for polymer components.

This module handles random co-polymer sequence generation, loading pre-built
polymer structures from SDF files, and optionally generating new polymers
using Polymerist when cached structures are not available.
"""

from __future__ import annotations

import logging
import os
import random
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union

from openff.toolkit import Molecule

if TYPE_CHECKING:
    from polyzymd.config.schema import MonomerSpec, PolymerConfig

LOGGER = logging.getLogger(__name__)


def canonical_sequence(sequence: str) -> str:
    """Get the canonical form of a polymer sequence.

    Since polymers can be read in either direction, we use the
    lexicographically smaller of the sequence and its reverse
    as the canonical form.

    Args:
        sequence: Polymer sequence string (e.g., "AABBA").

    Returns:
        Canonical form of the sequence.

    Example:
        >>> canonical_sequence("ABAAA")
        'AAABA'  # reverse is smaller
        >>> canonical_sequence("AABBA")
        'AABBA'  # original is smaller
    """
    return min(sequence, sequence[::-1])


def generate_random_sequence(
    length: int,
    characters: List[str],
    probabilities: List[float],
) -> str:
    """Generate a random polymer sequence.

    Args:
        length: Number of monomers in the sequence.
        characters: List of monomer labels (e.g., ["A", "B"]).
        probabilities: Probability weights for each monomer.

    Returns:
        Random sequence string.

    Example:
        >>> generate_random_sequence(5, ["A", "B"], [0.7, 0.3])
        'AABAA'  # example output
    """
    return "".join(random.choices(characters, weights=probabilities, k=length))


class PolymerBuilder:
    """Builder for loading and generating polymer structures.

    This class implements a cache-first, generate-fallback approach:
    1. First checks for pre-built SDF files in the cache/SDF directory
    2. If not found, can optionally generate using Polymerist

    Example:
        >>> builder = PolymerBuilder(
        ...     characters=["A", "B"],
        ...     probabilities=[0.7, 0.3],
        ...     length=5,
        ...     sdf_directory="polymers/",
        ...     type_prefix="SBMA-EGPMA"
        ... )
        >>> molecules, sequences = builder.build(count=10)
        >>> print(f"Generated {len(molecules)} unique polymer types")
    """

    def __init__(
        self,
        characters: List[str],
        probabilities: List[float],
        length: int,
        type_prefix: str,
        sdf_directory: Optional[Union[str, Path]] = None,
        cache_directory: Optional[Union[str, Path]] = None,
        allow_generation: bool = False,
    ) -> None:
        """Initialize the PolymerBuilder.

        Args:
            characters: Monomer unit labels (e.g., ["A", "B"]).
            probabilities: Selection probability for each monomer.
            length: Number of monomers per polymer chain.
            type_prefix: Prefix for filenames (e.g., "SBMA-EGPMA").
            sdf_directory: Directory containing pre-built polymer SDFs.
            cache_directory: Directory for caching generated polymers.
            allow_generation: If True, generate missing polymers with Polymerist.

        Raises:
            ValueError: If probabilities don't sum to 1.0 or lengths mismatch.
        """
        if len(characters) != len(probabilities):
            raise ValueError("Characters and probabilities must have same length")

        if abs(sum(probabilities) - 1.0) > 1e-6:
            raise ValueError(f"Probabilities must sum to 1.0, got {sum(probabilities)}")

        self._characters = characters
        self._probabilities = probabilities
        self._length = length
        self._type_prefix = type_prefix
        self._sdf_directory = Path(sdf_directory) if sdf_directory else None
        self._cache_directory = Path(cache_directory) if cache_directory else Path(".polymer_cache")
        self._allow_generation = allow_generation

        # State
        self._loaded_molecules: Dict[str, Molecule] = {}
        self._sequence_counts: Optional[Counter] = None
        self._generated_sequences: Optional[List[str]] = None

    @property
    def characters(self) -> List[str]:
        """Get monomer characters."""
        return self._characters

    @property
    def probabilities(self) -> List[float]:
        """Get monomer probabilities."""
        return self._probabilities

    @property
    def length(self) -> int:
        """Get polymer length."""
        return self._length

    @property
    def loaded_molecules(self) -> Dict[str, Molecule]:
        """Get loaded polymer molecules keyed by canonical sequence."""
        return self._loaded_molecules

    @property
    def sequence_counts(self) -> Optional[Counter]:
        """Get counts of each unique sequence."""
        return self._sequence_counts

    def build(self, count: int, seed: Optional[int] = None) -> Tuple[List[Molecule], List[str]]:
        """Generate random polymer sequences and load/create corresponding molecules.

        Args:
            count: Number of polymer chains to generate.
            seed: Random seed for reproducibility.

        Returns:
            Tuple of (list of molecules for packing, list of canonical sequences).

        Raises:
            FileNotFoundError: If SDF file not found and generation not allowed.
        """
        if seed is not None:
            random.seed(seed)

        LOGGER.info(
            f"Generating {count} polymer chains with length {self._length}, "
            f"monomers: {dict(zip(self._characters, self._probabilities))}"
        )

        # Generate random sequences
        raw_sequences = [
            generate_random_sequence(self._length, self._characters, self._probabilities)
            for _ in range(count)
        ]

        # Convert to canonical form
        canonical_sequences = [canonical_sequence(seq) for seq in raw_sequences]
        self._generated_sequences = canonical_sequences

        # Count unique sequences
        self._sequence_counts = Counter(canonical_sequences)

        LOGGER.info(
            f"Generated {len(self._sequence_counts)} unique sequences: {dict(self._sequence_counts)}"
        )

        # Load molecules for each unique sequence
        molecules_for_packing = []
        sequences_for_packing = []

        for sequence, count in self._sequence_counts.items():
            mol = self._get_or_create_molecule(sequence)
            self._loaded_molecules[sequence] = mol

            # Add to packing lists (molecule repeated by count)
            molecules_for_packing.append(mol)
            sequences_for_packing.append(sequence)

        return molecules_for_packing, list(self._sequence_counts.values())

    def build_from_config(
        self, config: "PolymerConfig", seed: Optional[int] = None
    ) -> Tuple[List[Molecule], List[int]]:
        """Build polymers from a configuration object.

        Args:
            config: PolymerConfig with polymer settings.
            seed: Random seed for reproducibility.

        Returns:
            Tuple of (list of unique molecules, list of counts).
        """
        if not config.enabled:
            LOGGER.info("Polymers disabled in config, returning empty lists")
            return [], []

        # Extract characters and probabilities from config
        characters = [m.label for m in config.monomers]
        probabilities = [m.probability for m in config.monomers]

        # Update builder state
        self._characters = characters
        self._probabilities = probabilities
        self._length = config.length
        self._type_prefix = config.type_prefix

        if config.sdf_directory:
            self._sdf_directory = config.sdf_directory
        if config.cache_directory:
            self._cache_directory = config.cache_directory

        LOGGER.info(f"Building polymers: {config.type_prefix}")
        return self.build(config.count, seed=seed)

    def _get_or_create_molecule(self, sequence: str) -> Molecule:
        """Get a molecule for a sequence, loading from SDF or generating.

        Args:
            sequence: Canonical polymer sequence.

        Returns:
            OpenFF Molecule for the sequence.

        Raises:
            FileNotFoundError: If SDF not found and generation not allowed.
        """
        # Check if already loaded
        if sequence in self._loaded_molecules:
            return self._loaded_molecules[sequence]

        # Try to load from SDF directory
        if self._sdf_directory:
            sdf_path = self._get_sdf_path(sequence, self._sdf_directory)
            if sdf_path.exists():
                return self._load_from_sdf(sdf_path)

        # Try to load from cache directory
        cache_path = self._get_sdf_path(sequence, self._cache_directory)
        if cache_path.exists():
            return self._load_from_sdf(cache_path)

        # Generate if allowed
        if self._allow_generation:
            return self._generate_polymer(sequence)

        # No options left - raise error
        searched_paths = []
        if self._sdf_directory:
            searched_paths.append(str(self._get_sdf_path(sequence, self._sdf_directory)))
        searched_paths.append(str(cache_path))

        raise FileNotFoundError(
            f"Could not find SDF file for sequence '{sequence}'. "
            f"Searched: {searched_paths}. "
            f"Set allow_generation=True to generate missing polymers."
        )

    def _get_sdf_path(self, sequence: str, directory: Path) -> Path:
        """Get the expected SDF file path for a sequence.

        Args:
            sequence: Canonical polymer sequence.
            directory: Directory containing SDF files.

        Returns:
            Path to the expected SDF file.
        """
        # Format: {type_prefix}_{sequence}_{length}-mer_charged.sdf
        filename = f"{self._type_prefix}_{sequence}_{self._length}-mer_charged.sdf"
        return directory / filename

    def _load_from_sdf(self, sdf_path: Path) -> Molecule:
        """Load a polymer molecule from an SDF file.

        Args:
            sdf_path: Path to the SDF file.

        Returns:
            OpenFF Molecule loaded from the SDF.
        """
        from polyzymd.utils import (
            get_largest_offmol,
            topology_from_sdf,
        )

        LOGGER.info(f"Loading polymer from {sdf_path}")

        topology = topology_from_sdf(str(sdf_path))
        molecule = get_largest_offmol(topology)

        return molecule

    def _generate_polymer(self, sequence: str) -> Molecule:
        """Generate a polymer molecule using Polymerist.

        This is a placeholder - actual implementation would use
        Polymerist's polymer building capabilities.

        Args:
            sequence: Canonical polymer sequence.

        Returns:
            OpenFF Molecule for the polymer.

        Raises:
            NotImplementedError: Generation not yet implemented.
        """
        # TODO: Implement polymer generation using Polymerist
        # This would involve:
        # 1. Creating monomer Molecule objects
        # 2. Building the polymer chain
        # 3. Assigning partial charges
        # 4. Saving to cache directory
        raise NotImplementedError(
            f"Polymer generation not yet implemented for sequence '{sequence}'. "
            f"Please provide pre-built SDF files."
        )

    def get_packing_info(self) -> Tuple[List[Molecule], List[int]]:
        """Get molecules and counts for PACKMOL packing.

        Returns:
            Tuple of (list of unique molecules, list of counts).

        Raises:
            RuntimeError: If build() has not been called.
        """
        if self._sequence_counts is None:
            raise RuntimeError("No polymers generated. Call build() first.")

        molecules = []
        counts = []

        for sequence, count in self._sequence_counts.items():
            molecules.append(self._loaded_molecules[sequence])
            counts.append(count)

        return molecules, counts

    def validate(self) -> bool:
        """Validate the loaded polymers.

        Returns:
            True if validation passes.

        Raises:
            RuntimeError: If no polymers have been loaded.
            ValueError: If validation fails.
        """
        if not self._loaded_molecules:
            raise RuntimeError("No polymers loaded. Call build() first.")

        for sequence, mol in self._loaded_molecules.items():
            if mol.n_atoms == 0:
                raise ValueError(f"Polymer {sequence} has no atoms")

        LOGGER.info(f"Polymer validation passed for {len(self._loaded_molecules)} sequences")
        return True
