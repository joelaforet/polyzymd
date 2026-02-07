"""
Builder for polymer components.

This module handles random co-polymer sequence generation, loading pre-built
polymer structures from SDF files, and optionally generating new polymers
using Polymerist when cached structures are not available.

Supports two generation modes:
- Cached: Load pre-built SDF files from disk
- Dynamic: Generate polymers on-the-fly using Polymerist from raw monomer SMILES

Made by PolyzyMD, by Joseph R. Laforet Jr.
"""

from __future__ import annotations

import logging
import random
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union

from openff.toolkit import Molecule

if TYPE_CHECKING:
    from polyzymd.config.schema import MonomerSpec, PolymerConfig, ReactionConfig

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

    This class supports two generation modes:

    1. Cached mode (legacy): Load pre-built SDF files from disk
       - Requires sdf_directory with pre-built polymer files
       - Filenames: {type_prefix}_{sequence}_{length}-mer_charged.sdf

    2. Dynamic mode: Generate polymers on-the-fly using Polymerist
       - Requires monomer SMILES and ATRP reaction templates
       - Automatically generates fragments, builds chains, assigns charges
       - Caches results for subsequent runs

    Example (cached mode):
        >>> builder = PolymerBuilder(
        ...     characters=["A", "B"],
        ...     probabilities=[0.7, 0.3],
        ...     length=5,
        ...     sdf_directory="polymers/",
        ...     type_prefix="SBMA-EGPMA"
        ... )
        >>> molecules, counts = builder.build(count=10)

    Example (dynamic mode):
        >>> builder = PolymerBuilder(
        ...     characters=["A", "B"],
        ...     probabilities=[0.7, 0.3],
        ...     length=5,
        ...     type_prefix="SBMA-EGPMA",
        ...     generation_mode="dynamic",
        ...     monomer_smiles={"SBMA": "...", "EGPMA": "..."},
        ...     monomer_names={"A": "SBMA", "B": "EGPMA"},
        ...     reactions=reaction_config,
        ... )
        >>> molecules, counts = builder.build(count=10)
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
        # New parameters for dynamic generation
        generation_mode: str = "cached",
        monomer_smiles: Optional[Dict[str, str]] = None,
        monomer_names: Optional[Dict[str, str]] = None,
        residue_names: Optional[Dict[str, str]] = None,
        reactions: Optional["ReactionConfig"] = None,
        charger_type: str = "nagl",
        max_retries: int = 10,
    ) -> None:
        """Initialize the PolymerBuilder.

        Args:
            characters: Monomer unit labels (e.g., ["A", "B"]).
            probabilities: Selection probability for each monomer.
            length: Number of monomers per polymer chain.
            type_prefix: Prefix for filenames (e.g., "SBMA-EGPMA").
            sdf_directory: Directory containing pre-built polymer SDFs (cached mode).
            cache_directory: Directory for caching generated polymers.
            allow_generation: If True, generate missing polymers (for cached mode fallback).
            generation_mode: "cached" for pre-built SDFs, "dynamic" for on-the-fly generation.
            monomer_smiles: Dictionary of monomer name -> raw SMILES (dynamic mode).
            monomer_names: Dictionary of label -> monomer name (dynamic mode).
            residue_names: Dictionary of monomer name -> 3-char PDB residue name.
            reactions: ReactionConfig with paths to ATRP .rxn files (dynamic mode).
            charger_type: Charge method ("nagl", "espaloma", "am1bcc") for dynamic mode.
            max_retries: Maximum retries for polymer generation (ring-piercing failures).

        Raises:
            ValueError: If probabilities don't sum to 1.0 or lengths mismatch.
            ValueError: If dynamic mode but missing required parameters.
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

        # Dynamic generation parameters
        self._generation_mode = generation_mode.lower()
        self._monomer_smiles = monomer_smiles or {}
        self._monomer_names = monomer_names or {}
        self._residue_names = residue_names or {}
        self._reactions = reactions
        self._charger_type = charger_type.lower()
        self._max_retries = max_retries

        # Validate dynamic mode requirements
        if self._generation_mode == "dynamic":
            if not self._monomer_smiles:
                raise ValueError("Dynamic generation mode requires monomer_smiles")
            if not self._monomer_names:
                raise ValueError("Dynamic generation mode requires monomer_names")
            if not self._reactions:
                raise ValueError("Dynamic generation mode requires reactions (ReactionConfig)")

        # Lazy-initialized generators for dynamic mode
        self._fragment_generator = None
        self._polymer_generator = None
        self._monomer_group = None

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

        This method extracts all configuration values and updates the builder
        state accordingly, then calls build() to generate the polymers.

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

        # Update builder state - basic parameters
        self._characters = characters
        self._probabilities = probabilities
        self._length = config.length
        self._type_prefix = config.type_prefix

        if config.sdf_directory:
            self._sdf_directory = config.sdf_directory
        if config.cache_directory:
            self._cache_directory = config.cache_directory

        # Update dynamic generation parameters
        self._generation_mode = config.generation_mode.value
        self._charger_type = config.charger.value
        self._max_retries = config.max_retries

        # Extract monomer-related mappings for dynamic mode
        if config.generation_mode.value == "dynamic":
            # Build monomer_smiles: name -> SMILES
            self._monomer_smiles = {
                m.name: m.smiles for m in config.monomers if m.smiles is not None
            }

            # Build monomer_names: label -> name
            self._monomer_names = {m.label: m.name for m in config.monomers}

            # Build residue_names: name -> 3-char residue name
            self._residue_names = {
                m.name: m.residue_name for m in config.monomers if m.residue_name is not None
            }

            # Store reaction config
            self._reactions = config.reactions

            # Reset generators to force re-initialization with new config
            self._fragment_generator = None
            self._polymer_generator = None
            self._monomer_group = None

            LOGGER.info(
                f"Dynamic mode configured with monomers: {list(self._monomer_smiles.keys())}"
            )

        LOGGER.info(f"Building polymers: {config.type_prefix} (mode: {self._generation_mode})")
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

        # Dynamic mode: Generate polymers using FragmentGenerator and PolymerGenerator
        if self._generation_mode == "dynamic":
            return self._generate_polymer(sequence)

        # Cached mode: Try to load from SDF directory
        if self._sdf_directory:
            sdf_path = self._get_sdf_path(sequence, self._sdf_directory)
            if sdf_path.exists():
                return self._load_from_sdf(sdf_path)

        # Try to load from cache directory
        cache_path = self._get_sdf_path(sequence, self._cache_directory)
        if cache_path.exists():
            return self._load_from_sdf(cache_path)

        # Generate if allowed (cached mode fallback)
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

        For dynamic mode, this uses FragmentGenerator and PolymerGenerator
        to build polymer structures from raw monomer SMILES.

        Args:
            sequence: Canonical polymer sequence.

        Returns:
            OpenFF Molecule for the polymer.

        Raises:
            RuntimeError: If not in dynamic mode and generation was attempted.
            PolymerGenerationError: If generation fails after all retries.
        """
        if self._generation_mode != "dynamic":
            raise RuntimeError(
                f"Polymer generation not available in cached mode for sequence '{sequence}'. "
                f"Either provide pre-built SDF files or switch to dynamic generation mode."
            )

        # Ensure generators are initialized
        self._ensure_generators_initialized()

        # Generate the polymer using PolymerGenerator
        LOGGER.info(f"Generating polymer for sequence: {sequence}")
        mol = self._polymer_generator.generate_polymer(
            sequence=sequence,
            monomer_names=self._monomer_names,
            residue_names=self._residue_names if self._residue_names else None,
        )

        return mol

    def _ensure_generators_initialized(self) -> None:
        """Lazy-initialize FragmentGenerator and PolymerGenerator.

        This method ensures that the fragment and polymer generators are
        created only when first needed, avoiding overhead when loading
        pre-built polymers.
        """
        if self._fragment_generator is not None and self._polymer_generator is not None:
            return

        from polyzymd.builders.fragment_generator import FragmentGenerator
        from polyzymd.builders.polymer_generator import PolymerGenerator

        LOGGER.info("Initializing dynamic polymer generation pipeline...")

        # Create cache directory
        self._cache_directory.mkdir(parents=True, exist_ok=True)

        # Initialize fragment generator
        self._fragment_generator = FragmentGenerator(
            initiation_rxn_path=self._reactions.initiation,
            polymerization_rxn_path=self._reactions.polymerization,
            termination_rxn_path=self._reactions.termination,
            cache_directory=self._cache_directory,
        )

        # Load or generate MonomerGroup
        self._monomer_group = self._fragment_generator.load_or_generate(
            monomer_smiles=self._monomer_smiles,
            type_prefix=self._type_prefix,
        )

        LOGGER.info(
            f"MonomerGroup ready with {len(self._monomer_group.monomers)} fragments: "
            f"{list(self._monomer_group.monomers.keys())}"
        )

        # Initialize polymer generator
        self._polymer_generator = PolymerGenerator(
            monomer_group=self._monomer_group,
            cache_directory=self._cache_directory,
            max_retries=self._max_retries,
            charger_type=self._charger_type,
        )

        LOGGER.info("Dynamic polymer generation pipeline initialized")

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
