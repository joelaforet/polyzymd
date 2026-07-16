"""
Builder for polymer components.

This module handles random co-polymer sequence generation, loading legacy
sequence-derived polymer structures from SDF files, assembling native explicit
linear fragments, adding user-provided charged SDF molecules, and optionally
generating new polymers using the native bundled methacrylate generator.

Supports four free-polymer modes:

- Cached: Load legacy sequence-derived SDF files from disk
- Dynamic: Generate default methacrylate polymers natively from raw monomer SMILES
- Fragments: Assemble explicit terminal/middle fragments natively with mBuild ports
- Provided: Pack only explicit user-provided charged SDF molecule pools

In cached, dynamic, and fragments modes, ``provided_molecules`` remains additive
and merges opaque pre-generated charged SDF molecules with sequence-derived
chemistry before packing.

Made by PolyzyMD, by Joseph R. Laforet Jr.
"""

from __future__ import annotations

import logging
import random
import warnings
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

if TYPE_CHECKING:
    from openff.toolkit import Molecule

    from polyzymd.config.schema import MonomerSpec, PolymerConfig, ReactionConfig

LOGGER = logging.getLogger(__name__)
_SDF_DIRECTORY_DEPRECATION_WARNED = False

SDF_DIRECTORY_DEPRECATION_WARNING = (
    "polymers.sdf_directory is deprecated for pre-generated polymer inventories. Use "
    "polymers.provided_molecules for explicit charged SDF molecule pools. The legacy "
    "sdf_directory path derives filenames from generated polymer sequences and remains "
    "available for historical cached and dynamic workflows."
)


def _warn_sdf_directory_deprecated() -> None:
    """Emit the legacy sdf_directory deprecation warning once per process."""
    global _SDF_DIRECTORY_DEPRECATION_WARNED  # noqa: PLW0603
    if _SDF_DIRECTORY_DEPRECATION_WARNED:
        return
    warnings.warn(SDF_DIRECTORY_DEPRECATION_WARNING, UserWarning, stacklevel=3)
    _SDF_DIRECTORY_DEPRECATION_WARNED = True


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

    This class supports three sequence-derived generation modes plus one
    provided-only mode:

    - Cached mode loads legacy sequence-derived SDF files from disk using
      deprecated ``sdf_directory`` inventories with filenames such as
      ``{type_prefix}_seq={sequence}_{length}-mer_charged.sdf``.
    - Dynamic mode generates polymers on-the-fly from monomer SMILES using the
      native bundled methacrylate chemistry selected by default reaction
      markers. It builds chains, assigns charges, and caches results.
    - Fragments mode assembles explicit linear fragments natively from
      terminal/middle fragment specifications with mBuild ``Port`` and
      ``force_overlap`` stitching.
    - Provided mode skips sequence generation and packs only explicit charged
      SDF molecule pools supplied through ``provided_molecules``.

    In cached, dynamic, and fragments modes, ``provided_molecules`` supplies
    additional opaque charged SDF molecules that are merged with
    sequence-derived molecules before packing.

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
        characters: List[str] | None = None,
        probabilities: List[float] | None = None,
        length: int | None = None,
        type_prefix: str | None = None,
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
        fragments: dict[str, Any] | None = None,
        provided_molecules: list[Any] | None = None,
        polymer_random_seed: int | None = None,
    ) -> None:
        """Initialize the PolymerBuilder.

        Args:
            characters: Monomer unit labels (e.g., ["A", "B"]).
            probabilities: Selection probability for each monomer.
            length: Number of monomers per polymer chain.
            type_prefix: Prefix for filenames (e.g., "SBMA-EGPMA").
            sdf_directory: Deprecated directory containing sequence-derived SDFs
                for cached and historical dynamic modes.
            cache_directory: Directory for caching generated polymers.
            allow_generation: If True, generate missing polymers (for cached mode fallback).
            generation_mode: "cached", "dynamic", "fragments", or "provided".
            monomer_smiles: Dictionary of monomer name -> raw SMILES (dynamic mode).
            monomer_names: Dictionary of label -> monomer name (dynamic mode).
            residue_names: Dictionary of monomer name -> 3-char PDB residue name.
            reactions: ReactionConfig with literal "default" markers for native dynamic mode;
                custom .rxn files are unsupported.
            charger_type: Charge method ("nagl", "espaloma", "am1bcc") for dynamic mode.
            max_retries: Maximum retries for polymer generation (ring-piercing failures).
            fragments: Explicit terminal/middle fragment specs for fragments mode.
            provided_molecules: Additive opaque user-provided charged SDF molecule pools.
            polymer_random_seed: Polymer-level seed for provided molecule selection precedence.

        Raises:
            ValueError: If probabilities don't sum to 1.0 or lengths mismatch.
            ValueError: If dynamic or fragments mode is missing required parameters.
        """
        self._generation_mode = generation_mode.lower()
        characters = characters or []
        probabilities = probabilities or []
        if len(characters) != len(probabilities):
            raise ValueError("Characters and probabilities must have same length")

        if self._generation_mode != "provided" and abs(sum(probabilities) - 1.0) > 1e-6:
            raise ValueError(f"Probabilities must sum to 1.0, got {sum(probabilities)}")

        self._characters = characters
        self._probabilities = probabilities
        self._length = length
        self._type_prefix = type_prefix
        self._sdf_directory = Path(sdf_directory) if sdf_directory else None
        self._cache_directory = Path(cache_directory) if cache_directory else Path(".polymer_cache")
        self._allow_generation = allow_generation

        # Dynamic generation parameters
        self._monomer_smiles = monomer_smiles or {}
        self._monomer_names = monomer_names or {}
        self._residue_names = residue_names or {}
        self._reactions = reactions
        self._charger_type = charger_type.lower()
        self._max_retries = max_retries
        self._fragments = fragments or {}
        self._provided_molecules = provided_molecules or []
        self._polymer_random_seed = polymer_random_seed

        if self._generation_mode == "provided" and not self._provided_molecules:
            raise ValueError("Provided generation mode requires provided_molecules")
        # Validate dynamic mode requirements
        if self._generation_mode == "dynamic":
            if not self._monomer_smiles:
                raise ValueError("Dynamic generation mode requires monomer_smiles")
            if not self._monomer_names:
                raise ValueError("Dynamic generation mode requires monomer_names")
            if not self._reactions:
                raise ValueError("Dynamic generation mode requires reactions (ReactionConfig)")
        if self._generation_mode == "fragments" and not self._fragments:
            raise ValueError("Fragments generation mode requires fragment specifications")

        # State
        self._loaded_molecules: Dict[str, Molecule] = {}
        self._sequence_counts: Optional[Counter] = None
        self._generated_sequences: Optional[List[str]] = None
        self._packing_molecules: list[Molecule] | None = None
        self._packing_counts: list[int] | None = None

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

    def build(self, count: int, seed: Optional[int] = None) -> Tuple[List[Molecule], List[int]]:
        """Generate random polymer sequences and load/create corresponding molecules.

        Args:
            count: Number of polymer chains to generate.
            seed: Random seed for reproducibility.

        Returns:
            Tuple of (list of molecules for packing, list of counts).

        Raises:
            FileNotFoundError: If SDF file not found and generation not allowed.
        """
        rng = random.Random(seed)
        if self._generation_mode == "provided":
            molecules_for_packing, counts = self._build_provided_molecules(seed)
            self._generated_sequences = []
            self._sequence_counts = Counter()
            self._loaded_molecules = {}
            self._packing_molecules = molecules_for_packing
            self._packing_counts = counts
            return molecules_for_packing, counts

        if self._sdf_directory is not None:
            _warn_sdf_directory_deprecated()

        LOGGER.info(
            f"Generating {count} polymer chains with length {self._length}, "
            f"monomers: {dict(zip(self._characters, self._probabilities))}"
        )

        # Generate random sequences
        raw_sequences = [
            "".join(rng.choices(self._characters, weights=self._probabilities, k=self._length))
            for _ in range(count)
        ]

        # Fragments can be asymmetric, so exact sequence direction is preserved
        if self._generation_mode == "fragments":
            canonical_sequences = raw_sequences
        else:
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

        counts = list(self._sequence_counts.values())
        if self._provided_molecules:
            pool_molecules, pool_counts = self._build_provided_molecules(seed)
            molecules_for_packing.extend(pool_molecules)
            counts.extend(pool_counts)

        self._packing_molecules = molecules_for_packing
        self._packing_counts = counts
        return molecules_for_packing, counts

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
        self._fragments = config.fragments or {}
        self._provided_molecules = list(config.provided_molecules)
        self._polymer_random_seed = config.random_seed

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

            LOGGER.info(
                f"Dynamic mode configured with monomers: {list(self._monomer_smiles.keys())}"
            )

        LOGGER.info(f"Building polymers: {config.type_prefix} (mode: {self._generation_mode})")
        return self.build(config.count or 0, seed=seed)

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

        if self._generation_mode == "fragments":
            return self._generate_native_fragment_polymer(sequence)

        # Try to load from SDF directory (checked first in ALL modes)
        if self._sdf_directory:
            sdf_path = self._get_sdf_path(sequence, self._sdf_directory)
            if sdf_path.exists():
                return self._load_from_sdf(sdf_path)

        if self._generation_mode == "dynamic" and self._uses_native_methacrylate_backend():
            native_cache_path = self._native_artifact_paths(sequence).charged_sdf_path
            if native_cache_path.exists():
                return self._load_from_sdf(native_cache_path)
            return self._generate_polymer(sequence)

        # Try to load from cache directory for cached inventories
        cache_path = self._get_sdf_path(sequence, self._cache_directory)
        if cache_path.exists():
            return self._load_from_sdf(cache_path)

        if self._generation_mode == "dynamic":
            return self._generate_polymer(sequence)

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
        # Format: {type_prefix}_seq={sequence}_{length}-mer_charged.sdf
        filename = f"{self._type_prefix}_seq={sequence}_{self._length}-mer_charged.sdf"
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
        """Generate a polymer molecule using native dynamic routing.

        Args:
            sequence: Canonical polymer sequence.

        Returns:
            OpenFF Molecule for the polymer.

        Raises:
            RuntimeError: If not in dynamic mode and generation was attempted.
            PolymerGenerationError: If generation fails after all retries.
        """
        if self._generation_mode not in {"dynamic", "fragments"}:
            raise RuntimeError(
                f"Polymer generation not available in cached mode for sequence '{sequence}'. "
                f"Either provide pre-built SDF files or switch to dynamic generation mode."
            )

        if self._uses_native_methacrylate_backend():
            return self._generate_native_methacrylate_polymer(sequence)

        if self._generation_mode == "fragments":
            return self._generate_native_fragment_polymer(sequence)

        raise ValueError(
            "Dynamic polymer generation supports only native default methacrylate reactions. "
            "Use reactions: {initiation: default, polymerization: default, termination: default}; "
            "migrate custom chemistry to polymers.fragments, the CGSmiles notebook from "
            "`polyzymd init`, or polymers.provided_molecules."
        )

    def _generate_native_methacrylate_polymer(self, sequence: str) -> Molecule:
        """Generate a dynamic default polymer with the native methacrylate recipe.

        Args:
            sequence: Canonical polymer sequence.

        Returns:
            OpenFF Molecule for the charged polymer.
        """
        from polyzymd.builders.conjugation.polymer.native import (
            generate_native_methacrylate_polymer,
        )

        recipe = self._native_recipe(sequence)
        result = generate_native_methacrylate_polymer(
            recipe,
            self._cache_directory,
            sequence=sequence,
            charger_type=self._charger_type,
        )
        return result.charged_molecule

    def _generate_native_fragment_polymer(self, sequence: str) -> Molecule:
        """Generate a dynamic polymer with native explicit fragments.

        Parameters
        ----------
        sequence : str
            Exact sequence to assemble.

        Returns
        -------
        openff.toolkit.Molecule
            Charged OpenFF molecule for the assembled fragment polymer.
        """
        from polyzymd.builders.conjugation.polymer.fragments_native import (
            generate_native_fragment_polymer,
        )

        result = generate_native_fragment_polymer(
            self._type_prefix,
            self._fragments,
            sequence,
            self._cache_directory,
            charger_type=self._charger_type,
        )
        return result.charged_molecule

    def _native_recipe(self, sequence: str):
        """Return a native PolymerRecipe for the current dynamic builder state."""
        from polyzymd.builders.conjugation.polymer.recipe import (
            PolymerMonomerRecipe,
            PolymerRecipe,
        )

        monomers = []
        for label in self._characters:
            monomer_name = self._monomer_names[label]
            monomers.append(
                PolymerMonomerRecipe(
                    label=label,
                    name=monomer_name,
                    residue_name=self._residue_names.get(monomer_name, monomer_name[:3].upper()),
                    smiles=self._monomer_smiles[monomer_name],
                    probability=self._probabilities[self._characters.index(label)],
                )
            )
        return PolymerRecipe(
            name=self._type_prefix,
            monomers=tuple(monomers),
            length=self._length,
            fixed_sequence=sequence,
        )

    def _uses_native_methacrylate_backend(self) -> bool:
        """Return whether dynamic generation should use the bundled native backend."""
        if self._generation_mode != "dynamic" or self._reactions is None:
            return False
        from polyzymd.data.reactions import is_default_atrp_reaction_set

        return is_default_atrp_reaction_set(
            self._reactions.initiation,
            self._reactions.polymerization,
            self._reactions.termination,
        )

    def _native_artifact_paths(self, sequence: str):
        """Return centralized native artifact paths for a dynamic sequence."""
        from polyzymd.builders.conjugation.polymer.native import native_artifact_paths

        return native_artifact_paths(
            self._native_recipe(sequence),
            sequence,
            self._cache_directory,
            charger_type=self._charger_type,
        )

    def _native_fragment_artifact_paths(self, sequence: str):
        """Return centralized native fragment artifact paths for a sequence."""
        from polyzymd.builders.conjugation.polymer.fragments_native import (
            native_fragment_artifact_paths,
        )

        return native_fragment_artifact_paths(
            self._type_prefix,
            self._fragments,
            sequence,
            self._cache_directory,
            charger_type=self._charger_type,
        )

    def _build_provided_molecules(self, seed: int | None) -> tuple[list[Molecule], list[int]]:
        """Build additive provided molecules for the current config."""
        from polyzymd.builders.conjugation.polymer.provided_molecules import (
            build_provided_molecule_pool,
        )

        return build_provided_molecule_pool(
            self._provided_molecules,
            base_seed=self._polymer_random_seed,
            caller_seed=seed,
        )

    def get_packing_info(self) -> Tuple[List[Molecule], List[int]]:
        """Get molecules and counts for PACKMOL packing.

        Returns:
            Tuple of (list of unique molecules, list of counts).

        Raises:
            RuntimeError: If build() has not been called.
        """
        if self._packing_molecules is None or self._packing_counts is None:
            raise RuntimeError("No polymers generated. Call build() first.")
        return list(self._packing_molecules), list(self._packing_counts)

    def validate(self) -> bool:
        """Validate the loaded polymers.

        Returns:
            True if validation passes.

        Raises:
            RuntimeError: If no polymers have been loaded.
            ValueError: If validation fails.
        """
        if self._generation_mode == "provided":
            return self._validate_provided_packing()

        if not self._loaded_molecules:
            raise RuntimeError("No polymers loaded. Call build() first.")

        for sequence, mol in self._loaded_molecules.items():
            if mol.n_atoms == 0:
                raise ValueError(f"Polymer {sequence} has no atoms")

        LOGGER.info(f"Polymer validation passed for {len(self._loaded_molecules)} sequences")
        return True

    def _validate_provided_packing(self) -> bool:
        """Validate provided-only packing molecules and counts.

        Returns
        -------
        bool
            True if provided-only packing state is valid.

        Raises
        ------
        RuntimeError
            If provided-only build has not populated packing state.
        ValueError
            If molecule/count state is malformed or contains empty molecules.
        """
        if self._packing_molecules is None or self._packing_counts is None:
            raise RuntimeError("No provided molecules loaded. Call build() first.")
        if not self._packing_molecules:
            raise RuntimeError("No provided molecules loaded. Call build() first.")
        if len(self._packing_molecules) != len(self._packing_counts):
            raise ValueError("Provided molecule and count lists must have matching lengths")
        for index, (mol, count) in enumerate(zip(self._packing_molecules, self._packing_counts)):
            if count < 1:
                raise ValueError(f"Provided molecule count at index {index} must be positive")
            if getattr(mol, "n_atoms", 0) == 0:
                raise ValueError(f"Provided molecule at index {index} has no atoms")

        LOGGER.info(
            "Provided-only polymer validation passed for "
            f"{len(self._packing_molecules)} molecule entries"
        )
        return True
