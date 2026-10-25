"""
Builder for polymer components.

This module handles random co-polymer sequence generation, loading pre-built
polymer structures from SDF files, and building new chains with mBuild when
cached structures are not available.

Supports two generation modes:
- Cached: Load pre-built SDF files from disk
- Dynamic: Build chains with mBuild from star-sited monomer and cap SMILES

Made by PolyzyMD, by Joseph R. Laforet Jr.
"""

from __future__ import annotations

import logging
import random
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union

if TYPE_CHECKING:
    from openff.toolkit import Molecule

    from polyzymd.config.schema import MonomerSpec, PolymerSolute

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

    1. Cached mode: Load pre-built SDF files from disk
       - Requires sdf_directory with pre-built polymer files
       - Filenames: {type_prefix}_seq={sequence}_{length}-mer_charged.sdf

    2. Dynamic mode: Build chains with mBuild
       - Requires each monomer's star-sited repeat unit and the star-sited caps
       - Builds each chain, assigns charges, relaxes and checks it
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
        ...     monomer_smiles={"SBMA": "[*]CC([*])(C)C(=O)O...", "EGPMA": "..."},
        ...     monomer_names={"A": "SBMA", "B": "EGPMA"},
        ...     end_groups={
        ...         "head": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
        ...         "tail": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
        ...     },
        ...     force_field=["openff-2.3.0.offxml"],
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
        # New parameters for dynamic generation
        generation_mode: str = "cached",
        monomer_smiles: Optional[Dict[str, str]] = None,
        monomer_names: Optional[Dict[str, str]] = None,
        residue_names: Optional[Dict[str, str]] = None,
        end_groups: Optional[Dict[str, List[dict]]] = None,
        force_field: Optional[List[str]] = None,
        charger_type: str = "nagl",
        max_retries: int = 10,
        nagl_model: Optional[str] = None,
    ) -> None:
        """Initialize the PolymerBuilder.

        Args:
            characters: Monomer unit labels (e.g., ["A", "B"]).
            probabilities: Selection probability for each monomer.
            length: Number of monomers per polymer chain.
            type_prefix: Prefix for filenames (e.g., "SBMA-EGPMA").
            sdf_directory: Directory containing pre-built polymer SDFs (cached mode).
            cache_directory: Directory for caching generated polymers.
            generation_mode: "cached" for pre-built SDFs, "dynamic" for on-the-fly generation.
            monomer_smiles: Dictionary of monomer name -> star-sited repeat-unit SMILES
                (dynamic mode).
            monomer_names: Dictionary of label -> monomer name (dynamic mode).
            residue_names: Dictionary of monomer name -> 3-char PDB residue name.
            end_groups: {"head": [...], "tail": [...]}, the cap options of each end
                (mappings with label, smiles, probability and optionally name and
                residue_name), drawn per chain (dynamic mode).
            force_field: Force field files the chains are relaxed with (dynamic mode).
            charger_type: Charge method ("nagl", "espaloma", "am1bcc") for dynamic mode.
            max_retries: Builds tried before a chain with a pierced ring or
                overlapping atoms is refused.

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

        # Dynamic generation parameters
        self._generation_mode = generation_mode.lower()
        self._monomer_smiles = monomer_smiles or {}
        self._monomer_names = monomer_names or {}
        self._residue_names = residue_names or {}
        self._end_groups = end_groups or {}
        self._force_field = list(force_field or [])
        self._charger_type = charger_type.lower()
        self._max_retries = max_retries
        self._nagl_model = nagl_model

        # Validate dynamic mode requirements
        if self._generation_mode == "dynamic":
            if not self._monomer_smiles:
                raise ValueError("Dynamic generation mode requires monomer_smiles")
            if not self._monomer_names:
                raise ValueError("Dynamic generation mode requires monomer_names")
            if not self._end_groups:
                raise ValueError("Dynamic generation mode requires end_groups (head and tail)")
            if self._length < 3:
                raise ValueError(
                    "Dynamic generation mode requires length >= 3: the length counts the two "
                    "ends and at least one repeat unit"
                )
            if not self._force_field:
                raise ValueError("Dynamic generation mode requires the polymer force_field")

        # Lazy-initialized generator for dynamic mode
        self._polymer_generator = None

        # State
        self._loaded_molecules: Dict[str, Molecule] = {}
        # Build metadata of each cached chain that dynamic generation made.
        self._cached_sidecars: Dict[str, dict] = {}
        self._warned_unchecked_chemistry = False
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

    def molecule_source(self, sequence: str) -> Dict[str, object]:
        """Return where the molecule of ``sequence`` came from, for the build's provenance.

        Returns a dict with ``sdf_path`` (the file holding the molecule and its
        charges, or None) and ``charge_method``: ``preset`` when the charges were
        read from that file, otherwise the charger that computed them (dynamic
        generation), with ``nagl_model`` for NAGL.
        """
        candidates = []
        if self._sdf_directory:
            candidates.append(self._get_sdf_path(sequence, self._sdf_directory))
        candidates.append(self._get_sdf_path(sequence, self._cache_directory))
        sdf_path = next((path for path in candidates if path.exists()), None)
        # Mirrors _get_or_create_molecule: in dynamic mode a user SDF that does not
        # look like a cache entry is loaded as given, charges included.
        user_sdf = self._sdf_directory and sdf_path == candidates[0]
        if user_sdf and not self._is_dynamic_cache_like_sdf(sdf_path, sequence):
            return {"sdf_path": sdf_path, "charge_method": "preset", "nagl_model": None}
        if self._generation_mode == "dynamic":
            return {
                "sdf_path": sdf_path,
                "charge_method": self._charger_type,
                "nagl_model": self._nagl_model if self._charger_type == "nagl" else None,
            }
        return {"sdf_path": sdf_path, "charge_method": "preset", "nagl_model": None}

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

        # Generate random sequences; in dynamic mode each chain also draws the
        # cap of each end, and its sequence is head cap + repeat units + tail cap.
        raw_sequences = [self._draw_sequence() for _ in range(count)]

        # A chain built from star-sited SMILES has a direction: its first repeat
        # unit carries the head cap, so "AAB" and "BAA" are different molecules
        # and dynamic mode keeps the drawn sequence. Cached SDF libraries are
        # named by the canonical form of each sequence.
        if self._generation_mode == "dynamic":
            canonical_sequences = list(raw_sequences)
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

        return molecules_for_packing, list(self._sequence_counts.values())

    def build_from_config(
        self, config: "PolymerSolute", seed: Optional[int] = None
    ) -> Tuple[List[Molecule], List[int]]:
        """Build polymers from a configuration object.

        This method extracts all configuration values and updates the builder
        state accordingly, then calls build() to generate the polymers.

        Args:
            config: PolymerSolute with polymer settings.
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

            self._end_groups = {
                end: [option.model_dump() for option in getattr(config.end_groups, end)]
                for end in ("head", "tail")
            }
            self._force_field = [str(name) for name in config.force_field]

            # Reset the generator to force re-initialization with new config
            self._polymer_generator = None

            LOGGER.info(
                f"Dynamic mode configured with monomers: {list(self._monomer_smiles.keys())}"
            )

        LOGGER.info(f"Building polymers: {config.type_prefix} (mode: {self._generation_mode})")
        return self.build(config.count, seed=seed)

    def _draw_sequence(self) -> str:
        """Draw one chain's sequence, with its two cap labels when the config gives caps."""
        if not self._end_groups:
            return generate_random_sequence(self._length, self._characters, self._probabilities)
        # With caps the length counts the two ends: a 5-mer is a head cap,
        # three repeat units and a tail cap.
        head = self._draw_cap("head")
        middle = generate_random_sequence(self._length - 2, self._characters, self._probabilities)
        return head + middle + self._draw_cap("tail")

    def _draw_cap(self, end: str) -> str:
        """Draw the label of one end's cap from its options' probabilities."""
        options = self._end_groups[end]
        return random.choices(
            [option["label"] for option in options],
            weights=[option["probability"] for option in options],
            k=1,
        )[0]

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

        # User-provided dynamic SDFs are allowed only when they do not resemble cache artifacts
        if self._sdf_directory:
            sdf_path = self._get_sdf_path(sequence, self._sdf_directory)
            if sdf_path.exists():
                if self._generation_mode == "dynamic" and self._is_dynamic_cache_like_sdf(
                    sdf_path, sequence
                ):
                    if self._dynamic_sdf_metadata_matches(sdf_path, sequence):
                        return self._load_from_sdf(sdf_path)
                    LOGGER.info(
                        "Regenerating dynamic SDF with missing or stale metadata: %s", sdf_path
                    )
                else:
                    return self._load_from_sdf(sdf_path)

        # Dynamic cache entries need metadata validation before reuse
        if self._generation_mode == "dynamic":
            return self._generate_polymer(sequence)

        # Try to load from cache directory in cached mode
        cache_path = self._get_sdf_path(sequence, self._cache_directory)
        if cache_path.exists():
            return self._load_from_sdf(cache_path)

        searched_paths = []
        if self._sdf_directory:
            searched_paths.append(str(self._get_sdf_path(sequence, self._sdf_directory)))
        searched_paths.append(str(cache_path))

        raise FileNotFoundError(
            f"Could not find SDF file for sequence '{sequence}'. "
            f"Searched: {searched_paths}. "
            "Add the file, or set polymers.generation_mode: dynamic (with monomer smiles) "
            "to generate missing chains."
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
        filename = f"{self._type_prefix}_seq={sequence}_{len(sequence)}-mer_charged.sdf"
        path = directory / filename
        if self._generation_mode != "dynamic" and not path.exists():
            # A chain dynamic mode built is named in its drawn direction.
            reverse = sequence[::-1]
            reversed_path = directory / (
                f"{self._type_prefix}_seq={reverse}_{len(reverse)}-mer_charged.sdf"
            )
            if reversed_path.exists():
                LOGGER.info(
                    "Sequence %s: loading the chain of the reversed sequence, %s",
                    sequence,
                    reversed_path.name,
                )
                return reversed_path
        return path

    def _is_dynamic_cache_like_sdf(self, sdf_path: Path, sequence: str) -> bool:
        """Return whether an SDF path looks like a generated dynamic cache artifact.

        Parameters
        ----------
        sdf_path : Path
            Candidate SDF path.
        sequence : str
            Canonical polymer sequence.

        Returns
        -------
        bool
            True when the path should be routed through dynamic cache validation.
        """
        if self._same_directory(sdf_path.parent, self._cache_directory):
            return True

        if self._get_dynamic_cache_metadata_path(sdf_path).exists():
            return True

        return self._has_generated_dynamic_cache_name(sdf_path, sequence)

    def _same_directory(self, path_a: Path, path_b: Path) -> bool:
        """Return whether two paths identify the same directory.

        Parameters
        ----------
        path_a : Path
            First path.
        path_b : Path
            Second path.

        Returns
        -------
        bool
            True when both paths resolve to the same directory.
        """
        return path_a.resolve(strict=False) == path_b.resolve(strict=False)

    def _dynamic_sdf_metadata_matches(self, sdf_path: Path, sequence: str) -> bool:
        """Return whether a dynamic SDF sidecar matches current metadata.

        Parameters
        ----------
        sdf_path : Path
            Candidate generated-style SDF path.
        sequence : str
            Canonical polymer sequence.

        Returns
        -------
        bool
            True when the SDF may be loaded directly.
        """
        self._ensure_generators_initialized()
        return self._polymer_generator.cache_metadata_matches(sdf_path, sequence)

    def _get_dynamic_cache_metadata_path(self, sdf_path: Path) -> Path:
        """Get the dynamic cache metadata sidecar path for an SDF.

        Parameters
        ----------
        sdf_path : Path
            SDF path.

        Returns
        -------
        Path
            Expected metadata sidecar path.
        """
        return sdf_path.with_name(f"{sdf_path.name}.metadata.json")

    def _has_generated_dynamic_cache_name(self, sdf_path: Path, sequence: str) -> bool:
        """Return whether a filename follows the generated dynamic cache pattern.

        Parameters
        ----------
        sdf_path : Path
            Candidate SDF path.
        sequence : str
            Canonical polymer sequence.

        Returns
        -------
        bool
            True when the filename has PolyzyMD's generated dynamic cache shape.
        """
        suffix = f"_seq={sequence}_{len(sequence)}-mer_charged.sdf"
        return sdf_path.name.endswith(suffix) and len(sdf_path.name) > len(suffix)

    def _load_from_sdf(self, sdf_path: Path) -> Molecule:
        """Load a polymer molecule from an SDF file.

        In cached mode, a file that dynamic generation wrote has a metadata
        sidecar; the monomer SMILES and caps the config gives must match it.

        Args:
            sdf_path: Path to the SDF file.

        Returns:
            OpenFF Molecule loaded from the SDF.

        Raises:
            ValueError: If the config's monomer SMILES or caps differ from the
                ones the cached chain was built from.
        """
        from polyzymd.utils import (
            get_largest_offmol,
            topology_from_sdf,
        )

        LOGGER.info(f"Loading polymer from {sdf_path}")
        if self._generation_mode != "dynamic":
            self._check_cached_chemistry(Path(sdf_path))

        topology = topology_from_sdf(str(sdf_path))
        molecule = get_largest_offmol(topology)

        return molecule

    def _check_cached_chemistry(self, sdf_path: Path) -> None:
        """Compare the config's monomer SMILES and caps with a cached chain's build metadata.

        Cached mode takes the chemistry from the SDF. When the config also
        names the SMILES or caps, a chain built from other ones would be
        simulated while the config says otherwise, so a mismatch is refused.
        A file without metadata from this builder cannot be checked; a
        warning says so once.
        """
        given_smiles = {
            label: self._monomer_smiles.get(name)
            for label, name in self._monomer_names.items()
            if self._monomer_smiles.get(name)
        }
        given = bool(given_smiles or self._end_groups)
        metadata = None
        sidecar = sdf_path.with_name(f"{sdf_path.name}.metadata.json")
        if sidecar.exists():
            import json

            try:
                metadata = json.loads(sidecar.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                metadata = None
        if not isinstance(metadata, dict) or "monomers" not in metadata:
            if given and not self._warned_unchecked_chemistry:
                self._warned_unchecked_chemistry = True
                LOGGER.warning(
                    "Cached mode reads each chain's chemistry from its SDF; the monomer smiles "
                    "and end_groups in the config are not checked against %s, which has no "
                    "build metadata.",
                    sdf_path,
                )
            return
        self._cached_sidecars[sdf_path.name] = metadata

        def same(a: str, b: str) -> bool:
            from rdkit import Chem, RDLogger

            RDLogger.DisableLog("rdApp.*")
            try:
                return Chem.CanonSmiles(a) == Chem.CanonSmiles(b)
            except Exception:  # unreadable SMILES: compare as written
                return a == b
            finally:
                RDLogger.EnableLog("rdApp.*")

        problems = []
        for label, smiles in sorted(given_smiles.items()):
            built = (metadata.get("monomers") or {}).get(label, {}).get("smiles")
            if built and not same(built, smiles):
                problems.append(f"monomer {label!r}: config {smiles!r}, chain built from {built!r}")
        for end, options in sorted(self._end_groups.items()):
            built = (metadata.get("end_groups") or {}).get(end)
            if not isinstance(built, dict):
                continue
            given = {option["label"]: option["smiles"] for option in options}
            smiles = given.get(built.get("label"))
            if smiles is None:
                problems.append(
                    f"end_groups.{end}: the chain's cap {built.get('label')!r} "
                    f"({built.get('smiles')!r}) is not one of the config's caps"
                )
            elif built.get("smiles") and not same(built["smiles"], smiles):
                problems.append(
                    f"end_groups.{end} {built['label']!r}: config {smiles!r}, chain built from "
                    f"{built['smiles']!r}"
                )
        if problems:
            raise ValueError(
                f"the cached chain {sdf_path} was built from other chemistry than the config "
                "names:\n  "
                + "\n  ".join(problems)
                + "\nfix: correct the config's smiles and end_groups, or point sdf_directory "
                "to chains built from them."
            )

    def chain_builders(self) -> List[dict]:
        """Return the builder record of each chain source this build used.

        Chains built now in dynamic mode give the installed mBuild; cached
        chains that dynamic generation made give the builder named in their
        metadata. Other cached chains give nothing.
        """
        records = []
        if self._generation_mode == "dynamic":
            if any(
                self.molecule_source(sequence).get("charge_method") != "preset"
                for sequence in self._sequence_counts or {}
            ):
                from polyzymd.builders.polymer_generator import mbuild_provenance

                records.append(mbuild_provenance())
        for metadata in self._cached_sidecars.values():
            builder = metadata.get("builder")
            if isinstance(builder, dict) and builder not in records:
                records.append(builder)
        return records

    def chain_chemistry(self) -> Dict[str, object]:
        """Return the monomer SMILES and caps the chains were built from, when known."""
        if self._generation_mode == "dynamic":
            return {
                "monomers": {
                    label: {
                        "name": name,
                        "residue_name": self._residue_names.get(name),
                        "smiles": self._monomer_smiles.get(name),
                    }
                    for label, name in sorted(self._monomer_names.items())
                },
                "end_groups": {
                    end: [
                        {key: value for key, value in option.items() if value is not None}
                        for option in options
                    ]
                    for end, options in self._end_groups.items()
                },
            }
        if self._cached_sidecars:
            return {
                "chains": {
                    name: {
                        "sequence": metadata.get("sequence"),
                        "monomers": metadata.get("monomers"),
                        "end_groups": metadata.get("end_groups"),
                    }
                    for name, metadata in sorted(self._cached_sidecars.items())
                }
            }
        return {}

    def _generate_polymer(self, sequence: str) -> Molecule:
        """Build the chain of ``sequence`` with mBuild (dynamic mode).

        Args:
            sequence: Canonical polymer sequence.

        Returns:
            OpenFF Molecule for the polymer, with partial charges.

        Raises:
            RuntimeError: If not in dynamic mode and generation was attempted.
            PolymerGenerationError: If the chain cannot be built.
        """
        if self._generation_mode != "dynamic":
            raise RuntimeError(
                f"Polymer generation not available in cached mode for sequence '{sequence}'. "
                f"Either provide pre-built SDF files or switch to dynamic generation mode."
            )

        self._ensure_generators_initialized()
        LOGGER.info(f"Generating polymer for sequence: {sequence}")
        return self._polymer_generator.generate_polymer(sequence)

    def _ensure_generators_initialized(self) -> None:
        """Create the mBuild polymer generator when first needed."""
        if self._polymer_generator is not None:
            return

        from polyzymd.builders.polymer_generator import PolymerGenerator

        monomers = []
        for label, name in self._monomer_names.items():
            residue_name = self._residue_names.get(name) or name[:3].upper()
            monomers.append((label, name, residue_name, self._monomer_smiles[name]))
        self._polymer_generator = PolymerGenerator.from_smiles(
            monomers,
            head=self._end_groups["head"],
            tail=self._end_groups["tail"],
            cache_directory=self._cache_directory,
            force_field=self._force_field,
            max_retries=self._max_retries,
            charger_type=self._charger_type,
            type_prefix=self._type_prefix,
            nagl_model=self._nagl_model,
        )

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
