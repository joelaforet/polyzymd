"""Polymer generator for dynamic polymer generation.

This module builds complete polymer chains from MonomerGroup fragments,
including 3D structure generation, ring-piercing validation, partial charge
assignment, and metadata-validated dynamic caching.

The workflow mirrors the notebook process:

- Set terminal orientations based on sequence head and tail
- Call ``build_linear_polymer()`` with the middle sequence
- Validate no ring-piercing and retry up to ``max_retries`` times
- Create an OpenFF topology and partition it
- Assign partial charges with NAGL, Espaloma, or AM1BCC
- Cache the charged SDF with a validation sidecar

"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from openff.toolkit import Molecule as OFFMolecule
    from openff.toolkit import Topology as OFFTopology
    from polymerist.polymers.monomers import MonomerGroup

logger = logging.getLogger(__name__)

DYNAMIC_CACHE_SCHEMA_VERSION = 3


def _build_linear_polymer(**kwargs: Any) -> Any:
    """Build a linear polymer through Polymerist.

    Parameters
    ----------
    **kwargs : Any
        Keyword arguments forwarded to Polymerist.

    Returns
    -------
    Any
        Polymerist chain object.
    """
    from polymerist.polymers.building import build_linear_polymer

    return build_linear_polymer(**kwargs)


def _make_monomer_group(monomers: dict[str, Any]) -> "MonomerGroup":
    """Create a Polymerist MonomerGroup lazily.

    Parameters
    ----------
    monomers : dict[str, Any]
        Mapping of fragment names to fragment chemistry objects.

    Returns
    -------
    MonomerGroup
        Polymerist monomer group containing the supplied fragments.
    """
    from polymerist.polymers.monomers import MonomerGroup

    return MonomerGroup(monomers=monomers)


def _mbmol_to_openmm_pdb(pdb_path: Path, chain: Any, resname_map: dict[str, str]) -> None:
    """Write a Polymerist chain to an OpenMM-compatible PDB lazily.

    Parameters
    ----------
    pdb_path : Path
        Destination PDB path.
    chain : Any
        Polymerist chain object.
    resname_map : dict[str, str]
        Fragment-to-residue-name mapping.
    """
    from polymerist.polymers.building import mbmol_to_openmm_pdb

    mbmol_to_openmm_pdb(pdb_path, chain, resname_map=resname_map)


def _mbmol_to_rdmol(chain: Any) -> Any:
    """Convert a Polymerist chain to an RDKit molecule lazily.

    Parameters
    ----------
    chain : Any
        Polymerist chain object.

    Returns
    -------
    Any
        RDKit molecule-like object.
    """
    from polymerist.polymers.building import mbmol_to_rdmol

    return mbmol_to_rdmol(chain)


def _summarize_ring_piercing(mol: Any) -> dict[str, Any]:
    """Summarize ring-piercing events through Polymerist lazily.

    Parameters
    ----------
    mol : Any
        RDKit molecule-like object.

    Returns
    -------
    dict[str, Any]
        Ring-piercing summary from Polymerist.
    """
    from polymerist.rdutils.rdcoords.piercing import summarize_ring_piercing

    return summarize_ring_piercing(mol)


def _topology_from_sdf(sdf_path: Path) -> Any:
    """Load an OpenFF topology from SDF through Polymerist lazily.

    Parameters
    ----------
    sdf_path : Path
        Source SDF path.

    Returns
    -------
    Any
        OpenFF topology-like object.
    """
    from polymerist.mdtools.openfftools.topology import topology_from_sdf

    return topology_from_sdf(sdf_path)


def _topology_to_sdf(sdf_path: Path, topology: Any) -> None:
    """Write an OpenFF topology to SDF through Polymerist lazily.

    Parameters
    ----------
    sdf_path : Path
        Destination SDF path.
    topology : Any
        OpenFF topology-like object.
    """
    from polymerist.mdtools.openfftools.topology import topology_to_sdf

    topology_to_sdf(sdf_path, topology)


def _get_largest_offmol(topology: Any) -> "OFFMolecule":
    """Extract the largest OpenFF molecule through Polymerist lazily.

    Parameters
    ----------
    topology : Any
        OpenFF topology-like object.

    Returns
    -------
    OFFMolecule
        Largest molecule in the topology.
    """
    from polymerist.mdtools.openfftools.topology import get_largest_offmol

    return get_largest_offmol(topology)


def _partition(topology: "OFFTopology") -> bool:
    """Partition an OpenFF topology through Polymerist lazily.

    Parameters
    ----------
    topology : OFFTopology
        OpenFF topology to partition.

    Returns
    -------
    bool
        True when partitioning succeeds.
    """
    from polymerist.mdtools.openfftools.partition import partition

    return partition(topology)


class PolymerGenerationError(Exception):
    """Raised when polymer generation fails after all retries."""


def _validate_dynamic_sequence_labels(sequence: str, monomer_names: dict[str, str]) -> list[str]:
    """Validate labels used by a dynamic polymer sequence.

    Parameters
    ----------
    sequence : str
        Polymer sequence string to validate.
    monomer_names : dict[str, str]
        Mapping of sequence labels to monomer names.

    Returns
    -------
    list[str]
        Unique sequence labels in first-seen order.

    Raises
    ------
    ValueError
        If any sequence label does not have a configured monomer name.
    """
    labels_used = []
    for label in sequence:
        if label not in monomer_names:
            raise ValueError(f"No monomer name configured for sequence label: {label}")
        if label not in labels_used:
            labels_used.append(label)

    return labels_used


class PolymerGenerator:
    """Generates complete polymer chains from MonomerGroup fragments.

    This class handles the full polymer building workflow:

    - Building 3D structures using Polymerist's ``build_linear_polymer``
    - Validating no ring-piercing events occur
    - Creating OpenFF topologies
    - Assigning partial charges
    - Caching results as SDF files guarded by deterministic metadata
    """

    def __init__(
        self,
        monomer_group: "MonomerGroup",
        cache_directory: Path,
        max_retries: int = 10,
        charger_type: str = "nagl",
    ):
        """Initialize the polymer generator.

        Parameters
        ----------
        monomer_group : MonomerGroup
            MonomerGroup containing all named fragments.
        cache_directory : Path
            Directory for caching polymer SDF files.
        max_retries : int, optional
            Maximum attempts for building after ring-piercing failures, by default 10.
        charger_type : str, optional
            Charge method, by default "nagl".
        """
        self.monomer_group = monomer_group
        self.cache_directory = Path(cache_directory)
        self.max_retries = max_retries
        self.charger_type = charger_type.lower()

        # Ensure cache directory exists
        self.cache_directory.mkdir(parents=True, exist_ok=True)

        # Lazy-loaded charger
        self._charger = None

    @property
    def charger(self):
        """Get the molecular charger lazily.

        Returns
        -------
        Any
            Molecular charger for the configured charge method.
        """
        if self._charger is None:
            self._charger = self._create_charger()
        return self._charger

    def _create_charger(self):
        """Create the appropriate charger based on ``charger_type``.

        Returns
        -------
        Any
            Polymerist molecular charger.

        Raises
        ------
        ValueError
            If ``charger_type`` is unknown.
        """
        if self.charger_type == "nagl":
            from polymerist.mdtools.openfftools.partialcharge.molchargers import (
                NAGLCharger,
            )

            return NAGLCharger()
        elif self.charger_type == "espaloma":
            from polymerist.mdtools.openfftools.partialcharge.molchargers import (
                EspalomaCharger,
            )

            return EspalomaCharger()
        elif self.charger_type == "am1bcc":
            from polymerist.mdtools.openfftools.partialcharge.molchargers import (
                AM1BCCCharger,
            )

            return AM1BCCCharger()
        else:
            raise ValueError(f"Unknown charger type: {self.charger_type}")

    def _get_terminal_fragment_name(self, monomer_name: str) -> str:
        """Get the 1-site fragment name for a monomer (terminal position).

        Parameters
        ----------
        monomer_name : str
            Base monomer name, such as "SBMA".

        Returns
        -------
        str
            Fragment name for 1-site, such as "SBMA_1-site".

        Raises
        ------
        ValueError
            If the terminal fragment is unavailable.
        """
        frag_name = f"{monomer_name}_1-site"
        if frag_name in self.monomer_group.monomers:
            return frag_name

        raise ValueError(f"No 1-site terminal fragment found for monomer: {monomer_name}")

    def _get_middle_fragment_name(self, monomer_name: str) -> str:
        """Get the 2-site fragment name for a monomer (middle position).

        Parameters
        ----------
        monomer_name : str
            Base monomer name, such as "SBMA".

        Returns
        -------
        str
            Fragment name for 2-site, such as "SBMA_2-site".

        Raises
        ------
        ValueError
            If the middle fragment is unavailable.
        """
        frag_name = f"{monomer_name}_2-site"
        if frag_name in self.monomer_group.monomers:
            return frag_name

        raise ValueError(f"No 2-site middle fragment found for monomer: {monomer_name}")

    def _build_middle_sequence_map(
        self,
        middle_sequence: str,
        monomer_names: dict[str, str],
    ) -> dict[str, str]:
        """Build an explicit Polymerist sequence map for middle monomers.

        Parameters
        ----------
        middle_sequence : str
            Sequence labels excluding terminal monomers.
        monomer_names : dict[str, str]
            Mapping of sequence labels to monomer names.

        Returns
        -------
        dict[str, str]
            Mapping of middle sequence labels to exact 2-site fragment names.

        Raises
        ------
        ValueError
            If a sequence label or fragment is unavailable.
        """
        labels_used = _validate_dynamic_sequence_labels(middle_sequence, monomer_names)
        sequence_map = {}
        for label in labels_used:
            sequence_map[label] = self._get_middle_fragment_name(monomer_names[label])

        return sequence_map

    def _build_dynamic_cache_metadata(
        self,
        sequence: str,
        monomer_names: dict[str, str],
        residue_names: dict[str, str] | None = None,
    ) -> dict[str, Any]:
        """Build dynamic-cache metadata for a generated charged polymer.

        Parameters
        ----------
        sequence : str
            Polymer sequence string.
        monomer_names : dict[str, str]
            Mapping of sequence labels to monomer names.
        residue_names : dict[str, str] | None, optional
            Optional mapping of monomer names to residue names.

        Returns
        -------
        dict[str, Any]
            Metadata payload used to validate dynamic cached SDF files.
        """
        self._validate_dynamic_sequence_length(sequence)
        head_label = sequence[0]
        tail_label = sequence[-1]
        middle_sequence = sequence[1:-1] if len(sequence) > 2 else ""
        labels_used = _validate_dynamic_sequence_labels(sequence, monomer_names)

        head_monomer = monomer_names[head_label]
        tail_monomer = monomer_names[tail_label]
        sequence_map = self._build_middle_sequence_map(middle_sequence, monomer_names)
        terminal_fragments = {
            "head": self._get_terminal_fragment_name(head_monomer),
            "tail": self._get_terminal_fragment_name(tail_monomer),
        }
        used_fragment_names = sorted(set(terminal_fragments.values()) | set(sequence_map.values()))
        monomer_names_used = {label: monomer_names[label] for label in labels_used}
        residue_names_used = self._build_used_residue_names(monomer_names_used, residue_names)

        return {
            "schema_version": DYNAMIC_CACHE_SCHEMA_VERSION,
            "charger_type": self.charger_type,
            "sequence": sequence,
            "middle_sequence": middle_sequence,
            "monomer_names_used": monomer_names_used,
            "residue_names_used": residue_names_used,
            "terminal_fragments": terminal_fragments,
            "sequence_map": sequence_map,
            "used_fragment_names": used_fragment_names,
            "monomer_group_fingerprint": self._build_monomer_group_fingerprint(),
        }

    def _validate_dynamic_sequence_length(self, sequence: str) -> None:
        """Validate that a sequence can be generated through Polymerist.

        Parameters
        ----------
        sequence : str
            Polymer sequence string to validate.

        Raises
        ------
        ValueError
            If ``sequence`` is too short for dynamic generation.
        """
        if len(sequence) < 3:
            raise ValueError(
                "Dynamic polymer generation requires sequence length >= 3 because Polymerist "
                "requires a non-empty middle sequence"
            )

    def _build_used_residue_names(
        self,
        monomer_names_used: dict[str, str],
        residue_names: dict[str, str] | None,
    ) -> dict[str, str]:
        """Build a metadata subset of residue names used by a sequence.

        Parameters
        ----------
        monomer_names_used : dict[str, str]
            Sequence-label-to-monomer-name mapping for labels used by a sequence.
        residue_names : dict[str, str] | None
            Optional monomer-name-to-residue-name mapping.

        Returns
        -------
        dict[str, str]
            Residue names keyed by monomer name for used monomers.
        """
        if residue_names is None:
            return {}

        used_monomers = set(monomer_names_used.values())
        return {
            monomer_name: residue_name
            for monomer_name, residue_name in sorted(residue_names.items())
            if monomer_name in used_monomers
        }

    def _build_monomer_group_fingerprint(self) -> dict[str, Any]:
        """Build a deterministic fingerprint for MonomerGroup contents.

        The fingerprint is based on sorted fragment names and deterministic
        chemistry representations where available, rather than object identity
        or Python's process-randomized ``hash()``.

        Returns
        -------
        dict[str, Any]
            Fingerprint payload with algorithm, digest, and fragment entries.
        """
        fragment_entries = {
            fragment_name: self._fingerprint_fragment(fragment)
            for fragment_name, fragment in sorted(self.monomer_group.monomers.items())
        }
        digest_payload = json.dumps(fragment_entries, sort_keys=True, separators=(",", ":"))

        return {
            "algorithm": "polyzymd-monomer-group-sha256-v2",
            "digest": hashlib.sha256(digest_payload.encode("utf-8")).hexdigest(),
            "fragments": fragment_entries,
        }

    def _fingerprint_fragment(self, fragment: Any) -> dict[str, Any]:
        """Build a deterministic fingerprint entry for one fragment.

        Parameters
        ----------
        fragment : Any
            Fragment chemistry object stored in a MonomerGroup.

        Returns
        -------
        dict[str, Any]
            Fragment fingerprint entry.
        """
        chemistry = self._fragment_chemistry_payload(fragment)
        digest_payload = json.dumps(chemistry, sort_keys=True, separators=(",", ":"))
        return {
            "algorithm": "polyzymd-fragment-sha256-v2",
            "digest": hashlib.sha256(digest_payload.encode("utf-8")).hexdigest(),
            "chemistry": chemistry,
        }

    def _fragment_chemistry_payload(self, fragment: Any) -> dict[str, Any]:
        """Convert a fragment object into deterministic chemistry metadata.

        Parameters
        ----------
        fragment : Any
            Fragment chemistry object stored in a MonomerGroup.

        Returns
        -------
        dict[str, Any]
            Deterministic chemistry payload suitable for JSON hashing.
        """
        if self._is_json_fingerprint_value(fragment):
            return {
                "kind": "json_value",
                "value": self._to_json_stable(fragment),
            }

        rdkit_payload = self._rdkit_mol_payload(fragment)
        if rdkit_payload is not None:
            return rdkit_payload

        for method_name in ("to_smiles", "to_smiles_explicit_hydrogens"):
            method = getattr(fragment, method_name, None)
            if method is None:
                continue
            smiles = self._call_smiles_method(method)
            if smiles:
                return {
                    "kind": method_name,
                    "value": self._canonicalize_smiles(smiles) or smiles,
                }

        to_rdkit = getattr(fragment, "to_rdkit", None)
        if to_rdkit is not None:
            try:
                rdkit_payload = self._rdkit_mol_payload(to_rdkit())
            except (TypeError, ValueError, RuntimeError):
                rdkit_payload = None
            if rdkit_payload is not None:
                return rdkit_payload

        public_state = self._stable_public_state(fragment)
        return {
            "kind": "public_state",
            "class": f"{type(fragment).__module__}.{type(fragment).__qualname__}",
            "value": public_state,
        }

    def _is_json_fingerprint_value(self, value: Any) -> bool:
        """Return whether a value can be fingerprinted as deterministic JSON.

        Parameters
        ----------
        value : Any
            Candidate fragment value.

        Returns
        -------
        bool
            True for recursively supported JSON-like values.
        """
        if value is None or isinstance(value, (bool, int, float, str)):
            return True
        if isinstance(value, (list, tuple)):
            return all(self._is_json_fingerprint_value(item) for item in value)
        if isinstance(value, dict):
            return all(
                self._is_json_fingerprint_value(key) and self._is_json_fingerprint_value(item)
                for key, item in value.items()
            )
        return False

    def _canonicalize_smiles(self, smiles: str) -> str | None:
        """Canonicalize a SMILES string when RDKit is available.

        Parameters
        ----------
        smiles : str
            SMILES string to canonicalize.

        Returns
        -------
        str | None
            Canonical isomeric SMILES, or None if canonicalization fails.
        """
        try:
            from rdkit import Chem
        except ImportError:
            return None

        try:
            mol = Chem.MolFromSmiles(smiles)
        except AttributeError:
            return None
        if mol is None:
            return None
        try:
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        except AttributeError:
            return None

    def _rdkit_mol_payload(self, fragment: Any) -> dict[str, str] | None:
        """Build a payload from an RDKit molecule-like object.

        Parameters
        ----------
        fragment : Any
            Potential RDKit molecule-like object.

        Returns
        -------
        dict[str, str] | None
            RDKit payload, or None when the object is not supported.
        """
        try:
            from rdkit import Chem
        except ImportError:
            return None

        if not hasattr(fragment, "GetNumAtoms"):
            return None

        try:
            smiles = Chem.MolToSmiles(fragment, canonical=True, isomericSmiles=True)
        except (AttributeError, TypeError, ValueError, RuntimeError):
            return None

        return {"kind": "rdkit_mol", "value": smiles}

    def _call_smiles_method(self, method: Any) -> str | None:
        """Call a fragment SMILES method using common OpenFF signatures.

        Parameters
        ----------
        method : Any
            Bound SMILES-producing method.

        Returns
        -------
        str | None
            SMILES string, or None if no supported call succeeds.
        """
        call_kwargs = (
            {"mapped": False, "isomeric": True, "explicit_hydrogens": True},
            {"isomeric": True, "explicit_hydrogens": True},
            {},
        )
        for kwargs in call_kwargs:
            try:
                smiles = method(**kwargs)
            except (TypeError, ValueError, RuntimeError):
                continue
            if isinstance(smiles, str) and smiles:
                return smiles
        return None

    def _stable_public_state(self, fragment: Any) -> Any:
        """Extract JSON-stable public state from a fragment object.

        Parameters
        ----------
        fragment : Any
            Fragment object without an available chemistry serializer.

        Returns
        -------
        Any
            JSON-stable public state used as a deterministic fallback.
        """
        try:
            state = vars(fragment)
        except TypeError:
            return str(type(fragment))

        public_state = {
            key: self._to_json_stable(value)
            for key, value in sorted(state.items())
            if not key.startswith("_")
        }
        return public_state

    def _to_json_stable(self, value: Any) -> Any:
        """Convert a value into a JSON-stable representation.

        Parameters
        ----------
        value : Any
            Value from public fragment state.

        Returns
        -------
        Any
            JSON-serializable representation without memory addresses.
        """
        if value is None or isinstance(value, (bool, int, float, str)):
            return value
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, dict):
            return {
                str(key): self._to_json_stable(item)
                for key, item in sorted(value.items(), key=lambda entry: str(entry[0]))
            }
        if isinstance(value, (list, tuple)):
            return [self._to_json_stable(item) for item in value]
        if isinstance(value, set):
            return sorted(self._to_json_stable(item) for item in value)
        return f"{type(value).__module__}.{type(value).__qualname__}"

    def _get_dynamic_cache_metadata_path(self, charged_sdf_path: Path) -> Path:
        """Get the metadata sidecar path for a dynamic cached charged SDF.

        Parameters
        ----------
        charged_sdf_path : Path
            Path to a charged dynamic-cache SDF file.

        Returns
        -------
        Path
            Path to the JSON metadata sidecar.
        """
        return charged_sdf_path.with_name(f"{charged_sdf_path.name}.metadata.json")

    def _dynamic_cache_metadata_matches(
        self,
        charged_sdf_path: Path,
        sequence: str,
        monomer_names: dict[str, str],
        residue_names: dict[str, str] | None = None,
    ) -> bool:
        """Return whether a dynamic charged SDF has matching metadata.

        Parameters
        ----------
        charged_sdf_path : Path
            Path to a charged dynamic-cache SDF file.
        sequence : str
            Polymer sequence string.
        monomer_names : dict[str, str]
            Mapping of sequence labels to monomer names.
        residue_names : dict[str, str] | None, optional
            Optional mapping of monomer names to residue names.

        Returns
        -------
        bool
            True if the metadata sidecar exists and exactly matches expectations.
        """
        metadata_path = self._get_dynamic_cache_metadata_path(charged_sdf_path)
        if not metadata_path.exists():
            logger.warning(
                "Ignoring dynamic polymer cache without metadata sidecar: %s",
                charged_sdf_path,
            )
            return False

        try:
            cached_metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            logger.warning(
                "Ignoring dynamic polymer cache with invalid metadata: %s", metadata_path
            )
            return False
        if not isinstance(cached_metadata, dict):
            logger.warning(
                "Ignoring dynamic polymer cache with non-object metadata: %s", metadata_path
            )
            return False

        expected_metadata = self._build_dynamic_cache_metadata(
            sequence,
            monomer_names,
            residue_names,
        )
        if cached_metadata != expected_metadata:
            mismatched_fields = sorted(
                field
                for field in set(cached_metadata) | set(expected_metadata)
                if cached_metadata.get(field) != expected_metadata.get(field)
            )
            logger.warning(
                "Ignoring dynamic polymer cache with mismatched metadata fields %s: %s",
                mismatched_fields,
                metadata_path,
            )
            return False

        return True

    def _write_dynamic_cache_metadata(
        self,
        charged_sdf_path: Path,
        sequence: str,
        monomer_names: dict[str, str],
        residue_names: dict[str, str] | None = None,
    ) -> None:
        """Write metadata for a dynamic cached charged SDF.

        Parameters
        ----------
        charged_sdf_path : Path
            Path to a charged dynamic-cache SDF file.
        sequence : str
            Polymer sequence string.
        monomer_names : dict[str, str]
            Mapping of sequence labels to monomer names.
        residue_names : dict[str, str] | None, optional
            Optional mapping of monomer names to residue names.
        """
        metadata_path = self._get_dynamic_cache_metadata_path(charged_sdf_path)
        metadata = self._build_dynamic_cache_metadata(sequence, monomer_names, residue_names)
        metadata_path.write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def _build_polymer_structure(
        self,
        sequence: str,
        monomer_names: dict[str, str],
        residue_names: dict[str, str] | None = None,
    ) -> tuple[Any, Path]:
        """Build a polymer 3D structure from sequence.

        The sequence parameter should be a string of block identifiers (e.g., "ABCAB")
        where each character maps to a monomer via the monomer_names dict.

        Polymerist's build_linear_polymer expects:

        - Terminal fragments set via monogrp.term_orient (head/tail)
        - A middle sequence of block identifiers (omitting head/tail)
        - An explicit sequence_map that maps middle labels to 2-site fragments

        Parameters
        ----------
        sequence : str
            Polymer sequence string, such as "ABCAB".
        monomer_names : dict[str, str]
            Mapping of sequence labels to monomer names.
        residue_names : dict[str, str] | None, optional
            Optional mapping of monomer names to 3-char residue names.

        Returns
        -------
        tuple[Any, Path]
            Tuple of the mBuild compound and path to the PDB file.

        Raises
        ------
        PolymerGenerationError
            If building fails after ``max_retries``.
        """
        self._validate_dynamic_sequence_length(sequence)
        _validate_dynamic_sequence_labels(sequence, monomer_names)

        # Parse sequence - head and tail are terminals, middle is for repeating
        head_label = sequence[0]
        tail_label = sequence[-1]
        middle_sequence = sequence[1:-1] if len(sequence) > 2 else ""

        head_monomer = monomer_names[head_label]
        tail_monomer = monomer_names[tail_label]

        # Set terminal orientations - these determine the 1-site end groups
        monogrp_local = _make_monomer_group(self.monomer_group.monomers)
        monogrp_local.term_orient = {
            "head": self._get_terminal_fragment_name(head_monomer),
            "tail": self._get_terminal_fragment_name(tail_monomer),
        }
        sequence_map = self._build_middle_sequence_map(middle_sequence, monomer_names)

        logger.debug(
            "Terminal orientations: head=%s, tail=%s",
            monogrp_local.term_orient["head"],
            monogrp_local.term_orient["tail"],
        )
        logger.debug(f"Middle sequence (block identifiers): {middle_sequence}")
        logger.debug(f"Middle sequence map: {sequence_map}")

        # Attempt building with retries for ring-piercing
        for attempt in range(self.max_retries):
            logger.debug(f"Building polymer attempt {attempt + 1}/{self.max_retries}")

            # Pin label-to-fragment mapping because Polymerist defaults are order-dependent
            chain = _build_linear_polymer(
                monomers=monogrp_local,
                n_monomers=len(sequence),
                sequence=middle_sequence,
                sequence_map=sequence_map,
                energy_minimize=True,
                allow_partial_sequences=True,
            )

            # Check for ring-piercing
            poly_mol = _mbmol_to_rdmol(chain)
            piercing_summary = _summarize_ring_piercing(poly_mol)

            if not piercing_summary:
                logger.info(f"Polymer built successfully on attempt {attempt + 1}")
                break
            else:
                logger.warning(
                    f"Ring-piercing detected on attempt {attempt + 1}: {piercing_summary}"
                )
        else:
            raise PolymerGenerationError(
                f"Failed to build polymer after {self.max_retries} attempts due to ring-piercing"
            )

        # Save PDB
        pdb_filename = self._make_polymer_filename(sequence, monomer_names, charged=False)
        pdb_path = self.cache_directory / f"{pdb_filename}.pdb"

        resname_map = self._build_resname_map(monomer_names, residue_names)
        _mbmol_to_openmm_pdb(pdb_path, chain, resname_map=resname_map)
        logger.info(f"Saved polymer PDB: {pdb_path}")

        return chain, pdb_path

    def _build_resname_map(
        self,
        monomer_names: dict[str, str],
        residue_names: dict[str, str] | None = None,
    ) -> dict[str, str]:
        """Build residue name mapping for PDB output.

        Parameters
        ----------
        monomer_names : dict[str, str]
            Mapping of labels to monomer names.
        residue_names : dict[str, str] | None, optional
            Optional custom residue names.

        Returns
        -------
        dict[str, str]
            Mapping of fragment names to 3-char residue names.
        """
        resname_map = {}
        for label, monomer_name in monomer_names.items():
            # Get 3-char residue name
            if residue_names and monomer_name in residue_names:
                base_resname = residue_names[monomer_name]
            else:
                base_resname = monomer_name[:3].upper()

            # Map both 1-site and 2-site fragments
            for frag_name in self.monomer_group.monomers.keys():
                if frag_name.startswith(monomer_name):
                    if "_1-site" in frag_name:
                        resname_map[frag_name] = f"{base_resname[:2]}1"
                    elif "_2-site" in frag_name:
                        resname_map[frag_name] = f"{base_resname[:2]}2"

        return resname_map

    def _make_polymer_filename(
        self,
        sequence: str,
        monomer_names: dict[str, str],
        charged: bool = True,
    ) -> str:
        """Create a descriptive filename for a polymer.

        Parameters
        ----------
        sequence : str
            Polymer sequence string.
        monomer_names : dict[str, str]
            Mapping of labels to monomer names.
        charged : bool, optional
            Whether this is a charged structure, by default True.

        Returns
        -------
        str
            Filename without extension.
        """
        unique_labels = _validate_dynamic_sequence_labels(sequence, monomer_names)

        # Build monomer prefix
        monomers_used = [monomer_names[label] for label in unique_labels]
        monomer_prefix = "-".join(monomers_used) if monomers_used else "NO_MONOMERS"

        length = len(sequence)
        filename = f"{monomer_prefix}_seq={sequence}_{length}-mer"
        if charged:
            filename += "_charged"

        return filename

    def generate_polymer(
        self,
        sequence: str,
        monomer_names: dict[str, str],
        residue_names: dict[str, str] | None = None,
        force_regenerate: bool = False,
    ) -> "OFFMolecule":
        """Generate a complete, charged polymer molecule.

        This method handles the full workflow:

        - Check cache for an existing charged SDF
        - Build a 3D structure if needed
        - Create an OpenFF topology
        - Assign partial charges
        - Cache the result

        Parameters
        ----------
        sequence : str
            Polymer sequence string, such as "ABCAB".
        monomer_names : dict[str, str]
            Mapping of sequence labels to monomer names.
        residue_names : dict[str, str] | None, optional
            Optional mapping of monomer names to 3-char residue names.
        force_regenerate : bool, optional
            If True, regenerate even if cached, by default False.

        Returns
        -------
        OFFMolecule
            OpenFF Molecule with partial charges assigned.

        Raises
        ------
        PolymerGenerationError
            If generation fails.
        """
        self._validate_dynamic_sequence_length(sequence)
        _validate_dynamic_sequence_labels(sequence, monomer_names)

        from openff.toolkit import Topology as OFFTopology

        # Check cache
        charged_filename = self._make_polymer_filename(sequence, monomer_names, charged=True)
        charged_sdf_path = self.cache_directory / f"{charged_filename}.sdf"

        if (
            charged_sdf_path.exists()
            and not force_regenerate
            and self._dynamic_cache_metadata_matches(
                charged_sdf_path,
                sequence,
                monomer_names,
                residue_names,
            )
        ):
            logger.info(f"Loading cached polymer from {charged_sdf_path}")
            off_top = _topology_from_sdf(charged_sdf_path)
            return _get_largest_offmol(off_top)

        # Build structure
        chain, pdb_path = self._build_polymer_structure(sequence, monomer_names, residue_names)

        # Create topology
        uncharged_filename = self._make_polymer_filename(sequence, monomer_names, charged=False)
        uncharged_sdf_path = self.cache_directory / f"{uncharged_filename}.sdf"

        logger.info("Creating OpenFF topology...")
        off_top = OFFTopology.from_pdb(
            str(pdb_path), _custom_substructures=self.monomer_group.monomers
        )
        was_partitioned = _partition(off_top)
        if not was_partitioned:
            raise PolymerGenerationError("Failed to partition polymer topology")

        # Fix residue names (truncate to 3 chars)
        for mol in off_top.molecules:
            for atom in mol.atoms:
                if "residue_name" in atom.metadata:
                    atom.metadata["extended_name"] = atom.metadata["residue_name"]
                    atom.metadata["residue_name"] = atom.metadata["residue_name"][:3]

        # Save uncharged SDF
        _topology_to_sdf(uncharged_sdf_path, off_top)
        logger.info(f"Saved uncharged SDF: {uncharged_sdf_path}")

        # Assign partial charges
        logger.info(f"Assigning partial charges using {self.charger_type}...")
        off_mol = _get_largest_offmol(off_top)
        charged_mol = self.charger.charge_molecule(off_mol)

        # Save charged SDF
        charged_top = charged_mol.to_topology()
        _topology_to_sdf(charged_sdf_path, charged_top)
        self._write_dynamic_cache_metadata(charged_sdf_path, sequence, monomer_names, residue_names)
        logger.info(f"Saved charged SDF: {charged_sdf_path}")

        return charged_mol

    def generate_polymers_batch(
        self,
        sequences: list[str],
        monomer_names: dict[str, str],
        residue_names: dict[str, str] | None = None,
    ) -> dict[str, "OFFMolecule"]:
        """Generate multiple polymers (sequential processing).

        Parameters
        ----------
        sequences : list[str]
            List of polymer sequences to generate.
        monomer_names : dict[str, str]
            Mapping of sequence labels to monomer names.
        residue_names : dict[str, str] | None, optional
            Optional mapping of monomer names to 3-char residue names.

        Returns
        -------
        dict[str, OFFMolecule]
            Dictionary mapping sequence to OpenFF Molecule.
        """
        results = {}
        for i, sequence in enumerate(sequences):
            logger.info(f"Generating polymer {i + 1}/{len(sequences)}: {sequence}")
            try:
                mol = self.generate_polymer(sequence, monomer_names, residue_names)
                results[sequence] = mol
            except Exception as e:
                logger.error(f"Failed to generate polymer {sequence}: {e}")
                raise

        return results

    def get_cached_polymer(
        self,
        sequence: str,
        monomer_names: dict[str, str],
        residue_names: dict[str, str] | None = None,
    ) -> "OFFMolecule" | None:
        """Load a polymer from cache if it exists.

        Parameters
        ----------
        sequence : str
            Polymer sequence string.
        monomer_names : dict[str, str]
            Mapping of labels to monomer names.
        residue_names : dict[str, str] | None, optional
            Optional mapping of monomer names to residue names.

        Returns
        -------
        OFFMolecule | None
            OpenFF Molecule if cached, otherwise None.
        """
        self._validate_dynamic_sequence_length(sequence)
        _validate_dynamic_sequence_labels(sequence, monomer_names)

        charged_filename = self._make_polymer_filename(sequence, monomer_names, charged=True)
        charged_sdf_path = self.cache_directory / f"{charged_filename}.sdf"

        if charged_sdf_path.exists() and self._dynamic_cache_metadata_matches(
            charged_sdf_path,
            sequence,
            monomer_names,
            residue_names,
        ):
            off_top = _topology_from_sdf(charged_sdf_path)
            return _get_largest_offmol(off_top)

        return None
