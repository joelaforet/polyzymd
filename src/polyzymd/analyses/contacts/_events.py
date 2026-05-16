"""Sparse contact-event utilities for the contacts analysis."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class ResidueIdentity:
    """Stable residue identity for contact sidecars.

    Parameters
    ----------
    key : str
        Stable chain-aware residue key.
    resid : int
        Residue identifier from the topology.
    resname : str
        Residue name from the topology.
    chain_id : str
        Chain or segment identifier, or ``"?"`` when unavailable.
    residue_index : int
        MDAnalysis residue index when available.
    group : str | None, optional
        Contacts grouping label for protein residues.
    """

    key: str
    resid: int
    resname: str
    chain_id: str
    residue_index: int
    group: str | None = None


class _IdentityResidueGrouping:
    """Group residues by residue name without biochemical coarsening."""

    @property
    def available_groups(self) -> list[str]:
        """Return an empty group list for open-ended residue-name labels."""

        return []

    def classify(self, resname: str) -> str:
        """Return the normalized residue name as the group label."""

        normalized = str(resname).upper().strip()
        return normalized or "unknown"


class _SecondaryStructureGrouping:
    """Group residues by attached secondary-structure annotations."""

    _DSSP_LABELS = {
        "H": "helix",
        "G": "helix",
        "I": "helix",
        "E": "sheet",
        "B": "sheet",
        "T": "turn",
        "S": "turn",
        "C": "coil",
        "-": "coil",
        "": "coil",
    }

    @property
    def available_groups(self) -> list[str]:
        """Return supported secondary-structure group labels."""

        return ["helix", "sheet", "turn", "coil", "unknown"]

    def classify(self, resname: str) -> str:
        """Classify residue names without structural annotations."""

        del resname
        return "unknown"

    def classify_residue(self, residue: Any) -> str:
        """Classify one residue from common secondary-structure attributes."""

        raw_label = self._extract_secondary_structure_label(residue)
        if raw_label is None:
            return "unknown"
        label = str(raw_label).strip()
        if not label:
            return "coil"
        return self._DSSP_LABELS.get(label.upper(), label.lower())

    @staticmethod
    def _extract_secondary_structure_label(residue: Any) -> Any | None:
        """Extract a secondary-structure label from common residue attributes."""

        for attr in ("secondary_structure", "secstruct", "ss", "dssp"):
            value = getattr(residue, attr, None)
            if value is not None:
                return value

        atoms = getattr(residue, "atoms", None)
        if atoms is None:
            return None
        for attr in ("secondary_structures", "secstructs", "ss", "dssp"):
            values = getattr(atoms, attr, None)
            if values is None:
                continue
            try:
                if len(values) == 0:
                    continue
                return values[0]
            except TypeError:
                return values
        return None


def build_contact_grouping(grouping: str) -> Any:
    """Build a residue grouping strategy from contacts settings."""

    from polyzymd.analyses.shared.groupings import ProteinAAClassification

    if grouping == "aa_class":
        return ProteinAAClassification()
    if grouping == "secondary_structure":
        return _SecondaryStructureGrouping()
    if grouping == "none":
        return _IdentityResidueGrouping()
    raise ValueError(f"Unsupported contacts grouping mode: {grouping}")


def classify_residue(residue: Any, grouping: Any) -> str:
    """Return the contacts group label for one protein residue."""

    classify_residue_method = getattr(grouping, "classify_residue", None)
    if callable(classify_residue_method):
        return str(classify_residue_method(residue))
    return str(grouping.classify(getattr(residue, "resname", "")))


def residue_identity(residue: Any, *, group: str | None = None) -> ResidueIdentity:
    """Build a stable chain-aware identity row for one residue."""

    residue_index = _residue_index(residue)
    resid = int(getattr(residue, "resid", residue_index))
    resname = str(getattr(residue, "resname", "UNK"))
    chain_id = _chain_id(residue)
    return ResidueIdentity(
        key=f"{chain_id}:{residue_index}:{resid}:{resname}",
        resid=resid,
        resname=resname,
        chain_id=chain_id,
        residue_index=residue_index,
        group=group,
    )


def unique_residue_key(residue: Any) -> int | tuple[str, int | Any, str]:
    """Build a collision-safe key for mapping residue-like objects."""

    residue_ix = getattr(residue, "ix", None)
    if residue_ix is not None:
        try:
            return int(residue_ix)
        except (TypeError, ValueError):
            pass
    return (_chain_id(residue), getattr(residue, "resid", None), getattr(residue, "resname", ""))


def build_atom_to_residue_map(atoms: Any, residues: Any) -> NDArray[np.int64]:
    """Map each atom in an atom group to its selected residue index."""

    residue_lookup = {unique_residue_key(residue): index for index, residue in enumerate(residues)}
    atom_to_residue = np.zeros(len(atoms), dtype=np.int64)
    for atom_index, atom in enumerate(atoms):
        residue_key = unique_residue_key(atom.residue)
        atom_to_residue[atom_index] = residue_lookup[residue_key]
    return atom_to_residue


def identify_polymer_chains(query_atoms: Any) -> tuple[NDArray[np.int64], list[str]]:
    """Assign selected polymer residues to chain indices using fragments.

    Parameters
    ----------
    query_atoms : Any
        MDAnalysis atom group containing polymer atoms.

    Returns
    -------
    tuple of ndarray and list of str
        Polymer residue chain indices aligned to ``query_atoms.residues`` and
        warning messages emitted while resolving fragments.
    """

    query_residues = query_atoms.residues
    chain_indices = np.zeros(len(query_residues), dtype=np.int64)
    warnings: list[str] = []
    fragments = fragments_or_single(
        query_atoms.atoms,
        context="contacts polymer chain detection",
        warnings=warnings,
    )
    residue_lookup = {
        unique_residue_key(residue): idx for idx, residue in enumerate(query_residues)
    }
    for chain_idx, fragment in enumerate(fragments):
        for residue in getattr(fragment, "residues", []):
            query_idx = residue_lookup.get(unique_residue_key(residue))
            if query_idx is not None:
                chain_indices[query_idx] = chain_idx
    return chain_indices, warnings


def fragments_or_single(
    atom_group: Any, *, context: str, warnings: list[str] | None = None
) -> list[Any]:
    """Return MDAnalysis fragments with a no-bond topology fallback."""

    from MDAnalysis.exceptions import NoDataError

    try:
        fragments = atom_group.fragments
    except NoDataError:
        message = (
            f"{context}: topology has no bond information; treating the selected polymer as "
            "one fragment"
        )
        LOGGER.warning("%s", message)
        if warnings is not None:
            warnings.append(message)
        return [atom_group]
    return list(fragments) if fragments else [atom_group]


def _residue_index(residue: Any) -> int:
    """Return a stable integer residue index for one residue-like object."""

    for attr in ("ix", "resindex"):
        value = getattr(residue, attr, None)
        if value is not None:
            try:
                return int(value)
            except (TypeError, ValueError):
                continue
    return int(getattr(residue, "resid", 0))


def _chain_id(residue: Any) -> str:
    """Return the best available chain or segment identifier for a residue."""

    for attr in ("chainID", "chainid", "segid"):
        value = getattr(residue, attr, None)
        if value is not None and str(value).strip():
            return str(value).strip()
    atoms = getattr(residue, "atoms", None)
    if atoms is not None:
        for attr in ("chainIDs", "chainids", "segids"):
            values = getattr(atoms, attr, None)
            if values is None:
                continue
            try:
                if len(values) > 0 and str(values[0]).strip():
                    return str(values[0]).strip()
            except TypeError:
                if str(values).strip():
                    return str(values).strip()
    return "?"
