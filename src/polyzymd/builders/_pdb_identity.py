"""Internal helpers for canonical PDB identity metadata."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

CHAIN_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
PROTEIN_CHAIN = "A"
SUBSTRATE_CHAIN = "B"
POLYMER_CHAIN = "C"
SOLVENT_START_CHAIN_INDEX = 3
PDB_MAX_RESIDUE_ID = 9999
PDB_MAX_ATOM_SERIAL = 99999
WATER_RESIDUE_ALIASES = {"HOH", "WAT", "H2O", "SOL", "TIP3", "TIP3P"}
ION_RESIDUE_ALIASES = {
    "NA": "NA",
    "NA+": "NA",
    "SOD": "NA",
    "K": "K",
    "K+": "K",
    "POT": "K",
    "CL": "CL",
    "CL-": "CL",
    "CLA": "CL",
    "MG": "MG",
    "MG2": "MG",
    "MG2+": "MG",
    "CA": "CA",
    "CA2": "CA",
    "CA2+": "CA",
}
_ORIGINAL_RESIDUE_TOKEN_KEY = "_polyzymd_original_residue_token"


def normalize_topology_pdb_identifiers(
    topology: Any,
    *,
    n_enzyme_molecules: int = 0,
    n_substrate_molecules: int = 0,
    n_polymer_chains: int = 0,
    preserve_enzyme_chain_ids: bool = False,
) -> None:
    """Normalize OpenFF atom metadata to the PolyzyMD chain convention.

    Parameters
    ----------
    topology : Any
        OpenFF-like topology whose atom metadata should be normalized in place.
    n_enzyme_molecules : int, optional
        Number of leading protein molecules, by default 0.
    n_substrate_molecules : int, optional
        Number of substrate molecules after protein molecules, by default 0.
    n_polymer_chains : int, optional
        Number of polymer or glycan molecules after substrates, by default 0.
    preserve_enzyme_chain_ids : bool, optional
        Preserve existing protein chain IDs when available, by default False.
    """

    if not hasattr(topology, "n_molecules") or not hasattr(topology, "molecule"):
        return

    mol_idx = 0
    protein_residue_num = 1
    for _ in range(n_enzyme_molecules):
        mol = topology.molecule(mol_idx)
        protein_residue_num = _assign_protein_metadata(
            mol,
            preserve_enzyme_chain_ids=preserve_enzyme_chain_ids,
            start_residue_number=protein_residue_num,
        )
        mol_idx += 1

    for _ in range(n_substrate_molecules):
        mol = topology.molecule(mol_idx)
        _assign_single_residue_metadata(mol, chain_id=SUBSTRATE_CHAIN, residue_number=1)
        mol_idx += 1

    polymer_residue_num = _next_residue_number_for_chain(
        topology,
        POLYMER_CHAIN,
        end_mol_idx=mol_idx,
        preserve_enzyme_chain_ids=preserve_enzyme_chain_ids,
    )
    for _ in range(n_polymer_chains):
        mol = topology.molecule(mol_idx)
        polymer_residue_num = _assign_polymer_metadata(
            mol, start_residue_number=polymer_residue_num
        )
        mol_idx += 1

    _assign_solvent_metadata(topology, start_mol_idx=mol_idx)


def require_classic_pdb_atom_capacity(topology: Any) -> None:
    """Raise if classic PDB output cannot represent the topology atom count.

    Parameters
    ----------
    topology : Any
        Topology-like object with atoms available through ``atoms()`` or
        ``n_atoms``.

    Raises
    ------
    ValueError
        If the topology exceeds the classic PDB atom serial limit.
    """

    atom_count = _topology_atom_count(topology)
    if atom_count > PDB_MAX_ATOM_SERIAL:
        raise ValueError(
            "Classic PDB output supports at most 99999 atoms without invalid serial wrapping; "
            f"this topology has {atom_count} atoms. Write mmCIF for OpenMM workflows or GRO "
            "for GROMACS workflows instead of classic PDB."
        )


def _assign_protein_metadata(
    molecule: Any,
    *,
    preserve_enzyme_chain_ids: bool,
    start_residue_number: int,
) -> int:
    """Assign protein metadata and return the next available residue number."""

    residue_map: dict[str, int] = {}
    next_residue_number = start_residue_number
    for atom_index, atom in enumerate(molecule.atoms):
        if preserve_enzyme_chain_ids:
            atom.metadata.setdefault("chain_id", PROTEIN_CHAIN)
            if "residue_number" in atom.metadata:
                atom.metadata["residue_number"] = str(atom.metadata["residue_number"])
        else:
            residue_token = _original_residue_token(atom, atom_index)
            atom.metadata["chain_id"] = PROTEIN_CHAIN
            if residue_token not in residue_map:
                residue_map[residue_token] = next_residue_number
                next_residue_number += 1
            atom.metadata["residue_number"] = str(residue_map[residue_token])
        atom.metadata["insertion_code"] = str(atom.metadata.get("insertion_code", "") or "")
    return next_residue_number


def _original_residue_token(atom: Any, atom_index: int) -> str:
    """Return an idempotent grouping token for one source protein residue."""

    metadata = atom.metadata
    stored = metadata.get(_ORIGINAL_RESIDUE_TOKEN_KEY)
    if stored is not None:
        return str(stored)
    residue_number = metadata.get("residue_number")
    residue_name = metadata.get("residue_name") or metadata.get("residue_name_3") or ""
    insertion_code = metadata.get("insertion_code") or ""
    chain_id = metadata.get("chain_id") or ""
    if residue_number is None:
        token = f"{chain_id!s}|missing-residue-metadata:{atom_index}"
    else:
        token = f"{chain_id!s}|{residue_number!s}|{residue_name!s}|{insertion_code!s}"
    metadata[_ORIGINAL_RESIDUE_TOKEN_KEY] = token
    return token


def _assign_single_residue_metadata(molecule: Any, *, chain_id: str, residue_number: int) -> None:
    """Assign a single residue identity to every atom in a molecule."""

    for atom in molecule.atoms:
        atom.metadata["chain_id"] = chain_id
        atom.metadata["residue_number"] = str(residue_number)
        atom.metadata["insertion_code"] = ""


def _assign_polymer_metadata(molecule: Any, *, start_residue_number: int) -> int:
    """Assign chain C sequential residue numbers to a polymer-like molecule."""

    residue_num = start_residue_number
    current_monomer_residue = None
    for atom in molecule.atoms:
        atom.metadata["chain_id"] = POLYMER_CHAIN
        atom.metadata["insertion_code"] = str(atom.metadata.get("insertion_code", "") or "")
        atom_residue = atom.metadata.get("residue_number", "0")
        if atom_residue != current_monomer_residue:
            if current_monomer_residue is not None:
                residue_num += 1
            current_monomer_residue = atom_residue
        atom.metadata["residue_number"] = str(residue_num)
    return residue_num + 1


def _assign_solvent_metadata(topology: Any, *, start_mol_idx: int) -> None:
    """Assign unique D-Z solvent molecule identities."""

    chain_idx = SOLVENT_START_CHAIN_INDEX
    residue_num = 1
    for mol_idx in range(start_mol_idx, topology.n_molecules):
        if residue_num > PDB_MAX_RESIDUE_ID:
            chain_idx += 1
            residue_num = 1
        if chain_idx >= len(CHAIN_LETTERS):
            capacity = (len(CHAIN_LETTERS) - SOLVENT_START_CHAIN_INDEX) * PDB_MAX_RESIDUE_ID
            raise ValueError(
                "Classic PDB chain/residue capacity exceeded for solvent/ion molecules: "
                f"requires more than {capacity} molecules across chains D-Z. Write mmCIF for "
                "OpenMM workflows or GRO for GROMACS workflows instead of classic PDB."
            )
        mol = topology.molecule(mol_idx)
        chain_id = CHAIN_LETTERS[chain_idx]
        _assign_solvent_molecule_metadata(mol, chain_id=chain_id, residue_number=residue_num)
        residue_num += 1


def _assign_solvent_molecule_metadata(molecule: Any, *, chain_id: str, residue_number: int) -> None:
    """Assign canonical solvent, co-solvent, or ion metadata to one molecule."""

    residue_name = _canonical_solvent_residue_name(molecule)
    atom_names = _canonical_solvent_atom_names(molecule, residue_name)
    for atom, atom_name in zip(molecule.atoms, atom_names, strict=True):
        atom.metadata["chain_id"] = chain_id
        atom.metadata["residue_number"] = str(residue_number)
        atom.metadata["residue_name"] = residue_name
        atom.metadata["insertion_code"] = ""
        atom.metadata["atom_name"] = atom_name


def _canonical_solvent_residue_name(molecule: Any) -> str:
    """Return the canonical residue name for a solvent-like molecule."""

    atoms = tuple(molecule.atoms)
    residue_names = {
        str(atom.metadata.get("residue_name", atom.metadata.get("residue", "")) or "")
        .strip()
        .upper()
        for atom in atoms
    }
    residue_names.discard("")
    if _is_water_molecule(atoms, residue_names):
        return "HOH"
    if len(atoms) == 1:
        alias = next(iter(residue_names), "")
        element_symbol = _atom_element_symbol(atoms[0])
        return ION_RESIDUE_ALIASES.get(
            alias, ION_RESIDUE_ALIASES.get(element_symbol, element_symbol)
        )
    return (next(iter(residue_names), getattr(molecule, "name", "MOL")) or "MOL")[:3]


def _canonical_solvent_atom_names(molecule: Any, residue_name: str) -> tuple[str, ...]:
    """Return canonical atom names for a solvent-like molecule without reordering atoms."""

    atoms = tuple(molecule.atoms)
    if residue_name == "HOH" and _is_water_molecule(atoms, {residue_name}):
        hydrogen_seen = 0
        names = []
        for atom in atoms:
            if int(atom.atomic_number) == 8:
                names.append("O")
            elif int(atom.atomic_number) == 1:
                hydrogen_seen += 1
                names.append(f"H{hydrogen_seen}")
            else:
                names.append(_metadata_atom_name(atom))
        return tuple(names)
    if len(atoms) == 1 and residue_name in set(ION_RESIDUE_ALIASES.values()):
        return (residue_name,)
    return tuple(_metadata_atom_name(atom) for atom in atoms)


def _is_water_molecule(atoms: tuple[Any, ...], residue_names: set[str]) -> bool:
    """Return whether atom composition and metadata indicate a water molecule."""

    atomic_numbers = sorted(int(atom.atomic_number) for atom in atoms)
    return atomic_numbers == [1, 1, 8] and (
        not residue_names or bool(residue_names & WATER_RESIDUE_ALIASES)
    )


def _metadata_atom_name(atom: Any) -> str:
    """Return a nonblank atom name from metadata, atom name, or element symbol."""

    for value in (
        atom.metadata.get("atom_name"),
        getattr(atom, "name", None),
        _atom_element_symbol(atom),
    ):
        name = str(value or "").strip().upper()
        if name:
            return name
    return "X"


def _atom_element_symbol(atom: Any) -> str:
    """Return a conservative element symbol for an atom-like object."""

    symbol = str(getattr(atom, "symbol", "") or "").strip().upper()
    if symbol:
        return symbol
    atomic_number = int(atom.atomic_number)
    return {1: "H", 6: "C", 7: "N", 8: "O", 11: "NA", 12: "MG", 17: "CL", 19: "K", 20: "CA"}.get(
        atomic_number, "X"
    )


def _next_residue_number_for_chain(
    topology: Any, chain_id: str, *, end_mol_idx: int, preserve_enzyme_chain_ids: bool
) -> int:
    """Return the next residue number after preserved residues on a chain."""

    if not preserve_enzyme_chain_ids:
        return 1
    max_residue_number = 0
    for mol_idx in range(end_mol_idx):
        mol = topology.molecule(mol_idx)
        for atom in mol.atoms:
            if str(atom.metadata.get("chain_id", "")).upper() != chain_id:
                continue
            residue_number = _metadata_residue_number(atom.metadata)
            if residue_number is not None:
                max_residue_number = max(max_residue_number, residue_number)
    return max_residue_number + 1


def _metadata_residue_number(metadata: dict[str, Any]) -> int | None:
    """Return an integer residue number from OpenFF metadata when possible."""

    try:
        return int(str(metadata.get("residue_number", "")).strip())
    except ValueError:
        return None


def _topology_atom_count(topology: Any) -> int:
    """Return the atom count for an OpenFF or OpenMM topology-like object."""

    n_atoms = getattr(topology, "n_atoms", None)
    if n_atoms is not None:
        return int(n_atoms)
    atoms = getattr(topology, "atoms", None)
    if callable(atoms):
        return sum(1 for _ in atoms())
    if isinstance(atoms, Iterable):
        return sum(1 for _ in atoms)
    raise ValueError("Cannot determine topology atom count for classic PDB capacity check")
