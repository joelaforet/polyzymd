"""Source-to-canonical identity mapping for prepared protein PDB files."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.structure.parsing import (
    ATOM_RECORD_PREFIXES,
    PdbAtomRecord,
    parse_pdb_atom_records,
    parse_pdb_conect_pairs,
)


class PreparedProteinAtomIdentity(BaseModel):
    """Stable PDB identity for one atom at one preparation stage."""

    chain_id: str
    residue_name: str
    residue_number: int
    insertion_code: str = ""
    atom_name: str
    serial: int
    atom_index: int = Field(..., ge=0)


class PreparedProteinIdentityRecord(BaseModel):
    """Mapping record for a surviving, missing, or preparation-added atom."""

    status: Literal["survived", "missing_after_preparation", "prepared_added"]
    source: PreparedProteinAtomIdentity | None = None
    prepared: PreparedProteinAtomIdentity | None = None
    canonical: PreparedProteinAtomIdentity | None = None


class PreparedProteinIdentityResult(BaseModel):
    """Canonical construction PDB and its complete atom identity mapping."""

    source_path: Path
    prepared_path: Path
    canonical_path: Path
    mapping_path: Path
    records: tuple[PreparedProteinIdentityRecord, ...]
    canonical_conect_pair_count: int

    def canonical_for_source(self, selector: Any) -> PreparedProteinAtomIdentity:
        """Return the unique canonical survivor for a source selector."""
        matches = [
            record
            for record in self.records
            if record.source is not None and _source_identity_matches(record.source, selector)
        ]
        if len(matches) != 1:
            raise ValueError(
                "Source protein selector must resolve exactly one atom in the identity map; "
                f"found {len(matches)}"
            )
        canonical = matches[0].canonical
        if canonical is None:
            raise ValueError("Selected source protein atom did not survive protein preparation")
        return canonical


def canonicalize_prepared_protein_identity(
    source_path: Path | str,
    prepared_path: Path | str,
    canonical_path: Path | str,
) -> PreparedProteinIdentityResult:
    """Write canonical chain-A protein identity and a source mapping sidecar."""
    source = Path(source_path)
    prepared = Path(prepared_path)
    destination = Path(canonical_path)
    sidecar = destination.with_suffix(".identity.json")
    source_atoms = tuple(parse_pdb_atom_records(source, require_atoms=True))
    prepared_atoms = tuple(parse_pdb_atom_records(prepared, require_atoms=True))
    _require_unambiguous_atoms(source_atoms, label="source")
    _require_unambiguous_atoms(prepared_atoms, label="prepared")

    source_boundary_residues = _boundary_residue_keys(source, source_atoms)
    prepared_boundary_indices = _boundary_atom_indices(prepared)
    prepared_atom_lines = _atom_lines(prepared)

    canonical_by_prepared_index: dict[int, PreparedProteinAtomIdentity] = {}
    canonical_lines: list[str] = []
    residue_number_by_key: dict[tuple[str, str, int, str], int] = {}
    next_residue_number = 1
    previous_residue_key = None
    for atom in prepared_atoms:
        prepared_index = _required_atom_index(atom, prepared)
        residue_key = _residue_key(atom)
        if residue_key not in residue_number_by_key:
            residue_number_by_key[residue_key] = next_residue_number
            next_residue_number += 1
        if prepared_index > 0 and (
            prepared_index in prepared_boundary_indices
            or (residue_key != previous_residue_key and residue_key in source_boundary_residues)
        ):
            if not canonical_lines or canonical_lines[-1] != "TER\n":
                canonical_lines.append("TER\n")
        canonical_identity = PreparedProteinAtomIdentity(
            chain_id="A",
            residue_name=atom.residue_name,
            residue_number=residue_number_by_key[residue_key],
            insertion_code="",
            atom_name=atom.atom_name,
            serial=prepared_index + 1,
            atom_index=prepared_index,
        )
        canonical_by_prepared_index[prepared_index] = canonical_identity
        previous_residue_key = residue_key
        canonical_lines.append(
            _rewrite_atom_line(prepared_atom_lines[prepared_index], canonical_identity)
        )

    if canonical_lines[-1] != "TER\n":
        canonical_lines.append("TER\n")

    records = _mapping_records(
        source_atoms,
        prepared_atoms,
        canonical_by_prepared_index=canonical_by_prepared_index,
    )
    canonical_pairs = _canonical_conect_pairs(
        source,
        prepared,
        records,
    )
    canonical_lines.extend(_format_conect_lines(canonical_pairs))
    canonical_lines.append("END\n")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(canonical_lines), encoding="utf-8")

    result = PreparedProteinIdentityResult(
        source_path=source,
        prepared_path=prepared,
        canonical_path=destination,
        mapping_path=sidecar,
        records=records,
        canonical_conect_pair_count=len(canonical_pairs),
    )
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(json.dumps(result.model_dump(mode="json"), indent=2) + "\n")
    return result


def _mapping_records(
    source_atoms: tuple[PdbAtomRecord, ...],
    prepared_atoms: tuple[PdbAtomRecord, ...],
    *,
    canonical_by_prepared_index: dict[int, PreparedProteinAtomIdentity],
) -> tuple[PreparedProteinIdentityRecord, ...]:
    source_by_key = {_atom_key(atom): atom for atom in source_atoms}
    prepared_by_key = {_atom_key(atom): atom for atom in prepared_atoms}
    records: list[PreparedProteinIdentityRecord] = []
    for source_atom in source_atoms:
        prepared_atom = prepared_by_key.get(_atom_key(source_atom))
        canonical = (
            canonical_by_prepared_index[_required_atom_index(prepared_atom, Path("prepared"))]
            if prepared_atom is not None
            else None
        )
        records.append(
            PreparedProteinIdentityRecord(
                status="survived" if prepared_atom is not None else "missing_after_preparation",
                source=_identity(source_atom),
                prepared=_identity(prepared_atom) if prepared_atom is not None else None,
                canonical=canonical,
            )
        )
    for prepared_atom in prepared_atoms:
        if _atom_key(prepared_atom) in source_by_key:
            continue
        prepared_index = _required_atom_index(prepared_atom, Path("prepared"))
        records.append(
            PreparedProteinIdentityRecord(
                status="prepared_added",
                prepared=_identity(prepared_atom),
                canonical=canonical_by_prepared_index[prepared_index],
            )
        )
    return tuple(records)


def _canonical_conect_pairs(
    source: Path,
    prepared: Path,
    records: tuple[PreparedProteinIdentityRecord, ...],
) -> set[tuple[int, int]]:
    source_serial_to_canonical = {
        record.source.serial: record.canonical.serial
        for record in records
        if record.source is not None and record.canonical is not None
    }
    prepared_serial_to_canonical = {
        record.prepared.serial: record.canonical.serial
        for record in records
        if record.prepared is not None and record.canonical is not None
    }
    pairs: set[tuple[int, int]] = set()
    for path, lookup in (
        (source, source_serial_to_canonical),
        (prepared, prepared_serial_to_canonical),
    ):
        for left, right in _undirected_pairs(parse_pdb_conect_pairs(path)):
            if left in lookup and right in lookup:
                pairs.add(tuple(sorted((lookup[left], lookup[right]))))
    return pairs


def _require_unambiguous_atoms(atoms: tuple[PdbAtomRecord, ...], *, label: str) -> None:
    atom_keys: set[tuple[str, str, int, str, str]] = set()
    seen_residues: set[tuple[str, str, int, str]] = set()
    previous_residue = None
    for atom in atoms:
        residue_key = _residue_key(atom)
        if residue_key != previous_residue:
            if residue_key in seen_residues:
                raise ValueError(f"Ambiguous {label} PDB residue identity: {residue_key}")
            seen_residues.add(residue_key)
            previous_residue = residue_key
        atom_key = _atom_key(atom)
        if atom_key in atom_keys:
            raise ValueError(f"Ambiguous {label} PDB atom identity: {atom_key}")
        atom_keys.add(atom_key)


def _boundary_atom_indices(path: Path) -> set[int]:
    boundaries: set[int] = set()
    atom_index = 0
    previous_chain: str | None = None
    pending = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("TER"):
            pending = atom_index > 0
        elif line.startswith(ATOM_RECORD_PREFIXES):
            chain = line[21:22].strip()
            if atom_index > 0 and (pending or chain != previous_chain):
                boundaries.add(atom_index)
            previous_chain = chain
            pending = False
            atom_index += 1
    return boundaries


def _boundary_residue_keys(
    path: Path, atoms: tuple[PdbAtomRecord, ...]
) -> set[tuple[str, str, int, str]]:
    return {_residue_key(atoms[index]) for index in _boundary_atom_indices(path)}


def _atom_lines(path: Path) -> tuple[str, ...]:
    return tuple(
        line
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
        if line.startswith(ATOM_RECORD_PREFIXES)
    )


def _rewrite_atom_line(line: str, identity: PreparedProteinAtomIdentity) -> str:
    body = line.rstrip("\r\n").ljust(80)
    return (
        f"{body[:6]}{identity.serial:5d}{body[11:21]}A"
        f"{identity.residue_number:4d} {body[27:]}\n"
    )


def _format_conect_lines(pairs: set[tuple[int, int]]) -> list[str]:
    neighbors: dict[int, set[int]] = defaultdict(set)
    for left, right in pairs:
        neighbors[left].add(right)
        neighbors[right].add(left)
    return [
        f"CONECT{serial:5d}{''.join(f'{neighbor:5d}' for neighbor in sorted(entries))}\n"
        for serial, entries in sorted(neighbors.items())
    ]


def _undirected_pairs(pairs: Any) -> set[tuple[int, int]]:
    return {tuple(sorted((int(left), int(right)))) for left, right in pairs if left != right}


def _identity(atom: PdbAtomRecord) -> PreparedProteinAtomIdentity:
    if atom.serial is None:
        raise ValueError(f"PDB atom {atom.atom_name} is missing a serial")
    return PreparedProteinAtomIdentity(
        chain_id=atom.chain_id,
        residue_name=atom.residue_name,
        residue_number=atom.residue_number,
        insertion_code=atom.insertion_code or "",
        atom_name=atom.atom_name,
        serial=atom.serial,
        atom_index=_required_atom_index(atom, Path("PDB")),
    )


def _required_atom_index(atom: PdbAtomRecord, path: Path) -> int:
    if atom.atom_index is None:
        raise ValueError(f"PDB atom {atom.atom_name} is missing atom_index in {path}")
    return atom.atom_index


def _residue_key(atom: PdbAtomRecord) -> tuple[str, str, int, str]:
    return (atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code or "")


def _atom_key(atom: PdbAtomRecord) -> tuple[str, str, int, str, str]:
    return (*_residue_key(atom), atom.atom_name)


def _source_identity_matches(identity: PreparedProteinAtomIdentity, selector: Any) -> bool:
    return all(
        expected is None or observed == expected
        for observed, expected in (
            (identity.chain_id, getattr(selector, "chain_id", None)),
            (identity.residue_name, getattr(selector, "residue_name", None)),
            (identity.residue_number, getattr(selector, "residue_number", None)),
            (identity.insertion_code, getattr(selector, "insertion_code", None)),
            (identity.atom_name, getattr(selector, "atom_name", None)),
            (identity.serial, getattr(selector, "atom_serial", None)),
            (identity.atom_index, getattr(selector, "atom_index", None)),
        )
    )
