"""Dependency-light PDB atom and connectivity parsing utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field

ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")


class PdbAtomRecord(BaseModel):
    """Fixed-width PDB atom record data used by conjugation builders."""

    serial: int | None = Field(None, description="Input or output PDB atom serial")
    atom_index: int | None = Field(None, description="Zero-based source atom index")
    atom_name: str = Field(..., min_length=1, description="PDB atom name")
    residue_name: str = Field(..., min_length=1, description="PDB residue name")
    chain_id: str = Field("", max_length=1, description="PDB chain identifier")
    residue_number: int = Field(..., description="PDB residue sequence number")
    insertion_code: str = Field("", max_length=1, description="PDB insertion code")
    x: float = Field(..., description="X coordinate in angstrom")
    y: float = Field(..., description="Y coordinate in angstrom")
    z: float = Field(..., description="Z coordinate in angstrom")
    occupancy: float = Field(1.0, description="PDB occupancy")
    temp_factor: float = Field(0.0, description="PDB temperature factor")
    element: str = Field("", description="Element symbol")
    charge: str = Field("", max_length=2, description="PDB formal charge field")
    alt_loc: str = Field("", max_length=1, description="PDB alternate location code")
    record_name: Literal["ATOM", "HETATM"] = Field("ATOM", description="PDB record type")

    @classmethod
    def from_pdb_line(cls, line: str, *, atom_index: int | None = None) -> PdbAtomRecord:
        """Parse one fixed-width PDB ATOM/HETATM line."""
        return parse_pdb_atom_line(line, atom_index=atom_index)


PdbAtomRecord.__module__ = "polyzymd.builders.conjugation.structure.pdb"


class PdbLinkSide(BaseModel):
    """One residue side parsed from a PDB LINK record."""

    atom_name: str | None = None
    residue_name: str = ""
    chain_id: str = ""
    residue_number: int | None = None
    res_seq: str | None = None
    insertion_code: str | None = None


def parse_pdb_atom_line(line: str, *, atom_index: int | None = None) -> PdbAtomRecord:
    """Parse one fixed-width PDB ATOM/HETATM line into a canonical atom record."""
    if not is_pdb_atom_record(line):
        raise ValueError("PDB atom parsing requires an ATOM or HETATM record")

    atom_name = safe_pdb_slice(line, 12, 16).strip()
    residue_name = safe_pdb_slice(line, 17, 20).strip()
    residue_number = parse_int(safe_pdb_slice(line, 22, 26).strip())
    x = parse_float(safe_pdb_slice(line, 30, 38).strip())
    y = parse_float(safe_pdb_slice(line, 38, 46).strip())
    z = parse_float(safe_pdb_slice(line, 46, 54).strip())
    if not atom_name or not residue_name or residue_number is None:
        raise ValueError(f"Invalid PDB atom record fields: {line.rstrip()}")
    if x is None or y is None or z is None:
        raise ValueError(f"Invalid PDB coordinate fields: {line.rstrip()}")

    record_name = "HETATM" if line.startswith("HETATM") else "ATOM"
    occupancy = parse_float(safe_pdb_slice(line, 54, 60).strip())
    temp_factor = parse_float(safe_pdb_slice(line, 60, 66).strip())
    return PdbAtomRecord(
        serial=parse_int(safe_pdb_slice(line, 6, 11).strip()),
        atom_index=atom_index,
        atom_name=atom_name,
        residue_name=residue_name.upper(),
        chain_id=safe_pdb_slice(line, 21, 22).strip(),
        residue_number=residue_number,
        insertion_code=safe_pdb_slice(line, 26, 27).strip(),
        x=x,
        y=y,
        z=z,
        occupancy=1.0 if occupancy is None else occupancy,
        temp_factor=0.0 if temp_factor is None else temp_factor,
        element=parse_pdb_element(line),
        charge=safe_pdb_slice(line, 78, 80).strip(),
        alt_loc=safe_pdb_slice(line, 16, 17).strip(),
        record_name=record_name,
    )


def parse_pdb_atom_records(
    path: Path | str | None,
    *,
    require_atoms: bool = False,
) -> tuple[PdbAtomRecord, ...]:
    """Parse all PDB ATOM/HETATM records from a file in source order."""
    if path is None:
        if require_atoms:
            raise ValueError("PDB atom records require a path")
        return ()
    pdb_path = Path(path)
    atoms: list[PdbAtomRecord] = []
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if is_pdb_atom_record(line):
                atoms.append(parse_pdb_atom_line(line, atom_index=len(atoms)))
    if require_atoms and not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {pdb_path}")
    return tuple(atoms)


def parse_pdb_atom_lines(lines: Any) -> tuple[PdbAtomRecord, ...]:
    """Parse ATOM/HETATM records from an iterable of PDB lines in source order."""
    atoms: list[PdbAtomRecord] = []
    for line in lines:
        if is_pdb_atom_record(line):
            atoms.append(parse_pdb_atom_line(line, atom_index=len(atoms)))
    return tuple(atoms)


def pdb_coordinates(
    path: Path | str, *, require_atoms: bool = True
) -> tuple[tuple[float, float, float], ...]:
    """Parse PDB atom coordinates in source order."""
    coords = tuple((atom.x, atom.y, atom.z) for atom in parse_pdb_atom_records(path))
    if require_atoms and not coords:
        raise ValueError(f"No ATOM/HETATM coordinates found in {path}")
    return coords


def parse_pdb_conect_serials(line: str) -> tuple[int, ...]:
    """Parse source and target atom serials from one CONECT record."""
    if not line.startswith("CONECT"):
        return ()
    serials: list[int] = []
    for start in range(6, len(line), 5):
        serial = parse_int(line[start : start + 5].strip())
        if serial is not None:
            serials.append(serial)
    return tuple(serials)


def parse_pdb_conect_target_serials(line: str) -> tuple[int, ...]:
    """Parse target atom serials from one CONECT record."""
    serials = parse_pdb_conect_serials(line)
    return serials[1:] if len(serials) > 1 else ()


def parse_pdb_conect_pairs(path: Path | str | None) -> tuple[tuple[int, int], ...]:
    """Parse unique PDB CONECT bonds as normalized serial-number pairs."""
    if path is None:
        return ()
    pdb_path = Path(path)
    if not pdb_path.exists():
        return ()
    bonds: set[tuple[int, int]] = set()
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            serials = parse_pdb_conect_serials(line)
            if len(serials) < 2:
                continue
            source = serials[0]
            for target in serials[1:]:
                if target != source:
                    bonds.add(tuple(sorted((source, target))))
    return tuple(sorted(bonds))


def parse_pdb_conect_adjacency(path: Path | str | None) -> dict[int, tuple[int, ...]]:
    """Parse CONECT records as an undirected adjacency mapping keyed by serial."""
    adjacency: dict[int, set[int]] = {}
    for left, right in parse_pdb_conect_pairs(path):
        adjacency.setdefault(left, set()).add(right)
        adjacency.setdefault(right, set()).add(left)
    return {serial: tuple(sorted(neighbors)) for serial, neighbors in sorted(adjacency.items())}


def parse_pdb_conect_records(path: Path | str | None) -> dict[int, tuple[int, ...]]:
    """Parse CONECT records as directed source-to-target serial lists."""
    if path is None:
        return {}
    records: dict[int, list[int]] = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            serials = parse_pdb_conect_serials(line)
            if len(serials) < 2:
                continue
            records.setdefault(serials[0], []).extend(serials[1:])
    return {source: tuple(targets) for source, targets in records.items()}


def pdb_has_conect_records(path: Path | str) -> bool:
    """Return whether a PDB file contains one or more CONECT records."""
    with Path(path).open("r", encoding="utf-8") as handle:
        return any(line.startswith("CONECT") for line in handle)


def parse_pdb_link_sides(
    line: str, *, whitespace_fallback: bool = True
) -> tuple[PdbLinkSide, PdbLinkSide]:
    """Parse both residue sides from a PDB LINK record.

    The optional whitespace fallback preserves diagnostics for hand-written LINK
    records that do not strictly follow fixed-width PDB columns.
    """
    fixed_sides = _parse_fixed_width_link_sides(line)
    if not whitespace_fallback or _link_sides_have_required_identity(fixed_sides):
        return fixed_sides

    parts = line.split()
    if len(parts) >= 9 and parts[0] == "LINK":
        return _parse_whitespace_link_sides(parts)
    return fixed_sides


def _parse_fixed_width_link_sides(line: str) -> tuple[PdbLinkSide, PdbLinkSide]:
    """Parse LINK sides using the PDB fixed-width columns."""
    return (
        PdbLinkSide(
            atom_name=safe_pdb_slice(line, 12, 16).strip() or None,
            residue_name=(safe_pdb_slice(line, 17, 20).strip() or "").upper(),
            chain_id=safe_pdb_slice(line, 21, 22).strip(),
            residue_number=parse_int(safe_pdb_slice(line, 22, 26).strip()),
            res_seq=safe_pdb_slice(line, 22, 26).strip() or None,
            insertion_code=safe_pdb_slice(line, 26, 27).strip() or None,
        ),
        PdbLinkSide(
            atom_name=safe_pdb_slice(line, 42, 46).strip() or None,
            residue_name=(safe_pdb_slice(line, 47, 50).strip() or "").upper(),
            chain_id=safe_pdb_slice(line, 51, 52).strip(),
            residue_number=parse_int(safe_pdb_slice(line, 52, 56).strip()),
            res_seq=safe_pdb_slice(line, 52, 56).strip() or None,
            insertion_code=safe_pdb_slice(line, 56, 57).strip() or None,
        ),
    )


def _link_sides_have_required_identity(sides: tuple[PdbLinkSide, PdbLinkSide]) -> bool:
    """Return whether fixed-width LINK sides contain enough residue identity."""
    return all(
        side.atom_name is not None
        and bool(side.residue_name)
        and side.residue_number is not None
        and side.res_seq is not None
        for side in sides
    )


def _parse_whitespace_link_sides(parts: list[str]) -> tuple[PdbLinkSide, PdbLinkSide]:
    """Parse malformed hand-written LINK sides from whitespace tokens."""
    return (
        PdbLinkSide(
            atom_name=parts[1],
            residue_name=parts[2].upper(),
            chain_id=parts[3],
            residue_number=parse_int(parts[4]),
            res_seq=parts[4],
        ),
        PdbLinkSide(
            atom_name=parts[5],
            residue_name=parts[6].upper(),
            chain_id=parts[7],
            residue_number=parse_int(parts[8]),
            res_seq=parts[8],
        ),
    )


def pdb_link_side_dicts(
    line: str, *, whitespace_fallback: bool = True
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Parse LINK sides and return dictionaries for legacy diagnostic code."""
    side_1, side_2 = parse_pdb_link_sides(line, whitespace_fallback=whitespace_fallback)
    return side_1.model_dump(mode="python"), side_2.model_dump(mode="python")


def parse_pdb_atom_metadata_line(line: str, *, atom_index: int, line_number: int) -> dict[str, Any]:
    """Parse tolerant atom metadata for diagnostics and preflight paths."""
    residue_number = parse_int(safe_pdb_slice(line, 22, 26).strip())
    atom_serial = parse_int(safe_pdb_slice(line, 6, 11).strip())
    return {
        "atom_index": atom_index,
        "atom_serial": atom_serial,
        "atom_name": safe_pdb_slice(line, 12, 16).strip() or None,
        "residue_name": (safe_pdb_slice(line, 17, 20).strip() or "").upper(),
        "chain_id": safe_pdb_slice(line, 21, 22).strip(),
        "residue_number": residue_number,
        "res_seq": safe_pdb_slice(line, 22, 26).strip() or None,
        "residue_index": None,
        "insertion_code": safe_pdb_slice(line, 26, 27).strip() or None,
        "record_name": safe_pdb_slice(line, 0, 6).strip(),
        "line_number": line_number,
        "element": (
            safe_pdb_slice(line, 76, 78).strip() or _element_from_atom_name(line[12:16])
        ).upper()
        or None,
        "x": parse_float(safe_pdb_slice(line, 30, 38).strip()),
        "y": parse_float(safe_pdb_slice(line, 38, 46).strip()),
        "z": parse_float(safe_pdb_slice(line, 46, 54).strip()),
    }


def is_pdb_atom_record(line: str) -> bool:
    """Return whether a line is an ATOM or HETATM PDB record."""
    return line.startswith(ATOM_RECORD_PREFIXES)


def safe_pdb_slice(value: str, start: int, stop: int) -> str:
    """Return a safe fixed-width slice from a possibly short line."""
    if len(value) < stop:
        value = value.ljust(stop)
    return value[start:stop]


def parse_int(value: Any) -> int | None:
    """Parse an integer-like value, returning ``None`` for blanks or invalid text."""
    if value is None:
        return None
    if isinstance(value, int):
        return value
    try:
        text = str(value).strip()
        return int(text) if text else None
    except (TypeError, ValueError):
        return None


def parse_float(value: Any) -> float | None:
    """Parse a float-like value, returning ``None`` for blanks or invalid text."""
    if value is None:
        return None
    if isinstance(value, float):
        return value
    try:
        text = str(value).strip()
        return float(text) if text else None
    except (TypeError, ValueError):
        return None


def parse_pdb_element(line: str) -> str:
    """Parse or infer a PDB element symbol."""
    element = safe_pdb_slice(line, 76, 78).strip()
    if element:
        return format_pdb_element(element, safe_pdb_slice(line, 12, 16).strip())
    return format_pdb_element("", safe_pdb_slice(line, 12, 16).strip())


def format_pdb_element(element: str, atom_name: str) -> str:
    """Return a two-character PDB element field."""
    normalized = (element or "").strip()
    if normalized:
        return normalized[:1].upper() if len(normalized) == 1 else normalized[:2].title()
    inferred = "".join(ch for ch in atom_name if ch.isalpha())
    if not inferred:
        return ""
    if len(inferred) >= 2 and inferred[:2].title() in {"Cl", "Na", "Mg", "Ca", "Zn"}:
        return inferred[:2].title()
    return inferred[0].upper()


def _element_from_atom_name(atom_name: str) -> str:
    """Infer a one-character element from a PDB atom name field."""
    stripped = atom_name.strip()
    if not stripped:
        return ""
    if stripped[0].isdigit() and len(stripped) > 1:
        return stripped[1]
    return stripped[0]
