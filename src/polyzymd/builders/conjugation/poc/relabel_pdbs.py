#!/usr/bin/env python
"""Relabel POC PDB outputs to the PolyzyMD chain convention

This ad-hoc script rewrites existing PDB files in place and stores a .bak
backup before writing the corrected records
"""

from __future__ import annotations

import shutil
from pathlib import Path

_POC_DIR = Path(__file__).resolve().parent

_ASSEMBLED_REFERENCE_PDB = _POC_DIR / "output" / "conjugate_assembled.pdb"

_TARGET_PDBS = [
    _POC_DIR / "output" / "simulation_output" / "solvated_minimized.pdb",
    _POC_DIR / "output" / "simulation_output" / "equilibrated.pdb",
    _POC_DIR / "output" / "simulation_output" / "production_final.pdb",
    _POC_DIR / "output" / "conjugate_assembled.pdb",
]

_N_PROTEIN_ATOMS = 2717
_N_CONJUGATED_POLYMER_ATOMS = 718
_N_FREE_POLYMER_ATOMS = 945
_POLYMER_END = _N_PROTEIN_ATOMS + _N_CONJUGATED_POLYMER_ATOMS + _N_FREE_POLYMER_ATOMS

_CROSSLINK_RESIDS = {5, 15}


def _parse_element(line: str) -> str:
    """Extract element from PDB line with robust fallback

    Parameters
    ----------
    line : str
        ATOM/HETATM line from PDB

    Returns
    -------
    str
        Element symbol
    """

    elem = line[76:78].strip()
    if elem:
        if len(elem) == 1:
            return elem.upper()
        return elem[0].upper() + elem[1:].lower()

    atom_name = line[12:16].strip()
    if not atom_name:
        return ""

    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if not letters:
        return ""

    if len(letters) >= 2 and letters[:2].title() in {"Na", "Cl"}:
        return letters[:2].title()
    return letters[0].upper()


def _extract_protein_records(
    assembled_atoms: list[dict[str, str]],
) -> list[tuple[str, str, str]]:
    """Extract canonical protein atom naming from assembled reference PDB

    The assembled PDB already contains the correct 2717 protein atoms
    (with 4 HZ atoms removed during conjugation), so we use it directly
    instead of the crystal PDB which has 2721

    Parameters
    ----------
    assembled_atoms : list[dict[str, str]]
        Parsed atom records from the assembled reference PDB

    Returns
    -------
    list[tuple[str, str, str]]
        Ordered tuples of ``(atom_name, residue_name, residue_number)``
        for the first ``_N_PROTEIN_ATOMS`` entries
    """

    if len(assembled_atoms) < _N_PROTEIN_ATOMS:
        raise ValueError(
            f"Assembled PDB too short for protein extraction: "
            f"need {_N_PROTEIN_ATOMS}, got {len(assembled_atoms)}"
        )
    return [
        (rec["atom_name"], rec["residue_name"], rec["residue_number"])
        for rec in assembled_atoms[:_N_PROTEIN_ATOMS]
    ]


def _load_atom_records(path: Path) -> list[dict[str, str]]:
    """Load ATOM/HETATM records from a PDB file

    Parameters
    ----------
    path : Path
        PDB file path

    Returns
    -------
    list[dict[str, str]]
        Parsed ATOM/HETATM records in file order
    """

    atom_records: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_records.append(
                {
                    "record_name": line[0:6],
                    "atom_name": line[12:16].strip(),
                    "serial": line[6:11],
                    "residue_name": line[17:20].strip(),
                    "residue_number": line[22:26].strip(),
                    "x": line[30:38],
                    "y": line[38:46],
                    "z": line[46:54],
                    "occ": line[54:60],
                    "temp": line[60:66],
                    "element": _parse_element(line),
                    "charge": line[78:80] if len(line) >= 80 else "  ",
                }
            )
    return atom_records


def _build_polymer_reference(
    assembled_atoms: list[dict[str, str]],
) -> tuple[dict[int, tuple[str, int]], dict[int, int]]:
    """Build polymer residue mappings from assembled reference PDB

    Parameters
    ----------
    assembled_atoms : list[dict[str, str]]
        Parsed atom records from assembled PDB

    Returns
    -------
    tuple[dict[int, tuple[str, int]], dict[int, int]]
        Mapping of atom index to ``(resname, old_resid)`` and atom index to
        new sequential polymer residue number
    """

    if len(assembled_atoms) < _POLYMER_END:
        raise ValueError(
            "Assembled PDB is too short for expected polymer reference atoms: "
            f"need {_POLYMER_END}, got {len(assembled_atoms)}"
        )

    atom_to_residue: dict[int, tuple[str, int]] = {}
    atom_to_new_resid: dict[int, int] = {}
    prev_old_resid: int | None = None
    new_resid = 0

    for atom_index in range(_N_PROTEIN_ATOMS, _POLYMER_END):
        record = assembled_atoms[atom_index]
        old_resid = int(record["residue_number"])
        old_resname = record["residue_name"]
        atom_to_residue[atom_index] = (old_resname, old_resid)

        if old_resid != prev_old_resid:
            new_resid += 1
            prev_old_resid = old_resid
        atom_to_new_resid[atom_index] = new_resid

    return atom_to_residue, atom_to_new_resid


def _canonicalize_conjugated_resname(resname: str, assembled_resid: int) -> str:
    """Map connectivity names to canonical conjugated polymer names

    Parameters
    ----------
    resname : str
        Original residue name from assembled PDB reference
    assembled_resid : int
        Original assembled-PDB residue number

    Returns
    -------
    str
        Canonical residue name
    """

    if resname in {"SB1", "SB2"}:
        return "SBM"
    if resname in {"EG1", "EG2"}:
        return "EGP"
    if resname in {"NH1", "NH2"}:
        return "NHX" if assembled_resid in _CROSSLINK_RESIDS else "NHS"
    return resname


def _chain_from_solvent_resid(resid: int) -> tuple[str, int]:
    """Convert global solvent residue count to chain and local residue number

    Parameters
    ----------
    resid : int
        One-based global solvent residue counter

    Returns
    -------
    tuple[str, int]
        Chain ID and residue number within that chain
    """

    chain_offset = (resid - 1) // 9999
    chain_code = ord("D") + chain_offset
    if chain_code > ord("Z"):
        raise ValueError("Solvent chain overflow beyond chain Z")
    local_resid = ((resid - 1) % 9999) + 1
    return chr(chain_code), local_resid


def _format_pdb_line(
    record_name: str,
    serial: str,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: str,
    y: str,
    z: str,
    occ: str,
    temp: str,
    element: str,
    charge: str,
) -> str:
    """Format a PDB ATOM/HETATM line preserving coordinates and serial

    Parameters
    ----------
    record_name : str
        Record field including spacing (ATOM/HETATM)
    serial : str
        Atom serial field
    atom_name : str
        Atom name
    residue_name : str
        Residue name
    chain_id : str
        Chain identifier
    residue_number : int
        Residue number
    x : str
        X coordinate field
    y : str
        Y coordinate field
    z : str
        Z coordinate field
    occ : str
        Occupancy field
    temp : str
        Temperature factor field
    element : str
        Element symbol
    charge : str
        Charge field

    Returns
    -------
    str
        Formatted PDB record line
    """

    if len(atom_name) < 4:
        atom_field = f" {atom_name:<3}"
    else:
        atom_field = f"{atom_name:<4}"

    return (
        f"{record_name:<6}{serial:>5} {atom_field} {residue_name:>3} {chain_id:1}"
        f"{residue_number:>4d}    {x:>8}{y:>8}{z:>8}{occ:>6}{temp:>6}          "
        f"{element:>2}{charge:>2}\n"
    )


def _rewrite_pdb(
    path: Path,
    protein_records: list[tuple[str, str, str]],
    polymer_ref: dict[int, tuple[str, int]],
    polymer_new_resid: dict[int, int],
) -> dict[str, int]:
    """Rewrite one PDB file in place to the PolyzyMD convention

    Parameters
    ----------
    path : Path
        Target PDB path to rewrite
    protein_records : list[tuple[str, str, str]]
        Canonical protein atom records from assembled reference PDB
    polymer_ref : dict[int, tuple[str, int]]
        Atom-index mapping to assembled polymer ``(resname, resid)``
    polymer_new_resid : dict[int, int]
        Atom-index mapping to sequential polymer residue numbers

    Returns
    -------
    dict[str, int]
        Summary counters for rewritten components
    """

    backup_path = path.with_suffix(path.suffix + ".bak")
    shutil.copy2(path, backup_path)

    with path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()

    rewritten: list[str] = []
    atom_index = 0
    solvent_resid_counter = 0

    n_protein = 0
    n_conjugated = 0
    n_free = 0
    n_water = 0
    n_na = 0
    n_cl = 0

    while_index = 0
    while while_index < len(lines):
        line = lines[while_index]
        if not line.startswith(("ATOM", "HETATM")):
            # Strip stale TER/END records — we regenerate them at the end
            # Keep other metadata lines (CRYST1, REMARK, etc.) as-is
            if not line.startswith(("TER", "END")):
                rewritten.append(line)
            while_index += 1
            continue

        if atom_index < _N_PROTEIN_ATOMS:
            atom_name, residue_name, residue_number = protein_records[atom_index]
            rec = {
                "record_name": line[0:6],
                "serial": line[6:11],
                "x": line[30:38],
                "y": line[38:46],
                "z": line[46:54],
                "occ": line[54:60],
                "temp": line[60:66],
                "element": _parse_element(line),
                "charge": line[78:80] if len(line) >= 80 else "  ",
            }
            rewritten.append(
                _format_pdb_line(
                    record_name=rec["record_name"],
                    serial=rec["serial"],
                    atom_name=atom_name,
                    residue_name=residue_name,
                    chain_id="A",
                    residue_number=int(residue_number),
                    x=rec["x"],
                    y=rec["y"],
                    z=rec["z"],
                    occ=rec["occ"],
                    temp=rec["temp"],
                    element=rec["element"],
                    charge=rec["charge"],
                )
            )
            atom_index += 1
            n_protein += 1
            while_index += 1
            continue

        if atom_index < _N_PROTEIN_ATOMS + _N_CONJUGATED_POLYMER_ATOMS:
            old_resname, old_resid = polymer_ref[atom_index]
            rec = {
                "record_name": line[0:6],
                "serial": line[6:11],
                "x": line[30:38],
                "y": line[38:46],
                "z": line[46:54],
                "occ": line[54:60],
                "temp": line[60:66],
                "element": _parse_element(line),
                "charge": line[78:80] if len(line) >= 80 else "  ",
            }
            rewritten.append(
                _format_pdb_line(
                    record_name=rec["record_name"],
                    serial=rec["serial"],
                    atom_name=line[12:16].strip(),
                    residue_name=_canonicalize_conjugated_resname(old_resname, old_resid),
                    chain_id="C",
                    residue_number=polymer_new_resid[atom_index],
                    x=rec["x"],
                    y=rec["y"],
                    z=rec["z"],
                    occ=rec["occ"],
                    temp=rec["temp"],
                    element=rec["element"],
                    charge=rec["charge"],
                )
            )
            atom_index += 1
            n_conjugated += 1
            while_index += 1
            continue

        if atom_index < _POLYMER_END:
            old_resname, _old_resid = polymer_ref[atom_index]
            rec = {
                "record_name": line[0:6],
                "serial": line[6:11],
                "x": line[30:38],
                "y": line[38:46],
                "z": line[46:54],
                "occ": line[54:60],
                "temp": line[60:66],
                "element": _parse_element(line),
                "charge": line[78:80] if len(line) >= 80 else "  ",
            }
            rewritten.append(
                _format_pdb_line(
                    record_name=rec["record_name"],
                    serial=rec["serial"],
                    atom_name=line[12:16].strip(),
                    residue_name=old_resname,
                    chain_id="C",
                    residue_number=polymer_new_resid[atom_index],
                    x=rec["x"],
                    y=rec["y"],
                    z=rec["z"],
                    occ=rec["occ"],
                    temp=rec["temp"],
                    element=rec["element"],
                    charge=rec["charge"],
                )
            )
            atom_index += 1
            n_free += 1
            while_index += 1
            continue

        solvent_records = []
        while while_index < len(lines):
            candidate = lines[while_index]
            if not candidate.startswith(("ATOM", "HETATM")):
                break
            solvent_records.append(candidate)
            break

        if not solvent_records:
            rewritten.append(line)
            while_index += 1
            continue

        first_line = solvent_records[0]
        element_0 = _parse_element(first_line)

        if element_0 == "O" and while_index + 2 < len(lines):
            next_1 = lines[while_index + 1]
            next_2 = lines[while_index + 2]
            if next_1.startswith(("ATOM", "HETATM")) and next_2.startswith(("ATOM", "HETATM")):
                e1 = _parse_element(next_1)
                e2 = _parse_element(next_2)
                if e1 == "H" and e2 == "H":
                    solvent_resid_counter += 1
                    chain_id, residue_number = _chain_from_solvent_resid(solvent_resid_counter)

                    for offset, atom_name in enumerate(["O", "H1", "H2"]):
                        src_line = lines[while_index + offset]
                        rewritten.append(
                            _format_pdb_line(
                                record_name=src_line[0:6],
                                serial=src_line[6:11],
                                atom_name=atom_name,
                                residue_name="HOH",
                                chain_id=chain_id,
                                residue_number=residue_number,
                                x=src_line[30:38],
                                y=src_line[38:46],
                                z=src_line[46:54],
                                occ=src_line[54:60],
                                temp=src_line[60:66],
                                element=_parse_element(src_line),
                                charge=src_line[78:80] if len(src_line) >= 80 else "  ",
                            )
                        )
                        atom_index += 1
                    n_water += 1
                    while_index += 3
                    continue

        if element_0 in {"Na", "Cl"}:
            solvent_resid_counter += 1
            chain_id, residue_number = _chain_from_solvent_resid(solvent_resid_counter)
            residue_name = "NA" if element_0 == "Na" else "CL"
            atom_name = residue_name
            rewritten.append(
                _format_pdb_line(
                    record_name=first_line[0:6],
                    serial=first_line[6:11],
                    atom_name=atom_name,
                    residue_name=residue_name,
                    chain_id=chain_id,
                    residue_number=residue_number,
                    x=first_line[30:38],
                    y=first_line[38:46],
                    z=first_line[46:54],
                    occ=first_line[54:60],
                    temp=first_line[60:66],
                    element=element_0,
                    charge=first_line[78:80] if len(first_line) >= 80 else "  ",
                )
            )
            atom_index += 1
            if residue_name == "NA":
                n_na += 1
            else:
                n_cl += 1
            while_index += 1
            continue

        # Fallback: unknown solvent species treated as single-atom residue.
        # This is correct for this system (only water, Na+, Cl-) but would
        # need grouping logic for multi-atom co-solvents.
        solvent_resid_counter += 1
        chain_id, residue_number = _chain_from_solvent_resid(solvent_resid_counter)
        rewritten.append(
            _format_pdb_line(
                record_name=first_line[0:6],
                serial=first_line[6:11],
                atom_name=first_line[12:16].strip(),
                residue_name=first_line[17:20].strip() or "UNK",
                chain_id=chain_id,
                residue_number=residue_number,
                x=first_line[30:38],
                y=first_line[38:46],
                z=first_line[46:54],
                occ=first_line[54:60],
                temp=first_line[60:66],
                element=element_0,
                charge=first_line[78:80] if len(first_line) >= 80 else "  ",
            )
        )
        atom_index += 1
        while_index += 1

    # Insert TER records at chain boundaries and trailing END
    final_lines: list[str] = []
    prev_chain: str | None = None
    for out_line in rewritten:
        if out_line.startswith(("ATOM", "HETATM")):
            current_chain = out_line[21]
            if prev_chain is not None and current_chain != prev_chain:
                final_lines.append("TER\n")
            prev_chain = current_chain
        final_lines.append(out_line)
    final_lines.append("TER\n")
    final_lines.append("END\n")

    with path.open("w", encoding="utf-8") as handle:
        handle.writelines(final_lines)

    return {
        "n_protein_atoms": n_protein,
        "n_conjugated_polymer_atoms": n_conjugated,
        "n_free_polymer_atoms": n_free,
        "n_water": n_water,
        "n_na": n_na,
        "n_cl": n_cl,
    }


def main() -> int:
    """Run ad-hoc relabeling for all POC target PDB files

    Returns
    -------
    int
        Process exit code
    """

    assembled_atoms = _load_atom_records(_ASSEMBLED_REFERENCE_PDB)
    protein_records = _extract_protein_records(assembled_atoms)
    polymer_ref, polymer_new_resid = _build_polymer_reference(assembled_atoms)

    print("Relabeling PDB files to PolyzyMD chain convention")
    for pdb_path in _TARGET_PDBS:
        if not pdb_path.exists():
            print(f"  Skipping missing file: {pdb_path}")
            continue

        summary = _rewrite_pdb(
            path=pdb_path,
            protein_records=protein_records,
            polymer_ref=polymer_ref,
            polymer_new_resid=polymer_new_resid,
        )
        print(f"  Rewrote {pdb_path}")
        print(f"    Backup: {pdb_path.with_suffix(pdb_path.suffix + '.bak')}")
        print(f"    Summary: {summary}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
