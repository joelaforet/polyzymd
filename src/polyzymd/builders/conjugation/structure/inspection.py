"""Pure-Python PDB compatibility inspection for Pablo ingestion preflight."""

from __future__ import annotations

from collections import Counter, defaultdict
from math import dist
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.structure.parsing import (
    parse_int,
    parse_pdb_atom_metadata_line,
    parse_pdb_conect_target_serials,
    pdb_link_side_dicts,
)

POLYZMD_PROTEIN_CHAIN = "A"
POLYZMD_MOIETY_CHAIN = "C"
STRUCTURAL_ATTACHMENT_CUTOFF_ANGSTROM = 1.9

STANDARD_PROTEIN_RESIDUES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
)
MODIFIED_PROTEIN_RESIDUES = frozenset({"LYX"})
WATER_RESIDUES = frozenset({"HOH", "WAT", "H2O", "SOL", "TIP", "TIP3", "TIP3P"})
COMMON_ION_RESIDUES = frozenset(
    {
        "NA",
        "K",
        "CL",
        "CA",
        "MG",
        "ZN",
        "FE",
        "MN",
        "CU",
        "CO",
        "NI",
        "BR",
        "IOD",
        "LI",
        "RB",
        "CS",
        "F",
        "SR",
        "BA",
    }
)
COMMON_SOLVENT_RESIDUES = frozenset(
    {
        "ACE",
        "ACT",
        "ACN",
        "DMS",
        "DMSO",
        "EDO",
        "EOH",
        "GOL",
        "IPA",
        "MOH",
        "MPD",
        "SO4",
    }
)
GLYCAN_RESIDUES = frozenset(
    {
        "A2G",
        "BMA",
        "BGC",
        "FUC",
        "GAL",
        "GLC",
        "MAN",
        "NAG",
        "NDG",
        "SIA",
    }
)
POLYMER_LIKE_RESIDUES = frozenset(
    {
        "PEG",
        "PEO",
        "PGE",
        "PG4",
        "P6G",
        "PME",
        "PLL",
    }
)


class _DiagnosticAtomRecord(BaseModel):
    """Lightweight atom metadata parsed from a PDB ATOM/HETATM record."""

    atom_index: int
    atom_serial: int | None = None
    atom_name: str | None = None
    residue_name: str
    chain_id: str = ""
    residue_number: int | None = None
    res_seq: str | None = None
    insertion_code: str | None = None
    record_name: str
    line_number: int
    element: str | None = None
    x: float | None = None
    y: float | None = None
    z: float | None = None


class PDBResidueInspection(BaseModel):
    """Residue-level classification used for diagnostics and compatibility warnings."""

    residue_id: str
    chain_id: str = ""
    residue_name: str
    residue_number: int | None = None
    res_seq: str | None = None
    insertion_code: str | None = None
    atom_count: int
    record_names: list[str] = Field(default_factory=list)
    category: str
    is_standard_protein: bool = False
    is_noncanonical: bool = False
    is_polymer_ptm_candidate: bool = False
    nearest_protein_distance_angstrom: float | None = None


class PDBCovalentAttachmentCandidate(BaseModel):
    """Potential protein-to-noncanonical covalent attachment evidence."""

    source: str
    protein_residue_id: str
    candidate_residue_id: str
    protein_chain_id: str = ""
    candidate_chain_id: str = ""
    protein_residue_name: str
    candidate_residue_name: str
    candidate_category: str
    protein_atom_name: str | None = None
    candidate_atom_name: str | None = None
    protein_atom_serial: int | None = None
    candidate_atom_serial: int | None = None
    line_number: int | None = None
    distance_angstrom: float | None = None
    details: dict[str, Any] = Field(default_factory=dict)


class PDBStructureInspection(BaseModel):
    """Summary of PDB/Pablo compatibility diagnostics.

    The inspection is diagnostics-only. It does not rewrite residue names, chain IDs,
    connectivity records, or coordinates.
    """

    path: Path
    atom_count: int = 0
    residue_count: int = 0
    chain_ids: list[str] = Field(default_factory=list)
    blank_chain_atom_count: int = 0
    blank_chain_residue_count: int = 0
    residue_name_counts: dict[str, int] = Field(default_factory=dict)
    protein_like_canonical_residue_count: int = 0
    water_ion_solvent_like_residue_count: int = 0
    ligand_cocrystal_like_residue_count: int = 0
    noncanonical_residue_candidates: list[PDBResidueInspection] = Field(default_factory=list)
    polymer_ptm_candidates: list[PDBResidueInspection] = Field(default_factory=list)
    covalent_attachment_candidates: list[PDBCovalentAttachmentCandidate] = Field(
        default_factory=list
    )
    ssbond_count: int = 0
    convention_warnings: list[str] = Field(default_factory=list)
    compatibility_warnings: list[str] = Field(default_factory=list)


def inspect_pdb_structure(path: Path | str) -> PDBStructureInspection:
    """Inspect a PDB file for diagnostics-only Pablo compatibility signals.

    Parameters
    ----------
    path : Path or str
        PDB file to inspect.

    Returns
    -------
    PDBStructureInspection
        JSON-safe structure summary, residue classifications, attachment evidence,
        and warning messages.
    """
    pdb_path = Path(path)
    lines = pdb_path.read_text(errors="replace").splitlines()
    atoms = _parse_atom_records(lines)
    residues = _summarize_residues(atoms)
    residue_lookup = {residue.residue_id: residue for residue in residues}
    serial_lookup = {atom.atom_serial: atom for atom in atoms if atom.atom_serial is not None}
    residue_atoms = _group_atoms_by_residue(atoms)
    _annotate_nearest_protein_distances(residues, residue_atoms)

    covalent_candidates = _parse_link_candidates(lines, residue_lookup)
    covalent_candidates.extend(_parse_conect_candidates(lines, serial_lookup))
    covalent_candidates = _deduplicate_candidates(covalent_candidates)

    residue_name_counts = Counter(atom.residue_name for atom in atoms)
    chain_ids = sorted({atom.chain_id for atom in atoms if atom.chain_id})
    blank_residues = [residue for residue in residues if not residue.chain_id]
    noncanonical = [residue for residue in residues if residue.is_noncanonical]
    polymer_ptm = [residue for residue in noncanonical if residue.is_polymer_ptm_candidate]
    convention_warnings = _build_convention_warnings(residues, covalent_candidates)
    compatibility_warnings = _build_compatibility_warnings(residues, covalent_candidates)

    return PDBStructureInspection(
        path=pdb_path,
        atom_count=len(atoms),
        residue_count=len(residues),
        chain_ids=chain_ids,
        blank_chain_atom_count=sum(1 for atom in atoms if not atom.chain_id),
        blank_chain_residue_count=len(blank_residues),
        residue_name_counts=dict(sorted(residue_name_counts.items())),
        protein_like_canonical_residue_count=sum(
            1 for residue in residues if residue.is_standard_protein
        ),
        water_ion_solvent_like_residue_count=sum(
            1 for residue in residues if residue.category in {"water", "ion", "solvent"}
        ),
        ligand_cocrystal_like_residue_count=sum(
            1 for residue in residues if residue.category == "ligand_cocrystal"
        ),
        noncanonical_residue_candidates=noncanonical,
        polymer_ptm_candidates=polymer_ptm,
        covalent_attachment_candidates=covalent_candidates,
        ssbond_count=sum(1 for line in lines if line.startswith("SSBOND")),
        convention_warnings=convention_warnings,
        compatibility_warnings=compatibility_warnings,
    )


def pdb_atom_records_as_dicts(inspection: PDBStructureInspection) -> list[dict[str, Any]]:
    """Return PDB atom records in Pablo adapter metadata format.

    Parameters
    ----------
    inspection : PDBStructureInspection
        Structure inspection created by :func:`inspect_pdb_structure`.

    Returns
    -------
    list of dict
        Atom dictionaries compatible with the existing Pablo adapter metadata path.
    """
    atoms = _parse_atom_records(inspection.path.read_text(errors="replace").splitlines())
    return [atom.model_dump(mode="python") for atom in atoms]


def _parse_atom_records(lines: list[str]) -> list[_DiagnosticAtomRecord]:
    """Parse fixed-width ATOM/HETATM records from PDB lines."""
    atoms: list[_DiagnosticAtomRecord] = []
    for line_number, line in enumerate(lines, start=1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms.append(
            _DiagnosticAtomRecord(
                **parse_pdb_atom_metadata_line(line, atom_index=len(atoms), line_number=line_number)
            )
        )
    return atoms


def _summarize_residues(atoms: list[_DiagnosticAtomRecord]) -> list[PDBResidueInspection]:
    """Group atom records into residue diagnostics."""
    grouped = _group_atoms_by_residue(atoms)
    residues: list[PDBResidueInspection] = []
    for key, residue_atoms in grouped.items():
        chain_id, residue_name, residue_number, res_seq, insertion_code = key
        category = _classify_residue(residue_name)
        residue_id = _format_residue_id(
            chain_id, residue_name, residue_number, res_seq, insertion_code
        )
        residues.append(
            PDBResidueInspection(
                residue_id=residue_id,
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=residue_number,
                res_seq=res_seq,
                insertion_code=insertion_code,
                atom_count=len(residue_atoms),
                record_names=sorted({atom.record_name for atom in residue_atoms}),
                category=category,
                is_standard_protein=_is_protein_like_residue(residue_name),
                is_noncanonical=_is_noncanonical_residue(residue_name),
                is_polymer_ptm_candidate=category in {"glycan", "polymer_ptm"},
            )
        )
    return residues


def _group_atoms_by_residue(
    atoms: list[_DiagnosticAtomRecord],
) -> dict[tuple[str, str, int | None, str | None, str | None], list[_DiagnosticAtomRecord]]:
    """Group atoms by PDB residue identity."""
    grouped: dict[
        tuple[str, str, int | None, str | None, str | None], list[_DiagnosticAtomRecord]
    ] = defaultdict(list)
    for atom in atoms:
        grouped[
            (
                atom.chain_id,
                atom.residue_name,
                atom.residue_number,
                atom.res_seq,
                atom.insertion_code,
            )
        ].append(atom)
    return dict(grouped)


def _annotate_nearest_protein_distances(
    residues: list[PDBResidueInspection],
    residue_atoms: dict[
        tuple[str, str, int | None, str | None, str | None], list[_DiagnosticAtomRecord]
    ],
) -> None:
    """Annotate noncanonical residues with nearest heavy-atom protein distances."""
    protein_atoms = [
        atom
        for atoms in residue_atoms.values()
        for atom in atoms
        if _is_protein_like_residue(atom.residue_name) and _has_heavy_atom_coordinates(atom)
    ]
    if not protein_atoms:
        return

    for residue in residues:
        if not residue.is_noncanonical:
            continue
        atoms = residue_atoms.get(
            (
                residue.chain_id,
                residue.residue_name,
                residue.residue_number,
                residue.res_seq,
                residue.insertion_code,
            ),
            [],
        )
        candidate_atoms = [atom for atom in atoms if _has_heavy_atom_coordinates(atom)]
        if not candidate_atoms:
            continue
        residue.nearest_protein_distance_angstrom = min(
            dist((atom.x, atom.y, atom.z), (protein_atom.x, protein_atom.y, protein_atom.z))
            for atom in candidate_atoms
            for protein_atom in protein_atoms
        )


def _parse_link_candidates(
    lines: list[str],
    residue_lookup: dict[str, PDBResidueInspection],
) -> list[PDBCovalentAttachmentCandidate]:
    """Parse protein-to-noncanonical candidates from LINK records."""
    candidates: list[PDBCovalentAttachmentCandidate] = []
    for line_number, line in enumerate(lines, start=1):
        if not line.startswith("LINK"):
            continue
        side_1, side_2 = _parse_link_sides(line)
        candidate = _candidate_from_link_sides(
            source="LINK",
            side_1=side_1,
            side_2=side_2,
            line_number=line_number,
            residue_lookup=residue_lookup,
            details={"raw_record": line.rstrip()},
        )
        if candidate is not None:
            candidates.append(candidate)
    return candidates


def _parse_link_sides(line: str) -> tuple[dict[str, Any], dict[str, Any]]:
    """Parse LINK residue sides with a whitespace fallback for hand-written records."""
    return pdb_link_side_dicts(line)


def _parse_conect_candidates(
    lines: list[str],
    serial_lookup: dict[int, _DiagnosticAtomRecord],
) -> list[PDBCovalentAttachmentCandidate]:
    """Parse protein-to-noncanonical candidates from CONECT records."""
    candidates: list[PDBCovalentAttachmentCandidate] = []
    seen_pairs: set[tuple[int, int]] = set()
    for line_number, line in enumerate(lines, start=1):
        if not line.startswith("CONECT"):
            continue
        source_serial = _parse_int(line[6:11].strip())
        if source_serial is None:
            continue
        for target_serial in parse_pdb_conect_target_serials(line):
            pair = tuple(sorted((source_serial, target_serial)))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            atom_1 = serial_lookup.get(source_serial)
            atom_2 = serial_lookup.get(target_serial)
            if atom_1 is None or atom_2 is None:
                continue
            candidate = _candidate_from_atoms(atom_1, atom_2, line_number)
            if candidate is not None:
                candidates.append(candidate)
    return candidates


def _candidate_from_link_sides(
    *,
    source: str,
    side_1: dict[str, Any],
    side_2: dict[str, Any],
    line_number: int,
    residue_lookup: dict[str, PDBResidueInspection],
    details: dict[str, Any] | None = None,
) -> PDBCovalentAttachmentCandidate | None:
    """Build a covalent candidate when one LINK side is protein and the other is noncanonical."""
    residue_1 = _residue_from_link_side(side_1, residue_lookup)
    residue_2 = _residue_from_link_side(side_2, residue_lookup)
    if residue_1 is None or residue_2 is None:
        return None
    return _candidate_from_residues(
        source=source,
        residue_1=residue_1,
        residue_2=residue_2,
        atom_name_1=side_1["atom_name"],
        atom_name_2=side_2["atom_name"],
        atom_serial_1=None,
        atom_serial_2=None,
        line_number=line_number,
        details=details,
    )


def _candidate_from_atoms(
    atom_1: _DiagnosticAtomRecord,
    atom_2: _DiagnosticAtomRecord,
    line_number: int,
) -> PDBCovalentAttachmentCandidate | None:
    """Build a covalent candidate from two atom records when roles match."""
    residue_1 = _residue_from_atom(atom_1)
    residue_2 = _residue_from_atom(atom_2)
    distance_angstrom = None
    if _has_coordinates(atom_1) and _has_coordinates(atom_2):
        distance_angstrom = dist((atom_1.x, atom_1.y, atom_1.z), (atom_2.x, atom_2.y, atom_2.z))
    return _candidate_from_residues(
        source="CONECT",
        residue_1=residue_1,
        residue_2=residue_2,
        atom_name_1=atom_1.atom_name,
        atom_name_2=atom_2.atom_name,
        atom_serial_1=atom_1.atom_serial,
        atom_serial_2=atom_2.atom_serial,
        line_number=line_number,
        distance_angstrom=distance_angstrom,
        details={"atom_serials": [atom_1.atom_serial, atom_2.atom_serial]},
    )


def _candidate_from_residues(
    *,
    source: str,
    residue_1: PDBResidueInspection,
    residue_2: PDBResidueInspection,
    atom_name_1: str | None,
    atom_name_2: str | None,
    atom_serial_1: int | None,
    atom_serial_2: int | None,
    line_number: int,
    distance_angstrom: float | None = None,
    details: dict[str, Any] | None = None,
) -> PDBCovalentAttachmentCandidate | None:
    """Orient a candidate as protein side then noncanonical side."""
    if residue_1.is_standard_protein and _is_attachment_candidate_residue(residue_2):
        protein = residue_1
        candidate = residue_2
        protein_atom_name = atom_name_1
        candidate_atom_name = atom_name_2
        protein_atom_serial = atom_serial_1
        candidate_atom_serial = atom_serial_2
    elif residue_2.is_standard_protein and _is_attachment_candidate_residue(residue_1):
        protein = residue_2
        candidate = residue_1
        protein_atom_name = atom_name_2
        candidate_atom_name = atom_name_1
        protein_atom_serial = atom_serial_2
        candidate_atom_serial = atom_serial_1
    else:
        return None

    return PDBCovalentAttachmentCandidate(
        source=source,
        protein_residue_id=protein.residue_id,
        candidate_residue_id=candidate.residue_id,
        protein_chain_id=protein.chain_id,
        candidate_chain_id=candidate.chain_id,
        protein_residue_name=protein.residue_name,
        candidate_residue_name=candidate.residue_name,
        candidate_category=candidate.category,
        protein_atom_name=protein_atom_name,
        candidate_atom_name=candidate_atom_name,
        protein_atom_serial=protein_atom_serial,
        candidate_atom_serial=candidate_atom_serial,
        line_number=line_number,
        distance_angstrom=distance_angstrom,
        details=details or {},
    )


def _deduplicate_candidates(
    candidates: list[PDBCovalentAttachmentCandidate],
) -> list[PDBCovalentAttachmentCandidate]:
    """Deduplicate connectivity records that describe the same residue attachment."""
    seen: set[tuple[str, str, str, str | None, str | None]] = set()
    unique: list[PDBCovalentAttachmentCandidate] = []
    for candidate in candidates:
        key = (
            candidate.source,
            candidate.protein_residue_id,
            candidate.candidate_residue_id,
            candidate.protein_atom_name,
            candidate.candidate_atom_name,
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(candidate)
    return unique


def _residue_from_link_side(
    side: dict[str, Any],
    residue_lookup: dict[str, PDBResidueInspection],
) -> PDBResidueInspection | None:
    """Resolve a LINK side into an inspected residue."""
    residue_id = _format_residue_id(
        side["chain_id"],
        side["residue_name"],
        side["residue_number"],
        side["res_seq"],
        side["insertion_code"],
    )
    residue = residue_lookup.get(residue_id)
    if residue is not None:
        return residue
    category = _classify_residue(side["residue_name"])
    if not side["residue_name"]:
        return None
    return PDBResidueInspection(
        residue_id=residue_id,
        chain_id=side["chain_id"],
        residue_name=side["residue_name"],
        residue_number=side["residue_number"],
        res_seq=side["res_seq"],
        insertion_code=side["insertion_code"],
        atom_count=0,
        category=category,
        is_standard_protein=_is_protein_like_residue(side["residue_name"]),
        is_noncanonical=_is_noncanonical_residue(side["residue_name"]),
        is_polymer_ptm_candidate=category in {"glycan", "polymer_ptm"},
    )


def _residue_from_atom(atom: _DiagnosticAtomRecord) -> PDBResidueInspection:
    """Create a residue classification from one atom record."""
    category = _classify_residue(atom.residue_name)
    return PDBResidueInspection(
        residue_id=_format_residue_id(
            atom.chain_id,
            atom.residue_name,
            atom.residue_number,
            atom.res_seq,
            atom.insertion_code,
        ),
        chain_id=atom.chain_id,
        residue_name=atom.residue_name,
        residue_number=atom.residue_number,
        res_seq=atom.res_seq,
        insertion_code=atom.insertion_code,
        atom_count=1,
        record_names=[atom.record_name],
        category=category,
        is_standard_protein=_is_protein_like_residue(atom.residue_name),
        is_noncanonical=_is_noncanonical_residue(atom.residue_name),
        is_polymer_ptm_candidate=category in {"glycan", "polymer_ptm"},
    )


def _build_convention_warnings(
    residues: list[PDBResidueInspection],
    covalent_candidates: list[PDBCovalentAttachmentCandidate],
) -> list[str]:
    """Build warning-only PolyzyMD chain convention messages."""
    warnings: list[str] = []
    protein_off_chain = sorted(
        {
            residue.chain_id or "blank"
            for residue in residues
            if residue.is_standard_protein and residue.chain_id != POLYZMD_PROTEIN_CHAIN
        }
    )
    if protein_off_chain:
        warnings.append(
            "PolyzyMD convention expects protein residues on chain A; found canonical protein "
            f"residues on chain(s): {', '.join(protein_off_chain)}"
        )

    blank_atom_count = sum(residue.atom_count for residue in residues if not residue.chain_id)
    blank_residue_count = sum(1 for residue in residues if not residue.chain_id)
    if blank_residue_count:
        warnings.append(
            "Blank chain IDs were found in the PDB; PolyzyMD convention expects explicit chain "
            f"IDs, including chain A for protein and chain C for attached moieties "
            f"({blank_residue_count} residues, {blank_atom_count} atoms)"
        )

    off_chain_candidates = sorted(
        {
            candidate.candidate_residue_id
            for candidate in covalent_candidates
            if candidate.candidate_chain_id != POLYZMD_MOIETY_CHAIN
        }
    )
    if off_chain_candidates:
        warnings.append(
            "PolyzyMD convention expects covalently attached PTM/glycan/polymer moieties on "
            f"chain C; linked candidates outside chain C: {', '.join(off_chain_candidates)}"
        )
    return warnings


def _build_compatibility_warnings(
    residues: list[PDBResidueInspection],
    covalent_candidates: list[PDBCovalentAttachmentCandidate],
) -> list[str]:
    """Build Pablo/CCD naming and connectivity compatibility warnings."""
    warnings: list[str] = []
    if residues:
        warnings.append(
            "Residue names should follow PDB CCD conventions where possible; diagnostics are "
            "warning-only and do not normalize names or chain IDs"
        )

    linked_ids = {candidate.candidate_residue_id for candidate in covalent_candidates}
    likely_missing_evidence = [
        residue
        for residue in residues
        if residue.is_polymer_ptm_candidate
        and residue.residue_id not in linked_ids
        and residue.nearest_protein_distance_angstrom is not None
        and residue.nearest_protein_distance_angstrom <= STRUCTURAL_ATTACHMENT_CUTOFF_ANGSTROM
    ]
    if likely_missing_evidence:
        examples = ", ".join(residue.residue_id for residue in likely_missing_evidence[:8])
        warnings.append(
            "Likely modified glycan/polymer/PTM residues are close to protein heavy atoms but "
            "lack LINK/CONECT evidence; add CCD-compatible connectivity records before "
            f"production ingestion. Examples: {examples}"
        )

    polymer_names = sorted(
        {residue.residue_name for residue in residues if residue.category == "polymer_ptm"}
    )
    if polymer_names:
        warnings.append(
            "Polymer/PTM-like residue names detected for diagnostics-only handling: "
            f"{', '.join(polymer_names)}"
        )
    return warnings


def _classify_residue(residue_name: str) -> str:
    """Classify a residue name for preflight diagnostics."""
    normalized = (residue_name or "").upper()
    if _is_protein_like_residue(normalized):
        return "protein"
    if normalized in WATER_RESIDUES:
        return "water"
    if normalized in COMMON_ION_RESIDUES:
        return "ion"
    if normalized in COMMON_SOLVENT_RESIDUES:
        return "solvent"
    if normalized in GLYCAN_RESIDUES:
        return "glycan"
    if normalized in POLYMER_LIKE_RESIDUES:
        return "polymer_ptm"
    return "ligand_cocrystal"


def _is_noncanonical_residue(residue_name: str) -> bool:
    """Return whether a residue should be considered noncanonical for protein ingestion."""
    category = _classify_residue(residue_name)
    return category not in {"protein", "water", "ion", "solvent"}


def _is_protein_like_residue(residue_name: str) -> bool:
    """Return whether a residue should be handled as protein-like."""
    normalized = (residue_name or "").upper()
    return normalized in STANDARD_PROTEIN_RESIDUES or normalized in MODIFIED_PROTEIN_RESIDUES


def _is_attachment_candidate_residue(residue: PDBResidueInspection) -> bool:
    """Return whether a residue can be a linked PTM/glycan/polymer candidate."""
    if not residue.is_noncanonical:
        return False
    return residue.category not in {"water", "ion", "solvent"}


def _format_residue_id(
    chain_id: str,
    residue_name: str,
    residue_number: int | None,
    res_seq: str | None,
    insertion_code: str | None,
) -> str:
    """Format a stable residue identifier for diagnostics."""
    number = residue_number if residue_number is not None else res_seq
    if number in (None, ""):
        number = "?"
    return f"{chain_id or '_'}:{residue_name}{number}{insertion_code or ''}"


def _has_coordinates(atom: _DiagnosticAtomRecord) -> bool:
    """Return whether an atom has parsed Cartesian coordinates."""
    return atom.x is not None and atom.y is not None and atom.z is not None


def _has_heavy_atom_coordinates(atom: _DiagnosticAtomRecord) -> bool:
    """Return whether an atom has non-hydrogen coordinates."""
    return _has_coordinates(atom) and (atom.element or "").upper() != "H"


def _parse_int(value: Any) -> int | None:
    """Parse an integer-like value and return ``None`` on failure."""
    return parse_int(value)
