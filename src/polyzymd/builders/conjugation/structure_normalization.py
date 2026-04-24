"""Clean-PDB chain normalization planning and writer utilities."""

from __future__ import annotations

from collections import defaultdict, deque
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.structure_inspection import (
    POLYZMD_MOIETY_CHAIN,
    POLYZMD_PROTEIN_CHAIN,
    PDBAtomRecord,
    PDBResidueInspection,
    PDBStructureInspection,
    _format_residue_id,
    _parse_atom_records,
    _parse_conect_targets,
    _parse_int,
    _parse_link_sides,
    _summarize_residues,
    inspect_pdb_structure,
)

ERROR_SEVERITY = "error"


class PDBCleanlinessIssue(BaseModel):
    """Clean-PDB validation issue for normalization preflight.

    The model is JSON-safe and intentionally keeps future opt-in support for
    intentional cofactors or metals as message details rather than behavior.
    """

    code: str
    severity: str = ERROR_SEVERITY
    message: str
    residue_id: str | None = None
    residue_name: str | None = None
    chain_id: str = ""
    category: str | None = None
    details: dict[str, Any] = Field(default_factory=dict)


class PDBChainNormalizationAction(BaseModel):
    """Residue-level chain reassignment planned for a normalized PDB copy."""

    residue_id: str
    residue_name: str
    residue_number: int | None = None
    res_seq: str | None = None
    insertion_code: str | None = None
    source_chain: str = ""
    target_chain: str
    target_residue_number: int
    target_res_seq: str
    target_insertion_code: str | None = None
    category: str
    reason: str


class PDBNormalizationPlan(BaseModel):
    """Clean-PDB normalization plan for PolyzyMD chain convention.

    A valid plan may be written to a separate normalized PDB copy. Invalid
    plans are diagnostics-only and should prompt users to strip waters, ions,
    solvents, or unsupported free components before ingestion.
    """

    path: Path
    clean: bool
    valid: bool
    protein_residue_count: int = 0
    moiety_residue_count: int = 0
    issues: list[PDBCleanlinessIssue] = Field(default_factory=list)
    actions: list[PDBChainNormalizationAction] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    source_chains_collapsed: dict[str, list[str]] = Field(default_factory=dict)
    chain_id_change_count: int = 0
    residue_number_change_count: int = 0
    default_output_path: Path | None = None
    output_recommendations: list[str] = Field(default_factory=list)


def default_cleaned_pdb_path(path: Path | str) -> Path:
    """Return the default normalized PDB output path.

    Parameters
    ----------
    path : Path or str
        Input PDB path.

    Returns
    -------
    Path
        Path beside the input named ``<stem>_cleaned.pdb``.
    """
    pdb_path = Path(path)
    return pdb_path.with_name(f"{pdb_path.stem}_cleaned{pdb_path.suffix or '.pdb'}")


def plan_pdb_chain_normalization(
    path_or_inspection: Path | str | PDBStructureInspection,
) -> PDBNormalizationPlan:
    """Plan clean-PDB chain normalization without mutating input files.

    Parameters
    ----------
    path_or_inspection : Path, str, or PDBStructureInspection
        PDB path or a precomputed pure-Python structure inspection.

    Returns
    -------
    PDBNormalizationPlan
        JSON-safe plan containing cleanliness issues, warnings, and residue
        actions that would normalize accepted residues to chains A and C.
    """
    inspection = _coerce_inspection(path_or_inspection)
    lines = inspection.path.read_text(errors="replace").splitlines()
    atoms = _parse_atom_records(lines)
    residues = _summarize_residues(atoms)
    residue_lookup = {residue.residue_id: residue for residue in residues}
    graph = _build_covalent_residue_graph(lines, atoms, residue_lookup)
    protein_residue_ids = {
        residue.residue_id for residue in residues if residue.is_standard_protein
    }
    reachable_residue_ids = _reachable_residue_ids(protein_residue_ids, graph)

    issues: list[PDBCleanlinessIssue] = []
    actions: list[PDBChainNormalizationAction] = []
    if not protein_residue_ids:
        issues.append(
            PDBCleanlinessIssue(
                code="no_protein_residues",
                message=(
                    "Clean-PDB normalization requires at least one canonical protein residue; "
                    "no protein-like residues were found"
                ),
                details={"action": "Provide the protein entity PDB before conjugation preflight"},
            )
        )

    for residue in residues:
        if residue.is_standard_protein:
            actions.append(_protein_action(residue))
            continue
        if residue.category in {"water", "ion", "solvent"}:
            issues.append(_dirty_component_issue(residue))
            continue
        if residue.residue_id in reachable_residue_ids:
            actions.append(_moiety_action(residue))
            continue
        issues.append(_unlinked_noncanonical_issue(residue))

    _assign_target_residue_numbers(actions)
    error_issues = [issue for issue in issues if issue.severity == ERROR_SEVERITY]
    warnings = _build_normalization_warnings(actions)
    clean = not error_issues
    return PDBNormalizationPlan(
        path=inspection.path,
        clean=clean,
        valid=clean,
        protein_residue_count=sum(
            1 for action in actions if action.target_chain == POLYZMD_PROTEIN_CHAIN
        ),
        moiety_residue_count=sum(
            1 for action in actions if action.target_chain == POLYZMD_MOIETY_CHAIN
        ),
        issues=issues,
        actions=actions,
        warnings=warnings,
        source_chains_collapsed=_source_chains_by_target(actions),
        chain_id_change_count=sum(1 for action in actions if _chain_id_changed(action)),
        residue_number_change_count=sum(
            1 for action in actions if _residue_number_or_insertion_changed(action)
        ),
        default_output_path=default_cleaned_pdb_path(inspection.path),
        output_recommendations=_build_output_recommendations(inspection.path, clean),
    )


def write_normalized_pdb(
    input_path: Path | str,
    output_path: Path | str,
    plan: PDBNormalizationPlan | None = None,
) -> Path:
    """Write a normalized PDB copy using a validated normalization plan.

    Parameters
    ----------
    input_path : Path or str
        Source PDB path. The file is never modified in place.
    output_path : Path or str
        Explicit destination path for the normalized copy.
    plan : PDBNormalizationPlan or None, optional
        Precomputed normalization plan, by default ``None``.

    Returns
    -------
    Path
        Destination path written.

    Raises
    ------
    ValueError
        If the output path is omitted, points at the input file, or the plan
        contains clean-PDB errors.
    """
    source_path = Path(input_path)
    destination_path = Path(output_path) if output_path is not None else None
    if destination_path is None:
        raise ValueError("write_normalized_pdb requires an explicit output path")
    if source_path.resolve() == destination_path.resolve():
        raise ValueError("write_normalized_pdb never mutates input; choose a separate output path")

    normalization_plan = plan or plan_pdb_chain_normalization(source_path)
    if not normalization_plan.valid:
        issue_summary = "; ".join(issue.message for issue in normalization_plan.issues[:5])
        raise ValueError(
            f"Refusing to write normalized PDB because clean-PDB validation failed: {issue_summary}"
        )

    action_lookup = {_action_key(action): action for action in normalization_plan.actions}
    output_lines = [
        _rewrite_pdb_line(line, action_lookup)
        for line in source_path.read_text(errors="replace").splitlines(keepends=True)
    ]
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    destination_path.write_text("".join(output_lines))
    return destination_path


def _coerce_inspection(
    path_or_inspection: Path | str | PDBStructureInspection,
) -> PDBStructureInspection:
    """Return a PDB inspection from a path or existing inspection object."""
    if isinstance(path_or_inspection, PDBStructureInspection):
        return path_or_inspection
    return inspect_pdb_structure(path_or_inspection)


def _build_covalent_residue_graph(
    lines: list[str],
    atoms: list[PDBAtomRecord],
    residue_lookup: dict[str, PDBResidueInspection],
) -> dict[str, set[str]]:
    """Build residue connectivity from LINK and CONECT records."""
    graph: dict[str, set[str]] = defaultdict(set)
    serial_lookup = {atom.atom_serial: atom for atom in atoms if atom.atom_serial is not None}
    for line in lines:
        if line.startswith("LINK"):
            residue_ids = _residue_ids_from_link(line, residue_lookup)
        elif line.startswith("CONECT"):
            residue_ids = _residue_id_pairs_from_conect(line, serial_lookup)
        else:
            continue
        for residue_id_1, residue_id_2 in residue_ids:
            if residue_id_1 == residue_id_2:
                continue
            graph[residue_id_1].add(residue_id_2)
            graph[residue_id_2].add(residue_id_1)
    return dict(graph)


def _residue_ids_from_link(
    line: str,
    residue_lookup: dict[str, PDBResidueInspection],
) -> list[tuple[str, str]]:
    """Return existing residue IDs connected by one LINK record."""
    side_1, side_2 = _parse_link_sides(line)
    residue_id_1 = _residue_id_from_link_side(side_1)
    residue_id_2 = _residue_id_from_link_side(side_2)
    if residue_id_1 not in residue_lookup or residue_id_2 not in residue_lookup:
        return []
    return [(residue_id_1, residue_id_2)]


def _residue_id_pairs_from_conect(
    line: str,
    serial_lookup: dict[int, PDBAtomRecord],
) -> list[tuple[str, str]]:
    """Return residue ID pairs connected by one CONECT record."""
    source_serial = _parse_int(line[6:11].strip())
    if source_serial is None:
        return []
    source_atom = serial_lookup.get(source_serial)
    if source_atom is None:
        return []
    source_residue_id = _residue_id_from_atom(source_atom)
    pairs: list[tuple[str, str]] = []
    for target_serial in _parse_conect_targets(line):
        target_atom = serial_lookup.get(target_serial)
        if target_atom is None:
            continue
        pairs.append((source_residue_id, _residue_id_from_atom(target_atom)))
    return pairs


def _reachable_residue_ids(
    starting_residue_ids: set[str],
    graph: dict[str, set[str]],
) -> set[str]:
    """Return residues with explicit covalent connectivity to protein residues."""
    reachable = set(starting_residue_ids)
    queue = deque(starting_residue_ids)
    while queue:
        residue_id = queue.popleft()
        for neighbor in graph.get(residue_id, set()):
            if neighbor in reachable:
                continue
            reachable.add(neighbor)
            queue.append(neighbor)
    return reachable


def _protein_action(residue: PDBResidueInspection) -> PDBChainNormalizationAction:
    """Build a chain-A action for a canonical protein residue."""
    return _action_for_residue(
        residue,
        target_chain=POLYZMD_PROTEIN_CHAIN,
        reason="Canonical protein residue normalized to PolyzyMD chain A",
    )


def _moiety_action(residue: PDBResidueInspection) -> PDBChainNormalizationAction:
    """Build a chain-C action for a covalently attached noncanonical residue."""
    return _action_for_residue(
        residue,
        target_chain=POLYZMD_MOIETY_CHAIN,
        reason=(
            "Covalently connected glycan/PTM/polymer/noncanonical residue normalized to "
            "PolyzyMD chain C"
        ),
    )


def _action_for_residue(
    residue: PDBResidueInspection,
    *,
    target_chain: str,
    reason: str,
) -> PDBChainNormalizationAction:
    """Build a residue-level normalization action."""
    return PDBChainNormalizationAction(
        residue_id=residue.residue_id,
        residue_name=residue.residue_name,
        residue_number=residue.residue_number,
        res_seq=residue.res_seq,
        insertion_code=residue.insertion_code,
        source_chain=residue.chain_id,
        target_chain=target_chain,
        target_residue_number=0,
        target_res_seq="0",
        target_insertion_code=None,
        category=residue.category,
        reason=reason,
    )


def _assign_target_residue_numbers(actions: list[PDBChainNormalizationAction]) -> None:
    """Assign first-seen continuous residue numbers within output chains."""
    next_residue_number = {
        POLYZMD_PROTEIN_CHAIN: 1,
        POLYZMD_MOIETY_CHAIN: 1,
    }
    for action in actions:
        target_number = next_residue_number[action.target_chain]
        action.target_residue_number = target_number
        action.target_res_seq = str(target_number)
        action.target_insertion_code = None
        next_residue_number[action.target_chain] = target_number + 1


def _source_chains_by_target(actions: list[PDBChainNormalizationAction]) -> dict[str, list[str]]:
    """Return source chains collapsed into each target chain in first-seen order."""
    collapsed: dict[str, list[str]] = {}
    for action in actions:
        chains = collapsed.setdefault(action.target_chain, [])
        source_chain = action.source_chain or "blank"
        if source_chain not in chains:
            chains.append(source_chain)
    return collapsed


def _chain_id_changed(action: PDBChainNormalizationAction) -> bool:
    """Return whether a residue changes chain ID in the normalized copy."""
    return action.source_chain != action.target_chain


def _residue_number_or_insertion_changed(action: PDBChainNormalizationAction) -> bool:
    """Return whether a residue changes residue number or insertion code."""
    return action.res_seq != action.target_res_seq or (action.insertion_code or None) != (
        action.target_insertion_code or None
    )


def _dirty_component_issue(residue: PDBResidueInspection) -> PDBCleanlinessIssue:
    """Build an issue for waters, ions, and common solvent residues."""
    return PDBCleanlinessIssue(
        code="dirty_component",
        message=(
            "Clean-PDB input must not contain crystallographic waters, ions, or common "
            f"solvents; strip residue {residue.residue_id} before normalization"
        ),
        residue_id=residue.residue_id,
        residue_name=residue.residue_name,
        chain_id=residue.chain_id,
        category=residue.category,
        details={
            "action": "Remove waters, ions, and solvent records before this ingestion path",
            "future_support": "Intentional bound cofactors or metals will require explicit opt-in",
        },
    )


def _unlinked_noncanonical_issue(residue: PDBResidueInspection) -> PDBCleanlinessIssue:
    """Build an issue for unsupported free or unlinked noncanonical residues."""
    if residue.is_polymer_ptm_candidate:
        code = "missing_covalent_evidence"
        message = (
            "Glycan/PTM/polymer-like residue lacks LINK/CONECT evidence to the protein entity; "
            f"add covalent connectivity or remove residue {residue.residue_id}"
        )
    else:
        code = "free_ligand_or_unlinked_component"
        message = (
            "Clean-PDB input currently supports only protein plus covalently attached "
            f"noncanonical residues; remove unsupported free component {residue.residue_id}"
        )
    return PDBCleanlinessIssue(
        code=code,
        message=message,
        residue_id=residue.residue_id,
        residue_name=residue.residue_name,
        chain_id=residue.chain_id,
        category=residue.category,
        details={
            "action": "Provide one covalently connected protein entity for automatic chain reassignment",
            "future_support": "Intentional bound cofactors or metals will require explicit opt-in",
        },
    )


def _build_normalization_warnings(actions: list[PDBChainNormalizationAction]) -> list[str]:
    """Build warning-only messages for chain convention changes."""
    warnings: list[str] = []
    blank_actions = [action for action in actions if not action.source_chain]
    if blank_actions:
        warnings.append(
            "Blank chain IDs are valid only for a clean single covalent entity; normalized copy "
            "will assign accepted residues to chains A and C"
        )

    protein_source_chains = sorted(
        {
            action.source_chain or "blank"
            for action in actions
            if action.target_chain == POLYZMD_PROTEIN_CHAIN
            and action.source_chain != POLYZMD_PROTEIN_CHAIN
        }
    )
    if protein_source_chains:
        warnings.append(
            "Canonical protein residues will be normalized to chain A from chain(s): "
            f"{', '.join(protein_source_chains)}"
        )

    moiety_source_chains = sorted(
        {
            action.source_chain or "blank"
            for action in actions
            if action.target_chain == POLYZMD_MOIETY_CHAIN
            and action.source_chain != POLYZMD_MOIETY_CHAIN
        }
    )
    if moiety_source_chains:
        warnings.append(
            "Covalently attached noncanonical residues will be normalized to chain C from "
            f"chain(s): {', '.join(moiety_source_chains)}"
        )
    residue_number_changes = sum(
        1 for action in actions if _residue_number_or_insertion_changed(action)
    )
    if residue_number_changes:
        warnings.append(
            "Residue numbers and insertion codes will be renumbered continuously within "
            f"normalized chains ({residue_number_changes} residue mappings changed)"
        )
    return warnings


def _build_output_recommendations(path: Path, clean: bool) -> list[str]:
    """Build actionable output recommendations for a normalization plan."""
    cleaned_name = default_cleaned_pdb_path(path).name
    if clean:
        return [
            "Plan is valid for an explicit normalized copy; use write_normalized_pdb with an "
            f"output path such as {cleaned_name}",
            "Input PDB is never mutated by the normalization writer",
        ]
    return [
        "Do not write a normalized copy until clean-PDB errors are resolved",
        "Strip waters, ions, solvents, and unsupported free ligands before this ingestion path",
    ]


def _rewrite_pdb_line(
    line: str,
    action_lookup: dict[tuple[str, str, str | None, str | None], PDBChainNormalizationAction],
) -> str:
    """Rewrite supported PDB records according to residue normalization actions."""
    if line.startswith(("ATOM  ", "HETATM")):
        return _rewrite_atom_residue_identity(line, action_lookup)
    if line.startswith("LINK"):
        return _rewrite_link_residue_identity(line, action_lookup)
    return line


def _rewrite_atom_residue_identity(
    line: str,
    action_lookup: dict[tuple[str, str, str | None, str | None], PDBChainNormalizationAction],
) -> str:
    """Rewrite ATOM/HETATM chain, residue number, and insertion-code fields."""
    body = line.rstrip("\r\n")
    newline = line[len(body) :]
    action = action_lookup.get(_line_residue_key(body))
    if action is None:
        return line
    padded = body.ljust(27)
    target_res_seq = _format_target_res_seq(action)
    target_insertion_code = action.target_insertion_code or " "
    return (
        f"{padded[:21]}{action.target_chain}{target_res_seq}"
        f"{target_insertion_code}{padded[27:]}{newline}"
    )


def _rewrite_link_residue_identity(
    line: str,
    action_lookup: dict[tuple[str, str, str | None, str | None], PDBChainNormalizationAction],
) -> str:
    """Rewrite LINK residue chain, residue number, and insertion-code fields."""
    side_1, side_2 = _parse_link_sides(line.rstrip("\r\n"))
    action_1 = action_lookup.get(_link_side_key(side_1))
    action_2 = action_lookup.get(_link_side_key(side_2))
    if action_1 is None and action_2 is None:
        return line

    body = line.rstrip("\r\n")
    newline = line[len(body) :]
    padded = body.ljust(57)
    if action_1 is not None:
        padded = _replace_link_side(
            padded, action_1, chain_index=21, res_seq_start=22, icode_index=26
        )
    if action_2 is not None:
        padded = _replace_link_side(
            padded, action_2, chain_index=51, res_seq_start=52, icode_index=56
        )
    return f"{padded.rstrip()}{newline}"


def _replace_link_side(
    line: str,
    action: PDBChainNormalizationAction,
    *,
    chain_index: int,
    res_seq_start: int,
    icode_index: int,
) -> str:
    """Return a LINK record with one residue side updated in fixed-width columns."""
    target_res_seq = _format_target_res_seq(action)
    target_insertion_code = action.target_insertion_code or " "
    return (
        f"{line[:chain_index]}{action.target_chain}"
        f"{line[chain_index + 1 : res_seq_start]}{target_res_seq}"
        f"{target_insertion_code}{line[icode_index + 1 :]}"
    )


def _format_target_res_seq(action: PDBChainNormalizationAction) -> str:
    """Format a target residue number for a four-character PDB residue field."""
    if action.target_residue_number < -999 or action.target_residue_number > 9999:
        raise ValueError(
            "PDB residue numbers must fit the four-character resSeq field; "
            f"target {action.target_residue_number} for {action.residue_id} is unsupported"
        )
    return f"{action.target_residue_number:4d}"


def _action_key(action: PDBChainNormalizationAction) -> tuple[str, str, str | None, str | None]:
    """Return a residue identity key from a normalization action."""
    return (
        action.source_chain,
        action.residue_name,
        action.res_seq,
        action.insertion_code,
    )


def _line_residue_key(line: str) -> tuple[str, str, str | None, str | None]:
    """Return a residue identity key from an ATOM/HETATM record line."""
    return (
        line[21:22].strip(),
        (line[17:20].strip() or "").upper(),
        line[22:26].strip() or None,
        line[26:27].strip() or None,
    )


def _link_side_key(side: dict[str, Any]) -> tuple[str, str, str | None, str | None]:
    """Return a residue identity key from a parsed LINK side."""
    return (
        side["chain_id"],
        side["residue_name"],
        side["res_seq"],
        side["insertion_code"],
    )


def _residue_id_from_link_side(side: dict[str, Any]) -> str:
    """Format a residue ID from a parsed LINK side."""
    return _format_residue_id(
        side["chain_id"],
        side["residue_name"],
        side["residue_number"],
        side["res_seq"],
        side["insertion_code"],
    )


def _residue_id_from_atom(atom: PDBAtomRecord) -> str:
    """Format a residue ID from an atom record."""
    return _format_residue_id(
        atom.chain_id,
        atom.residue_name,
        atom.residue_number,
        atom.res_seq,
        atom.insertion_code,
    )
