"""RDKit graph helpers retained for NHS-Lys reaction tests and diagnostics."""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from pydantic import BaseModel, Field


class AddedBond(BaseModel):
    """Atom-index record for a bond created by graph surgery."""

    begin_atom_index: int = Field(..., ge=0, description="Begin atom index in the product molecule")
    end_atom_index: int = Field(..., ge=0, description="End atom index in the product molecule")
    order: int = Field(1, ge=1, description="Integer bond order used by the primitive")


class RdkitGraphEditResult(BaseModel):
    """Structured result from an RDKit graph edit primitive."""

    model_config = {"arbitrary_types_allowed": True}

    product_mol: Any = Field(..., exclude=True, description="Product RDKit molecule")
    removed_protein_atom_indices: tuple[int, ...] = Field(
        default_factory=tuple,
        description="Original protein atom indices removed during graph surgery",
    )
    removed_moiety_atom_indices: tuple[int, ...] = Field(
        default_factory=tuple,
        description="Original moiety atom indices removed during graph surgery",
    )
    added_bond: AddedBond = Field(..., description="Product atom indices for the created bond")
    protein_atom_index_map: dict[int, int] = Field(
        default_factory=dict,
        description="Original-to-product atom index map for retained protein atoms",
    )
    moiety_atom_index_map: dict[int, int] = Field(
        default_factory=dict,
        description="Original-to-product atom index map for retained moiety atoms",
    )
    warnings: tuple[str, ...] = Field(
        default_factory=tuple,
        description="Non-fatal warnings emitted by the conservative executor",
    )


class NhsLysAttachmentSite(BaseModel):
    """Normalized explicit lysine NZ attachment site for NHS-Lys graph helpers."""

    chain_id: str = Field(..., min_length=1, max_length=1)
    residue_name: str = Field(..., min_length=1)
    residue_number: int
    atom_name: str = Field(..., min_length=1)


class LysineReactiveSite(BaseModel):
    """Explicit lysine NZ reactive site resolved to topology atom indices."""

    chain_id: str = Field(..., min_length=1, description="Chain identifier")
    residue_name: str = Field("LYS", description="Residue name")
    residue_number: int = Field(..., description="Residue number in the input topology")
    atom_name: str = Field("NZ", description="Reactive atom name")
    nz_atom_index: int = Field(..., ge=0, description="Topology atom index for lysine NZ")
    nz_hydrogen_indices: tuple[int, ...] = Field(
        default_factory=tuple,
        description="Topology atom indices for hydrogens bonded to NZ",
    )
    ce_atom_index: int | None = Field(None, ge=0, description="Optional CE atom index")
    nz_position: tuple[float, float, float] | None = Field(
        None,
        description="Optional NZ Cartesian position",
    )
    ce_position: tuple[float, float, float] | None = Field(
        None,
        description="Optional CE Cartesian position",
    )
    evidence: dict[str, Any] = Field(
        default_factory=dict,
        description="Metadata evidence used to resolve the site",
    )
    warnings: tuple[str, ...] = Field(
        default_factory=tuple,
        description="Non-fatal site extraction warnings",
    )


class NhsReactiveGroup(BaseModel):
    """Reactive carbon and leaving group indices for an NHS ester moiety."""

    reactive_carbon_index: int = Field(..., ge=0, description="Ester acyl carbon atom index")
    bridging_oxygen_index: int = Field(
        ..., ge=0, description="Oxygen connecting acyl carbon to NHS"
    )
    leaving_group_atom_indices: tuple[int, ...] = Field(
        ..., min_length=1, description="Atom indices removed as the NHS leaving group"
    )
    retained_atom_indices: tuple[int, ...] = Field(
        ..., min_length=1, description="Atom indices retained after leaving group removal"
    )
    candidate_atom_indices: tuple[int, ...] = Field(
        default_factory=tuple,
        description="Atom indices searched when identifying the reactive group",
    )
    evidence: dict[str, Any] = Field(
        default_factory=dict,
        description="Topology evidence supporting the assignment",
    )
    diagnostics: tuple[str, ...] = Field(
        default_factory=tuple,
        description="Human-readable detection diagnostics",
    )


class NhsLysGraphEditPlan(BaseModel):
    """Typed graph edit plan for the ``nhs_lys_amide`` primitive."""

    mechanism: str = Field("nhs_lys_amide", description="Mechanism identifier")
    site: LysineReactiveSite = Field(..., description="Resolved lysine site")
    reactive_group: NhsReactiveGroup = Field(..., description="Resolved NHS reactive group")
    site_atom_index: int = Field(..., ge=0, description="Protein NZ atom index before graph edit")
    reactive_carbon_index: int = Field(..., ge=0, description="Moiety acyl carbon index")
    site_hydrogen_indices_to_remove: tuple[int, ...] = Field(
        default_factory=tuple,
        description="Protein hydrogen atom indices selected for removal",
    )
    leaving_group_atom_indices: tuple[int, ...] = Field(
        ..., min_length=1, description="Moiety leaving group atom indices selected for removal"
    )
    add_bond: tuple[int, int] = Field(
        ..., description="Pre-edit local atom indices for the NZ-carbon bond"
    )
    executable_with_rdkit: bool = Field(
        True,
        description="Whether the plan has enough explicit atom indices for RDKit graph surgery",
    )
    warnings: tuple[str, ...] = Field(
        default_factory=tuple,
        description="Non-fatal planning warnings",
    )


def extract_lysine_reactive_site(
    site: NhsLysAttachmentSite,
    atoms: Iterable[Any],
    *,
    bonds: Iterable[Any] | None = None,
    positions: Mapping[int, Sequence[float]] | Sequence[Sequence[float]] | None = None,
    require_hydrogens: bool = False,
) -> LysineReactiveSite:
    """Resolve an explicit lysine NZ attachment site from topology atom metadata.

    Parameters
    ----------
    site : NhsLysAttachmentSite
        Normalized explicit attachment site to resolve.
    atoms : Iterable[Any]
        Atom records from a topology-like object. The resolver supports dict
        records, OpenFF-like atom metadata, and RDKit-like atom objects without
        importing those packages at module load.
    bonds : Iterable[Any] or None, optional
        Optional topology bond records used to identify NZ-bound hydrogens and
        CE, by default ``None``.
    positions : Mapping[int, Sequence[float]] or Sequence[Sequence[float]] or None, optional
        Optional Cartesian positions keyed by atom index or ordered by atom
        index, by default ``None``.
    require_hydrogens : bool, optional
        Require at least one NZ-bound hydrogen when bond information is
        available, by default ``False``.

    Returns
    -------
    LysineReactiveSite
        Resolved lysine site with atom indices and optional bonded hydrogens.

    Raises
    ------
    ValueError
        If the attachment site is not an explicit lysine NZ site or cannot be
        matched unambiguously.
    """
    if site.residue_name != "LYS" or site.atom_name != "NZ":
        raise ValueError(
            "NHS-Lys reactive site extraction requires an explicit LYS:NZ attachment site"
        )

    atom_records = list(atoms)
    indexed_atoms = [
        (_atom_index(atom, fallback), atom) for fallback, atom in enumerate(atom_records)
    ]
    matching_atoms = [
        (index, atom)
        for index, atom in indexed_atoms
        if _atom_name(atom) == "NZ"
        and _residue_name(atom) == "LYS"
        and _residue_number(atom) == site.residue_number
        and _chain_id(atom) == site.chain_id
    ]

    if not matching_atoms:
        raise ValueError(
            f"No lysine NZ atom matched site {site.chain_id}:LYS{site.residue_number}:NZ"
        )
    if len(matching_atoms) > 1:
        indices = ", ".join(str(index) for index, _atom in matching_atoms)
        raise ValueError(
            f"Ambiguous lysine NZ site {site.chain_id}:LYS{site.residue_number}:NZ matched "
            f"multiple atoms: {indices}"
        )

    nz_index, nz_atom = matching_atoms[0]
    warnings: list[str] = []
    bond_pairs = _collect_bond_pairs(atom_records, bonds=bonds)
    neighbor_indices = _neighbor_indices(nz_index, bond_pairs)
    atom_by_index = dict(indexed_atoms)

    hydrogen_indices = tuple(
        sorted(
            index
            for index in neighbor_indices
            if index in atom_by_index and _atomic_number(atom_by_index[index]) == 1
        )
    )
    ce_indices = sorted(
        index
        for index in neighbor_indices
        if index in atom_by_index and _atom_name(atom_by_index[index]) == "CE"
    )
    ce_index = ce_indices[0] if ce_indices else None

    if not bond_pairs:
        warnings.append("No topology bonds were available; NZ-bound hydrogens were not resolved")
    elif require_hydrogens and not hydrogen_indices:
        raise ValueError(
            f"Lysine site {site.chain_id}:LYS{site.residue_number}:NZ has no NZ-bound hydrogens"
        )
    if bond_pairs and ce_index is None:
        warnings.append("Topology bonds were available, but CE was not bonded to NZ")

    return LysineReactiveSite(
        chain_id=site.chain_id,
        residue_name="LYS",
        residue_number=site.residue_number,
        atom_name="NZ",
        nz_atom_index=nz_index,
        nz_hydrogen_indices=hydrogen_indices,
        ce_atom_index=ce_index,
        nz_position=_position_for_index(nz_index, positions),
        ce_position=_position_for_index(ce_index, positions) if ce_index is not None else None,
        evidence={
            "matched_atom_index": nz_index,
            "matched_atom_name": _atom_name(nz_atom),
            "matched_chain_id": _chain_id(nz_atom),
            "matched_residue_number": _residue_number(nz_atom),
            "matched_residue_name": _residue_name(nz_atom),
            "bond_count": len(bond_pairs),
        },
        warnings=tuple(warnings),
    )


def detect_nhs_reactive_group(
    mol: Any,
    *,
    candidate_atom_indices: Iterable[int] | None = None,
    residue_name_prefixes: tuple[str, ...] = ("NHS", "NH"),
) -> NhsReactiveGroup:
    """Detect the reactive acyl carbon and leaving group of an NHS ester.

    Parameters
    ----------
    mol : Any
        RDKit molecule-like object with explicit graph connectivity.
    candidate_atom_indices : Iterable[int] or None, optional
        Optional atom indices limiting the search region. If omitted, PDB
        residue metadata matching ``residue_name_prefixes`` is used when
        present, otherwise all atoms are searched.
    residue_name_prefixes : tuple[str, ...], optional
        PDB residue name prefixes that indicate an NHS residue, by default
        ``("NHS", "NH")``.

    Returns
    -------
    NhsReactiveGroup
        Detected acyl carbon, bridging oxygen, leaving group, and retained atom
        indices.

    Raises
    ------
    ValueError
        If no NHS ester-like group is found or if multiple assignments are
        possible.
    """
    atom_count = mol.GetNumAtoms()
    search_indices = _resolve_nhs_search_indices(
        mol,
        candidate_atom_indices=candidate_atom_indices,
        residue_name_prefixes=residue_name_prefixes,
    )
    validate_atom_indices(mol, search_indices, label="NHS search")

    candidates: list[NhsReactiveGroup] = []
    diagnostics: list[str] = []
    for atom_index in sorted(search_indices):
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetSymbol() != "C":
            continue
        oxygen_neighbors = [
            neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == "O"
        ]
        if len(oxygen_neighbors) != 2:
            continue

        for oxygen_index in oxygen_neighbors:
            oxygen_atom = mol.GetAtomWithIdx(oxygen_index)
            non_carbon_neighbors = [
                neighbor.GetIdx()
                for neighbor in oxygen_atom.GetNeighbors()
                if neighbor.GetIdx() != atom_index
            ]
            has_nhs_nitrogen = any(
                mol.GetAtomWithIdx(neighbor_index).GetSymbol() == "N"
                for neighbor_index in non_carbon_neighbors
            )
            if not has_nhs_nitrogen:
                continue

            leaving_group = _traverse_leaving_group(mol, oxygen_index, excluded_index=atom_index)
            if not _leaving_group_has_nhs_evidence(mol, leaving_group):
                diagnostics.append(
                    f"Candidate carbon {atom_index} oxygen {oxygen_index} lacked NHS-like "
                    "leaving-group evidence"
                )
                continue

            retained = tuple(index for index in range(atom_count) if index not in leaving_group)
            candidates.append(
                NhsReactiveGroup(
                    reactive_carbon_index=atom_index,
                    bridging_oxygen_index=oxygen_index,
                    leaving_group_atom_indices=tuple(sorted(leaving_group)),
                    retained_atom_indices=retained,
                    candidate_atom_indices=tuple(sorted(search_indices)),
                    evidence={
                        "oxygen_neighbor_indices": tuple(sorted(oxygen_neighbors)),
                        "bridging_oxygen_neighbor_indices": tuple(sorted(non_carbon_neighbors)),
                        "leaving_group_contains_nitrogen": True,
                        "leaving_group_atom_count": len(leaving_group),
                    },
                    diagnostics=tuple(diagnostics),
                )
            )

    if not candidates:
        detail = f" Diagnostics: {'; '.join(diagnostics)}" if diagnostics else ""
        raise ValueError(
            "No unambiguous NHS ester reactive group was found. Expected an acyl carbon bonded "
            f"to two oxygens with one oxygen connected to an NHS nitrogen.{detail}"
        )
    unique_candidates = _deduplicate_reactive_groups(candidates)
    if len(unique_candidates) > 1:
        assignments = ", ".join(
            f"C{candidate.reactive_carbon_index}-O{candidate.bridging_oxygen_index}"
            for candidate in unique_candidates
        )
        raise ValueError(f"Ambiguous NHS ester reactive group assignment: {assignments}")
    return unique_candidates[0]


def plan_nhs_lys_amide(
    site: LysineReactiveSite,
    reactive_group: NhsReactiveGroup,
    *,
    site_hydrogen_indices_to_remove: Iterable[int] | None = None,
) -> NhsLysGraphEditPlan:
    """Create a typed graph edit plan for NHS-Lys amide formation.

    Parameters
    ----------
    site : LysineReactiveSite
        Resolved lysine NZ site.
    reactive_group : NhsReactiveGroup
        Resolved NHS ester reactive group on the moiety.
    site_hydrogen_indices_to_remove : Iterable[int] or None, optional
        Explicit NZ-bound hydrogens to remove. If omitted, the initial product
        proton policy keeps the first deterministically ordered hydrogen and
        removes one hydrogen for neutral lysine or two for protonated lysine, by
        default ``None``.

    Returns
    -------
    NhsLysGraphEditPlan
        Typed graph edit plan with explicit atom indices.
    """
    hydrogens = tuple(
        sorted(
            site_hydrogen_indices_to_remove
            if site_hydrogen_indices_to_remove is not None
            else _default_nz_hydrogens_to_remove(site.nz_hydrogen_indices)
        )
    )
    warnings: list[str] = []
    if not hydrogens:
        warnings.append("No site hydrogens were selected for removal; proton handling is deferred")

    return NhsLysGraphEditPlan(
        site=site,
        reactive_group=reactive_group,
        site_atom_index=site.nz_atom_index,
        reactive_carbon_index=reactive_group.reactive_carbon_index,
        site_hydrogen_indices_to_remove=hydrogens,
        leaving_group_atom_indices=reactive_group.leaving_group_atom_indices,
        add_bond=(site.nz_atom_index, reactive_group.reactive_carbon_index),
        warnings=tuple(warnings),
    )


def _default_nz_hydrogens_to_remove(hydrogen_indices: Iterable[int]) -> tuple[int, ...]:
    """Apply the initial NHS-Lys product proton policy to resolved hydrogens."""
    hydrogens = tuple(sorted(hydrogen_indices))
    if len(hydrogens) == 3:
        return hydrogens[1:]
    if len(hydrogens) == 2:
        return hydrogens[1:]
    raise ValueError(
        "NHS-Lys planning requires explicit lysine NZ hydrogens. Expected 2 neutral or 3 "
        f"protonated N-bound hydrogens, found {len(hydrogens)}. Automatic hydrogen addition "
        "or protonation normalization is not implemented."
    )


def execute_nhs_lys_amide_rdkit_graph_edit(
    *,
    protein_mol: Any,
    moiety_mol: Any,
    site_atom_index: int,
    reactive_carbon_index: int,
    leaving_atom_indices: Iterable[int],
    site_hydrogen_indices: Iterable[int] = (),
    sanitize: bool = True,
) -> RdkitGraphEditResult:
    """Execute conservative RDKit graph surgery for an NHS-Lys amide bond.

    Parameters
    ----------
    protein_mol : Any
        RDKit molecule representing the protein or protein fragment.
    moiety_mol : Any
        RDKit molecule representing the NHS ester moiety.
    site_atom_index : int
        Protein atom index for the lysine NZ atom before hydrogen removal.
    reactive_carbon_index : int
        Moiety atom index for the NHS ester acyl carbon before leaving group
        removal.
    leaving_atom_indices : Iterable[int]
        Moiety atom indices to remove as the NHS leaving group.
    site_hydrogen_indices : Iterable[int], optional
        Protein hydrogen atom indices to remove before bonding, by default
        ``()``.
    sanitize : bool, optional
        Attempt RDKit sanitization after graph surgery, by default ``True``.

    Returns
    -------
    RdkitGraphEditResult
        Product molecule without conformers, atom index maps, added bond, and
        non-fatal warnings.

    Raises
    ------
    ValueError
        If required atom indices are invalid or would be removed before bond
        creation.
    """
    from rdkit import Chem

    leaving_indices = set(leaving_atom_indices)
    hydrogen_indices = set(site_hydrogen_indices)
    validate_atom_indices(
        protein_mol,
        {site_atom_index, *hydrogen_indices},
        label="Protein",
    )
    validate_atom_indices(
        moiety_mol,
        {reactive_carbon_index, *leaving_indices},
        label="Moiety",
    )
    if site_atom_index in hydrogen_indices:
        raise ValueError("The lysine site atom cannot also be removed as a site hydrogen")
    if reactive_carbon_index in leaving_indices:
        raise ValueError("The NHS reactive carbon cannot be part of the leaving group")

    warnings: list[str] = []
    protein_atom_count = protein_mol.GetNumAtoms()
    moiety_atom_count = moiety_mol.GetNumAtoms()

    protein_rw = Chem.RWMol(protein_mol)
    for atom_index in sorted(hydrogen_indices, reverse=True):
        protein_rw.RemoveAtom(atom_index)
    protein_trimmed = protein_rw.GetMol()
    protein_old_to_new = build_old_to_new_atom_map(protein_atom_count, hydrogen_indices)
    site_atom_new = protein_old_to_new[site_atom_index]

    moiety_rw = Chem.RWMol(moiety_mol)
    for atom_index in sorted(leaving_indices, reverse=True):
        moiety_rw.RemoveAtom(atom_index)
    moiety_trimmed = moiety_rw.GetMol()
    moiety_old_to_new_local = build_old_to_new_atom_map(moiety_atom_count, leaving_indices)
    reactive_carbon_new_local = moiety_old_to_new_local[reactive_carbon_index]

    offset = protein_trimmed.GetNumAtoms()
    combined = Chem.CombineMols(protein_trimmed, moiety_trimmed)
    combined_rw = Chem.RWMol(combined)
    product_carbon_index = offset + reactive_carbon_new_local
    combined_rw.AddBond(site_atom_new, product_carbon_index, Chem.BondType.SINGLE)
    product = combined_rw.GetMol()

    if sanitize:
        try:
            Chem.SanitizeMol(product)
        except Exception as exc:  # noqa: BLE001 - preserve a conservative executor boundary
            warnings.append(f"RDKit sanitization warning: {exc}")

    if product.GetNumConformers():
        warnings.append("Removed stale RDKit conformers; product coordinates must be regenerated")
    product = clone_without_conformers(product)

    moiety_old_to_product = {
        old_index: offset + new_index for old_index, new_index in moiety_old_to_new_local.items()
    }

    return RdkitGraphEditResult(
        product_mol=product,
        removed_protein_atom_indices=tuple(sorted(hydrogen_indices)),
        removed_moiety_atom_indices=tuple(sorted(leaving_indices)),
        added_bond=AddedBond(
            begin_atom_index=site_atom_new,
            end_atom_index=product_carbon_index,
            order=1,
        ),
        protein_atom_index_map=protein_old_to_new,
        moiety_atom_index_map=moiety_old_to_product,
        warnings=tuple(warnings),
    )


def clone_without_conformers(mol: Any) -> Any:
    """Return an RDKit molecule copy with all conformers removed."""
    from rdkit import Chem

    rw_mol = Chem.RWMol(mol)
    rw_mol.RemoveAllConformers()
    return rw_mol.GetMol()


def validate_atom_indices(mol: Any, indices: set[int], *, label: str) -> None:
    """Validate that all atom indices exist in an RDKit molecule."""
    atom_count = mol.GetNumAtoms()
    invalid = sorted(index for index in indices if index < 0 or index >= atom_count)
    if invalid:
        invalid_text = ", ".join(str(index) for index in invalid)
        raise ValueError(
            f"{label} atom indices are outside the valid range 0..{atom_count - 1}: {invalid_text}"
        )


def build_old_to_new_atom_map(atom_count: int, removed_indices: set[int]) -> dict[int, int]:
    """Build an atom index map after deleting atoms."""
    old_to_new: dict[int, int] = {}
    cursor = 0
    for old_index in range(atom_count):
        if old_index in removed_indices:
            continue
        old_to_new[old_index] = cursor
        cursor += 1
    return old_to_new


def _resolve_nhs_search_indices(
    mol: Any,
    *,
    candidate_atom_indices: Iterable[int] | None,
    residue_name_prefixes: tuple[str, ...],
) -> set[int]:
    """Resolve the atom index set searched by NHS detection."""
    if candidate_atom_indices is not None:
        return set(candidate_atom_indices)

    metadata_indices: set[int] = set()
    normalized_prefixes = tuple(prefix.strip().upper() for prefix in residue_name_prefixes)
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        residue_name = info.GetResidueName().strip().upper()
        if residue_name.startswith(normalized_prefixes):
            metadata_indices.add(atom.GetIdx())

    return metadata_indices or set(range(mol.GetNumAtoms()))


def _traverse_leaving_group(mol: Any, start_index: int, *, excluded_index: int) -> set[int]:
    """Traverse the leaving group graph while excluding the reactive carbon."""
    leaving_group: set[int] = set()
    visited = {excluded_index}
    queue: deque[int] = deque([start_index])
    while queue:
        current_index = queue.popleft()
        if current_index in visited:
            continue
        visited.add(current_index)
        leaving_group.add(current_index)
        for neighbor in mol.GetAtomWithIdx(current_index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if neighbor_index not in visited:
                queue.append(neighbor_index)
    return leaving_group


def _leaving_group_has_nhs_evidence(mol: Any, leaving_group: set[int]) -> bool:
    """Return whether a leaving group has topology evidence consistent with NHS."""
    symbols = {mol.GetAtomWithIdx(index).GetSymbol() for index in leaving_group}
    return "N" in symbols and "O" in symbols and len(leaving_group) >= 3


def _deduplicate_reactive_groups(candidates: list[NhsReactiveGroup]) -> list[NhsReactiveGroup]:
    """Deduplicate reactive-group candidates by carbon, oxygen, and leaving group."""
    seen: set[tuple[int, int, tuple[int, ...]]] = set()
    unique: list[NhsReactiveGroup] = []
    for candidate in candidates:
        key = (
            candidate.reactive_carbon_index,
            candidate.bridging_oxygen_index,
            candidate.leaving_group_atom_indices,
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(candidate)
    return unique


def _atom_index(atom: Any, fallback: int) -> int:
    """Return an atom index from common topology atom representations."""
    if isinstance(atom, Mapping):
        value = atom.get("index", atom.get("molecule_atom_index", fallback))
        return int(value)
    if hasattr(atom, "molecule_atom_index"):
        return int(atom.molecule_atom_index)
    if hasattr(atom, "index"):
        return int(atom.index)
    if hasattr(atom, "GetIdx"):
        return int(atom.GetIdx())
    return fallback


def _atom_name(atom: Any) -> str | None:
    """Return an uppercase atom name from metadata or topology attributes."""
    value = _metadata_value(atom, "atom_name")
    if value is None and hasattr(atom, "name"):
        value = atom.name
    if value is None and hasattr(atom, "GetPDBResidueInfo"):
        info = atom.GetPDBResidueInfo()
        if info is not None:
            value = info.GetName()
    return str(value).strip().upper() if value is not None else None


def _chain_id(atom: Any) -> str | None:
    """Return an uppercase chain ID from metadata or PDB residue information."""
    value = _metadata_value(atom, "chain_id")
    if value is None and hasattr(atom, "GetPDBResidueInfo"):
        info = atom.GetPDBResidueInfo()
        if info is not None:
            value = info.GetChainId()
    return str(value).strip().upper() if value is not None else None


def _residue_name(atom: Any) -> str | None:
    """Return an uppercase residue name from metadata or PDB residue information."""
    value = _metadata_value(atom, "residue_name")
    if value is None and hasattr(atom, "GetPDBResidueInfo"):
        info = atom.GetPDBResidueInfo()
        if info is not None:
            value = info.GetResidueName()
    return str(value).strip().upper() if value is not None else None


def _residue_number(atom: Any) -> int | None:
    """Return a residue sequence number from metadata or PDB residue information."""
    value = _metadata_value(atom, "residue_number")
    if value is None and hasattr(atom, "GetPDBResidueInfo"):
        info = atom.GetPDBResidueInfo()
        if info is not None:
            value = info.GetResidueNumber()
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _metadata_value(atom: Any, key: str) -> Any:
    """Return a metadata value from dict or OpenFF-like atom records."""
    if isinstance(atom, Mapping):
        metadata = atom.get("metadata", {}) or {}
        return atom.get(key, metadata.get(key))
    metadata = getattr(atom, "metadata", {}) or {}
    if isinstance(metadata, Mapping):
        return metadata.get(key)
    return None


def _atomic_number(atom: Any) -> int | None:
    """Return the atomic number from common atom record types."""
    if isinstance(atom, Mapping):
        value = atom.get("atomic_number")
        return int(value) if value is not None else None
    if hasattr(atom, "atomic_number"):
        return int(atom.atomic_number)
    if hasattr(atom, "GetAtomicNum"):
        return int(atom.GetAtomicNum())
    return None


def _collect_bond_pairs(
    atom_records: Sequence[Any], *, bonds: Iterable[Any] | None
) -> set[tuple[int, int]]:
    """Collect topology bond pairs from explicit bonds or atom-local bonds."""
    bond_pairs: set[tuple[int, int]] = set()
    source_bonds: Iterable[Any]
    if bonds is not None:
        source_bonds = bonds
    else:
        gathered: list[Any] = []
        for atom in atom_records:
            gathered.extend(getattr(atom, "bonds", []) or [])
        source_bonds = gathered

    for bond in source_bonds:
        pair = _bond_pair(bond)
        if pair is not None:
            bond_pairs.add(pair)
    return bond_pairs


def _bond_pair(bond: Any) -> tuple[int, int] | None:
    """Return a sorted atom-index pair from common topology bond records."""
    if isinstance(bond, Sequence) and not isinstance(bond, str) and len(bond) == 2:
        first, second = int(bond[0]), int(bond[1])
    elif hasattr(bond, "atom1_index") and hasattr(bond, "atom2_index"):
        first, second = int(bond.atom1_index), int(bond.atom2_index)
    elif hasattr(bond, "GetBeginAtomIdx") and hasattr(bond, "GetEndAtomIdx"):
        first, second = int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())
    else:
        return None
    return tuple(sorted((first, second)))


def _neighbor_indices(atom_index: int, bond_pairs: set[tuple[int, int]]) -> set[int]:
    """Return neighbors of an atom from an undirected bond-pair set."""
    neighbors: set[int] = set()
    for first, second in bond_pairs:
        if first == atom_index:
            neighbors.add(second)
        elif second == atom_index:
            neighbors.add(first)
    return neighbors


def _position_for_index(
    atom_index: int | None,
    positions: Mapping[int, Sequence[float]] | Sequence[Sequence[float]] | None,
) -> tuple[float, float, float] | None:
    """Return an optional 3D position for an atom index."""
    if atom_index is None or positions is None:
        return None
    if isinstance(positions, Mapping):
        value = positions.get(atom_index)
    else:
        try:
            value = positions[atom_index]
        except IndexError:
            return None
    if value is None:
        return None
    if len(value) != 3:
        return None
    return (float(value[0]), float(value[1]), float(value[2]))
