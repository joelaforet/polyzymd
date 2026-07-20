"""Pure-Python PDB assembly for NHS-Lys crosslinked polymer fragments."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field, model_validator

from polyzymd.builders.conjugation.structure.parsing import (
    ATOM_RECORD_PREFIXES as _ATOM_RECORD_PREFIXES,
)
from polyzymd.builders.conjugation.structure.parsing import (
    PdbAtomRecord,
    parse_pdb_atom_lines,
    parse_pdb_atom_records,
)

_MAX_CONECT_TARGETS = 4

_POC_RESIDUE_NAME_MAP = {
    "SB1": "SBM",
    "SB2": "SBM",
    "SBMA": "SBM",
    "EG1": "EGP",
    "EG2": "EGP",
    "EGPMA": "EGP",
    "EGM": "EGM",
    "EGP": "EGP",
    "NH1": "NHS",
    "NH2": "NHS",
    "NHS": "NHS",
    "NHX": "NHX",
}


class AtomIdentity(BaseModel):
    """PDB atom identity associated with one RDKit atom index."""

    chain_id: str = Field("", max_length=1)
    residue_name: str = Field(..., min_length=1)
    residue_number: int
    insertion_code: str = Field("", max_length=1)
    atom_name: str = Field(..., min_length=1)
    atom_index: int = Field(..., ge=0)
    serial: int | None = Field(None, ge=1)
    element: str = Field(..., min_length=1)
    alt_loc: str = Field("", max_length=1)
    rdkit_index: int = Field(..., ge=0)

    def same_residue_as(self, atom: PdbAtomRecord) -> bool:
        """Return whether this identity belongs to the same PDB residue."""
        return (
            self.chain_id.upper() == atom.chain_id.upper()
            and self.residue_name.upper() == atom.residue_name.upper()
            and self.residue_number == atom.residue_number
            and self.insertion_code.upper() == atom.insertion_code.upper()
        )


class RdkitInputBundle(BaseModel):
    """RDKit molecule plus PDB identity lookups for each atom."""

    model_config = {"arbitrary_types_allowed": True}

    mol: Any
    source_path: Path
    source_kind: str = "pdb"
    atom_identities: tuple[AtomIdentity, ...]
    pdb_index_to_rdkit_index: dict[int, int]
    rdkit_index_to_pdb_index: dict[int, int]
    serial_to_rdkit_index: dict[int, int]

    def identity_for_rdkit_index(self, rdkit_index: int) -> AtomIdentity:
        """Return identity metadata for an RDKit atom index."""
        pdb_index = self.rdkit_index_to_pdb_index[rdkit_index]
        return self.atom_identities[pdb_index]


def load_pdb_as_rdkit_input(
    path: Path | str,
    *,
    sanitize: bool = False,
    proximity_bonding: bool = True,
) -> RdkitInputBundle:
    """Load a PDB file into RDKit while preserving PDB atom identity."""
    from rdkit import Chem

    pdb_path = Path(path)
    pdb_atoms = parse_pdb_atom_records(pdb_path)
    mol = Chem.MolFromPDBFile(
        str(pdb_path),
        sanitize=sanitize,
        removeHs=False,
        proximityBonding=proximity_bonding,
    )
    if mol is None:
        raise ValueError(f"RDKit could not load PDB file {pdb_path}")
    return build_rdkit_input_bundle(mol, pdb_atoms=pdb_atoms, source_path=pdb_path)


def build_rdkit_input_bundle(
    mol: Any,
    *,
    pdb_atoms: tuple[PdbAtomRecord, ...],
    source_path: Path | str,
    source_kind: str = "pdb",
) -> RdkitInputBundle:
    """Build identity maps for an existing RDKit molecule and PDB records."""
    atom_count = mol.GetNumAtoms()
    if atom_count != len(pdb_atoms):
        raise ValueError(
            "RDKit/PDB atom-count mismatch while loading identity map: "
            f"RDKit has {atom_count} atoms but PDB parser found {len(pdb_atoms)} atoms"
        )

    identities: list[AtomIdentity] = []
    serial_to_rdkit_index: dict[int, int] = {}
    for rdkit_index, pdb_atom in enumerate(pdb_atoms):
        rdkit_atom = mol.GetAtomWithIdx(rdkit_index)
        _validate_rdkit_pdb_identity(rdkit_atom, pdb_atom, rdkit_index)
        identity = _identity_from_pdb_atom(pdb_atom, rdkit_index=rdkit_index, rdkit_atom=rdkit_atom)
        identities.append(identity)
        if identity.serial is not None:
            if identity.serial in serial_to_rdkit_index:
                raise ValueError(f"Duplicate PDB atom serial {identity.serial} in {source_path}")
            serial_to_rdkit_index[identity.serial] = rdkit_index

    pdb_index_to_rdkit_index = {
        identity.atom_index: identity.rdkit_index for identity in identities
    }
    rdkit_index_to_pdb_index = {
        identity.rdkit_index: identity.atom_index for identity in identities
    }
    if sorted(pdb_index_to_rdkit_index) != list(range(len(identities))):
        raise ValueError("PDB atom indices must be contiguous and zero-based for RDKit mapping")

    return RdkitInputBundle(
        mol=mol,
        source_path=Path(source_path),
        source_kind=source_kind,
        atom_identities=tuple(identities),
        pdb_index_to_rdkit_index=pdb_index_to_rdkit_index,
        rdkit_index_to_pdb_index=rdkit_index_to_pdb_index,
        serial_to_rdkit_index=serial_to_rdkit_index,
    )


def _identity_from_pdb_atom(
    pdb_atom: PdbAtomRecord,
    *,
    rdkit_index: int,
    rdkit_atom: Any,
) -> AtomIdentity:
    """Create an atom identity after validating required metadata."""
    atom_index = pdb_atom.atom_index
    if atom_index is None:
        raise ValueError(f"PDB atom {pdb_atom.atom_name} is missing zero-based atom_index")
    element = (pdb_atom.element or rdkit_atom.GetSymbol()).strip()
    if not element:
        raise ValueError(
            f"PDB atom index {atom_index} serial {pdb_atom.serial} is missing element metadata"
        )
    return AtomIdentity(
        chain_id=pdb_atom.chain_id.strip(),
        residue_name=pdb_atom.residue_name.strip().upper(),
        residue_number=pdb_atom.residue_number,
        insertion_code=pdb_atom.insertion_code.strip().upper(),
        atom_name=pdb_atom.atom_name.strip().upper(),
        atom_index=atom_index,
        serial=pdb_atom.serial,
        element=element.upper(),
        alt_loc=pdb_atom.alt_loc.strip().upper(),
        rdkit_index=rdkit_index,
    )


def _validate_rdkit_pdb_identity(
    rdkit_atom: Any, pdb_atom: PdbAtomRecord, rdkit_index: int
) -> None:
    """Validate RDKit PDB residue metadata against parsed PDB order."""
    info = rdkit_atom.GetPDBResidueInfo()
    if info is None:
        return

    mismatches: list[str] = []
    _compare_text(mismatches, "atom name", info.GetName(), pdb_atom.atom_name)
    _compare_text(mismatches, "residue name", info.GetResidueName(), pdb_atom.residue_name)
    _compare_text(mismatches, "chain ID", info.GetChainId(), pdb_atom.chain_id)
    _compare_text(mismatches, "insertion code", info.GetInsertionCode(), pdb_atom.insertion_code)
    _compare_text(mismatches, "alt loc", info.GetAltLoc(), pdb_atom.alt_loc)
    if int(info.GetResidueNumber()) != pdb_atom.residue_number:
        mismatches.append(
            f"residue number RDKit={info.GetResidueNumber()} PDB={pdb_atom.residue_number}"
        )
    serial = int(info.GetSerialNumber())
    if pdb_atom.serial is not None and serial not in {0, pdb_atom.serial}:
        mismatches.append(f"serial RDKit={serial} PDB={pdb_atom.serial}")
    if mismatches:
        detail = "; ".join(mismatches)
        raise ValueError(f"RDKit/PDB atom-order mismatch at atom {rdkit_index}: {detail}")


def _compare_text(mismatches: list[str], label: str, rdkit_value: str, pdb_value: str) -> None:
    """Append a mismatch if normalized RDKit and PDB fields differ."""
    left = str(rdkit_value).strip().upper()
    right = str(pdb_value).strip().upper()
    if left != right:
        mismatches.append(f"{label} RDKit={left!r} PDB={right!r}")


class PlacedPolymerFragment(BaseModel):
    """Placed polymer fragment atoms and connectivity for PDB assembly."""

    atoms: tuple[PdbAtomRecord, ...] = Field(..., min_length=1)
    bonds: tuple[tuple[int | str, int | str], ...] = Field(default_factory=tuple)
    bond_orders: tuple[tuple[int | str, int | str, float], ...] = Field(default_factory=tuple)
    leaving_atom_serials: tuple[int, ...] = Field(default_factory=tuple)
    leaving_atom_indices: tuple[int, ...] = Field(default_factory=tuple)
    leaving_atom_names: tuple[str, ...] = Field(default_factory=tuple)
    reactive_atom_serial: int | None = None
    reactive_atom_index: int | None = None
    reactive_atom_name: str | None = None
    name: str = "polymer"

    @classmethod
    def from_pdb_lines(
        cls,
        lines: Iterable[str],
        *,
        bonds: Sequence[tuple[int | str, int | str]] = (),
        bond_orders: Sequence[tuple[int | str, int | str, float]] = (),
        leaving_atom_serials: Sequence[int] = (),
        leaving_atom_indices: Sequence[int] = (),
        leaving_atom_names: Sequence[str] = (),
        reactive_atom_serial: int | None = None,
        reactive_atom_index: int | None = None,
        reactive_atom_name: str | None = None,
        name: str = "polymer",
    ) -> PlacedPolymerFragment:
        """Build a placed polymer fragment from fixed-width PDB atom lines.

        Parameters
        ----------
        lines : iterable of str
            PDB lines containing ATOM/HETATM records.
        bonds : sequence of tuple, optional
            Internal polymer bonds using atom serials, indices, or names, by default ``()``.
        bond_orders : sequence of tuple, optional
            Internal polymer bonds with bond-order metadata, by default ``()``.
        leaving_atom_serials : sequence of int, optional
            Atom serials to omit as leaving-group atoms, by default ``()``.
        leaving_atom_indices : sequence of int, optional
            Zero-based atom indices to omit as leaving-group atoms, by default ``()``.
        leaving_atom_names : sequence of str, optional
            Atom names to omit as leaving-group atoms, by default ``()``.
        reactive_atom_serial : int or None, optional
            Reactive acyl-carbon atom serial, by default ``None``.
        reactive_atom_index : int or None, optional
            Reactive acyl-carbon atom index, by default ``None``.
        reactive_atom_name : str or None, optional
            Reactive acyl-carbon atom name, by default ``None``.
        name : str, optional
            Fragment label for diagnostics, by default ``"polymer"``.

        Returns
        -------
        PlacedPolymerFragment
            Parsed placed polymer fragment.
        """
        atoms = parse_pdb_atom_lines(lines)
        return cls(
            atoms=atoms,
            bonds=tuple(bonds),
            bond_orders=tuple(bond_orders),
            leaving_atom_serials=tuple(leaving_atom_serials),
            leaving_atom_indices=tuple(leaving_atom_indices),
            leaving_atom_names=tuple(leaving_atom_names),
            reactive_atom_serial=reactive_atom_serial,
            reactive_atom_index=reactive_atom_index,
            reactive_atom_name=reactive_atom_name,
            name=name,
        )

    @model_validator(mode="after")
    def validate_reactive_reference(self) -> PlacedPolymerFragment:
        """Validate that a reactive atom selector was provided.

        Returns
        -------
        PlacedPolymerFragment
            Validated model.

        Raises
        ------
        ValueError
            If no reactive atom selector is present.
        """
        if (
            self.reactive_atom_serial is None
            and self.reactive_atom_index is None
            and self.reactive_atom_name is None
        ):
            raise ValueError("Placed polymer fragments require a reactive atom selector")
        return self


class NhsLysPdbAttachment(BaseModel):
    """Explicit PDB-level NHS-Lys attachment data for assembly."""

    target_chain: str = Field("A", max_length=1)
    target_residue_name: str | None = None
    target_residue_number: int
    target_insertion_code: str = Field("", max_length=1)
    target_atom_name: str = "NZ"
    nz_hydrogen_atom_serials_to_remove: tuple[int, ...] = Field(default_factory=tuple)
    nz_hydrogen_atom_indices_to_remove: tuple[int, ...] = Field(default_factory=tuple)
    nz_hydrogen_atom_names_to_remove: tuple[str, ...] = Field(default_factory=tuple)
    max_nz_hydrogens_to_remove: int = Field(2, ge=0)
    lysine_target_resname: str = "LYX"
    polymer_target_resname: str = "NHX"


class PdbLinkageAttachment(BaseModel):
    """Generic explicit PDB linkage attachment data for assembly."""

    target_chain: str = Field("A", max_length=1)
    target_residue_name: str | None = None
    target_residue_number: int
    target_insertion_code: str = Field("", max_length=1)
    target_atom_name: str
    protein_leaving_atoms_to_remove: tuple[PdbAtomRecord, ...] = Field(default_factory=tuple)
    protein_target_resname: str
    modifier_target_resname: str


PdbAssemblyAttachment = NhsLysPdbAttachment | PdbLinkageAttachment


class CrosslinkedPdbAssemblyOptions(BaseModel):
    """Options controlling deterministic crosslinked PDB output."""

    protein_chain: str = Field("A", max_length=1)
    polymer_chain: str = Field("C", max_length=1)
    preserve_serials: bool = False
    include_link_records: bool = False
    append_ter_records: bool = True


class CrosslinkedPdbAssemblyResult(BaseModel):
    """Summary of a crosslinked PDB assembly artifact."""

    output_path: Path
    protein_atom_count: int
    polymer_atom_count: int
    removed_protein_atom_count: int
    removed_polymer_atom_count: int
    removed_atom_serials: tuple[int, ...]
    removed_atom_names: tuple[str, ...]
    added_conect_pair: tuple[int, int]
    added_conect_pairs: tuple[tuple[int, int], ...] = Field(default_factory=tuple)
    warnings: tuple[str, ...] = Field(default_factory=tuple)
    residue_mappings: dict[str, dict[str, int | str]] = Field(default_factory=dict)
    atom_mappings: dict[str, dict[str, int | str]] = Field(default_factory=dict)
    attachment_endpoint_records: tuple[dict[str, int | str | dict[str, int | str]], ...] = Field(
        default_factory=tuple
    )


def canonicalize_poc_residue_name(raw_residue_name: str, *, crosslinked: bool = False) -> str:
    """Canonicalize common POC polymer residue names for PDB output.

    Parameters
    ----------
    raw_residue_name : str
        Input residue name.
    crosslinked : bool, optional
        Whether this residue is the linked NHS product residue, by default ``False``.

    Returns
    -------
    str
        Three-character PDB-safe residue name.
    """
    normalized = (raw_residue_name or "POL").strip().upper()
    if crosslinked and normalized in {"NH1", "NH2", "NHS", "NHX"}:
        return "NHX"
    canonical = _POC_RESIDUE_NAME_MAP.get(normalized, normalized)
    return _pdb_safe_residue_name(canonical)


def write_crosslinked_pdb(
    protein_pdb_path: Path | str,
    polymer_fragment: PlacedPolymerFragment | Sequence[PlacedPolymerFragment],
    attachment: PdbAssemblyAttachment | Sequence[PdbAssemblyAttachment],
    output_path: Path | str,
    options: CrosslinkedPdbAssemblyOptions | None = None,
) -> CrosslinkedPdbAssemblyResult:
    """Write a deterministic crosslinked protein-polymer PDB with CONECT records.

    Parameters
    ----------
    protein_pdb_path : Path or str
        Source protein PDB. The file is read but never modified.
    polymer_fragment : PlacedPolymerFragment or sequence of PlacedPolymerFragment
        Already placed polymer fragment or fragments.
    attachment : NhsLysPdbAttachment or PdbLinkageAttachment
        Explicit attachment selectors. The NHS-Lys model preserves the legacy
        hydrogen-only behavior while the generic model removes resolved protein
        leaving atoms exactly.
    output_path : Path or str
        Destination PDB path.
    options : CrosslinkedPdbAssemblyOptions or None, optional
        Writer options, by default ``None``.

    Returns
    -------
    CrosslinkedPdbAssemblyResult
        Summary of written atoms, removed atoms, connectivity, and residue mappings.

    Raises
    ------
    ValueError
        If the attachment atoms cannot be resolved.
    """
    writer_options = options or CrosslinkedPdbAssemblyOptions()
    protein_path = Path(protein_pdb_path)
    destination = Path(output_path)
    fragments = _as_fragment_tuple(polymer_fragment)
    attachments = _as_attachment_tuple(attachment)
    if len(attachments) == 1 and len(fragments) > 1:
        attachments = attachments * len(fragments)
    if len(attachments) != len(fragments):
        raise ValueError(
            "Crosslinked PDB assembly requires one attachment per fragment "
            f"(attachments={len(attachments)}, fragments={len(fragments)})"
        )
    warnings: list[str] = []

    protein_atoms = _parse_pdb_atoms(protein_path)
    kept_protein_atoms, removed_protein_atoms = _prepare_protein_atoms(
        protein_atoms,
        attachments,
        writer_options,
        warnings,
    )

    output_atoms: list[PdbAtomRecord] = []
    atom_lines: list[str] = []
    conect_map: dict[int, set[int]] = {}
    residue_mappings: dict[str, dict[str, int | str]] = {}
    atom_mappings: dict[str, dict[str, int | str]] = {}
    removed_polymer_atoms: list[PdbAtomRecord] = []
    crosslink_pair: tuple[int, int] | None = None
    crosslink_pairs: list[tuple[int, int]] = []

    next_serial = 1
    used_serials: set[int] = set()
    protein_serial_by_index: dict[int, int] = {}
    protein_serial_by_input_serial: dict[int, int] = {}

    for atom in kept_protein_atoms:
        new_serial, next_serial = _next_atom_serial(atom, next_serial, used_serials, writer_options)
        updated = atom.model_copy(
            update={"serial": new_serial, "chain_id": writer_options.protein_chain}
        )
        output_atoms.append(updated)
        atom_lines.append(_format_pdb_atom_line(updated))
        if atom.atom_index is not None:
            protein_serial_by_index[atom.atom_index] = new_serial
        if atom.serial is not None:
            protein_serial_by_input_serial[atom.serial] = new_serial

    if writer_options.append_ter_records and output_atoms:
        atom_lines.append(_format_ter_line(output_atoms[-1], next_serial))
        next_serial += 1

    next_polymer_residue_number = 1
    for fragment_index, (fragment, fragment_attachment) in enumerate(
        zip(fragments, attachments, strict=True), start=1
    ):
        target_atom = _resolve_attachment_atom(kept_protein_atoms, fragment_attachment)
        target_serial = _output_serial_for_atom(
            target_atom,
            serial_by_index=protein_serial_by_index,
            serial_by_input_serial=protein_serial_by_input_serial,
        )
        fragment_result = _append_polymer_fragment(
            fragment,
            fragment_index=fragment_index,
            attachment=fragment_attachment,
            options=writer_options,
            starting_residue_number=next_polymer_residue_number,
            next_serial=next_serial,
            used_serials=used_serials,
            conect_map=conect_map,
            warnings=warnings,
        )
        next_serial = fragment_result.next_serial
        next_polymer_residue_number = fragment_result.next_residue_number
        atom_lines.extend(fragment_result.lines)
        output_atoms.extend(fragment_result.atoms)
        removed_polymer_atoms.extend(fragment_result.removed_atoms)
        residue_mappings.update(fragment_result.residue_mappings)
        atom_mappings.update(fragment_result.atom_mappings)

        conect_map.setdefault(target_serial, set()).add(fragment_result.reactive_serial)
        conect_map.setdefault(fragment_result.reactive_serial, set()).add(target_serial)
        crosslink_pair = (target_serial, fragment_result.reactive_serial)
        crosslink_pairs.append(crosslink_pair)

        if writer_options.append_ter_records and fragment_result.atoms:
            atom_lines.append(_format_ter_line(fragment_result.atoms[-1], next_serial))
            next_serial += 1

    if crosslink_pair is None:
        raise ValueError("No polymer fragment was available for crosslinked PDB assembly")

    link_lines = []
    if writer_options.include_link_records:
        for fragment_attachment, pair in zip(attachments, crosslink_pairs, strict=True):
            target_atom = _resolve_attachment_atom(kept_protein_atoms, fragment_attachment)
            link_lines.append(_format_link_line(target_atom, pair[1], output_atoms))

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        handle.writelines(link_lines)
        handle.writelines(atom_lines)
        for serial in sorted(conect_map):
            bonded_serials = sorted(conect_map[serial])
            handle.write(_format_conect_records(serial, bonded_serials))
        handle.write("END\n")

    removed_atoms = [*removed_protein_atoms, *removed_polymer_atoms]
    endpoint_records = tuple(
        _attachment_endpoint_record(
            fragment_index=index,
            protein_atom=_product_atom_by_serial(output_atoms, pair[0]),
            modifier_atom=_product_atom_by_serial(output_atoms, pair[1]),
            conect_pair=pair,
        )
        for index, pair in enumerate(crosslink_pairs, start=1)
    )
    return CrosslinkedPdbAssemblyResult(
        output_path=destination,
        protein_atom_count=len(kept_protein_atoms),
        polymer_atom_count=sum(
            1 for atom in output_atoms if atom.chain_id == writer_options.polymer_chain
        ),
        removed_protein_atom_count=len(removed_protein_atoms),
        removed_polymer_atom_count=len(removed_polymer_atoms),
        removed_atom_serials=tuple(
            atom.serial for atom in removed_atoms if atom.serial is not None
        ),
        removed_atom_names=tuple(atom.atom_name for atom in removed_atoms),
        added_conect_pair=crosslink_pair,
        added_conect_pairs=tuple(crosslink_pairs),
        warnings=tuple(warnings),
        residue_mappings=residue_mappings,
        atom_mappings=atom_mappings,
        attachment_endpoint_records=endpoint_records,
    )


class _PolymerAppendResult(BaseModel):
    """Internal result from appending one polymer fragment."""

    atoms: list[PdbAtomRecord]
    lines: list[str]
    removed_atoms: list[PdbAtomRecord]
    reactive_serial: int
    next_serial: int
    next_residue_number: int
    residue_mappings: dict[str, dict[str, int | str]]
    atom_mappings: dict[str, dict[str, int | str]]


def _append_polymer_fragment(
    fragment: PlacedPolymerFragment,
    *,
    fragment_index: int,
    attachment: PdbAssemblyAttachment,
    options: CrosslinkedPdbAssemblyOptions,
    starting_residue_number: int,
    next_serial: int,
    used_serials: set[int],
    conect_map: dict[int, set[int]],
    warnings: list[str],
) -> _PolymerAppendResult:
    """Append retained polymer atoms and internal connectivity for one fragment."""
    _validate_polymer_leaving_selectors(fragment)
    removed_atoms = [atom for atom in fragment.atoms if _is_removed_polymer_atom(atom, fragment)]
    removed_keys = {_atom_identity(atom) for atom in removed_atoms}
    kept_atoms = [atom for atom in fragment.atoms if _atom_identity(atom) not in removed_keys]
    reactive_atom = _resolve_reactive_polymer_atom(fragment)
    if _atom_identity(reactive_atom) in removed_keys:
        raise ValueError("The polymer reactive atom cannot be listed as a leaving-group atom")
    kept_atoms = _order_polymer_atoms_by_connectivity(
        fragment,
        kept_atoms,
        warnings,
    )

    residue_key_to_number: dict[tuple[int, str], int] = {}
    residue_mappings: dict[str, dict[str, int | str]] = {}
    atom_mappings: dict[str, dict[str, int | str]] = {}
    reactive_residue_key = _polymer_residue_key(reactive_atom)
    atom_serial_by_index: dict[int, int] = {}
    atom_serial_by_input_serial: dict[int, int] = {}
    atom_serial_by_name: dict[str, int] = {}
    lines: list[str] = []
    output_atoms: list[PdbAtomRecord] = []
    residue_cursor = starting_residue_number

    for atom in kept_atoms:
        residue_key = _polymer_residue_key(atom)
        if residue_key not in residue_key_to_number:
            residue_key_to_number[residue_key] = residue_cursor
            residue_mappings[f"fragment_{fragment_index}:{residue_key[0]}{residue_key[1]}"] = {
                "source_residue_number": residue_key[0],
                "source_insertion_code": residue_key[1] or "",
                "target_residue_number": residue_cursor,
                "target_insertion_code": "",
                "target_chain": options.polymer_chain,
            }
            residue_cursor += 1

        is_crosslinked_residue = residue_key == reactive_residue_key
        residue_name = (
            _pdb_safe_residue_name(_attachment_modifier_product_resname(attachment))
            if is_crosslinked_residue
            else _non_crosslinked_modifier_residue_name(atom, attachment)
        )
        new_serial, next_serial = _next_atom_serial(atom, next_serial, used_serials, options)
        updated = atom.model_copy(
            update={
                "serial": new_serial,
                "chain_id": options.polymer_chain,
                "residue_number": residue_key_to_number[residue_key],
                "residue_name": residue_name,
                "record_name": "HETATM",
            }
        )
        output_atoms.append(updated)
        lines.append(_format_pdb_atom_line(updated))
        if atom.atom_index is not None:
            atom_serial_by_index[atom.atom_index] = new_serial
        if atom.serial is not None:
            atom_serial_by_input_serial[atom.serial] = new_serial
        atom_serial_by_name.setdefault(atom.atom_name, new_serial)
        atom_mappings[
            f"fragment_{fragment_index}:{residue_key[0]}{residue_key[1]}:{atom.atom_name}"
        ] = {
            "source_serial": int(atom.serial) if atom.serial is not None else "",
            "source_atom_index": int(atom.atom_index) if atom.atom_index is not None else "",
            "source_atom_name": atom.atom_name,
            "source_residue_name": atom.residue_name,
            "source_residue_number": residue_key[0],
            "source_insertion_code": residue_key[1] or "",
            "target_serial": new_serial,
            "target_atom_name": updated.atom_name,
            "target_residue_name": updated.residue_name,
            "target_residue_number": updated.residue_number,
            "target_insertion_code": updated.insertion_code or "",
            "target_chain": updated.chain_id,
        }

    retained_keys = {_atom_identity(atom) for atom in kept_atoms}
    for atom_1_ref, atom_2_ref in fragment.bonds:
        atom_1 = _resolve_polymer_ref(fragment, atom_1_ref)
        atom_2 = _resolve_polymer_ref(fragment, atom_2_ref)
        if atom_1 is None or atom_2 is None:
            warnings.append(f"Skipped unresolved polymer bond reference in {fragment.name}")
            continue
        if (
            _atom_identity(atom_1) not in retained_keys
            or _atom_identity(atom_2) not in retained_keys
        ):
            warnings.append(
                f"Skipped polymer bond to omitted leaving-group atom in {fragment.name}"
            )
            continue
        serial_1 = _output_serial_for_atom(
            atom_1,
            serial_by_index=atom_serial_by_index,
            serial_by_input_serial=atom_serial_by_input_serial,
            serial_by_name=atom_serial_by_name,
        )
        serial_2 = _output_serial_for_atom(
            atom_2,
            serial_by_index=atom_serial_by_index,
            serial_by_input_serial=atom_serial_by_input_serial,
            serial_by_name=atom_serial_by_name,
        )
        conect_map.setdefault(serial_1, set()).add(serial_2)
        conect_map.setdefault(serial_2, set()).add(serial_1)

    reactive_serial = _output_serial_for_atom(
        reactive_atom,
        serial_by_index=atom_serial_by_index,
        serial_by_input_serial=atom_serial_by_input_serial,
        serial_by_name=atom_serial_by_name,
    )
    return _PolymerAppendResult(
        atoms=output_atoms,
        lines=lines,
        removed_atoms=removed_atoms,
        reactive_serial=reactive_serial,
        next_serial=next_serial,
        next_residue_number=residue_cursor,
        residue_mappings=residue_mappings,
        atom_mappings=atom_mappings,
    )


def _product_atom_by_serial(atoms: list[PdbAtomRecord], serial: int) -> PdbAtomRecord:
    """Return an emitted atom by exact output serial."""
    for atom in atoms:
        if atom.serial == serial:
            return atom
    raise ValueError(f"Crosslinked PDB assembly recorded missing output serial {serial}")


def _attachment_endpoint_record(
    *,
    fragment_index: int,
    protein_atom: PdbAtomRecord,
    modifier_atom: PdbAtomRecord,
    conect_pair: tuple[int, int],
) -> dict[str, int | str | dict[str, int | str]]:
    """Build serial-first endpoint provenance for one emitted attachment."""
    return {
        "attachment_index": fragment_index,
        "conect_pair": {"protein_serial": conect_pair[0], "modifier_serial": conect_pair[1]},
        "protein_endpoint": _endpoint_atom_payload(protein_atom),
        "modifier_endpoint": _endpoint_atom_payload(modifier_atom),
    }


def _endpoint_atom_payload(atom: PdbAtomRecord) -> dict[str, int | str]:
    """Return a JSON-safe atom endpoint record."""
    return {
        "serial": int(atom.serial) if atom.serial is not None else "",
        "atom_name": atom.atom_name,
        "residue_name": atom.residue_name,
        "chain_id": atom.chain_id,
        "residue_number": atom.residue_number,
        "insertion_code": atom.insertion_code or "",
    }


def _validate_polymer_leaving_selectors(fragment: PlacedPolymerFragment) -> None:
    """Validate name-only polymer leaving selectors before atom removal."""
    if fragment.leaving_atom_serials or fragment.leaving_atom_indices:
        return
    for name in fragment.leaving_atom_names:
        matches = [atom for atom in fragment.atoms if atom.atom_name.upper() == name.upper()]
        unique = _unique_atoms(matches)
        if len(unique) != 1:
            raise ValueError(
                f"Expected exactly one leaving atom named {name} in {fragment.name}, "
                f"found {len(unique)}"
            )


def _order_polymer_atoms_by_connectivity(
    fragment: PlacedPolymerFragment,
    kept_atoms: list[PdbAtomRecord],
    warnings: list[str],
) -> list[PdbAtomRecord]:
    """Return atoms grouped by linear polymer connectivity where available."""
    residue_order = list(dict.fromkeys(_polymer_residue_key(atom) for atom in kept_atoms))
    if len(residue_order) <= 2:
        return kept_atoms

    residue_rank = {key: index for index, key in enumerate(residue_order)}
    retained_keys = {_atom_identity(atom) for atom in kept_atoms}
    adjacency: dict[tuple[int, str], set[tuple[int, str]]] = defaultdict(set)
    for atom_1_ref, atom_2_ref in fragment.bonds:
        atom_1 = _resolve_polymer_ref(fragment, atom_1_ref)
        atom_2 = _resolve_polymer_ref(fragment, atom_2_ref)
        if atom_1 is None or atom_2 is None:
            continue
        if (
            _atom_identity(atom_1) not in retained_keys
            or _atom_identity(atom_2) not in retained_keys
        ):
            continue
        key_1 = _polymer_residue_key(atom_1)
        key_2 = _polymer_residue_key(atom_2)
        if key_1 == key_2:
            continue
        adjacency[key_1].add(key_2)
        adjacency[key_2].add(key_1)

    if not adjacency:
        return kept_atoms

    ordered_residues: list[tuple[int, str]] = []
    visited: set[tuple[int, str]] = set()
    for seed in residue_order:
        if seed in visited:
            continue
        if seed not in adjacency:
            ordered_residues.append(seed)
            visited.add(seed)
            continue
        component = _polymer_residue_component(seed, adjacency)
        visited.update(component)
        if any(len(adjacency[key]) > 2 for key in component):
            warnings.append(
                "Branched polymer residue connectivity was preserved in source PDB order"
            )
            ordered_residues.extend(
                key for key in residue_order if key in component and key not in ordered_residues
            )
            continue
        ordered_residues.extend(_linear_polymer_residue_order(component, adjacency, residue_rank))

    groups: dict[tuple[int, str], list[PdbAtomRecord]] = defaultdict(list)
    for atom in kept_atoms:
        groups[_polymer_residue_key(atom)].append(atom)
    return [atom for key in ordered_residues for atom in groups.get(key, ())]


def _polymer_residue_component(
    seed: tuple[int, str],
    adjacency: dict[tuple[int, str], set[tuple[int, str]]],
) -> set[tuple[int, str]]:
    """Return one connected residue component."""
    component: set[tuple[int, str]] = set()
    stack = [seed]
    while stack:
        key = stack.pop()
        if key in component:
            continue
        component.add(key)
        stack.extend(neighbor for neighbor in adjacency[key] if neighbor not in component)
    return component


def _linear_polymer_residue_order(
    component: set[tuple[int, str]],
    adjacency: dict[tuple[int, str], set[tuple[int, str]]],
    residue_rank: dict[tuple[int, str], int],
) -> list[tuple[int, str]]:
    """Walk a path-like residue component from the earliest source-order endpoint."""
    endpoints = [key for key in component if len(adjacency[key]) <= 1]
    start_candidates = endpoints or list(component)
    current = min(start_candidates, key=lambda key: residue_rank.get(key, len(residue_rank)))
    previous: tuple[int, str] | None = None
    ordered: list[tuple[int, str]] = []
    while current not in ordered:
        ordered.append(current)
        next_entries = [neighbor for neighbor in adjacency[current] if neighbor != previous]
        if not next_entries:
            break
        previous, current = current, min(
            next_entries,
            key=lambda key: residue_rank.get(key, len(residue_rank)),
        )
    ordered.extend(key for key in component if key not in ordered)
    return ordered


def _parse_pdb_atoms(path: Path) -> list[PdbAtomRecord]:
    """Parse all ATOM/HETATM records from a PDB file."""
    atoms = list(parse_pdb_atom_records(path, require_atoms=True))
    return atoms


def _prepare_protein_atoms(
    atoms: Sequence[PdbAtomRecord],
    attachment: PdbAssemblyAttachment | Sequence[PdbAssemblyAttachment],
    options: CrosslinkedPdbAssemblyOptions,
    warnings: list[str],
) -> tuple[list[PdbAtomRecord], list[PdbAtomRecord]]:
    """Rename the linked residue and remove selected leaving atoms."""
    attachments = _as_attachment_tuple(attachment)
    atom_indices_to_remove: set[int] = set()
    for item in attachments:
        atom_indices_to_remove.update(_select_protein_atoms_to_remove(atoms, item, warnings))
    kept_atoms: list[PdbAtomRecord] = []
    removed_atoms: list[PdbAtomRecord] = []
    for atom in atoms:
        if atom.atom_index in atom_indices_to_remove:
            removed_atoms.append(atom)
            continue
        update = {"chain_id": options.protein_chain}
        for item in attachments:
            if _matches_attachment_residue(atom, item):
                update["residue_name"] = _pdb_safe_residue_name(
                    _attachment_protein_product_resname(item)
                )
                break
        kept_atoms.append(atom.model_copy(update=update))
    return kept_atoms, removed_atoms


def _select_protein_atoms_to_remove(
    atoms: Sequence[PdbAtomRecord],
    attachment: PdbAssemblyAttachment,
    warnings: list[str],
) -> set[int]:
    """Select protein atom indices to remove for the attachment model."""
    if isinstance(attachment, NhsLysPdbAttachment):
        return _select_nz_hydrogens(atoms, attachment, warnings)
    return _select_generic_protein_leaving_atoms(atoms, attachment)


def _select_nz_hydrogens(
    atoms: Sequence[PdbAtomRecord],
    attachment: NhsLysPdbAttachment,
    warnings: list[str],
) -> set[int]:
    """Select at most two linked lysine NZ hydrogens for removal."""
    name_selectors = {name.upper() for name in attachment.nz_hydrogen_atom_names_to_remove}
    serial_selectors = set(attachment.nz_hydrogen_atom_serials_to_remove)
    index_selectors = set(attachment.nz_hydrogen_atom_indices_to_remove)
    residue_hydrogens = [
        atom
        for atom in atoms
        if _matches_attachment_residue(atom, attachment) and _is_hydrogen_atom(atom)
    ]

    if name_selectors or serial_selectors or index_selectors:
        selected = [
            atom
            for atom in residue_hydrogens
            if (atom.serial in serial_selectors)
            or (atom.atom_index in index_selectors)
            or (atom.atom_name.upper() in name_selectors)
        ]
    else:
        selected = sorted(residue_hydrogens, key=lambda atom: atom.atom_name)[1:]

    selected = sorted(
        selected, key=lambda atom: atom.atom_index if atom.atom_index is not None else -1
    )
    if len(selected) > attachment.max_nz_hydrogens_to_remove:
        warnings.append(
            "More than two linked lysine hydrogens were selected; only the first two were removed"
        )
        selected = selected[: attachment.max_nz_hydrogens_to_remove]
    return {atom.atom_index for atom in selected if atom.atom_index is not None}


def _select_generic_protein_leaving_atoms(
    atoms: Sequence[PdbAtomRecord], attachment: PdbLinkageAttachment
) -> set[int]:
    """Select exact resolved generic protein leaving atoms for removal."""
    selected_indices: set[int] = set()
    for leaving_atom in attachment.protein_leaving_atoms_to_remove:
        matches = [
            atom for atom in atoms if _matches_generic_leaving_atom(atom, leaving_atom, attachment)
        ]
        if len(matches) != 1:
            raise ValueError(
                "Expected exactly one resolved protein leaving atom during PDB assembly "
                f"for {leaving_atom.chain_id}:{leaving_atom.residue_name}:"
                f"{leaving_atom.residue_number}{leaving_atom.insertion_code}:"
                f"{leaving_atom.atom_name}, found {len(matches)}"
            )
        matched_atom = matches[0]
        if matched_atom.atom_index is None:
            raise ValueError(
                "Resolved protein leaving atom lacks an atom index during PDB assembly: "
                f"{matched_atom.atom_name}"
            )
        selected_indices.add(matched_atom.atom_index)
    return selected_indices


def _resolve_attachment_atom(
    atoms: Sequence[PdbAtomRecord], attachment: PdbAssemblyAttachment
) -> PdbAtomRecord:
    """Resolve the linked protein atom after protein rewriting."""
    matches = [
        atom
        for atom in atoms
        if _matches_rewritten_attachment_residue(atom, attachment)
        and atom.atom_name.upper() == attachment.target_atom_name.upper()
    ]
    if len(matches) != 1:
        raise ValueError(
            "Expected exactly one linked protein target atom "
            f"{attachment.target_chain}:{attachment.target_residue_number}:{attachment.target_atom_name}, "
            f"found {len(matches)}"
        )
    return matches[0]


def _matches_attachment_residue(atom: PdbAtomRecord, attachment: PdbAssemblyAttachment) -> bool:
    """Return whether an atom belongs to the attachment residue."""
    if atom.chain_id != attachment.target_chain:
        return False
    if atom.residue_number != attachment.target_residue_number:
        return False
    if (atom.insertion_code or "").upper() != attachment.target_insertion_code.upper():
        return False
    if attachment.target_residue_name is not None:
        return atom.residue_name.upper() == attachment.target_residue_name.upper()
    return True


def _matches_rewritten_attachment_residue(
    atom: PdbAtomRecord, attachment: PdbAssemblyAttachment
) -> bool:
    """Return whether an atom belongs to the rewritten attachment residue."""
    if atom.residue_number != attachment.target_residue_number:
        return False
    if (atom.insertion_code or "").upper() != attachment.target_insertion_code.upper():
        return False
    if atom.chain_id != attachment.target_chain and atom.residue_name != _pdb_safe_residue_name(
        _attachment_protein_product_resname(attachment)
    ):
        return False
    if attachment.target_residue_name is not None:
        allowed_residue_names = {
            attachment.target_residue_name.upper(),
            _pdb_safe_residue_name(_attachment_protein_product_resname(attachment)),
        }
        return atom.residue_name.upper() in allowed_residue_names
    return True


def _matches_generic_leaving_atom(
    atom: PdbAtomRecord,
    leaving_atom: PdbAtomRecord,
    attachment: PdbLinkageAttachment,
) -> bool:
    """Return whether an assembly atom is the resolved generic leaving atom."""
    if not _matches_attachment_residue(atom, attachment):
        return False
    if atom.atom_name.upper() != leaving_atom.atom_name.upper():
        return False
    if leaving_atom.atom_index is not None and atom.atom_index != leaving_atom.atom_index:
        return False
    if leaving_atom.serial is not None and atom.serial != leaving_atom.serial:
        return False
    return (
        atom.chain_id.upper() == leaving_atom.chain_id.upper()
        and atom.residue_name.upper() == leaving_atom.residue_name.upper()
        and atom.residue_number == leaving_atom.residue_number
        and (atom.insertion_code or "").upper() == (leaving_atom.insertion_code or "").upper()
    )


def _attachment_protein_product_resname(attachment: PdbAssemblyAttachment) -> str:
    """Return the protein-side product residue name for an attachment."""
    if isinstance(attachment, NhsLysPdbAttachment):
        return attachment.lysine_target_resname
    return attachment.protein_target_resname


def _attachment_modifier_product_resname(attachment: PdbAssemblyAttachment) -> str:
    """Return the modifier-side product residue name for an attachment."""
    if isinstance(attachment, NhsLysPdbAttachment):
        return attachment.polymer_target_resname
    return attachment.modifier_target_resname


def _non_crosslinked_modifier_residue_name(
    atom: PdbAtomRecord, attachment: PdbAssemblyAttachment
) -> str:
    """Return the output residue name for an unlinked modifier atom."""
    if isinstance(attachment, NhsLysPdbAttachment):
        return canonicalize_poc_residue_name(atom.residue_name)
    return _pdb_safe_residue_name(atom.residue_name)


def _is_hydrogen_atom(atom: PdbAtomRecord) -> bool:
    """Return whether an atom record represents hydrogen."""
    return (atom.element or "").upper() == "H" or atom.atom_name.upper().startswith("H")


def _is_removed_polymer_atom(atom: PdbAtomRecord, fragment: PlacedPolymerFragment) -> bool:
    """Return whether a polymer atom is selected as leaving group."""
    if fragment.leaving_atom_serials:
        return atom.serial is not None and atom.serial in fragment.leaving_atom_serials
    if fragment.leaving_atom_indices:
        return atom.atom_index is not None and atom.atom_index in fragment.leaving_atom_indices
    if fragment.leaving_atom_names:
        return atom.atom_name.upper() in {name.upper() for name in fragment.leaving_atom_names}
    return False


def _resolve_reactive_polymer_atom(fragment: PlacedPolymerFragment) -> PdbAtomRecord:
    """Resolve the polymer reactive atom from explicit selectors."""
    if fragment.reactive_atom_serial is not None:
        matches = [atom for atom in fragment.atoms if atom.serial == fragment.reactive_atom_serial]
    elif fragment.reactive_atom_index is not None:
        matches = [
            atom for atom in fragment.atoms if atom.atom_index == fragment.reactive_atom_index
        ]
    elif fragment.reactive_atom_name is not None:
        matches = [
            atom
            for atom in fragment.atoms
            if atom.atom_name.upper() == fragment.reactive_atom_name.upper()
        ]
    else:
        matches = []
    unique = _unique_atoms(matches)
    if len(unique) != 1:
        raise ValueError(
            f"Expected exactly one reactive atom in {fragment.name}, found {len(unique)}"
        )
    return unique[0]


def _resolve_polymer_ref(
    fragment: PlacedPolymerFragment, atom_ref: int | str
) -> PdbAtomRecord | None:
    """Resolve a polymer bond endpoint by serial, index, or atom name."""
    if isinstance(atom_ref, str):
        matches = [atom for atom in fragment.atoms if atom.atom_name.upper() == atom_ref.upper()]
    else:
        serial_matches = [atom for atom in fragment.atoms if atom.serial == atom_ref]
        index_matches = [atom for atom in fragment.atoms if atom.atom_index == atom_ref]
        matches = serial_matches or index_matches
    unique = _unique_atoms(matches)
    if len(unique) != 1:
        return None
    return unique[0]


def _unique_atoms(atoms: Iterable[PdbAtomRecord]) -> list[PdbAtomRecord]:
    """Return atom records with duplicate source identities removed."""
    unique: list[PdbAtomRecord] = []
    seen: set[tuple[int | None, int | None, str]] = set()
    for atom in atoms:
        key = (atom.atom_index, atom.serial, atom.atom_name)
        if key in seen:
            continue
        seen.add(key)
        unique.append(atom)
    return unique


def _atom_identity(atom: PdbAtomRecord) -> tuple[int | None, int | None, str, int, str]:
    """Return a stable source identity for an atom record."""
    return (
        atom.atom_index,
        atom.serial,
        atom.atom_name,
        atom.residue_number,
        atom.insertion_code,
    )


def _polymer_residue_key(atom: PdbAtomRecord) -> tuple[int, str]:
    """Return a stable key for polymer residue renumbering."""
    return (atom.residue_number, atom.insertion_code or "")


def _as_fragment_tuple(
    polymer_fragment: PlacedPolymerFragment | Sequence[PlacedPolymerFragment],
) -> tuple[PlacedPolymerFragment, ...]:
    """Normalize a single fragment or sequence to a tuple."""
    if isinstance(polymer_fragment, PlacedPolymerFragment):
        return (polymer_fragment,)
    return tuple(polymer_fragment)


def _as_attachment_tuple(
    attachment: PdbAssemblyAttachment | Sequence[PdbAssemblyAttachment],
) -> tuple[PdbAssemblyAttachment, ...]:
    """Normalize a single attachment or sequence to a tuple."""
    if isinstance(attachment, (NhsLysPdbAttachment, PdbLinkageAttachment)):
        return (attachment,)
    return tuple(attachment)


def _next_atom_serial(
    atom: PdbAtomRecord,
    next_serial: int,
    used_serials: set[int],
    options: CrosslinkedPdbAssemblyOptions,
) -> tuple[int, int]:
    """Choose the next output atom serial."""
    if options.preserve_serials and atom.serial is not None and atom.serial not in used_serials:
        used_serials.add(atom.serial)
        return atom.serial, next_serial
    while next_serial in used_serials:
        next_serial += 1
    serial = next_serial
    used_serials.add(serial)
    return serial, next_serial + 1


def _output_serial_for_atom(
    atom: PdbAtomRecord,
    *,
    serial_by_index: dict[int, int],
    serial_by_input_serial: dict[int, int],
    serial_by_name: dict[str, int] | None = None,
) -> int:
    """Resolve an atom's output serial from source-to-output maps."""
    if atom.atom_index is not None and atom.atom_index in serial_by_index:
        return serial_by_index[atom.atom_index]
    if atom.serial is not None and atom.serial in serial_by_input_serial:
        return serial_by_input_serial[atom.serial]
    if serial_by_name is not None and atom.atom_name in serial_by_name:
        return serial_by_name[atom.atom_name]
    raise ValueError(f"No output serial was assigned for atom {atom.atom_name}")


def _format_pdb_atom_line(atom: PdbAtomRecord) -> str:
    """Format one PDB ATOM/HETATM line with fixed-width fields."""
    serial = atom.serial if atom.serial is not None else 0
    atom_name = _format_atom_name(atom.atom_name)
    residue_name = _pdb_safe_residue_name(atom.residue_name)
    chain_id = (atom.chain_id or " ")[:1]
    insertion_code = (atom.insertion_code or " ")[:1]
    charge = (atom.charge or "")[:2]
    return (
        f"{atom.record_name:<6}{serial:5d} {atom_name}{atom.alt_loc[:1]:1}"
        f"{residue_name:>3} {chain_id}{atom.residue_number:4d}{insertion_code}   "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{atom.occupancy:6.2f}{atom.temp_factor:6.2f}          "
        f"{_format_element(atom.element, atom.atom_name):>2}{charge:>2}\n"
    )


def _format_ter_line(atom: PdbAtomRecord, serial: int) -> str:
    """Format a PDB TER line with residue identity preserved."""
    residue_name = _pdb_safe_residue_name(atom.residue_name)
    chain_id = (atom.chain_id or " ")[:1]
    insertion_code = (atom.insertion_code or " ")[:1]
    return f"TER   {serial:5d}      {residue_name:>3} {chain_id}{atom.residue_number:4d}{insertion_code}\n"


def _format_atom_name(atom_name: str) -> str:
    """Format a PDB atom name using the POC alignment convention."""
    stripped = atom_name.strip()
    if len(stripped) < 4:
        return f" {stripped:<3}"
    return f"{stripped:<4}"[:4]


def _format_element(element: str, atom_name: str) -> str:
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


def _format_conect_records(serial: int, bonded_serials: Sequence[int]) -> str:
    """Format symmetric CONECT rows for one source serial."""
    lines: list[str] = []
    for start in range(0, len(bonded_serials), _MAX_CONECT_TARGETS):
        chunk = bonded_serials[start : start + _MAX_CONECT_TARGETS]
        fields = "".join(f"{entry:5d}" for entry in chunk)
        lines.append(f"CONECT{serial:5d}{fields}\n")
    return "".join(lines)


def _format_link_line(
    protein_atom: PdbAtomRecord, polymer_serial: int, output_atoms: Sequence[PdbAtomRecord]
) -> str:
    """Format an optional LINK record for the crosslink."""
    polymer_atom = next(atom for atom in output_atoms if atom.serial == polymer_serial)
    return (
        f"LINK        {protein_atom.atom_name:<4} {protein_atom.residue_name:>3} "
        f"{protein_atom.chain_id:1}{protein_atom.residue_number:4d}{'':16}"
        f"{polymer_atom.atom_name:<4} {polymer_atom.residue_name:>3} "
        f"{polymer_atom.chain_id:1}{polymer_atom.residue_number:4d}\n"
    )


def _pdb_safe_residue_name(residue_name: str) -> str:
    """Return a three-character uppercase PDB residue name."""
    cleaned = "".join(ch for ch in (residue_name or "POL").upper() if ch.isalnum())
    return (cleaned or "POL")[:3]
