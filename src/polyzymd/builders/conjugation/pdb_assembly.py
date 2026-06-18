"""Pure-Python PDB assembly for NHS-Lys crosslinked polymer fragments."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field, model_validator

_ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")
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


class PdbAtomRecord(BaseModel):
    """Fixed-width PDB atom record data used by the crosslinked writer."""

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
        """Parse one fixed-width PDB ATOM/HETATM line.

        Parameters
        ----------
        line : str
            Input PDB line.
        atom_index : int or None, optional
            Optional zero-based source atom index, by default ``None``.

        Returns
        -------
        PdbAtomRecord
            Parsed atom record.

        Raises
        ------
        ValueError
            If the line is not an ATOM/HETATM record or lacks required fields.
        """
        if not line.startswith(_ATOM_RECORD_PREFIXES):
            raise ValueError("PDB atom parsing requires an ATOM or HETATM record")

        atom_name = _slice(line, 12, 16).strip()
        residue_name = _slice(line, 17, 20).strip()
        residue_number = _parse_int(_slice(line, 22, 26).strip())
        x = _parse_float(_slice(line, 30, 38).strip())
        y = _parse_float(_slice(line, 38, 46).strip())
        z = _parse_float(_slice(line, 46, 54).strip())
        if not atom_name or not residue_name or residue_number is None:
            raise ValueError(f"Invalid PDB atom record fields: {line.rstrip()}")
        if x is None or y is None or z is None:
            raise ValueError(f"Invalid PDB coordinate fields: {line.rstrip()}")

        record_name = "HETATM" if line.startswith("HETATM") else "ATOM"
        return cls(
            serial=_parse_int(_slice(line, 6, 11).strip()),
            atom_index=atom_index,
            atom_name=atom_name,
            residue_name=residue_name.upper(),
            chain_id=_slice(line, 21, 22).strip(),
            residue_number=residue_number,
            insertion_code=_slice(line, 26, 27).strip(),
            x=x,
            y=y,
            z=z,
            occupancy=_parse_float(_slice(line, 54, 60).strip()) or 1.0,
            temp_factor=_parse_float(_slice(line, 60, 66).strip()) or 0.0,
            element=_parse_element(line),
            charge=_slice(line, 78, 80).strip(),
            alt_loc=_slice(line, 16, 17).strip(),
            record_name=record_name,
        )


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
        atoms_list: list[PdbAtomRecord] = []
        for line in lines:
            if line.startswith(_ATOM_RECORD_PREFIXES):
                atoms_list.append(PdbAtomRecord.from_pdb_line(line, atom_index=len(atoms_list)))
        atoms = tuple(atoms_list)
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
    warnings: tuple[str, ...] = Field(default_factory=tuple)
    residue_mappings: dict[str, dict[str, int | str]] = Field(default_factory=dict)


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
    attachment: PdbAssemblyAttachment,
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
    warnings: list[str] = []

    protein_atoms = _parse_pdb_atoms(protein_path)
    kept_protein_atoms, removed_protein_atoms = _prepare_protein_atoms(
        protein_atoms,
        attachment,
        writer_options,
        warnings,
    )
    target_nz_atom = _resolve_attachment_atom(kept_protein_atoms, attachment)

    output_atoms: list[PdbAtomRecord] = []
    atom_lines: list[str] = []
    conect_map: dict[int, set[int]] = {}
    residue_mappings: dict[str, dict[str, int | str]] = {}
    removed_polymer_atoms: list[PdbAtomRecord] = []
    crosslink_pair: tuple[int, int] | None = None

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

    nz_serial = _output_serial_for_atom(
        target_nz_atom,
        serial_by_index=protein_serial_by_index,
        serial_by_input_serial=protein_serial_by_input_serial,
    )

    if writer_options.append_ter_records and output_atoms:
        atom_lines.append(_format_ter_line(output_atoms[-1], next_serial))
        next_serial += 1

    next_polymer_residue_number = 1
    for fragment_index, fragment in enumerate(fragments, start=1):
        fragment_result = _append_polymer_fragment(
            fragment,
            fragment_index=fragment_index,
            attachment=attachment,
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

        conect_map.setdefault(nz_serial, set()).add(fragment_result.reactive_serial)
        conect_map.setdefault(fragment_result.reactive_serial, set()).add(nz_serial)
        crosslink_pair = (nz_serial, fragment_result.reactive_serial)

        if writer_options.append_ter_records and fragment_result.atoms:
            atom_lines.append(_format_ter_line(fragment_result.atoms[-1], next_serial))
            next_serial += 1

    if crosslink_pair is None:
        raise ValueError("No polymer fragment was available for crosslinked PDB assembly")

    link_lines = []
    if writer_options.include_link_records:
        link_lines.append(_format_link_line(target_nz_atom, crosslink_pair[1], output_atoms))

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        handle.writelines(link_lines)
        handle.writelines(atom_lines)
        for serial in sorted(conect_map):
            bonded_serials = sorted(conect_map[serial])
            handle.write(_format_conect_records(serial, bonded_serials))
        handle.write("END\n")

    removed_atoms = [*removed_protein_atoms, *removed_polymer_atoms]
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
        warnings=tuple(warnings),
        residue_mappings=residue_mappings,
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
    removed_atoms = [atom for atom in fragment.atoms if _is_removed_polymer_atom(atom, fragment)]
    removed_keys = {_atom_identity(atom) for atom in removed_atoms}
    kept_atoms = [atom for atom in fragment.atoms if _atom_identity(atom) not in removed_keys]
    reactive_atom = _resolve_reactive_polymer_atom(fragment)
    if _atom_identity(reactive_atom) in removed_keys:
        raise ValueError("The polymer reactive atom cannot be listed as a leaving-group atom")

    residue_key_to_number: dict[tuple[int, str], int] = {}
    residue_mappings: dict[str, dict[str, int | str]] = {}
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
                "target_residue_number": residue_cursor,
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
    )


def _parse_pdb_atoms(path: Path) -> list[PdbAtomRecord]:
    """Parse all ATOM/HETATM records from a PDB file."""
    atoms: list[PdbAtomRecord] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(_ATOM_RECORD_PREFIXES):
                atoms.append(PdbAtomRecord.from_pdb_line(line, atom_index=len(atoms)))
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    return atoms


def _prepare_protein_atoms(
    atoms: Sequence[PdbAtomRecord],
    attachment: PdbAssemblyAttachment,
    options: CrosslinkedPdbAssemblyOptions,
    warnings: list[str],
) -> tuple[list[PdbAtomRecord], list[PdbAtomRecord]]:
    """Rename the linked residue and remove selected leaving atoms."""
    atom_indices_to_remove = _select_protein_atoms_to_remove(atoms, attachment, warnings)
    kept_atoms: list[PdbAtomRecord] = []
    removed_atoms: list[PdbAtomRecord] = []
    for atom in atoms:
        if atom.atom_index in atom_indices_to_remove:
            removed_atoms.append(atom)
            continue
        update = {"chain_id": options.protein_chain}
        if _matches_attachment_residue(atom, attachment):
            update["residue_name"] = _pdb_safe_residue_name(
                _attachment_protein_product_resname(attachment)
            )
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
    return (
        (atom.serial is not None and atom.serial in fragment.leaving_atom_serials)
        or (atom.atom_index is not None and atom.atom_index in fragment.leaving_atom_indices)
        or atom.atom_name.upper() in {name.upper() for name in fragment.leaving_atom_names}
    )


def _resolve_reactive_polymer_atom(fragment: PlacedPolymerFragment) -> PdbAtomRecord:
    """Resolve the polymer reactive atom from explicit selectors."""
    matches = [
        atom
        for atom in fragment.atoms
        if (
            fragment.reactive_atom_serial is not None
            and atom.serial == fragment.reactive_atom_serial
        )
        or (
            fragment.reactive_atom_index is not None
            and atom.atom_index == fragment.reactive_atom_index
        )
        or (
            fragment.reactive_atom_name is not None
            and atom.atom_name.upper() == fragment.reactive_atom_name.upper()
        )
    ]
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


def _parse_element(line: str) -> str:
    """Parse or infer a PDB element symbol."""
    element = _slice(line, 76, 78).strip()
    if element:
        return _format_element(element, _slice(line, 12, 16).strip())
    return _format_element("", _slice(line, 12, 16).strip())


def _slice(value: str, start: int, stop: int) -> str:
    """Return a safe fixed-width slice from a possibly short line."""
    if len(value) < stop:
        value = value.ljust(stop)
    return value[start:stop]


def _parse_int(value: str) -> int | None:
    """Parse an integer value, returning ``None`` for blanks or invalid text."""
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _parse_float(value: str) -> float | None:
    """Parse a float value, returning ``None`` for blanks or invalid text."""
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _pdb_safe_residue_name(residue_name: str) -> str:
    """Return a three-character uppercase PDB residue name."""
    cleaned = "".join(ch for ch in (residue_name or "POL").upper() if ch.isalnum())
    return (cleaned or "POL")[:3]
