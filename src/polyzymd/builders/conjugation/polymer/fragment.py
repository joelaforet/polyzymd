"""Generated polymer fragment records and PDB assembly adapter."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, model_validator

from polyzymd.builders.conjugation.structure.parsing import parse_pdb_atom_lines
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord, PlacedPolymerFragment


class PolymerFragmentResidue(BaseModel):
    """Residue-level metadata for a generated polymer fragment."""

    sequence_index: int = Field(..., ge=0)
    label: str | None = None
    name: str | None = None
    residue_name: str = Field(..., min_length=1, max_length=3)
    residue_number: int = Field(..., ge=1)
    insertion_code: str = Field("", max_length=1)


class PolymerFragmentAtom(BaseModel):
    """Atom record for a generated polymer before final PDB assembly."""

    atom_index: int = Field(..., ge=0)
    serial: int | None = Field(None, ge=1)
    atom_name: str = Field(..., min_length=1)
    residue_name: str = Field(..., min_length=1, max_length=3)
    residue_number: int = Field(..., ge=1)
    insertion_code: str = Field("", max_length=1)
    sequence_index: int | None = Field(None, ge=0)
    chain_id: str = Field("", max_length=1)
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    temp_factor: float = 0.0
    element: str = ""
    charge: str = Field("", max_length=2)
    formal_charge: int | None = None
    alt_loc: str = Field("", max_length=1)
    record_name: Literal["ATOM", "HETATM"] = "HETATM"

    @classmethod
    def from_pdb_atom(
        cls, atom: PdbAtomRecord, *, sequence_index: int | None = None
    ) -> PolymerFragmentAtom:
        """Build a fragment atom from an existing PDB atom record.

        Parameters
        ----------
        atom : PdbAtomRecord
            Source PDB atom record.
        sequence_index : int or None, optional
            Optional zero-based polymer sequence index, by default ``None``.

        Returns
        -------
        PolymerFragmentAtom
            Generated polymer atom record.
        """
        if atom.atom_index is None:
            raise ValueError("PDB atom records require atom_index for fragment conversion")
        return cls(
            atom_index=atom.atom_index,
            serial=atom.serial,
            atom_name=atom.atom_name,
            residue_name=atom.residue_name,
            residue_number=atom.residue_number,
            insertion_code=atom.insertion_code,
            sequence_index=sequence_index,
            chain_id=atom.chain_id,
            x=atom.x,
            y=atom.y,
            z=atom.z,
            occupancy=atom.occupancy,
            temp_factor=atom.temp_factor,
            element=atom.element,
            charge=atom.charge,
            formal_charge=_pdb_formal_charge(atom.charge),
            alt_loc=atom.alt_loc,
            record_name=atom.record_name,
        )

    def to_pdb_atom(self) -> PdbAtomRecord:
        """Convert this generated atom record to the PDB assembly atom model.

        Returns
        -------
        PdbAtomRecord
            PDB atom record preserving both serial and zero-based index selectors.
        """
        return PdbAtomRecord(
            serial=self.serial,
            atom_index=self.atom_index,
            atom_name=self.atom_name,
            residue_name=self.residue_name,
            chain_id=self.chain_id,
            residue_number=self.residue_number,
            insertion_code=self.insertion_code,
            x=self.x,
            y=self.y,
            z=self.z,
            occupancy=self.occupancy,
            temp_factor=self.temp_factor,
            element=self.element,
            charge=self.charge,
            alt_loc=self.alt_loc,
            record_name=self.record_name,
        )


def _pdb_formal_charge(value: str) -> int | None:
    """Parse a PDB formal-charge field into an integer when specified."""
    text = (value or "").strip()
    if not text:
        return None
    if len(text) == 2 and text[0].isdigit() and text[1] in "+-":
        magnitude = int(text[0])
        return magnitude if text[1] == "+" else -magnitude
    if len(text) == 2 and text[0] in "+-" and text[1].isdigit():
        magnitude = int(text[1])
        return magnitude if text[0] == "+" else -magnitude
    if text in {"+", "-"}:
        return 1 if text == "+" else -1
    return int(text)


class GeneratedPolymerFragment(BaseModel):
    """Generated multi-residue polymer fragment with assembly selectors."""

    atoms: tuple[PolymerFragmentAtom, ...] = Field(..., min_length=1)
    bonds: tuple[tuple[int | str, int | str], ...] = Field(default_factory=tuple)
    bond_orders: tuple[tuple[int | str, int | str, float], ...] = Field(default_factory=tuple)
    residues: tuple[PolymerFragmentResidue, ...] = Field(default_factory=tuple)
    sequence: str | None = None
    reactive_atom_serial: int | None = None
    reactive_atom_index: int | None = None
    reactive_atom_name: str | None = None
    leaving_atom_serials: tuple[int, ...] = Field(default_factory=tuple)
    leaving_atom_indices: tuple[int, ...] = Field(default_factory=tuple)
    leaving_atom_names: tuple[str, ...] = Field(default_factory=tuple)
    name: str = "polymer"

    @classmethod
    def from_pdb_lines(
        cls,
        lines: Iterable[str],
        *,
        bonds: Sequence[tuple[int | str, int | str]] = (),
        bond_orders: Sequence[tuple[int | str, int | str, float]] = (),
        reactive_atom_serial: int | None = None,
        reactive_atom_index: int | None = None,
        reactive_atom_name: str | None = None,
        leaving_atom_serials: Sequence[int] = (),
        leaving_atom_indices: Sequence[int] = (),
        leaving_atom_names: Sequence[str] = (),
        sequence: str | None = None,
        name: str = "polymer",
    ) -> GeneratedPolymerFragment:
        """Build a generated fragment from fixed-width PDB atom lines.

        Parameters
        ----------
        lines : iterable of str
            Input PDB lines containing ATOM/HETATM records.
        bonds : sequence of tuple, optional
            Internal polymer bonds by serial, index, or atom name, by default ``()``.
        bond_orders : sequence of tuple, optional
            Internal polymer bonds with bond-order metadata, by default ``()``.
        reactive_atom_serial : int or None, optional
            Reactive atom serial selector, by default ``None``.
        reactive_atom_index : int or None, optional
            Reactive atom zero-based index selector, by default ``None``.
        reactive_atom_name : str or None, optional
            Reactive atom-name selector, by default ``None``.
        leaving_atom_serials : sequence of int, optional
            Leaving-group atom serials, by default ``()``.
        leaving_atom_indices : sequence of int, optional
            Leaving-group atom indices, by default ``()``.
        leaving_atom_names : sequence of str, optional
            Leaving-group atom names, by default ``()``.
        sequence : str or None, optional
            Polymer monomer-label sequence, by default ``None``.
        name : str, optional
            Fragment name for diagnostics, by default ``"polymer"``.

        Returns
        -------
        GeneratedPolymerFragment
            Fragment ready for conversion to :class:`PlacedPolymerFragment`.
        """
        pdb_atoms = parse_pdb_atom_lines(lines)
        return cls.from_atom_records(
            pdb_atoms,
            bonds=bonds,
            bond_orders=bond_orders,
            reactive_atom_serial=reactive_atom_serial,
            reactive_atom_index=reactive_atom_index,
            reactive_atom_name=reactive_atom_name,
            leaving_atom_serials=leaving_atom_serials,
            leaving_atom_indices=leaving_atom_indices,
            leaving_atom_names=leaving_atom_names,
            sequence=sequence,
            name=name,
        )

    @classmethod
    def from_atom_records(
        cls,
        atoms: Sequence[PdbAtomRecord | PolymerFragmentAtom | Mapping[str, object]],
        *,
        bonds: Sequence[tuple[int | str, int | str]] = (),
        bond_orders: Sequence[tuple[int | str, int | str, float]] = (),
        residues: Sequence[PolymerFragmentResidue | Mapping[str, object]] = (),
        reactive_atom_serial: int | None = None,
        reactive_atom_index: int | None = None,
        reactive_atom_name: str | None = None,
        leaving_atom_serials: Sequence[int] = (),
        leaving_atom_indices: Sequence[int] = (),
        leaving_atom_names: Sequence[str] = (),
        sequence: str | None = None,
        name: str = "polymer",
    ) -> GeneratedPolymerFragment:
        """Build a generated fragment from typed or mapping atom records.

        Parameters
        ----------
        atoms : sequence
            Atom records from PDB parsing, Polymerist/RDKit adapters, or tests.
        bonds : sequence of tuple, optional
            Internal polymer bonds by serial, index, or atom name, by default ``()``.
        bond_orders : sequence of tuple, optional
            Internal polymer bonds with bond-order metadata, by default ``()``.
        residues : sequence, optional
            Optional explicit residue metadata, by default ``()``.
        reactive_atom_serial : int or None, optional
            Reactive atom serial selector, by default ``None``.
        reactive_atom_index : int or None, optional
            Reactive atom zero-based index selector, by default ``None``.
        reactive_atom_name : str or None, optional
            Reactive atom-name selector, by default ``None``.
        leaving_atom_serials : sequence of int, optional
            Leaving-group atom serials, by default ``()``.
        leaving_atom_indices : sequence of int, optional
            Leaving-group atom indices, by default ``()``.
        leaving_atom_names : sequence of str, optional
            Leaving-group atom names, by default ``()``.
        sequence : str or None, optional
            Polymer monomer-label sequence, by default ``None``.
        name : str, optional
            Fragment name for diagnostics, by default ``"polymer"``.

        Returns
        -------
        GeneratedPolymerFragment
            Validated generated polymer fragment.
        """
        fragment_atoms = tuple(_coerce_atom_record(atom) for atom in atoms)
        fragment_residues = (
            tuple(_coerce_residue_record(residue) for residue in residues)
            if residues
            else _infer_residues(fragment_atoms, sequence)
        )
        return cls(
            atoms=fragment_atoms,
            bonds=tuple(bonds),
            bond_orders=tuple(bond_orders),
            residues=fragment_residues,
            sequence=sequence,
            reactive_atom_serial=reactive_atom_serial,
            reactive_atom_index=reactive_atom_index,
            reactive_atom_name=reactive_atom_name,
            leaving_atom_serials=tuple(leaving_atom_serials),
            leaving_atom_indices=tuple(leaving_atom_indices),
            leaving_atom_names=tuple(leaving_atom_names),
            name=name,
        )

    @model_validator(mode="after")
    def validate_atom_references(self) -> GeneratedPolymerFragment:
        """Validate atom identity maps and reactive/leaving selectors."""
        atom_indices = [atom.atom_index for atom in self.atoms]
        if len(set(atom_indices)) != len(atom_indices):
            raise ValueError("Generated polymer atom indices must be unique")

        serials = [atom.serial for atom in self.atoms if atom.serial is not None]
        if len(set(serials)) != len(serials):
            raise ValueError("Generated polymer atom serials must be unique when provided")

        if (
            self.reactive_atom_serial is None
            and self.reactive_atom_index is None
            and self.reactive_atom_name is None
        ):
            raise ValueError("Generated polymer fragments require a reactive atom selector")

        if len(_matching_atoms(self, reactive=True)) != 1:
            raise ValueError("Reactive atom selector must resolve exactly one generated atom")

        leaving_atoms = _matching_atoms(self, reactive=False)
        reactive_atom = _matching_atoms(self, reactive=True)[0]
        if reactive_atom in leaving_atoms:
            raise ValueError("The reactive atom cannot also be a leaving atom")

        return self

    def to_placed_fragment(self, *, name: str | None = None) -> PlacedPolymerFragment:
        """Convert this generated fragment to the PDB assembly adapter model.

        Parameters
        ----------
        name : str or None, optional
            Override fragment name, by default ``None``.

        Returns
        -------
        PlacedPolymerFragment
            Existing PDB assembly fragment model.
        """
        return PlacedPolymerFragment(
            atoms=tuple(atom.to_pdb_atom() for atom in self.atoms),
            bonds=self.bonds,
            bond_orders=self.bond_orders,
            leaving_atom_serials=self.leaving_atom_serials,
            leaving_atom_indices=self.leaving_atom_indices,
            leaving_atom_names=self.leaving_atom_names,
            reactive_atom_serial=self.reactive_atom_serial,
            reactive_atom_index=self.reactive_atom_index,
            reactive_atom_name=self.reactive_atom_name,
            name=name or self.name,
        )


class PreparedFragment(GeneratedPolymerFragment):
    """Authoritative provider-neutral fragment used by conjugation construction."""

    model_config = {"arbitrary_types_allowed": True}

    source_identity: str = Field(..., min_length=1)
    source_kind: Literal["polymer", "smiles", "pdb_fragment"]
    sidecars: dict[str, Path] = Field(default_factory=dict)
    provenance: dict[str, Any] = Field(default_factory=dict)
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)

    @classmethod
    def from_generated_fragment(
        cls,
        fragment: GeneratedPolymerFragment,
        *,
        source_identity: str,
        source_kind: Literal["polymer", "smiles", "pdb_fragment"],
        sidecars: Mapping[str, Path] | None = None,
        provenance: Mapping[str, Any] | None = None,
        diagnostics: Sequence[str] = (),
    ) -> PreparedFragment:
        """Promote generation output into the construction contract once."""
        return cls(
            **fragment.model_dump(),
            source_identity=source_identity,
            source_kind=source_kind,
            sidecars=dict(sidecars or {}),
            provenance=dict(provenance or {}),
            diagnostics=tuple(diagnostics),
        )


def _coerce_atom_record(
    atom: PdbAtomRecord | PolymerFragmentAtom | Mapping[str, object],
) -> PolymerFragmentAtom:
    """Coerce supported atom record inputs to a generated polymer atom."""
    if isinstance(atom, PolymerFragmentAtom):
        return atom
    if isinstance(atom, PdbAtomRecord):
        return PolymerFragmentAtom.from_pdb_atom(atom)
    return PolymerFragmentAtom.model_validate(atom)


def _coerce_residue_record(
    residue: PolymerFragmentResidue | Mapping[str, object],
) -> PolymerFragmentResidue:
    """Coerce supported residue record inputs to generated residue metadata."""
    if isinstance(residue, PolymerFragmentResidue):
        return residue
    return PolymerFragmentResidue.model_validate(residue)


def _infer_residues(
    atoms: Sequence[PolymerFragmentAtom], sequence: str | None
) -> tuple[PolymerFragmentResidue, ...]:
    """Infer residue metadata from atom residue identifiers."""
    residues: list[PolymerFragmentResidue] = []
    seen: set[tuple[int, str, str]] = set()
    for atom in atoms:
        key = (atom.residue_number, atom.insertion_code, atom.residue_name)
        if key in seen:
            continue
        seen.add(key)
        sequence_index = len(residues)
        residues.append(
            PolymerFragmentResidue(
                sequence_index=sequence_index,
                label=(
                    sequence[sequence_index]
                    if sequence and sequence_index < len(sequence)
                    else None
                ),
                residue_name=atom.residue_name,
                residue_number=atom.residue_number,
                insertion_code=atom.insertion_code,
            )
        )
    return tuple(residues)


def _matching_atoms(
    fragment: GeneratedPolymerFragment, *, reactive: bool
) -> list[PolymerFragmentAtom]:
    """Return atoms matched by reactive or leaving-group selectors."""
    if reactive:
        serials = {fragment.reactive_atom_serial}
        indices = {fragment.reactive_atom_index}
        names = {fragment.reactive_atom_name.upper()} if fragment.reactive_atom_name else set()
    else:
        serials = set(fragment.leaving_atom_serials)
        indices = set(fragment.leaving_atom_indices)
        names = {name.upper() for name in fragment.leaving_atom_names}

    return [
        atom
        for atom in fragment.atoms
        if (atom.serial is not None and atom.serial in serials)
        or atom.atom_index in indices
        or atom.atom_name.upper() in names
    ]
