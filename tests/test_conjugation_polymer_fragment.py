"""Tests for generated polymer fragment conversion into PDB assembly."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.builders.conjugation.structure.pdb import NhsLysPdbAttachment, write_crosslinked_pdb
from polyzymd.builders.conjugation.polymer import (
    GeneratedPolymerFragment,
    PolymerFragmentAtom,
    PolymerFragmentResidue,
)


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    *,
    record_name: str = "ATOM",
    element: str = "C",
    insertion_code: str = "",
) -> str:
    """Format one fixed-width PDB atom line for assembly tests."""
    return (
        f"{record_name:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}{insertion_code[:1]:1}   {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )


def _protein_pdb(tmp_path: Path) -> Path:
    """Create a small protein PDB containing a reactive lysine."""
    path = tmp_path / "protein.pdb"
    path.write_text(
        _pdb_atom(1, "N", "LYS", "A", 23, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "CA", "LYS", "A", 23, 1.0, 0.0, 0.0)
        + _pdb_atom(3, "CE", "LYS", "A", 23, 1.5, 0.0, 0.0)
        + _pdb_atom(4, "NZ", "LYS", "A", 23, 2.0, 0.0, 0.0, element="N")
        + _pdb_atom(5, "HZ1", "LYS", "A", 23, 2.0, 0.7, 0.0, element="H")
        + _pdb_atom(6, "HZ2", "LYS", "A", 23, 2.0, -0.7, 0.0, element="H")
        + _pdb_atom(7, "HZ3", "LYS", "A", 23, 2.0, 0.0, 0.7, element="H")
        + "END\n"
    )
    return path


def _generated_fragment() -> GeneratedPolymerFragment:
    """Create a deterministic multi-residue SBMA/NHS/EGPMA fragment."""
    atoms = (
        PolymerFragmentAtom(
            atom_index=0,
            serial=101,
            atom_name="C1",
            residue_name="SBM",
            residue_number=10,
            sequence_index=0,
            x=5.0,
            y=0.0,
            z=0.0,
            element="C",
        ),
        PolymerFragmentAtom(
            atom_index=1,
            serial=102,
            atom_name="RC",
            residue_name="NHS",
            residue_number=11,
            sequence_index=1,
            x=3.3,
            y=0.0,
            z=0.0,
            element="C",
        ),
        PolymerFragmentAtom(
            atom_index=2,
            serial=103,
            atom_name="LG",
            residue_name="NHS",
            residue_number=11,
            sequence_index=1,
            x=4.2,
            y=1.0,
            z=0.0,
            element="O",
        ),
        PolymerFragmentAtom(
            atom_index=3,
            serial=104,
            atom_name="C2",
            residue_name="EGP",
            residue_number=12,
            sequence_index=2,
            x=5.8,
            y=0.0,
            z=0.0,
            element="C",
        ),
    )
    residues = (
        PolymerFragmentResidue(
            sequence_index=0,
            label="A",
            name="SBMA",
            residue_name="SBM",
            residue_number=10,
        ),
        PolymerFragmentResidue(
            sequence_index=1,
            label="C",
            name="NHS",
            residue_name="NHS",
            residue_number=11,
        ),
        PolymerFragmentResidue(
            sequence_index=2,
            label="B",
            name="EGPMA",
            residue_name="EGP",
            residue_number=12,
        ),
    )
    return GeneratedPolymerFragment(
        atoms=atoms,
        bonds=((0, 1), (1, 2), (1, 3)),
        residues=residues,
        sequence="ACB",
        reactive_atom_index=1,
        leaving_atom_indices=(2,),
        name="sbma_nhs_egpma_test",
    )


def test_generated_fragment_converts_to_placed_fragment_with_index_mapping():
    """Generated fragments should preserve serials, indices, and residue metadata."""
    generated = _generated_fragment()

    placed = generated.to_placed_fragment()

    assert placed.name == "sbma_nhs_egpma_test"
    assert placed.reactive_atom_index == 1
    assert placed.leaving_atom_indices == (2,)
    assert [atom.atom_index for atom in placed.atoms] == [0, 1, 2, 3]
    assert [atom.serial for atom in placed.atoms] == [101, 102, 103, 104]
    assert [residue.label for residue in generated.residues] == ["A", "C", "B"]


def test_generated_fragment_writes_crosslinked_nhx_lyx_pdb(tmp_path):
    """Converted fragments should feed the NHX/LYX crosslinked PDB writer."""
    output_path = tmp_path / "assembled.pdb"
    attachment = NhsLysPdbAttachment(
        target_chain="A",
        target_residue_number=23,
        nz_hydrogen_atom_names_to_remove=("HZ2", "HZ3"),
    )

    result = write_crosslinked_pdb(
        _protein_pdb(tmp_path),
        _generated_fragment().to_placed_fragment(),
        attachment,
        output_path,
    )

    lines = output_path.read_text().splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    polymer_lines = [line for line in atom_lines if line.startswith("HETATM") and line[21] == "C"]
    residue_names = [line[17:20].strip() for line in polymer_lines]
    atom_names = [line[12:16].strip() for line in polymer_lines]

    assert result.removed_polymer_atom_count == 1
    assert "LG" not in atom_names
    assert residue_names == ["SBM", "NHX", "EGP"]
    assert {line[17:20].strip() for line in atom_lines if line[21] == "A"} == {"LYX"}


def test_generated_fragment_rejects_ambiguous_reactive_selector():
    """Reactive selectors should resolve exactly one generated atom."""
    atoms = [atom.model_copy(update={"atom_name": "RC"}) for atom in _generated_fragment().atoms]

    with pytest.raises(ValueError, match="Reactive atom selector"):
        GeneratedPolymerFragment(atoms=tuple(atoms), reactive_atom_name="RC")


def test_generated_fragment_from_pdb_lines_infers_residues():
    """PDB-like atom sources should convert without Polymerist or RDKit."""
    lines = [
        _pdb_atom(201, "C1", "SBM", "C", 1, 0.0, 0.0, 0.0, record_name="HETATM"),
        _pdb_atom(202, "RC", "NHS", "C", 2, 1.0, 0.0, 0.0, record_name="HETATM"),
        _pdb_atom(203, "LG", "NHS", "C", 2, 1.5, 0.0, 0.0, record_name="HETATM"),
    ]

    fragment = GeneratedPolymerFragment.from_pdb_lines(
        lines,
        bonds=((201, 202), (202, 203)),
        reactive_atom_serial=202,
        leaving_atom_names=("LG",),
        sequence="AC",
    )

    assert [residue.residue_name for residue in fragment.residues] == ["SBM", "NHS"]
    assert fragment.to_placed_fragment().reactive_atom_serial == 202


def test_generated_fragment_preserves_insertion_codes_as_residue_keys():
    """Generated fragments should distinguish residues that share sequence numbers."""
    lines = [
        _pdb_atom(
            201,
            "C1",
            "SBM",
            "C",
            1,
            0.0,
            0.0,
            0.0,
            record_name="HETATM",
            insertion_code="A",
        ),
        _pdb_atom(
            202,
            "RC",
            "NHS",
            "C",
            1,
            1.0,
            0.0,
            0.0,
            record_name="HETATM",
            insertion_code="B",
        ),
        _pdb_atom(
            203,
            "LG",
            "NHS",
            "C",
            1,
            1.5,
            0.0,
            0.0,
            record_name="HETATM",
            insertion_code="B",
        ),
    ]

    fragment = GeneratedPolymerFragment.from_pdb_lines(
        lines,
        reactive_atom_serial=202,
        leaving_atom_serials=(203,),
        sequence="AC",
    )
    placed = fragment.to_placed_fragment()

    assert [residue.insertion_code for residue in fragment.residues] == ["A", "B"]
    assert [atom.insertion_code for atom in fragment.atoms] == ["A", "B", "B"]
    assert [atom.insertion_code for atom in placed.atoms] == ["A", "B", "B"]
