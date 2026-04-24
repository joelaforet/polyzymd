"""Tests for pure-Python crosslinked PDB assembly."""

from __future__ import annotations

from pathlib import Path

from polyzymd.builders.conjugation.pdb_assembly import (
    NhsLysPdbAttachment,
    PdbAtomRecord,
    PdbLinkageAttachment,
    PlacedPolymerFragment,
    write_crosslinked_pdb,
)
from polyzymd.builders.conjugation.structure_inspection import inspect_pdb_structure


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


def _protein_pdb(tmp_path: Path) -> tuple[Path, str]:
    """Create a small protein PDB containing a lysine NZ with three hydrogens."""
    original = (
        _pdb_atom(1, "N", "LYS", "A", 23, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "CA", "LYS", "A", 23, 1.0, 0.0, 0.0)
        + _pdb_atom(3, "CE", "LYS", "A", 23, 1.5, 0.0, 0.0)
        + _pdb_atom(4, "NZ", "LYS", "A", 23, 2.0, 0.0, 0.0, element="N")
        + _pdb_atom(5, "HZ1", "LYS", "A", 23, 2.0, 0.7, 0.0, element="H")
        + _pdb_atom(6, "HZ2", "LYS", "A", 23, 2.0, -0.7, 0.0, element="H")
        + _pdb_atom(7, "HZ3", "LYS", "A", 23, 2.0, 0.0, 0.7, element="H")
        + _pdb_atom(8, "N", "ALA", "A", 24, 4.0, 0.0, 0.0, element="N")
        + "END\n"
    )
    path = tmp_path / "protein.pdb"
    path.write_text(original)
    return path, original


def _polymer_fragment() -> PlacedPolymerFragment:
    """Create a placed polymer fragment with one leaving-group atom."""
    atoms = (
        PdbAtomRecord(
            serial=101,
            atom_index=0,
            atom_name="C1",
            residue_name="SB1",
            chain_id="Z",
            residue_number=10,
            x=5.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=102,
            atom_index=1,
            atom_name="RC",
            residue_name="NH1",
            chain_id="Z",
            residue_number=11,
            x=3.3,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=103,
            atom_index=2,
            atom_name="O1",
            residue_name="NH1",
            chain_id="Z",
            residue_number=11,
            x=3.8,
            y=0.5,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=104,
            atom_index=3,
            atom_name="LG",
            residue_name="NH1",
            chain_id="Z",
            residue_number=11,
            x=4.2,
            y=1.0,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
    )
    return PlacedPolymerFragment(
        atoms=atoms,
        bonds=((101, 102), (102, 103), (102, 104)),
        leaving_atom_serials=(104,),
        reactive_atom_serial=102,
        name="egpma_nhs_test",
    )


def _write_assembled(tmp_path: Path) -> tuple[Path, str]:
    """Write a crosslinked assembly used by multiple tests."""
    protein_path, original = _protein_pdb(tmp_path)
    output_path = tmp_path / "assembled.pdb"
    attachment = NhsLysPdbAttachment(
        target_chain="A",
        target_residue_number=23,
        nz_hydrogen_atom_names_to_remove=("HZ2", "HZ3"),
    )

    write_crosslinked_pdb(protein_path, _polymer_fragment(), attachment, output_path)

    return output_path, original


def test_crosslinked_writer_renames_lys_and_removes_selected_nz_hydrogens(tmp_path):
    """The linked lysine should become LYX while one NZ hydrogen is retained."""
    output_path, _original = _write_assembled(tmp_path)
    atom_lines = [
        line for line in output_path.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))
    ]

    lys_lines = [line for line in atom_lines if line[21] == "A" and line[22:26].strip() == "23"]
    atom_names = {line[12:16].strip() for line in lys_lines}

    assert {line[17:20].strip() for line in lys_lines} == {"LYX"}
    assert "HZ1" in atom_names
    assert "HZ2" not in atom_names
    assert "HZ3" not in atom_names


def test_crosslinked_writer_appends_polymer_chain_c_and_canonicalizes_names(tmp_path):
    """Polymer atoms should be appended on chain C with POC-safe residue names."""
    output_path, _original = _write_assembled(tmp_path)
    polymer_lines = [
        line
        for line in output_path.read_text().splitlines()
        if line.startswith("HETATM") and line[21] == "C"
    ]

    assert [line[12:16].strip() for line in polymer_lines] == ["C1", "RC", "O1"]
    assert [line[17:20].strip() for line in polymer_lines] == ["SBM", "NHX", "NHX"]
    assert {line[22:26].strip() for line in polymer_lines} == {"1", "2"}


def test_generic_writer_removes_resolved_non_hydrogen_leaving_atom(tmp_path):
    """Generic assembly should remove exact non-hydrogen leaving atoms only."""
    protein_path = tmp_path / "generic_protein.pdb"
    protein_path.write_text(
        _pdb_atom(1, "N", "CYS", "A", 10, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "CA", "CYS", "A", 10, 1.0, 0.0, 0.0)
        + _pdb_atom(3, "SG", "CYS", "A", 10, 2.0, 0.0, 0.0, element="S")
        + _pdb_atom(4, "OXT", "CYS", "A", 10, 2.5, 0.0, 0.0, element="O")
        + _pdb_atom(5, "N", "GLY", "A", 11, 4.0, 0.0, 0.0, element="N")
        + _pdb_atom(6, "OXT", "GLY", "A", 11, 5.0, 0.0, 0.0, element="O")
        + "END\n"
    )
    output_path = tmp_path / "generic_assembled.pdb"
    attachment = PdbLinkageAttachment(
        target_chain="A",
        target_residue_name="CYS",
        target_residue_number=10,
        target_atom_name="SG",
        protein_leaving_atoms_to_remove=(
            PdbAtomRecord(
                serial=4,
                atom_index=3,
                atom_name="OXT",
                residue_name="CYS",
                chain_id="A",
                residue_number=10,
                x=2.5,
                y=0.0,
                z=0.0,
                element="O",
            ),
        ),
        protein_target_resname="CYX",
        modifier_target_resname="MXL",
    )

    write_crosslinked_pdb(protein_path, _generic_polymer_fragment(), attachment, output_path)

    atom_lines = [
        line for line in output_path.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))
    ]
    selected_lines = [
        line for line in atom_lines if line[21] == "A" and line[22:26].strip() == "10"
    ]
    duplicate_residue_lines = [
        line for line in atom_lines if line[21] == "A" and line[22:26].strip() == "11"
    ]
    assert {line[17:20].strip() for line in selected_lines} == {"CYX"}
    assert "OXT" not in {line[12:16].strip() for line in selected_lines}
    assert "OXT" in {line[12:16].strip() for line in duplicate_residue_lines}


def test_generic_writer_preserves_unlinked_modifier_residue_names(tmp_path):
    """Generic assembly should avoid POC canonicalization for unlinked residues."""
    protein_path = tmp_path / "generic_protein.pdb"
    protein_path.write_text(
        _pdb_atom(1, "SG", "CYS", "A", 10, 0.0, 0.0, 0.0, element="S") + "END\n"
    )
    output_path = tmp_path / "generic_assembled.pdb"
    attachment = PdbLinkageAttachment(
        target_chain="A",
        target_residue_name="CYS",
        target_residue_number=10,
        target_atom_name="SG",
        protein_target_resname="CYX",
        modifier_target_resname="MXL",
    )

    write_crosslinked_pdb(protein_path, _generic_polymer_fragment(), attachment, output_path)

    polymer_lines = [
        line
        for line in output_path.read_text().splitlines()
        if line.startswith("HETATM") and line[21] == "C"
    ]

    assert [line[12:16].strip() for line in polymer_lines] == ["S1", "RC", "S2"]
    assert [line[17:20].strip() for line in polymer_lines] == ["SB1", "MXL", "SB2"]


def test_crosslinked_writer_emits_reciprocal_conect_without_removed_atoms(tmp_path):
    """Internal polymer bonds and the crosslink should be reciprocal and clean."""
    output_path, _original = _write_assembled(tmp_path)
    lines = output_path.read_text().splitlines()
    serial_by_name = {
        line[12:16].strip(): int(line[6:11])
        for line in lines
        if line.startswith(("ATOM", "HETATM"))
    }
    conect_pairs = _conect_pairs(lines)

    nz = serial_by_name["NZ"]
    reactive = serial_by_name["RC"]
    c1 = serial_by_name["C1"]
    o1 = serial_by_name["O1"]

    assert (nz, reactive) in conect_pairs
    assert (reactive, nz) in conect_pairs
    assert (c1, reactive) in conect_pairs
    assert (reactive, c1) in conect_pairs
    assert (reactive, o1) in conect_pairs
    assert (o1, reactive) in conect_pairs
    output_serials = set(serial_by_name.values())
    assert all(
        source in output_serials and target in output_serials for source, target in conect_pairs
    )


def test_crosslinked_writer_does_not_modify_input_protein(tmp_path):
    """The source protein file should be read-only for assembly."""
    protein_path, original = _protein_pdb(tmp_path)
    output_path = tmp_path / "assembled.pdb"
    attachment = NhsLysPdbAttachment(
        target_chain="A",
        target_residue_number=23,
        nz_hydrogen_atom_names_to_remove=("HZ2", "HZ3"),
    )

    write_crosslinked_pdb(protein_path, _polymer_fragment(), attachment, output_path)

    assert protein_path.read_text() == original


def test_crosslinked_writer_matches_residue_name_number_and_insertion_code(tmp_path):
    """Assembly should rewrite only the exact selected protein residue identity."""
    protein_path = tmp_path / "protein_insertions.pdb"
    protein_path.write_text(
        _pdb_atom(1, "NZ", "LYS", "A", 23, 0.0, 0.0, 0.0, element="N", insertion_code="A")
        + _pdb_atom(2, "HZ2", "LYS", "A", 23, 0.0, 0.6, 0.0, element="H", insertion_code="A")
        + _pdb_atom(3, "HZ3", "LYS", "A", 23, 0.0, 0.0, 0.6, element="H", insertion_code="A")
        + _pdb_atom(4, "NZ", "LYS", "A", 23, 2.0, 0.0, 0.0, element="N", insertion_code="B")
        + _pdb_atom(5, "HZ2", "LYS", "A", 23, 2.0, 0.6, 0.0, element="H", insertion_code="B")
        + _pdb_atom(6, "HZ3", "LYS", "A", 23, 2.0, 0.0, 0.6, element="H", insertion_code="B")
        + _pdb_atom(7, "NZ", "ALA", "A", 23, 4.0, 0.0, 0.0, element="N", insertion_code="B")
        + "END\n"
    )
    output_path = tmp_path / "assembled_insertions.pdb"
    attachment = NhsLysPdbAttachment(
        target_chain="A",
        target_residue_name="LYS",
        target_residue_number=23,
        target_insertion_code="B",
        nz_hydrogen_atom_names_to_remove=("HZ2", "HZ3"),
    )

    write_crosslinked_pdb(protein_path, _polymer_fragment(), attachment, output_path)

    atom_lines = [
        line for line in output_path.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))
    ]
    residue_a_lines = [line for line in atom_lines if line[21] == "A" and line[26] == "A"]
    residue_b_lines = [line for line in atom_lines if line[21] == "A" and line[26] == "B"]
    assert {line[17:20].strip() for line in residue_a_lines} == {"LYS"}
    assert any(line[17:20].strip() == "LYX" for line in residue_b_lines)
    assert any(line[17:20].strip() == "ALA" for line in residue_b_lines)
    assert "HZ2" in {line[12:16].strip() for line in residue_a_lines}
    assert "HZ2" not in {
        line[12:16].strip() for line in residue_b_lines if line[17:20].strip() == "LYX"
    }


def test_inspection_detects_crosslink_and_lyx_is_not_dirty_free_ligand(tmp_path):
    """Inspection should treat LYX as protein-like and detect the chain C crosslink."""
    output_path, _original = _write_assembled(tmp_path)

    inspection = inspect_pdb_structure(output_path)

    assert inspection.residue_name_counts["LYX"] == 5
    assert all(
        residue.residue_name != "LYX" for residue in inspection.noncanonical_residue_candidates
    )
    ligand_names = {
        residue.residue_name
        for residue in inspection.noncanonical_residue_candidates
        if residue.category == "ligand_cocrystal"
    }
    assert "LYX" not in ligand_names
    candidates = inspection.covalent_attachment_candidates
    assert len(candidates) == 1
    assert candidates[0].protein_residue_name == "LYX"
    assert candidates[0].candidate_residue_name == "NHX"
    assert candidates[0].candidate_chain_id == "C"


def _conect_pairs(lines: list[str]) -> set[tuple[int, int]]:
    """Parse directed CONECT pairs from PDB lines."""
    pairs: set[tuple[int, int]] = set()
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        source = int(line[6:11])
        for start in range(11, len(line), 5):
            target = line[start : start + 5].strip()
            if target:
                pairs.add((source, int(target)))
    return pairs


def _generic_polymer_fragment() -> PlacedPolymerFragment:
    """Create a generic placed modifier fragment with distinct residue names."""
    atoms = (
        PdbAtomRecord(
            serial=201,
            atom_index=0,
            atom_name="S1",
            residue_name="SB1",
            chain_id="Z",
            residue_number=20,
            x=5.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=202,
            atom_index=1,
            atom_name="RC",
            residue_name="LNK",
            chain_id="Z",
            residue_number=21,
            x=3.3,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=203,
            atom_index=2,
            atom_name="LG",
            residue_name="LNK",
            chain_id="Z",
            residue_number=21,
            x=3.8,
            y=0.5,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=204,
            atom_index=3,
            atom_name="S2",
            residue_name="SB2",
            chain_id="Z",
            residue_number=22,
            x=6.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
    )
    return PlacedPolymerFragment(
        atoms=atoms,
        bonds=((201, 202), (202, 203), (202, 204)),
        leaving_atom_serials=(203,),
        reactive_atom_serial=202,
        name="generic_modifier_test",
    )
