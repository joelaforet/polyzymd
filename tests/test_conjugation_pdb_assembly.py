"""Tests for pure-Python crosslinked PDB assembly."""

from __future__ import annotations

from pathlib import Path

from polyzymd.builders.conjugation.structure.inspection import inspect_pdb_structure
from polyzymd.builders.conjugation.structure.pdb import (
    NhsLysPdbAttachment,
    PdbAtomRecord,
    PdbLinkageAttachment,
    PlacedPolymerFragment,
    write_crosslinked_pdb,
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


def _three_mer_polymer_fragment_with_internal_reactive_residue() -> PlacedPolymerFragment:
    """Create raw Polymerist-like order: reactive middle residue first, terminals after."""
    atoms = (
        PdbAtomRecord(
            serial=201,
            atom_index=0,
            atom_name="N1",
            residue_name="NH1",
            chain_id="Z",
            residue_number=20,
            x=0.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=202,
            atom_index=1,
            atom_name="RC",
            residue_name="NH1",
            chain_id="Z",
            residue_number=20,
            x=1.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=203,
            atom_index=2,
            atom_name="A1",
            residue_name="SB1",
            chain_id="Z",
            residue_number=10,
            x=-1.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=204,
            atom_index=3,
            atom_name="B1",
            residue_name="EG1",
            chain_id="Z",
            residue_number=30,
            x=2.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
    )
    return PlacedPolymerFragment(
        atoms=atoms,
        bonds=((201, 202), (201, 203), (202, 204)),
        reactive_atom_serial=202,
        name="public_three_mer_like",
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


def test_crosslinked_writer_orders_linear_polymer_by_connectivity(tmp_path):
    """A graph-internal reactive residue should not be emitted as chain-C residue 1."""
    protein_path, _original = _protein_pdb(tmp_path)
    output_path = tmp_path / "assembled_three_mer.pdb"
    attachment = NhsLysPdbAttachment(
        target_chain="A",
        target_residue_number=23,
        nz_hydrogen_atom_names_to_remove=("HZ2", "HZ3"),
    )

    write_crosslinked_pdb(
        protein_path,
        _three_mer_polymer_fragment_with_internal_reactive_residue(),
        attachment,
        output_path,
    )

    polymer_lines = [
        line
        for line in output_path.read_text().splitlines()
        if line.startswith("HETATM") and line[21] == "C"
    ]
    residue_order = list(
        dict.fromkeys((line[17:20].strip(), line[22:26].strip()) for line in polymer_lines)
    )

    assert residue_order == [("SBM", "1"), ("NHX", "2"), ("EGP", "3")]


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
        modifier_link_atom=_generic_polymer_fragment().atoms[1],
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
        modifier_link_atom=_generic_polymer_fragment().atoms[1],
        modifier_leaving_atoms_to_remove=(_generic_polymer_fragment().atoms[2],),
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


def test_duplicate_modifier_names_use_serial_selectors_for_reactive_and_leaving(tmp_path):
    """Serial selectors should disambiguate repeated modifier reactive and leaving names."""
    protein_path = tmp_path / "generic_protein.pdb"
    protein_path.write_text(
        _pdb_atom(1, "SG", "CYS", "A", 10, 0.0, 0.0, 0.0, element="S") + "END\n"
    )
    fragment = _duplicate_name_polymer_fragment()
    output_path = tmp_path / "generic_assembled.pdb"
    attachment = PdbLinkageAttachment(
        target_chain="A",
        target_residue_name="CYS",
        target_residue_number=10,
        target_atom_name="SG",
        modifier_link_atom=fragment.atoms[2],
        modifier_leaving_atoms_to_remove=(fragment.atoms[1],),
        protein_target_resname="CYX",
        modifier_target_resname="MXL",
    )

    write_crosslinked_pdb(protein_path, fragment, attachment, output_path)

    polymer_lines = [
        line
        for line in output_path.read_text().splitlines()
        if line.startswith("HETATM") and line[21] == "C"
    ]
    rc_residues = [line[17:20].strip() for line in polymer_lines if line[12:16].strip() == "RC"]
    lg_x_coords = [line[30:38].strip() for line in polymer_lines if line[12:16].strip() == "LG"]

    assert rc_residues == ["SB1", "MXL"]
    assert lg_x_coords == ["6.000"]


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


def test_crosslinked_writer_preserves_source_components_and_chain_b_selector(tmp_path):
    """Canonical chain A output must retain source boundaries and chain-B selection."""
    from openmm.app import PDBFile

    protein_path = tmp_path / "multichain_protein.pdb"
    protein_path.write_text(
        _pdb_atom(1, "N", "CYS", "A", 1, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "CA", "CYS", "A", 1, 1.4, 0.0, 0.0)
        + _pdb_atom(3, "C", "CYS", "A", 1, 2.8, 0.0, 0.0)
        + _pdb_atom(4, "O", "CYS", "A", 1, 3.4, 1.0, 0.0, element="O")
        + _pdb_atom(5, "CB", "CYS", "A", 1, 1.4, -1.4, 0.0)
        + _pdb_atom(6, "SG", "CYS", "A", 1, 1.4, -2.8, 0.0, element="S")
        + "TER\n"
        + _pdb_atom(7, "N", "CYS", "B", 2, 4.1, 0.0, 0.0, element="N")
        + _pdb_atom(8, "CA", "CYS", "B", 2, 5.5, 0.0, 0.0)
        + _pdb_atom(9, "C", "CYS", "B", 2, 6.9, 0.0, 0.0)
        + _pdb_atom(10, "O", "CYS", "B", 2, 7.5, 1.0, 0.0, element="O")
        + _pdb_atom(11, "CB", "CYS", "B", 2, 5.5, -1.4, 0.0)
        + _pdb_atom(12, "SG", "CYS", "B", 2, 2.9, -2.8, 0.0, element="S")
        + _pdb_atom(13, "N", "ASN", "B", 3, 8.2, 0.0, 0.0, element="N")
        + _pdb_atom(14, "CA", "ASN", "B", 3, 9.6, 0.0, 0.0)
        + _pdb_atom(15, "C", "ASN", "B", 3, 11.0, 0.0, 0.0)
        + _pdb_atom(16, "O", "ASN", "B", 3, 11.6, 1.0, 0.0, element="O")
        + _pdb_atom(17, "CB", "ASN", "B", 3, 9.6, -1.4, 0.0)
        + _pdb_atom(18, "CG", "ASN", "B", 3, 10.8, -2.1, 0.0)
        + _pdb_atom(19, "OD1", "ASN", "B", 3, 11.9, -1.6, 0.0, element="O")
        + _pdb_atom(20, "ND2", "ASN", "B", 3, 10.7, -3.4, 0.0, element="N")
        + _pdb_atom(21, "HD21", "ASN", "B", 3, 11.5, -4.0, 0.0, element="H")
        + _pdb_atom(22, "HD22", "ASN", "B", 3, 9.8, -3.8, 0.0, element="H")
        + "CONECT    6   12\n"
        + "END\n"
    )
    fragment = _generic_polymer_fragment()
    output_path = tmp_path / "multichain_assembled.pdb"
    attachment = PdbLinkageAttachment(
        target_chain="B",
        target_residue_name="ASN",
        target_residue_number=3,
        target_atom_name="ND2",
        modifier_link_atom=fragment.atoms[1],
        protein_leaving_atoms_to_remove=(
            PdbAtomRecord(
                serial=21,
                atom_index=20,
                atom_name="HD21",
                residue_name="ASN",
                chain_id="B",
                residue_number=3,
                x=11.5,
                y=-4.0,
                z=0.0,
                element="H",
            ),
        ),
        protein_target_resname="NLN",
        modifier_target_resname="MXL",
    )

    write_crosslinked_pdb(protein_path, fragment, attachment, output_path)

    lines = output_path.read_text().splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    assert {line[21] for line in atom_lines if line.startswith("ATOM")} == {"A"}
    assert any(line.startswith("TER") and line[22:26].strip() == "1" for line in lines)
    nln_atoms = {
        line[12:16].strip()
        for line in atom_lines
        if line[17:20].strip() == "NLN" and line[22:26].strip() == "3"
    }
    assert "HD21" in nln_atoms
    assert "HD22" not in nln_atoms

    topology = PDBFile(str(output_path)).topology
    residue_by_id = {
        (residue.name, residue.id): residue
        for residue in topology.residues()
        if residue.chain.id == "A"
    }
    chain_a_c = next(atom for atom in residue_by_id[("CYS", "1")].atoms() if atom.name == "C")
    chain_b_n = next(atom for atom in residue_by_id[("CYS", "2")].atoms() if atom.name == "N")
    chain_a_sg = next(atom for atom in residue_by_id[("CYS", "1")].atoms() if atom.name == "SG")
    chain_b_sg = next(atom for atom in residue_by_id[("CYS", "2")].atoms() if atom.name == "SG")
    inferred_bonds = {frozenset((left, right)) for left, right in topology.bonds()}
    assert frozenset((chain_a_c, chain_b_n)) not in inferred_bonds
    assert frozenset((chain_a_sg, chain_b_sg)) in inferred_bonds


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


def _duplicate_name_polymer_fragment() -> PlacedPolymerFragment:
    """Create a modifier with duplicate reactive and leaving atom names."""
    atoms = (
        PdbAtomRecord(
            serial=301,
            atom_index=0,
            atom_name="RC",
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
            serial=302,
            atom_index=1,
            atom_name="LG",
            residue_name="LNK",
            chain_id="Z",
            residue_number=21,
            x=5.5,
            y=0.0,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=303,
            atom_index=2,
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
            serial=304,
            atom_index=3,
            atom_name="LG",
            residue_name="SB2",
            chain_id="Z",
            residue_number=22,
            x=6.0,
            y=0.0,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
    )
    return PlacedPolymerFragment(
        atoms=atoms,
        bonds=((301, 303), (302, 303), (303, 304)),
        leaving_atom_serials=(302,),
        reactive_atom_serial=303,
        reactive_atom_name="RC",
        leaving_atom_names=("LG",),
        name="duplicate_name_modifier_test",
    )
