"""Tests for clean-PDB chain normalization planning and writing."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.builders.conjugation.structure.inspection import inspect_pdb_structure
from polyzymd.builders.conjugation.structure.normalization import (
    default_cleaned_pdb_path,
    plan_pdb_chain_normalization,
    write_normalized_pdb,
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
) -> str:
    """Format one fixed-width PDB atom line for normalization tests."""
    return (
        f"{record_name:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )


def _pdb_link(
    atom_name_1: str,
    residue_name_1: str,
    chain_id_1: str,
    residue_number_1: int,
    atom_name_2: str,
    residue_name_2: str,
    chain_id_2: str,
    residue_number_2: int,
) -> str:
    """Format one fixed-width LINK line for normalization tests."""
    return (
        f"LINK        {atom_name_1:<4} {residue_name_1:>3} {chain_id_1:1}"
        f"{residue_number_1:4d}{'':16}{atom_name_2:<4} {residue_name_2:>3} "
        f"{chain_id_2:1}{residue_number_2:4d}\n"
    )


def _poc_path(name: str) -> Path:
    """Return a conjugation POC structure path."""
    return (
        Path(__file__).parents[1] / "src" / "polyzymd" / "builders" / "conjugation" / "poc" / name
    )


def test_clean_protein_only_plan_and_writer_normalize_to_chain_a(tmp_path):
    """A protein-only PDB should plan chain A and leave the source file unchanged."""
    structure = tmp_path / "protein_only.pdb"
    original = (
        _pdb_atom(1, "N", "ALA", "B", 1, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "CA", "ALA", "B", 1, 1.5, 0.0, 0.0)
        + "END\n"
    )
    structure.write_text(original)
    output = tmp_path / "protein_only_cleaned.pdb"

    plan = plan_pdb_chain_normalization(structure)
    written = write_normalized_pdb(structure, output, plan)

    assert plan.valid is True
    assert plan.clean is True
    assert plan.protein_residue_count == 1
    assert {action.target_chain for action in plan.actions} == {"A"}
    assert any("chain A" in warning for warning in plan.warnings)
    assert written == output
    assert structure.read_text() == original
    atom_lines = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]
    assert {line[21] for line in atom_lines} == {"A"}
    assert [line[22:26].strip() for line in atom_lines] == ["1", "1"]


def test_multichain_protein_collapse_renumbers_first_seen_residue_order(tmp_path):
    """Protein chains A and B should collapse to chain A with continuous numbering."""
    structure = tmp_path / "two_chain_protein.pdb"
    lines = []
    serial = 1
    for chain_id, residue_total in [("A", 5), ("B", 6)]:
        for residue_number in range(1, residue_total + 1):
            lines.append(
                _pdb_atom(serial, "CA", "ALA", chain_id, residue_number, float(serial), 0.0, 0.0)
            )
            serial += 1
    structure.write_text("".join(lines) + "END\n")
    output = tmp_path / "two_chain_protein_cleaned.pdb"

    plan = plan_pdb_chain_normalization(structure)
    write_normalized_pdb(structure, output, plan)
    atom_lines = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]

    assert plan.valid is True
    assert plan.protein_residue_count == 11
    assert plan.source_chains_collapsed["A"] == ["A", "B"]
    assert {line[21] for line in atom_lines} == {"A"}
    assert [int(line[22:26]) for line in atom_lines] == list(range(1, 12))


def test_noncontiguous_protein_numbering_renumbers_by_first_seen_order(tmp_path):
    """Non-contiguous source residue numbers should normalize to A:1, A:2, A:3."""
    structure = tmp_path / "noncontiguous.pdb"
    structure.write_text(
        _pdb_atom(1, "CA", "ALA", "A", 10, 0.0, 0.0, 0.0)
        + _pdb_atom(2, "CA", "GLY", "A", 20, 1.0, 0.0, 0.0)
        + _pdb_atom(3, "CA", "SER", "B", 5, 2.0, 0.0, 0.0)
        + "END\n"
    )
    output = tmp_path / "noncontiguous_cleaned.pdb"

    plan = plan_pdb_chain_normalization(structure)
    write_normalized_pdb(structure, output, plan)
    atom_lines = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]

    assert plan.valid is True
    assert [action.target_residue_number for action in plan.actions] == [1, 2, 3]
    assert [int(line[22:26]) for line in atom_lines] == [1, 2, 3]
    assert {line[21] for line in atom_lines} == {"A"}


def test_clean_covalent_glycan_plan_and_writer_update_link_records(tmp_path):
    """A covalently attached glycan should normalize to chain C and update LINK."""
    structure = tmp_path / "linked_glycan.pdb"
    link_line = _pdb_link("ND2", "ASN", "A", 10, "C1", "NAG", "B", 1)
    conect_line = "CONECT    1    2\n"
    structure.write_text(
        _pdb_atom(1, "ND2", "ASN", "A", 10, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "C1", "NAG", "B", 1, 1.4, 0.0, 0.0, record_name="HETATM")
        + link_line
        + conect_line
        + "END\n"
    )
    output = tmp_path / "linked_glycan_cleaned.pdb"

    plan = plan_pdb_chain_normalization(inspect_pdb_structure(structure))
    write_normalized_pdb(structure, output, plan)
    output_text = output.read_text()

    assert plan.valid is True
    assert plan.moiety_residue_count == 1
    assert any(
        action.residue_name == "NAG" and action.target_chain == "C" for action in plan.actions
    )
    assert link_line not in output_text
    link_output = next(line for line in output_text.splitlines() if line.startswith("LINK"))
    assert link_output[21] == "A"
    assert link_output[22:26].strip() == "1"
    assert link_output[51] == "C"
    assert link_output[52:56].strip() == "1"
    assert conect_line in output_text
    glycan_line = next(line for line in output_text.splitlines() if line.startswith("HETATM"))
    assert glycan_line[21] == "C"
    assert glycan_line[22:26].strip() == "1"


def test_protein_on_non_a_chain_warns_but_is_valid(tmp_path):
    """Protein residues outside chain A should be warning-only for clean input."""
    structure = tmp_path / "protein_chain_b.pdb"
    structure.write_text(_pdb_atom(1, "CA", "LYS", "B", 5, 0.0, 0.0, 0.0) + "END\n")

    plan = plan_pdb_chain_normalization(structure)

    assert plan.valid is True
    assert plan.issues == []
    assert plan.actions[0].source_chain == "B"
    assert plan.actions[0].target_chain == "A"
    assert any("chain A" in warning for warning in plan.warnings)


def test_blank_chains_are_valid_for_single_covalent_entity(tmp_path):
    """Blank chains should normalize when all residues belong to one covalent entity."""
    structure = tmp_path / "blank_entity.pdb"
    structure.write_text(
        _pdb_atom(1, "ND2", "ASN", "", 10, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "C1", "NAG", "", 1, 1.4, 0.0, 0.0, record_name="HETATM")
        + "CONECT    1    2\n"
        + "END\n"
    )
    output = tmp_path / "blank_entity_cleaned.pdb"

    plan = plan_pdb_chain_normalization(structure)
    write_normalized_pdb(structure, output, plan)
    chains_by_residue = {
        line[17:20].strip(): line[21]
        for line in output.read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    }

    assert plan.valid is True
    assert any("Blank chain IDs" in warning for warning in plan.warnings)
    assert chains_by_residue == {"ASN": "A", "NAG": "C"}


def test_dirty_pdb_with_water_ion_and_free_ligand_is_invalid_and_not_written(tmp_path):
    """Dirty extra components should invalidate the plan and block the writer."""
    structure = tmp_path / "dirty.pdb"
    output = tmp_path / "dirty_cleaned.pdb"
    structure.write_text(
        _pdb_atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)
        + _pdb_atom(2, "O", "HOH", "D", 1, 5.0, 0.0, 0.0, record_name="HETATM", element="O")
        + _pdb_atom(3, "NA", "NA", "D", 2, 6.0, 0.0, 0.0, record_name="HETATM", element="NA")
        + _pdb_atom(4, "C1", "CIT", "B", 3, 7.0, 0.0, 0.0, record_name="HETATM")
        + "END\n"
    )

    plan = plan_pdb_chain_normalization(structure)

    assert plan.valid is False
    assert {issue.code for issue in plan.issues} >= {
        "dirty_component",
        "free_ligand_or_unlinked_component",
    }
    with pytest.raises(ValueError, match="clean-PDB validation failed"):
        write_normalized_pdb(structure, output, plan)
    assert not output.exists()


def test_default_cleaned_pdb_path_uses_cleaned_suffix(tmp_path):
    """Default normalized output paths should use the approved cleaned suffix."""
    structure = tmp_path / "raw_model.pdb"

    assert default_cleaned_pdb_path(structure) == tmp_path / "raw_model_cleaned.pdb"


def test_glycan_like_residue_without_covalent_evidence_is_invalid(tmp_path):
    """A glycan-like residue requires explicit LINK/CONECT evidence to protein."""
    structure = tmp_path / "unlinked_glycan.pdb"
    structure.write_text(
        _pdb_atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)
        + _pdb_atom(2, "C1", "NAG", "C", 1, 8.0, 0.0, 0.0, record_name="HETATM")
        + "END\n"
    )

    plan = plan_pdb_chain_normalization(structure)

    assert plan.valid is False
    assert any(issue.code == "missing_covalent_evidence" for issue in plan.issues)
    assert any("LINK/CONECT" in issue.message for issue in plan.issues)


def test_poc_1lyz_normalization_plan_reports_actionable_dirty_or_blank_issues():
    """The 1LYZ conjugation POC should plan or fail with clear diagnostics."""
    structure = _poc_path("1LYZ_conj.pdb")
    if not structure.exists():
        pytest.skip("POC conjugation structure is not present")

    plan = plan_pdb_chain_normalization(structure)

    assert plan.protein_residue_count > 0 or any(
        issue.code == "no_protein_residues" for issue in plan.issues
    )
    assert plan.warnings or plan.issues
    combined_messages = " ".join([*plan.warnings, *(issue.message for issue in plan.issues)])
    assert any(
        term in combined_messages for term in ["Blank chain IDs", "Clean-PDB", "LINK/CONECT"]
    )


def test_poc_5fyj_normalization_plan_returns_actionable_result_without_crashing():
    """The glycosylated 5FYJ POC should return a valid or actionable invalid plan."""
    structure = _poc_path("5fyj-monomer-threeGlycans.pdb")
    if not structure.exists():
        pytest.skip("POC conjugation structure is not present")

    plan = plan_pdb_chain_normalization(structure)

    assert plan.protein_residue_count > 0
    assert plan.actions or plan.issues
    if not plan.valid:
        assert plan.issues
        assert plan.output_recommendations
