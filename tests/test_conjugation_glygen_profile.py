"""Tests for GlyGen/GlyCAM residue-resolved glycan PDB ingestion."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from polyzymd.builders.conjugation import ConjugatedPolymerSystemSettings as PublicSettings
from polyzymd.builders.conjugation._moiety_provider import validate_moiety_source_config
from polyzymd.builders.conjugation._pdb_fragment import load_pdb_fragment
from polyzymd.builders.conjugation.reactions.n_glycosylation import glygen_pdb_profile_from_fragment
from polyzymd.builders.conjugation.structure.parsing import (
    parse_pdb_atom_records,
    parse_pdb_conect_pairs,
)
from polyzymd.builders.conjugation.system_workflow import (
    ConjugatedPolymerSystemSettings,
    _build_attachment_spec,
    _build_pdb_fragment_coordinate_only_result,
    _is_protein_fragment_pair,
    _pdb_fragment_final_clash_graph_bonds,
    _summarize_pdb_fragment_true_nonbonded_contacts,
)
from polyzymd.builders.conjugation.validation import (
    classify_bond_path_length,
    summarize_nonbonded_heavy_clashes,
)


def test_glygen_loader_uses_conect_and_preserves_residues(tmp_path: Path) -> None:
    """GlyGen loader should validate CONECT graphs and strict ROH leaving atoms."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)

    load_result = load_pdb_fragment(glycan_path)
    result = glygen_pdb_profile_from_fragment(load_result)

    assert load_result.connectivity_provenance == "conect"
    assert result.reducing_c1_serial == 1
    assert result.fragment.leaving_atom_serials == (8, 9)
    assert [residue.residue_name for residue in result.fragment.residues] == ["NAG", "MAN", "ROH"]
    assert any(item.plausible_glycosidic for item in result.linkage_diagnostics)


def test_conjugation_facade_exports_workflow_settings() -> None:
    """Public conjugation facade should expose workflow settings."""
    assert PublicSettings is ConjugatedPolymerSystemSettings


def test_glygen_loader_infers_coordinates_without_conect(tmp_path: Path) -> None:
    """GlyGen loader should report coordinate-inferred provenance without CONECT."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_no_conect.pdb", include_conect=False)

    load_result = load_pdb_fragment(glycan_path)
    result = glygen_pdb_profile_from_fragment(load_result)

    assert load_result.connectivity_provenance == "coordinate_inferred"
    assert result.reducing_c1_serial == 1


def test_glygen_loader_rejects_missing_roh(tmp_path: Path) -> None:
    """GlyGen loader should fail actionably when ROH labels are absent."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad.pdb", include_conect=True)
    glycan_path.write_text(
        glycan_path.read_text(encoding="utf-8").replace("ROH", "BAD"), encoding="utf-8"
    )

    with pytest.raises(ValueError, match="ROH residue"):
        glygen_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_glygen_loader_rejects_wrong_roh_ho1_element(tmp_path: Path) -> None:
    """GlyGen loader should require ROH:HO1 to be a hydrogen atom."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_ho1.pdb", include_conect=True)
    _replace_atom_element(glycan_path, serial=9, element="O")

    with pytest.raises(ValueError, match="HO1 element H"):
        glygen_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_glygen_loader_rejects_wrong_roh_o1_element(tmp_path: Path) -> None:
    """GlyGen loader should require ROH:O1 to be an oxygen atom."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_o1.pdb", include_conect=True)
    _replace_atom_element(glycan_path, serial=8, element="C")

    with pytest.raises(ValueError, match="O1 element O"):
        glygen_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_glygen_loader_rejects_wrong_reducing_c1_element(tmp_path: Path) -> None:
    """GlyGen loader should require the reducing C1 candidate to be carbon."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_c1.pdb", include_conect=True)
    _replace_atom_element(glycan_path, serial=1, element="O")

    with pytest.raises(ValueError, match="unique reducing sugar C1"):
        glygen_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_glygen_loader_rejects_missing_exact_roh_bond(tmp_path: Path) -> None:
    """GlyGen loader should prove the exact ROH:HO1-ROH:O1-C1 path."""
    glycan_path = _write_glygen_fixture(
        tmp_path / "glycan_missing_roh_bond.pdb", include_conect=True
    )
    text = glycan_path.read_text(encoding="utf-8")
    text = text.replace("CONECT    8    9\n", "CONECT    1    9\n")
    glycan_path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="exact ROH:HO1-ROH:O1-C1"):
        glygen_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_glygen_loader_rejects_extra_roh_o1_carbon_bond(tmp_path: Path) -> None:
    """GlyGen loader should not accept unrelated C-O bonds around ROH:O1."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_extra_roh_bond.pdb", include_conect=True)
    text = glycan_path.read_text(encoding="utf-8")
    text = text.replace("CONECT    8    9\n", "CONECT    8    9    3\n")
    glycan_path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="bonded exactly"):
        glygen_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_coordinate_inferred_loader_rejects_false_hydrogen_proximity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Coordinate-inferred validation should reject impossible H-H bonds."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_false_hh.pdb", include_conect=False)
    _patch_coordinate_bonds(
        monkeypatch,
        atom_count=9,
        bonds=((1, 2), (1, 3), (1, 8), (3, 4), (4, 6), (6, 5), (6, 7), (8, 9), (9, 9)),
    )

    with pytest.raises(ValueError, match="unsafe element bond H-H"):
        glygen_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_coordinate_inferred_loader_rejects_overbonded_atoms(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Coordinate-inferred validation should reject unsafe extra valence."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_overbonded.pdb", include_conect=False)
    _patch_coordinate_bonds(
        monkeypatch,
        atom_count=9,
        bonds=((1, 2), (1, 3), (1, 8), (3, 4), (4, 6), (6, 5), (6, 7), (8, 9), (3, 8)),
    )

    with pytest.raises(ValueError, match="overbonded atoms"):
        glygen_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_coordinate_inferred_loader_accepts_glygen_branch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Coordinate-inferred validation should accept a plausible branched glycan graph."""
    glycan_path = _write_glygen_branch_fixture(tmp_path / "glycan_branch.pdb")
    _patch_coordinate_bonds(
        monkeypatch,
        atom_count=13,
        bonds=(
            (1, 2),
            (1, 3),
            (1, 8),
            (3, 4),
            (3, 10),
            (4, 6),
            (6, 5),
            (6, 7),
            (8, 9),
            (10, 11),
            (11, 12),
            (12, 13),
        ),
    )

    load_result = load_pdb_fragment(glycan_path)
    result = glygen_pdb_profile_from_fragment(load_result)

    assert load_result.connectivity_provenance == "coordinate_inferred"
    assert [residue.residue_name for residue in result.fragment.residues] == [
        "NAG",
        "MAN",
        "ROH",
        "GAL",
    ]


def test_input_path_preserves_generic_exactly_one_source_semantics() -> None:
    """Input-path moieties should be accepted as a generic source at config validation."""
    moiety = SimpleNamespace(
        input_path=Path("glycan.pdb"), smiles=None, residue_name=None, polymer_recipe=None
    )

    assert validate_moiety_source_config(moiety, mechanism_name="nhs_lys") == ["input_path"]
    assert validate_moiety_source_config(moiety, mechanism_name="n_glycosylation") == ["input_path"]


def test_coordinate_only_workflow_removes_roh_and_links_asn(tmp_path: Path) -> None:
    """Coordinate-only workflow should write residue-preserved Asn-glycan artifacts."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    attachment = _attachment(glycan_path)
    settings = ConjugatedPolymerSystemSettings(pdb_fragment_output_mode="coordinate_only")
    packmol_calls: list[Path] = []
    spec, *_ = _build_attachment_spec(
        attachment,
        attachment_index=1,
        protein_pdb_path=protein_path,
        artifact_dir=tmp_path,
        workflow_settings=settings,
    )

    result = _build_pdb_fragment_coordinate_only_result(
        protein_pdb_path=protein_path,
        specs=(spec,),
        output_dir=tmp_path,
        construction_dir=tmp_path / "construction",
        protein_canonicalization=None,
        placement_settings=settings.placement,
        run_packmol_func=_fake_packmol_executor(packmol_calls),
    )
    output = result.crosslinked_conjugate_pdb_path.read_text(encoding="utf-8")

    assert result.status == "coordinate_only"
    assert " ROH " not in output
    assert " NAG C" in output
    assert " MAN C" in output
    assert "CONECT" in output
    assert " ASX A" in output
    assert "HD21" not in output
    assert Path(result.artifact_paths["pdb_fragment_pdb_fragment_ingestion"]).exists()
    assert packmol_calls == [tmp_path / "construction" / "packmol_modifier_placement"]
    assert result.construction.placement.packmol_input_path.exists()
    assert result.construction.placement.packmol_input_text.count("inside sphere") == 1
    assert "atoms 1\n  inside sphere" in result.construction.placement.packmol_input_text
    assert "outside sphere" not in result.construction.placement.packmol_input_text


@pytest.mark.parametrize("external_like", [False, True])
def test_coordinate_only_places_glycan_at_n_glycosylation_bond_length(
    tmp_path: Path, external_like: bool
) -> None:
    """Coordinate-only GlyGen output should use Packmol-constrained placement."""
    glycan_path = _write_glygen_fixture(
        tmp_path / f"glycan_{external_like}.pdb",
        include_conect=True,
        coordinate_offset=(25.0, -7.0, 3.0) if external_like else (0.0, 0.0, 0.0),
    )
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    result, spec = _build_coordinate_only_fixture_result(tmp_path, protein_path, glycan_path)

    output_atoms = tuple(parse_pdb_atom_records(result.crosslinked_conjugate_pdb_path))
    nd2 = _single_atom(output_atoms, chain_id="A", residue_number=1, atom_name="ND2")
    c1 = _single_atom(output_atoms, chain_id="C", residue_number=1, atom_name="C1")
    bond_length = _distance(_xyz(nd2), _xyz(c1))

    assert 1.35 <= bond_length <= 1.65
    assert not np.allclose(_xyz(c1), _xyz(spec.resolved_plan.modifier_link_atom), atol=1.0e-3)


def test_coordinate_only_preserves_internal_glycan_distances(tmp_path: Path) -> None:
    """Rigid coordinate-only placement should not distort glycan internal geometry."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    result, _spec = _build_coordinate_only_fixture_result(tmp_path, protein_path, glycan_path)

    source_atoms = tuple(
        atom for atom in parse_pdb_atom_records(glycan_path) if atom.residue_name.upper() != "ROH"
    )
    output_atoms = tuple(
        atom
        for atom in parse_pdb_atom_records(result.crosslinked_conjugate_pdb_path)
        if atom.chain_id == "C"
    )
    output_by_source = _output_glycan_atoms_by_source(source_atoms, output_atoms)

    for index, left in enumerate(source_atoms):
        for right in source_atoms[index + 1 :]:
            source_distance = _distance(_xyz(left), _xyz(right))
            output_distance = _distance(
                _xyz(output_by_source[left.serial]), _xyz(output_by_source[right.serial])
            )
            assert output_distance == pytest.approx(source_distance, abs=1.0e-6)


def test_coordinate_only_has_no_heavy_nonbonded_clashes(tmp_path: Path) -> None:
    """Coordinate-only output should avoid severe protein-glycan heavy clashes."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    result, _spec = _build_coordinate_only_fixture_result(tmp_path, protein_path, glycan_path)

    output_atoms = tuple(parse_pdb_atom_records(result.crosslinked_conjugate_pdb_path))
    bonds = _pdb_fragment_final_clash_graph_bonds(
        result.crosslinked_conjugate_pdb_path, output_atoms
    )
    summary = summarize_nonbonded_heavy_clashes(
        output_atoms,
        bonds,
        cutoff_angstrom=2.0,
        excluded_bond_depth=3,
        include_pair=_is_protein_fragment_pair,
    )

    assert summary.contact_count == 0
    assert result.construction.placement.true_nonbonded_heavy_contact_count_below_2_angstrom == 0
    assert result.construction.placement.final_conect_graph_valid is True


def test_coordinate_only_retries_code_173_when_output_is_missing(tmp_path: Path) -> None:
    """Packmol code 173 without an output PDB should retry without FileNotFoundError."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    calls: list[Path] = []
    valid_executor = _fake_packmol_executor(calls)

    def fake_packmol(input_text: str, work_dir: Path | str) -> Path:
        """Skip the first best-effort output, then write a valid retry output."""
        working_directory = Path(work_dir)
        if not calls:
            calls.append(working_directory)
            (working_directory / "packmol_error.log").write_text(
                "Packmol ended with code 173\n", encoding="utf-8"
            )
            return working_directory / "packmol_output.pdb"
        return valid_executor(input_text, working_directory)

    result, _spec = _build_coordinate_only_fixture_result_with_executor(
        tmp_path,
        protein_path,
        glycan_path,
        fake_packmol,
    )

    assert calls == [
        tmp_path / "construction" / "packmol_modifier_placement",
        tmp_path / "construction" / "packmol_modifier_placement_attempt_02",
    ]
    assert result.construction.placement.final_conect_graph_valid is True


def test_coordinate_only_accepts_valid_code_173_packmol_output(tmp_path: Path) -> None:
    """Packmol code 173 should remain usable when the output PDB validates."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    calls: list[Path] = []
    valid_executor = _fake_packmol_executor(calls)

    def fake_packmol(input_text: str, work_dir: Path | str) -> Path:
        """Write a valid best-effort output and the retained Packmol error log."""
        output_path = valid_executor(input_text, work_dir)
        (Path(work_dir) / "packmol_error.log").write_text(
            "Packmol ended with code 173\n", encoding="utf-8"
        )
        return output_path

    result, _spec = _build_coordinate_only_fixture_result_with_executor(
        tmp_path,
        protein_path,
        glycan_path,
        fake_packmol,
    )

    assert calls == [tmp_path / "construction" / "packmol_modifier_placement"]
    assert result.construction.placement.packmol_exit_status == "173_imperfect_accepted"
    assert result.construction.placement.final_conect_graph_valid is True


def test_coordinate_only_classifies_glycosidic_neighbors_by_graph_distance(
    tmp_path: Path,
) -> None:
    """Close ND2-O5 and CG-O5 contacts should classify as near-neighbor geometry."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    result, _spec = _build_coordinate_only_fixture_result(tmp_path, protein_path, glycan_path)

    output_atoms = tuple(parse_pdb_atom_records(result.crosslinked_conjugate_pdb_path))
    bonds = _pdb_fragment_final_clash_graph_bonds(
        result.crosslinked_conjugate_pdb_path, output_atoms
    )
    nd2 = _single_atom(output_atoms, chain_id="A", residue_number=1, atom_name="ND2")
    cg = _single_atom(output_atoms, chain_id="A", residue_number=1, atom_name="CG")
    o5 = _single_atom(output_atoms, chain_id="C", residue_number=1, atom_name="O5")
    summary = _summarize_pdb_fragment_true_nonbonded_contacts(result.crosslinked_conjugate_pdb_path)

    assert classify_bond_path_length(nd2.serial, o5.serial, bonds) == 2
    assert classify_bond_path_length(cg.serial, o5.serial, bonds) == 3
    assert summary.contact_count == 0


def test_coordinate_only_conect_matches_retained_source_graph_and_crosslink(
    tmp_path: Path,
) -> None:
    """Output CONECT should equal retained source graph plus the Asn-glycan crosslink."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    result, _spec = _build_coordinate_only_fixture_result(tmp_path, protein_path, glycan_path)

    source_atoms = tuple(parse_pdb_atom_records(glycan_path))
    output_atoms = tuple(parse_pdb_atom_records(result.crosslinked_conjugate_pdb_path))
    retained_source_atoms = tuple(
        atom for atom in source_atoms if atom.residue_name.upper() != "ROH"
    )
    output_glycan_atoms = tuple(atom for atom in output_atoms if atom.chain_id == "C")
    output_by_source = _output_glycan_atoms_by_source(retained_source_atoms, output_glycan_atoms)
    source_to_output = {serial: atom.serial for serial, atom in output_by_source.items()}

    retained_source_serials = set(source_to_output)
    expected_edges = {
        frozenset((source_to_output[left], source_to_output[right]))
        for left, right in parse_pdb_conect_pairs(glycan_path)
        if left in retained_source_serials and right in retained_source_serials
    }
    nd2 = _single_atom(output_atoms, chain_id="A", residue_number=1, atom_name="ND2")
    c1 = _single_atom(output_atoms, chain_id="C", residue_number=1, atom_name="C1")
    expected_edges.add(frozenset((nd2.serial, c1.serial)))
    output_edges = {
        frozenset(edge) for edge in parse_pdb_conect_pairs(result.crosslinked_conjugate_pdb_path)
    }

    assert output_edges == expected_edges


def test_coordinate_only_emits_no_implausible_conect_lengths(tmp_path: Path) -> None:
    """Coordinate-only CONECT records should not contain long covalent edges."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    result, _spec = _build_coordinate_only_fixture_result(tmp_path, protein_path, glycan_path)

    atoms_by_serial = {
        atom.serial: atom for atom in parse_pdb_atom_records(result.crosslinked_conjugate_pdb_path)
    }
    long_edges = []
    for left, right in parse_pdb_conect_pairs(result.crosslinked_conjugate_pdb_path):
        atom_left = atoms_by_serial[left]
        atom_right = atoms_by_serial[right]
        distance = _distance(_xyz(atom_left), _xyz(atom_right))
        has_hydrogen = "H" in {_element(atom_left), _element(atom_right)}
        limit = 1.405 if has_hydrogen else 2.205
        if distance > limit:
            long_edges.append((left, right, distance))

    assert long_edges == []


def test_glygen_loader_converts_coordinate_inferred_indices_to_source_serials(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Coordinate-inferred GlyGen bonds should not retain raw RDKit atom indices."""
    glycan_path = _write_glygen_fixture(
        tmp_path / "glycan_high_serials.pdb",
        include_conect=False,
        serial_offset=100,
    )
    _patch_coordinate_bonds(
        monkeypatch,
        atom_count=9,
        bonds=(
            (1, 2),
            (1, 3),
            (1, 8),
            (3, 4),
            (4, 6),
            (6, 5),
            (6, 7),
            (8, 9),
        ),
    )

    result = glygen_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))

    assert result.fragment.bonds[0] == (101, 102)
    assert all(endpoint >= 101 for bond in result.fragment.bonds for endpoint in bond)


def test_coordinate_only_high_serial_graph_avoids_raw_index_endpoint_ambiguity(
    tmp_path: Path,
) -> None:
    """High source serials guard against writer serial-first/index-fallback ambiguity."""
    glycan_path = _write_glygen_fixture(
        tmp_path / "glycan_high_serials.pdb",
        include_conect=True,
        serial_offset=100,
    )
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    result, spec = _build_coordinate_only_fixture_result(tmp_path, protein_path, glycan_path)

    assert spec.generated_fragment.bonds[0] == (101, 102)
    assert all(endpoint >= 101 for bond in spec.generated_fragment.bonds for endpoint in bond)

    output_edges = {
        frozenset(edge) for edge in parse_pdb_conect_pairs(result.crosslinked_conjugate_pdb_path)
    }
    output_atoms = tuple(parse_pdb_atom_records(result.crosslinked_conjugate_pdb_path))
    output_glycan_atoms = tuple(atom for atom in output_atoms if atom.chain_id == "C")
    source_atoms = tuple(
        atom for atom in parse_pdb_atom_records(glycan_path) if atom.residue_name.upper() != "ROH"
    )
    output_by_source = _output_glycan_atoms_by_source(source_atoms, output_glycan_atoms)
    expected_c1_o5 = frozenset((output_by_source[101].serial, output_by_source[102].serial))

    assert expected_c1_o5 in output_edges


def _attachment(glycan_path: Path) -> SimpleNamespace:
    """Build a minimal n_glycosylation attachment namespace."""
    return SimpleNamespace(
        name="asn_glycan",
        enabled=True,
        site=SimpleNamespace(
            chain_id="A",
            residue_name="ASN",
            residue_number=1,
            atom_name="ND2",
            insertion_code="",
            atom_serial=None,
            atom_index=None,
        ),
        moiety=SimpleNamespace(
            name="glycan",
            input_path=glycan_path,
            smiles=None,
            residue_name=None,
            polymer_recipe=None,
        ),
        mechanism=SimpleNamespace(
            name="n_glycosylation",
            product_residues=SimpleNamespace(site=None, moiety=None),
            bond=SimpleNamespace(order=1, site_atom=None, target_bond_length_angstrom=1.45),
        ),
    )


def _build_coordinate_only_fixture_result(
    tmp_path: Path,
    protein_path: Path,
    glycan_path: Path,
) -> tuple[Any, Any]:
    """Build a coordinate-only GlyGen fixture result and its attachment spec."""
    packmol_calls: list[Path] = []
    result, spec = _build_coordinate_only_fixture_result_with_executor(
        tmp_path,
        protein_path,
        glycan_path,
        _fake_packmol_executor(packmol_calls),
    )
    assert packmol_calls == [tmp_path / "construction" / "packmol_modifier_placement"]
    return result, spec


def _build_coordinate_only_fixture_result_with_executor(
    tmp_path: Path,
    protein_path: Path,
    glycan_path: Path,
    run_packmol_func: Any,
) -> tuple[Any, Any]:
    """Build a coordinate-only GlyGen fixture result with a custom Packmol fake."""
    attachment = _attachment(glycan_path)
    settings = ConjugatedPolymerSystemSettings(pdb_fragment_output_mode="coordinate_only")
    spec, *_ = _build_attachment_spec(
        attachment,
        attachment_index=1,
        protein_pdb_path=protein_path,
        artifact_dir=tmp_path,
        workflow_settings=settings,
    )
    result = _build_pdb_fragment_coordinate_only_result(
        protein_pdb_path=protein_path,
        specs=(spec,),
        output_dir=tmp_path,
        construction_dir=tmp_path / "construction",
        protein_canonicalization=None,
        placement_settings=settings.placement,
        run_packmol_func=run_packmol_func,
    )
    return result, spec


def _fake_packmol_executor(calls: list[Path]) -> Any:
    """Return a Packmol fake that writes fixed protein plus translated modifier coordinates."""

    def fake_packmol(input_text: str, work_dir: Path | str) -> Path:
        """Write deterministic Packmol-like output and record invocation."""
        working_directory = Path(work_dir)
        calls.append(working_directory)
        structure_paths = _packmol_structure_paths(input_text)
        assert len(structure_paths) == 2
        protein_atoms = tuple(parse_pdb_atom_records(structure_paths[0]))
        modifier_atoms = tuple(parse_pdb_atom_records(structure_paths[1]))
        sphere_center = _packmol_inside_sphere_center(input_text)
        reactive_index = _packmol_reactive_atom_index(input_text)
        reactive_atom = modifier_atoms[reactive_index]
        translation = sphere_center + np.asarray([3.0, 0.0, 0.0]) - _xyz(reactive_atom)
        output_path = working_directory / "packmol_output.pdb"
        lines = []
        serial = 1
        for atom in protein_atoms:
            lines.append(
                _pdb_atom(serial, atom.atom_name, "PRO", "A", 1, atom.x, atom.y, atom.z, "C")
            )
            serial += 1
        for atom in modifier_atoms:
            xyz = _xyz(atom) + translation
            lines.append(
                _pdb_atom(
                    serial,
                    atom.atom_name,
                    atom.residue_name,
                    "C",
                    atom.residue_number,
                    float(xyz[0]),
                    float(xyz[1]),
                    float(xyz[2]),
                    _element(atom),
                )
            )
            serial += 1
        output_path.write_text("".join(lines) + "END\n", encoding="utf-8")
        return output_path

    return fake_packmol


def _packmol_structure_paths(input_text: str) -> tuple[Path, ...]:
    """Extract structure paths from Packmol input text."""
    return tuple(
        Path(line.split(maxsplit=1)[1])
        for line in input_text.splitlines()
        if line.startswith("structure ")
    )


def _packmol_inside_sphere_center(input_text: str) -> np.ndarray:
    """Extract the reactive-site sphere center from Packmol input text."""
    for line in input_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("inside sphere "):
            fields = stripped.split()
            return np.asarray([float(fields[2]), float(fields[3]), float(fields[4])])
    raise AssertionError("Packmol input did not contain an inside sphere constraint")


def _packmol_reactive_atom_index(input_text: str) -> int:
    """Extract the zero-based reactive atom index from Packmol input text."""
    for line in input_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("atoms "):
            return int(stripped.split()[1]) - 1
    raise AssertionError("Packmol input did not contain an atom constraint")


def _single_atom(
    atoms: tuple[Any, ...], *, chain_id: str, residue_number: int, atom_name: str
) -> Any:
    """Return one atom matching simple PDB identity fields."""
    matches = [
        atom
        for atom in atoms
        if atom.chain_id == chain_id
        and atom.residue_number == residue_number
        and atom.atom_name.strip().upper() == atom_name.upper()
    ]
    assert len(matches) == 1
    return matches[0]


def _output_glycan_atoms_by_source(
    source_atoms: tuple[Any, ...], output_atoms: tuple[Any, ...]
) -> dict[int, Any]:
    """Map retained source glycan serials to output atoms by residue and atom name."""
    output_by_key = {
        (atom.residue_name, atom.residue_number, atom.atom_name): atom for atom in output_atoms
    }
    mapping = {}
    for atom in source_atoms:
        key = (atom.residue_name, atom.residue_number, atom.atom_name)
        assert key in output_by_key
        mapping[atom.serial] = output_by_key[key]
    return mapping


def _xyz(atom: Any) -> np.ndarray:
    """Return one atom coordinate vector."""
    return np.asarray([atom.x, atom.y, atom.z], dtype=float)


def _distance(left: np.ndarray, right: np.ndarray) -> float:
    """Return Euclidean distance between two coordinate vectors."""
    return float(np.linalg.norm(left - right))


def _element(atom: Any) -> str:
    """Return a normalized atom element symbol."""
    return (atom.element or atom.atom_name[:1]).strip().upper()


def _write_asn_fixture(path: Path) -> Path:
    """Write a tiny explicit-hydrogen ASN residue fixture."""
    atoms = [
        (1, "N", "ASN", "A", 1, 0.0, 0.0, 0.0, "N"),
        (2, "CA", "ASN", "A", 1, 1.4, 0.0, 0.0, "C"),
        (3, "C", "ASN", "A", 1, 2.0, 1.2, 0.0, "C"),
        (4, "CB", "ASN", "A", 1, 1.8, -1.2, 0.0, "C"),
        (5, "CG", "ASN", "A", 1, 3.0, -1.2, 0.0, "C"),
        (6, "OD1", "ASN", "A", 1, 3.6, -2.2, 0.0, "O"),
        (7, "ND2", "ASN", "A", 1, 3.5, 0.0, 0.0, "N"),
        (8, "HD21", "ASN", "A", 1, 4.4, 0.0, 0.0, "H"),
        (9, "HD22", "ASN", "A", 1, 3.1, 0.8, 0.0, "H"),
    ]
    path.write_text("".join(_pdb_atom(*atom) for atom in atoms) + "END\n", encoding="utf-8")
    return path


def _write_glygen_fixture(
    path: Path,
    *,
    include_conect: bool,
    coordinate_offset: tuple[float, float, float] = (0.0, 0.0, 0.0),
    serial_offset: int = 0,
) -> Path:
    """Write a small labeled multi-residue GlyGen-style glycan fixture."""
    dx, dy, dz = coordinate_offset
    atoms = [
        (1, "C1", "NAG", "", 1, 0.0 + dx, 0.0 + dy, 0.0 + dz, "C"),
        (2, "O5", "NAG", "", 1, 0.0 + dx, 1.4 + dy, 0.0 + dz, "O"),
        (3, "C2", "NAG", "", 1, 1.5 + dx, 0.0 + dy, 0.0 + dz, "C"),
        (4, "O2", "NAG", "", 1, 2.8 + dx, 0.0 + dy, 0.0 + dz, "O"),
        (5, "O1", "MAN", "", 2, 5.2 + dx, 0.0 + dy, 0.0 + dz, "O"),
        (6, "C1", "MAN", "", 2, 4.1 + dx, 0.0 + dy, 0.0 + dz, "C"),
        (7, "O5", "MAN", "", 2, 4.1 + dx, 1.4 + dy, 0.0 + dz, "O"),
        (8, "O1", "ROH", "", 3, -1.3 + dx, 0.0 + dy, 0.0 + dz, "O"),
        (9, "HO1", "ROH", "", 3, -2.2 + dx, 0.0 + dy, 0.0 + dz, "H"),
    ]
    lines = [_pdb_atom(atom[0] + serial_offset, *atom[1:]) for atom in atoms]
    if include_conect:
        serial = {index: index + serial_offset for index in range(1, 10)}
        lines.extend(
            [
                _conect(serial[1], serial[2], serial[3], serial[8]),
                _conect(serial[3], serial[4]),
                _conect(serial[4], serial[6]),
                _conect(serial[6], serial[5]),
                _conect(serial[6], serial[7]),
                _conect(serial[8], serial[9]),
            ]
        )
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def _conect(source: int, *targets: int) -> str:
    """Format one fixed-width CONECT fixture line."""
    return f"CONECT{source:5d}{''.join(f'{target:5d}' for target in targets)}\n"


def _write_glygen_branch_fixture(path: Path) -> Path:
    """Write a labeled branched GlyGen-style glycan fixture without CONECT."""
    atoms = [
        (1, "C1", "NAG", "", 1, 0.0, 0.0, 0.0, "C"),
        (2, "O5", "NAG", "", 1, 0.0, 1.4, 0.0, "O"),
        (3, "C2", "NAG", "", 1, 1.5, 0.0, 0.0, "C"),
        (4, "O2", "NAG", "", 1, 2.8, 0.0, 0.0, "O"),
        (5, "O1", "MAN", "", 2, 5.2, 0.0, 0.0, "O"),
        (6, "C1", "MAN", "", 2, 4.1, 0.0, 0.0, "C"),
        (7, "O5", "MAN", "", 2, 4.1, 1.4, 0.0, "O"),
        (8, "O1", "ROH", "", 3, -1.3, 0.0, 0.0, "O"),
        (9, "HO1", "ROH", "", 3, -2.2, 0.0, 0.0, "H"),
        (10, "C3", "NAG", "", 1, 1.5, 1.4, 0.0, "C"),
        (11, "O3", "NAG", "", 1, 2.8, 1.4, 0.0, "O"),
        (12, "C1", "GAL", "", 4, 3.9, 1.4, 0.0, "C"),
        (13, "O5", "GAL", "", 4, 3.9, 2.8, 0.0, "O"),
    ]
    path.write_text("".join(_pdb_atom(*atom) for atom in atoms) + "END\n", encoding="utf-8")
    return path


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    """Format one fixed-width PDB atom line for fixtures."""
    return (
        f"HETATM{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )


def _replace_atom_element(path: Path, *, serial: int, element: str) -> None:
    """Replace one fixture atom element while preserving fixed-width columns."""
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")) and int(line[6:11]) == serial:
            line = f"{line[:76]}{element:>2}{line[78:]}"
        lines.append(line)
    path.write_text("".join(lines), encoding="utf-8")


def _patch_coordinate_bonds(
    monkeypatch: pytest.MonkeyPatch,
    *,
    atom_count: int,
    bonds: tuple[tuple[int, int], ...],
) -> None:
    """Patch RDKit proximity bonding to make graph validation deterministic."""
    from rdkit import Chem

    class FakeBond:
        """Minimal RDKit bond double for serial-order fixtures."""

        def __init__(self, begin_serial: int, end_serial: int) -> None:
            self._begin_index = begin_serial - 1
            self._end_index = end_serial - 1

        def GetBeginAtomIdx(self) -> int:
            """Return the zero-based begin atom index."""
            return self._begin_index

        def GetEndAtomIdx(self) -> int:
            """Return the zero-based end atom index."""
            return self._end_index

    class FakeMol:
        """Minimal RDKit molecule double for coordinate inference."""

        def GetNumAtoms(self) -> int:
            """Return the fixture atom count."""
            return atom_count

        def GetBonds(self) -> tuple[FakeBond, ...]:
            """Return deterministic fake bonds."""
            return tuple(FakeBond(left, right) for left, right in bonds)

    def fake_mol_from_pdb_file(*_args: Any, **_kwargs: Any) -> FakeMol:
        """Return the fake molecule for the loader's RDKit call."""
        return FakeMol()

    monkeypatch.setattr(Chem, "MolFromPDBFile", fake_mol_from_pdb_file)
