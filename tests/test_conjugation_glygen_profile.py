"""Tests for residue-resolved glycan PDB ingestion."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from polyzymd.builders.conjugation import ConjugatedPolymerSystemSettings as PublicSettings
from polyzymd.builders.conjugation._linkage import PdbAtomSelector
from polyzymd.builders.conjugation._moiety_provider import (
    resolve_moiety_source,
    validate_moiety_source_config,
)
from polyzymd.builders.conjugation._pdb_fragment import load_pdb_fragment
from polyzymd.builders.conjugation.reactions.n_glycosylation import (
    NGlycosylationReaction,
    _resolve_asn_nd2_hydrogen,
    residue_resolved_glycan_pdb_profile_from_fragment,
)
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

G42666_CONECT_PATH = (
    Path(__file__).resolve().parent / "data/conjugation/glygen/synthetic_g42666_style_conect.pdb"
)


def test_glygen_loader_uses_conect_and_preserves_residues(tmp_path: Path) -> None:
    """GlyGen loader should validate CONECT graphs and strict ROH leaving atoms."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)

    load_result = load_pdb_fragment(glycan_path)
    result = residue_resolved_glycan_pdb_profile_from_fragment(load_result)

    assert load_result.connectivity_provenance == "conect"
    assert result.reducing_c1_serial == 4
    assert result.fragment.leaving_atom_serials == (1, 2)
    assert [residue.residue_name for residue in result.fragment.residues] == ["ROH", "NAG", "MAN"]
    assert any(item.plausible_glycosidic for item in result.linkage_diagnostics)


def test_external_g42666_conect_exact_anomeric_group_and_plan() -> None:
    """Local external G42666 fixture should resolve the audited reducing-end serials."""
    if not G42666_CONECT_PATH.exists():
        pytest.skip(f"External G42666 fixture not present: {G42666_CONECT_PATH}")
    load_result = load_pdb_fragment(G42666_CONECT_PATH)
    profile = residue_resolved_glycan_pdb_profile_from_fragment(load_result)

    assert profile.reducing_c1_serial == 4
    assert profile.fragment.leaving_atom_serials == (1, 2)
    assert profile.fragment.atoms[profile.hydroxyl_oxygen_atom_index].serial == 1
    assert profile.fragment.atoms[profile.hydroxyl_hydrogen_atom_index].serial == 2
    assert profile.retained_ring_oxygen_serial == 14
    assert profile.fragment.reactive_atom_serial == 4


def test_external_g42666_bond_orders_preserve_graph_and_carbonyls() -> None:
    """Bond-order assignment should preserve G42666 CONECT and repair acetamide C=O."""
    if not G42666_CONECT_PATH.exists():
        pytest.skip(f"External G42666 fixture not present: {G42666_CONECT_PATH}")

    load_result = load_pdb_fragment(G42666_CONECT_PATH)
    bond_set = {frozenset(bond) for bond in load_result.serial_bonds}
    order_set = {frozenset((left, right)) for left, right, _order in load_result.serial_bond_orders}
    orders = {
        frozenset((left, right)): order for left, right, order in load_result.serial_bond_orders
    }

    assert order_set == bond_set
    assert orders[frozenset((6, 8))] == pytest.approx(2.0)
    assert orders[frozenset((33, 35))] == pytest.approx(2.0)
    assert len(load_result.source_atoms) == 57


def test_generated_fragment_equivalent_uses_graph_detector(tmp_path: Path) -> None:
    """Generated fragments should follow the same graph semantics as PDB input."""
    load_result = load_pdb_fragment(
        _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    )
    profile = residue_resolved_glycan_pdb_profile_from_fragment(load_result)

    contract = NGlycosylationReaction.build_contract(
        {"chain_id": "A", "residue_name": "ASN", "residue_number": 1},
        profile.fragment,
    )

    assert contract.modifier_endpoint.selector.atom_serial == 4
    assert tuple(
        selector.atom_serial for selector in contract.modifier_endpoint.leaving_atom_selectors
    ) == (
        1,
        2,
    )


def test_g42666_resolve_plan_removes_one_asn_h_and_keeps_link_atom(tmp_path: Path) -> None:
    """N-glycosylation plan should remove one ND2 H and retain the C1 link atom."""
    if not G42666_CONECT_PATH.exists():
        pytest.skip(f"External G42666 fixture not present: {G42666_CONECT_PATH}")
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    profile = residue_resolved_glycan_pdb_profile_from_fragment(
        load_pdb_fragment(G42666_CONECT_PATH)
    )

    plan = NGlycosylationReaction.resolve_plan(
        protein_path,
        {"chain_id": "A", "residue_name": "ASN", "residue_number": 1},
        profile.fragment,
    )

    assert plan.protein_link_atom.atom_name == "ND2"
    assert tuple(atom.atom_name for atom in plan.protein_leaving_atoms) == ("HD21",)
    assert plan.modifier_link_atom.serial == 4
    assert tuple(atom.serial for atom in plan.modifier_leaving_atoms) == (1, 2)
    assert plan.pablo_crosslink_requirement.linking_atoms == ("ND2", "C1")
    assert plan.pablo_crosslink_requirement.leaving_atoms[0] == ("HD21",)


def test_asn_nd2_hydrogen_resolution_uses_template_not_coordinates(tmp_path: Path) -> None:
    """Pablo ASN topology should select ND2 H even when coordinates are nonbondlike."""
    protein_path = _write_asn_fixture(tmp_path / "asn_far_h.pdb", hydrogens=("HD21",))

    hydrogen = _resolve_asn_nd2_hydrogen(protein_path, _asn_selector())

    assert hydrogen.atom_name == "HD21"


def test_asn_nd2_hydrogen_resolution_rejects_missing_template_h(tmp_path: Path) -> None:
    """Missing Pablo-template ND2 H names should fail without geometry fallback."""
    protein_path = _write_asn_fixture(tmp_path / "asn_no_h.pdb", hydrogens=())

    with pytest.raises(ValueError, match="exactly one explicit Asn ND2 hydrogen"):
        _resolve_asn_nd2_hydrogen(protein_path, _asn_selector())


def test_asn_nd2_hydrogen_resolution_selects_hd21_from_canonical_template_pair(
    tmp_path: Path,
) -> None:
    """Canonical Pablo-template ND2 H names should resolve deterministically to HD21."""
    protein_path = _write_asn_fixture(tmp_path / "asn_two_h.pdb", hydrogens=("HD22", "HD21"))

    hydrogen = _resolve_asn_nd2_hydrogen(protein_path, _asn_selector())

    assert hydrogen.atom_name == "HD21"


def test_g42666_product_pablo_definitions_retain_carbonyl_order_two(tmp_path: Path) -> None:
    """Product-state Pablo residue definitions should keep G42666 carbonyl orders."""
    if not G42666_CONECT_PATH.exists():
        pytest.skip(f"External G42666 fixture not present: {G42666_CONECT_PATH}")
    pytest.importorskip("openff.pablo")
    from polyzymd.builders.conjugation.pablo.product import (
        build_product_state_pablo_library_for_specs,
    )

    glycan_path = G42666_CONECT_PATH
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
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
        run_packmol_func=_fake_packmol_executor([]),
    )

    library = build_product_state_pablo_library_for_specs(
        result.crosslinked_conjugate_pdb_path,
        protein_path,
        (spec,),
    )
    carbonyl_orders = [
        getattr(bond, "order", None)
        for definition in library.definitions
        for bond in getattr(definition, "bonds", ())
        if {getattr(bond, "atom1", ""), getattr(bond, "atom2", "")} == {"C2N", "O2N"}
    ]

    assert carbonyl_orders.count(2) >= 2


def test_product_pablo_endpoint_prefers_assembly_provenance_over_source_serial() -> None:
    """Product-state endpoint lookup should distinguish source and product serials."""
    from polyzymd.builders.conjugation.pablo.product import _locate_residue_key
    from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord

    atoms = [
        _pdb_record(serial=4, atom_name="CB", residue_name="ASN", chain_id="A", residue_number=1),
        _pdb_record(serial=20, atom_name="C1", residue_name="NAG", chain_id="C", residue_number=1),
    ]
    source_link_atom = PdbAtomRecord(
        serial=4,
        atom_name="C1",
        residue_name="NAG",
        chain_id="C",
        residue_number=1,
        x=0.0,
        y=0.0,
        z=0.0,
        element="C",
    )

    key = _locate_residue_key(
        atoms,
        residue_name="NAG",
        atom_name="C1",
        resolved_atom=source_link_atom,
        endpoint_provenance={"conect_pair": {"protein_serial": 1, "modifier_serial": 20}},
        endpoint_role="modifier",
    )

    assert key == ("C", "NAG", 1, "")


def test_product_pablo_mapped_endpoint_rejects_ambiguous_fallback() -> None:
    """Mapped multi-attachment endpoint lookup should fail instead of guessing."""
    from polyzymd.builders.conjugation.pablo.product import _locate_residue_key

    atoms = [
        _pdb_record(serial=10, atom_name="C1", residue_name="NAG", chain_id="C", residue_number=1),
        _pdb_record(serial=40, atom_name="C1", residue_name="NAG", chain_id="C", residue_number=2),
    ]

    with pytest.raises(ValueError, match="requires exact assembly provenance"):
        _locate_residue_key(
            atoms,
            residue_name="NAG",
            atom_name="C1",
            resolved_atom=None,
            allow_legacy_serial_fallback=False,
        )


def test_conjugation_facade_exports_workflow_settings() -> None:
    """Public conjugation facade should expose workflow settings."""
    assert PublicSettings is ConjugatedPolymerSystemSettings


def test_glygen_loader_rejects_missing_conect_without_coordinate_inference(tmp_path: Path) -> None:
    """GlyGen loader should require explicit CONECT and never infer from coordinates."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_no_conect.pdb", include_conect=False)

    with pytest.raises(ValueError, match="requires complete CONECT records"):
        load_pdb_fragment(glycan_path)


def test_loader_is_not_roh_residue_name_dependent(tmp_path: Path) -> None:
    """PDB-fragment compatibility should use graph chemistry instead of ROH names."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad.pdb", include_conect=True)
    glycan_path.write_text(
        glycan_path.read_text(encoding="utf-8").replace("ROH", "BAD"), encoding="utf-8"
    )

    result = residue_resolved_glycan_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))

    assert result.fragment.leaving_atom_serials == (1, 2)


def test_glygen_loader_rejects_wrong_roh_ho1_element(tmp_path: Path) -> None:
    """GlyGen loader should require ROH:HO1 to be a hydrogen atom."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_ho1.pdb", include_conect=True)
    _replace_atom_element(glycan_path, serial=2, element="O")
    glycan_path.write_text(
        glycan_path.read_text(encoding="utf-8").replace("HO1", "OO1"), encoding="utf-8"
    )

    with pytest.raises(
        ValueError,
        match="above obvious upper valence|No glycan anomeric motif|bond orders could not be assigned",
    ):
        residue_resolved_glycan_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_glygen_loader_rejects_wrong_roh_o1_element(tmp_path: Path) -> None:
    """GlyGen loader should require ROH:O1 to be an oxygen atom."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_o1.pdb", include_conect=True)
    _replace_atom_element(glycan_path, serial=1, element="C")

    with pytest.raises(
        ValueError,
        match="above obvious upper valence|No glycan anomeric motif|bond orders could not be assigned",
    ):
        residue_resolved_glycan_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_glygen_loader_rejects_wrong_reducing_c1_element(tmp_path: Path) -> None:
    """GlyGen loader should require the reducing C1 candidate to be carbon."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_c1.pdb", include_conect=True)
    _replace_atom_element(glycan_path, serial=4, element="O")

    with pytest.raises(
        ValueError,
        match="above obvious upper valence|No glycan anomeric motif|bond orders could not be assigned",
    ):
        residue_resolved_glycan_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_glygen_loader_rejects_missing_exact_roh_bond(tmp_path: Path) -> None:
    """GlyGen loader should prove the exact ROH:HO1-ROH:O1-C1 path."""
    glycan_path = _write_glygen_fixture(
        tmp_path / "glycan_missing_roh_bond.pdb", include_conect=True
    )
    text = glycan_path.read_text(encoding="utf-8")
    text = text.replace("CONECT    1    2\n", "CONECT    4    2\n")
    glycan_path.write_text(text, encoding="utf-8")

    with pytest.raises(
        ValueError,
        match="above obvious upper valence|No glycan anomeric motif|bond orders could not be assigned",
    ):
        residue_resolved_glycan_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_glygen_loader_rejects_extra_roh_o1_carbon_bond(tmp_path: Path) -> None:
    """GlyGen loader should not accept unrelated C-O bonds around ROH:O1."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_extra_roh_bond.pdb", include_conect=True)
    text = glycan_path.read_text(encoding="utf-8")
    text = text.replace("CONECT    1    2\n", "CONECT    1    2    3\n")
    glycan_path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="above obvious upper valence|No glycan anomeric motif"):
        residue_resolved_glycan_pdb_profile_from_fragment(load_pdb_fragment(glycan_path))


def test_loader_rejects_unknown_conect_serial(tmp_path: Path) -> None:
    """Strict CONECT loading should reject references outside the atom table."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_unknown_serial.pdb", include_conect=True)
    glycan_path.write_text(
        glycan_path.read_text(encoding="utf-8").replace("CONECT    1    2", "CONECT    1   99"),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="unknown atom serials"):
        load_pdb_fragment(glycan_path)


def test_loader_rejects_self_conect_bond(tmp_path: Path) -> None:
    """Strict CONECT loading should reject self bonds instead of dropping them."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_self_bond.pdb", include_conect=True)
    glycan_path.write_text(
        glycan_path.read_text(encoding="utf-8") + "CONECT    9    9\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="self bonds"):
        load_pdb_fragment(glycan_path)


def test_loader_rejects_isolated_atoms_from_partial_conect(tmp_path: Path) -> None:
    """Strict CONECT loading should reject partial graphs with isolated atoms."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_partial_conect.pdb", include_conect=True)
    text = glycan_path.read_text(encoding="utf-8").replace("CONECT    1    2\n", "")
    glycan_path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="hydrogens must have degree 1|disconnected"):
        load_pdb_fragment(glycan_path)


def test_loader_rejects_explicit_hydrogen_degree_not_one(tmp_path: Path) -> None:
    """Strict CONECT loading should reject hydrogens with degree other than one."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_h_degree.pdb", include_conect=True)
    text = glycan_path.read_text(encoding="utf-8")
    text = text.replace("CONECT    1    2\n", "CONECT    1    2\nCONECT    2    4\n")
    glycan_path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="hydrogens must have degree 1"):
        load_pdb_fragment(glycan_path)


def test_loader_rejects_obvious_upper_valence(tmp_path: Path) -> None:
    """Strict CONECT loading should reject impossible valence before motif detection."""
    glycan_path = _write_glygen_fixture(
        tmp_path / "glycan_overvalent_oxygen.pdb", include_conect=True
    )
    text = glycan_path.read_text(encoding="utf-8").replace(
        "CONECT    1    2\n",
        "CONECT    1    2    3    4\n",
    )
    glycan_path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="above obvious upper valence"):
        load_pdb_fragment(glycan_path)


def test_input_path_preserves_generic_exactly_one_source_semantics() -> None:
    """Input-path moieties should be accepted as a generic source at config validation."""
    moiety = SimpleNamespace(
        input_path=Path("glycan.pdb"), smiles=None, residue_name=None, polymer_recipe=None
    )

    assert validate_moiety_source_config(moiety, mechanism_name="nhs_lys") == ["input_path"]
    assert validate_moiety_source_config(moiety, mechanism_name="n_glycosylation") == ["input_path"]


def test_input_path_provider_routes_pdb_fragment_compatibility_to_template(
    tmp_path: Path,
) -> None:
    """Provider should let the selected reaction template reject unsupported fragments."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    attachment = _attachment(glycan_path)
    attachment.mechanism.name = "nhs_lys_amide"
    attachment.site.residue_name = "LYS"
    attachment.site.atom_name = "NZ"

    with pytest.raises(ValueError, match="does not support PDB-fragment moiety sources"):
        resolve_moiety_source(
            attachment,
            attachment_index=1,
            output_dir=tmp_path,
            protein_pdb_path=tmp_path / "protein.pdb",
        )


def test_n_glycosylation_template_owns_pdb_fragment_profile_resolution(
    tmp_path: Path,
) -> None:
    """N-glycosylation should add its profile only through the reaction hook."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    source = resolve_moiety_source(
        _attachment(glycan_path),
        attachment_index=1,
        output_dir=tmp_path,
        protein_pdb_path=tmp_path / "protein.pdb",
    )

    sidecar = source.sidecars["pdb_fragment_ingestion"].read_text(encoding="utf-8")

    assert source.source_kind == "pdb_fragment"
    assert source.reactive_selector["atom_name"] == "C1"
    assert "n_glycosylation_profile" in sidecar


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
    assert "atoms 2\n  inside sphere" in result.construction.placement.packmol_input_text
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
    assert not np.allclose(_xyz(c1), _xyz(spec.modifier_link_atom), atol=1.0e-3)


def test_coordinate_only_accepts_generic_local_hydroxyl_glycan(tmp_path: Path) -> None:
    """Coordinate-only export should accept generic residue-resolved glycan labels."""
    glycan_path = _write_glygen_fixture(tmp_path / "external_like.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    result, _spec = _build_coordinate_only_fixture_result(tmp_path, protein_path, glycan_path)

    output = result.crosslinked_conjugate_pdb_path.read_text(encoding="utf-8")
    output_atoms = tuple(parse_pdb_atom_records(result.crosslinked_conjugate_pdb_path))
    bonds = parse_pdb_conect_pairs(result.crosslinked_conjugate_pdb_path)
    nd2 = _single_atom(output_atoms, chain_id="A", residue_number=1, atom_name="ND2")
    c1 = _single_atom(output_atoms, chain_id="C", residue_number=1, atom_name="C1")

    assert result.status == "coordinate_only"
    assert " NAG C" in output
    assert " MAN C" in output
    assert "HO1" not in output
    assert any(set(pair) == {nd2.serial, c1.serial} for pair in bonds)
    assert "CONECT" in output


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


def test_glygen_loader_rejects_high_serial_coordinate_only_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Strict GlyGen loading should reject high-serial input without CONECT records."""
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

    with pytest.raises(ValueError, match="requires complete CONECT"):
        load_pdb_fragment(glycan_path)


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

    assert spec.fragment.bonds[0] == (101, 102)
    assert all(endpoint >= 101 for bond in spec.fragment.bonds for endpoint in bond)

    output_edges = {
        frozenset(edge) for edge in parse_pdb_conect_pairs(result.crosslinked_conjugate_pdb_path)
    }
    output_atoms = tuple(parse_pdb_atom_records(result.crosslinked_conjugate_pdb_path))
    output_glycan_atoms = tuple(atom for atom in output_atoms if atom.chain_id == "C")
    source_atoms = tuple(
        atom for atom in parse_pdb_atom_records(glycan_path) if atom.residue_name.upper() != "ROH"
    )
    output_by_source = _output_glycan_atoms_by_source(source_atoms, output_glycan_atoms)
    expected_c1_o5 = frozenset((output_by_source[104].serial, output_by_source[114].serial))

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
            force_field="glycam06",
        ),
        mechanism=SimpleNamespace(
            name="n_glycosylation",
            product_residues=SimpleNamespace(site=None, moiety=None),
            bond=SimpleNamespace(order=1, site_atom=None, target_bond_length_angstrom=1.45),
        ),
    )


def _asn_selector() -> PdbAtomSelector:
    """Return the default ASN ND2 selector for tests."""
    return PdbAtomSelector(
        chain_id="A",
        residue_name="ASN",
        residue_number=1,
        atom_name="ND2",
        insertion_code="",
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


def _write_asn_fixture(path: Path, *, hydrogens: tuple[str, ...] = ("HD21",)) -> Path:
    """Write a tiny explicit-hydrogen ASN residue fixture."""
    atoms = [
        (1, "N", "ASN", "A", 1, 0.0, 0.0, 0.0, "N"),
        (2, "CA", "ASN", "A", 1, 1.4, 0.0, 0.0, "C"),
        (3, "C", "ASN", "A", 1, 2.0, 1.2, 0.0, "C"),
        (4, "CB", "ASN", "A", 1, 1.8, -1.2, 0.0, "C"),
        (5, "CG", "ASN", "A", 1, 3.0, -1.2, 0.0, "C"),
        (6, "OD1", "ASN", "A", 1, 3.6, -2.2, 0.0, "O"),
        (7, "ND2", "ASN", "A", 1, 3.5, 0.0, 0.0, "N"),
    ]
    hydrogen_coords = {
        "HD21": (8, "HD21", "ASN", "A", 1, 14.4, 0.0, 0.0, "H"),
        "HD22": (9, "HD22", "ASN", "A", 1, 13.1, 0.8, 0.0, "H"),
    }
    atoms.extend(hydrogen_coords[name] for name in hydrogens)
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
    source_atoms = parse_pdb_atom_records(G42666_CONECT_PATH)
    residue_map = {"4YB": ("NAG", 1), "0YB": ("MAN", 2), "ROH": ("ROH", 3)}
    lines = []
    for atom in source_atoms:
        residue_name, residue_number = residue_map[atom.residue_name]
        lines.append(
            _pdb_atom(
                atom.serial + serial_offset,
                atom.atom_name,
                residue_name,
                "",
                residue_number,
                atom.x + dx,
                atom.y + dy,
                atom.z + dz,
                atom.element,
            )
        )
    if include_conect:
        for left, right in parse_pdb_conect_pairs(G42666_CONECT_PATH):
            if left < right:
                lines.append(_conect(left + serial_offset, right + serial_offset))
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def _write_generic_local_hydroxyl_fixture(path: Path) -> Path:
    """Write a generic residue-resolved glycan with a residue-local reducing OH."""
    atoms = [
        (1, "A1", "QAA", "", 1, 0.0, 0.0, 0.0, "C"),
        (2, "B2", "QAA", "", 1, 0.0, 1.4, 0.0, "O"),
        (3, "C3", "QAA", "", 1, 1.5, 0.0, 0.0, "C"),
        (4, "D4", "QAA", "", 1, 2.8, 0.0, 0.0, "O"),
        (5, "E5", "Z9B", "", 2, 5.2, 0.0, 0.0, "O"),
        (6, "F6", "Z9B", "", 2, 4.1, 0.0, 0.0, "C"),
        (7, "G7", "Z9B", "", 2, 4.1, 1.4, 0.0, "O"),
        (8, "X7", "QAA", "", 1, -1.3, 0.0, 0.0, "O"),
        (9, "Y8", "QAA", "", 1, -2.2, 0.0, 0.0, "H"),
    ]
    lines = [_pdb_atom(*atom) for atom in atoms]
    lines.extend(
        [
            _conect(1, 2, 3, 8),
            _conect(3, 4),
            _conect(4, 6),
            _conect(6, 5, 7, 2),
            _conect(8, 9),
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


def _pdb_record(
    *,
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
) -> Any:
    """Return one parsed PDB atom record for endpoint-resolution tests."""
    from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord

    return PdbAtomRecord(
        serial=serial,
        atom_name=atom_name,
        residue_name=residue_name,
        chain_id=chain_id,
        residue_number=residue_number,
        x=0.0,
        y=0.0,
        z=0.0,
        element=atom_name[0],
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
