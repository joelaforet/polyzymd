"""Tests for peptide-capped product-state charge patches."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.pablo import charge_patch


def test_peptide_capped_records_apply_strict_cap_closure(monkeypatch):
    """Local residuals should close over sorted modified-residue atoms only."""
    spec = _spec()
    reference = _reference_build()
    molecule = SimpleNamespace(
        atoms=tuple(SimpleNamespace(metadata={"atom_map": index}) for index in range(1, 6)),
        partial_charges=[0.10, -0.05, 0.041, -0.10, 0.01],
    )

    monkeypatch.setattr(
        charge_patch, "_build_reference_product", lambda *_args, **_kwargs: reference
    )
    monkeypatch.setattr(charge_patch, "_charge_with_nagl", lambda *_args, **_kwargs: molecule)

    records, model_name = charge_patch.build_local_product_charge_patch_records(
        spec,
        product_atoms=_product_atoms(),
    )

    charges = {record.atom_name: record.charge_e for record in records}
    assert model_name == charge_patch.DEFAULT_PATCH_NAGL_MODEL
    assert charges["CA"] == pytest.approx(0.10 + 0.009 / 3.0)
    assert charges["CB"] == pytest.approx(-0.05 + 0.009 / 3.0)
    assert charges["NZ"] == pytest.approx(0.041 + 0.009 / 3.0)
    assert charges["C001"] == pytest.approx(-0.10)
    assert sum(charges.values()) == pytest.approx(0.0, abs=1.0e-12)


def test_peptide_capped_records_keep_unmapped_cap_charges_distinct(monkeypatch):
    """Unmapped cap atoms should not collide with mapped product atom numbers."""
    spec = _spec()
    reference = _reference_build()
    molecule = SimpleNamespace(
        atoms=(
            SimpleNamespace(metadata={}),
            SimpleNamespace(metadata={"atom_map": 1}),
            SimpleNamespace(metadata={"atom_map": 2}),
            SimpleNamespace(metadata={"atom_map": 3}),
            SimpleNamespace(metadata={"atom_map": 4}),
        ),
        partial_charges=[0.01, 0.10, -0.05, 0.041, -0.10],
    )

    monkeypatch.setattr(
        charge_patch, "_build_reference_product", lambda *_args, **_kwargs: reference
    )
    monkeypatch.setattr(charge_patch, "_charge_with_nagl", lambda *_args, **_kwargs: molecule)

    records, _model_name = charge_patch.build_local_product_charge_patch_records(
        spec,
        product_atoms=_product_atoms(),
    )

    charges = {record.atom_name: record.charge_e for record in records}
    assert charges["CA"] == pytest.approx(0.10 + 0.009 / 3.0)
    assert charges["CB"] == pytest.approx(-0.05 + 0.009 / 3.0)
    assert charges["NZ"] == pytest.approx(0.041 + 0.009 / 3.0)
    assert charges["C001"] == pytest.approx(-0.10)


def test_partial_charges_assign_negative_keys_to_unmapped_reference_atoms():
    """Unmapped cap and flank atoms should retain charge without stealing product maps."""
    molecule = SimpleNamespace(
        atoms=(
            SimpleNamespace(metadata={}),
            SimpleNamespace(metadata={"atom_map": 1}),
            SimpleNamespace(metadata={}),
            SimpleNamespace(metadata={"atom_map": 2}),
        ),
        partial_charges=[0.25, -0.10, 0.05, -0.20],
    )

    charges = charge_patch._partial_charges(molecule)

    assert charges == {-1: 0.25, 1: -0.10, -3: 0.05, 2: -0.20}


def test_peptide_capped_records_reject_large_cap_residual(monkeypatch):
    """Large local residuals should fail instead of silently normalizing."""
    molecule = SimpleNamespace(
        atoms=tuple(SimpleNamespace(metadata={"atom_map": index}) for index in range(1, 6)),
        partial_charges=[0.10, -0.05, 0.20, -0.10, 0.0],
    )
    monkeypatch.setattr(
        charge_patch, "_build_reference_product", lambda *_args, **_kwargs: _reference_build()
    )
    monkeypatch.setattr(charge_patch, "_charge_with_nagl", lambda *_args, **_kwargs: molecule)

    with pytest.raises(charge_patch.LocalChargePatchError, match="closure failed"):
        charge_patch.build_local_product_charge_patch_records(
            _spec(),
            product_atoms=_product_atoms(),
        )


def test_local_formal_projection_accepts_g42666_asx_residual():
    """The observed G42666 ASX residual should distribute below the per-atom bound."""
    reference = _reference_build_with_closure_domain(13)
    charges = dict.fromkeys(range(1, 14), -0.18544238 / 13.0)
    charges[14] = 0.15030887

    records, closure = charge_patch._records_with_cap_closure(reference, charges)

    assert closure["target_formal_charge_e"] == pytest.approx(0.0)
    assert closure["target_scope"] == "all_emitted_mapped_local_product_atoms"
    assert closure["correction_domain"] == "modified_protein_residue_closure_atoms"
    assert closure["residual_to_target_e"] == pytest.approx(0.03513351)
    assert closure["closure_atom_count"] == 13
    assert closure["per_atom_closure_e"] == pytest.approx(0.03513351 / 13.0)
    assert closure["max_per_atom_closure_e"] < 0.005
    assert closure["final_projected_total_e"] == pytest.approx(0.0, abs=1.0e-12)
    assert sum(record.charge_e for record in records) == pytest.approx(0.0, abs=1.0e-12)


def test_local_formal_projection_rejects_excessive_per_atom_correction():
    """Residuals above 0.005 e per closure atom should remain hard failures."""
    reference = _reference_build_with_closure_domain(13)
    charges = dict.fromkeys(range(1, 14), -0.216)
    charges[14] = 2.742

    with pytest.raises(charge_patch.LocalChargePatchError, match="exceeds 0.005 e"):
        charge_patch._records_with_cap_closure(reference, charges)


def test_peptide_capped_records_reject_terminal_site():
    """Terminal side-chain modifications are outside the first-release scope."""
    spec = _spec(residue_number=1)

    with pytest.raises(charge_patch.LocalChargePatchError, match="Terminal"):
        charge_patch.build_local_product_charge_patch_records(spec, product_atoms=_product_atoms(1))


def test_fragment_bonds_prefer_serial_when_index_collides():
    """One-based serial bonds should not be remapped to zero-based atom-index neighbors."""
    atom0 = _atom(0, 1, "C1", "MOD", 1, element="C", chain_id="C")
    atom1 = _atom(1, 2, "H1", "MOD", 1, element="H", chain_id="C")
    fragment = SimpleNamespace(atoms=(atom0, atom1), bonds=((1, 2),), bond_orders=())

    bonds = charge_patch._fragment_bonds(fragment)

    assert bonds == ((atom0, atom1, 1.0),)


def test_g42666_modifier_link_prefers_source_serial_and_index_with_duplicate_c1_names():
    """Glycan-like C1 duplicates should resolve through source serial and index first."""
    spec = _g42666_like_spec()

    resolved = charge_patch._resolve_retained_modifier_link_atom(spec, spec.fragment)
    retained = charge_patch._retained_fragment_atoms(spec, spec.fragment)
    mapped = charge_patch._retained_modifier_atoms(
        spec,
        spec.fragment,
        product_atoms=_g42666_product_atoms(),
    )

    assert resolved.serial == 4
    assert resolved.atom_index == 3
    assert {atom.atom_name for atom in retained} == {"C1", "C2"}
    assert {atom.atom_name for atom, _mapped in mapped} == {"C1", "C2"}
    assert all(atom.atom_name not in {"O1", "HO1"} for atom in retained)
    assert {mapped_atom.residue_number for _atom, mapped_atom in mapped} == {42666, 42667}


def test_g42666_modifier_link_rejects_ambiguous_name_only_reactive_atom():
    """Name fallback should fail when retained source atom names are not unique."""
    spec = _g42666_like_spec(reactive_atom_serial=None, reactive_atom_index=None)

    with pytest.raises(charge_patch.LocalChargePatchError, match="ambiguous"):
        charge_patch._resolve_retained_modifier_link_atom(spec, spec.fragment)


def test_g42666_modifier_link_rejects_explicit_mapping_mismatch():
    """Explicit source-to-product residue mappings should not use global name fallback."""
    spec = _g42666_like_spec(mapping_residue_number=99999)

    with pytest.raises(charge_patch.LocalChargePatchError, match="exact product atom identity"):
        charge_patch._retained_modifier_atoms(
            spec,
            spec.fragment,
            product_atoms=_g42666_product_atoms(),
        )


def test_g42666_modifier_link_respects_explicit_non_c_target_chain():
    """Explicit target_chain mappings should resolve modifier atoms outside chain C."""
    spec = _g42666_like_spec(mapping_chain_id="D")
    mapped = charge_patch._retained_modifier_atoms(
        spec,
        spec.fragment,
        product_atoms=_g42666_product_atoms(chain_id="D"),
    )

    assert {mapped_atom.chain_id for _atom, mapped_atom in mapped} == {"D"}
    assert {mapped_atom.residue_number for _atom, mapped_atom in mapped} == {42666, 42667}


def test_modifier_lookup_rejects_duplicate_exact_product_identity():
    """Duplicate canonical product identities should fail rather than overwrite."""
    atoms = (
        _atom(100, 100, "C1", "NAG", 42666, element="C", chain_id="C"),
        _atom(101, 101, "C1", "NAG", 42666, element="C", chain_id="C"),
    )

    with pytest.raises(charge_patch.LocalChargePatchError, match="Duplicate product atom identity"):
        charge_patch._modifier_product_atom_lookup(atoms)


def test_unmapped_modifier_name_fallback_does_not_match_protein_chain_a():
    """Unmapped modifier fallback should not steal unique atom names from protein chain A."""
    spec = _spec()
    spec.product_residue_mappings = {}
    protein_only_name_match = _product_atoms()[:-1]

    mapped = charge_patch._retained_modifier_atoms(
        spec,
        spec.fragment,
        product_atoms=protein_only_name_match
        + (_atom(30, 30, "C001", "ALA", 11, element="C", chain_id="A"),),
    )

    assert mapped[0][1].chain_id == "C"
    assert mapped[0][1].residue_name == "NHX"


def test_g42666_modifier_link_rejects_reactive_serial_index_mismatch():
    """Source reactive serial and index should identify the same retained atom."""
    spec = _g42666_like_spec(reactive_atom_index=4)

    with pytest.raises(charge_patch.LocalChargePatchError, match="different source atoms"):
        charge_patch._resolve_retained_modifier_link_atom(spec, spec.fragment)


def _spec(*, residue_number: int = 10) -> SimpleNamespace:
    """Build a generic supported Lys linkage spec."""
    protein = _atom(10, 10, "NZ", "LYS", residue_number, element="N", chain_id="A")
    modifier = _atom(0, 1, "C001", "MOD", 1, element="C", chain_id="C")
    leaving = _atom(1, 2, "O001", "MOD", 1, element="O", chain_id="C")
    fragment = SimpleNamespace(
        atoms=(modifier, leaving),
        bonds=((1, 2),),
        bond_orders=((1, 2, 1.0),),
        leaving_atom_names=("O001",),
    )
    plan = SimpleNamespace(
        protein_link_atom=protein,
        modifier_link_atom=modifier,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(leaving,),
        protein_product_residue_name="LYX",
        modifier_product_residue_name="NHX",
        pablo_crosslink_requirement=SimpleNamespace(bond_order=1),
        contract=SimpleNamespace(mechanism_name="nhs_lys_amide"),
    )
    return SimpleNamespace(
        name="lys_patch",
        resolved_plan=plan,
        fragment=fragment,
        product_residue_mappings={
            "1": {"target_chain": "C", "target_residue_name": "NHX", "target_residue_number": 1}
        },
    )


def _g42666_like_spec(
    *,
    reactive_atom_serial: int | None = 4,
    reactive_atom_index: int | None = 3,
    mapping_residue_number: int = 42666,
    mapping_chain_id: str = "C",
) -> SimpleNamespace:
    """Build a glycan-like Asn-linked spec with duplicate retained C1 source names."""
    protein = _atom(10, 10, "ND2", "ASN", 42, element="N", chain_id="A")
    selected_c1 = _atom(3, 4, "C1", "NAG", 1, element="C", chain_id="C")
    leaving_o1 = _atom(1, 2, "O1", "NAG", 1, element="O", chain_id="C")
    leaving_ho1 = _atom(2, 3, "HO1", "NAG", 1, element="H", chain_id="C")
    second_c1 = _atom(4, 5, "C1", "NAG", 2, element="C", chain_id="C")
    second_c2 = _atom(5, 6, "C2", "NAG", 2, element="C", chain_id="C")
    fragment = SimpleNamespace(
        atoms=(leaving_o1, leaving_ho1, selected_c1, second_c1, second_c2),
        bonds=((2, 3), (4, 5), (5, 6)),
        bond_orders=((4, 5, 1.0), (5, 6, 1.0)),
        reactive_atom_serial=reactive_atom_serial,
        reactive_atom_index=reactive_atom_index,
        reactive_atom_name="C1",
        leaving_atom_serials=(2, 3),
        leaving_atom_indices=(1, 2),
        leaving_atom_names=("O1", "HO1"),
    )
    plan = SimpleNamespace(
        protein_link_atom=protein,
        modifier_link_atom=selected_c1,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(leaving_o1, leaving_ho1),
        protein_product_residue_name="ASX",
        modifier_product_residue_name="NAG",
        pablo_crosslink_requirement=SimpleNamespace(bond_order=1),
        contract=SimpleNamespace(mechanism_name="n_glycosylation"),
    )
    return SimpleNamespace(
        name="g42666_patch",
        resolved_plan=plan,
        fragment=fragment,
        product_residue_mappings={
            "1": {
                "target_chain": mapping_chain_id,
                "target_residue_name": "NAG",
                "target_residue_number": mapping_residue_number,
            },
            "2": {
                "target_chain": mapping_chain_id,
                "target_residue_name": "NAG",
                "target_residue_number": 42667,
            },
        },
    )


def _g42666_product_atoms(*, chain_id: str = "C") -> tuple[SimpleNamespace, ...]:
    """Build final product atoms for G42666-like charge patch identity tests."""
    return (
        _atom(100, 100, "C1", "NAG", 42666, element="C", chain_id=chain_id),
        _atom(101, 101, "C1", "NAG", 42667, element="C", chain_id=chain_id),
        _atom(102, 102, "C2", "NAG", 42667, element="C", chain_id=chain_id),
    )


def _product_atoms(residue_number: int = 10) -> tuple[SimpleNamespace, ...]:
    """Build final product atoms for the supported test linkage."""
    return (
        _atom(1, 1, "N", "LYX", residue_number, element="N", chain_id="A"),
        _atom(2, 2, "CA", "LYX", residue_number, element="C", chain_id="A"),
        _atom(3, 3, "C", "LYX", residue_number, element="C", chain_id="A"),
        _atom(4, 4, "O", "LYX", residue_number, element="O", chain_id="A"),
        _atom(5, 5, "CB", "LYX", residue_number, element="C", chain_id="A"),
        _atom(6, 6, "CG", "LYX", residue_number, element="C", chain_id="A"),
        _atom(7, 7, "CD", "LYX", residue_number, element="C", chain_id="A"),
        _atom(8, 8, "CE", "LYX", residue_number, element="C", chain_id="A"),
        _atom(10, 10, "NZ", "LYX", residue_number, element="N", chain_id="A"),
        _atom(20, 1, "C001", "NHX", 1, element="C", chain_id="C"),
    )


def _reference_build() -> charge_patch._ReferenceBuild:
    """Build a minimal mapped reference for closure tests."""
    mapped = {
        1: charge_patch._MappedAtom(
            ("protein", 2, "CA"), "CA", "C", "A", "LYX", 10, "", "modified_protein_product"
        ),
        2: charge_patch._MappedAtom(
            ("protein", 5, "CB"), "CB", "C", "A", "LYX", 10, "", "modified_protein_product"
        ),
        3: charge_patch._MappedAtom(
            ("protein", 10, "NZ"), "NZ", "N", "A", "LYX", 10, "", "modified_protein_product"
        ),
        4: charge_patch._MappedAtom(
            ("modifier", 0, "C001"), "C001", "C", "C", "NHX", 1, "", "moiety_product"
        ),
    }
    return charge_patch._ReferenceBuild(
        molecule=object(),
        mapped_atoms=mapped,
        closure_map_numbers=(1, 2, 3),
        cap_map_numbers=(5,),
        product_atom_count=4,
        reference_atom_count=5,
    )


def _reference_build_with_closure_domain(closure_count: int) -> charge_patch._ReferenceBuild:
    """Build a mapped neutral reference with a configurable closure domain."""
    mapped = {
        index: charge_patch._MappedAtom(
            ("protein", index, f"C{index}"),
            f"C{index}",
            "C",
            "A",
            "ASX",
            10,
            "",
            "modified_protein_product",
        )
        for index in range(1, closure_count + 1)
    }
    mapped[closure_count + 1] = charge_patch._MappedAtom(
        ("modifier", closure_count + 1, "C001"),
        "C001",
        "C",
        "C",
        "NAG",
        1,
        "",
        "moiety_product",
    )
    return charge_patch._ReferenceBuild(
        molecule=object(),
        mapped_atoms=mapped,
        closure_map_numbers=tuple(range(1, closure_count + 1)),
        cap_map_numbers=(),
        product_atom_count=len(mapped),
        reference_atom_count=len(mapped),
    )


def _atom(
    atom_index: int,
    serial: int,
    atom_name: str,
    residue_name: str,
    residue_number: int,
    *,
    element: str,
    chain_id: str,
) -> SimpleNamespace:
    """Build a PDB-like atom double."""
    return SimpleNamespace(
        atom_index=atom_index,
        serial=serial,
        atom_name=atom_name,
        element=element,
        chain_id=chain_id,
        residue_name=residue_name,
        residue_number=residue_number,
        insertion_code="",
    )
