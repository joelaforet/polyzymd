"""Tests for peptide-capped product-state charge patches."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.pablo import charge_patch


def test_peptide_capped_records_apply_strict_cap_closure(monkeypatch):
    """Cap/flank residuals should close over sorted modified-residue atoms only."""
    spec = _spec()
    reference = _reference_build()
    molecule = SimpleNamespace(
        atoms=tuple(SimpleNamespace(metadata={"atom_map": index}) for index in range(1, 6)),
        partial_charges=[0.10, -0.05, 0.20, -0.10, 0.01],
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
    assert charges["CA"] == pytest.approx(0.10 - 0.01 / 3.0)
    assert charges["CB"] == pytest.approx(-0.05 - 0.01 / 3.0)
    assert charges["NZ"] == pytest.approx(0.20 - 0.01 / 3.0)
    assert charges["C001"] == pytest.approx(-0.10)


def test_peptide_capped_records_reject_large_cap_residual(monkeypatch):
    """Large cap/flank residuals should fail instead of silently normalizing."""
    molecule = SimpleNamespace(
        atoms=tuple(SimpleNamespace(metadata={"atom_map": index}) for index in range(1, 6)),
        partial_charges=[0.10, -0.05, 0.20, -0.10, 0.5],
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
        generated_fragment=fragment,
        product_residue_mappings={
            "1": {"target_chain": "C", "target_residue_name": "NHX", "target_residue_number": 1}
        },
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
