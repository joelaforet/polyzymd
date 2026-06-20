"""Tests for generic local product-state charge patches."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.pablo import charge_patch


def test_patch_graph_uses_metadata_without_nhs_atom_names():
    """Patch graph roots and leaving atoms should come from metadata only."""
    spec = _spec(
        protein_atom_name="QX",
        modifier_atom_name="M1",
        leaving_atom_name="LG",
    )

    graph = charge_patch._build_product_graph(spec, product_atoms=())
    selected = charge_patch._select_patch_keys(graph, radius=1)
    selected_names = {graph.atoms[key].atom_name for key in selected}

    assert "QX" in selected_names
    assert "M1" in selected_names
    assert "LG" not in selected_names
    assert not {"NZ", "C047", "O020"}.intersection(selected_names)


def test_local_patch_records_map_charges_to_real_atoms(monkeypatch):
    """Local patch records should map charged patch atoms back to product identities."""
    spec = _spec()
    molecule = SimpleNamespace(
        atoms=(
            SimpleNamespace(metadata={"atom_map": 1}),
            SimpleNamespace(metadata={"atom_map": 2}),
            SimpleNamespace(metadata={"atom_map": 3}),
        ),
        partial_charges=[0.2, -0.1, -0.3],
    )

    def fake_charge(molecule, *, model_name):
        return molecule

    def fake_build_patch(graph, selected):
        map_numbers = {}
        for key in selected:
            name = graph.atoms[key].atom_name
            map_numbers[key] = {"QX": 0, "M1": 1}.get(name, 2)
        return molecule, map_numbers

    monkeypatch.setattr(charge_patch, "_build_off_patch_molecule", fake_build_patch)
    monkeypatch.setattr(charge_patch, "_charge_with_nagl", fake_charge)

    records, model_name = charge_patch.build_local_product_charge_patch_records(spec)

    charges = {(record.residue_name, record.atom_name): record.charge_e for record in records}
    assert model_name == charge_patch.DEFAULT_PATCH_NAGL_MODEL
    assert charges[("PRX", "QX")] == pytest.approx(0.2)
    assert charges[("MRX", "M1")] == pytest.approx(-0.1)


def test_patch_fails_when_fragment_metadata_missing():
    """Insufficient mechanism metadata should fail clearly instead of hardcoding chemistry."""
    spec = SimpleNamespace(resolved_plan=SimpleNamespace(), generated_fragment=None)

    with pytest.raises(charge_patch.LocalChargePatchError, match="generated fragment metadata"):
        charge_patch.build_local_product_charge_patch_records(spec)


def test_fragment_bonds_prefer_serial_when_index_collides():
    """One-based serial bonds should not be remapped to zero-based atom-index neighbors."""
    atom0 = _atom(
        atom_index=0,
        atom_name="C1",
        element="C",
        chain_id="C",
        residue_name="MOD",
        residue_number=1,
    )
    atom1 = _atom(
        atom_index=1,
        atom_name="H1",
        element="H",
        chain_id="C",
        residue_name="MOD",
        residue_number=1,
    )
    fragment = SimpleNamespace(atoms=(atom0, atom1), bonds=((1, 2),), bond_orders=())

    bonds = charge_patch._fragment_bonds(fragment)

    assert bonds == ((atom0, atom1, 1.0),)


def _spec(
    *,
    protein_atom_name: str = "QX",
    modifier_atom_name: str = "M1",
    leaving_atom_name: str = "LG",
) -> SimpleNamespace:
    """Build a generic linkage spec without NHS-Lys names."""
    protein = _atom(
        atom_index=0,
        atom_name=protein_atom_name,
        element="N",
        chain_id="A",
        residue_name="SRC",
        residue_number=10,
    )
    modifier = _atom(
        atom_index=1,
        atom_name=modifier_atom_name,
        element="C",
        chain_id="C",
        residue_name="MOD",
        residue_number=1,
    )
    oxygen = _atom(
        atom_index=2,
        atom_name="MO",
        element="O",
        chain_id="C",
        residue_name="MOD",
        residue_number=1,
    )
    leaving = _atom(
        atom_index=3,
        atom_name=leaving_atom_name,
        element="O",
        chain_id="C",
        residue_name="MOD",
        residue_number=1,
    )
    fragment = SimpleNamespace(
        atoms=(modifier, oxygen, leaving),
        bonds=((1, 2), (1, 3)),
        bond_orders=((1, 2, 2.0), (1, 3, 1.0)),
        leaving_atom_names=(leaving_atom_name,),
    )
    plan = SimpleNamespace(
        protein_link_atom=protein,
        modifier_link_atom=modifier,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(leaving,),
        protein_product_residue_name="PRX",
        modifier_product_residue_name="MRX",
        pablo_crosslink_requirement=SimpleNamespace(bond_order=1),
    )
    return SimpleNamespace(
        resolved_plan=plan,
        generated_fragment=fragment,
        product_residue_mappings={
            "1": {"target_chain": "C", "target_residue_name": "MRX", "target_residue_number": 1}
        },
    )


def _atom(
    *,
    atom_index: int,
    atom_name: str,
    element: str,
    chain_id: str,
    residue_name: str,
    residue_number: int,
) -> SimpleNamespace:
    """Build a PDB-like atom double."""
    return SimpleNamespace(
        atom_index=atom_index,
        serial=atom_index + 1,
        atom_name=atom_name,
        element=element,
        chain_id=chain_id,
        residue_name=residue_name,
        residue_number=residue_number,
        insertion_code="",
    )
