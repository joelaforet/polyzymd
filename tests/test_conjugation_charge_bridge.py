"""Tests for product-state charge bridge records and mapping."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.pablo import charge_bridge
from polyzymd.builders.conjugation.pablo.charge_records import (
    AtomPartialChargeRecord,
    ResiduePartialChargeRecord,
    validate_unique_atom_records,
)


def test_residue_records_group_atom_records():
    """Atom charge records should group by residue and provenance."""
    records = (
        AtomPartialChargeRecord(
            chain_id="A",
            residue_name="LYX",
            residue_number=10,
            atom_name="NZ",
            charge_e=-0.2,
            source="production:ff14SB",
            source_role="protein_ff14sb",
        ),
        AtomPartialChargeRecord(
            chain_id="A",
            residue_name="LYX",
            residue_number=10,
            atom_name="CE",
            charge_e=0.2,
            source="production:ff14SB",
            source_role="protein_ff14sb",
        ),
    )

    grouped = ResiduePartialChargeRecord.from_atom_records(records)

    assert len(grouped) == 1
    assert grouped[0].atom_charges == {"NZ": -0.2, "CE": 0.2}


def test_validate_unique_atom_records_rejects_duplicate_identity():
    """Duplicate atom identities should fail before template assembly."""
    record = AtomPartialChargeRecord(
        chain_id="A",
        residue_name="LYX",
        residue_number=10,
        atom_name="NZ",
        charge_e=-0.2,
        source="production:ff14SB",
        source_role="protein_ff14sb",
    )

    with pytest.raises(ValueError, match="Duplicate partial-charge record"):
        validate_unique_atom_records((record, record))


def test_build_product_state_charge_bridge_combines_sources(monkeypatch, tmp_path):
    """Bridge should assemble ff14SB, polymer-template, and local patch records."""
    target = _molecule(
        [
            _atom("A", "LYX", 10, "NZ", 0),
            _atom("C", "NHX", 1, "C047", 0),
            _atom("C", "NHX", 1, "O020", 0),
            _atom("C", "NHX", 1, "C001", 0),
        ]
    )
    library = SimpleNamespace(residue_names=("LYX", "NHX"), definitions=())
    spec = SimpleNamespace()

    monkeypatch.setattr(
        charge_bridge,
        "_protein_ff14sb_records",
        lambda **_: (
            AtomPartialChargeRecord(
                chain_id="A",
                residue_name="LYX",
                residue_number=10,
                atom_name="NZ",
                charge_e=-0.1,
                source="production:ff14SB",
                source_role="protein_ff14sb",
            ),
        ),
    )
    monkeypatch.setattr(
        charge_bridge,
        "_polymer_template_records",
        lambda _: (
            AtomPartialChargeRecord(
                chain_id="C",
                residue_name="NHX",
                residue_number=1,
                atom_name="C047",
                charge_e=0.1,
                source="production:polymer",
                source_role="polymer_template",
            ),
            AtomPartialChargeRecord(
                chain_id="C",
                residue_name="NHX",
                residue_number=1,
                atom_name="O020",
                charge_e=-0.3,
                source="production:polymer",
                source_role="polymer_template",
            ),
            AtomPartialChargeRecord(
                chain_id="C",
                residue_name="NHX",
                residue_number=1,
                atom_name="C001",
                charge_e=0.3,
                source="production:polymer",
                source_role="polymer_template",
            ),
        ),
    )
    monkeypatch.setattr(
        charge_bridge,
        "_local_nagl_patch_records",
        lambda _: (
            (
                AtomPartialChargeRecord(
                    chain_id="A",
                    residue_name="LYX",
                    residue_number=10,
                    atom_name="NZ",
                    charge_e=-0.25,
                    source="production:nagl",
                    source_role="local_nagl_patch",
                ),
            ),
            "test-nagl.pt",
        ),
    )
    monkeypatch.setattr(charge_bridge, "parse_pdb_atom_records", lambda _: ())

    result = charge_bridge.build_product_state_charge_bridge(
        product_state_pablo_library=library,
        product_topology=SimpleNamespace(molecules=(target,)),
        product_pdb=tmp_path / "product.pdb",
        source_protein_pdb=tmp_path / "source.pdb",
        specs=(spec,),
        output_dir=tmp_path,
    )

    assert result.charge_bridge_report.local_nagl_patch_atom_count == 1
    assert result.charge_bridge_report.polymer_template_atom_count == 2
    assert (tmp_path / "product_state_charge_bridge.json").is_file()
    charges = {
        (record.residue_name, atom_name): charge
        for record in result.residue_partial_charges
        for atom_name, charge in record.atom_charges.items()
    }
    assert charges[("LYX", "NZ")] == pytest.approx(-0.25)


def _molecule(atoms: list[SimpleNamespace]) -> SimpleNamespace:
    """Build an OpenFF-like molecule double."""
    return SimpleNamespace(atoms=tuple(atoms), n_atoms=len(atoms), partial_charges=None)


def _atom(
    chain_id: str,
    residue_name: str,
    residue_number: int,
    atom_name: str,
    formal_charge: int,
) -> SimpleNamespace:
    """Build an OpenFF-like atom double."""
    return SimpleNamespace(
        name=atom_name,
        formal_charge=formal_charge,
        metadata={
            "chain_id": chain_id,
            "residue_name": residue_name,
            "residue_number": residue_number,
            "atom_name": atom_name,
        },
    )
