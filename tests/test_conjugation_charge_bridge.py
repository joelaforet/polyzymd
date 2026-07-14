"""Tests for product-state charge bridge records and mapping."""

from __future__ import annotations

import inspect
import json
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.pablo import charge_bridge
from polyzymd.builders.conjugation.pablo.charge_records import (
    AtomPartialChargeRecord,
    ResiduePartialChargeRecord,
    validate_unique_atom_records,
)
from polyzymd.builders.conjugation.pablo.product_state import target_identities_from_molecule


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


def test_residue_records_preserve_ordered_atom_records():
    """Ordered atom records should emit one residue record per atom."""
    records = (
        AtomPartialChargeRecord(
            chain_id="A",
            residue_name="LYX",
            residue_number=10,
            atom_name="NZ",
            charge_e=-0.4,
            source="production:ff14SB",
            source_role="protein_ff14sb",
        ),
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
            chain_id="A",
            residue_name="LYX",
            residue_number=10,
            atom_name="CE",
            charge_e=0.3,
            source="production:ff14SB",
            source_role="protein_ff14sb",
        ),
    )

    ordered = ResiduePartialChargeRecord.from_ordered_atom_records(records)

    assert len(ordered) == len(records)
    assert [tuple(record.atom_charges) for record in ordered] == [("NZ",), ("C047",), ("CE",)]
    assert [record.source_role for record in ordered] == [
        "protein_ff14sb",
        "polymer_template",
        "protein_ff14sb",
    ]
    assert [next(iter(record.atom_charges.values())) for record in ordered] == pytest.approx(
        (-0.4, 0.1, 0.3)
    )


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
        lambda _, **__: (
            AtomPartialChargeRecord(
                chain_id="C",
                residue_name="NHX",
                residue_number=1,
                atom_name="C047",
                charge_e=0.25,
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
        lambda _, **__: (
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
    assert result.charge_bridge_report.polymer_template_atom_count == 3
    assert result.charge_bridge_report.total_partial_charge_before_correction_e == pytest.approx(
        0.0
    )
    assert (
        result.charge_bridge_report.source
        == "production:product-state-peptide-capped-nagl-charge-bridge"
    )
    assert result.charge_bridge_report.order_preserving_atom_records is True
    assert (tmp_path / "product_state_charge_bridge.json").is_file()
    assert len(result.residue_partial_charges) == target.n_atoms
    assert [tuple(record.atom_charges) for record in result.residue_partial_charges] == [
        ("NZ",),
        ("C047",),
        ("O020",),
        ("C001",),
    ]
    charges = {
        (record.residue_name, atom_name): charge
        for record in result.residue_partial_charges
        for atom_name, charge in record.atom_charges.items()
    }
    assert charges[("LYX", "NZ")] == pytest.approx(-0.25)


def test_build_product_state_charge_bridge_has_no_public_total_charge_tolerance():
    """Strict bridge reconciliation should not expose a caller-relaxed tolerance."""
    signature = inspect.signature(charge_bridge.build_product_state_charge_bridge)

    assert "total_charge_tolerance" not in signature.parameters


def test_build_product_state_charge_bridge_rejects_caller_total_charge_tolerance(tmp_path):
    """Callers should not be able to bypass fixed bridge reconciliation limits."""
    with pytest.raises(TypeError, match="total_charge_tolerance"):
        charge_bridge.build_product_state_charge_bridge(
            product_state_pablo_library=SimpleNamespace(residue_names=("LYX",), definitions=()),
            product_topology=SimpleNamespace(molecules=()),
            product_pdb=tmp_path / "product.pdb",
            source_protein_pdb=tmp_path / "source.pdb",
            specs=(SimpleNamespace(),),
            output_dir=tmp_path,
            total_charge_tolerance=1.0,
        )


def test_target_identities_validate_product_pdb_metadata():
    """Product-state target identities should accept same-order PDB metadata."""
    product_atoms = (_product_atom(1, 0, "A", "LYX", 10, "NZ"),)
    target = _molecule([_atom_from_product_atom(product_atoms[0])])

    identities = target_identities_from_molecule(target, product_atoms=product_atoms)

    assert identities == (("A", "LYX", 10, "", "NZ"),)


def test_target_identities_reject_complete_stale_product_metadata():
    """Complete but stale atom metadata should not drive charge transfer."""
    product_atoms = (_product_atom(1, 0, "A", "LYX", 10, "NZ"),)
    stale_atom = _atom_from_product_atom(_product_atom(99, 0, "C", "NHX", 1, "C001"))
    target = _molecule([stale_atom])

    with pytest.raises(ValueError, match="stale or reordered metadata"):
        target_identities_from_molecule(target, product_atoms=product_atoms)


def test_target_identities_reject_product_atom_count_mismatch():
    """Topology and product PDB atom counts must agree before charge transfer."""
    product_atoms = (
        _product_atom(1, 0, "A", "LYX", 10, "NZ"),
        _product_atom(2, 1, "A", "LYX", 10, "CE"),
    )
    target = _molecule([_atom_from_product_atom(product_atoms[0])])

    with pytest.raises(ValueError, match="atom-count mismatch"):
        target_identities_from_molecule(target, product_atoms=product_atoms)


def test_target_identities_reject_product_atom_order_mismatch():
    """Product charge transfer should fail when metadata order is swapped."""
    product_atoms = (
        _product_atom(1, 0, "A", "LYX", 10, "NZ"),
        _product_atom(2, 1, "A", "LYX", 10, "CE"),
    )
    target = _molecule(
        [_atom_from_product_atom(product_atoms[1]), _atom_from_product_atom(product_atoms[0])]
    )

    with pytest.raises(ValueError, match="stale or reordered metadata"):
        target_identities_from_molecule(target, product_atoms=product_atoms)


def test_target_identities_reject_serial_renumbering_for_charge_transfer():
    """Serial renumbering should not be accepted during product charge transfer."""
    product_atom = _product_atom(1, 0, "A", "LYX", 10, "NZ")
    target_atom = _atom_from_product_atom(product_atom)
    target_atom.metadata["product_atom_serial"] = 101
    target = _molecule([target_atom])

    with pytest.raises(ValueError, match="product_atom_serial"):
        target_identities_from_molecule(target, product_atoms=(product_atom,))


def test_bridge_rejects_reconciliation_without_local_patch_atoms(monkeypatch, tmp_path):
    """Charge residuals should not be hidden on unrelated polymer atoms."""
    target = _molecule([_atom("C", "NHX", 1, "C001", 0)])
    library = SimpleNamespace(residue_names=("NHX",), definitions=())
    record = AtomPartialChargeRecord(
        chain_id="C",
        residue_name="NHX",
        residue_number=1,
        atom_name="C001",
        charge_e=-1.2,
        source="production:polymer",
        source_role="polymer_template",
    )
    monkeypatch.setattr(charge_bridge, "_protein_ff14sb_records", lambda **_: ())
    monkeypatch.setattr(charge_bridge, "_polymer_template_records", lambda _, **__: (record,))
    monkeypatch.setattr(charge_bridge, "_local_nagl_patch_records", lambda _, **__: ((), None))
    monkeypatch.setattr(charge_bridge, "parse_pdb_atom_records", lambda _: ())

    with pytest.raises(ValueError, match="no mapped modified-protein product residue atoms"):
        charge_bridge.build_product_state_charge_bridge(
            product_state_pablo_library=library,
            product_topology=SimpleNamespace(molecules=(target,)),
            product_pdb=tmp_path / "product.pdb",
            source_protein_pdb=tmp_path / "source.pdb",
            specs=(SimpleNamespace(),),
            output_dir=tmp_path,
        )

    assert not (tmp_path / "product_state_charge_bridge_local_reconciliation.json").exists()


def test_bridge_reconciles_small_residual_over_local_patch_atoms(monkeypatch, tmp_path):
    """Small residuals should be auditable and local to NAGL patch atoms."""
    target = _molecule([_atom("A", "LYX", 10, "NZ", 0), _atom("C", "NHX", 1, "C001", 0)])
    library = SimpleNamespace(residue_names=("NHX",), definitions=())
    polymer_record = AtomPartialChargeRecord(
        chain_id="C",
        residue_name="NHX",
        residue_number=1,
        atom_name="C001",
        charge_e=-0.10,
        source="production:polymer",
        source_role="polymer_template",
    )
    patch_record = AtomPartialChargeRecord(
        chain_id="A",
        residue_name="LYX",
        residue_number=10,
        atom_name="NZ",
        charge_e=0.096,
        source="production:nagl-patch",
        source_role="local_nagl_patch",
    )
    monkeypatch.setattr(charge_bridge, "_protein_ff14sb_records", lambda **_: ())
    monkeypatch.setattr(
        charge_bridge, "_polymer_template_records", lambda _, **__: (polymer_record,)
    )
    monkeypatch.setattr(
        charge_bridge,
        "_local_nagl_patch_records",
        lambda _, **__: ((patch_record,), "nagl-test"),
    )
    monkeypatch.setattr(charge_bridge, "parse_pdb_atom_records", lambda _: ())

    result = charge_bridge.build_product_state_charge_bridge(
        product_state_pablo_library=library,
        product_topology=SimpleNamespace(molecules=(target,)),
        product_pdb=tmp_path / "product.pdb",
        source_protein_pdb=tmp_path / "source.pdb",
        specs=(SimpleNamespace(),),
        output_dir=tmp_path,
    )

    report = result.charge_bridge_report
    assert report.normalization_correction_e == pytest.approx(0.004)
    assert report.max_per_atom_correction_e == pytest.approx(0.004)
    assert report.correction_atom_identities == ("chain A residue LYX 10 atom NZ",)
    diagnostic = json.loads(
        (tmp_path / "product_state_charge_bridge_local_reconciliation.json").read_text(
            encoding="utf-8"
        )
    )
    reconciliation = diagnostic["local_reconciliation"]
    assert reconciliation["success"] is True
    assert reconciliation["corrected_atom_count"] == 1
    assert reconciliation["per_atom_correction_e"] == pytest.approx(0.004)


def test_bridge_refuses_raw_sdf_as_production_charge_source(tmp_path):
    """Raw bond SDF sidecars should not be accepted as production charges."""
    raw_sdf = tmp_path / "polymer.sdf"
    raw_sdf.write_text("", encoding="utf-8")
    spec = SimpleNamespace(source_sidecars={"sdf": raw_sdf, "bond_sdf": raw_sdf})

    with pytest.raises(ValueError, match=r"requires source_sidecars\['charged_sdf'\]"):
        charge_bridge._source_sdf_path(spec)


def test_bridge_skips_raw_sdf_for_smiles_moiety_patch(tmp_path):
    """SMILES moieties should be charged by product patch, not raw SDF transfer."""
    raw_sdf = tmp_path / "glycan.sdf"
    raw_sdf.write_text("", encoding="utf-8")
    spec = SimpleNamespace(
        fragment=SimpleNamespace(source_kind="moiety"),
        source_sidecars={"sdf": raw_sdf, "bond_sdf": raw_sdf},
    )

    assert charge_bridge._source_sdf_path(spec) is None


def test_bridge_prefers_charged_sdf_source(tmp_path):
    """Charged SDF sidecars should be selected for production charge transfer."""
    raw_sdf = tmp_path / "polymer.sdf"
    charged_sdf = tmp_path / "polymer_charged.sdf"
    spec = SimpleNamespace(
        source_sidecars={"sdf": raw_sdf, "bond_sdf": raw_sdf, "charged_sdf": charged_sdf}
    )

    assert charge_bridge._source_sdf_path(spec) == charged_sdf


def test_polymer_records_refine_duplicate_atom_names_by_residue_mapping(monkeypatch, tmp_path):
    """Polymer charge transfer should resolve duplicate atom names by mapped product residue."""
    charged_sdf = tmp_path / "polymer_charged.sdf"
    charged_sdf.write_text("", encoding="utf-8")
    fragment = SimpleNamespace(
        atoms=(
            SimpleNamespace(
                atom_index=0,
                atom_name="C060",
                element="C",
                residue_name="NH2",
                residue_number=5,
                insertion_code="",
            ),
        ),
        leaving_atom_names=(),
    )
    spec = SimpleNamespace(
        source_sidecars={"charged_sdf": charged_sdf},
        generated_fragment=fragment,
        product_residue_mappings={
            "5": {"target_chain": "C", "target_residue_number": 6},
        },
    )
    product_atoms = (
        SimpleNamespace(chain_id="C", residue_name="NHX", residue_number=6, atom_name="C060"),
        SimpleNamespace(chain_id="C", residue_name="PE2", residue_number=15, atom_name="C060"),
    )
    monkeypatch.setattr(
        charge_bridge,
        "_charged_sdf_atom_charges",
        lambda *_args, **_kwargs: (0.125,),
    )

    records = charge_bridge._polymer_template_records((spec,), product_atoms=product_atoms)

    assert len(records) == 1
    assert records[0].chain_id == "C"
    assert records[0].residue_name == "NHX"
    assert records[0].residue_number == 6
    assert records[0].atom_name == "C060"
    assert records[0].charge_e == pytest.approx(0.125)


def test_polymer_records_assign_charged_sdf_charges_by_atom_index(monkeypatch, tmp_path):
    """Polymer template charges should follow atom_index order, not tuple order."""
    charged_sdf = tmp_path / "polymer_charged.sdf"
    charged_sdf.write_text("", encoding="utf-8")
    fragment = SimpleNamespace(
        atoms=(
            SimpleNamespace(
                atom_index=1,
                atom_name="O001",
                element="O",
                residue_name="NHX",
                residue_number=1,
                insertion_code="",
            ),
            SimpleNamespace(
                atom_index=0,
                atom_name="C001",
                element="C",
                residue_name="NHX",
                residue_number=1,
                insertion_code="",
            ),
        ),
        leaving_atom_names=(),
    )
    spec = SimpleNamespace(
        source_sidecars={"charged_sdf": charged_sdf},
        generated_fragment=fragment,
        product_residue_mappings={},
    )
    monkeypatch.setattr(
        charge_bridge,
        "_charged_sdf_atom_charges",
        lambda *_args, **_kwargs: (0.25, -0.5),
    )

    records = charge_bridge._polymer_template_records((spec,))

    charges_by_atom_name = {record.atom_name: record.charge_e for record in records}
    assert charges_by_atom_name == pytest.approx({"C001": 0.25, "O001": -0.5})


def test_bridge_validates_charged_sdf_atom_order(tmp_path):
    """Charged SDF sources must match generated-fragment atom order."""
    charged_sdf = tmp_path / "polymer_charged.sdf"
    charged_sdf.write_text(_mini_sdf(("C", "O")), encoding="utf-8")
    fragment = SimpleNamespace(
        atoms=(
            SimpleNamespace(atom_index=0, atom_name="C1", element="C"),
            SimpleNamespace(atom_index=1, atom_name="O1", element="O"),
        )
    )

    charge_bridge._validate_charged_sdf_matches_fragment(
        charged_sdf,
        generated_fragment=fragment,
    )

    mismatched_sdf = tmp_path / "polymer_charged_mismatch.sdf"
    mismatched_sdf.write_text(_mini_sdf(("O", "C")), encoding="utf-8")
    with pytest.raises(ValueError, match="element/order does not match"):
        charge_bridge._validate_charged_sdf_matches_fragment(
            mismatched_sdf,
            generated_fragment=fragment,
        )


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


def _product_atom(
    serial: int,
    atom_index: int,
    chain_id: str,
    residue_name: str,
    residue_number: int,
    atom_name: str,
) -> SimpleNamespace:
    """Build a product PDB atom-record double."""
    return SimpleNamespace(
        serial=serial,
        atom_index=atom_index,
        chain_id=chain_id,
        residue_name=residue_name,
        residue_number=residue_number,
        insertion_code="",
        atom_name=atom_name,
    )


def _atom_from_product_atom(product_atom: SimpleNamespace) -> SimpleNamespace:
    """Build an OpenFF-like atom double with product PDB metadata."""
    atom = _atom(
        product_atom.chain_id,
        product_atom.residue_name,
        product_atom.residue_number,
        product_atom.atom_name,
        0,
    )
    atom.metadata.update(
        {
            "product_identity_source": "product_pdb",
            "product_atom_index": product_atom.atom_index,
            "serial": product_atom.serial,
            "product_atom_serial": product_atom.serial,
            "insertion_code": product_atom.insertion_code,
        }
    )
    return atom


def _mini_sdf(elements: tuple[str, ...]) -> str:
    """Build a minimal V2000 SDF for validation tests."""
    lines = [
        "\n",
        "  PolyzyMD test fixture\n",
        "\n",
        f"{len(elements):3d}{max(len(elements) - 1, 0):3d}  0  0  0  0  0  0  0  0999 V2000\n",
    ]
    for index, element in enumerate(elements):
        lines.append(
            f"{float(index):10.4f}{0.0:10.4f}{0.0:10.4f} {element:<3} 0  0"
            "  0  0  0  0  0  0  0  0  0  0\n"
        )
    for index in range(1, len(elements)):
        lines.append(f"{index:3d}{index + 1:3d}{1:3d}  0\n")
    lines.append("M  END\n$$$$\n")
    return "".join(lines)
