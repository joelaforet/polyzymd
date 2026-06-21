"""Tests for conjugate validation reports."""

from __future__ import annotations

import json
import math
from pathlib import Path
from types import SimpleNamespace

from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord
from polyzymd.builders.conjugation.validation import (
    ConjugateValidationReport,
    ValidationStatus,
    audit_charge_reports,
    audit_linkage_geometry,
    audit_openmm_smoke_reports,
    audit_parameter_coverage,
    build_conjugate_validation_report,
    validate_atom_presence,
    validate_product_bond_graph,
)


class FakeSystem:
    """Fake OpenMM System exposing a particle count."""

    def __init__(self, particle_count: int) -> None:
        """Initialize the fake system."""
        self._particle_count = particle_count

    def getNumParticles(self) -> int:  # noqa: N802 - OpenMM API compatibility
        """Return the configured fake particle count."""
        return self._particle_count


class FakeInterchange:
    """Fake Interchange exposing OpenMM system conversion."""

    def __init__(self, particle_count: int) -> None:
        """Initialize the fake interchange."""
        self._particle_count = particle_count

    def to_openmm_system(self) -> FakeSystem:
        """Return a fake OpenMM system."""
        return FakeSystem(self._particle_count)


class FailingBackendInterchange:
    """Fake Interchange that raises an expected backend conversion error."""

    def to_openmm_system(self) -> FakeSystem:
        """Raise a backend-like conversion error."""
        raise ValueError("missing parameter")


class BuggyInterchange:
    """Fake Interchange that raises an unexpected programming error."""

    def to_openmm_system(self) -> FakeSystem:
        """Raise an unexpected conversion implementation error."""
        raise AttributeError("buggy object")


def test_validation_report_serializes_json(tmp_path):
    """Validation reports should write canonical JSON artifacts."""
    pdb_path = tmp_path / "product.pdb"
    _write_product_pdb(pdb_path, include_link=True)
    report = build_conjugate_validation_report(
        product_pdb_path=pdb_path,
        assembly=SimpleNamespace(added_conect_pairs=((1, 2),)),
        output_dir=tmp_path,
        interchange=FakeInterchange(2),
        expected_particle_count=2,
        write=True,
    )

    report_path = tmp_path / "conjugate_validation_report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert report_path.exists()
    assert payload["status"] in {"pass", "skipped"}
    assert payload["bond_graph"]["status"] == "pass"
    assert report.report_path == report_path


def test_bond_graph_reports_pass_and_fail():
    """Bond graph validation should fail only when expected bonds are absent."""
    passed = validate_product_bond_graph(((1, 2),), ((1, 2), (2, 3)))
    failed = validate_product_bond_graph(((1, 2),), ((2, 3),))

    assert passed.status == ValidationStatus.PASS
    assert failed.status == ValidationStatus.FAIL
    assert failed.missing_bonds == ((1, 2),)


def test_atom_presence_detects_lingering_leaving_atom():
    """Atom presence validation should fail for lingering leaving atoms."""
    link_atom = _atom(serial=1, atom_name="NZ", residue_name="LYS", residue_number=10)
    modifier_atom = _atom(serial=2, atom_name="C1", residue_name="LIG", residue_number=1)
    leaving_atom = _atom(serial=3, atom_name="HZ1", residue_name="LYS", residue_number=10)
    plan = SimpleNamespace(
        protein_link_atom=link_atom,
        modifier_link_atom=modifier_atom,
        protein_leaving_atoms=(leaving_atom,),
        modifier_leaving_atoms=(),
    )

    report = validate_atom_presence(
        (link_atom, modifier_atom, leaving_atom), resolved_plans=(plan,)
    )

    assert report.status == ValidationStatus.FAIL
    assert len(report.lingering_leaving_atoms) == 1


def test_atom_presence_allows_remapped_modifier_link_atom():
    """Modifier link atom matching should tolerate product PDB identity remapping."""
    protein_link_atom = _atom(
        serial=1,
        atom_name="NZ",
        residue_name="LYS",
        residue_number=10,
    )
    modifier_link_atom = _atom(
        serial=50,
        atom_name="C7",
        residue_name="MOD",
        residue_number=4,
    )
    remapped_modifier_atom = _atom(
        serial=2,
        atom_name="C7",
        residue_name="PRD",
        residue_number=99,
        chain_id="B",
    )
    plan = SimpleNamespace(
        protein_link_atom=protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(),
    )

    report = validate_atom_presence(
        (protein_link_atom, remapped_modifier_atom),
        resolved_plans=(plan,),
    )

    assert report.status == ValidationStatus.PASS
    assert report.missing_atoms == ()


def test_atom_presence_detects_remapped_modifier_leaving_atom():
    """Modifier leaving atoms should be found after product residue remapping."""
    protein_link_atom = _atom(
        serial=1,
        atom_name="NZ",
        residue_name="LYS",
        residue_number=10,
    )
    modifier_link_atom = _atom(
        serial=50,
        atom_name="C7",
        residue_name="MOD",
        residue_number=4,
    )
    modifier_leaving_atom = _atom(
        serial=51,
        atom_name="H7",
        residue_name="MOD",
        residue_number=4,
    )
    remapped_modifier_link_atom = _atom(
        serial=2,
        atom_name="C7",
        residue_name="PRD",
        residue_number=99,
        chain_id="C",
    )
    remapped_leaving_atom = _atom(
        serial=3,
        atom_name="H7",
        residue_name="PRD",
        residue_number=99,
        chain_id="C",
    )
    plan = SimpleNamespace(
        protein_link_atom=protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(modifier_leaving_atom,),
        modifier_product_residue_name="PRD",
    )

    report = validate_atom_presence(
        (protein_link_atom, remapped_modifier_link_atom, remapped_leaving_atom),
        resolved_plans=(plan,),
    )

    assert report.status == ValidationStatus.FAIL
    assert len(report.lingering_leaving_atoms) == 1
    assert report.lingering_leaving_atoms[0].atom_name == "H7"


def test_atom_presence_rejects_ambiguous_remapped_modifier_link_atom():
    """Modifier link matching should not pass duplicate name/element candidates."""
    protein_link_atom = _atom(
        serial=1,
        atom_name="NZ",
        residue_name="LYS",
        residue_number=10,
    )
    modifier_link_atom = _atom(
        serial=50,
        atom_name="C7",
        residue_name="MOD",
        residue_number=4,
    )
    first_candidate = _atom(
        serial=2,
        atom_name="C7",
        residue_name="PR1",
        residue_number=99,
        chain_id="C",
    )
    second_candidate = _atom(
        serial=3,
        atom_name="C7",
        residue_name="PR2",
        residue_number=100,
        chain_id="C",
    )
    plan = SimpleNamespace(
        protein_link_atom=protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(),
    )

    report = validate_atom_presence(
        (protein_link_atom, first_candidate, second_candidate),
        resolved_plans=(plan,),
    )

    assert report.status == ValidationStatus.FAIL
    assert len(report.missing_atoms) == 1
    assert report.missing_atoms[0].atom_name == "C7"


def test_charge_audit_pass_warn_and_fail(tmp_path):
    """Charge audit should summarize fake bridge payloads."""
    bridge_path = tmp_path / "product_state_charge_bridge.json"
    bridge_path.write_text(
        json.dumps({"success": True, "total_charge_e": 0.0, "formal_charge_e": 0.0}),
        encoding="utf-8",
    )
    assert audit_charge_reports(tmp_path).status == ValidationStatus.PASS

    bridge_path.write_text(
        json.dumps({"success": True, "total_charge_e": 0.1, "formal_charge_e": 0.0}),
        encoding="utf-8",
    )
    assert audit_charge_reports(tmp_path).status == ValidationStatus.WARN

    bridge_path.write_text(
        json.dumps({"success": False, "total_charge_e": 0.0, "formal_charge_e": 0.0}),
        encoding="utf-8",
    )
    assert audit_charge_reports(tmp_path).status == ValidationStatus.FAIL


def test_charge_audit_checks_all_numeric_charge_fields(tmp_path):
    """Charge audit should reject non-finite correction fields as well as totals."""
    bridge_path = tmp_path / "product_state_charge_bridge.json"
    bridge_path.write_text(
        json.dumps(
            {
                "success": True,
                "total_charge_e": 0.0,
                "formal_charge_e": 0.0,
                "normalization_correction_e": math.inf,
                "max_per_atom_correction_e": "nan",
            }
        ),
        encoding="utf-8",
    )

    report = audit_charge_reports(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert report.normalization_correction_e is None
    assert report.max_per_atom_correction_e is None
    assert report.checks[0].evidence["fields"] == (
        "normalization_correction_e",
        "max_per_atom_correction_e",
    )


def test_validation_report_writes_strict_json_for_nonfinite_values(tmp_path):
    """Validation report JSON should not emit non-standard NaN tokens."""
    pdb_path = tmp_path / "product.pdb"
    _write_product_pdb(pdb_path, include_link=True)
    (tmp_path / "product_state_charge_bridge.json").write_text(
        json.dumps(
            {
                "success": True,
                "total_charge_e": math.nan,
                "formal_charge_e": 0.0,
            }
        ),
        encoding="utf-8",
    )

    build_conjugate_validation_report(
        product_pdb_path=pdb_path,
        assembly=SimpleNamespace(added_conect_pairs=((1, 2),)),
        output_dir=tmp_path,
        write=True,
    )

    report_text = (tmp_path / "conjugate_validation_report.json").read_text(encoding="utf-8")
    payload = json.loads(report_text)
    assert "NaN" not in report_text
    assert payload["charge_audit"]["total_charge_e"] is None
    assert payload["charge_audit"]["status"] == "fail"


def test_parameter_coverage_pass_and_fail_with_fakes():
    """Parameter coverage should use fake OpenMM objects without heavy imports."""
    passed = audit_parameter_coverage(FakeInterchange(5), expected_particle_count=5)
    failed = audit_parameter_coverage(FakeInterchange(4), expected_particle_count=5)

    assert passed.status == ValidationStatus.PASS
    assert failed.status == ValidationStatus.FAIL
    assert failed.observed_particle_count == 4


def test_parameter_coverage_classifies_expected_backend_errors():
    """Expected conversion errors should be reported as parameter coverage failures."""
    report = audit_parameter_coverage(FailingBackendInterchange(), expected_particle_count=5)

    assert report.status == ValidationStatus.FAIL
    assert report.checks[0].evidence["error"] == "missing parameter"


def test_parameter_coverage_reraises_unexpected_errors():
    """Unexpected programming errors should not be hidden as validation evidence."""
    try:
        audit_parameter_coverage(BuggyInterchange(), expected_particle_count=5)
    except AttributeError as exc:
        assert str(exc) == "buggy object"
    else:  # pragma: no cover - explicit failure branch improves assertion message
        raise AssertionError("Unexpected conversion errors should be re-raised")


def test_linkage_geometry_warns_for_nonbonded_close_contacts():
    """Close contacts below the clash threshold should affect geometry status."""
    atoms = (
        _atom(serial=1, atom_name="N1", residue_name="LYS", residue_number=10),
        _atom(serial=2, atom_name="C1", residue_name="MOD", residue_number=1),
        _atom(serial=3, atom_name="C2", residue_name="MOD", residue_number=1, x=0.2),
    )

    report = audit_linkage_geometry(atoms, ((1, 2),), ((1, 2),))

    assert report.status == ValidationStatus.WARN
    assert report.close_contact_count == 1
    assert any(check.name == "nonbonded_close_contacts" for check in report.checks)


def test_smoke_audit_pass_fail_and_skipped(tmp_path):
    """Smoke audit should consume fake JSON payloads."""
    assert audit_openmm_smoke_reports(tmp_path).status == ValidationStatus.SKIPPED

    smoke_path = tmp_path / "vacuum_smoke.json"
    smoke_path.write_text(
        json.dumps({"success": True, "energy_after_min_kj_mol": -1.0}),
        encoding="utf-8",
    )
    assert audit_openmm_smoke_reports(tmp_path).status == ValidationStatus.PASS

    smoke_path.write_text(
        json.dumps({"success": True, "energy_after_min_kj_mol": "nan"}),
        encoding="utf-8",
    )
    assert audit_openmm_smoke_reports(tmp_path).status == ValidationStatus.FAIL


def test_validation_report_status_aggregation(tmp_path):
    """Top-level report status should aggregate child failures."""
    pdb_path = tmp_path / "product.pdb"
    _write_product_pdb(pdb_path, include_link=False)

    report = build_conjugate_validation_report(
        product_pdb_path=pdb_path,
        assembly=SimpleNamespace(added_conect_pairs=((1, 2),)),
        output_dir=tmp_path,
        write=False,
    )

    assert isinstance(report, ConjugateValidationReport)
    assert report.status == ValidationStatus.FAIL
    assert report.bond_graph.status == ValidationStatus.FAIL


def _atom(
    *,
    serial: int,
    atom_name: str,
    residue_name: str,
    residue_number: int,
    chain_id: str | None = None,
    x: float | None = None,
) -> PdbAtomRecord:
    """Build a minimal PDB atom record for tests."""
    return PdbAtomRecord(
        serial=serial,
        atom_index=serial - 1,
        atom_name=atom_name,
        residue_name=residue_name,
        chain_id=chain_id if chain_id is not None else "A" if residue_name == "LYS" else "C",
        residue_number=residue_number,
        x=float(serial - 1) if x is None else x,
        y=0.0,
        z=0.0,
        element=atom_name[0],
        record_name="ATOM" if residue_name == "LYS" else "HETATM",
    )


def _write_product_pdb(path: Path, *, include_link: bool) -> None:
    """Write a tiny product PDB with optional CONECT linkage."""
    lines = [
        "ATOM      1  NZ  LYS A  10       0.000   0.000   0.000  1.00  0.00           N\n",
        "HETATM    2  C1  LIG C   1       1.400   0.000   0.000  1.00  0.00           C\n",
    ]
    if include_link:
        lines.extend(["CONECT    1    2\n", "CONECT    2    1\n"])
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")
