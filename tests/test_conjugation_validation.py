"""Tests for conjugate validation reports."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord
from polyzymd.builders.conjugation.validation import (
    ConjugateValidationReport,
    ValidationStatus,
    audit_charge_reports,
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


def test_parameter_coverage_pass_and_fail_with_fakes():
    """Parameter coverage should use fake OpenMM objects without heavy imports."""
    passed = audit_parameter_coverage(FakeInterchange(5), expected_particle_count=5)
    failed = audit_parameter_coverage(FakeInterchange(4), expected_particle_count=5)

    assert passed.status == ValidationStatus.PASS
    assert failed.status == ValidationStatus.FAIL
    assert failed.observed_particle_count == 4


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
) -> PdbAtomRecord:
    """Build a minimal PDB atom record for tests."""
    return PdbAtomRecord(
        serial=serial,
        atom_index=serial - 1,
        atom_name=atom_name,
        residue_name=residue_name,
        chain_id="A" if residue_name == "LYS" else "C",
        residue_number=residue_number,
        x=float(serial - 1),
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
