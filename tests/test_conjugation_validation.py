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
    audit_parameter_coverage,
    audit_relaxation_evidence,
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


class RaisingInterchange:
    """Fake Interchange that must not be converted during validation."""

    def to_openmm_system(self) -> FakeSystem:
        """Raise if validation attempts a validation-only OpenMM conversion."""
        raise AssertionError("validation must not call to_openmm_system")


def test_validation_report_serializes_json(tmp_path):
    """Validation reports should write canonical JSON artifacts."""
    pdb_path = tmp_path / "product.pdb"
    _write_product_pdb(pdb_path, include_link=True)
    report = build_conjugate_validation_report(
        product_pdb_path=pdb_path,
        assembly=SimpleNamespace(added_conect_pairs=((1, 2),)),
        output_dir=tmp_path,
        openmm_system=FakeSystem(2),
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


def test_atom_presence_detects_multiple_remapped_modifier_leaving_candidates():
    """Any plausible forbidden remap candidate should count as lingering."""
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
    first_lingering_candidate = _atom(
        serial=2,
        atom_name="H7",
        residue_name="PRD",
        residue_number=99,
        chain_id="C",
    )
    second_lingering_candidate = _atom(
        serial=3,
        atom_name="H7",
        residue_name="PRD",
        residue_number=100,
        chain_id="C",
    )
    plan = SimpleNamespace(
        protein_link_atom=protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(modifier_leaving_atom,),
        modifier_product_residue_name="PRD",
    )
    assembly = SimpleNamespace(
        residue_mappings={
            "fragment_1:4a": {
                "source_residue_number": 4,
                "target_residue_number": 99,
                "target_chain": "C",
            },
            "fragment_1:4b": {
                "source_residue_number": 4,
                "target_residue_number": 100,
                "target_chain": "C",
            },
        }
    )

    report = validate_atom_presence(
        (
            protein_link_atom,
            modifier_link_atom,
            first_lingering_candidate,
            second_lingering_candidate,
        ),
        resolved_plans=(plan,),
        assembly=assembly,
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


def test_atom_presence_scopes_remapped_modifier_link_to_attachment():
    """Identical modifier fragments should match link atoms using plan metadata."""
    first_protein_link_atom = _atom(
        serial=1,
        atom_name="NZ",
        residue_name="LYS",
        residue_number=10,
    )
    second_protein_link_atom = _atom(
        serial=3,
        atom_name="NZ",
        residue_name="LYS",
        residue_number=11,
    )
    modifier_link_atom = _atom(
        serial=50,
        atom_name="C7",
        residue_name="MOD",
        residue_number=4,
    )
    first_product_modifier_link = _atom(
        serial=2,
        atom_name="C7",
        residue_name="PRD",
        residue_number=99,
        chain_id="C",
    )
    second_product_modifier_link = _atom(
        serial=4,
        atom_name="C7",
        residue_name="PRD",
        residue_number=100,
        chain_id="C",
    )
    first_plan = SimpleNamespace(
        protein_link_atom=first_protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(),
        modifier_product_residue_name="PRD",
    )
    second_plan = SimpleNamespace(
        protein_link_atom=second_protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(),
        modifier_product_residue_name="PRD",
    )
    assembly = SimpleNamespace(
        added_conect_pairs=((1, 2), (3, 4)),
        residue_mappings={
            "fragment_1:4": {
                "source_residue_number": 4,
                "target_residue_number": 99,
                "target_chain": "C",
            },
            "fragment_2:4": {
                "source_residue_number": 4,
                "target_residue_number": 100,
                "target_chain": "C",
            },
        },
    )

    report = validate_atom_presence(
        (
            first_protein_link_atom,
            first_product_modifier_link,
            second_protein_link_atom,
            second_product_modifier_link,
        ),
        resolved_plans=(first_plan, second_plan),
        assembly=assembly,
    )

    assert report.status == ValidationStatus.PASS
    assert report.missing_atoms == ()


def test_atom_presence_keeps_duplicate_modifier_identities_plan_scoped():
    """A remapped match for one attachment should not satisfy another."""
    first_protein_link_atom = _atom(
        serial=1,
        atom_name="NZ",
        residue_name="LYS",
        residue_number=10,
    )
    second_protein_link_atom = _atom(
        serial=3,
        atom_name="NZ",
        residue_name="LYS",
        residue_number=11,
    )
    modifier_link_atom = _atom(
        serial=50,
        atom_name="C7",
        residue_name="MOD",
        residue_number=4,
    )
    first_product_modifier_link = _atom(
        serial=2,
        atom_name="C7",
        residue_name="PRD",
        residue_number=99,
        chain_id="C",
    )
    first_plan = SimpleNamespace(
        protein_link_atom=first_protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(),
        modifier_product_residue_name="PRD",
    )
    second_plan = SimpleNamespace(
        protein_link_atom=second_protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=(),
        modifier_leaving_atoms=(),
        modifier_product_residue_name="PRD",
    )
    assembly = SimpleNamespace(
        added_conect_pairs=((1, 2), (3, 999)),
        residue_mappings={
            "fragment_1:4": {
                "source_residue_number": 4,
                "target_residue_number": 99,
                "target_chain": "C",
            },
            "fragment_2:4": {
                "source_residue_number": 4,
                "target_residue_number": 100,
                "target_chain": "C",
            },
        },
    )

    report = validate_atom_presence(
        (first_protein_link_atom, first_product_modifier_link, second_protein_link_atom),
        resolved_plans=(first_plan, second_plan),
        assembly=assembly,
    )

    assert report.status == ValidationStatus.FAIL
    assert len(report.missing_atoms) == 1
    missing_modifier_identity = report.missing_atoms[0]
    assert report.present_atoms.count(missing_modifier_identity) == 1
    assert missing_modifier_identity.atom_name == "C7"


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
    """Parameter coverage should use production OpenMM evidence without heavy imports."""
    passed = audit_parameter_coverage(openmm_system=FakeSystem(5), expected_particle_count=5)
    failed = audit_parameter_coverage(openmm_system=FakeSystem(4), expected_particle_count=5)

    assert passed.status == ValidationStatus.PASS
    assert failed.status == ValidationStatus.FAIL
    assert failed.observed_particle_count == 4


def test_parameter_coverage_skips_without_production_system_evidence():
    """Parameter coverage should skip when production OpenMM evidence is unavailable."""
    report = audit_parameter_coverage(expected_particle_count=5)

    assert report.status == ValidationStatus.SKIPPED
    assert "No production OpenMM System evidence" in report.checks[0].message


def test_parameter_coverage_does_not_convert_interchange_like_evidence():
    """Validation reporting must not perform validation-only Interchange conversion."""
    report = audit_parameter_coverage(
        openmm_system=RaisingInterchange(),
        expected_particle_count=5,
    )

    assert report.status == ValidationStatus.SKIPPED
    assert "getNumParticles" in report.checks[0].message


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


def test_relaxation_evidence_audit_pass_fail_and_skipped(tmp_path):
    """OpenMM relaxation audit should consume relaxation JSON payloads."""
    assert audit_relaxation_evidence(tmp_path).status == ValidationStatus.SKIPPED

    relaxation_path = tmp_path / "conjugate_relaxation.json"
    relaxation_path.write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": True,
                "stage_b_success": True,
                "barostat_used": False,
                "stage_a_energy_before_min_kj_mol": -0.5,
                "stage_a_energy_after_min_kj_mol": -1.0,
                "stage_b_energy_before_md_kj_mol": -1.5,
                "stage_b_energy_after_md_kj_mol": -2.0,
                "stage_b_protein_rmsd_from_stage_a_angstrom": 0.0,
                "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.0,
                "stage_b_linkage_distance_errors_angstrom": [0.1],
                "settings": {"md_steps": 10},
            }
        ),
        encoding="utf-8",
    )
    assert audit_relaxation_evidence(tmp_path).status == ValidationStatus.PASS

    relaxation_path.write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": True,
                "stage_b_success": True,
                "barostat_used": False,
                "stage_a_energy_after_min_kj_mol": "nan",
                "stage_b_energy_after_md_kj_mol": -2.0,
                "settings": {"md_steps": 10},
            }
        ),
        encoding="utf-8",
    )
    assert audit_relaxation_evidence(tmp_path).status == ValidationStatus.FAIL


def test_relaxation_evidence_audit_ignores_legacy_validation_artifacts(tmp_path):
    """Legacy validation artifacts should not affect the new relaxation audit."""
    (tmp_path / "legacy_validation_evidence.json").write_text(
        json.dumps({"success": False}),
        encoding="utf-8",
    )
    diagnostics_path = tmp_path / "conjugate_relaxation.json"
    diagnostics_path.write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": True,
                "stage_b_success": True,
                "temporary_anchor_count": 2,
                "barostat_used": False,
                "stage_a_energy_before_min_kj_mol": -9.0,
                "stage_a_energy_after_min_kj_mol": -10.0,
                "stage_b_energy_before_md_kj_mol": -10.5,
                "stage_b_energy_after_md_kj_mol": -11.0,
                "stage_b_protein_rmsd_from_stage_a_angstrom": 0.0,
                "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.0,
                "stage_b_linkage_distance_errors_angstrom": [0.1, 0.2],
                "settings": {
                    "md_steps": 10,
                    "max_protein_rmsd_angstrom": 0.05,
                    "max_protein_displacement_angstrom": 0.25,
                    "max_linkage_distance_error_angstrom": 0.35,
                },
            }
        ),
        encoding="utf-8",
    )

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.PASS
    assert report.relaxation_diagnostics_json_path == diagnostics_path


def test_relaxation_evidence_audit_accepts_passing_conjugate_relaxation(tmp_path):
    """Validation should pass when all present OpenMM relaxation evidence passes."""
    diagnostics_path = tmp_path / "conjugate_relaxation.json"
    diagnostics_path.write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": True,
                "stage_b_success": True,
                "temporary_anchor_count": 2,
                "barostat_used": False,
                "stage_a_energy_before_min_kj_mol": -9.0,
                "stage_a_energy_after_min_kj_mol": -10.0,
                "stage_b_energy_before_md_kj_mol": -10.5,
                "stage_b_energy_after_md_kj_mol": -11.0,
                "stage_b_protein_rmsd_from_stage_a_angstrom": 0.0,
                "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.0,
                "stage_b_linkage_distance_errors_angstrom": [0.1, 0.2],
                "settings": {
                    "md_steps": 10,
                    "max_protein_rmsd_angstrom": 0.05,
                    "max_protein_displacement_angstrom": 0.25,
                    "max_linkage_distance_error_angstrom": 0.35,
                },
            }
        ),
        encoding="utf-8",
    )

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.PASS
    assert report.relaxation_diagnostics_json_path == diagnostics_path


def test_relaxation_evidence_audit_accepts_force_group_energy_maps(tmp_path):
    """Validation should pass finite scalar and force-group relaxation energies."""
    diagnostics_path = tmp_path / "conjugate_relaxation.json"
    diagnostics_path.write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": True,
                "stage_b_success": True,
                "temporary_anchor_count": 2,
                "barostat_used": False,
                "stage_a_energy_before_min_kj_mol": -9.0,
                "stage_a_energy_after_min_kj_mol": -10.0,
                "stage_b_energy_before_md_kj_mol": -10.5,
                "stage_b_energy_after_md_kj_mol": -11.0,
                "stage_a_force_group_energies_before_min_kj_mol": {
                    "0": -9.0,
                    "1": 0.5,
                    "2": 1.5,
                },
                "stage_a_force_group_energies_after_min_kj_mol": {
                    "bonded": -10.5,
                    "nonbonded": -0.25,
                },
                "stage_b_force_group_energies_after_md_kj_mol": {
                    "bonded": -11.25,
                    "nonbonded": 0.25,
                },
                "relaxation_note": "finite force-group energies should not fail audit",
                "stage_b_protein_rmsd_from_stage_a_angstrom": 0.0,
                "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.0,
                "stage_b_linkage_distance_errors_angstrom": [0.1, 0.2],
                "settings": {
                    "md_steps": 10,
                    "max_protein_rmsd_angstrom": 0.05,
                    "max_protein_displacement_angstrom": 0.25,
                    "max_linkage_distance_error_angstrom": 0.35,
                },
            }
        ),
        encoding="utf-8",
    )

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.PASS
    assert report.relaxation_diagnostics_json_path == diagnostics_path


def test_relaxation_evidence_audit_requires_stage_b_rmsd(tmp_path):
    """Validation should reject successful diagnostics missing Stage B protein RMSD."""
    payload = _canonical_relaxation_payload()
    payload.pop("stage_b_protein_rmsd_from_stage_a_angstrom")
    _write_relaxation_diagnostics(tmp_path, payload)

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(
        check.name == "conjugate_relaxation_required_protein_immobilization"
        and "stage_b_protein_rmsd_from_stage_a_angstrom" in check.evidence["fields"]
        for check in report.checks
    )


def test_relaxation_evidence_audit_requires_stage_b_max_displacement(tmp_path):
    """Validation should reject successful diagnostics missing Stage B max displacement."""
    payload = _canonical_relaxation_payload()
    payload.pop("stage_b_protein_max_displacement_from_stage_a_angstrom")
    _write_relaxation_diagnostics(tmp_path, payload)

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(
        check.name == "conjugate_relaxation_required_protein_immobilization"
        and "stage_b_protein_max_displacement_from_stage_a_angstrom" in check.evidence["fields"]
        for check in report.checks
    )


def test_relaxation_evidence_audit_requires_stage_b_linkage_distance_errors(tmp_path):
    """Validation should reject successful diagnostics missing linkage distance errors."""
    payload = _canonical_relaxation_payload()
    payload.pop("stage_b_linkage_distance_errors_angstrom")
    _write_relaxation_diagnostics(tmp_path, payload)

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(
        check.name == "conjugate_relaxation_required_linkage_distances"
        and check.evidence["field"] == "stage_b_linkage_distance_errors_angstrom"
        for check in report.checks
    )


def test_relaxation_evidence_audit_requires_canonical_energy_fields(tmp_path):
    """Validation should reject successful diagnostics missing canonical energy fields."""
    required_fields = (
        "stage_a_energy_before_min_kj_mol",
        "stage_a_energy_after_min_kj_mol",
        "stage_b_energy_before_md_kj_mol",
        "stage_b_energy_after_md_kj_mol",
    )
    for field in required_fields:
        payload = _canonical_relaxation_payload()
        payload.pop(field)
        _write_relaxation_diagnostics(tmp_path, payload)

        report = audit_relaxation_evidence(tmp_path)

        assert report.status == ValidationStatus.FAIL
        assert any(
            check.name == "conjugate_relaxation_required_energies"
            and field in check.evidence["fields"]
            for check in report.checks
        )


def test_relaxation_evidence_audit_rejects_stale_zero_step_relaxation(tmp_path):
    """Validation should reject stale frozen diagnostics with no Stage B MD."""
    (tmp_path / "conjugate_relaxation.json").write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": True,
                "stage_b_success": True,
                "barostat_used": False,
                "stage_a_energy_before_min_kj_mol": -9.0,
                "stage_a_energy_after_min_kj_mol": -10.0,
                "stage_b_energy_before_md_kj_mol": -10.5,
                "stage_b_energy_after_md_kj_mol": -11.0,
                "stage_b_protein_rmsd_from_stage_a_angstrom": 0.0,
                "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.0,
                "stage_b_linkage_distance_errors_angstrom": [0.1],
                "settings": {"md_steps": 0},
            }
        ),
        encoding="utf-8",
    )

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(check.name == "conjugate_relaxation_md_steps" for check in report.checks)


def test_relaxation_evidence_audit_requires_conjugate_relaxation_md_steps(tmp_path):
    """Validation should reject otherwise-passing frozen diagnostics without MD steps."""
    (tmp_path / "conjugate_relaxation.json").write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": True,
                "stage_b_success": True,
                "barostat_used": False,
                "stage_a_energy_before_min_kj_mol": -9.0,
                "stage_a_energy_after_min_kj_mol": -10.0,
                "stage_b_energy_before_md_kj_mol": -10.5,
                "stage_b_energy_after_md_kj_mol": -11.0,
                "stage_b_protein_rmsd_from_stage_a_angstrom": 0.0,
                "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.0,
                "stage_b_linkage_distance_errors_angstrom": [0.1],
                "settings": {
                    "max_protein_rmsd_angstrom": 0.05,
                    "max_protein_displacement_angstrom": 0.25,
                    "max_linkage_distance_error_angstrom": 0.35,
                },
            }
        ),
        encoding="utf-8",
    )

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(check.name == "conjugate_relaxation_md_steps" for check in report.checks)


def test_relaxation_evidence_audit_accepts_top_level_relaxation_md_steps(tmp_path):
    """Validation should accept top-level MD steps when settings omit them."""
    (tmp_path / "conjugate_relaxation.json").write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": True,
                "stage_b_success": True,
                "barostat_used": False,
                "stage_a_energy_before_min_kj_mol": -9.0,
                "stage_a_energy_after_min_kj_mol": -10.0,
                "stage_b_energy_before_md_kj_mol": -10.5,
                "stage_b_energy_after_md_kj_mol": -11.0,
                "stage_b_protein_rmsd_from_stage_a_angstrom": 0.0,
                "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.0,
                "stage_b_linkage_distance_errors_angstrom": [0.1],
                "md_steps": 10,
                "settings": {
                    "max_protein_rmsd_angstrom": 0.05,
                    "max_protein_displacement_angstrom": 0.25,
                    "max_linkage_distance_error_angstrom": 0.35,
                },
            }
        ),
        encoding="utf-8",
    )

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.PASS


def test_relaxation_evidence_audit_fails_non_numeric_tolerance_setting(tmp_path):
    """Validation should fail invalid tolerance settings without raising."""
    payload = _canonical_relaxation_payload()
    settings = payload["settings"]
    assert isinstance(settings, dict)
    settings["max_protein_rmsd_angstrom"] = "not-a-number"
    _write_relaxation_diagnostics(tmp_path, payload)

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(
        check.name == "conjugate_relaxation_tolerance_settings"
        and check.evidence["field"] == "max_protein_rmsd_angstrom"
        for check in report.checks
    )


def test_relaxation_evidence_audit_fails_negative_rmsd(tmp_path):
    """Validation should reject negative Stage B protein RMSD evidence."""
    payload = _canonical_relaxation_payload()
    payload["stage_b_protein_rmsd_from_stage_a_angstrom"] = -0.01
    _write_relaxation_diagnostics(tmp_path, payload)

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(
        check.name == "conjugate_relaxation_required_protein_immobilization"
        and "stage_b_protein_rmsd_from_stage_a_angstrom" in check.evidence["fields"]
        for check in report.checks
    )


def test_relaxation_evidence_audit_fails_negative_max_displacement(tmp_path):
    """Validation should reject negative Stage B protein max displacement evidence."""
    payload = _canonical_relaxation_payload()
    payload["stage_b_protein_max_displacement_from_stage_a_angstrom"] = -0.01
    _write_relaxation_diagnostics(tmp_path, payload)

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(
        check.name == "conjugate_relaxation_required_protein_immobilization"
        and "stage_b_protein_max_displacement_from_stage_a_angstrom" in check.evidence["fields"]
        for check in report.checks
    )


def test_relaxation_evidence_audit_fails_negative_linkage_error(tmp_path):
    """Validation should reject negative linkage distance error evidence."""
    payload = _canonical_relaxation_payload()
    payload["stage_b_linkage_distance_errors_angstrom"] = [0.1, -0.01]
    _write_relaxation_diagnostics(tmp_path, payload)

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(
        check.name == "conjugate_relaxation_required_linkage_distances"
        and check.evidence["field"] == "stage_b_linkage_distance_errors_angstrom"
        for check in report.checks
    )


def test_relaxation_evidence_audit_fails_unfixed_frozen_stage_b(tmp_path):
    """Validation should fail when Stage B does not keep protein fixed."""
    diagnostics_path = tmp_path / "conjugate_relaxation.json"
    diagnostics_path.write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": True,
                "stage_b_success": True,
                "barostat_used": False,
                "stage_a_energy_before_min_kj_mol": -9.0,
                "stage_a_energy_after_min_kj_mol": -10.0,
                "stage_b_energy_before_md_kj_mol": -10.5,
                "stage_b_energy_after_md_kj_mol": -11.0,
                "stage_b_protein_rmsd_from_stage_a_angstrom": 0.2,
                "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.5,
                "stage_b_linkage_distance_errors_angstrom": [0.1],
                "settings": {
                    "max_protein_rmsd_angstrom": 0.05,
                    "max_protein_displacement_angstrom": 0.25,
                    "max_linkage_distance_error_angstrom": 0.35,
                },
            }
        ),
        encoding="utf-8",
    )

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(check.name == "conjugate_relaxation_protein_rmsd" for check in report.checks)


def test_relaxation_evidence_audit_requires_frozen_stage_a_success(tmp_path):
    """Validation should fail when Stage A minimization is missing."""
    diagnostics_path = tmp_path / "conjugate_relaxation.json"
    diagnostics_path.write_text(
        json.dumps(
            {
                "success": True,
                "stage_a_success": False,
                "stage_b_success": True,
                "barostat_used": False,
                "stage_a_energy_before_min_kj_mol": -9.0,
                "stage_a_energy_after_min_kj_mol": -10.0,
                "stage_b_energy_before_md_kj_mol": -10.5,
                "stage_b_energy_after_md_kj_mol": -11.0,
                "stage_b_protein_rmsd_from_stage_a_angstrom": 0.0,
                "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.0,
                "stage_b_linkage_distance_errors_angstrom": [0.1],
            }
        ),
        encoding="utf-8",
    )

    report = audit_relaxation_evidence(tmp_path)

    assert report.status == ValidationStatus.FAIL
    assert any(check.name == "conjugate_relaxation_stage_a" for check in report.checks)


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


def _canonical_relaxation_payload() -> dict[str, object]:
    """Return otherwise-passing conjugate relaxation diagnostics."""
    return {
        "success": True,
        "stage_a_success": True,
        "stage_b_success": True,
        "barostat_used": False,
        "stage_a_energy_before_min_kj_mol": -9.0,
        "stage_a_energy_after_min_kj_mol": -10.0,
        "stage_b_energy_before_md_kj_mol": -10.5,
        "stage_b_energy_after_md_kj_mol": -11.0,
        "stage_b_protein_rmsd_from_stage_a_angstrom": 0.0,
        "stage_b_protein_max_displacement_from_stage_a_angstrom": 0.0,
        "stage_b_linkage_distance_errors_angstrom": [0.1, 0.2],
        "settings": {
            "md_steps": 10,
            "max_protein_rmsd_angstrom": 0.05,
            "max_protein_displacement_angstrom": 0.25,
            "max_linkage_distance_error_angstrom": 0.35,
        },
    }


def _write_relaxation_diagnostics(tmp_path: Path, payload: dict[str, object]) -> None:
    """Write conjugate relaxation diagnostics for validation tests."""
    (tmp_path / "conjugate_relaxation.json").write_text(
        json.dumps(payload),
        encoding="utf-8",
    )
