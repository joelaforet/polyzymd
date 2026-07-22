"""Opt-in fast mixed polymer/glycan conjugation E2E checkpoint."""

from __future__ import annotations

import os

import pytest

from tests._support.conjugation_fast_e2e import (
    run_fast_mixed_build_export,
    validate_fast_mixed_o_summary,
    validate_fast_mixed_summary,
)

RUN_FAST_E2E = "POLYZYMD_RUN_CONJUGATION_FAST_E2E"
RUN_FAST_E2E_CUDA = "POLYZYMD_RUN_CONJUGATION_FAST_E2E_CUDA"
RUN_FAST_O_E2E = "POLYZYMD_RUN_CONJUGATION_FAST_O_E2E"
RUN_FAST_O_THR_E2E = "POLYZYMD_RUN_CONJUGATION_FAST_O_THR_E2E"


@pytest.mark.slow
def test_fast_mixed_conjugation_build_export_checkpoint(tmp_path):
    """Build and export the fast mixed 1UBQ checkpoint when opted in."""
    if os.environ.get(RUN_FAST_E2E) != "1":
        pytest.skip(f"Set {RUN_FAST_E2E}=1 to run the fast mixed conjugation E2E checkpoint")

    summary = run_fast_mixed_build_export(tmp_path)

    validate_fast_mixed_summary(summary)


@pytest.mark.slow
def test_fast_mixed_conjugation_cuda_checkpoint(tmp_path):
    """Run the separate 100000-step CUDA checkpoint only when explicitly opted in."""
    if os.environ.get(RUN_FAST_E2E_CUDA) != "1":
        pytest.skip(f"Set {RUN_FAST_E2E_CUDA}=1 to run the CUDA fast mixed checkpoint")

    summary = run_fast_mixed_build_export(tmp_path, run_cuda=True)

    validate_fast_mixed_summary(summary, require_cuda=True)


@pytest.mark.slow
def test_fast_mixed_o_glycosylation_build_export_checkpoint(tmp_path):
    """Build the mixed NHS-Lys plus Ser O-glycosylation checkpoint when opted in."""
    if os.environ.get(RUN_FAST_O_E2E) != "1":
        pytest.skip(f"Set {RUN_FAST_O_E2E}=1 to run the O-glycosylation E2E checkpoint")

    summary = run_fast_mixed_build_export(tmp_path, glycosylation="o_ser")

    validate_fast_mixed_o_summary(summary)


@pytest.mark.slow
def test_fast_mixed_thr_o_glycosylation_build_export_checkpoint(tmp_path):
    """Build the mixed NHS-Lys plus Thr O-glycosylation checkpoint when opted in."""
    if os.environ.get(RUN_FAST_O_THR_E2E) != "1":
        pytest.skip(f"Set {RUN_FAST_O_THR_E2E}=1 to run the Thr O-glycosylation E2E")

    summary = run_fast_mixed_build_export(tmp_path, glycosylation="o_thr")

    validate_fast_mixed_o_summary(
        summary, product_residue="OLT", residue_number=60, site_atom="OG1"
    )


def test_fast_mixed_summary_validation_accepts_expected_payload():
    """Checkpoint validation should accept a compact passing payload."""
    validate_fast_mixed_summary(_summary_payload())


def test_fast_mixed_summary_validation_rejects_mono_nag_surrogate():
    """Checkpoint validation should reject mono-NAG-like glycan evidence."""
    summary = _summary_payload()
    summary["fixture_glycan"]["heavy_atom_count"] = 14

    with pytest.raises(AssertionError):
        validate_fast_mixed_summary(summary)


def test_fast_mixed_summary_validation_checks_cuda_when_requested():
    """CUDA validation should be separate from the CPU build/export contract."""
    summary = _summary_payload()
    summary["cuda"] = {
        "platform": "CUDA",
        "steps": 100_000,
        "checkpoint_exists": True,
        "potential_energy_kj_mol": -100.0,
    }

    validate_fast_mixed_summary(summary, require_cuda=True)


def test_fast_mixed_summary_validation_rejects_missing_polymer_charge_coverage():
    """Checkpoint validation should require authoritative polymer charge coverage."""
    summary = _summary_payload()
    summary["charge_bridge"]["patch_owned_polymer_atom_count"] = 299

    with pytest.raises(AssertionError):
        validate_fast_mixed_summary(summary)


def _summary_payload() -> dict[str, object]:
    """Return a minimal passing fast mixed checkpoint summary."""
    return {
        "fixture_glycan": {
            "heavy_atom_count": 117,
            "hydrogenated_atom_count": 225,
            "formula": "C64H108N2O51",
            "anomeric_carbon_index": 9,
            "leaving_atom_indices": [10, 126],
            "leaving_atom_count": 2,
        },
        "authoritative_system_created": True,
        "resolved_attachment_count": 2,
        "solvated_atom_count": 2_500,
        "crosslinked_atom_count": 1_500,
        "finite_coordinates": True,
        "residue_counts": {"LYX": 10, "ASX": 8, "NAG": 223},
        "product_glycan": {
            "residue_name": "NAG",
            "source_atom_count": 225,
            "leaving_atom_count": 2,
            "retained_atom_count": 223,
        },
        "n_glycosylation_linkage": {
            "mechanism_name": "n_glycosylation",
            "site_residue_number": 60,
            "protein_product_residue_name": "ASX",
            "modifier_product_residue_name": "NAG",
            "protein_link_atom_name": "ND2",
            "product_site_atom_present": True,
            "protein_leaving_atom_names": ["HD21"],
            "modifier_leaving_atom_count": 2,
            "crosslink_residues": ["ASX", "NAG"],
            "crosslink_linking_atoms": ["ND2", "C1"],
        },
        "export_exists": {"gro": True, "top": True, "em_mdp": True, "prod_mdp": True},
        "validation": {
            "bond_graph_status": "pass",
            "atom_presence_status": "pass",
            "charge_audit_status": "pass",
            "relaxation_evidence_status": "pass",
            "linkage_distances_angstrom": [1.33, 1.45],
            "close_contact_count": 0,
        },
        "charge_bridge": {
            "success": True,
            "ff14sb_atom_count": 1_198,
            "polymer_template_atom_count": 0,
            "local_nagl_patch_atom_count": 302,
            "source_role_atom_count": 1_500,
            "polymer_atom_count": 300,
            "polymer_source_coverage_count": 300,
            "patch_owned_polymer_atom_count": 300,
            "polymer_template_mapped_product_atom_count": 3,
            "polymer_template_overwrite_count": 3,
            "polymer_template_coverage_count": 3,
            "unmodified_polymer_template_atom_count": 0,
            "invalid_polymer_template_source_count": 0,
            "total_charge_e": 0.0,
            "formal_charge_e": 0.0,
            "normalization_correction_e": 0.0,
            "max_per_atom_correction_e": 0.0,
        },
        "relaxation": {
            "stage_a_energy_before_min_kj_mol": -10.0,
            "stage_a_energy_after_min_kj_mol": -11.0,
            "stage_a_linkage_distances_angstrom": [1.33, 1.45],
        },
    }
