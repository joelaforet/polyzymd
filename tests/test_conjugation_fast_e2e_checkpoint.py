"""Opt-in fast mixed polymer/glycan conjugation E2E checkpoint."""

from __future__ import annotations

import os

import pytest

from tests._support.conjugation_fast_e2e import (
    run_fast_mixed_build_export,
    validate_fast_mixed_summary,
)

RUN_FAST_E2E = "POLYZYMD_RUN_CONJUGATION_FAST_E2E"
RUN_FAST_E2E_CUDA = "POLYZYMD_RUN_CONJUGATION_FAST_E2E_CUDA"


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


def _summary_payload() -> dict[str, object]:
    """Return a minimal passing fast mixed checkpoint summary."""
    return {
        "fixture_glycan": {
            "heavy_atom_count": 117,
            "hydrogenated_atom_count": 225,
            "formula": "C64H108N2O51",
            "anomeric_carbon_index": 9,
            "leaving_atom_indices": [10, 126],
        },
        "final_interchange_created": True,
        "resolved_attachment_count": 2,
        "solvated_atom_count": 2_500,
        "crosslinked_atom_count": 1_500,
        "finite_coordinates": True,
        "residue_counts": {"LYX": 10, "ASX": 8, "NAG": 225},
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
            "ff14sb_atom_count": 1_100,
            "polymer_template_atom_count": 3,
            "local_nagl_patch_atom_count": 2,
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
