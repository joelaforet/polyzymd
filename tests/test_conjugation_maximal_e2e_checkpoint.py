"""Focused tests and opt-in maximal mixed-conjugation acceptance."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from polyzymd.config.loader import load_config
from tests._support.conjugation_maximal_e2e import (
    MaximalMdProtocol,
    composition_evidence,
    maximal_config_payload,
    run_maximal_mixed_acceptance,
    validate_maximal_manifest,
    write_maximal_runtime_config,
)

RUN_MAXIMAL_E2E = "POLYZYMD_RUN_CONJUGATION_MAXIMAL_E2E"


@pytest.mark.slow
def test_maximal_mixed_conjugation_cuda_acceptance(tmp_path: Path) -> None:
    """Run only with the opt-in flag using the approved maximal protocol."""
    if os.environ.get(RUN_MAXIMAL_E2E) != "1":
        pytest.skip(f"Set {RUN_MAXIMAL_E2E}=1 to run the maximal acceptance case")
    protocol = MaximalMdProtocol()

    manifest = run_maximal_mixed_acceptance(tmp_path, protocol=protocol)

    assert manifest["cuda"]["frame_count"] == protocol.expected_frames


def test_maximal_protocol_requires_explicit_divisible_steps() -> None:
    """Frame count must be exact and derived from explicit caller inputs."""
    protocol = MaximalMdProtocol()

    assert protocol.expected_frames == 150
    with pytest.raises(ValueError, match="divisible"):
        MaximalMdProtocol(production_steps=150_001)


def test_maximal_protocol_rejects_superseded_conditions() -> None:
    """The harness must fail closed on the superseded 1 fs/300 K protocol."""
    with pytest.raises(ValueError, match="fixed and approved"):
        MaximalMdProtocol(
            time_step_fs=1.0,
            temperature_k=300.0,
        )


def test_maximal_config_payload_is_explicit_and_simultaneous(tmp_path: Path) -> None:
    """Payload generation should encode every maximal build-domain choice."""
    protocol = MaximalMdProtocol()

    payload = maximal_config_payload(tmp_path, protocol)

    assert payload["force_field"] == {
        "protein": "ff14sb_off_impropers_0.0.4.offxml",
        "small_molecule": "openff-2.3.0.offxml",
    }
    assert [item["moiety"]["force_field"] for item in payload["conjugation"]["attachments"]] == [
        "openff-2.3.0.offxml",
        "glycam06",
        "glycam06",
    ]
    assert [item["site"]["atom_name"] for item in payload["conjugation"]["attachments"]] == [
        "NZ",
        "ND2",
        "OG",
    ]
    assert payload["polymers"]["count"] == 3
    assert payload["polymers"]["length"] == 3
    assert payload["solvent"]["co_solvents"] == [{"name": "dmso", "mole_fraction": 0.2}]
    assert payload["solvent"]["ions"] == {"neutralize": True, "nacl_concentration": 0.0}
    assert payload["thermodynamics"] == {"temperature": 293.0, "pressure": 1.0}
    assert payload["simulation_phases"]["production"]["time_step"] == 2.0
    assert payload["simulation_phases"]["production"]["report_interval"] == 1_000
    assert payload["simulation_phases"]["production"]["duration"] == pytest.approx(0.3)
    stages = payload["simulation_phases"]["equilibration_stages"]
    assert [(stage["duration"], stage["ensemble"]) for stage in stages] == [
        (0.02, "NVT"),
        (0.02, "NPT"),
    ]
    assert stages[0]["position_restraints"] == [
        {"group": "protein_backbone", "force_constant": 4_184.0}
    ]
    assert stages[1]["position_restraints"] == []
    assert stages[1]["barostat_frequency"] == 25


def test_maximal_runtime_config_loads_without_fallback(tmp_path: Path) -> None:
    """The generated payload should pass the real config loader with all domains explicit."""
    protocol = MaximalMdProtocol()

    config = load_config(write_maximal_runtime_config(tmp_path, protocol))

    assert config.force_field.protein == "ff14sb_off_impropers_0.0.4.offxml"
    assert config.force_field.small_molecule == "openff-2.3.0.offxml"
    assert [attachment.moiety.force_field for attachment in config.conjugation.attachments] == [
        "openff-2.3.0.offxml",
        "glycam06",
        "glycam06",
    ]


def test_maximal_manifest_validator_checks_protocol_and_products(tmp_path: Path) -> None:
    """Summary validation should enforce composition, chemistry, export, and CUDA evidence."""
    protocol = MaximalMdProtocol()
    manifest = _manifest_payload(tmp_path, protocol)

    validate_maximal_manifest(manifest, protocol=protocol, require_cuda=True)
    manifest["products"]["link_atoms"] = ["NZ", "ND2"]
    with pytest.raises(AssertionError):
        validate_maximal_manifest(manifest, protocol=protocol, require_cuda=True)


def test_composition_evidence_audits_integer_rounding() -> None:
    """Composition summary should accept the nearest-integer 20 mol% count."""
    atoms = []
    for residue_number in range(1, 9):
        atoms.append(_atom("HOH", residue_number))
    for residue_number in range(9, 11):
        atoms.append(_atom("DMS", residue_number))
    atoms.append(_atom("NA", 11))
    payload = {
        "polymers": {"count": 3},
        "solvent": {
            "co_solvents": [{"mole_fraction": 0.2}],
            "ions": {"neutralize": True, "nacl_concentration": 0.0},
        },
    }

    evidence = composition_evidence(atoms, payload)

    assert evidence["water_count"] == 8
    assert evidence["dmso_count"] == 2
    assert evidence["ion_count"] == 1
    assert evidence["dmso_achieved_mole_fraction"] == pytest.approx(0.2)
    assert evidence["dmso_within_integer_rounding"] is True


def _atom(residue_name: str, residue_number: int) -> dict[str, object]:
    """Return one residue-identifying atom record."""
    return {
        "chain_id": "D",
        "residue_number": residue_number,
        "residue_name": residue_name,
    }


def _manifest_payload(tmp_path: Path, protocol: MaximalMdProtocol) -> dict[str, object]:
    """Return a compact passing maximal manifest for validator unit tests."""
    artifact = tmp_path / "artifact"
    artifact.touch()
    return {
        "effective_config": maximal_config_payload(tmp_path, protocol),
        "resolved_attachment_count": 3,
        "products": {
            "link_atoms": ["NZ", "ND2", "OG"],
            "residue_counts": {"LYX": 1, "ASX": 1, "OLS": 1},
        },
        "validation": {
            "linkage_distances_angstrom": [1.3, 1.4, 1.5],
            "bond_graph_status": "pass",
            "atom_presence_status": "pass",
            "charge_audit_status": "pass",
        },
        "charge_ownership": {"success": True},
        "composition": {
            "free_polymer_count_requested": 3,
            "water_count": 80,
            "dmso_count": 20,
            "ion_count": 1,
            "dmso_within_integer_rounding": True,
        },
        "artifacts": {"exports_exist": {"gro": True, "top": True}},
        "cuda": {
            "platform": "CUDA",
            "precision": "mixed",
            "production_steps": protocol.production_steps,
            "restrained_nvt_steps": protocol.restrained_nvt_steps,
            "unrestrained_npt_steps": protocol.unrestrained_npt_steps,
            "report_interval": protocol.report_interval,
            "checkpoint_interval": protocol.checkpoint_interval,
            "frame_count": protocol.expected_frames,
            "integrator": "LangevinMiddleIntegrator",
            "time_step_fs": 2.0,
            "temperature_k": 293.0,
            "pressure_atm": 1.0,
            "friction_per_ps": 1.0,
            "barostat_frequency": 25,
            "restrained_backbone_heavy_atom_count": 300,
            "artifacts_exist": {"trajectory": artifact.exists()},
            "minimized_potential_energy_kj_mol": -10.0,
            "start_potential_energy_kj_mol": -9.5,
            "final_potential_energy_kj_mol": -9.0,
        },
    }
