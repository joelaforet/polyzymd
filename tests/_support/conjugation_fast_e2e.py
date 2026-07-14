"""Opt-in fast mixed conjugation E2E helper."""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any

import pytest
import yaml

from polyzymd.builders.conjugation.api import build_conjugate_from_config
from polyzymd.builders.conjugation.pablo.parameterization import InterchangeParameterizationSettings
from polyzymd.builders.conjugation.relaxation import ConjugateRelaxationSettings
from polyzymd.builders.conjugation.system_workflow import ConjugatedPolymerSystemSettings
from polyzymd.config.loader import load_config
from polyzymd.exporters.gromacs import GromacsExporter

FIXTURE_DIR = Path(__file__).resolve().parents[1] / "data/conjugation/fast_mixed_e2e_1ubq"
FIXTURE_CONFIG = FIXTURE_DIR / "config.yaml"
SUMMARY_NAME = "fast_mixed_e2e_summary.json"


def run_fast_mixed_build_export(
    output_dir: Path,
    *,
    run_cuda: bool = False,
    cuda_steps: int = 100_000,
) -> dict[str, Any]:
    """Run the fast mixed conjugation checkpoint and return a compact summary.

    Parameters
    ----------
    output_dir : pathlib.Path
        Temporary output directory for generated artifacts.
    run_cuda : bool, optional
        Whether to run the separate CUDA dynamics smoke test, by default ``False``.
    cuda_steps : int, optional
        CUDA production-step count, by default 100000.

    Returns
    -------
    dict[str, Any]
        JSON-serializable validation summary.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    runtime_config = _runtime_config_path(output_dir)
    config = load_config(runtime_config)
    build_dir = output_dir / "build"
    settings = ConjugatedPolymerSystemSettings(
        create_final_interchange=True,
        force_regenerate_conjugate_polymer=True,
        conjugate_polymerist_max_retries=3,
        conjugate_polymerist_energy_minimize=True,
        conjugate_parameterization=InterchangeParameterizationSettings(
            small_molecule_force_field="openff-2.3.0.offxml",
        ),
        relaxation=ConjugateRelaxationSettings(platform_name="CPU"),
    )

    start = time.perf_counter()
    result = build_conjugate_from_config(
        config,
        output_dir=build_dir,
        settings=settings,
        free_polymer_seed=202610,
    )
    build_seconds = time.perf_counter() - start

    if result.final_interchange is None:
        raise AssertionError("Fast mixed E2E build did not create a final Interchange")

    export_dir = output_dir / "gromacs"
    exported = GromacsExporter(result.final_interchange, config).export(
        export_dir,
        prefix="fast_mixed_e2e_1ubq",
    )

    summary = _build_summary(
        result,
        exported=exported,
        build_seconds=build_seconds,
    )
    if run_cuda:
        summary["cuda"] = _run_cuda_smoke(result.final_interchange, output_dir / "cuda", cuda_steps)
    _write_summary(output_dir / SUMMARY_NAME, summary)
    return summary


def validate_fast_mixed_summary(summary: dict[str, Any], *, require_cuda: bool = False) -> None:
    """Assert the compact summary satisfies the fast mixed checkpoint contract.

    Parameters
    ----------
    summary : dict[str, Any]
        Summary returned by :func:`run_fast_mixed_build_export`.
    require_cuda : bool, optional
        Whether to assert the optional CUDA smoke-test fields, by default ``False``.
    """
    assert summary["final_interchange_created"] is True
    assert summary["resolved_attachment_count"] == 2
    assert summary["solvated_atom_count"] > summary["crosslinked_atom_count"] > 1_000
    assert summary["finite_coordinates"] is True
    assert summary["residue_counts"]["LYX"] > 0
    assert summary["residue_counts"]["ASX"] > 0
    assert summary["residue_counts"]["NAG"] > 0
    assert all(summary["export_exists"].values())

    validation = summary["validation"]
    assert validation["bond_graph_status"] == "pass"
    assert validation["atom_presence_status"] == "pass"
    assert validation["charge_audit_status"] == "pass"
    assert validation["relaxation_evidence_status"] == "pass"
    assert len(validation["linkage_distances_angstrom"]) == 2
    assert validation["linkage_distances_angstrom"][0] == pytest.approx(1.33, abs=0.10)
    assert validation["linkage_distances_angstrom"][1] == pytest.approx(1.45, abs=0.10)
    assert validation["close_contact_count"] <= 10

    charge_bridge = summary["charge_bridge"]
    assert charge_bridge["success"] is True
    assert charge_bridge["ff14sb_atom_count"] > 1_000
    assert charge_bridge["polymer_template_atom_count"] > 0
    assert charge_bridge["local_nagl_patch_atom_count"] > 0
    assert charge_bridge["total_charge_e"] == pytest.approx(charge_bridge["formal_charge_e"])

    relaxation = summary["relaxation"]
    assert math.isfinite(relaxation["stage_a_energy_before_min_kj_mol"])
    assert math.isfinite(relaxation["stage_a_energy_after_min_kj_mol"])
    assert len(relaxation["stage_a_linkage_distances_angstrom"]) == 2

    if require_cuda:
        cuda = summary["cuda"]
        assert cuda["platform"] == "CUDA"
        assert cuda["steps"] == 100_000
        assert cuda["checkpoint_exists"] is True
        assert math.isfinite(cuda["potential_energy_kj_mol"])


def _runtime_config_path(output_dir: Path) -> Path:
    """Write a runtime-local config so fixture directories stay read-only."""
    payload = yaml.safe_load(FIXTURE_CONFIG.read_text(encoding="utf-8"))
    runtime_dir = output_dir / "runtime_inputs"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    payload["enzyme"]["pdb_path"] = str(FIXTURE_DIR / "1ubq.pdb")
    payload["polymers"]["cache_directory"] = str(runtime_dir / "polymer_cache")
    payload["output"]["projects_directory"] = str(output_dir / "projects")
    payload["output"]["scratch_directory"] = str(output_dir / "scratch")
    runtime_config = runtime_dir / "config.yaml"
    runtime_config.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return runtime_config


def _build_summary(
    result: Any, *, exported: dict[str, Any], build_seconds: float
) -> dict[str, Any]:
    """Return a compact JSON-compatible build/export validation summary."""
    solvated_pdb = Path(result.solvated_pdb_path)
    crosslinked_pdb = Path(result.construction.crosslinked_pdb_path)
    validation = _load_json(getattr(result.construction, "validation_report_path", None))
    charge_bridge = _load_json(
        Path(result.construction.output_dir) / "product_state_charge_bridge.json"
    )
    relaxation = _load_json(Path(result.construction.output_dir) / "conjugate_relaxation.json")
    atoms = _pdb_atom_records(solvated_pdb)
    crosslinked_atoms = _pdb_atom_records(crosslinked_pdb)
    finite = all(
        math.isfinite(coord) for atom in atoms for coord in (atom["x"], atom["y"], atom["z"])
    )
    residue_counts = _residue_counts(crosslinked_atoms)
    attachment_specs = tuple(getattr(result, "attachment_specs", ()) or ())
    resolved = [spec.resolved_plan.model_dump(mode="json") for spec in attachment_specs]
    export_paths = _flatten_paths(exported)
    return {
        "fixture_config": str(FIXTURE_CONFIG),
        "build_seconds": build_seconds,
        "solvated_pdb": str(solvated_pdb),
        "crosslinked_pdb": str(crosslinked_pdb),
        "solvated_atom_count": len(atoms),
        "crosslinked_atom_count": len(crosslinked_atoms),
        "finite_coordinates": finite,
        "final_interchange_created": bool(result.final_interchange_created),
        "resolved_attachment_count": len(resolved),
        "resolved_plans": resolved,
        "residue_counts": residue_counts,
        "validation": _validation_summary(validation),
        "charge_bridge": _charge_bridge_summary(charge_bridge),
        "relaxation": _relaxation_summary(relaxation),
        "exports": {name: str(path) for name, path in export_paths.items()},
        "export_exists": {name: path.exists() for name, path in export_paths.items()},
    }


def _run_cuda_smoke(interchange: Any, output_dir: Path, steps: int) -> dict[str, Any]:
    """Run the optional CUDA smoke test without platform fallback."""
    import openmm
    from openmm import unit
    from openmm.app import CheckpointReporter, Simulation, StateDataReporter

    output_dir.mkdir(parents=True, exist_ok=True)
    platform = openmm.Platform.getPlatformByName("CUDA")
    system = interchange.to_openmm(combine_nonbonded_forces=True)
    topology = interchange.to_openmm_topology()
    positions = interchange.positions.to_openmm()
    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        2.0 * unit.femtosecond,
    )
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy(maxIterations=200)
    log_path = output_dir / "cuda_state.csv"
    checkpoint_path = output_dir / "cuda.chk"
    simulation.reporters.append(
        StateDataReporter(
            str(log_path),
            max(1, steps // 10),
            step=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            totalSteps=steps,
            separator=",",
        )
    )
    simulation.reporters.append(CheckpointReporter(str(checkpoint_path), max(1, steps // 2)))
    start = time.perf_counter()
    simulation.step(steps)
    state = simulation.context.getState(getEnergy=True)
    elapsed = time.perf_counter() - start
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(float(energy)):
        raise AssertionError("CUDA smoke test produced a non-finite potential energy")
    return {
        "platform": platform.getName(),
        "steps": steps,
        "elapsed_seconds": elapsed,
        "potential_energy_kj_mol": float(energy),
        "state_log": str(log_path),
        "checkpoint": str(checkpoint_path),
        "checkpoint_exists": checkpoint_path.exists(),
    }


def _validation_summary(payload: dict[str, Any]) -> dict[str, Any]:
    """Extract lightweight validation report fields for assertions."""
    if not payload:
        return {}
    linkage = payload.get("linkage_geometry", {}) or {}
    return {
        "status": payload.get("status"),
        "bond_graph_status": (payload.get("bond_graph", {}) or {}).get("status"),
        "atom_presence_status": (payload.get("atom_presence", {}) or {}).get("status"),
        "charge_audit_status": (payload.get("charge_audit", {}) or {}).get("status"),
        "parameter_coverage_status": (payload.get("parameter_coverage", {}) or {}).get("status"),
        "relaxation_evidence_status": (payload.get("relaxation_evidence", {}) or {}).get("status"),
        "linkage_geometry_status": linkage.get("status"),
        "linkage_distances_angstrom": linkage.get("linkage_distances_angstrom", []),
        "close_contact_count": linkage.get("close_contact_count"),
    }


def _charge_bridge_summary(payload: dict[str, Any]) -> dict[str, Any]:
    """Extract lightweight charge bridge fields for assertions."""
    if not payload:
        return {}
    return {
        "success": payload.get("success"),
        "ff14sb_atom_count": payload.get("ff14sb_atom_count"),
        "polymer_template_atom_count": payload.get("polymer_template_atom_count"),
        "local_nagl_patch_atom_count": payload.get("local_nagl_patch_atom_count"),
        "total_charge_e": payload.get("total_charge_e"),
        "formal_charge_e": payload.get("formal_charge_e"),
        "normalization_correction_e": payload.get("normalization_correction_e"),
        "max_per_atom_correction_e": payload.get("max_per_atom_correction_e"),
    }


def _relaxation_summary(payload: dict[str, Any]) -> dict[str, Any]:
    """Extract Stage A and energy fields from relaxation diagnostics."""
    if not payload:
        return {}
    keys = (
        "stage_a_energy_before_min_kj_mol",
        "stage_a_energy_after_min_kj_mol",
        "stage_a_linkage_distance_errors_angstrom",
        "stage_a_linkage_distances_angstrom",
        "stage_b_energy_after_md_kj_mol",
    )
    return {key: payload.get(key) for key in keys if key in payload}


def _pdb_atom_records(path: Path) -> list[dict[str, Any]]:
    """Parse ATOM and HETATM records needed by the E2E assertions."""
    records = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        records.append(
            {
                "atom_name": line[12:16].strip(),
                "residue_name": line[17:20].strip(),
                "chain_id": line[21].strip(),
                "residue_number": int(line[22:26]),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            }
        )
    return records


def _residue_counts(atoms: list[dict[str, Any]]) -> dict[str, int]:
    """Count atom records by residue name."""
    counts: dict[str, int] = {}
    for atom in atoms:
        residue_name = str(atom["residue_name"])
        counts[residue_name] = counts.get(residue_name, 0) + 1
    return counts


def _flatten_paths(exported: dict[str, Any]) -> dict[str, Path]:
    """Flatten exporter path payloads to the core runtime files."""
    paths: dict[str, Path] = {}
    for key in ("gro", "top", "em_mdp", "prod_mdp", "run_script"):
        if key in exported:
            paths[key] = Path(exported[key])
    for index, path in enumerate(exported.get("eq_mdps", []) or [], start=1):
        paths[f"eq_mdp_{index}"] = Path(path)
    return paths


def _load_json(path: Path | str | None) -> dict[str, Any]:
    """Load a JSON object if present, otherwise return an empty dict."""
    if path is None:
        return {}
    json_path = Path(path)
    if not json_path.exists():
        return {}
    return json.loads(json_path.read_text(encoding="utf-8"))


def _write_summary(path: Path, summary: dict[str, Any]) -> None:
    """Write the compact summary for local debugging."""
    path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
