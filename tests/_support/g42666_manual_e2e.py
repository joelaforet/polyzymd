"""Manual G42666 N-glycosylation build and production-MD acceptance runner."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from pathlib import Path
from typing import Any

import yaml

from polyzymd.builders.conjugation.api import build_conjugate_from_config
from polyzymd.builders.conjugation.pablo.parameterization import InterchangeParameterizationSettings
from polyzymd.builders.conjugation.relaxation import ConjugateRelaxationSettings
from polyzymd.builders.conjugation.system_workflow import ConjugatedPolymerSystemSettings
from polyzymd.config.loader import load_config

FIXTURE_DIR = Path(__file__).resolve().parents[1] / "data/conjugation/fast_mixed_e2e_1ubq"
FIXTURE_CONFIG = FIXTURE_DIR / "config.yaml"
G42666_ENV_VAR = "POLYZYMD_G42666_CONECT_PDB"


def run_g42666_manual_e2e(
    output_dir: Path,
    *,
    glycan_path: Path,
    production_steps: int = 150_000,
    report_interval: int = 1_000,
    platform_name: str = "CPU",
) -> dict[str, Any]:
    """Run the guarded local G42666 build and OpenMM production acceptance check.

    Parameters
    ----------
    output_dir : pathlib.Path
        Directory for local build and simulation artifacts.
    glycan_path : pathlib.Path, optional
        Strict-CONECT G42666 PDB fragment, by default the local acceptance path.
    production_steps : int, optional
        OpenMM production steps to run, by default 150000.
    report_interval : int, optional
        Reporter interval in steps, by default 1000.
    platform_name : str, optional
        OpenMM platform name, by default ``"CPU"``.

    Returns
    -------
    dict[str, Any]
        JSON-compatible diagnostics for acceptance review.
    """
    if not glycan_path.exists():
        raise FileNotFoundError(f"G42666 CONECT glycan fixture not found: {glycan_path}")
    output_dir.mkdir(parents=True, exist_ok=True)
    runtime_config = _runtime_config_path(output_dir, glycan_path=glycan_path)
    config = load_config(runtime_config)
    settings = ConjugatedPolymerSystemSettings(
        create_final_interchange=True,
        pdb_fragment_output_mode="experimental_pablo",
        force_regenerate_conjugate_polymer=True,
        conjugate_polymerist_max_retries=3,
        conjugate_polymerist_energy_minimize=True,
        conjugate_parameterization=InterchangeParameterizationSettings(
            small_molecule_force_field="openff-2.3.0.offxml",
        ),
        relaxation=ConjugateRelaxationSettings(platform_name=platform_name),
    )
    start = time.perf_counter()
    result = build_conjugate_from_config(
        config,
        output_dir=output_dir / "build",
        settings=settings,
        free_polymer_seed=202610,
    )
    if result.final_interchange is None:
        raise AssertionError("G42666 E2E build did not retain a final Interchange")
    simulation_summary = _run_openmm_production(
        result.final_interchange,
        output_dir / "production_md",
        production_steps=production_steps,
        report_interval=report_interval,
        platform_name=platform_name,
    )
    summary = {
        "runtime_config": str(runtime_config),
        "glycan_path": str(glycan_path),
        "build_seconds": time.perf_counter() - start,
        "crosslinked_pdb": str(result.crosslinked_conjugate_pdb_path),
        "relaxed_conjugate_pdb": str(result.relaxed_conjugate_pdb_path),
        "solvated_pdb": str(result.solvated_pdb_path),
        "final_interchange_created": bool(result.final_interchange_created),
        "attachment_diagnostics": _attachment_diagnostics(result.attachment_specs),
        "simulation": simulation_summary,
    }
    _validate_summary(summary)
    summary_path = output_dir / "g42666_manual_e2e_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def _runtime_config_path(output_dir: Path, *, glycan_path: Path) -> Path:
    """Write a runtime config that swaps the glycan source to G42666 input_path."""
    payload = yaml.safe_load(FIXTURE_CONFIG.read_text(encoding="utf-8"))
    runtime_dir = output_dir / "runtime_inputs"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    payload["enzyme"]["pdb_path"] = str(FIXTURE_DIR / "1ubq.pdb")
    payload["polymers"]["cache_directory"] = str(runtime_dir / "polymer_cache")
    payload["polymers"]["enabled"] = False
    payload["conjugation"]["attachments"] = [payload["conjugation"]["attachments"][1]]
    glycan = payload["conjugation"]["attachments"][0]["moiety"]
    glycan.pop("smiles", None)
    glycan.pop("residue_name", None)
    glycan["input_path"] = str(glycan_path)
    payload["simulation_phases"]["production"].update(
        {"duration": 300.0, "samples": 150, "report_interval": 1000, "time_step": 2.0}
    )
    payload["output"]["projects_directory"] = str(output_dir / "projects")
    payload["output"]["scratch_directory"] = str(output_dir / "scratch")
    runtime_config = runtime_dir / "config.yaml"
    runtime_config.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return runtime_config


def _run_openmm_production(
    interchange: Any,
    output_dir: Path,
    *,
    production_steps: int,
    report_interval: int,
    platform_name: str,
) -> dict[str, Any]:
    """Run production MD and return finite-energy and trajectory diagnostics."""
    import openmm
    from openmm import unit
    from openmm.app import DCDReporter, Simulation, StateDataReporter

    output_dir.mkdir(parents=True, exist_ok=True)
    platform = openmm.Platform.getPlatformByName(platform_name)
    system = interchange.to_openmm(combine_nonbonded_forces=True)
    topology = interchange.to_openmm_topology()
    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        2.0 * unit.femtosecond,
    )
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(interchange.positions.to_openmm())
    initial_energy = _potential_energy(simulation)
    trajectory_path = output_dir / "g42666_production.dcd"
    state_path = output_dir / "g42666_state.csv"
    simulation.reporters.append(DCDReporter(str(trajectory_path), report_interval))
    simulation.reporters.append(
        StateDataReporter(
            str(state_path),
            report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
        )
    )
    simulation.step(production_steps)
    final_energy = _potential_energy(simulation)
    state_rows = _read_state_rows(state_path)
    frame_count = _trajectory_frame_count(trajectory_path, topology)
    return {
        "platform": platform.getName(),
        "steps": production_steps,
        "report_interval": report_interval,
        "requested_frames": production_steps // report_interval,
        "trajectory": str(trajectory_path),
        "state_csv": str(state_path),
        "state_rows": len(state_rows),
        "trajectory_frames": frame_count,
        "initial_potential_energy_kj_mol": initial_energy,
        "final_potential_energy_kj_mol": final_energy,
        "state_energies_finite": all(math.isfinite(row["potential_energy"]) for row in state_rows),
    }


def _potential_energy(simulation: Any) -> float:
    """Return the current potential energy in kJ/mol."""
    from openmm import unit

    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    return float(energy.value_in_unit(unit.kilojoule_per_mole))


def _read_state_rows(path: Path) -> list[dict[str, float]]:
    """Read OpenMM state data and normalize potential-energy values."""
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(_state_csv_lines(handle))
        for row in reader:
            energy_key = _potential_energy_key(row)
            rows.append({"potential_energy": float(row[energy_key])})
    return rows


def _state_csv_lines(lines: Any) -> Any:
    """Yield parseable OpenMM state CSV lines while preserving the commented header."""
    header_seen = False
    for line in lines:
        if not line.startswith("#"):
            yield line
            continue
        uncommented = line[1:]
        if _is_state_csv_header(uncommented):
            if header_seen:
                raise ValueError("OpenMM state CSV contains multiple StateDataReporter headers")
            header_seen = True
            yield uncommented


def _is_state_csv_header(line: str) -> bool:
    """Return whether an uncommented line is the OpenMM StateDataReporter header."""
    try:
        fields = next(csv.reader([line]))
    except csv.Error:
        return False
    normalized = [field.strip().strip('"') for field in fields]
    return bool(
        len(normalized) >= 2
        and normalized[0] == "Step"
        and any(field.startswith("Potential Energy") for field in normalized[1:])
    )


def _potential_energy_key(row: dict[str, str]) -> str:
    """Return the CSV column key that contains potential energy values."""
    try:
        return next(key for key in row if key is not None and "Potential Energy" in key)
    except StopIteration as exc:
        raise ValueError("OpenMM state CSV is missing potential-energy data") from exc


def _trajectory_frame_count(path: Path, topology: Any) -> int | None:
    """Return saved trajectory frame count when MDAnalysis is available."""
    try:
        import MDAnalysis as mda
    except ImportError:
        return None
    universe = mda.Universe(topology, str(path))
    return len(universe.trajectory)


def _attachment_diagnostics(specs: tuple[Any, ...]) -> list[dict[str, Any]]:
    """Return compact resolved-plan diagnostics for attachment assertions."""
    diagnostics = []
    for spec in specs:
        plan = spec.resolved_plan
        diagnostics.append(
            {
                "mechanism_name": plan.contract.mechanism_name,
                "protein_link_atom": plan.protein_link_atom.atom_name,
                "protein_leaving_atoms": [atom.atom_name for atom in plan.protein_leaving_atoms],
                "modifier_link_serial": plan.modifier_link_atom.serial,
                "modifier_leaving_serials": [atom.serial for atom in plan.modifier_leaving_atoms],
                "crosslink": list(plan.pablo_crosslink_requirement.linking_atoms),
            }
        )
    return diagnostics


def _validate_summary(summary: dict[str, Any]) -> None:
    """Validate the manual E2E acceptance summary before writing it."""
    n_glycan = next(
        item
        for item in summary["attachment_diagnostics"]
        if item["mechanism_name"] == "n_glycosylation"
    )
    assert n_glycan["protein_link_atom"] == "ND2"
    assert n_glycan["protein_leaving_atoms"] == ["HD21"]
    assert n_glycan["modifier_link_serial"] == 4
    assert n_glycan["modifier_leaving_serials"] == [1, 2]
    simulation = summary["simulation"]
    assert simulation["state_rows"] == simulation["requested_frames"]
    if simulation["trajectory_frames"] is not None:
        assert simulation["trajectory_frames"] == simulation["requested_frames"]
    assert math.isfinite(simulation["initial_potential_energy_kj_mol"])
    assert math.isfinite(simulation["final_potential_energy_kj_mol"])
    assert simulation["state_energies_finite"] is True


def main() -> None:
    """Run the manual E2E from the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("/tmp/opencode/g42666_manual_e2e"))
    parser.add_argument("--glycan-path", type=Path, default=None)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--production-steps", type=int, default=150_000)
    parser.add_argument("--report-interval", type=int, default=1_000)
    args = parser.parse_args()
    glycan_path = args.glycan_path or _env_glycan_path()
    if glycan_path is None:
        raise SystemExit(f"Provide --glycan-path or set {G42666_ENV_VAR}")
    summary = run_g42666_manual_e2e(
        args.output_dir,
        glycan_path=glycan_path,
        production_steps=args.production_steps,
        report_interval=args.report_interval,
        platform_name=args.platform,
    )
    print(json.dumps(summary, indent=2, allow_nan=False))


def _env_glycan_path() -> Path | None:
    """Return the optional G42666 fixture path from the environment."""
    value = os.environ.get(G42666_ENV_VAR)
    return Path(value) if value else None


if __name__ == "__main__":
    main()
