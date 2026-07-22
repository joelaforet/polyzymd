"""Opt-in maximal mixed-conjugation acceptance helper."""

from __future__ import annotations

import importlib.metadata
import json
import math
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from polyzymd.builders.conjugation.api import build_conjugate_from_config
from polyzymd.builders.conjugation.pablo.parameterization import InterchangeParameterizationSettings
from polyzymd.builders.conjugation.relaxation import ConjugateRelaxationSettings
from polyzymd.builders.conjugation.system_workflow import ConjugatedPolymerSystemSettings
from polyzymd.config.loader import load_config
from polyzymd.exporters.gromacs import GromacsExporter
from tests._support.conjugation_fast_e2e import (
    FIXTURE_CONFIG as FAST_CONFIG,
)
from tests._support.conjugation_fast_e2e import (
    FIXTURE_DIR as FAST_FIXTURE_DIR,
)
from tests._support.conjugation_fast_e2e import (
    _charge_bridge_summary,
    _flatten_paths,
    _glycosylation_linkage_evidence,
    _load_json,
    _pdb_atom_records,
    _residue_counts,
    _validation_summary,
)

FIXTURE_DIR = Path(__file__).resolve().parents[1] / "data/conjugation/maximal_mixed_e2e_1ubq"
SUMMARY_NAME = "maximal_mixed_e2e_manifest.json"
INTEGRATOR_SEED = 202_611
VELOCITY_SEED = 202_612
FREE_POLYMER_SEED = 202_613


@dataclass(frozen=True)
class MaximalMdProtocol:
    """Approved maximal-acceptance protocol, not a universal production default."""

    restrained_nvt_steps: int = 10_000
    unrestrained_npt_steps: int = 10_000
    production_steps: int = 150_000
    report_interval: int = 1_000
    checkpoint_interval: int = 25_000
    time_step_fs: float = 2.0
    temperature_k: float = 293.0
    pressure_atm: float = 1.0
    friction_per_ps: float = 1.0
    barostat_frequency: int = 25
    restraint_k_kj_mol_nm2: float = 4_184.0

    def __post_init__(self) -> None:
        positive = (
            self.restrained_nvt_steps,
            self.unrestrained_npt_steps,
            self.production_steps,
            self.report_interval,
            self.checkpoint_interval,
            self.barostat_frequency,
        )
        if any(value <= 0 for value in positive):
            raise ValueError("All maximal protocol step counts and intervals must be positive")
        if self.production_steps % self.report_interval:
            raise ValueError("production_steps must be divisible by report_interval")
        if self.production_steps % self.checkpoint_interval:
            raise ValueError("production_steps must be divisible by checkpoint_interval")
        if (
            self.restrained_nvt_steps,
            self.unrestrained_npt_steps,
            self.production_steps,
            self.report_interval,
            self.checkpoint_interval,
            self.time_step_fs,
            self.temperature_k,
            self.pressure_atm,
            self.friction_per_ps,
            self.barostat_frequency,
        ) != (10_000, 10_000, 150_000, 1_000, 25_000, 2.0, 293.0, 1.0, 1.0, 25):
            raise ValueError("Maximal acceptance protocol values are fixed and approved")

    @property
    def expected_frames(self) -> int:
        """Return the exact expected reporter frame count."""
        return self.production_steps // self.report_interval


def maximal_config_payload(output_dir: Path, protocol: MaximalMdProtocol) -> dict[str, Any]:
    """Return the explicit effective configuration payload for the maximal case."""
    payload = yaml.safe_load(FAST_CONFIG.read_text(encoding="utf-8"))
    runtime_inputs = output_dir / "runtime_inputs"
    payload["name"] = "maximal_mixed_e2e_1ubq"
    payload["description"] = "Opt-in maximal simultaneous mixed-conjugation acceptance"
    payload["enzyme"]["pdb_path"] = str(FAST_FIXTURE_DIR / "1ubq.pdb")
    payload["polymers"]["count"] = 3
    payload["polymers"]["length"] = 3
    payload["polymers"]["cache_directory"] = str(runtime_inputs / "polymer_cache")
    attachments = payload["conjugation"]["attachments"]
    attachments[1]["name"] = "asn60_g42666ht_nglycan"
    attachments[1]["moiety"] = {
        "name": "glycam_G42666HT",
        "force_field": "glycam06",
        "input_path": str(FIXTURE_DIR / "glycam_G42666HT_CONECT.pdb"),
    }
    o_attachment = yaml.safe_load(yaml.safe_dump(attachments[1]))
    o_attachment["name"] = "ser20_g57321fi_oglycan"
    o_attachment["site"] = {
        "chain_id": "A",
        "residue_name": "SER",
        "residue_number": 20,
        "atom_name": "OG",
    }
    o_attachment["moiety"] = {
        "name": "glycam_G57321FI",
        "force_field": "glycam06",
        "input_path": str(FAST_FIXTURE_DIR / "glycam_G57321FI_CONECT.pdb"),
    }
    o_attachment["mechanism"] = {"name": "o_glycosylation"}
    attachments.append(o_attachment)
    payload["solvent"]["primary"] = {"type": "water", "model": "tip3p"}
    payload["solvent"]["co_solvents"] = [{"name": "dmso", "mole_fraction": 0.2}]
    payload["solvent"]["ions"] = {"neutralize": True, "nacl_concentration": 0.0}
    payload["force_field"] = {
        "protein": "ff14sb_off_impropers_0.0.4.offxml",
        "small_molecule": "openff-2.3.0.offxml",
    }
    payload["thermodynamics"] = {
        "temperature": protocol.temperature_k,
        "pressure": protocol.pressure_atm,
    }
    payload["simulation_phases"]["production"].update(
        {
            "duration": protocol.production_steps * protocol.time_step_fs / 1_000_000.0,
            "samples": protocol.expected_frames,
            "report_interval": protocol.report_interval,
            "time_step": protocol.time_step_fs,
            "thermostat": "LangevinMiddle",
            "thermostat_timescale": 1.0 / protocol.friction_per_ps,
            "barostat": "MC",
            "barostat_frequency": protocol.barostat_frequency,
            "checkpoint_interval": (
                protocol.checkpoint_interval * protocol.time_step_fs / 1_000_000.0
            ),
        }
    )
    payload["simulation_phases"]["equilibration_stages"] = [
        {
            "name": "maximal_restrained_backbone_nvt",
            "duration": protocol.restrained_nvt_steps * protocol.time_step_fs / 1_000_000.0,
            "samples": 20,
            "ensemble": "NVT",
            "temperature": protocol.temperature_k,
            "time_step": protocol.time_step_fs,
            "thermostat": "LangevinMiddle",
            "thermostat_timescale": 1.0 / protocol.friction_per_ps,
            "position_restraints": [
                {
                    "group": "protein_backbone",
                    "force_constant": protocol.restraint_k_kj_mol_nm2,
                }
            ],
        },
        {
            "name": "maximal_unrestrained_npt",
            "duration": protocol.unrestrained_npt_steps * protocol.time_step_fs / 1_000_000.0,
            "samples": 20,
            "ensemble": "NPT",
            "temperature": protocol.temperature_k,
            "time_step": protocol.time_step_fs,
            "thermostat": "LangevinMiddle",
            "thermostat_timescale": 1.0 / protocol.friction_per_ps,
            "barostat": "MC",
            "barostat_frequency": protocol.barostat_frequency,
            "position_restraints": [],
        },
    ]
    payload["output"].update(
        {
            "projects_directory": str(output_dir / "projects"),
            "scratch_directory": str(output_dir / "scratch"),
            "save_checkpoint": True,
            "save_state_data": True,
            "trajectory_format": "dcd",
        }
    )
    return payload


def write_maximal_runtime_config(output_dir: Path, protocol: MaximalMdProtocol) -> Path:
    """Write and return the runtime-local maximal configuration."""
    runtime_dir = output_dir / "runtime_inputs"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    path = runtime_dir / "config.yaml"
    path.write_text(
        yaml.safe_dump(maximal_config_payload(output_dir, protocol), sort_keys=False),
        encoding="utf-8",
    )
    return path


def run_maximal_mixed_acceptance(
    output_dir: Path, *, protocol: MaximalMdProtocol
) -> dict[str, Any]:
    """Build, exactly export, and run the explicit CUDA acceptance case."""
    _require_maximal_cuda_environment()
    output_dir.mkdir(parents=True, exist_ok=True)
    runtime_config = write_maximal_runtime_config(output_dir, protocol)
    config = load_config(runtime_config)
    settings = ConjugatedPolymerSystemSettings(
        create_final_interchange=True,
        force_regenerate_conjugate_polymer=True,
        conjugate_polymerist_max_retries=3,
        conjugate_polymerist_energy_minimize=True,
        conjugate_parameterization=InterchangeParameterizationSettings(
            small_molecule_force_field="openff-2.3.0.offxml"
        ),
        relaxation=ConjugateRelaxationSettings(platform_name="CPU"),
    )
    started = time.perf_counter()
    result = build_conjugate_from_config(
        config,
        output_dir=output_dir / "build",
        settings=settings,
        free_polymer_seed=FREE_POLYMER_SEED,
    )
    bundle = result.exact_export_bundle
    if bundle is None or not getattr(bundle, "is_exact_export_bundle", False):
        raise AssertionError("Maximal acceptance requires the exact native export bundle")
    exported = GromacsExporter(bundle, config).export(
        output_dir / "gromacs", prefix="maximal_mixed_e2e_1ubq"
    )
    manifest = build_maximal_manifest(
        result,
        config_payload=maximal_config_payload(output_dir, protocol),
        runtime_config=runtime_config,
        exported=exported,
        build_seconds=time.perf_counter() - started,
    )
    manifest["cuda"] = _run_explicit_cuda(bundle, output_dir / "cuda", protocol)
    _persist_and_validate_manifest(output_dir, manifest, protocol=protocol)
    return manifest


def _persist_and_validate_manifest(
    output_dir: Path, manifest: dict[str, Any], *, protocol: MaximalMdProtocol
) -> None:
    """Persist completed-run evidence before applying final acceptance assertions."""
    (output_dir / SUMMARY_NAME).write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    validate_maximal_manifest(manifest, protocol=protocol, require_cuda=True)


def build_maximal_manifest(
    result: Any,
    *,
    config_payload: dict[str, Any],
    runtime_config: Path,
    exported: dict[str, Any],
    build_seconds: float,
) -> dict[str, Any]:
    """Build compact machine-readable construction and ownership evidence."""
    crosslinked_atoms = _pdb_atom_records(Path(result.construction.crosslinked_pdb_path))
    solvated_atoms = _pdb_atom_records(Path(result.solvated_pdb_path))
    resolved = [item.model_dump(mode="json") for item in result.attachment_specs]
    validation = _load_json(result.construction.validation_report_path)
    charge_path = Path(result.construction.output_dir) / "product_state_charge_bridge.json"
    export_paths = _flatten_paths(exported)
    composition = composition_evidence(solvated_atoms, config_payload)
    return {
        "schema_version": 1,
        "source": _source_provenance(runtime_config),
        "build_seconds": build_seconds,
        "effective_config": config_payload,
        "resolved_attachment_count": len(resolved),
        "resolved_plans": resolved,
        "products": {
            "residue_counts": _residue_counts(crosslinked_atoms),
            "link_atoms": [plan["protein_link_atom"]["atom_name"] for plan in resolved],
            "leaving_groups": [
                {
                    "protein": [a["atom_name"] for a in plan.get("protein_leaving_atoms", [])],
                    "modifier_count": len(plan.get("modifier_leaving_atoms", [])),
                }
                for plan in resolved
            ],
            "n_glycosylation": _glycosylation_linkage_evidence(
                resolved,
                crosslinked_atoms,
                mechanism_name="n_glycosylation",
                product_site_residue_name="ASX",
                site_residue_number=60,
                site_atom_name="ND2",
            ),
            "o_glycosylation": _glycosylation_linkage_evidence(
                resolved,
                crosslinked_atoms,
                mechanism_name="o_glycosylation",
                product_site_residue_name="OLS",
                site_residue_number=20,
                site_atom_name="OG",
            ),
        },
        "composition": composition,
        "validation": _validation_summary(validation),
        "charge_ownership": _charge_bridge_summary(_load_json(charge_path), crosslinked_atoms),
        "artifacts": {
            "runtime_config": str(runtime_config),
            "crosslinked_pdb": str(result.construction.crosslinked_pdb_path),
            "solvated_pdb": str(result.solvated_pdb_path),
            "exports": {name: str(path) for name, path in export_paths.items()},
            "exports_exist": {name: path.exists() for name, path in export_paths.items()},
        },
    }


def composition_evidence(
    atoms: list[dict[str, Any]], config_payload: dict[str, Any]
) -> dict[str, Any]:
    """Summarize packed components without trusting reused PDB residue identities."""
    residues: dict[tuple[str, int, str], str] = {}
    for atom in atoms:
        key = (str(atom["chain_id"]), int(atom["residue_number"]), str(atom["residue_name"]))
        residues[key] = str(atom["residue_name"])
    counts: dict[str, int] = {}
    for residue_name in residues.values():
        counts[residue_name] = counts.get(residue_name, 0) + 1
    water_count = _packed_anchor_count(
        atoms,
        residue_names={"HOH", "WAT", "SOL"},
        anchor_atom_names={"O", "OW", "O1", "O1x"},
    )
    dmso_count = _packed_anchor_count(
        atoms,
        residue_names={"DMS", "DMSO"},
        anchor_atom_names={"S", "S1", "S1x"},
    )
    polymer_evidence = _free_polymer_multiplicity_evidence(atoms)
    neutral_solvent_count = water_count + dmso_count
    achieved = dmso_count / neutral_solvent_count if neutral_solvent_count else None
    target = float(config_payload["solvent"]["co_solvents"][0]["mole_fraction"])
    max_rounding_error = 0.5 / neutral_solvent_count if neutral_solvent_count else None
    ion_count = sum(1 for atom in atoms if str(atom["residue_name"]) in {"NA", "Na+", "CL", "Cl-"})
    return {
        "generic_residue_identity_counts": counts,
        "water_count": water_count,
        "dmso_count": dmso_count,
        "ion_count": ion_count,
        "free_polymer_count_requested": int(config_payload["polymers"]["count"]),
        **polymer_evidence,
        "dmso_target_mole_fraction": target,
        "dmso_achieved_mole_fraction": achieved,
        "dmso_integer_rounding_bound": max_rounding_error,
        "dmso_within_integer_rounding": achieved is not None
        and max_rounding_error is not None
        and abs(achieved - target) <= max_rounding_error + 1e-12,
        "neutralize": bool(config_payload["solvent"]["ions"]["neutralize"]),
        "nacl_concentration": float(config_payload["solvent"]["ions"]["nacl_concentration"]),
    }


def _packed_anchor_count(
    atoms: list[dict[str, Any]],
    *,
    residue_names: set[str],
    anchor_atom_names: set[str],
) -> int:
    """Count packed molecules by one chemically unique atom-name anchor."""
    return sum(
        1
        for atom in atoms
        if str(atom["residue_name"]) in residue_names
        and str(atom["atom_name"]) in anchor_atom_names
    )


def _free_polymer_multiplicity_evidence(atoms: list[dict[str, Any]]) -> dict[str, Any]:
    """Infer free SBM count from repeated per-molecule unique atom names."""
    atom_name_counts: dict[str, int] = {}
    for atom in atoms:
        if str(atom["residue_name"]) == "SBM":
            atom_name = str(atom["atom_name"])
            atom_name_counts[atom_name] = atom_name_counts.get(atom_name, 0) + 1
    if not atom_name_counts:
        return {
            "free_polymer_count_achieved": 0,
            "free_polymer_anchor_atom_name": None,
            "free_polymer_unique_atom_name_count": 0,
            "free_polymer_atom_name_multiplicity_consistent": False,
        }
    multiplicities = set(atom_name_counts.values())
    anchor = min(atom_name_counts)
    return {
        "free_polymer_count_achieved": (
            atom_name_counts[anchor] if len(multiplicities) == 1 else None
        ),
        "free_polymer_anchor_atom_name": anchor,
        "free_polymer_unique_atom_name_count": len(atom_name_counts),
        "free_polymer_atom_name_multiplicity_consistent": len(multiplicities) == 1,
        "free_polymer_atom_name_multiplicity_min": min(multiplicities),
        "free_polymer_atom_name_multiplicity_max": max(multiplicities),
    }


def validate_maximal_manifest(
    manifest: dict[str, Any], *, protocol: MaximalMdProtocol, require_cuda: bool
) -> None:
    """Validate the maximal acceptance manifest without silent weakening."""
    config = manifest["effective_config"]
    assert config["force_field"]["protein"] == "ff14sb_off_impropers_0.0.4.offxml"
    assert config["force_field"]["small_molecule"] == "openff-2.3.0.offxml"
    assert [a["moiety"]["force_field"] for a in config["conjugation"]["attachments"]] == [
        "openff-2.3.0.offxml",
        "glycam06",
        "glycam06",
    ]
    assert config["polymers"]["count"] == 3
    assert config["polymers"]["length"] == 3
    assert config["solvent"]["primary"] == {"type": "water", "model": "tip3p"}
    assert config["solvent"]["co_solvents"] == [{"name": "dmso", "mole_fraction": 0.2}]
    assert config["solvent"]["ions"]["neutralize"] is True
    assert config["solvent"]["ions"]["nacl_concentration"] == 0.0
    assert manifest["resolved_attachment_count"] == 3
    assert set(manifest["products"]["link_atoms"]) == {"NZ", "ND2", "OG"}
    residue_counts = manifest["products"]["residue_counts"]
    assert all(residue_counts[name] > 0 for name in ("LYX", "ASX", "OLS"))
    assert len(manifest["validation"]["linkage_distances_angstrom"]) == 3
    assert manifest["validation"]["bond_graph_status"] == "pass"
    assert manifest["validation"]["atom_presence_status"] == "pass"
    assert manifest["validation"]["charge_audit_status"] == "pass"
    assert manifest["charge_ownership"]["success"] is True
    composition = manifest["composition"]
    assert composition["free_polymer_count_requested"] == 3
    assert composition["free_polymer_atom_name_multiplicity_consistent"] is True
    assert composition["free_polymer_count_achieved"] == composition["free_polymer_count_requested"]
    assert composition["water_count"] > 0 and composition["dmso_count"] > 0
    assert composition["ion_count"] > 0
    assert composition["dmso_within_integer_rounding"] is True
    assert all(manifest["artifacts"]["exports_exist"].values())
    if require_cuda:
        cuda = manifest["cuda"]
        assert cuda["platform"] == "CUDA"
        assert cuda["precision"] == "mixed"
        assert cuda["production_steps"] == protocol.production_steps
        assert cuda["restrained_nvt_steps"] == protocol.restrained_nvt_steps
        assert cuda["unrestrained_npt_steps"] == protocol.unrestrained_npt_steps
        assert cuda["report_interval"] == protocol.report_interval
        assert cuda["checkpoint_interval"] == protocol.checkpoint_interval
        assert cuda["frame_count"] == protocol.expected_frames
        assert cuda["integrator"] == "LangevinMiddleIntegrator"
        assert cuda["time_step_fs"] == 2.0
        assert cuda["temperature_k"] == 293.0
        assert cuda["pressure_atm"] == 1.0
        assert cuda["friction_per_ps"] == 1.0
        assert cuda["barostat_frequency"] == 25
        assert cuda["restrained_backbone_heavy_atom_count"] > 0
        assert all(cuda["artifacts_exist"].values())
        assert math.isfinite(cuda["minimized_potential_energy_kj_mol"])
        assert math.isfinite(cuda["start_potential_energy_kj_mol"])
        assert math.isfinite(cuda["final_potential_energy_kj_mol"])


def _require_maximal_cuda_environment() -> None:
    """Fail before construction unless the PTX-compatible pixi environment is active."""
    pixi_environment = os.environ.get("PIXI_ENVIRONMENT_NAME")
    conda_prefix_name = Path(os.environ.get("CONDA_PREFIX", "")).name
    if "build-cuda-12-6" not in {pixi_environment, conda_prefix_name}:
        raise RuntimeError(
            "Maximal CUDA acceptance must run via "
            "'/home/joelaforet/.pixi/bin/pixi run -e build-cuda-12-6'. "
            "The ordinary build environment reproduced "
            "CUDA_ERROR_UNSUPPORTED_PTX_VERSION."
        )


def _run_explicit_cuda(
    bundle: Any, output_dir: Path, protocol: MaximalMdProtocol
) -> dict[str, Any]:
    """Run the explicit CUDA protocol and verify every expected DCD frame."""
    import MDAnalysis as mda
    import openmm
    from openmm import unit
    from openmm.app import (
        CheckpointReporter,
        DCDReporter,
        PDBFile,
        Simulation,
        StateDataReporter,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    platform = openmm.Platform.getPlatformByName("CUDA")
    properties = {"Precision": "mixed"}
    system = bundle.to_openmm()
    if any(isinstance(force, openmm.MonteCarloBarostat) for force in system.getForces()):
        raise AssertionError("Exact bundle unexpectedly contains a pre-existing barostat")
    integrator = openmm.LangevinMiddleIntegrator(
        protocol.temperature_k * unit.kelvin,
        protocol.friction_per_ps / unit.picosecond,
        protocol.time_step_fs * unit.femtosecond,
    )
    integrator.setRandomNumberSeed(INTEGRATOR_SEED)
    simulation = Simulation(bundle.to_openmm_topology(), system, integrator, platform, properties)
    simulation.context.setPositions(bundle.to_openmm_positions())
    simulation.minimizeEnergy(maxIterations=0)
    minimized_state = simulation.context.getState(getPositions=True, getEnergy=True)
    minimized_pdb = output_dir / "minimized.pdb"
    with minimized_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(
            bundle.to_openmm_topology(), minimized_state.getPositions(), handle, keepIds=True
        )
    restraint = openmm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    restraint.addGlobalParameter(
        "k", protocol.restraint_k_kj_mol_nm2 * unit.kilojoule_per_mole / unit.nanometer**2
    )
    for name in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(name)
    minimized_positions = minimized_state.getPositions(asNumpy=True)
    restrained_indices = []
    for atom in bundle.to_openmm_topology().atoms():
        if atom.residue.chain.id == "A" and atom.name in {"N", "CA", "C", "O"}:
            position = minimized_positions[atom.index].value_in_unit(unit.nanometer)
            restraint.addParticle(atom.index, [position[0], position[1], position[2]])
            restrained_indices.append(atom.index)
    if not restrained_indices:
        raise AssertionError("No protein backbone heavy atoms selected for NVT restraints")
    restraint_index = system.addForce(restraint)
    simulation.context.reinitialize(preserveState=True)
    simulation.context.setVelocitiesToTemperature(
        protocol.temperature_k * unit.kelvin, VELOCITY_SEED
    )
    simulation.step(protocol.restrained_nvt_steps)
    system.removeForce(restraint_index)
    barostat = openmm.MonteCarloBarostat(
        protocol.pressure_atm * unit.atmosphere,
        protocol.temperature_k * unit.kelvin,
        protocol.barostat_frequency,
    )
    barostat.setRandomNumberSeed(INTEGRATOR_SEED + 1)
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)
    simulation.step(protocol.unrestrained_npt_steps)
    start_state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    start_pdb = output_dir / "start.pdb"
    with start_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(
            bundle.to_openmm_topology(), start_state.getPositions(), handle, keepIds=True
        )
    dcd = output_dir / "production.dcd"
    state_csv = output_dir / "state.csv"
    checkpoint = output_dir / "production.chk"
    simulation.reporters.extend(
        [
            DCDReporter(str(dcd), protocol.report_interval),
            StateDataReporter(
                str(state_csv),
                protocol.report_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                volume=True,
            ),
            CheckpointReporter(str(checkpoint), protocol.checkpoint_interval),
        ]
    )
    simulation.step(protocol.production_steps)
    final_state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    final_pdb = output_dir / "final.pdb"
    with final_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(
            bundle.to_openmm_topology(), final_state.getPositions(), handle, keepIds=True
        )
    system_xml = output_dir / "system.xml"
    state_xml = output_dir / "final_state.xml"
    system_xml.write_text(openmm.XmlSerializer.serialize(system), encoding="utf-8")
    state_xml.write_text(openmm.XmlSerializer.serialize(final_state), encoding="utf-8")
    frame_count = len(mda.Universe(str(minimized_pdb), str(dcd)).trajectory)
    if frame_count != protocol.expected_frames:
        raise AssertionError(
            f"Expected exactly {protocol.expected_frames} DCD frames, got {frame_count}"
        )
    artifacts = {
        "trajectory": dcd,
        "state_data": state_csv,
        "checkpoint": checkpoint,
        "minimized_pdb": minimized_pdb,
        "start_pdb": start_pdb,
        "final_pdb": final_pdb,
        "system_xml": system_xml,
        "state_xml": state_xml,
    }
    energy_unit = unit.kilojoule_per_mole
    return {
        "platform": platform.getName(),
        "device": platform.getPropertyValue(simulation.context, "DeviceIndex"),
        "precision": platform.getPropertyValue(simulation.context, "Precision"),
        "integrator": type(integrator).__name__,
        "temperature_k": protocol.temperature_k,
        "pressure_atm": protocol.pressure_atm,
        "time_step_fs": protocol.time_step_fs,
        "production_steps": protocol.production_steps,
        "restrained_nvt_steps": protocol.restrained_nvt_steps,
        "unrestrained_npt_steps": protocol.unrestrained_npt_steps,
        "report_interval": protocol.report_interval,
        "checkpoint_interval": protocol.checkpoint_interval,
        "friction_per_ps": protocol.friction_per_ps,
        "barostat_frequency": protocol.barostat_frequency,
        "restraint_k_kj_mol_nm2": protocol.restraint_k_kj_mol_nm2,
        "restrained_backbone_heavy_atom_count": len(restrained_indices),
        "frame_count": frame_count,
        "integrator_seed": INTEGRATOR_SEED,
        "velocity_seed": VELOCITY_SEED,
        "minimized_potential_energy_kj_mol": float(
            minimized_state.getPotentialEnergy().value_in_unit(energy_unit)
        ),
        "start_potential_energy_kj_mol": float(
            start_state.getPotentialEnergy().value_in_unit(energy_unit)
        ),
        "final_potential_energy_kj_mol": float(
            final_state.getPotentialEnergy().value_in_unit(energy_unit)
        ),
        "artifacts": {name: str(path) for name, path in artifacts.items()},
        "artifacts_exist": {name: path.exists() for name, path in artifacts.items()},
    }


def _source_provenance(runtime_config: Path) -> dict[str, Any]:
    """Return repository, import, and package provenance for the manifest."""
    import polyzymd
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    root = Path(__file__).resolve().parents[2]

    def git(*args: str) -> str:
        return subprocess.check_output(["git", *args], cwd=root, text=True).strip()

    packages = {}
    for name in ("polyzymd", "openmm", "openff-toolkit", "MDAnalysis"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "worktree": str(root),
        "branch": git("branch", "--show-current"),
        "head": git("rev-parse", "HEAD"),
        "dirty_status": git("status", "--short"),
        "runtime_config": str(runtime_config),
        "imports": {
            "polyzymd": str(polyzymd.__file__),
            "system_workflow": str(workflow_module.__file__),
        },
        "packages": packages,
        "environment": {"CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES")},
        "seeds": {
            "free_polymer": FREE_POLYMER_SEED,
            "integrator": INTEGRATOR_SEED,
            "velocities": VELOCITY_SEED,
        },
    }
