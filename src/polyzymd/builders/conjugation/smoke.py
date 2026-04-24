"""Restrained vacuum OpenMM smoke runner for conjugation construction."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from pydantic import BaseModel, Field


class VacuumSmokeSettings(BaseModel):
    """Settings for the short restrained vacuum smoke run."""

    minimization_max_iterations: int = Field(100, ge=0)
    minimization_tolerance_kj_mol_nm: float = Field(10.0, gt=0)
    nvt_steps: int = Field(10, ge=0)
    temperature_kelvin: float = Field(310.0, gt=0)
    timestep_femtoseconds: float = Field(1.0, gt=0)
    friction_per_picosecond: float = Field(1.0, gt=0)
    protein_heavy_restraint_k_kj_mol_nm2: float = Field(5000.0, gt=0)
    platform_name: str | None = None
    smoke_json_name: str = "vacuum_smoke.json"
    minimized_pdb_name: str = "assembled_minimized.pdb"
    equilibrated_pdb_name: str = "assembled_equilibrated.pdb"


class VacuumSmokeResult(BaseModel):
    """Summary and artifacts from a restrained vacuum smoke run."""

    success: bool
    output_dir: Path
    smoke_json_path: Path
    minimized_pdb_path: Path | None = None
    equilibrated_pdb_path: Path | None = None
    platform_name: str
    restrained_atom_count: int
    energy_before_min_kj_mol: float
    energy_after_min_kj_mol: float
    energy_before_nvt_kj_mol: float
    energy_after_nvt_kj_mol: float
    nvt_steps: int
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)


def run_restrained_vacuum_smoke(
    interchange: Any,
    output_dir: Path | str,
    *,
    protein_heavy_atom_indices: tuple[int, ...] | None = None,
    settings: VacuumSmokeSettings | None = None,
) -> VacuumSmokeResult:
    """Run restrained minimization and very short vacuum NVT MD.

    Parameters
    ----------
    interchange : Any
        OpenFF Interchange-like object exposing OpenMM conversion methods.
    output_dir : pathlib.Path or str
        Directory for smoke artifacts.
    protein_heavy_atom_indices : tuple of int or None, optional
        Atom indices to restrain. When omitted, chain-A non-hydrogen atoms are
        inferred from the OpenMM topology where possible.
    settings : VacuumSmokeSettings or None, optional
        Smoke settings, by default ``None``.

    Returns
    -------
    VacuumSmokeResult
        Energies and artifact paths for the smoke run.

    Raises
    ------
    RuntimeError
        If OpenMM is unavailable, energies/positions are non-finite, or the
        simulation stack raises an exception.
    """
    try:
        import openmm
        from openmm import app as openmm_app
        from openmm import unit as openmm_unit
    except ImportError as exc:
        raise RuntimeError(
            "OpenMM is required for the conjugation vacuum smoke. Run under a suitable pixi "
            "environment such as cuda-12-4 when GPU resources are allocated."
        ) from exc

    smoke_settings = settings or VacuumSmokeSettings()
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    topology = interchange.to_openmm_topology()
    initial_positions = _openmm_positions_from_interchange(interchange, openmm_unit)
    restrained_indices = protein_heavy_atom_indices or _infer_chain_a_heavy_indices(topology)
    if not restrained_indices:
        raise RuntimeError("No protein heavy atoms were available for restrained vacuum smoke")

    platform = _select_platform(openmm, smoke_settings.platform_name)
    platform_name = platform.getName() if hasattr(platform, "getName") else str(platform)

    system_min = interchange.to_openmm_system()
    _add_positional_restraints(
        system_min,
        initial_positions,
        restrained_indices,
        smoke_settings.protein_heavy_restraint_k_kj_mol_nm2,
        openmm,
        openmm_unit,
    )
    integrator_min = openmm.VerletIntegrator(0.001 * openmm_unit.picoseconds)
    simulation_min = openmm_app.Simulation(topology, system_min, integrator_min, platform)
    simulation_min.context.setPositions(initial_positions)

    energy_before_min = _state_energy_kj_mol(
        simulation_min.context.getState(getEnergy=True), openmm_unit
    )
    validate_finite_energy(energy_before_min, label="energy_before_min_kj_mol")
    openmm.LocalEnergyMinimizer.minimize(
        simulation_min.context,
        tolerance=smoke_settings.minimization_tolerance_kj_mol_nm
        * openmm_unit.kilojoule_per_mole
        / openmm_unit.nanometer,
        maxIterations=smoke_settings.minimization_max_iterations,
    )
    state_after_min = simulation_min.context.getState(getEnergy=True, getPositions=True)
    energy_after_min = _state_energy_kj_mol(state_after_min, openmm_unit)
    minimized_positions = state_after_min.getPositions(asNumpy=True)
    validate_finite_energy(energy_after_min, label="energy_after_min_kj_mol")
    validate_finite_positions(minimized_positions, openmm_unit, label="minimized_positions")

    minimized_pdb = artifact_dir / smoke_settings.minimized_pdb_name
    _write_openmm_pdb(openmm_app, topology, minimized_positions, minimized_pdb)

    system_nvt = interchange.to_openmm_system()
    _add_positional_restraints(
        system_nvt,
        minimized_positions,
        restrained_indices,
        smoke_settings.protein_heavy_restraint_k_kj_mol_nm2,
        openmm,
        openmm_unit,
    )
    integrator_nvt = openmm.LangevinMiddleIntegrator(
        smoke_settings.temperature_kelvin * openmm_unit.kelvin,
        smoke_settings.friction_per_picosecond / openmm_unit.picosecond,
        smoke_settings.timestep_femtoseconds * openmm_unit.femtosecond,
    )
    simulation_nvt = openmm_app.Simulation(topology, system_nvt, integrator_nvt, platform)
    simulation_nvt.context.setPositions(minimized_positions)
    simulation_nvt.context.setVelocitiesToTemperature(
        smoke_settings.temperature_kelvin * openmm_unit.kelvin
    )

    energy_before_nvt = _state_energy_kj_mol(
        simulation_nvt.context.getState(getEnergy=True), openmm_unit
    )
    validate_finite_energy(energy_before_nvt, label="energy_before_nvt_kj_mol")
    if smoke_settings.nvt_steps:
        simulation_nvt.step(smoke_settings.nvt_steps)
    state_after_nvt = simulation_nvt.context.getState(getEnergy=True, getPositions=True)
    energy_after_nvt = _state_energy_kj_mol(state_after_nvt, openmm_unit)
    equilibrated_positions = state_after_nvt.getPositions(asNumpy=True)
    validate_finite_energy(energy_after_nvt, label="energy_after_nvt_kj_mol")
    validate_finite_positions(equilibrated_positions, openmm_unit, label="equilibrated_positions")

    equilibrated_pdb = artifact_dir / smoke_settings.equilibrated_pdb_name
    _write_openmm_pdb(openmm_app, topology, equilibrated_positions, equilibrated_pdb)

    result = VacuumSmokeResult(
        success=True,
        output_dir=artifact_dir,
        smoke_json_path=artifact_dir / smoke_settings.smoke_json_name,
        minimized_pdb_path=minimized_pdb,
        equilibrated_pdb_path=equilibrated_pdb,
        platform_name=platform_name,
        restrained_atom_count=len(restrained_indices),
        energy_before_min_kj_mol=float(energy_before_min),
        energy_after_min_kj_mol=float(energy_after_min),
        energy_before_nvt_kj_mol=float(energy_before_nvt),
        energy_after_nvt_kj_mol=float(energy_after_nvt),
        nvt_steps=smoke_settings.nvt_steps,
        diagnostics=("Restrained vacuum OpenMM smoke completed with finite state values",),
    )
    result.smoke_json_path.write_text(json.dumps(result.model_dump(mode="json"), indent=2) + "\n")
    return result


def validate_finite_energy(value: float, *, label: str = "energy") -> None:
    """Validate that an energy-like scalar is finite."""
    if not math.isfinite(float(value)):
        raise RuntimeError(f"OpenMM smoke produced non-finite {label}: {value}")


def validate_finite_positions(
    positions: Any, unit_module: Any | None = None, *, label: str
) -> None:
    """Validate that a coordinate array or OpenMM quantity has finite values."""
    array = _positions_to_numpy(positions, unit_module)
    if array.size == 0:
        raise RuntimeError(f"OpenMM smoke produced empty {label}")
    if not np.all(np.isfinite(array)):
        raise RuntimeError(f"OpenMM smoke produced non-finite {label}")


def _openmm_positions_from_interchange(interchange: Any, openmm_unit: Any) -> Any:
    """Extract OpenMM-compatible positions from an Interchange-like object."""
    for attr_name in ("positions",):
        positions = getattr(interchange, attr_name, None)
        if positions is not None:
            return _to_openmm_positions(positions, openmm_unit)
    topology = getattr(interchange, "topology", None)
    if topology is not None and hasattr(topology, "get_positions"):
        positions = topology.get_positions()
        if positions is not None:
            return _to_openmm_positions(positions, openmm_unit)
    raise RuntimeError("Interchange object does not expose positions for OpenMM smoke")


def _to_openmm_positions(positions: Any, openmm_unit: Any) -> Any:
    """Convert position containers to OpenMM nanometer quantities."""
    if hasattr(positions, "to_openmm"):
        return positions.to_openmm()
    if hasattr(positions, "m_as"):
        return positions.m_as("nanometer") * openmm_unit.nanometer
    return positions


def _positions_to_numpy(positions: Any, unit_module: Any | None) -> np.ndarray:
    """Convert coordinate containers to a NumPy array for finite checks."""
    if unit_module is not None and hasattr(positions, "value_in_unit"):
        try:
            return np.asarray(positions.value_in_unit(unit_module.nanometer), dtype=float)
        except Exception:  # noqa: BLE001 - mock and OpenMM quantity APIs vary
            pass
    if hasattr(positions, "m_as"):
        return np.asarray(positions.m_as("nanometer"), dtype=float)
    return np.asarray(positions, dtype=float)


def _state_energy_kj_mol(state: Any, openmm_unit: Any) -> float:
    """Return potential energy from an OpenMM state in kJ/mol."""
    energy = state.getPotentialEnergy()
    if hasattr(energy, "value_in_unit"):
        return float(energy.value_in_unit(openmm_unit.kilojoule_per_mole))
    return float(energy)


def _add_positional_restraints(
    system: Any,
    positions: Any,
    atom_indices: tuple[int, ...],
    restraint_k: float,
    openmm: Any,
    openmm_unit: Any,
) -> None:
    """Add harmonic positional restraints for selected atoms."""
    reference_nm = _positions_to_numpy(positions, openmm_unit)
    restraint = openmm.CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    restraint.addGlobalParameter(
        "k", restraint_k * openmm_unit.kilojoule_per_mole / openmm_unit.nanometer**2
    )
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for atom_index in atom_indices:
        x_coord, y_coord, z_coord = reference_nm[atom_index]
        restraint.addParticle(atom_index, [float(x_coord), float(y_coord), float(z_coord)])
    system.addForce(restraint)


def _infer_chain_a_heavy_indices(topology: Any) -> tuple[int, ...]:
    """Infer protein heavy atoms from chain A in an OpenMM topology."""
    indices: list[int] = []
    for atom in topology.atoms():
        chain_id = getattr(getattr(atom, "residue", None), "chain", None)
        chain_id = getattr(chain_id, "id", None)
        element_symbol = getattr(getattr(atom, "element", None), "symbol", "")
        if chain_id == "A" and element_symbol.upper() != "H":
            indices.append(atom.index)
    return tuple(indices)


def _select_platform(openmm: Any, requested_platform: str | None) -> Any:
    """Select an OpenMM platform, preferring accelerators when not specified."""
    names = (requested_platform,) if requested_platform else ("CUDA", "OpenCL", "CPU", "Reference")
    errors: list[str] = []
    for name in names:
        if name is None:
            continue
        try:
            return openmm.Platform.getPlatformByName(name)
        except Exception as exc:  # noqa: BLE001 - OpenMM platform errors vary
            errors.append(f"{name}: {exc}")
    raise RuntimeError(
        "No suitable OpenMM platform found for conjugation smoke: " + "; ".join(errors)
    )


def _write_openmm_pdb(openmm_app: Any, topology: Any, positions: Any, output_path: Path) -> None:
    """Write an OpenMM PDB artifact."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        openmm_app.PDBFile.writeFile(topology, positions, handle)
