"""Restrained vacuum OpenMM smoke runner for conjugation construction."""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path
from typing import Any

import numpy as np
from pydantic import BaseModel, Field

LOGGER = logging.getLogger(__name__)

_POSITION_CONVERSION_ERRORS = (AttributeError, TypeError, ValueError)


class VacuumSmokeSettings(BaseModel):
    """Settings for the short restrained vacuum smoke run."""

    minimization_max_iterations: int = Field(100, ge=0)
    minimization_tolerance_kj_mol_nm: float = Field(10.0, gt=0)
    nvt_steps: int = Field(10, ge=1)
    temperature_kelvin: float = Field(50.0, gt=0)
    timestep_femtoseconds: float = Field(0.25, gt=0)
    friction_per_picosecond: float = Field(10.0, gt=0)
    protein_heavy_restraint_k_kj_mol_nm2: float = Field(50000.0, gt=0)
    restrain_all_heavy_atoms: bool = True
    max_position_span_nm: float = Field(50.0, gt=0)
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
        Legacy chain-A atom indices to restrain when all-heavy restraints are
        disabled. When all-heavy restraints are enabled, the complete non-hydrogen
        atom set is inferred from the OpenMM topology and this selection cannot
        restrict the smoke to protein-only restraints.
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
    restrained_indices = _resolve_restrained_indices(
        topology,
        protein_heavy_atom_indices=protein_heavy_atom_indices,
        restrain_all_heavy_atoms=smoke_settings.restrain_all_heavy_atoms,
    )
    if not restrained_indices:
        raise RuntimeError("No heavy atoms were available for restrained vacuum smoke")

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
    minimized_span_nm = validate_finite_positions(
        minimized_positions,
        openmm_unit,
        label="minimized_positions",
        max_span_nm=smoke_settings.max_position_span_nm,
    )

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
    equilibrated_span_nm = validate_finite_positions(
        equilibrated_positions,
        openmm_unit,
        label="equilibrated_positions",
        max_span_nm=smoke_settings.max_position_span_nm,
    )

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
        diagnostics=(
            "Restrained vacuum OpenMM smoke completed with finite state values",
            f"Maximum minimized coordinate span was {minimized_span_nm:.3f} nm",
            f"Maximum equilibrated coordinate span was {equilibrated_span_nm:.3f} nm",
        ),
    )
    result.smoke_json_path.write_text(json.dumps(result.model_dump(mode="json"), indent=2) + "\n")
    return result


def validate_finite_energy(value: float, *, label: str = "energy") -> None:
    """Validate that an energy-like scalar is finite."""
    if not math.isfinite(float(value)):
        raise RuntimeError(f"OpenMM smoke produced non-finite {label}: {value}")


def validate_finite_positions(
    positions: Any,
    unit_module: Any | None = None,
    *,
    label: str,
    max_span_nm: float | None = 50.0,
) -> float:
    """Validate that coordinates are finite and contained in a realistic span.

    Parameters
    ----------
    positions : Any
        Coordinate array or OpenMM-like quantity in nanometers.
    unit_module : Any or None, optional
        OpenMM unit module used to extract nanometer values, by default ``None``.
    label : str
        Human-readable label included in validation errors.
    max_span_nm : float or None, optional
        Maximum allowed per-axis coordinate span in nanometers. Set to ``None``
        to skip the span check, by default 50.0.

    Returns
    -------
    float
        Maximum per-axis coordinate span in nanometers.
    """
    array = _positions_to_numpy(positions, unit_module)
    if array.size == 0:
        raise RuntimeError(f"OpenMM smoke produced empty {label}")
    if not np.all(np.isfinite(array)):
        raise RuntimeError(f"OpenMM smoke produced non-finite {label}")
    if array.ndim != 2 or array.shape[1] != 3:
        raise RuntimeError(f"OpenMM smoke produced invalid {label} shape: {array.shape}")
    span_nm = float(np.max(np.ptp(array, axis=0)))
    if max_span_nm is not None and span_nm > max_span_nm:
        raise RuntimeError(
            f"OpenMM smoke produced unrealistic coordinate span for {label}: "
            f"{span_nm:.3f} nm exceeds {max_span_nm:.3f} nm. "
            "The restrained vacuum MD likely became unstable; refusing to solvate blown-up "
            "post-MD coordinates."
        )
    return span_nm


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
    conversion_error: Exception | None = None
    if unit_module is not None and hasattr(positions, "value_in_unit"):
        try:
            return np.asarray(positions.value_in_unit(unit_module.nanometer), dtype=float)
        except _POSITION_CONVERSION_ERRORS as exc:
            conversion_error = exc
    if hasattr(positions, "m_as"):
        try:
            return np.asarray(positions.m_as("nanometer"), dtype=float)
        except _POSITION_CONVERSION_ERRORS as exc:
            conversion_error = exc
    if conversion_error is not None:
        LOGGER.warning(
            "Falling back to raw np.asarray() for positions of type %s after unit-aware "
            "coordinate conversion failed: %s",
            type(positions).__name__,
            conversion_error,
        )
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


def _resolve_restrained_indices(
    topology: Any,
    *,
    protein_heavy_atom_indices: tuple[int, ...] | None,
    restrain_all_heavy_atoms: bool,
) -> tuple[int, ...]:
    """Resolve restrained atom indices for a stable vacuum smoke."""
    if restrain_all_heavy_atoms:
        return _infer_heavy_indices(topology)
    selected = set(protein_heavy_atom_indices or _infer_chain_a_heavy_indices(topology))
    return tuple(sorted(selected))


def _infer_heavy_indices(topology: Any) -> tuple[int, ...]:
    """Infer all non-hydrogen atom indices from an OpenMM topology."""
    indices: list[int] = []
    for atom in topology.atoms():
        if _is_heavy_atom(atom):
            indices.append(atom.index)
    return tuple(indices)


def _infer_chain_a_heavy_indices(topology: Any) -> tuple[int, ...]:
    """Infer protein heavy atoms from chain A in an OpenMM topology."""
    indices: list[int] = []
    for atom in topology.atoms():
        chain_id = getattr(getattr(atom, "residue", None), "chain", None)
        chain_id = getattr(chain_id, "id", None)
        if chain_id == "A" and _is_heavy_atom(atom):
            indices.append(atom.index)
    return tuple(indices)


def _is_heavy_atom(atom: Any) -> bool:
    """Return whether an OpenMM atom-like object is not hydrogen."""
    element_symbol = getattr(getattr(atom, "element", None), "symbol", "")
    atom_name = getattr(atom, "name", "")
    return element_symbol.upper() != "H" and not str(atom_name).upper().startswith("H")


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
