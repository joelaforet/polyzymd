"""OpenMM system helpers for conjugate relaxation."""

from __future__ import annotations

import logging
import math
import os
from pathlib import Path
from typing import Any

import numpy as np

from polyzymd.builders.conjugation.relaxation._diagnostics import (
    validate_finite_energy,
    validate_finite_positions,
)
from polyzymd.builders.conjugation.relaxation.models import (
    ConjugateRelaxationSettings,
    ProductLinkage,
)

LOGGER = logging.getLogger(__name__)


def _freeze_protein_chain_masses(
    system: Any,
    topology: Any,
    openmm_unit: Any,
    *,
    chain_ids: tuple[str, ...] = ("A",),
) -> tuple[tuple[int, ...], dict[int, Any]]:
    """Set selected protein-chain particle masses to zero.

    Parameters
    ----------
    system : Any
        OpenMM-like system with particle mass accessors.
    topology : Any
        OpenMM-like topology exposing atoms with residue chain IDs.
    openmm_unit : Any
        OpenMM unit module or compatible double.
    chain_ids : tuple of str, optional
        Protein chain IDs to freeze, by default ``("A",)``.

    Returns
    -------
    tuple[tuple[int, ...], dict[int, Any]]
        Fixed particle indices and their original masses.
    """
    indices = _protein_chain_indices(topology, chain_ids)
    original = {index: system.getParticleMass(index) for index in indices}
    zero_mass = 0.0 * openmm_unit.dalton
    for index in indices:
        system.setParticleMass(index, zero_mass)
    return indices, original


def _protein_chain_indices(topology: Any, chain_ids: tuple[str, ...]) -> tuple[int, ...]:
    """Return topology atom indices belonging to selected protein chains."""
    selected = {chain_id.strip() for chain_id in chain_ids if chain_id.strip()}
    return tuple(
        int(atom.index)
        for atom in topology.atoms()
        if str(getattr(getattr(getattr(atom, "residue", None), "chain", None), "id", ""))
        in selected
    )


def _restore_particle_masses(system: Any, original_masses: dict[int, Any]) -> None:
    """Restore particle masses captured before temporary freezing."""
    for index, mass in original_masses.items():
        system.setParticleMass(index, mass)


def _set_zero_initial_velocities(context: Any, topology: Any, openmm_unit: Any) -> None:
    """Set finite initial velocities before the thermostat heats mobile atoms."""
    atom_count = len(tuple(topology.atoms()))
    velocities = (
        np.zeros((atom_count, 3), dtype=float) * openmm_unit.nanometer / openmm_unit.picosecond
    )
    context.setVelocities(velocities)


def _create_relaxation_verlet_integrator(openmm: Any, openmm_unit: Any) -> Any:
    """Create the Verlet integrator used for transient relaxation setup.

    Parameters
    ----------
    openmm : Any
        Imported OpenMM module or compatible test double.
    openmm_unit : Any
        Imported OpenMM unit module or compatible test double.

    Returns
    -------
    Any
        OpenMM Verlet integrator with the established relaxation timestep.
    """
    return openmm.VerletIntegrator(0.001 * openmm_unit.picoseconds)


def _create_simulation(
    topology: Any,
    system: Any,
    integrator: Any,
    openmm_app: Any,
    platform: Any,
) -> Any:
    """Create an OpenMM simulation with the selected platform.

    Parameters
    ----------
    topology : Any
        OpenMM-compatible topology.
    system : Any
        OpenMM-compatible system.
    integrator : Any
        OpenMM-compatible integrator.
    openmm_app : Any
        Imported ``openmm.app`` module or compatible test double.
    platform : Any
        Selected OpenMM platform.

    Returns
    -------
    Any
        OpenMM simulation object.
    """
    return openmm_app.Simulation(topology, system, integrator, platform)


def _run_fixed_product_md(
    topology: Any,
    system: Any,
    positions: Any,
    settings: ConjugateRelaxationSettings,
    openmm: Any,
    openmm_app: Any,
    openmm_unit: Any,
    platform: Any,
    warnings: list[str],
) -> tuple[float, float, Any]:
    """Run finite fixed-protein product MD, retrying with safer timesteps if needed."""
    timestep_schedule = _md_timestep_retry_schedule(settings.timestep_femtoseconds)
    last_error: Exception | None = None
    for timestep_fs in timestep_schedule:
        try:
            integrator = openmm.LangevinMiddleIntegrator(
                settings.temperature_kelvin * openmm_unit.kelvin,
                settings.friction_per_picosecond / openmm_unit.picosecond,
                timestep_fs * openmm_unit.femtosecond,
            )
            simulation = _create_simulation(topology, system, integrator, openmm_app, platform)
            simulation.context.setPositions(positions)
            _set_zero_initial_velocities(simulation.context, topology, openmm_unit)
            energy_before = _state_energy_kj_mol(
                simulation.context.getState(getEnergy=True),
                openmm_unit,
            )
            validate_finite_energy(energy_before, label="energy_before_md_kj_mol")
            simulation.step(settings.md_steps)
            state_after = simulation.context.getState(getEnergy=True, getPositions=True)
            energy_after = _state_energy_kj_mol(state_after, openmm_unit)
            validate_finite_energy(energy_after, label="energy_after_md_kj_mol")
            final_positions = state_after.getPositions(asNumpy=True)
            validate_finite_positions(final_positions, openmm_unit, label="relaxed_positions")
            if timestep_fs != settings.timestep_femtoseconds:
                warnings.append(
                    "Conjugate MD retried with a smaller timestep "
                    f"({timestep_fs:.3f} fs) after instability at the requested timestep"
                )
            return energy_before, energy_after, final_positions
        except _openmm_runtime_exceptions(openmm) as exc:
            last_error = exc
            if timestep_fs == timestep_schedule[-1]:
                break
            warnings.append(
                f"Conjugate MD was unstable at {timestep_fs:.3f} fs; retrying with a "
                "smaller timestep"
            )
    if last_error is None:
        raise RuntimeError("Conjugate MD did not run")
    raise last_error


def _openmm_runtime_exceptions(openmm: Any) -> tuple[type[Exception], ...]:
    """Return expected OpenMM runtime instability exception classes.

    Parameters
    ----------
    openmm : Any
        Imported OpenMM module or compatible test double.

    Returns
    -------
    tuple[type[Exception], ...]
        Exception classes that represent runtime or platform instability.
    """
    exceptions: list[type[Exception]] = [RuntimeError, OSError]
    openmm_exception = getattr(openmm, "OpenMMException", None)
    if isinstance(openmm_exception, type) and issubclass(openmm_exception, Exception):
        exceptions.append(openmm_exception)
    return tuple(exceptions)


def _md_timestep_retry_schedule(requested_femtoseconds: float) -> tuple[float, ...]:
    """Return requested timestep followed by conservative relaxation fallbacks."""
    candidates = (requested_femtoseconds, 1.0, 0.5, 0.25)
    schedule: list[float] = []
    for value in candidates:
        if value <= requested_femtoseconds and value not in schedule:
            schedule.append(value)
    return tuple(schedule)


def _openmm_positions_from_interchange_or_pdb(
    interchange: Any,
    openmm_app: Any,
    openmm_unit: Any,
    product_pdb_path: Path | str,
) -> Any:
    """Extract positions from Interchange, falling back to the product PDB."""
    try:
        return _openmm_positions_from_interchange(interchange, openmm_unit)
    except RuntimeError:
        pdb = openmm_app.PDBFile(str(product_pdb_path))
        return pdb.positions


def _remove_barostats(system: Any) -> int:
    """Remove barostat-like forces from a transient relaxation system."""
    removed = 0
    for index in reversed(range(int(system.getNumForces()))):
        force = system.getForce(index)
        if "Barostat" in type(force).__name__:
            system.removeForce(index)
            removed += 1
    return removed


def _add_linkage_anchor_restraints(
    system: Any,
    pairs: tuple[ProductLinkage, ...],
    force_constant_kj_mol_nm2: float,
    openmm: Any,
    openmm_unit: Any,
) -> int:
    """Add temporary harmonic distance anchors for resolved product linkages."""
    if not pairs:
        return 0
    force = openmm.CustomBondForce("0.5*k*(r-r0)^2")
    force.addPerBondParameter("r0")
    force.addGlobalParameter(
        "k",
        force_constant_kj_mol_nm2
        * openmm_unit.kilojoule_per_mole
        / (openmm_unit.nanometer * openmm_unit.nanometer),
    )
    for pair in pairs:
        force.addBond(
            pair.protein_atom_index,
            pair.modifier_atom_index,
            [pair.target_bond_length_angstrom * 0.1 * openmm_unit.nanometer],
        )
    system.addForce(force)
    return len(pairs)


def _protein_displacements_angstrom(
    initial_nm: np.ndarray,
    final_nm: np.ndarray,
    indices: tuple[int, ...],
) -> tuple[float, float]:
    """Return RMSD and maximum displacement for fixed protein atoms."""
    if not indices:
        return 0.0, 0.0
    delta_angstrom = (final_nm[list(indices)] - initial_nm[list(indices)]) * 10.0
    distances = np.linalg.norm(delta_angstrom, axis=1)
    return float(np.sqrt(np.mean(distances**2))), float(np.max(distances))


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
    raise RuntimeError("Interchange object does not expose positions for OpenMM relaxation")


def _to_openmm_positions(positions: Any, openmm_unit: Any) -> Any:
    """Convert position containers to OpenMM nanometer quantities."""
    if hasattr(positions, "to_openmm"):
        return positions.to_openmm()
    if hasattr(positions, "m_as"):
        return positions.m_as("nanometer") * openmm_unit.nanometer
    return positions


def _assign_force_groups(system: Any) -> dict[int, str]:
    """Assign one force group per force when supported."""
    labels: dict[int, str] = {}
    try:
        force_count = int(system.getNumForces())
    except (AttributeError, TypeError, ValueError):
        return labels
    if force_count > 31:
        return labels
    for index in range(force_count):
        force = system.getForce(index)
        try:
            force.setForceGroup(index)
            labels[index] = _force_label(force, index)
        except (AttributeError, ValueError):
            continue
    return labels


def _force_group_labels(system: Any, *, existing_labels: dict[int, str]) -> dict[int, str]:
    """Return force group labels after adding restraints."""
    labels = dict(existing_labels)
    try:
        force_count = int(system.getNumForces())
    except (AttributeError, TypeError, ValueError):
        return labels
    if force_count > 31:
        return labels
    for index in range(force_count):
        force = system.getForce(index)
        try:
            group = int(force.getForceGroup())
        except (AttributeError, TypeError, ValueError):
            group = index
            try:
                force.setForceGroup(group)
            except (AttributeError, ValueError):
                continue
        if group in labels and labels[group] != _force_label(force, index):
            group = index
            try:
                force.setForceGroup(group)
            except (AttributeError, ValueError):
                continue
        labels.setdefault(group, _force_label(force, index))
    return labels


def _force_label(force: Any, index: int) -> str:
    """Return a stable diagnostic label for an OpenMM force."""
    class_name = type(force).__name__
    if "Bond" in class_name and "Nonbonded" not in class_name:
        prefix = "bond"
    elif "Angle" in class_name:
        prefix = "angle"
    elif "Torsion" in class_name:
        prefix = "torsion"
    elif "Nonbonded" in class_name:
        prefix = "nonbonded"
    elif "External" in class_name:
        prefix = "restraint"
    else:
        prefix = "other"
    return f"{prefix}:{index}:{class_name}"


def _force_group_energies(
    context: Any,
    group_labels: dict[int, str],
    openmm_unit: Any,
) -> dict[str, float]:
    """Return potential energies per OpenMM force group with graceful fallback."""
    energies: dict[str, float] = {}
    for group, label in sorted(group_labels.items()):
        try:
            state = context.getState(getEnergy=True, groups={group})
            energies[label] = _state_energy_kj_mol(state, openmm_unit)
        except (AttributeError, TypeError, ValueError, RuntimeError) as exc:
            LOGGER.debug("Could not collect force-group energy for %s: %s", label, exc)
    return energies


def _state_energy_kj_mol(state: Any, openmm_unit: Any) -> float:
    """Return potential energy from an OpenMM state in kJ/mol."""
    energy = state.getPotentialEnergy()
    if hasattr(energy, "value_in_unit"):
        return float(energy.value_in_unit(openmm_unit.kilojoule_per_mole))
    return float(energy)


def _select_platform(openmm: Any, requested_platform: str | None) -> Any:
    """Select a usable OpenMM platform, preferring accelerators when available."""
    requested = requested_platform or os.environ.get("POLYZYMD_CONJUGATION_VALIDATION_PLATFORM")
    names = (requested,) if requested else ("CUDA", "OpenCL", "CPU", "Reference")
    errors: list[str] = []
    for name in names:
        if name is None:
            continue
        try:
            platform = openmm.Platform.getPlatformByName(name)
            _validate_platform_context(openmm, platform)
            return platform
        except (AttributeError, TypeError, ValueError, RuntimeError) as exc:
            errors.append(f"{name}: {exc}")
    raise RuntimeError(
        "No suitable OpenMM platform found for conjugation validation: " + "; ".join(errors)
    )


def _validate_platform_context(openmm: Any, platform: Any) -> None:
    """Create a tiny context to confirm the registered platform is usable."""
    from openmm import unit as openmm_unit

    system = openmm.System()
    system.addParticle(39.9)
    integrator = _create_relaxation_verlet_integrator(openmm, openmm_unit)
    context = openmm.Context(system, integrator, platform)
    context.setPositions([[0.0, 0.0, 0.0]] * openmm_unit.nanometer)
    context.getState(getEnergy=True)
    del context
    del integrator


def _write_openmm_pdb(openmm_app: Any, topology: Any, positions: Any, output_path: Path) -> None:
    """Write an OpenMM PDB artifact."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        openmm_app.PDBFile.writeFile(topology, positions, handle)
