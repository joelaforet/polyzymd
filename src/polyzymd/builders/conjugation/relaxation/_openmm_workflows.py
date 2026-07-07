"""Private OpenMM workflow execution for conjugate relaxation and validation."""

from __future__ import annotations

import json
import logging
import math
import traceback
from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation.relaxation._diagnostics import (
    _first_invalid_phase,
    _invalid_phase_reason,
    _positions_to_numpy,
    _state_max_force_kj_mol_nm,
    _validation_phase_diagnostics,
    analyze_product_geometry,
    validate_finite_energy,
    validate_finite_positions,
)
from polyzymd.builders.conjugation.relaxation._linkages import (
    _linkage_distances_angstrom,
    resolve_product_linkage_pairs,
)
from polyzymd.builders.conjugation.relaxation._openmm_system import (
    _add_linkage_anchor_restraints,
    _add_positional_restraints,
    _assign_force_groups,
    _force_group_energies,
    _force_group_labels,
    _freeze_protein_chain_masses,
    _openmm_positions_from_interchange,
    _openmm_positions_from_interchange_or_pdb,
    _protein_chain_indices,
    _protein_displacements_angstrom,
    _remove_barostats,
    _resolve_restrained_indices,
    _restore_particle_masses,
    _run_fixed_product_md,
    _select_platform,
    _state_energy_kj_mol,
    _validation_settings_from_environment,
    _write_openmm_pdb,
)
from polyzymd.builders.conjugation.relaxation.models import (
    ConjugateRelaxationDiagnostics,
    ConjugateRelaxationResult,
    ConjugateRelaxationSettings,
    OpenMMValidationPhaseDiagnostics,
    OpenMMValidationResult,
    OpenMMValidationSettings,
)

LOGGER = logging.getLogger(__name__)


def _openmm_stage_errors(openmm: Any) -> tuple[type[BaseException], ...]:
    """Return expected exception types raised by OpenMM workflow stages."""
    openmm_error = getattr(openmm, "OpenMMException", RuntimeError)
    return (RuntimeError, ValueError, ArithmeticError, openmm_error)


def _finite_or_none(value: float) -> float | None:
    """Return a finite float or ``None`` for incomplete failed stages."""
    return float(value) if math.isfinite(float(value)) else None


def _write_validation_result(result: OpenMMValidationResult) -> OpenMMValidationResult:
    """Persist the canonical OpenMM validation result once and return it."""
    result.write_json(result.validation_json_path)
    return result


def run_openmm_validation_workflow(
    interchange: Any,
    output_dir: Path | str,
    *,
    protein_heavy_atom_indices: tuple[int, ...] | None = None,
    settings: OpenMMValidationSettings | None = None,
    crosslinked_pdb_path: Path | str | None = None,
    attachment_specs: tuple[Any, ...] = (),
) -> OpenMMValidationResult:
    """Run restrained minimization and short OpenMM product validation.

    Parameters
    ----------
    interchange : Any
        OpenFF Interchange-like object exposing OpenMM conversion methods.
    output_dir : pathlib.Path or str
        Directory for validation artifacts.
    protein_heavy_atom_indices : tuple of int or None, optional
        Protein heavy-atom indices to restrain when all-heavy restraints are
        disabled. When all-heavy restraints are enabled, the complete
        non-hydrogen atom set is inferred from the OpenMM topology.
    settings : OpenMMValidationSettings or None, optional
        Validation settings, by default ``None``.
    crosslinked_pdb_path : pathlib.Path, str, or None, optional
        Product PDB used to measure resolved crosslink lengths, by default
        ``None``.
    attachment_specs : tuple of Any, optional
        Resolved attachment build specs used for crosslink-specific diagnostics,
        by default ``()``.

    Returns
    -------
    OpenMMValidationResult
        Energies and artifact paths for the OpenMM validation.

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
            "OpenMM is required for the conjugation OpenMM validation. Run under a suitable pixi "
            "environment such as cuda-12-4 when GPU resources are allocated."
        ) from exc

    validation_settings = _validation_settings_from_environment(
        settings or OpenMMValidationSettings()
    )
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    topology = interchange.to_openmm_topology()
    initial_positions = _openmm_positions_from_interchange(interchange, openmm_unit)
    restrained_indices = _resolve_restrained_indices(
        topology,
        protein_heavy_atom_indices=protein_heavy_atom_indices,
        restrain_all_heavy_atoms=validation_settings.restrain_all_heavy_atoms,
    )
    if not restrained_indices:
        raise RuntimeError("No heavy atoms were available for restrained OpenMM validation")

    platform = _select_platform(openmm, validation_settings.platform_name)
    platform_name = platform.getName() if hasattr(platform, "getName") else str(platform)
    diagnostics_path = artifact_dir / validation_settings.validation_json_name
    phase_diagnostics: list[OpenMMValidationPhaseDiagnostics] = []
    validation_segment_diagnostics: list[OpenMMValidationPhaseDiagnostics] = []
    first_invalid_phase: str | None = None

    energy_before_min = math.nan
    energy_after_min = math.nan
    energy_before_nvt = math.nan
    energy_after_nvt = math.nan
    failure_json_path: Path | None = None
    force_energies_before: dict[str, float] = {}
    force_energies_after: dict[str, float] = {}
    product_geometry = analyze_product_geometry(
        topology,
        initial_positions,
        openmm_unit,
        crosslinked_pdb_path=crosslinked_pdb_path,
        attachment_specs=attachment_specs,
    )
    product_geometry_path = artifact_dir / validation_settings.product_geometry_json_name
    product_geometry.write_json(product_geometry_path)
    phase_diagnostics.append(
        _validation_phase_diagnostics(
            "before_minimization",
            topology,
            initial_positions,
            openmm_unit,
            potential_energy_kj_mol=None,
            attachment_specs=attachment_specs,
        )
    )

    try:
        system_min = interchange.to_openmm_system()
        group_labels_min = _assign_force_groups(system_min)
        _add_positional_restraints(
            system_min,
            initial_positions,
            restrained_indices,
            validation_settings.protein_heavy_restraint_k_kj_mol_nm2,
            openmm,
            openmm_unit,
        )
        group_labels_min = _force_group_labels(system_min, existing_labels=group_labels_min)
        integrator_min = openmm.VerletIntegrator(0.001 * openmm_unit.picoseconds)
        simulation_min = openmm_app.Simulation(topology, system_min, integrator_min, platform)
        simulation_min.context.setPositions(initial_positions)

        energy_before_min = _state_energy_kj_mol(
            simulation_min.context.getState(getEnergy=True), openmm_unit
        )
        validate_finite_energy(energy_before_min, label="energy_before_min_kj_mol")
        phase_diagnostics[-1].potential_energy_kj_mol = float(energy_before_min)
        force_energies_before = _force_group_energies(
            simulation_min.context,
            group_labels_min,
            openmm_unit,
        )
        openmm.LocalEnergyMinimizer.minimize(
            simulation_min.context,
            tolerance=validation_settings.minimization_tolerance_kj_mol_nm
            * openmm_unit.kilojoule_per_mole
            / openmm_unit.nanometer,
            maxIterations=validation_settings.minimization_max_iterations,
        )
        state_after_min = simulation_min.context.getState(
            getEnergy=True, getPositions=True, getForces=True
        )
        energy_after_min = _state_energy_kj_mol(state_after_min, openmm_unit)
        minimized_positions = state_after_min.getPositions(asNumpy=True)
        validate_finite_energy(energy_after_min, label="energy_after_min_kj_mol")
        phase_diagnostics.append(
            _validation_phase_diagnostics(
                "after_minimization",
                topology,
                minimized_positions,
                openmm_unit,
                potential_energy_kj_mol=energy_after_min,
                max_force_kj_mol_nm=_state_max_force_kj_mol_nm(state_after_min, openmm_unit),
                attachment_specs=attachment_specs,
            )
        )
        minimized_span_nm = validate_finite_positions(
            minimized_positions,
            openmm_unit,
            label="minimized_positions",
            max_span_nm=validation_settings.max_position_span_nm,
        )
        force_energies_after = _force_group_energies(
            simulation_min.context,
            group_labels_min,
            openmm_unit,
        )

        if validation_settings.nvt_steps == 0:
            energy_before_nvt = energy_after_min
            energy_after_nvt = energy_after_min
            equilibrated_positions = minimized_positions
            equilibrated_span_nm = minimized_span_nm
        else:
            system_nvt = interchange.to_openmm_system()
            _add_positional_restraints(
                system_nvt,
                minimized_positions,
                restrained_indices,
                validation_settings.protein_heavy_restraint_k_kj_mol_nm2,
                openmm,
                openmm_unit,
            )
            integrator_nvt = openmm.LangevinMiddleIntegrator(
                validation_settings.temperature_kelvin * openmm_unit.kelvin,
                validation_settings.friction_per_picosecond / openmm_unit.picosecond,
                validation_settings.timestep_femtoseconds * openmm_unit.femtosecond,
            )
            simulation_nvt = openmm_app.Simulation(topology, system_nvt, integrator_nvt, platform)
            simulation_nvt.context.setPositions(minimized_positions)
            simulation_nvt.context.setVelocitiesToTemperature(
                validation_settings.temperature_kelvin * openmm_unit.kelvin
            )

            state_before_nvt = simulation_nvt.context.getState(getEnergy=True, getPositions=True)
            energy_before_nvt = _state_energy_kj_mol(state_before_nvt, openmm_unit)
            validate_finite_energy(energy_before_nvt, label="energy_before_nvt_kj_mol")
            phase_diagnostics.append(
                _validation_phase_diagnostics(
                    "before_md",
                    topology,
                    state_before_nvt.getPositions(asNumpy=True),
                    openmm_unit,
                    potential_energy_kj_mol=energy_before_nvt,
                    attachment_specs=attachment_specs,
                )
            )
            for step_index in range(validation_settings.nvt_steps):
                simulation_nvt.step(1)
                segment_state = simulation_nvt.context.getState(getEnergy=True, getPositions=True)
                segment_energy = _state_energy_kj_mol(segment_state, openmm_unit)
                segment_positions = segment_state.getPositions(asNumpy=True)
                segment = _validation_phase_diagnostics(
                    f"after_md_step_{step_index + 1}",
                    topology,
                    segment_positions,
                    openmm_unit,
                    potential_energy_kj_mol=segment_energy,
                    attachment_specs=attachment_specs,
                )
                validation_segment_diagnostics.append(segment)
                if first_invalid_phase is None:
                    first_invalid_phase = _invalid_phase_reason(
                        segment,
                        max_span_nm=validation_settings.max_position_span_nm,
                    )
            state_after_nvt = simulation_nvt.context.getState(getEnergy=True, getPositions=True)
            energy_after_nvt = _state_energy_kj_mol(state_after_nvt, openmm_unit)
            equilibrated_positions = state_after_nvt.getPositions(asNumpy=True)
            validate_finite_energy(energy_after_nvt, label="energy_after_nvt_kj_mol")
            phase_diagnostics.append(
                _validation_phase_diagnostics(
                    "after_md",
                    topology,
                    equilibrated_positions,
                    openmm_unit,
                    potential_energy_kj_mol=energy_after_nvt,
                    attachment_specs=attachment_specs,
                )
            )
            if first_invalid_phase is None:
                first_invalid_phase = _invalid_phase_reason(
                    phase_diagnostics[-1],
                    max_span_nm=validation_settings.max_position_span_nm,
                )
            equilibrated_span_nm = validate_finite_positions(
                equilibrated_positions,
                openmm_unit,
                label="equilibrated_positions",
                max_span_nm=validation_settings.max_position_span_nm,
            )
    except _openmm_stage_errors(openmm) as exc:
        if first_invalid_phase is None:
            first_invalid_phase = _first_invalid_phase(
                (*phase_diagnostics, *validation_segment_diagnostics),
                max_span_nm=validation_settings.max_position_span_nm,
            )
        result = OpenMMValidationResult(
            success=False,
            output_dir=artifact_dir,
            validation_json_path=diagnostics_path,
            diagnostics_json_path=diagnostics_path,
            product_geometry_json_path=product_geometry_path,
            failure_json_path=failure_json_path,
            platform_name=platform_name,
            restrained_atom_count=len(restrained_indices),
            energy_before_min_kj_mol=_finite_or_none(energy_before_min),
            energy_after_min_kj_mol=_finite_or_none(energy_after_min),
            energy_before_nvt_kj_mol=_finite_or_none(energy_before_nvt),
            energy_after_nvt_kj_mol=_finite_or_none(energy_after_nvt),
            force_group_energies_before_min_kj_mol=force_energies_before,
            force_group_energies_after_min_kj_mol=force_energies_after,
            nvt_steps=validation_settings.nvt_steps,
            settings=validation_settings.model_dump(mode="json"),
            phases=tuple(phase_diagnostics),
            validation_segments=tuple(validation_segment_diagnostics),
            first_invalid_phase=first_invalid_phase,
            error_type=type(exc).__name__,
            error_message=str(exc),
            traceback="".join(traceback.format_exception(type(exc), exc, exc.__traceback__)),
            diagnostics=("OpenMM validation failed before completion",),
        )
        _write_validation_result(result)
        raise

    result = OpenMMValidationResult(
        success=True,
        output_dir=artifact_dir,
        validation_json_path=artifact_dir / validation_settings.validation_json_name,
        diagnostics_json_path=diagnostics_path,
        product_geometry_json_path=product_geometry_path,
        failure_json_path=failure_json_path,
        platform_name=platform_name,
        restrained_atom_count=len(restrained_indices),
        energy_before_min_kj_mol=float(energy_before_min),
        energy_after_min_kj_mol=float(energy_after_min),
        energy_before_nvt_kj_mol=float(energy_before_nvt),
        energy_after_nvt_kj_mol=float(energy_after_nvt),
        minimized_coordinate_span_nm=minimized_span_nm,
        equilibrated_coordinate_span_nm=equilibrated_span_nm,
        force_group_energies_before_min_kj_mol=force_energies_before,
        force_group_energies_after_min_kj_mol=force_energies_after,
        nvt_steps=validation_settings.nvt_steps,
        settings=validation_settings.model_dump(mode="json"),
        phases=tuple(phase_diagnostics),
        validation_segments=tuple(validation_segment_diagnostics),
        first_invalid_phase=None,
        diagnostics=(
            "Short restrained OpenMM validation completed with finite state values",
            f"Maximum minimized coordinate span was {minimized_span_nm:.3f} nm",
            (
                "NVT dynamics were skipped because nvt_steps=0"
                if validation_settings.nvt_steps == 0
                else f"Maximum equilibrated coordinate span was {equilibrated_span_nm:.3f} nm"
            ),
        ),
    )
    return _write_validation_result(result)


def run_conjugate_relaxation_workflow(
    interchange: Any,
    output_dir: Path | str,
    *,
    product_pdb_path: Path | str,
    attachment_specs: tuple[Any, ...],
    assembly: Any | None = None,
    settings: ConjugateRelaxationSettings | None = None,
) -> ConjugateRelaxationResult:
    """Relax a product-state conjugate with staged minimization then fixed-protein MD.

    Parameters
    ----------
    interchange : Any
        OpenFF Interchange-like object exposing OpenMM conversion methods.
    output_dir : pathlib.Path or str
        Directory for relaxation artifacts.
    product_pdb_path : pathlib.Path or str
        Product PDB carrying assembly serial and residue mapping metadata.
    attachment_specs : tuple of Any
        Resolved attachment build specs used to resolve generic product linkage atoms.
    assembly : Any or None, optional
        Assembly result with ``added_conect_pairs`` and residue mappings, by default ``None``.
    settings : ConjugateRelaxationSettings or None, optional
        Relaxation settings, by default ``None``.

    Returns
    -------
    ConjugateRelaxationResult
        Energies, validation metrics, and artifact paths.
    """
    try:
        import openmm
        from openmm import app as openmm_app
        from openmm import unit as openmm_unit
    except ImportError as exc:
        raise RuntimeError(
            "OpenMM is required for conjugate conjugation relaxation. Run under a "
            "suitable PolyzyMD pixi environment."
        ) from exc

    relaxation_settings = settings or ConjugateRelaxationSettings()
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_path = artifact_dir / relaxation_settings.diagnostics_json_name
    failure_path: Path | None = None
    platform = _select_platform(openmm, relaxation_settings.platform_name)
    platform_name = platform.getName() if hasattr(platform, "getName") else str(platform)
    warnings: list[str] = []

    topology = interchange.to_openmm_topology()
    initial_positions = _openmm_positions_from_interchange_or_pdb(
        interchange,
        openmm_app,
        openmm_unit,
        product_pdb_path,
    )
    initial_coords_nm = _positions_to_numpy(initial_positions, openmm_unit)
    linkage_pairs = resolve_product_linkage_pairs(
        topology,
        product_pdb_path=product_pdb_path,
        attachment_specs=attachment_specs,
        assembly=assembly,
        fallback_bond_length_angstrom=relaxation_settings.fallback_bond_length_angstrom,
        warnings=warnings,
    )

    protein_indices = _protein_chain_indices(topology, relaxation_settings.protein_chain_ids)
    if not protein_indices:
        chains = ", ".join(relaxation_settings.protein_chain_ids)
        raise RuntimeError(
            "Conjugate relaxation could not identify protein chain atoms to freeze "
            f"for chain IDs: {chains}"
        )

    energy_before_min = math.nan
    energy_after_min = math.nan
    energy_before_md = math.nan
    energy_after_md = math.nan
    force_energies_before_min: dict[str, float] = {}
    force_energies_after_min: dict[str, float] = {}
    stage_a_protein_rmsd = math.nan
    stage_a_protein_max = math.nan
    stage_b_protein_rmsd = math.nan
    stage_b_protein_max = math.nan
    initial_to_final_protein_rmsd = math.nan
    initial_to_final_protein_max = math.nan
    stage_a_distances: tuple[float, ...] = ()
    distances: tuple[float, ...] = ()
    errors: tuple[float, ...] = ()
    removed_barostats = 0
    fixed_indices: tuple[int, ...] = ()
    anchor_count = 0
    relaxed_pdb = artifact_dir / relaxation_settings.relaxed_pdb_name
    try:
        LOGGER.info("Running Stage A normal-mass full-system minimization for conjugate product")
        system_min = interchange.to_openmm_system()
        removed_barostats += _remove_barostats(system_min)
        group_labels_min = _assign_force_groups(system_min)
        anchor_count = _add_linkage_anchor_restraints(
            system_min,
            linkage_pairs,
            relaxation_settings.anchor_k_kj_mol_nm2,
            openmm,
            openmm_unit,
        )
        group_labels_min = _force_group_labels(system_min, existing_labels=group_labels_min)
        integrator_min = openmm.VerletIntegrator(0.001 * openmm_unit.picoseconds)
        simulation = openmm_app.Simulation(topology, system_min, integrator_min, platform)
        simulation.context.setPositions(initial_positions)
        energy_before_min = _state_energy_kj_mol(
            simulation.context.getState(getEnergy=True),
            openmm_unit,
        )
        validate_finite_energy(energy_before_min, label="energy_before_min_kj_mol")
        force_energies_before_min = _force_group_energies(
            simulation.context,
            group_labels_min,
            openmm_unit,
        )
        openmm.LocalEnergyMinimizer.minimize(
            simulation.context,
            tolerance=relaxation_settings.minimization_tolerance_kj_mol_nm
            * openmm_unit.kilojoule_per_mole
            / openmm_unit.nanometer,
            maxIterations=relaxation_settings.minimization_max_iterations,
        )
        min_state = simulation.context.getState(getEnergy=True, getPositions=True)
        energy_after_min = _state_energy_kj_mol(min_state, openmm_unit)
        validate_finite_energy(energy_after_min, label="energy_after_min_kj_mol")
        minimized_positions = min_state.getPositions(asNumpy=True)
        validate_finite_positions(minimized_positions, openmm_unit, label="minimized_positions")
        force_energies_after_min = _force_group_energies(
            simulation.context,
            group_labels_min,
            openmm_unit,
        )
        minimized_coords_nm = _positions_to_numpy(minimized_positions, openmm_unit)
        stage_a_protein_rmsd, stage_a_protein_max = _protein_displacements_angstrom(
            initial_coords_nm,
            minimized_coords_nm,
            protein_indices,
        )
        stage_a_distances = _linkage_distances_angstrom(minimized_coords_nm, linkage_pairs)

        LOGGER.info("Running Stage B restrained OpenMM relaxation for conjugate product")
        system_md = interchange.to_openmm_system()
        removed_barostats += _remove_barostats(system_md)
        fixed_indices, _original_masses = _freeze_protein_chain_masses(
            system_md,
            topology,
            openmm_unit,
            chain_ids=relaxation_settings.protein_chain_ids,
        )
        anchor_count = _add_linkage_anchor_restraints(
            system_md,
            linkage_pairs,
            relaxation_settings.anchor_k_kj_mol_nm2,
            openmm,
            openmm_unit,
        )
        try:
            energy_before_md, energy_after_md, final_positions = _run_fixed_product_md(
                topology,
                system_md,
                minimized_positions,
                relaxation_settings,
                openmm,
                openmm_app,
                openmm_unit,
                platform,
                warnings,
            )
        finally:
            _restore_particle_masses(system_md, _original_masses)

        _write_openmm_pdb(openmm_app, topology, final_positions, relaxed_pdb)
        final_coords_nm = _positions_to_numpy(final_positions, openmm_unit)
        stage_b_protein_rmsd, stage_b_protein_max = _protein_displacements_angstrom(
            minimized_coords_nm,
            final_coords_nm,
            fixed_indices,
        )
        initial_to_final_protein_rmsd, initial_to_final_protein_max = (
            _protein_displacements_angstrom(
                initial_coords_nm,
                final_coords_nm,
                fixed_indices,
            )
        )
        distances = _linkage_distances_angstrom(final_coords_nm, linkage_pairs)
        errors = tuple(
            abs(distance - pair.target_bond_length_angstrom)
            for distance, pair in zip(distances, linkage_pairs, strict=True)
        )
        if stage_b_protein_rmsd > relaxation_settings.max_protein_rmsd_angstrom:
            warnings.append(
                f"Stage B protein RMSD {stage_b_protein_rmsd:.4f} A relative to Stage A exceeds "
                f"{relaxation_settings.max_protein_rmsd_angstrom:.4f} A"
            )
        if stage_b_protein_max > relaxation_settings.max_protein_displacement_angstrom:
            warnings.append(
                f"Stage B protein max displacement {stage_b_protein_max:.4f} A relative to Stage A "
                "exceeds "
                f"{relaxation_settings.max_protein_displacement_angstrom:.4f} A"
            )
        if any(error > relaxation_settings.max_linkage_distance_error_angstrom for error in errors):
            warnings.append("One or more linkage distances deviates from its target")
    except _openmm_stage_errors(openmm) as exc:
        failure_path = artifact_dir / relaxation_settings.failure_json_name
        failure_payload = {
            "success": False,
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "traceback": traceback.format_exc(),
            "energy_before_min_kj_mol": energy_before_min,
            "energy_after_min_kj_mol": energy_after_min,
            "energy_before_md_kj_mol": energy_before_md,
            "energy_after_md_kj_mol": energy_after_md,
        }
        failure_path.write_text(json.dumps(failure_payload, indent=2) + "\n", encoding="utf-8")
        diagnostics = ConjugateRelaxationDiagnostics(
            success=False,
            platform_name=platform_name,
            settings=relaxation_settings.model_dump(mode="json"),
            md_steps=relaxation_settings.md_steps,
            fixed_atom_count=len(fixed_indices),
            temporary_anchor_count=anchor_count,
            removed_barostat_count=removed_barostats,
            barostat_used=False,
            stage_a_energy_before_min_kj_mol=(
                float(energy_before_min) if math.isfinite(energy_before_min) else None
            ),
            stage_a_energy_after_min_kj_mol=(
                float(energy_after_min) if math.isfinite(energy_after_min) else None
            ),
            stage_a_force_group_energies_before_min_kj_mol=force_energies_before_min,
            stage_a_force_group_energies_after_min_kj_mol=force_energies_after_min,
            stage_a_linkage_distances_angstrom=stage_a_distances,
            linkage_pairs=tuple(linkage_pairs),
            warnings=tuple(warnings),
            error_type=type(exc).__name__,
            error_message=str(exc),
            traceback=traceback.format_exc(),
        )
        diagnostics.write_json(diagnostics_path)
        raise

    diagnostics = ConjugateRelaxationDiagnostics(
        success=True,
        stage_a_success=True,
        stage_b_success=True,
        platform_name=platform_name,
        settings=relaxation_settings.model_dump(mode="json"),
        md_steps=relaxation_settings.md_steps,
        fixed_atom_count=len(fixed_indices),
        temporary_anchor_count=anchor_count,
        removed_barostat_count=removed_barostats,
        barostat_used=False,
        stage_a_energy_before_min_kj_mol=float(energy_before_min),
        stage_a_energy_after_min_kj_mol=float(energy_after_min),
        stage_a_force_group_energies_before_min_kj_mol=force_energies_before_min,
        stage_a_force_group_energies_after_min_kj_mol=force_energies_after_min,
        stage_a_protein_rmsd_from_initial_angstrom=stage_a_protein_rmsd,
        stage_a_protein_max_displacement_from_initial_angstrom=stage_a_protein_max,
        stage_a_linkage_distances_angstrom=stage_a_distances,
        stage_b_energy_before_md_kj_mol=float(energy_before_md),
        stage_b_energy_after_md_kj_mol=float(energy_after_md),
        stage_b_protein_rmsd_from_stage_a_angstrom=stage_b_protein_rmsd,
        stage_b_protein_max_displacement_from_stage_a_angstrom=stage_b_protein_max,
        stage_b_protein_rmsd_from_initial_angstrom=initial_to_final_protein_rmsd,
        stage_b_protein_max_displacement_from_initial_angstrom=initial_to_final_protein_max,
        stage_b_linkage_distances_angstrom=distances,
        stage_b_linkage_distance_errors_angstrom=errors,
        final_relaxed_pdb_path=relaxed_pdb,
        protein_rmsd_angstrom=stage_b_protein_rmsd,
        protein_max_displacement_angstrom=stage_b_protein_max,
        linkage_distances_angstrom=distances,
        linkage_distance_errors_angstrom=errors,
        linkage_pairs=tuple(linkage_pairs),
        warnings=tuple(warnings),
    )
    diagnostics.write_json(diagnostics_path)
    return ConjugateRelaxationResult(
        success=True,
        output_dir=artifact_dir,
        diagnostics_json_path=diagnostics_path,
        relaxed_pdb_path=relaxed_pdb,
        failure_json_path=failure_path,
        platform_name=platform_name,
        fixed_atom_count=len(fixed_indices),
        temporary_anchor_count=anchor_count,
        removed_barostat_count=removed_barostats,
        energy_before_min_kj_mol=float(energy_before_min),
        energy_after_min_kj_mol=float(energy_after_min),
        energy_before_md_kj_mol=float(energy_before_md),
        energy_after_md_kj_mol=float(energy_after_md),
        md_steps=relaxation_settings.md_steps,
        stage_a_protein_rmsd_angstrom=stage_a_protein_rmsd,
        stage_a_protein_max_displacement_angstrom=stage_a_protein_max,
        protein_rmsd_angstrom=stage_b_protein_rmsd,
        protein_max_displacement_angstrom=stage_b_protein_max,
        linkage_distances_angstrom=distances,
        diagnostics=(
            "Stage A normal-mass minimization completed",
            "Stage B restrained OpenMM relaxation completed",
            *tuple(warnings),
        ),
    )
