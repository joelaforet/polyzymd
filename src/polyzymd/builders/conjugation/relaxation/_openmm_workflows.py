"""Private OpenMM workflow execution for conjugate relaxation."""

from __future__ import annotations

import json
import logging
import math
import traceback
from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation.relaxation._diagnostics import (
    validate_finite_energy,
    validate_finite_positions,
)
from polyzymd.builders.conjugation.relaxation._linkages import (
    _linkage_distances_angstrom,
    resolve_product_linkage_pairs,
)
from polyzymd.builders.conjugation.relaxation._openmm_system import (
    _add_linkage_anchor_restraints,
    _assign_force_groups,
    _create_relaxation_verlet_integrator,
    _create_simulation,
    _force_group_energies,
    _force_group_labels,
    _freeze_protein_chain_masses,
    _openmm_positions_from_interchange_or_pdb,
    _protein_chain_indices,
    _protein_displacements_angstrom,
    _remove_barostats,
    _restore_particle_masses,
    _run_fixed_product_md,
    _select_platform,
    _state_energy_kj_mol,
    _write_relaxed_product_pdb,
)
from polyzymd.builders.conjugation.relaxation.geometry import positions_to_numpy
from polyzymd.builders.conjugation.relaxation.models import (
    ConjugateRelaxationDiagnostics,
    ConjugateRelaxationResult,
    ConjugateRelaxationSettings,
)

LOGGER = logging.getLogger(__name__)


def _openmm_stage_errors(openmm: Any) -> tuple[type[BaseException], ...]:
    """Return expected exception types raised by OpenMM workflow stages."""
    openmm_error = getattr(openmm, "OpenMMException", RuntimeError)
    return (RuntimeError, ValueError, ArithmeticError, openmm_error)


def _finite_or_none(value: float) -> float | None:
    """Return a finite diagnostic value or ``None`` for JSON-safe output.

    Parameters
    ----------
    value : float
        Numeric diagnostic value that may be non-finite.

    Returns
    -------
    float or None
        ``float(value)`` when finite, otherwise ``None``.
    """
    return float(value) if math.isfinite(value) else None


def _build_relaxation_diagnostics(
    *,
    success: bool,
    platform_name: str,
    relaxation_settings: ConjugateRelaxationSettings,
    fixed_indices: tuple[int, ...],
    anchor_count: int,
    removed_barostats: int,
    energy_before_min: float,
    energy_after_min: float,
    energy_before_md: float,
    energy_after_md: float,
    force_energies_before_min: dict[str, float],
    force_energies_after_min: dict[str, float],
    stage_a_protein_rmsd: float,
    stage_a_protein_max: float,
    stage_a_distances: tuple[float, ...],
    stage_b_protein_rmsd: float,
    stage_b_protein_max: float,
    initial_to_final_protein_rmsd: float,
    initial_to_final_protein_max: float,
    distances: tuple[float, ...],
    errors: tuple[float, ...],
    relaxed_pdb: Path,
    linkage_pairs: tuple[Any, ...],
    warnings: list[str],
    error_type: str | None = None,
    error_message: str | None = None,
    error_traceback: str | None = None,
) -> ConjugateRelaxationDiagnostics:
    """Build the relaxation diagnostics payload without changing schema.

    Parameters
    ----------
    success : bool
        Whether both relaxation stages completed successfully.
    platform_name : str
        Name of the OpenMM platform used for the workflow.
    relaxation_settings : ConjugateRelaxationSettings
        Settings used to run relaxation.
    fixed_indices : tuple of int
        Fixed protein atom indices from Stage B.
    anchor_count : int
        Number of temporary linkage anchor restraints.
    removed_barostats : int
        Number of removed barostat forces.
    energy_before_min, energy_after_min, energy_before_md, energy_after_md : float
        Stage energy diagnostics in kJ/mol.
    force_energies_before_min, force_energies_after_min : dict of str to float
        Stage A force-group energy diagnostics.
    stage_a_protein_rmsd, stage_a_protein_max : float
        Stage A protein displacement diagnostics.
    stage_a_distances : tuple of float
        Stage A linkage distances.
    stage_b_protein_rmsd, stage_b_protein_max : float
        Stage B protein displacement diagnostics relative to Stage A.
    initial_to_final_protein_rmsd, initial_to_final_protein_max : float
        Stage B protein displacement diagnostics relative to initial coordinates.
    distances : tuple of float
        Final linkage distances.
    errors : tuple of float
        Final linkage distance errors.
    relaxed_pdb : pathlib.Path
        Final relaxed PDB path.
    linkage_pairs : tuple of Any
        Resolved product linkage pairs.
    warnings : list of str
        Warning messages accumulated during relaxation.
    error_type, error_message, error_traceback : str or None, optional
        Failure diagnostics, by default ``None``.

    Returns
    -------
    ConjugateRelaxationDiagnostics
        Diagnostics model matching the existing JSON schema.
    """
    if success:
        stage_a_energy_before = float(energy_before_min)
        stage_a_energy_after = float(energy_after_min)
        stage_b_energy_before = float(energy_before_md)
        stage_b_energy_after = float(energy_after_md)
        stage_a_rmsd = stage_a_protein_rmsd
        stage_a_max = stage_a_protein_max
        stage_b_rmsd = stage_b_protein_rmsd
        stage_b_max = stage_b_protein_max
        initial_to_final_rmsd = initial_to_final_protein_rmsd
        initial_to_final_max = initial_to_final_protein_max
        final_relaxed_pdb_path = relaxed_pdb
    else:
        stage_a_energy_before = _finite_or_none(energy_before_min)
        stage_a_energy_after = _finite_or_none(energy_after_min)
        stage_b_energy_before = _finite_or_none(energy_before_md)
        stage_b_energy_after = _finite_or_none(energy_after_md)
        stage_a_rmsd = _finite_or_none(stage_a_protein_rmsd)
        stage_a_max = _finite_or_none(stage_a_protein_max)
        stage_b_rmsd = _finite_or_none(stage_b_protein_rmsd)
        stage_b_max = _finite_or_none(stage_b_protein_max)
        initial_to_final_rmsd = _finite_or_none(initial_to_final_protein_rmsd)
        initial_to_final_max = _finite_or_none(initial_to_final_protein_max)
        final_relaxed_pdb_path = relaxed_pdb if relaxed_pdb.exists() else None

    return ConjugateRelaxationDiagnostics(
        success=success,
        stage_a_success=success or math.isfinite(energy_after_min),
        stage_b_success=success,
        platform_name=platform_name,
        settings=relaxation_settings.model_dump(mode="json"),
        md_steps=relaxation_settings.md_steps,
        fixed_atom_count=len(fixed_indices),
        temporary_anchor_count=anchor_count,
        removed_barostat_count=removed_barostats,
        barostat_used=False,
        stage_a_energy_before_min_kj_mol=stage_a_energy_before,
        stage_a_energy_after_min_kj_mol=stage_a_energy_after,
        stage_a_force_group_energies_before_min_kj_mol=force_energies_before_min,
        stage_a_force_group_energies_after_min_kj_mol=force_energies_after_min,
        stage_a_protein_rmsd_from_initial_angstrom=stage_a_rmsd,
        stage_a_protein_max_displacement_from_initial_angstrom=stage_a_max,
        stage_a_linkage_distances_angstrom=stage_a_distances,
        stage_b_energy_before_md_kj_mol=stage_b_energy_before,
        stage_b_energy_after_md_kj_mol=stage_b_energy_after,
        stage_b_protein_rmsd_from_stage_a_angstrom=stage_b_rmsd,
        stage_b_protein_max_displacement_from_stage_a_angstrom=stage_b_max,
        stage_b_protein_rmsd_from_initial_angstrom=initial_to_final_rmsd,
        stage_b_protein_max_displacement_from_initial_angstrom=initial_to_final_max,
        stage_b_linkage_distances_angstrom=distances,
        stage_b_linkage_distance_errors_angstrom=errors,
        final_relaxed_pdb_path=final_relaxed_pdb_path,
        protein_rmsd_angstrom=stage_b_rmsd,
        protein_max_displacement_angstrom=stage_b_max,
        linkage_distances_angstrom=distances,
        linkage_distance_errors_angstrom=errors,
        linkage_pairs=tuple(linkage_pairs),
        warnings=tuple(warnings),
        error_type=error_type,
        error_message=error_message,
        traceback=error_traceback,
    )


def _stage_b_tolerance_violations(
    *,
    stage_b_protein_rmsd: float,
    stage_b_protein_max: float,
    linkage_errors: tuple[float, ...],
    settings: ConjugateRelaxationSettings,
) -> tuple[str, ...]:
    """Return Stage B evidence failures that must block production builds.

    Parameters
    ----------
    stage_b_protein_rmsd : float
        Protein RMSD relative to the Stage A minimized coordinates.
    stage_b_protein_max : float
        Maximum protein displacement relative to the Stage A minimized coordinates.
    linkage_errors : tuple of float
        Linkage distance errors relative to target bond lengths.
    settings : ConjugateRelaxationSettings
        Relaxation tolerance settings used for the run.

    Returns
    -------
    tuple of str
        Human-readable tolerance violations.
    """
    violations: list[str] = []
    if not math.isfinite(stage_b_protein_rmsd):
        violations.append("Stage B protein RMSD relative to Stage A is non-finite")
    elif stage_b_protein_rmsd > settings.max_protein_rmsd_angstrom:
        violations.append(
            f"Stage B protein RMSD {stage_b_protein_rmsd:.4f} A relative to Stage A exceeds "
            f"{settings.max_protein_rmsd_angstrom:.4f} A"
        )
    if not math.isfinite(stage_b_protein_max):
        violations.append("Stage B protein max displacement relative to Stage A is non-finite")
    elif stage_b_protein_max > settings.max_protein_displacement_angstrom:
        violations.append(
            f"Stage B protein max displacement {stage_b_protein_max:.4f} A relative to Stage A "
            f"exceeds {settings.max_protein_displacement_angstrom:.4f} A"
        )
    nonfinite_errors = [error for error in linkage_errors if not math.isfinite(error)]
    if nonfinite_errors:
        violations.append("One or more linkage distance errors is non-finite")
    elif any(error > settings.max_linkage_distance_error_angstrom for error in linkage_errors):
        max_error = max(linkage_errors)
        violations.append(
            f"One or more linkage distances deviates from its target; max error "
            f"{max_error:.4f} A exceeds {settings.max_linkage_distance_error_angstrom:.4f} A"
        )
    return tuple(violations)


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
    initial_coords_nm = positions_to_numpy(initial_positions, openmm_unit)
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
        integrator_min = _create_relaxation_verlet_integrator(openmm, openmm_unit)
        simulation = _create_simulation(topology, system_min, integrator_min, openmm_app, platform)
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
        minimized_coords_nm = positions_to_numpy(minimized_positions, openmm_unit)
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

        _write_relaxed_product_pdb(
            product_pdb_path,
            final_positions,
            relaxed_pdb,
            openmm_unit,
            expected_atom_count=len(tuple(topology.atoms())),
        )
        final_coords_nm = positions_to_numpy(final_positions, openmm_unit)
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
        violations = _stage_b_tolerance_violations(
            stage_b_protein_rmsd=stage_b_protein_rmsd,
            stage_b_protein_max=stage_b_protein_max,
            linkage_errors=errors,
            settings=relaxation_settings,
        )
        warnings.extend(violations)
        if violations:
            violation_text = "; ".join(violations)
            raise RuntimeError(f"Conjugate relaxation Stage B evidence failed: {violation_text}")
    except _openmm_stage_errors(openmm) as exc:
        failure_path = artifact_dir / relaxation_settings.failure_json_name
        exc_traceback = traceback.format_exc()
        failure_payload = {
            "success": False,
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "traceback": exc_traceback,
            "energy_before_min_kj_mol": energy_before_min,
            "energy_after_min_kj_mol": energy_after_min,
            "energy_before_md_kj_mol": energy_before_md,
            "energy_after_md_kj_mol": energy_after_md,
        }
        failure_path.write_text(json.dumps(failure_payload, indent=2) + "\n", encoding="utf-8")
        diagnostics = _build_relaxation_diagnostics(
            success=False,
            platform_name=platform_name,
            relaxation_settings=relaxation_settings,
            fixed_indices=fixed_indices,
            anchor_count=anchor_count,
            removed_barostats=removed_barostats,
            energy_before_min=energy_before_min,
            energy_after_min=energy_after_min,
            energy_before_md=energy_before_md,
            energy_after_md=energy_after_md,
            force_energies_before_min=force_energies_before_min,
            force_energies_after_min=force_energies_after_min,
            stage_a_protein_rmsd=stage_a_protein_rmsd,
            stage_a_protein_max=stage_a_protein_max,
            stage_a_distances=stage_a_distances,
            stage_b_protein_rmsd=stage_b_protein_rmsd,
            stage_b_protein_max=stage_b_protein_max,
            initial_to_final_protein_rmsd=initial_to_final_protein_rmsd,
            initial_to_final_protein_max=initial_to_final_protein_max,
            distances=distances,
            errors=errors,
            relaxed_pdb=relaxed_pdb,
            linkage_pairs=tuple(linkage_pairs),
            warnings=warnings,
            error_type=type(exc).__name__,
            error_message=str(exc),
            error_traceback=exc_traceback,
        )
        diagnostics.write_json(diagnostics_path)
        raise

    diagnostics = _build_relaxation_diagnostics(
        success=True,
        platform_name=platform_name,
        relaxation_settings=relaxation_settings,
        fixed_indices=fixed_indices,
        anchor_count=anchor_count,
        removed_barostats=removed_barostats,
        energy_before_min=energy_before_min,
        energy_after_min=energy_after_min,
        energy_before_md=energy_before_md,
        energy_after_md=energy_after_md,
        force_energies_before_min=force_energies_before_min,
        force_energies_after_min=force_energies_after_min,
        stage_a_protein_rmsd=stage_a_protein_rmsd,
        stage_a_protein_max=stage_a_protein_max,
        stage_a_distances=stage_a_distances,
        stage_b_protein_rmsd=stage_b_protein_rmsd,
        stage_b_protein_max=stage_b_protein_max,
        initial_to_final_protein_rmsd=initial_to_final_protein_rmsd,
        initial_to_final_protein_max=initial_to_final_protein_max,
        distances=distances,
        errors=errors,
        relaxed_pdb=relaxed_pdb,
        linkage_pairs=tuple(linkage_pairs),
        warnings=warnings,
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
