"""Internal relaxation helpers for conjugation workflows."""

from __future__ import annotations

import json
import logging
import math
import os
import traceback
from pathlib import Path
from typing import Any

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.builders.conjugation._linkage import PabloCrosslinkRequirement, parse_pdb_atom_records
from polyzymd.builders.conjugation.pablo.ingestion import PabloIngestor
from polyzymd.builders.conjugation.pablo.parameterization import (
    build_formal_charge_smoke_template,
    create_interchange_from_pablo_topology,
)
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord
from polyzymd.config.schema import ConjugationCcdCrosslinkConfig, ConjugationCcdPabloPolicyConfig

LOGGER = logging.getLogger(__name__)

_POSITION_CONVERSION_ERRORS = (AttributeError, TypeError, ValueError)


class VacuumSmokeSettings(BaseModel):
    """Settings for the short restrained vacuum smoke run."""

    minimization_max_iterations: int = Field(100, ge=0)
    minimization_tolerance_kj_mol_nm: float = Field(10.0, gt=0)
    nvt_steps: int = Field(10, ge=0)
    temperature_kelvin: float = Field(50.0, gt=0)
    timestep_femtoseconds: float = Field(0.25, gt=0)
    friction_per_picosecond: float = Field(10.0, gt=0)
    protein_heavy_restraint_k_kj_mol_nm2: float = Field(50000.0, gt=0)
    restrain_all_heavy_atoms: bool = True
    max_position_span_nm: float = Field(50.0, gt=0)
    platform_name: str | None = None
    smoke_json_name: str = "vacuum_smoke.json"
    diagnostics_json_name: str = "restrained_smoke_diagnostics.json"
    pre_smoke_geometry_json_name: str = "pre_smoke_geometry.json"
    failure_json_name: str = "vacuum_smoke_failure.json"
    failed_pdb_name: str = "assembled_failed_smoke.pdb"
    minimized_pdb_name: str = "assembled_minimized.pdb"
    equilibrated_pdb_name: str = "assembled_equilibrated.pdb"


class GeometryPairDiagnostic(BaseModel):
    """Diagnostic for a close contact or suspicious bonded distance."""

    atom_i: int
    atom_j: int
    distance_nm: float
    distance_angstrom: float
    atom_i_identity: str | None = None
    atom_j_identity: str | None = None
    category: str


class CrosslinkBondDiagnostic(BaseModel):
    """Measured crosslink bond length for one resolved attachment."""

    attachment_id: str | None = None
    attachment_index: int | None = None
    reaction_name: str | None = None
    protein_atom: str | None = None
    modifier_atom: str | None = None
    distance_angstrom: float | None = None
    target_distance_angstrom: float | None = None
    status: str = "measured"


class SmokePhaseDiagnostics(BaseModel):
    """Geometry and energy diagnostics for one restrained smoke phase."""

    phase: str
    coordinate_span_nm: float | None = None
    coordinates_are_finite: bool
    has_nan: bool
    has_inf: bool
    potential_energy_kj_mol: float | None = None
    max_force_kj_mol_nm: float | None = None
    min_heavy_heavy_distance_nm: float | None = None
    min_h_heavy_distance_nm: float | None = None
    close_contacts: tuple[GeometryPairDiagnostic, ...] = Field(default_factory=tuple)
    bonded_distance_outliers: tuple[GeometryPairDiagnostic, ...] = Field(default_factory=tuple)
    crosslink_bonds: tuple[CrosslinkBondDiagnostic, ...] = Field(default_factory=tuple)


class RestrainedSmokeDiagnostics(BaseModel):
    """Full diagnostic report for restrained vacuum minimization and smoke MD."""

    success: bool
    first_invalid_phase: str | None = None
    platform_name: str | None = None
    restrained_atom_count: int = 0
    settings: dict[str, Any] = Field(default_factory=dict)
    phases: tuple[SmokePhaseDiagnostics, ...] = Field(default_factory=tuple)
    smoke_segments: tuple[SmokePhaseDiagnostics, ...] = Field(default_factory=tuple)
    error_type: str | None = None
    error_message: str | None = None
    traceback: str | None = None
    json_path: Path | None = None

    def write_json(self, path: Path | str) -> Path:
        """Write diagnostics as JSON and return the output path."""
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.model_copy(update={"json_path": target}).model_dump(mode="json")
        target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.json_path = target
        return target


class PreSmokeGeometryDiagnostics(BaseModel):
    """Geometry diagnostics collected immediately before OpenMM smoke relaxation."""

    atom_count: int
    coordinate_span_nm: float
    min_heavy_heavy_distance_nm: float | None = None
    min_h_heavy_distance_nm: float | None = None
    close_contacts: tuple[GeometryPairDiagnostic, ...] = Field(default_factory=tuple)
    bonded_distance_outliers: tuple[GeometryPairDiagnostic, ...] = Field(default_factory=tuple)
    crosslink_bonds: tuple[CrosslinkBondDiagnostic, ...] = Field(default_factory=tuple)
    json_path: Path | None = None

    def write_json(self, path: Path | str) -> Path:
        """Write diagnostics as JSON and return the output path."""
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.model_copy(update={"json_path": target}).model_dump(mode="json")
        target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.json_path = target
        return target


class VacuumSmokeResult(BaseModel):
    """Summary and artifacts from a restrained vacuum smoke run."""

    success: bool
    output_dir: Path
    smoke_json_path: Path
    diagnostics_json_path: Path | None = None
    pre_smoke_geometry_json_path: Path | None = None
    failure_json_path: Path | None = None
    minimized_pdb_path: Path | None = None
    equilibrated_pdb_path: Path | None = None
    failed_pdb_path: Path | None = None
    platform_name: str
    restrained_atom_count: int
    energy_before_min_kj_mol: float
    energy_after_min_kj_mol: float
    energy_before_nvt_kj_mol: float
    energy_after_nvt_kj_mol: float
    minimized_coordinate_span_nm: float | None = None
    equilibrated_coordinate_span_nm: float | None = None
    force_group_energies_before_min_kj_mol: dict[str, float] = Field(default_factory=dict)
    force_group_energies_after_min_kj_mol: dict[str, float] = Field(default_factory=dict)
    nvt_steps: int
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)


class FrozenProteinRelaxationSettings(BaseModel):
    """Settings for generic product relaxation with a frozen protein."""

    minimization_max_iterations: int = Field(0, ge=0)
    minimization_tolerance_kj_mol_nm: float = Field(10.0, gt=0)
    md_steps: int = Field(10_000, gt=0)
    temperature_kelvin: float = Field(300.0, gt=0)
    timestep_femtoseconds: float = Field(2.0, gt=0)
    friction_per_picosecond: float = Field(1.0, gt=0)
    anchor_k_kj_mol_nm2: float = Field(10_000.0, gt=0)
    fallback_bond_length_angstrom: float = Field(1.5, gt=0)
    protein_chain_ids: tuple[str, ...] = ("A",)
    max_protein_rmsd_angstrom: float = Field(0.05, ge=0)
    max_protein_displacement_angstrom: float = Field(0.25, ge=0)
    max_linkage_distance_error_angstrom: float = Field(0.35, ge=0)
    platform_name: str | None = None
    diagnostics_json_name: str = "frozen_protein_relaxation_diagnostics.json"
    relaxed_pdb_name: str = "assembled_frozen_protein_relaxed.pdb"
    minimized_pdb_name: str = "assembled_frozen_protein_minimized.pdb"
    failure_json_name: str = "frozen_protein_relaxation_failure.json"


class ProductLinkagePair(BaseModel):
    """Product topology atom-index pair for one resolved attachment linkage."""

    attachment_id: str | None = None
    attachment_index: int | None = None
    protein_atom_index: int
    modifier_atom_index: int
    protein_serial: int | None = None
    modifier_serial: int | None = None
    target_bond_length_angstrom: float
    used_fallback_target: bool = False


class FrozenProteinRelaxationDiagnostics(BaseModel):
    """Diagnostics from generic frozen-protein product relaxation."""

    success: bool
    stage_a_success: bool = False
    stage_b_success: bool = False
    platform_name: str | None = None
    settings: dict[str, Any] = Field(default_factory=dict)
    frozen_atom_count: int = 0
    temporary_anchor_count: int = 0
    removed_barostat_count: int = 0
    barostat_used: bool = False
    stage_a_energy_before_min_kj_mol: float | None = None
    stage_a_energy_after_min_kj_mol: float | None = None
    stage_a_force_group_energies_before_min_kj_mol: dict[str, float] = Field(default_factory=dict)
    stage_a_force_group_energies_after_min_kj_mol: dict[str, float] = Field(default_factory=dict)
    stage_a_protein_rmsd_from_initial_angstrom: float | None = None
    stage_a_protein_max_displacement_from_initial_angstrom: float | None = None
    stage_a_linkage_distances_angstrom: tuple[float, ...] = Field(default_factory=tuple)
    stage_b_energy_before_md_kj_mol: float | None = None
    stage_b_energy_after_md_kj_mol: float | None = None
    stage_b_protein_rmsd_from_stage_a_angstrom: float | None = None
    stage_b_protein_max_displacement_from_stage_a_angstrom: float | None = None
    stage_b_protein_rmsd_from_initial_angstrom: float | None = None
    stage_b_protein_max_displacement_from_initial_angstrom: float | None = None
    stage_b_linkage_distances_angstrom: tuple[float, ...] = Field(default_factory=tuple)
    stage_b_linkage_distance_errors_angstrom: tuple[float, ...] = Field(default_factory=tuple)
    final_relaxed_pdb_path: Path | None = None
    protein_rmsd_angstrom: float | None = None
    protein_max_displacement_angstrom: float | None = None
    linkage_distances_angstrom: tuple[float, ...] = Field(default_factory=tuple)
    linkage_distance_errors_angstrom: tuple[float, ...] = Field(default_factory=tuple)
    linkage_pairs: tuple[ProductLinkagePair, ...] = Field(default_factory=tuple)
    warnings: tuple[str, ...] = Field(default_factory=tuple)
    error_type: str | None = None
    error_message: str | None = None
    traceback: str | None = None
    json_path: Path | None = None

    def write_json(self, path: Path | str) -> Path:
        """Write diagnostics as JSON and return the output path."""
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.model_copy(update={"json_path": target}).model_dump(mode="json")
        target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.json_path = target
        return target


class FrozenProteinRelaxationResult(BaseModel):
    """Summary and artifacts from generic frozen-protein product relaxation."""

    success: bool
    output_dir: Path
    diagnostics_json_path: Path
    relaxed_pdb_path: Path | None = None
    equilibrated_pdb_path: Path | None = None
    minimized_pdb_path: Path | None = None
    failure_json_path: Path | None = None
    platform_name: str
    frozen_atom_count: int
    temporary_anchor_count: int
    removed_barostat_count: int
    energy_before_min_kj_mol: float
    energy_after_min_kj_mol: float
    energy_before_md_kj_mol: float
    energy_after_md_kj_mol: float
    md_steps: int
    stage_a_protein_rmsd_angstrom: float | None = None
    stage_a_protein_max_displacement_angstrom: float | None = None
    protein_rmsd_angstrom: float | None = None
    protein_max_displacement_angstrom: float | None = None
    linkage_distances_angstrom: tuple[float, ...] = Field(default_factory=tuple)
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)


def run_restrained_vacuum_smoke(
    interchange: Any,
    output_dir: Path | str,
    *,
    protein_heavy_atom_indices: tuple[int, ...] | None = None,
    settings: VacuumSmokeSettings | None = None,
    crosslinked_pdb_path: Path | str | None = None,
    attachment_specs: tuple[Any, ...] = (),
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
    crosslinked_pdb_path : pathlib.Path, str, or None, optional
        Product PDB used to measure resolved crosslink lengths, by default
        ``None``.
    attachment_specs : tuple of Any, optional
        Resolved attachment build specs used for crosslink-specific diagnostics,
        by default ``()``.

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

    smoke_settings = _smoke_settings_from_environment(settings or VacuumSmokeSettings())
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
    diagnostics_path = artifact_dir / smoke_settings.diagnostics_json_name
    phase_diagnostics: list[SmokePhaseDiagnostics] = []
    smoke_segment_diagnostics: list[SmokePhaseDiagnostics] = []
    first_invalid_phase: str | None = None

    current_positions = initial_positions
    energy_before_min = math.nan
    energy_after_min = math.nan
    energy_before_nvt = math.nan
    energy_after_nvt = math.nan
    minimized_pdb: Path | None = None
    equilibrated_pdb: Path | None = None
    failure_json_path: Path | None = None
    failed_pdb_path: Path | None = None
    force_energies_before: dict[str, float] = {}
    force_energies_after: dict[str, float] = {}
    pre_smoke = analyze_pre_smoke_geometry(
        topology,
        initial_positions,
        openmm_unit,
        crosslinked_pdb_path=crosslinked_pdb_path,
        attachment_specs=attachment_specs,
    )
    pre_smoke_path = artifact_dir / smoke_settings.pre_smoke_geometry_json_name
    pre_smoke.write_json(pre_smoke_path)
    phase_diagnostics.append(
        _smoke_phase_diagnostics(
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
            smoke_settings.protein_heavy_restraint_k_kj_mol_nm2,
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
            tolerance=smoke_settings.minimization_tolerance_kj_mol_nm
            * openmm_unit.kilojoule_per_mole
            / openmm_unit.nanometer,
            maxIterations=smoke_settings.minimization_max_iterations,
        )
        state_after_min = simulation_min.context.getState(
            getEnergy=True, getPositions=True, getForces=True
        )
        energy_after_min = _state_energy_kj_mol(state_after_min, openmm_unit)
        minimized_positions = state_after_min.getPositions(asNumpy=True)
        current_positions = minimized_positions
        validate_finite_energy(energy_after_min, label="energy_after_min_kj_mol")
        phase_diagnostics.append(
            _smoke_phase_diagnostics(
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
            max_span_nm=smoke_settings.max_position_span_nm,
        )
        force_energies_after = _force_group_energies(
            simulation_min.context,
            group_labels_min,
            openmm_unit,
        )

        minimized_pdb = artifact_dir / smoke_settings.minimized_pdb_name
        _write_openmm_pdb(openmm_app, topology, minimized_positions, minimized_pdb)

        if smoke_settings.nvt_steps == 0:
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

            state_before_nvt = simulation_nvt.context.getState(getEnergy=True, getPositions=True)
            energy_before_nvt = _state_energy_kj_mol(state_before_nvt, openmm_unit)
            validate_finite_energy(energy_before_nvt, label="energy_before_nvt_kj_mol")
            phase_diagnostics.append(
                _smoke_phase_diagnostics(
                    "before_md",
                    topology,
                    state_before_nvt.getPositions(asNumpy=True),
                    openmm_unit,
                    potential_energy_kj_mol=energy_before_nvt,
                    attachment_specs=attachment_specs,
                )
            )
            for step_index in range(smoke_settings.nvt_steps):
                simulation_nvt.step(1)
                segment_state = simulation_nvt.context.getState(getEnergy=True, getPositions=True)
                segment_energy = _state_energy_kj_mol(segment_state, openmm_unit)
                segment_positions = segment_state.getPositions(asNumpy=True)
                segment = _smoke_phase_diagnostics(
                    f"after_md_step_{step_index + 1}",
                    topology,
                    segment_positions,
                    openmm_unit,
                    potential_energy_kj_mol=segment_energy,
                    attachment_specs=attachment_specs,
                )
                smoke_segment_diagnostics.append(segment)
                if first_invalid_phase is None:
                    first_invalid_phase = _invalid_phase_reason(
                        segment,
                        max_span_nm=smoke_settings.max_position_span_nm,
                    )
            state_after_nvt = simulation_nvt.context.getState(getEnergy=True, getPositions=True)
            energy_after_nvt = _state_energy_kj_mol(state_after_nvt, openmm_unit)
            equilibrated_positions = state_after_nvt.getPositions(asNumpy=True)
            current_positions = equilibrated_positions
            validate_finite_energy(energy_after_nvt, label="energy_after_nvt_kj_mol")
            phase_diagnostics.append(
                _smoke_phase_diagnostics(
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
                    max_span_nm=smoke_settings.max_position_span_nm,
                )
            equilibrated_span_nm = validate_finite_positions(
                equilibrated_positions,
                openmm_unit,
                label="equilibrated_positions",
                max_span_nm=smoke_settings.max_position_span_nm,
            )
    except Exception as exc:
        if first_invalid_phase is None:
            first_invalid_phase = _first_invalid_phase(
                (*phase_diagnostics, *smoke_segment_diagnostics),
                max_span_nm=smoke_settings.max_position_span_nm,
            )
        failed_pdb_path = artifact_dir / smoke_settings.failed_pdb_name
        _safe_write_failed_pdb(openmm_app, topology, current_positions, failed_pdb_path)
        _write_restrained_smoke_diagnostics(
            diagnostics_path,
            success=False,
            settings=smoke_settings,
            platform_name=platform_name,
            restrained_atom_count=len(restrained_indices),
            phases=tuple(phase_diagnostics),
            smoke_segments=tuple(smoke_segment_diagnostics),
            first_invalid_phase=first_invalid_phase,
            exc=exc,
        )
        failure_json_path = _write_vacuum_smoke_failure(
            artifact_dir / smoke_settings.failure_json_name,
            exc=exc,
            settings=smoke_settings,
            pre_smoke=pre_smoke,
            energy_before_min=energy_before_min,
            energy_after_min=energy_after_min,
            energy_before_nvt=energy_before_nvt,
            energy_after_nvt=energy_after_nvt,
            failed_pdb_path=failed_pdb_path,
            diagnostics_path=diagnostics_path,
            first_invalid_phase=first_invalid_phase,
        )
        raise

    if smoke_settings.nvt_steps == 0:
        equilibrated_pdb = minimized_pdb
    else:
        equilibrated_pdb = artifact_dir / smoke_settings.equilibrated_pdb_name
        _write_openmm_pdb(openmm_app, topology, equilibrated_positions, equilibrated_pdb)

    result = VacuumSmokeResult(
        success=True,
        output_dir=artifact_dir,
        smoke_json_path=artifact_dir / smoke_settings.smoke_json_name,
        diagnostics_json_path=diagnostics_path,
        pre_smoke_geometry_json_path=pre_smoke_path,
        failure_json_path=failure_json_path,
        minimized_pdb_path=minimized_pdb,
        equilibrated_pdb_path=equilibrated_pdb,
        failed_pdb_path=failed_pdb_path,
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
        nvt_steps=smoke_settings.nvt_steps,
        diagnostics=(
            "Restrained vacuum OpenMM smoke completed with finite state values",
            f"Maximum minimized coordinate span was {minimized_span_nm:.3f} nm",
            (
                "NVT dynamics were skipped because nvt_steps=0"
                if smoke_settings.nvt_steps == 0
                else f"Maximum equilibrated coordinate span was {equilibrated_span_nm:.3f} nm"
            ),
        ),
    )
    result.smoke_json_path.write_text(json.dumps(result.model_dump(mode="json"), indent=2) + "\n")
    _write_restrained_smoke_diagnostics(
        diagnostics_path,
        success=True,
        settings=smoke_settings,
        platform_name=platform_name,
        restrained_atom_count=len(restrained_indices),
        phases=tuple(phase_diagnostics),
        smoke_segments=tuple(smoke_segment_diagnostics),
        first_invalid_phase=None,
        exc=None,
    )
    return result


def run_frozen_protein_product_relaxation(
    interchange: Any,
    output_dir: Path | str,
    *,
    product_pdb_path: Path | str,
    attachment_specs: tuple[Any, ...],
    assembly: Any | None = None,
    settings: FrozenProteinRelaxationSettings | None = None,
) -> FrozenProteinRelaxationResult:
    """Relax a product-state conjugate with staged minimization then frozen MD.

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
    settings : FrozenProteinRelaxationSettings or None, optional
        Relaxation settings, by default ``None``.

    Returns
    -------
    FrozenProteinRelaxationResult
        Energies, validation metrics, and artifact paths.
    """
    try:
        import openmm
        from openmm import app as openmm_app
        from openmm import unit as openmm_unit
    except ImportError as exc:
        raise RuntimeError(
            "OpenMM is required for frozen-protein conjugation relaxation. Run under a "
            "suitable PolyzyMD pixi environment."
        ) from exc

    relaxation_settings = settings or FrozenProteinRelaxationSettings()
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
            "Frozen-protein relaxation could not identify protein chain atoms to freeze "
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
    frozen_indices: tuple[int, ...] = ()
    anchor_count = 0
    minimized_pdb = artifact_dir / relaxation_settings.minimized_pdb_name
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
        _write_openmm_pdb(openmm_app, topology, minimized_positions, minimized_pdb)
        minimized_coords_nm = _positions_to_numpy(minimized_positions, openmm_unit)
        stage_a_protein_rmsd, stage_a_protein_max = _protein_displacements_angstrom(
            initial_coords_nm,
            minimized_coords_nm,
            protein_indices,
        )
        stage_a_distances = _linkage_distances_angstrom(minimized_coords_nm, linkage_pairs)

        LOGGER.info("Running Stage B frozen-protein vacuum MD for conjugate product")
        system_md = interchange.to_openmm_system()
        removed_barostats += _remove_barostats(system_md)
        frozen_indices, _original_masses = freeze_protein_chain_masses(
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
            energy_before_md, energy_after_md, final_positions = _run_frozen_product_md(
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
            restore_particle_masses(system_md, _original_masses)

        _write_openmm_pdb(openmm_app, topology, final_positions, relaxed_pdb)
        final_coords_nm = _positions_to_numpy(final_positions, openmm_unit)
        stage_b_protein_rmsd, stage_b_protein_max = _protein_displacements_angstrom(
            minimized_coords_nm,
            final_coords_nm,
            frozen_indices,
        )
        initial_to_final_protein_rmsd, initial_to_final_protein_max = (
            _protein_displacements_angstrom(
                initial_coords_nm,
                final_coords_nm,
                frozen_indices,
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
    except Exception as exc:
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
        diagnostics = FrozenProteinRelaxationDiagnostics(
            success=False,
            platform_name=platform_name,
            settings=relaxation_settings.model_dump(mode="json"),
            frozen_atom_count=len(frozen_indices),
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

    diagnostics = FrozenProteinRelaxationDiagnostics(
        success=True,
        stage_a_success=True,
        stage_b_success=True,
        platform_name=platform_name,
        settings=relaxation_settings.model_dump(mode="json"),
        frozen_atom_count=len(frozen_indices),
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
    return FrozenProteinRelaxationResult(
        success=True,
        output_dir=artifact_dir,
        diagnostics_json_path=diagnostics_path,
        relaxed_pdb_path=relaxed_pdb,
        equilibrated_pdb_path=relaxed_pdb,
        minimized_pdb_path=minimized_pdb,
        failure_json_path=failure_path,
        platform_name=platform_name,
        frozen_atom_count=len(frozen_indices),
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
            "Stage B frozen-protein vacuum MD completed",
            *tuple(warnings),
        ),
    )


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


def analyze_pre_smoke_geometry(
    topology: Any,
    positions: Any,
    unit_module: Any | None = None,
    *,
    crosslinked_pdb_path: Path | str | None = None,
    attachment_specs: tuple[Any, ...] = (),
    heavy_heavy_close_nm: float = 0.12,
    h_heavy_close_nm: float = 0.08,
    max_pairs: int = 50,
) -> PreSmokeGeometryDiagnostics:
    """Measure geometry diagnostics before restrained vacuum smoke.

    Parameters
    ----------
    topology : Any
        OpenMM topology-like object with atoms and bonds.
    positions : Any
        Coordinates in nanometers or an OpenMM-compatible quantity.
    unit_module : Any or None, optional
        OpenMM unit module used for coordinate conversion, by default ``None``.
    crosslinked_pdb_path : pathlib.Path, str, or None, optional
        Product PDB used to measure crosslink lengths, by default ``None``.
    attachment_specs : tuple of Any, optional
        Attachment specs carrying resolved plans, by default ``()``.
    heavy_heavy_close_nm : float, optional
        Heavy-heavy close-contact threshold, by default 0.12.
    h_heavy_close_nm : float, optional
        Hydrogen-heavy close-contact threshold, by default 0.08.
    max_pairs : int, optional
        Maximum number of pair diagnostics to retain per category, by default 50.

    Returns
    -------
    PreSmokeGeometryDiagnostics
        Structured geometry report.
    """
    coords = _positions_to_numpy(positions, unit_module)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise RuntimeError(f"Pre-smoke coordinates have invalid shape: {coords.shape}")
    if not np.all(np.isfinite(coords)):
        raise RuntimeError("Pre-smoke coordinates contain non-finite values")
    atom_records = tuple(topology.atoms())
    identities = tuple(_topology_atom_identity(atom) for atom in atom_records)
    heavy = tuple(_is_heavy_atom(atom) for atom in atom_records)
    bonded_pairs = _topology_bond_index_pairs(topology)
    bonded_set = {tuple(sorted(pair)) for pair in bonded_pairs}
    min_heavy_heavy, min_h_heavy, contacts = _close_contact_diagnostics(
        coords,
        heavy,
        identities,
        bonded_set=bonded_set,
        heavy_heavy_close_nm=heavy_heavy_close_nm,
        h_heavy_close_nm=h_heavy_close_nm,
        max_pairs=max_pairs,
    )
    outliers = _bonded_distance_outliers(coords, identities, bonded_pairs, max_pairs=max_pairs)
    crosslinks = _crosslink_bond_diagnostics(crosslinked_pdb_path, attachment_specs)
    return PreSmokeGeometryDiagnostics(
        atom_count=int(coords.shape[0]),
        coordinate_span_nm=float(np.max(np.ptp(coords, axis=0))) if coords.size else 0.0,
        min_heavy_heavy_distance_nm=min_heavy_heavy,
        min_h_heavy_distance_nm=min_h_heavy,
        close_contacts=tuple(contacts),
        bonded_distance_outliers=tuple(outliers),
        crosslink_bonds=tuple(crosslinks),
    )


def _smoke_phase_diagnostics(
    phase: str,
    topology: Any,
    positions: Any,
    unit_module: Any | None,
    *,
    potential_energy_kj_mol: float | None,
    max_force_kj_mol_nm: float | None = None,
    attachment_specs: tuple[Any, ...] = (),
) -> SmokePhaseDiagnostics:
    """Collect finite, span, contact, bond, and linkage diagnostics for a phase."""
    coords = _positions_to_numpy(positions, unit_module)
    has_nan = bool(np.any(np.isnan(coords))) if coords.size else False
    has_inf = bool(np.any(np.isinf(coords))) if coords.size else False
    finite = bool(np.all(np.isfinite(coords))) if coords.size else False
    span_nm = None
    min_heavy_heavy = None
    min_h_heavy = None
    contacts: list[GeometryPairDiagnostic] = []
    outliers: list[GeometryPairDiagnostic] = []
    crosslinks: list[CrosslinkBondDiagnostic] = []
    if finite and coords.ndim == 2 and coords.shape[1] == 3:
        atom_records = tuple(topology.atoms())
        identities = tuple(_topology_atom_identity(atom) for atom in atom_records)
        heavy = tuple(_is_heavy_atom(atom) for atom in atom_records)
        bonded_pairs = _topology_bond_index_pairs(topology)
        bonded_set = {tuple(sorted(pair)) for pair in bonded_pairs}
        span_nm = float(np.max(np.ptp(coords, axis=0))) if coords.size else 0.0
        min_heavy_heavy, min_h_heavy, contacts = _close_contact_diagnostics(
            coords,
            heavy,
            identities,
            bonded_set=bonded_set,
            heavy_heavy_close_nm=0.12,
            h_heavy_close_nm=0.08,
            max_pairs=50,
        )
        outliers = _bonded_distance_outliers(coords, identities, bonded_pairs, max_pairs=50)
        crosslinks = _crosslink_bond_diagnostics_from_topology(topology, coords, attachment_specs)
    return SmokePhaseDiagnostics(
        phase=phase,
        coordinate_span_nm=span_nm,
        coordinates_are_finite=finite,
        has_nan=has_nan,
        has_inf=has_inf,
        potential_energy_kj_mol=potential_energy_kj_mol,
        max_force_kj_mol_nm=max_force_kj_mol_nm,
        min_heavy_heavy_distance_nm=min_heavy_heavy,
        min_h_heavy_distance_nm=min_h_heavy,
        close_contacts=tuple(contacts),
        bonded_distance_outliers=tuple(outliers),
        crosslink_bonds=tuple(crosslinks),
    )


def _invalid_phase_reason(
    phase: SmokePhaseDiagnostics,
    *,
    max_span_nm: float,
) -> str | None:
    """Return the phase name when coordinates or energies first become invalid."""
    if not phase.coordinates_are_finite or phase.has_nan or phase.has_inf:
        return phase.phase
    if phase.coordinate_span_nm is not None and phase.coordinate_span_nm > max_span_nm:
        return phase.phase
    energy = phase.potential_energy_kj_mol
    if energy is not None and not math.isfinite(float(energy)):
        return phase.phase
    return None


def _first_invalid_phase(
    phases: tuple[SmokePhaseDiagnostics, ...],
    *,
    max_span_nm: float,
) -> str | None:
    """Return the first invalid phase name from ordered diagnostics."""
    for phase in phases:
        reason = _invalid_phase_reason(phase, max_span_nm=max_span_nm)
        if reason is not None:
            return reason
    return None


def _close_contact_diagnostics(
    coords: np.ndarray,
    heavy: tuple[bool, ...],
    identities: tuple[str, ...],
    *,
    bonded_set: set[tuple[int, int]],
    heavy_heavy_close_nm: float,
    h_heavy_close_nm: float,
    max_pairs: int,
) -> tuple[float | None, float | None, list[GeometryPairDiagnostic]]:
    """Return minimum nonbonded distances and close-contact pairs."""
    min_heavy_heavy: float | None = None
    min_h_heavy: float | None = None
    contacts: list[GeometryPairDiagnostic] = []
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            if (i, j) in bonded_set:
                continue
            distance_nm = float(np.linalg.norm(coords[i] - coords[j]))
            both_heavy = heavy[i] and heavy[j]
            one_h_one_heavy = heavy[i] != heavy[j]
            if both_heavy:
                min_heavy_heavy = (
                    distance_nm if min_heavy_heavy is None else min(min_heavy_heavy, distance_nm)
                )
            if one_h_one_heavy:
                min_h_heavy = distance_nm if min_h_heavy is None else min(min_h_heavy, distance_nm)
            category = None
            if both_heavy and distance_nm < heavy_heavy_close_nm:
                category = "heavy-heavy-close-contact"
            elif one_h_one_heavy and distance_nm < h_heavy_close_nm:
                category = "h-heavy-close-contact"
            if category is not None and len(contacts) < max_pairs:
                contacts.append(_pair_diagnostic(i, j, distance_nm, identities, category))
    contacts.sort(key=lambda item: item.distance_nm)
    return min_heavy_heavy, min_h_heavy, contacts


def _bonded_distance_outliers(
    coords: np.ndarray,
    identities: tuple[str, ...],
    bonded_pairs: tuple[tuple[int, int], ...],
    *,
    max_pairs: int,
) -> list[GeometryPairDiagnostic]:
    """Return topology bond distances outside broad covalent bounds."""
    outliers: list[GeometryPairDiagnostic] = []
    for atom_i, atom_j in bonded_pairs:
        if atom_i >= len(coords) or atom_j >= len(coords):
            continue
        distance_nm = float(np.linalg.norm(coords[atom_i] - coords[atom_j]))
        if distance_nm < 0.06 or distance_nm > 0.25:
            outliers.append(
                _pair_diagnostic(atom_i, atom_j, distance_nm, identities, "bonded-distance-outlier")
            )
    outliers.sort(key=lambda item: abs(item.distance_nm - 0.15), reverse=True)
    return outliers[:max_pairs]


def _crosslink_bond_diagnostics(
    crosslinked_pdb_path: Path | str | None,
    attachment_specs: tuple[Any, ...],
) -> list[CrosslinkBondDiagnostic]:
    """Measure resolved crosslink bonds from the emitted product PDB."""
    if crosslinked_pdb_path is None or not attachment_specs:
        return []
    try:
        atoms = parse_pdb_atom_records(Path(crosslinked_pdb_path))
    except (OSError, ValueError) as exc:
        return [CrosslinkBondDiagnostic(status=f"unavailable: {exc}")]
    diagnostics: list[CrosslinkBondDiagnostic] = []
    for spec in attachment_specs:
        plan = getattr(spec, "resolved_plan", spec)
        protein_atom = getattr(plan, "protein_link_atom", None)
        modifier_atom = getattr(plan, "modifier_link_atom", None)
        product_protein = _matching_product_atom(atoms, protein_atom, plan, role="protein")
        product_modifier = _matching_product_atom(atoms, modifier_atom, plan, role="modifier")
        distance = None
        status = "measured"
        if product_protein is None or product_modifier is None:
            status = "missing product atom"
        else:
            distance = _distance(_atom_position(product_protein), _atom_position(product_modifier))
        diagnostics.append(
            CrosslinkBondDiagnostic(
                attachment_id=getattr(spec, "attachment_id", None),
                attachment_index=getattr(spec, "attachment_index", None),
                reaction_name=getattr(spec, "reaction_name", None),
                protein_atom=_pdb_atom_identity(product_protein or protein_atom),
                modifier_atom=_pdb_atom_identity(product_modifier or modifier_atom),
                distance_angstrom=distance,
                target_distance_angstrom=getattr(plan, "target_bond_length_angstrom", None),
                status=status,
            )
        )
    return diagnostics


def _crosslink_bond_diagnostics_from_topology(
    topology: Any,
    coords_nm: np.ndarray,
    attachment_specs: tuple[Any, ...],
) -> list[CrosslinkBondDiagnostic]:
    """Measure resolved crosslink bond lengths from topology-ordered coordinates."""
    if not attachment_specs:
        return []
    topology_atoms = tuple(topology.atoms())
    diagnostics: list[CrosslinkBondDiagnostic] = []
    for spec in attachment_specs:
        plan = getattr(spec, "resolved_plan", spec)
        protein_atom = getattr(plan, "protein_link_atom", None)
        modifier_atom = getattr(plan, "modifier_link_atom", None)
        protein_index = _matching_topology_atom_index(
            topology_atoms,
            protein_atom,
            plan,
            role="protein",
        )
        modifier_index = _matching_topology_atom_index(
            topology_atoms,
            modifier_atom,
            plan,
            role="modifier",
        )
        distance = None
        status = "measured"
        if protein_index is None or modifier_index is None:
            status = "missing topology atom"
        elif protein_index >= len(coords_nm) or modifier_index >= len(coords_nm):
            status = "coordinate index out of range"
        else:
            distance = (
                float(np.linalg.norm(coords_nm[protein_index] - coords_nm[modifier_index])) * 10.0
            )
        diagnostics.append(
            CrosslinkBondDiagnostic(
                attachment_id=getattr(spec, "attachment_id", None),
                attachment_index=getattr(spec, "attachment_index", None),
                reaction_name=getattr(spec, "reaction_name", None),
                protein_atom=(
                    _topology_atom_identity(topology_atoms[protein_index])
                    if protein_index is not None and protein_index < len(topology_atoms)
                    else _pdb_atom_identity(protein_atom)
                ),
                modifier_atom=(
                    _topology_atom_identity(topology_atoms[modifier_index])
                    if modifier_index is not None and modifier_index < len(topology_atoms)
                    else _pdb_atom_identity(modifier_atom)
                ),
                distance_angstrom=distance,
                target_distance_angstrom=getattr(plan, "target_bond_length_angstrom", None),
                status=status,
            )
        )
    return diagnostics


def _matching_topology_atom_index(
    topology_atoms: tuple[Any, ...],
    source_atom: Any,
    plan: Any,
    *,
    role: str,
) -> int | None:
    """Find a topology atom index corresponding to a resolved product link atom."""
    if source_atom is None:
        return None
    target_resname = (
        getattr(plan, "protein_product_residue_name", None)
        if role == "protein"
        else getattr(plan, "modifier_product_residue_name", None)
    )
    source_chain = str(getattr(source_atom, "chain_id", "") or "").strip()
    source_number = getattr(source_atom, "residue_number", None)
    source_insertion = str(getattr(source_atom, "insertion_code", "") or "").strip()
    source_name = str(getattr(source_atom, "atom_name", "") or "").strip()
    for atom in topology_atoms:
        residue = getattr(atom, "residue", None)
        chain = getattr(residue, "chain", None)
        if source_chain and str(getattr(chain, "id", "") or "").strip() != source_chain:
            continue
        residue_id = str(getattr(residue, "id", "") or "").strip()
        if source_number is not None and residue_id != str(source_number):
            continue
        if source_insertion and not residue_id.endswith(source_insertion):
            continue
        if (
            target_resname
            and str(getattr(residue, "name", "") or "").upper() != str(target_resname).upper()
        ):
            continue
        if str(getattr(atom, "name", "") or "").strip().upper() == source_name.upper():
            return int(atom.index)
    return None


def resolve_product_linkage_pairs(
    topology: Any,
    *,
    product_pdb_path: Path | str,
    attachment_specs: tuple[Any, ...],
    assembly: Any | None = None,
    fallback_bond_length_angstrom: float = 1.5,
    warnings: list[str] | None = None,
) -> tuple[ProductLinkagePair, ...]:
    """Resolve generic product linkage atom-index pairs for relaxation anchors.

    Resolution first uses assembly ``added_conect_pairs`` serial metadata, then
    falls back to remapped atoms from resolved attachment plans. It does not rely
    on residue names, atom names, or chemistry-specific assumptions beyond the
    metadata already carried by the resolved build plan.
    """
    if not attachment_specs:
        return ()
    product_atoms = parse_pdb_atom_records(Path(product_pdb_path))
    topology_atoms = tuple(topology.atoms())
    serial_to_index = _product_serial_to_topology_index(product_atoms, topology_atoms)
    assembly_pairs = tuple(getattr(assembly, "added_conect_pairs", ()) or ())
    if not assembly_pairs:
        pair = getattr(assembly, "added_conect_pair", None)
        assembly_pairs = (pair,) if pair is not None else ()

    resolved: list[ProductLinkagePair] = []
    for plan_index, spec in enumerate(attachment_specs, start=1):
        plan = getattr(spec, "resolved_plan", spec)
        serial_pair = _serial_pair_for_attachment(
            plan,
            spec,
            plan_index=plan_index,
            product_atoms=product_atoms,
            assembly_pairs=assembly_pairs,
        )
        if serial_pair is None:
            raise RuntimeError(
                f"Could not resolve product linkage atoms for attachment {plan_index}"
            )
        protein_serial, modifier_serial = serial_pair
        if protein_serial not in serial_to_index or modifier_serial not in serial_to_index:
            raise RuntimeError(
                "Product linkage serials could not be mapped to OpenMM topology indices: "
                f"{protein_serial}, {modifier_serial}"
            )
        target = getattr(plan, "target_bond_length_angstrom", None)
        used_fallback = False
        if target is None:
            used_fallback = True
            target = fallback_bond_length_angstrom
            if warnings is not None:
                warnings.append(
                    f"Attachment {plan_index} has no target bond length; using generic "
                    f"fallback {fallback_bond_length_angstrom:.3f} A"
                )
        resolved.append(
            ProductLinkagePair(
                attachment_id=getattr(spec, "attachment_id", None),
                attachment_index=getattr(spec, "attachment_index", plan_index),
                protein_atom_index=serial_to_index[protein_serial],
                modifier_atom_index=serial_to_index[modifier_serial],
                protein_serial=protein_serial,
                modifier_serial=modifier_serial,
                target_bond_length_angstrom=float(target),
                used_fallback_target=used_fallback,
            )
        )
    return tuple(resolved)


def freeze_chain_a_masses(
    system: Any,
    topology: Any,
    openmm_unit: Any,
) -> tuple[tuple[int, ...], dict[int, Any]]:
    """Set chain A particle masses to zero and return original masses."""
    return freeze_protein_chain_masses(system, topology, openmm_unit, chain_ids=("A",))


def freeze_protein_chain_masses(
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
        Frozen particle indices and their original masses.
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


def restore_particle_masses(system: Any, original_masses: dict[int, Any]) -> None:
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


def _run_frozen_product_md(
    topology: Any,
    system: Any,
    positions: Any,
    settings: FrozenProteinRelaxationSettings,
    openmm: Any,
    openmm_app: Any,
    openmm_unit: Any,
    platform: Any,
    warnings: list[str],
) -> tuple[float, float, Any]:
    """Run finite frozen-product MD, retrying with safer timesteps if needed."""
    timestep_schedule = _md_timestep_retry_schedule(settings.timestep_femtoseconds)
    last_error: Exception | None = None
    for timestep_fs in timestep_schedule:
        try:
            integrator = openmm.LangevinMiddleIntegrator(
                settings.temperature_kelvin * openmm_unit.kelvin,
                settings.friction_per_picosecond / openmm_unit.picosecond,
                timestep_fs * openmm_unit.femtosecond,
            )
            simulation = openmm_app.Simulation(topology, system, integrator, platform)
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
                    "Frozen-protein MD retried with a smaller timestep "
                    f"({timestep_fs:.3f} fs) after instability at the requested timestep"
                )
            return energy_before, energy_after, final_positions
        except _openmm_runtime_exceptions(openmm) as exc:
            last_error = exc
            if timestep_fs == timestep_schedule[-1]:
                break
            warnings.append(
                f"Frozen-protein MD was unstable at {timestep_fs:.3f} fs; retrying with a "
                "smaller timestep"
            )
    if last_error is None:
        raise RuntimeError("Frozen-protein MD did not run")
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
    """Return requested timestep followed by conservative vacuum-MD fallbacks."""
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
    """Remove barostat-like forces from a transient vacuum relaxation system."""
    removed = 0
    for index in reversed(range(int(system.getNumForces()))):
        force = system.getForce(index)
        if "Barostat" in type(force).__name__:
            system.removeForce(index)
            removed += 1
    return removed


def _add_linkage_anchor_restraints(
    system: Any,
    pairs: tuple[ProductLinkagePair, ...],
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


def _product_serial_to_topology_index(
    product_atoms: tuple[PdbAtomRecord, ...],
    topology_atoms: tuple[Any, ...],
) -> dict[int, int]:
    """Map product PDB atom serials to OpenMM topology atom indices."""
    identity_to_serials: dict[tuple[str, str, str, str], list[int]] = {}
    for atom in product_atoms:
        if atom.serial is None:
            continue
        identity_to_serials.setdefault(_pdb_product_identity(atom), []).append(atom.serial)
    serial_to_index: dict[int, int] = {}
    consumed: dict[tuple[str, str, str, str], int] = {}
    for atom in topology_atoms:
        residue = getattr(atom, "residue", None)
        chain = getattr(residue, "chain", None)
        identity = (
            str(getattr(chain, "id", "") or "").strip(),
            str(getattr(residue, "name", "") or "").strip(),
            str(getattr(residue, "id", "") or "").strip(),
            str(getattr(atom, "name", "") or "").strip(),
        )
        serials = identity_to_serials.get(identity, [])
        offset = consumed.get(identity, 0)
        if offset < len(serials):
            serial_to_index[serials[offset]] = int(atom.index)
            consumed[identity] = offset + 1
    return serial_to_index


def _pdb_product_identity(atom: PdbAtomRecord) -> tuple[str, str, str, str]:
    """Return identity fields shared by PDB atoms and OpenMM topology atoms."""
    residue_id = f"{atom.residue_number}{atom.insertion_code.strip()}".strip()
    return (
        atom.chain_id.strip(),
        atom.residue_name.strip(),
        residue_id,
        atom.atom_name.strip(),
    )


def _serial_pair_for_attachment(
    plan: Any,
    spec: Any,
    *,
    plan_index: int,
    product_atoms: tuple[PdbAtomRecord, ...],
    assembly_pairs: tuple[Any, ...],
) -> tuple[int, int] | None:
    """Resolve product serials for one attachment from assembly or plan metadata."""
    if 1 <= plan_index <= len(assembly_pairs):
        normalized = _normalize_serial_pair(assembly_pairs[plan_index - 1])
        if normalized is not None:
            return normalized
    protein = _matching_product_atom(
        product_atoms, getattr(plan, "protein_link_atom", None), plan, role="protein"
    )
    modifier = _matching_product_atom(
        product_atoms,
        getattr(plan, "modifier_link_atom", None),
        plan,
        role="modifier",
    )
    if protein is not None and modifier is not None and protein.serial and modifier.serial:
        return int(protein.serial), int(modifier.serial)
    mappings = getattr(spec, "product_residue_mappings", {}) or {}
    if mappings:
        # Preserve generic metadata for future richer disambiguation without guessing chemistry
        return None
    return None


def _normalize_serial_pair(pair: Any) -> tuple[int, int] | None:
    """Normalize a two-item serial pair from assembly metadata."""
    try:
        first, second = pair
        return int(first), int(second)
    except (TypeError, ValueError):
        return None


def _protein_displacements_angstrom(
    initial_nm: np.ndarray,
    final_nm: np.ndarray,
    indices: tuple[int, ...],
) -> tuple[float, float]:
    """Return RMSD and maximum displacement for frozen protein atoms."""
    if not indices:
        return 0.0, 0.0
    delta_angstrom = (final_nm[list(indices)] - initial_nm[list(indices)]) * 10.0
    distances = np.linalg.norm(delta_angstrom, axis=1)
    return float(np.sqrt(np.mean(distances**2))), float(np.max(distances))


def _linkage_distances_angstrom(
    coords_nm: np.ndarray,
    pairs: tuple[ProductLinkagePair, ...],
) -> tuple[float, ...]:
    """Measure product linkage distances from topology-order coordinates."""
    return tuple(
        float(
            np.linalg.norm(coords_nm[pair.protein_atom_index] - coords_nm[pair.modifier_atom_index])
        )
        * 10.0
        for pair in pairs
    )


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


def _topology_bond_index_pairs(topology: Any) -> tuple[tuple[int, int], ...]:
    """Return atom-index pairs from an OpenMM topology-like object."""
    pairs: list[tuple[int, int]] = []
    bonds = getattr(topology, "bonds", None)
    if bonds is None:
        return ()
    for bond in bonds():
        try:
            atom_i, atom_j = bond
            pairs.append((int(atom_i.index), int(atom_j.index)))
        except (AttributeError, TypeError, ValueError):
            continue
    return tuple(pairs)


def _topology_atom_identity(atom: Any) -> str:
    """Format an OpenMM atom-like object for diagnostics."""
    residue = getattr(atom, "residue", None)
    chain = getattr(residue, "chain", None)
    chain_id = str(getattr(chain, "id", "") or "").strip()
    residue_name = str(getattr(residue, "name", "") or "").strip()
    residue_id = str(getattr(residue, "id", "") or "").strip()
    atom_name = str(getattr(atom, "name", "") or "").strip()
    atom_index = getattr(atom, "index", None)
    return f"{atom_index}:{chain_id}:{residue_name}{residue_id}:{atom_name}"


def _pair_diagnostic(
    atom_i: int,
    atom_j: int,
    distance_nm: float,
    identities: tuple[str, ...],
    category: str,
) -> GeometryPairDiagnostic:
    """Build a pair diagnostic with atom identity text."""
    return GeometryPairDiagnostic(
        atom_i=atom_i,
        atom_j=atom_j,
        distance_nm=distance_nm,
        distance_angstrom=distance_nm * 10.0,
        atom_i_identity=identities[atom_i] if atom_i < len(identities) else None,
        atom_j_identity=identities[atom_j] if atom_j < len(identities) else None,
        category=category,
    )


def _matching_product_atom(
    atoms: tuple[PdbAtomRecord, ...],
    source_atom: Any,
    plan: Any,
    *,
    role: str,
) -> PdbAtomRecord | None:
    """Find a product PDB atom corresponding to a resolved link atom."""
    if source_atom is None:
        return None
    target_resname = (
        getattr(plan, "protein_product_residue_name", None)
        if role == "protein"
        else getattr(plan, "modifier_product_residue_name", None)
    )
    source_chain = str(getattr(source_atom, "chain_id", "") or "").strip()
    source_number = getattr(source_atom, "residue_number", None)
    source_insertion = str(getattr(source_atom, "insertion_code", "") or "").strip()
    source_name = str(getattr(source_atom, "atom_name", "") or "").strip()
    for atom in atoms:
        if source_chain and atom.chain_id.strip() != source_chain:
            continue
        if source_number is not None and atom.residue_number != source_number:
            continue
        if source_insertion and atom.insertion_code.strip() != source_insertion:
            continue
        if target_resname and atom.residue_name.strip().upper() != str(target_resname).upper():
            continue
        if atom.atom_name.strip().upper() == source_name.upper():
            return atom
    return None


def _pdb_atom_identity(atom: Any | None) -> str | None:
    """Format a PDB atom-like record for diagnostics."""
    if atom is None:
        return None
    chain_id = str(getattr(atom, "chain_id", "") or "").strip()
    residue_name = str(getattr(atom, "residue_name", "") or "").strip()
    residue_number = getattr(atom, "residue_number", None)
    atom_name = str(getattr(atom, "atom_name", "") or "").strip()
    serial = getattr(atom, "serial", None)
    return f"{serial}:{chain_id}:{residue_name}{residue_number}:{atom_name}"


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
        except Exception as exc:  # noqa: BLE001 - OpenMM errors vary by platform
            LOGGER.debug("Could not collect force-group energy for %s: %s", label, exc)
    return energies


def _safe_write_failed_pdb(
    openmm_app: Any,
    topology: Any,
    positions: Any,
    output_path: Path,
) -> None:
    """Best-effort PDB writer used while handling smoke failures."""
    try:
        _write_openmm_pdb(openmm_app, topology, positions, output_path)
    except Exception as exc:  # noqa: BLE001 - failure artifact writing is best effort
        LOGGER.warning("Could not write failed smoke PDB %s: %s", output_path, exc)


def _write_vacuum_smoke_failure(
    output_path: Path,
    *,
    exc: BaseException,
    settings: VacuumSmokeSettings,
    pre_smoke: PreSmokeGeometryDiagnostics,
    energy_before_min: float,
    energy_after_min: float,
    energy_before_nvt: float,
    energy_after_nvt: float,
    failed_pdb_path: Path | None,
    diagnostics_path: Path | None = None,
    first_invalid_phase: str | None = None,
) -> Path:
    """Write structured failure diagnostics for the vacuum smoke stage."""
    payload = {
        "success": False,
        "error_type": type(exc).__name__,
        "error_message": str(exc),
        "traceback": "".join(traceback.format_exception(type(exc), exc, exc.__traceback__)),
        "settings": settings.model_dump(mode="json"),
        "energies_kj_mol": {
            "before_min": energy_before_min if math.isfinite(energy_before_min) else None,
            "after_min": energy_after_min if math.isfinite(energy_after_min) else None,
            "before_nvt": energy_before_nvt if math.isfinite(energy_before_nvt) else None,
            "after_nvt": energy_after_nvt if math.isfinite(energy_after_nvt) else None,
        },
        "pre_smoke_geometry": pre_smoke.model_dump(mode="json"),
        "failed_pdb_path": str(failed_pdb_path) if failed_pdb_path is not None else None,
        "diagnostics_json_path": str(diagnostics_path) if diagnostics_path is not None else None,
        "first_invalid_phase": first_invalid_phase,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return output_path


def _write_restrained_smoke_diagnostics(
    output_path: Path,
    *,
    success: bool,
    settings: VacuumSmokeSettings,
    platform_name: str,
    restrained_atom_count: int,
    phases: tuple[SmokePhaseDiagnostics, ...],
    smoke_segments: tuple[SmokePhaseDiagnostics, ...],
    first_invalid_phase: str | None,
    exc: BaseException | None,
) -> Path:
    """Write the detailed restrained smoke diagnostic JSON sidecar."""
    diagnostics = RestrainedSmokeDiagnostics(
        success=success,
        first_invalid_phase=first_invalid_phase,
        platform_name=platform_name,
        restrained_atom_count=restrained_atom_count,
        settings=settings.model_dump(mode="json"),
        phases=phases,
        smoke_segments=smoke_segments,
        error_type=type(exc).__name__ if exc is not None else None,
        error_message=str(exc) if exc is not None else None,
        traceback=(
            "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
            if exc is not None
            else None
        ),
    )
    return diagnostics.write_json(output_path)


def _state_energy_kj_mol(state: Any, openmm_unit: Any) -> float:
    """Return potential energy from an OpenMM state in kJ/mol."""
    energy = state.getPotentialEnergy()
    if hasattr(energy, "value_in_unit"):
        return float(energy.value_in_unit(openmm_unit.kilojoule_per_mole))
    return float(energy)


def _state_max_force_kj_mol_nm(state: Any, openmm_unit: Any) -> float | None:
    """Return the maximum force norm from an OpenMM state when available."""
    try:
        forces = state.getForces(asNumpy=True)
    except Exception as exc:  # noqa: BLE001 - force retrieval varies by OpenMM state
        LOGGER.debug("Could not collect max force diagnostic: %s", exc)
        return None
    force_array = np.asarray(
        forces.value_in_unit(openmm_unit.kilojoule_per_mole / openmm_unit.nanometer),
        dtype=float,
    )
    if force_array.size == 0 or not np.all(np.isfinite(force_array)):
        return None
    return float(np.max(np.linalg.norm(force_array, axis=1)))


def _smoke_settings_from_environment(settings: VacuumSmokeSettings) -> VacuumSmokeSettings:
    """Apply internal environment-variable diagnostic overrides to smoke settings."""
    updates: dict[str, Any] = {}
    if _truthy_environment("POLYZYMD_CONJUGATION_SMOKE_MIN_ONLY"):
        updates["nvt_steps"] = 0
    nvt_steps = os.environ.get("POLYZYMD_CONJUGATION_SMOKE_NVT_STEPS")
    if nvt_steps not in (None, ""):
        updates["nvt_steps"] = int(nvt_steps)
    timestep_fs = os.environ.get("POLYZYMD_CONJUGATION_SMOKE_TIMESTEP_FS")
    if timestep_fs not in (None, ""):
        updates["timestep_femtoseconds"] = float(timestep_fs)
    temperature_kelvin = os.environ.get("POLYZYMD_CONJUGATION_SMOKE_TEMPERATURE_K")
    if temperature_kelvin not in (None, ""):
        updates["temperature_kelvin"] = float(temperature_kelvin)
    if not updates:
        return settings
    LOGGER.info("Applying restrained smoke diagnostic environment overrides: %s", updates)
    return settings.model_copy(update=updates)


def _truthy_environment(name: str) -> bool:
    """Return whether an environment variable contains a truthy diagnostic flag."""
    value = os.environ.get(name)
    return value is not None and value.strip().lower() not in {"", "0", "false", "no", "off"}


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
    """Select a usable OpenMM platform, preferring accelerators when available."""
    requested = requested_platform or os.environ.get("POLYZYMD_CONJUGATION_SMOKE_PLATFORM")
    names = (requested,) if requested else ("CUDA", "OpenCL", "CPU", "Reference")
    errors: list[str] = []
    for name in names:
        if name is None:
            continue
        try:
            platform = openmm.Platform.getPlatformByName(name)
            _validate_platform_context(openmm, platform)
            return platform
        except Exception as exc:  # noqa: BLE001 - OpenMM platform errors vary
            errors.append(f"{name}: {exc}")
    raise RuntimeError(
        "No suitable OpenMM platform found for conjugation smoke: " + "; ".join(errors)
    )


def _validate_platform_context(openmm: Any, platform: Any) -> None:
    """Create a tiny context to confirm the registered platform is usable."""
    from openmm import unit as openmm_unit

    system = openmm.System()
    system.addParticle(39.9)
    integrator = openmm.VerletIntegrator(0.001 * openmm_unit.picoseconds)
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


__all__ = [
    "CrosslinkAtomSelector",
    "FrozenProteinRelaxationDiagnostics",
    "FrozenProteinRelaxationResult",
    "FrozenProteinRelaxationSettings",
    "LocalGeometryMetrics",
    "LocalMinimizationResult",
    "LocalMinimizationSettings",
    "analyze_crosslink_geometry",
    "analyze_pre_smoke_geometry",
    "build_product_state_pablo_policy",
    "freeze_chain_a_masses",
    "freeze_protein_chain_masses",
    "product_state_pablo_crosslink_requirement",
    "resolve_product_linkage_pairs",
    "restore_particle_masses",
    "run_frozen_protein_product_relaxation",
    "run_post_crosslink_local_minimization",
    "write_pdb_with_replaced_coordinates",
]


class CrosslinkAtomSelector(BaseModel):
    """Fixed PDB atom selector for the local crosslink geometry check."""

    serial: int | None = Field(None, ge=1)
    chain_id: str
    residue_name: str
    residue_number: int
    atom_name: str


class LocalGeometryMetrics(BaseModel):
    """Geometry and connectivity metrics around the NHS-Lys product amide."""

    nz_c047_distance_angstrom: float
    c047_o020_distance_angstrom: float
    nz_o020_distance_angstrom: float
    nz_c047_o020_angle_degrees: float
    reciprocal_nz_c047_conect: bool
    nz_o020_conect_present: bool
    passes: bool
    failures: tuple[str, ...] = Field(default_factory=tuple)


class LocalMinimizationSettings(BaseModel):
    """Settings for restrained post-crosslink OpenMM/OpenFF/Pablo minimization."""

    nz_selector: CrosslinkAtomSelector = Field(
        default_factory=lambda: CrosslinkAtomSelector(
            serial=324,
            chain_id="A",
            residue_name="LYX",
            residue_number=23,
            atom_name="NZ",
        )
    )
    c047_selector: CrosslinkAtomSelector = Field(
        default_factory=lambda: CrosslinkAtomSelector(
            serial=2883,
            chain_id="C",
            residue_name="NHX",
            residue_number=5,
            atom_name="C047",
        )
    )
    o020_selector: CrosslinkAtomSelector = Field(
        default_factory=lambda: CrosslinkAtomSelector(
            serial=2884,
            chain_id="C",
            residue_name="NHX",
            residue_number=5,
            atom_name="O020",
        )
    )
    polymer_mobile_residue_window: int = Field(1, ge=0)
    restraint_k_kj_mol_nm2: float = Field(1_000_000.0, gt=0)
    minimization_tolerance_kj_mol_nm: float = Field(1.0, gt=0)
    minimization_max_iterations: int = Field(500, ge=0)
    platform_name: str | None = None
    relaxed_pdb_name: str = "seeded_random_10mer_crosslinked_relaxed.pdb"
    result_json_name: str = "local_minimization_result.json"
    use_default_pablo_crosslink: bool = True
    use_formal_charge_templates: bool = True


class LocalMinimizationResult(BaseModel):
    """Result record for a local minimization attempt."""

    success: bool
    input_pdb: Path
    output_dir: Path
    result_json_path: Path | None = None
    relaxed_pdb_path: Path | None = None
    before_geometry: LocalGeometryMetrics
    after_geometry: LocalGeometryMetrics | None = None
    platform_name: str | None = None
    mobile_atom_count: int = 0
    restrained_atom_count: int = 0
    energy_before_min_kj_mol: float | None = None
    energy_after_min_kj_mol: float | None = None
    max_restrained_protein_displacement_angstrom: float | None = None
    mean_restrained_protein_displacement_angstrom: float | None = None
    blocker: str | None = None
    blocker_traceback: str | None = None
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)

    def save(self, path: Path | str | None = None) -> Path:
        """Write this result as JSON and return the output path."""
        output_path = (
            Path(path) if path is not None else self.output_dir / "local_minimization_result.json"
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self.result_json_path = output_path
        output_path.write_text(json.dumps(self.model_dump(mode="json"), indent=2) + "\n")
        return output_path


def analyze_crosslink_geometry(
    pdb_path: Path | str,
    *,
    settings: LocalMinimizationSettings | None = None,
) -> LocalGeometryMetrics:
    """Measure local amide geometry and required ``CONECT`` records."""
    minimization_settings = settings or LocalMinimizationSettings()
    atoms = parse_pdb_atom_records(pdb_path)
    nz_atom = _select_atom(atoms, minimization_settings.nz_selector)
    c047_atom = _select_atom(atoms, minimization_settings.c047_selector)
    o020_atom = _select_atom(atoms, minimization_settings.o020_selector)
    conect = _parse_conect_records(Path(pdb_path))

    nz_position = _atom_position(nz_atom)
    c047_position = _atom_position(c047_atom)
    o020_position = _atom_position(o020_atom)

    nz_c047 = _distance(nz_position, c047_position)
    c047_o020 = _distance(c047_position, o020_position)
    nz_o020 = _distance(nz_position, o020_position)
    angle = _angle_degrees(nz_position, c047_position, o020_position)
    reciprocal = _has_conect(conect, nz_atom.serial, c047_atom.serial) and _has_conect(
        conect, c047_atom.serial, nz_atom.serial
    )
    nz_o020_conect = _has_conect(conect, nz_atom.serial, o020_atom.serial) or _has_conect(
        conect, o020_atom.serial, nz_atom.serial
    )

    failures: list[str] = []
    if not 1.25 <= nz_c047 <= 1.65:
        failures.append(f"NZ-C047 distance {nz_c047:.3f} A outside 1.25-1.65 A")
    if not 1.15 <= c047_o020 <= 1.35:
        failures.append(f"C047-O020 distance {c047_o020:.3f} A outside 1.15-1.35 A")
    if nz_o020 < 1.8:
        failures.append(f"NZ-O020 distance {nz_o020:.3f} A below hard 1.8 A minimum")
    if not 105.0 <= angle <= 135.0:
        failures.append(f"NZ-C047-O020 angle {angle:.2f} deg outside 105-135 deg")
    if not reciprocal:
        failures.append("NZ-C047 CONECT is not reciprocal")
    if nz_o020_conect:
        failures.append("Unexpected NZ-O020 CONECT is present")

    return LocalGeometryMetrics(
        nz_c047_distance_angstrom=nz_c047,
        c047_o020_distance_angstrom=c047_o020,
        nz_o020_distance_angstrom=nz_o020,
        nz_c047_o020_angle_degrees=angle,
        reciprocal_nz_c047_conect=reciprocal,
        nz_o020_conect_present=nz_o020_conect,
        passes=not failures,
        failures=tuple(failures),
    )


def run_post_crosslink_local_minimization(
    pdb_path: Path | str,
    output_dir: Path | str | None = None,
    *,
    settings: LocalMinimizationSettings | None = None,
    pablo_policy: Any | None = None,
    pablo_crosslink_requirement: Any | None = None,
    product_state_pablo_library: Any | None = None,
    resolved_plan: Any | None = None,
) -> LocalMinimizationResult:
    """Run restrained local minimization through OpenMM/OpenFF/Pablo.

    The function never falls back to RDKit/UFF. If Pablo ingestion or OpenFF
    parameterization fails, the returned result records the blocker and traceback.
    """
    minimization_settings = settings or LocalMinimizationSettings()
    input_pdb = Path(pdb_path)
    artifact_dir = Path(output_dir) if output_dir is not None else input_pdb.parent
    artifact_dir.mkdir(parents=True, exist_ok=True)
    before_geometry = analyze_crosslink_geometry(input_pdb, settings=minimization_settings)

    try:
        import openmm
        from openmm import app as openmm_app
        from openmm import unit as openmm_unit
    except ImportError as exc:
        return _blocked_result(
            input_pdb,
            artifact_dir,
            before_geometry,
            exc,
            "OpenMM is required for local minimization; run in the build pixi environment.",
            minimization_settings,
        )

    try:
        custom_residue_library = _extract_residue_library(product_state_pablo_library)
        policy = pablo_policy
        if policy is None and custom_residue_library is None:
            policy = _product_state_pablo_policy(
                input_pdb,
                minimization_settings,
                pablo_crosslink_requirement=pablo_crosslink_requirement,
                resolved_plan=resolved_plan,
            )
        ingestion = PabloIngestor(policy).ingest_structure(
            input_pdb,
            output_dir=artifact_dir,
            residue_library=custom_residue_library,
        )
        if not ingestion.success or ingestion.topology is None:
            raise RuntimeError(_pablo_failure_message(ingestion))

        charge_templates = ()
        if minimization_settings.use_formal_charge_templates:
            charge_templates = _formal_charge_templates_from_topology(ingestion.topology)
        interchange_result = create_interchange_from_pablo_topology(
            ingestion.topology,
            charge_from_molecules=charge_templates,
        )
        interchange = interchange_result.interchange
        topology = interchange.to_openmm_topology()
        initial_positions = _openmm_positions_from_interchange(interchange, openmm_unit)
        initial_positions_array = _positions_to_numpy(initial_positions, openmm_unit)

        atoms = parse_pdb_atom_records(input_pdb)
        if len(atoms) != len(initial_positions_array):
            raise RuntimeError(
                "Pablo/OpenFF atom count does not match the input PDB order: "
                f"PDB has {len(atoms)} atoms but positions have {len(initial_positions_array)}"
            )

        mobile_indices = _mobile_atom_indices(atoms, minimization_settings)
        restrained_indices = _restrained_protein_indices(atoms, mobile_indices)
        if not restrained_indices:
            raise RuntimeError("No nonparticipating protein atoms were selected for restraints")

        system = interchange.to_openmm_system()
        _add_positional_restraints(
            system,
            initial_positions,
            restrained_indices,
            minimization_settings.restraint_k_kj_mol_nm2,
            openmm,
            openmm_unit,
        )
        simulation, platform_name = _build_simulation_with_platform_fallback(
            openmm,
            openmm_app,
            topology,
            system,
            minimization_settings,
            openmm_unit,
        )
        simulation.context.setPositions(initial_positions)

        energy_before = _state_energy_kj_mol(
            simulation.context.getState(getEnergy=True), openmm_unit
        )
        validate_finite_energy(energy_before, label="energy_before_min_kj_mol")
        openmm.LocalEnergyMinimizer.minimize(
            simulation.context,
            tolerance=minimization_settings.minimization_tolerance_kj_mol_nm
            * openmm_unit.kilojoule_per_mole
            / openmm_unit.nanometer,
            maxIterations=minimization_settings.minimization_max_iterations,
        )
        state_after = simulation.context.getState(getEnergy=True, getPositions=True)
        energy_after = _state_energy_kj_mol(state_after, openmm_unit)
        minimized_positions = state_after.getPositions(asNumpy=True)
        validate_finite_energy(energy_after, label="energy_after_min_kj_mol")
        validate_finite_positions(minimized_positions, openmm_unit, label="minimized_positions")

        minimized_angstrom = _positions_to_numpy(minimized_positions, openmm_unit) * 10.0
        relaxed_path = artifact_dir / minimization_settings.relaxed_pdb_name
        write_pdb_with_replaced_coordinates(input_pdb, minimized_angstrom, relaxed_path)
        after_geometry = analyze_crosslink_geometry(relaxed_path, settings=minimization_settings)
        max_displacement, mean_displacement = _protein_restraint_displacements_angstrom(
            initial_positions_array * 10.0,
            minimized_angstrom,
            restrained_indices,
        )

        result = LocalMinimizationResult(
            success=after_geometry.passes and (max_displacement < 0.05),
            input_pdb=input_pdb,
            output_dir=artifact_dir,
            relaxed_pdb_path=relaxed_path,
            before_geometry=before_geometry,
            after_geometry=after_geometry,
            platform_name=platform_name,
            mobile_atom_count=len(mobile_indices),
            restrained_atom_count=len(restrained_indices),
            energy_before_min_kj_mol=energy_before,
            energy_after_min_kj_mol=energy_after,
            max_restrained_protein_displacement_angstrom=max_displacement,
            mean_restrained_protein_displacement_angstrom=mean_displacement,
            diagnostics=tuple(_success_diagnostics(after_geometry, max_displacement)),
        )
    except Exception as exc:  # noqa: BLE001 - third-party chemistry errors need capture
        result = _blocked_result(
            input_pdb,
            artifact_dir,
            before_geometry,
            exc,
            "OpenMM/OpenFF/Pablo local minimization could not be completed.",
            minimization_settings,
        )

    result.save(artifact_dir / minimization_settings.result_json_name)
    return result


def _formal_charge_templates_from_topology(topology: Any) -> tuple[Any, ...]:
    """Build smoke-only charge templates for Pablo/OpenFF local minimization."""
    molecules = tuple(getattr(topology, "molecules", ()) or ())
    return tuple(build_formal_charge_smoke_template(molecule) for molecule in molecules)


def _build_simulation_with_platform_fallback(
    openmm: Any,
    openmm_app: Any,
    topology: Any,
    system: Any,
    settings: LocalMinimizationSettings,
    openmm_unit: Any,
) -> tuple[Any, str]:
    """Create an OpenMM Simulation, retrying CPU/Reference after accelerator failures."""
    names = (
        (settings.platform_name,)
        if settings.platform_name is not None
        else ("CUDA", "OpenCL", "CPU", "Reference")
    )
    errors: list[str] = []
    for name in names:
        if name is None:
            continue
        try:
            platform = openmm.Platform.getPlatformByName(name)
            integrator = openmm.VerletIntegrator(0.001 * openmm_unit.picoseconds)
            simulation = openmm_app.Simulation(topology, system, integrator, platform)
            platform_name = platform.getName() if hasattr(platform, "getName") else str(platform)
            return simulation, platform_name
        except Exception as exc:  # noqa: BLE001 - OpenMM platform/context errors vary
            errors.append(f"{name}: {exc}")
            if settings.platform_name is not None:
                break
    raise RuntimeError(
        "No suitable OpenMM platform found for local minimization: " + "; ".join(errors)
    )


def write_pdb_with_replaced_coordinates(
    template_pdb_path: Path | str,
    coordinates_angstrom: Any,
    output_path: Path | str,
) -> Path:
    """Write a PDB by replacing only ATOM/HETATM coordinate columns."""
    template_path = Path(template_pdb_path)
    output = Path(output_path)
    coords = np.asarray(coordinates_angstrom, dtype=float)
    atom_count = sum(
        1
        for line in template_path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )
    if coords.shape != (atom_count, 3):
        raise ValueError(
            "Coordinate replacement requires one XYZ triplet per PDB atom: "
            f"expected {(atom_count, 3)}, got {coords.shape}"
        )
    if not np.all(np.isfinite(coords)):
        raise ValueError("Replacement coordinates contain non-finite values")

    output.parent.mkdir(parents=True, exist_ok=True)
    coord_index = 0
    out_lines: list[str] = []
    with template_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(("ATOM", "HETATM")):
                x_coord, y_coord, z_coord = coords[coord_index]
                padded = f"{line:<80}"
                line = f"{padded[:30]}{x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}" f"{padded[54:]}"
                coord_index += 1
            out_lines.append(line)
    output.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return output


def product_state_pablo_crosslink_requirement(
    product_pdb_path: Path | str,
    *,
    settings: LocalMinimizationSettings | None = None,
    pablo_crosslink_requirement: Any | None = None,
    resolved_plan: Any | None = None,
    leaving_atoms: tuple[tuple[str, ...], tuple[str, ...]] = ((), ()),
) -> PabloCrosslinkRequirement:
    """Return a Pablo crosslink requirement for an already-modified product PDB.

    PolyzyMD removes reactant-side leaving atoms before writing the local
    minimization input PDB. Pablo should therefore be told to match the emitted
    product residues and crosslink atom names, not to remove reactant-state atoms
    such as lysine NZ hydrogens again.
    """
    minimization_settings = settings or LocalMinimizationSettings()
    source = _coerce_pablo_crosslink_requirement(
        pablo_crosslink_requirement or resolved_plan,
        default_settings=minimization_settings,
    )
    requirement = PabloCrosslinkRequirement(
        residues=source.residues,
        linking_atoms=source.linking_atoms,
        leaving_atoms=leaving_atoms,
        bond_order=source.bond_order,
    )
    _validate_product_state_crosslink_atoms(product_pdb_path, requirement)
    return requirement


def build_product_state_pablo_policy(
    product_pdb_path: Path | str,
    *,
    settings: LocalMinimizationSettings | None = None,
    pablo_crosslink_requirement: Any | None = None,
    resolved_plan: Any | None = None,
) -> ConjugationCcdPabloPolicyConfig | None:
    """Build a Pablo policy for product-state local minimization ingestion."""
    minimization_settings = settings or LocalMinimizationSettings()
    if not minimization_settings.use_default_pablo_crosslink and not (
        pablo_crosslink_requirement or resolved_plan
    ):
        return None
    requirement = product_state_pablo_crosslink_requirement(
        product_pdb_path,
        settings=minimization_settings,
        pablo_crosslink_requirement=pablo_crosslink_requirement,
        resolved_plan=resolved_plan,
    )
    return _policy_from_crosslink_requirement(requirement)


def _blocked_result(
    input_pdb: Path,
    output_dir: Path,
    before_geometry: LocalGeometryMetrics,
    exc: BaseException,
    message: str,
    settings: LocalMinimizationSettings,
) -> LocalMinimizationResult:
    """Build a failed result preserving the full exception traceback."""
    return LocalMinimizationResult(
        success=False,
        input_pdb=input_pdb,
        output_dir=output_dir,
        before_geometry=before_geometry,
        blocker=f"{message} {type(exc).__name__}: {exc}",
        blocker_traceback="".join(traceback.format_exception(type(exc), exc, exc.__traceback__)),
        diagnostics=(
            "No RDKit/UFF fallback was attempted; minimization requires the OpenMM/OpenFF/Pablo route.",
            "The input PDB metadata and CONECT records were preserved; only geometry was measured.",
            "Product-state Pablo crosslinks use empty leaving atom groups only when no custom product library is supplied.",
            f"Pablo crosslink policy enabled: {settings.use_default_pablo_crosslink}",
        ),
    )


def _extract_residue_library(product_state_pablo_library: Any | None) -> Any | None:
    """Return a Pablo cache from a product-state library wrapper or raw cache."""
    if product_state_pablo_library is None:
        return None
    return getattr(product_state_pablo_library, "residue_library", product_state_pablo_library)


def _product_state_pablo_policy(
    product_pdb_path: Path | str,
    settings: LocalMinimizationSettings,
    *,
    pablo_crosslink_requirement: Any | None = None,
    resolved_plan: Any | None = None,
) -> ConjugationCcdPabloPolicyConfig | None:
    """Return the Pablo policy used by the POC local minimizer."""
    return build_product_state_pablo_policy(
        product_pdb_path,
        settings=settings,
        pablo_crosslink_requirement=pablo_crosslink_requirement,
        resolved_plan=resolved_plan,
    )


def _policy_from_crosslink_requirement(
    requirement: PabloCrosslinkRequirement,
) -> ConjugationCcdPabloPolicyConfig:
    """Build a Pablo policy from a normalized crosslink requirement."""
    return ConjugationCcdPabloPolicyConfig(
        crosslinks=[
            ConjugationCcdCrosslinkConfig(
                residues=requirement.residues,
                linking_atoms=requirement.linking_atoms,
                leaving_atoms=requirement.leaving_atoms,
                bond_order=requirement.bond_order,
            )
        ]
    )


def _coerce_pablo_crosslink_requirement(
    value: Any | None,
    *,
    default_settings: LocalMinimizationSettings,
) -> PabloCrosslinkRequirement:
    """Normalize resolved plans, requirements, dicts, or objects to a requirement."""
    if value is None:
        return PabloCrosslinkRequirement(
            residues=(
                default_settings.nz_selector.residue_name,
                default_settings.c047_selector.residue_name,
            ),
            linking_atoms=(
                default_settings.nz_selector.atom_name,
                default_settings.c047_selector.atom_name,
            ),
            leaving_atoms=((), ()),
            bond_order=1,
        )
    if isinstance(value, PabloCrosslinkRequirement):
        return value
    plan_requirement = getattr(value, "pablo_crosslink_requirement", None)
    if plan_requirement is not None:
        return _coerce_pablo_crosslink_requirement(
            plan_requirement,
            default_settings=default_settings,
        )
    if isinstance(value, dict):
        return PabloCrosslinkRequirement.model_validate(value)
    return PabloCrosslinkRequirement(
        residues=value.residues,
        linking_atoms=value.linking_atoms,
        leaving_atoms=getattr(value, "leaving_atoms", ((), ())),
        bond_order=getattr(value, "bond_order", 1),
    )


def _validate_product_state_crosslink_atoms(
    product_pdb_path: Path | str,
    requirement: PabloCrosslinkRequirement,
) -> None:
    """Validate that product-state residue and crosslink atom names exist in the PDB."""
    atoms = parse_pdb_atom_records(product_pdb_path)
    missing: list[str] = []
    for residue_name, atom_name in zip(requirement.residues, requirement.linking_atoms):
        if not any(
            atom.residue_name.upper() == residue_name and atom.atom_name.upper() == atom_name
            for atom in atoms
        ):
            missing.append(f"{residue_name}:{atom_name}")
    if missing:
        raise ValueError(
            "Product-state Pablo crosslink atoms were not found in the emitted PDB: "
            + ", ".join(missing)
        )


def _pablo_failure_message(ingestion: Any) -> str:
    """Return a compact Pablo failure message from an ingestion result."""
    errors = [
        diagnostic
        for diagnostic in getattr(ingestion, "diagnostics", ())
        if getattr(getattr(diagnostic, "severity", None), "value", None) == "error"
    ]
    if not errors:
        return "Pablo ingestion failed without a detailed diagnostic"
    details = getattr(errors[-1], "details", {})
    return str(details.get("error") or getattr(errors[-1], "message", "Pablo ingestion failed"))


def _success_diagnostics(
    geometry: LocalGeometryMetrics,
    max_displacement_angstrom: float,
) -> list[str]:
    """Return human-readable success/failure diagnostics after minimization."""
    diagnostics = [
        "OpenMM/OpenFF/Pablo local minimization completed without RDKit/UFF fallback",
        f"Maximum restrained protein displacement was {max_displacement_angstrom:.4f} A",
    ]
    if geometry.passes:
        diagnostics.append("Relaxed geometry passed local amide validation")
    else:
        diagnostics.extend(geometry.failures)
    return diagnostics


def _select_atom(
    atoms: tuple[PdbAtomRecord, ...], selector: CrosslinkAtomSelector
) -> PdbAtomRecord:
    """Return exactly one atom matching a fixed selector."""
    matches = [atom for atom in atoms if _selector_matches(atom, selector)]
    if len(matches) != 1:
        raise ValueError(
            "Expected exactly one atom for selector "
            f"{selector.model_dump(mode='json')}, found {len(matches)}"
        )
    return matches[0]


def _selector_matches(atom: PdbAtomRecord, selector: CrosslinkAtomSelector) -> bool:
    """Return whether an atom matches a crosslink selector."""
    return (
        (selector.serial is None or atom.serial == selector.serial)
        and atom.chain_id.upper() == selector.chain_id.upper()
        and atom.residue_name.upper() == selector.residue_name.upper()
        and atom.residue_number == selector.residue_number
        and atom.atom_name.upper() == selector.atom_name.upper()
    )


def _atom_position(atom: PdbAtomRecord) -> np.ndarray:
    """Return an atom position in Angstrom."""
    return np.asarray([atom.x, atom.y, atom.z], dtype=float)


def _distance(first: np.ndarray, second: np.ndarray) -> float:
    """Return a Cartesian distance in Angstrom."""
    return float(np.linalg.norm(first - second))


def _angle_degrees(first: np.ndarray, vertex: np.ndarray, third: np.ndarray) -> float:
    """Return the angle first-vertex-third in degrees."""
    vector_a = first - vertex
    vector_b = third - vertex
    norm_product = float(np.linalg.norm(vector_a) * np.linalg.norm(vector_b))
    if norm_product == 0.0:
        raise ValueError("Cannot compute angle with coincident atoms")
    cosine = float(np.dot(vector_a, vector_b) / norm_product)
    return float(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))


def _parse_conect_records(path: Path) -> dict[int, set[int]]:
    """Parse PDB ``CONECT`` records into a serial adjacency map."""
    conect: dict[int, set[int]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("CONECT"):
                continue
            try:
                source = int(line[6:11])
            except ValueError:
                continue
            targets = conect.setdefault(source, set())
            for start in range(11, len(line), 5):
                field = line[start : start + 5].strip()
                if not field:
                    continue
                try:
                    targets.add(int(field))
                except ValueError:
                    continue
    return conect


def _has_conect(conect: dict[int, set[int]], source: int | None, target: int | None) -> bool:
    """Return whether a directed ``CONECT`` entry exists."""
    if source is None or target is None:
        return False
    return target in conect.get(source, set())


def _mobile_atom_indices(
    atoms: tuple[PdbAtomRecord, ...],
    settings: LocalMinimizationSettings,
) -> tuple[int, ...]:
    """Select lysine side chain and local polymer atoms allowed to move."""
    mobile: set[int] = set()
    lysine_sidechain_names = {
        "CB",
        "CG",
        "CD",
        "CE",
        "NZ",
        "H02",
        "H03",
        "H04",
        "H05",
        "H06",
        "H07",
        "H08",
        "H09",
        "H10",
    }
    polymer_residue_min = (
        settings.c047_selector.residue_number - settings.polymer_mobile_residue_window
    )
    polymer_residue_max = (
        settings.c047_selector.residue_number + settings.polymer_mobile_residue_window
    )
    for atom in atoms:
        if atom.atom_index is None:
            continue
        if (
            atom.chain_id.upper() == settings.nz_selector.chain_id.upper()
            and atom.residue_name.upper() == settings.nz_selector.residue_name.upper()
            and atom.residue_number == settings.nz_selector.residue_number
            and atom.atom_name.upper() in lysine_sidechain_names
        ):
            mobile.add(atom.atom_index)
        if (
            atom.chain_id.upper() == settings.c047_selector.chain_id.upper()
            and polymer_residue_min <= atom.residue_number <= polymer_residue_max
        ):
            mobile.add(atom.atom_index)
    return tuple(sorted(mobile))


def _restrained_protein_indices(
    atoms: tuple[PdbAtomRecord, ...],
    mobile_indices: tuple[int, ...],
) -> tuple[int, ...]:
    """Select all non-mobile protein atoms for strong positional restraints."""
    mobile = set(mobile_indices)
    indices = [
        atom.atom_index
        for atom in atoms
        if atom.atom_index is not None
        and atom.chain_id.upper() == "A"
        and atom.atom_index not in mobile
    ]
    return tuple(sorted(indices))


def _protein_restraint_displacements_angstrom(
    initial_angstrom: np.ndarray,
    final_angstrom: np.ndarray,
    restrained_indices: tuple[int, ...],
) -> tuple[float, float]:
    """Return max and mean displacement for restrained protein atoms."""
    if not restrained_indices:
        return 0.0, 0.0
    displacement = np.linalg.norm(
        final_angstrom[list(restrained_indices)] - initial_angstrom[list(restrained_indices)],
        axis=1,
    )
    return float(np.max(displacement)), float(np.mean(displacement))
