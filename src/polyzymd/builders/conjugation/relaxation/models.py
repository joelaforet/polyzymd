"""Data models for conjugate relaxation and OpenMM validation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field


class OpenMMValidationSettings(BaseModel):
    """Settings for the short restrained OpenMM validation run."""

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
    validation_json_name: str = "openmm_validation.json"
    product_geometry_json_name: str = "product_geometry_diagnostics.json"


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


class OpenMMValidationPhaseDiagnostics(BaseModel):
    """Geometry and energy diagnostics for one restrained validation phase."""

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


class ProductGeometryDiagnostics(BaseModel):
    """Geometry diagnostics collected immediately before OpenMM validation relaxation."""

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


class OpenMMValidationResult(BaseModel):
    """Summary and artifacts from a restrained OpenMM validation run."""

    success: bool
    output_dir: Path
    validation_json_path: Path
    diagnostics_json_path: Path | None = None
    product_geometry_json_path: Path | None = None
    failure_json_path: Path | None = None
    platform_name: str | None = None
    restrained_atom_count: int = 0
    energy_before_min_kj_mol: float | None = None
    energy_after_min_kj_mol: float | None = None
    energy_before_nvt_kj_mol: float | None = None
    energy_after_nvt_kj_mol: float | None = None
    minimized_coordinate_span_nm: float | None = None
    equilibrated_coordinate_span_nm: float | None = None
    force_group_energies_before_min_kj_mol: dict[str, float] = Field(default_factory=dict)
    force_group_energies_after_min_kj_mol: dict[str, float] = Field(default_factory=dict)
    nvt_steps: int = 0
    settings: dict[str, Any] = Field(default_factory=dict)
    phases: tuple[OpenMMValidationPhaseDiagnostics, ...] = Field(default_factory=tuple)
    validation_segments: tuple[OpenMMValidationPhaseDiagnostics, ...] = Field(default_factory=tuple)
    first_invalid_phase: str | None = None
    error_type: str | None = None
    error_message: str | None = None
    traceback: str | None = None
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)

    def write_json(self, path: Path | str) -> Path:
        """Write the canonical validation result as JSON.

        Parameters
        ----------
        path : pathlib.Path or str
            Destination path.

        Returns
        -------
        pathlib.Path
            Written path.
        """
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.model_copy(update={"validation_json_path": target}).model_dump(mode="json")
        target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.validation_json_path = target
        return target


class ConjugateRelaxationSettings(BaseModel):
    """Settings for conjugate relaxation."""

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
    diagnostics_json_name: str = "conjugate_relaxation.json"
    relaxed_pdb_name: str = "conjugate_relaxed.pdb"
    failure_json_name: str = "conjugate_relaxation_failure.json"


class ProductLinkage(BaseModel):
    """Product topology atom-index pair for one resolved attachment linkage."""

    attachment_id: str | None = None
    attachment_index: int | None = None
    protein_atom_index: int
    modifier_atom_index: int
    protein_serial: int | None = None
    modifier_serial: int | None = None
    target_bond_length_angstrom: float
    used_fallback_target: bool = False


class ConjugateRelaxationDiagnostics(BaseModel):
    """Diagnostics from generic conjugate product relaxation."""

    success: bool
    stage_a_success: bool = False
    stage_b_success: bool = False
    platform_name: str | None = None
    settings: dict[str, Any] = Field(default_factory=dict)
    md_steps: int | None = None
    fixed_atom_count: int = 0
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
    linkage_pairs: tuple[ProductLinkage, ...] = Field(default_factory=tuple)
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


class ConjugateRelaxationResult(BaseModel):
    """Summary and artifacts from generic conjugate product relaxation."""

    success: bool
    output_dir: Path
    diagnostics_json_path: Path
    relaxed_pdb_path: Path | None = None
    failure_json_path: Path | None = None
    platform_name: str
    fixed_atom_count: int
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
