"""Data models for conjugate relaxation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field


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
