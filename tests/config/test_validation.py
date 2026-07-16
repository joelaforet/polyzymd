"""Tests for runtime configuration reference validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.config.schema import SimulationConfig
from polyzymd.config.validation import collect_reference_warnings


def _minimal_config_data(tmp_path: Path) -> dict[str, object]:
    """Return minimal simulation config data for reference-validation tests."""

    return {
        "name": "reference_validation",
        "engine": "openmm",
        "enzyme": {"name": "Enz", "pdb_path": tmp_path / "missing.pdb"},
        "thermodynamics": {"temperature": 300.0},
        "simulation_phases": {
            "equilibration_stages": [
                {
                    "name": "eq",
                    "duration": 0.1,
                    "temperature": 300.0,
                    "ensemble": "NVT",
                }
            ],
            "production": {
                "ensemble": "NPT",
                "duration": 1.0,
                "samples": 10,
                "report_interval": 50000,
                "checkpoint_interval": 60.0,
            },
        },
    }


def test_collect_reference_warnings_reports_missing_enzyme_pdb(tmp_path: Path) -> None:
    """Missing enzyme PDB paths should be warnings outside schema validation."""

    config = SimulationConfig(**_minimal_config_data(tmp_path))

    warnings = collect_reference_warnings(config)

    assert any("Missing enzyme PDB" in warning for warning in warnings)


def test_collect_reference_warnings_accepts_existing_pdb_and_sdf(tmp_path: Path) -> None:
    """Existing enzyme and substrate structures should not emit warnings."""

    pdb_path = tmp_path / "enzyme.pdb"
    sdf_path = tmp_path / "substrate.sdf"
    pdb_path.write_text("HEADER test\n", encoding="utf-8")
    sdf_path.write_text("substrate\n$$$$\n", encoding="utf-8")
    data = _minimal_config_data(tmp_path)
    data["enzyme"] = {"name": "Enz", "pdb_path": pdb_path}
    data["substrate"] = {"name": "Lig", "sdf_path": sdf_path}
    config = SimulationConfig(**data)

    warnings = collect_reference_warnings(config)

    assert warnings == []


def test_collect_reference_warnings_reports_cached_polymer_references(tmp_path: Path) -> None:
    """Cached polymer configs should warn for missing directories and SDFs."""

    data = _minimal_config_data(tmp_path)
    data["polymers"] = {
        "enabled": True,
        "generation_mode": "cached",
        "type_prefix": "PEG",
        "monomers": [{"label": "A", "probability": 1.0}],
        "length": 4,
        "count": 1,
        "sdf_directory": tmp_path / "missing_polymers",
    }
    config = SimulationConfig(**data)

    warnings = collect_reference_warnings(config)

    assert any("Missing polymer SDF directory" in warning for warning in warnings)


def test_collect_reference_warnings_reports_empty_cached_polymer_directory(
    tmp_path: Path,
) -> None:
    """Cached polymer directories should contain matching charged SDF files."""

    polymer_dir = tmp_path / "polymers"
    polymer_dir.mkdir()
    data = _minimal_config_data(tmp_path)
    data["polymers"] = {
        "enabled": True,
        "generation_mode": "cached",
        "type_prefix": "PEG",
        "monomers": [{"label": "A", "probability": 1.0}],
        "length": 4,
        "count": 1,
        "sdf_directory": polymer_dir,
    }
    config = SimulationConfig(**data)

    warnings = collect_reference_warnings(config)

    assert any("Missing cached polymer SDF files" in warning for warning in warnings)


def test_custom_reaction_templates_reject_at_schema_boundary(tmp_path: Path) -> None:
    """Dynamic polymer custom reaction paths should be rejected immediately."""

    data = _minimal_config_data(tmp_path)
    data["polymers"] = {
        "enabled": True,
        "generation_mode": "dynamic",
        "type_prefix": "PEG",
        "monomers": [{"label": "A", "probability": 1.0, "smiles": "C=C"}],
        "length": 4,
        "count": 1,
        "reactions": {
            "initiation": tmp_path / "missing_init.rxn",
            "polymerization": tmp_path / "missing_poly.rxn",
            "termination": tmp_path / "missing_term.rxn",
        },
    }
    with pytest.raises(ValueError, match="provided_molecules"):
        SimulationConfig(**data)
