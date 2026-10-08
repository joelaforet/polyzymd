"""Tests for runtime configuration reference validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.config.schema import SimulationConfig
from polyzymd.config.validation import collect_reference_warnings, require_inputs


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


def test_collect_reference_warnings_reports_missing_reaction_templates(tmp_path: Path) -> None:
    """Dynamic polymer custom reaction paths should be checked at runtime."""

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
    config = SimulationConfig(**data)

    warnings = collect_reference_warnings(config)

    assert any("Missing polymer initiation reaction template" in warning for warning in warnings)
    assert any(
        "Missing polymer polymerization reaction template" in warning for warning in warnings
    )
    assert any("Missing polymer termination reaction template" in warning for warning in warnings)


_ONE_ATOM_PDB = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
)


def _buildable(tmp_path: Path, **sections) -> SimulationConfig:
    """A config with a one-atom enzyme PDB and the given extra sections."""
    (tmp_path / "enz.pdb").write_text(_ONE_ATOM_PDB)
    data = _minimal_config_data(tmp_path)
    data["enzyme"] = {"name": "Enz", "pdb_path": tmp_path / "enz.pdb"}
    data.update(sections)
    return SimulationConfig(**data)


def test_inputs_of_a_buildable_config_pass(tmp_path: Path) -> None:
    require_inputs(_buildable(tmp_path))


def test_an_enzyme_pdb_without_atoms_is_refused(tmp_path: Path) -> None:
    config = _buildable(tmp_path)
    (tmp_path / "enz.pdb").write_text("garbage\n")
    with pytest.raises(ValueError, match="no ATOM or HETATM records"):
        require_inputs(config)


def test_a_substrate_conformer_the_sdf_lacks_is_refused(tmp_path: Path) -> None:
    from rdkit import Chem

    (tmp_path / "lig.sdf").write_text(Chem.MolToMolBlock(Chem.MolFromSmiles("CCO")) + "$$$$\n")
    substrate = {"name": "lig", "sdf_path": tmp_path / "lig.sdf", "conformer_index": 1}
    with pytest.raises(ValueError, match="holds 1 conformer"):
        require_inputs(_buildable(tmp_path, substrate=substrate))
    (tmp_path / "lig.sdf").write_text("garbage\n")
    with pytest.raises(ValueError, match="cannot read conformer 0"):
        require_inputs(_buildable(tmp_path, substrate={**substrate, "conformer_index": 0}))


@pytest.mark.parametrize(
    ("smiles", "message"),
    [
        ("C1CC", "cannot read the SMILES"),
        ("[K+]", "the build adds only Na\\+ and Cl- ions"),
        ("C[Si](C)(C)O", "Si, which the NAGL charge model does not cover"),
    ],
)
def test_cosolvents_the_build_cannot_make_are_refused(tmp_path: Path, smiles, message) -> None:
    """Unreadable SMILES, ions other than NaCl and elements NAGL cannot charge are refused."""
    solvent = {"co_solvents": [{"name": "x", "smiles": smiles, "concentration": 0.1}]}
    with pytest.raises(ValueError, match=message):
        require_inputs(_buildable(tmp_path, solvent=solvent))


def test_a_monomer_without_a_polymerisable_group_is_refused(tmp_path: Path) -> None:
    """The default ATRP initiation reacts with methacrylates only."""
    polymers = {
        "generation_mode": "dynamic",
        "type_prefix": "P",
        "length": 5,
        "count": 1,
        "reactions": {
            "initiation": "default",
            "polymerization": "default",
            "termination": "default",
        },
        "monomers": [{"label": "A", "probability": 1.0, "smiles": "CC(=C)C(=O)OCCO"}],
    }
    require_inputs(_buildable(tmp_path, polymers=polymers))
    polymers["monomers"] = [{"label": "A", "probability": 1.0, "smiles": "CCO"}]
    with pytest.raises(ValueError, match="no group the initiation reaction"):
        require_inputs(_buildable(tmp_path, polymers=polymers))


def test_a_force_field_that_is_not_installed_is_refused(tmp_path: Path) -> None:
    config = _buildable(tmp_path, force_field={"protein": "ff14sb_typo.offxml"})
    with pytest.raises(ValueError, match="not an installed force field"):
        require_inputs(config)
