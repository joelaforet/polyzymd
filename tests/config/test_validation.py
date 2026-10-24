"""Tests for runtime configuration reference validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.config.schema import SimulationConfig
from polyzymd.config.validation import collect_reference_warnings, require_inputs

SAGE = ["openff-2.3.0.offxml"]

#: NAGL charges with the model the config must name.
NAGL = {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"}

#: The solvent section's required keys: the small-molecule force field and the water model.
SOLVENT = {"force_field": SAGE, "water": {"model": "tip3p"}}


def _protein(pdb_path: Path, **keys: object) -> dict[str, object]:
    """A protein solutes entry."""
    return {
        "type": "protein",
        "name": "Enz",
        "pdb_path": pdb_path,
        "force_field": ["ff14sb_off_impropers_0.0.4.offxml"],
        "charges": {"method": "library"},
        **keys,
    }


def _solute(type_: str, **keys: object) -> dict[str, object]:
    """A non-protein solutes entry with the solvent's force field."""
    return {"type": type_, "force_field": SAGE, **keys}


def _minimal_config_data(tmp_path: Path) -> dict[str, object]:
    """Return minimal simulation config data for reference-validation tests."""

    return {
        "name": "reference_validation",
        "engine": "openmm",
        "solutes": [_protein(tmp_path / "missing.pdb")],
        "solvent": dict(SOLVENT),
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
    data["solutes"] = [
        _protein(pdb_path),
        _solute("small_molecule", name="Lig", sdf_path=sdf_path, charges=NAGL),
    ]
    config = SimulationConfig(**data)

    warnings = collect_reference_warnings(config)

    assert warnings == []


def test_collect_reference_warnings_reports_cached_polymer_references(tmp_path: Path) -> None:
    """Cached polymer configs should warn for missing directories and SDFs."""

    data = _minimal_config_data(tmp_path)
    data["solutes"].append(
        {
            "type": "polymer",
            "force_field": SAGE,
            "charges": {"method": "preset"},
            "generation_mode": "cached",
            "type_prefix": "PEG",
            "monomers": [{"label": "A", "probability": 1.0}],
            "length": 4,
            "count": 1,
            "sdf_directory": tmp_path / "missing_polymers",
        }
    )
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
    data["solutes"].append(
        {
            "type": "polymer",
            "force_field": SAGE,
            "charges": {"method": "preset"},
            "generation_mode": "cached",
            "type_prefix": "PEG",
            "monomers": [{"label": "A", "probability": 1.0}],
            "length": 4,
            "count": 1,
            "sdf_directory": polymer_dir,
        }
    )
    config = SimulationConfig(**data)

    warnings = collect_reference_warnings(config)

    assert any("Missing cached polymer SDF files" in warning for warning in warnings)


def test_collect_reference_warnings_has_nothing_to_check_for_dynamic_polymers(
    tmp_path: Path,
) -> None:
    """Dynamic polymers are built from SMILES, so they reference no files."""

    data = _minimal_config_data(tmp_path)
    data["solutes"].append(
        {
            "type": "polymer",
            "force_field": SAGE,
            "charges": NAGL,
            "generation_mode": "dynamic",
            "type_prefix": "PEG",
            "monomers": [{"label": "A", "probability": 1.0, "smiles": "[*]CC([*])(C)C(=O)OC"}],
            "end_groups": {
                "head": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
                "tail": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
            },
            "length": 4,
            "count": 1,
        }
    )
    config = SimulationConfig(**data)

    warnings = collect_reference_warnings(config)

    assert not any("polymer" in warning.lower() for warning in warnings)


def test_am1bcc_without_ambertools_is_reported_before_the_build(
    tmp_path: Path, monkeypatch
) -> None:
    """A config that asks for AM1-BCC warns when no toolkit here can compute it."""
    from polyzymd.utils import charging

    monkeypatch.setattr(charging, "am1bcc_available", lambda: False)
    data = _minimal_config_data(tmp_path)
    data["solvent"] = {
        **SOLVENT,
        "co_solvents": [
            {"name": "etac", "smiles": "CCOC(C)=O", "count": 5, "charges": {"method": "am1bcc"}}
        ],
    }
    warnings = collect_reference_warnings(SimulationConfig(**data))
    assert any("AmberTools" in warning for warning in warnings)

    monkeypatch.setattr(charging, "am1bcc_available", lambda: True)
    assert not any("AmberTools" in w for w in collect_reference_warnings(SimulationConfig(**data)))


_ONE_ATOM_PDB = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
)


def _buildable(tmp_path: Path, **sections) -> SimulationConfig:
    """A config with a one-atom enzyme PDB and the given extra sections.

    ``protein``, ``small_molecule`` and ``polymer`` hold the keys of those solutes.
    """
    (tmp_path / "enz.pdb").write_text(_ONE_ATOM_PDB)
    data = _minimal_config_data(tmp_path)
    data["solutes"] = [_protein(tmp_path / "enz.pdb", **sections.pop("protein", {}))]
    for type_ in ("small_molecule", "polymer"):
        if type_ in sections:
            data["solutes"].append(_solute(type_, **sections.pop(type_)))
    if "solvent" in sections:
        sections["solvent"] = {**SOLVENT, **sections["solvent"]}
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
    substrate = {
        "name": "lig",
        "sdf_path": tmp_path / "lig.sdf",
        "conformer_index": 1,
        "charges": NAGL,
    }
    with pytest.raises(ValueError, match="holds 1 conformer"):
        require_inputs(_buildable(tmp_path, small_molecule=substrate))
    (tmp_path / "lig.sdf").write_text("garbage\n")
    with pytest.raises(ValueError, match="cannot read conformer 0"):
        require_inputs(_buildable(tmp_path, small_molecule={**substrate, "conformer_index": 0}))


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
    solvent = {
        "co_solvents": [{"name": "x", "smiles": smiles, "concentration": 0.1, "charges": NAGL}]
    }
    with pytest.raises(ValueError, match=message):
        require_inputs(_buildable(tmp_path, solvent=solvent))


def _dynamic_polymer(smiles: str, head: str = "[*][H]", tail: str = "[*][H]") -> dict:
    """The keys of a dynamic one-monomer polymer solute."""
    return {
        "generation_mode": "dynamic",
        "type_prefix": "P",
        "length": 5,
        "count": 1,
        "monomers": [{"label": "A", "probability": 1.0, "smiles": smiles}],
        "end_groups": {
            "head": [{"label": "H", "smiles": head, "probability": 1.0}],
            "tail": [{"label": "H", "smiles": tail, "probability": 1.0}],
        },
        "charges": NAGL,
    }


def test_star_sited_monomers_and_caps_are_accepted(tmp_path: Path) -> None:
    """A repeat unit with two '*' and caps with one pass the input checks."""
    polymer = _dynamic_polymer("[*]CC([*])(C)C(=O)OCCO", head="[*]C")
    require_inputs(_buildable(tmp_path, polymer=polymer))


@pytest.mark.parametrize(
    ("polymer", "message"),
    [
        (_dynamic_polymer("[*]CC([*])(C)C(=O"), "cannot read the SMILES"),
        (_dynamic_polymer("[*]=CC([*])C"), "bonded to one atom by a single bond"),
        (_dynamic_polymer("[*]C([*])C(=O)OC"), "both '\\*' .* are on one atom"),
        (_dynamic_polymer("[*]CC([*])(C)C(=O)OC", tail="[*]=C"), "end_groups.tail 'H'"),
        (_dynamic_polymer("[*]CC([*])[Si](C)(C)C"), "Si, which the NAGL charge model"),
    ],
    ids=["unreadable", "double_bond", "one_atom", "cap", "element"],
)
def test_star_sites_the_build_cannot_use_are_refused(tmp_path: Path, polymer, message) -> None:
    """Bad repeat units and caps are refused before the build, with a fix."""
    with pytest.raises(ValueError, match=message) as error:
        require_inputs(_buildable(tmp_path, polymer=polymer))
    assert "fix:" in str(error.value)


@pytest.mark.parametrize("residue_name", ["ALA", "SOL", "NA"])
def test_a_monomer_residue_name_the_selections_take_is_refused(
    tmp_path: Path, residue_name: str
) -> None:
    """A monomer named like a protein, water or ion residue would join those selections."""
    polymer = _dynamic_polymer("[*]CC([*])(C)C(=O)OC")
    polymer["monomers"][0]["residue_name"] = "MMA"
    require_inputs(_buildable(tmp_path, polymer=polymer))
    polymer["monomers"][0]["residue_name"] = residue_name
    config = _buildable(tmp_path, polymer=polymer)
    config.require_buildable(inputs=False)
    with pytest.raises(ValueError, match=f"residue_name '{residue_name}'") as error:
        require_inputs(config)
    assert "fix:" in str(error.value)


def test_every_cap_option_is_parsed(tmp_path: Path) -> None:
    """The second option of an end is checked like the first."""
    polymer = _dynamic_polymer("[*]CC([*])(C)C(=O)OC")
    polymer["end_groups"]["head"] = [
        {"label": "H", "smiles": "[*][H]", "probability": 0.5},
        {"label": "S", "smiles": "[*][Si](C)(C)C", "probability": 0.5},
    ]
    with pytest.raises(ValueError, match="end_groups.head 'S'") as error:
        require_inputs(_buildable(tmp_path, polymer=polymer))
    assert "Si" in str(error.value)


@pytest.mark.parametrize("residue_name", ["ALA", "HOH"])
def test_a_cap_residue_name_the_selections_take_is_refused(
    tmp_path: Path, residue_name: str
) -> None:
    """A cap that is its own residue may not take a protein or water residue name."""
    polymer = _dynamic_polymer("[*]CC([*])(C)C(=O)OC")
    polymer["end_groups"]["tail"] = [
        {"label": "V", "smiles": "[*]C=C(C)C(=O)OC", "probability": 1.0, "residue_name": "MMV"}
    ]
    require_inputs(_buildable(tmp_path, polymer=polymer))
    polymer["end_groups"]["tail"][0]["residue_name"] = residue_name
    with pytest.raises(ValueError, match="end_groups.tail 'V' has residue_name") as error:
        require_inputs(_buildable(tmp_path, polymer=polymer))
    assert "fix:" in str(error.value)


def test_a_force_field_that_is_not_installed_is_refused(tmp_path: Path) -> None:
    config = _buildable(tmp_path, protein={"force_field": ["ff14sb_typo.offxml"]})
    with pytest.raises(ValueError, match="not an installed force field"):
        require_inputs(config)


@pytest.mark.parametrize("smiles", ["CC(=O)[O-]", "CC(=O)[O-].[Na+]"])
def test_a_charged_cosolvent_with_ionic_strength_is_accepted(tmp_path: Path, smiles) -> None:
    """Charged co-solvents count toward ionic_strength, so validate accepts them."""
    solvent = {
        "co_solvents": [{"name": "x", "smiles": smiles, "concentration": 0.1, "charges": NAGL}],
        "ions": {"ionic_strength": 0.13},
    }
    require_inputs(_buildable(tmp_path, solvent=solvent))
