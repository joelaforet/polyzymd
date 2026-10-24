"""Tests for runtime configuration reference validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.config.schema import SimulationConfig
from polyzymd.config.validation import (
    collect_protocol_warnings,
    collect_reference_warnings,
    require_inputs,
)

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


def test_collect_reference_warnings_reports_missing_reaction_templates(tmp_path: Path) -> None:
    """Dynamic polymer custom reaction paths should be checked at runtime."""

    data = _minimal_config_data(tmp_path)
    data["solutes"].append(
        {
            "type": "polymer",
            "force_field": SAGE,
            "charges": NAGL,
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
    )
    config = SimulationConfig(**data)

    warnings = collect_reference_warnings(config)

    assert any("Missing polymer initiation reaction template" in warning for warning in warnings)
    assert any(
        "Missing polymer polymerization reaction template" in warning for warning in warnings
    )
    assert any("Missing polymer termination reaction template" in warning for warning in warnings)


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
        "charges": NAGL,
    }
    require_inputs(_buildable(tmp_path, polymer=polymers))
    polymers["monomers"] = [{"label": "A", "probability": 1.0, "smiles": "CCO"}]
    with pytest.raises(ValueError, match="no group the initiation reaction"):
        require_inputs(_buildable(tmp_path, polymer=polymers))


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


def _benzene_sdf(tmp_path: Path) -> Path:
    """An SDF of benzene with hydrogens, as an SDF from a crystal ligand: no atom names."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(molecule, randomSeed=1)
    path = tmp_path / "benzene.sdf"
    path.write_text(Chem.MolToMolBlock(molecule) + "$$$$\n")
    return path


def _ligand_restraint(selection: str, **keys: object) -> dict[str, object]:
    return {
        "type": "flat_bottom",
        "name": "in_cavity",
        "atom1": {"component": "Enz", "selection": "name CA"},
        "atom2": {"component": "benzene", "selection": selection},
        "distance": 7.0,
        **keys,
    }


def test_small_molecule_atom_names_follow_openff() -> None:
    """The names validate predicts are the names OpenFF gives the built topology."""
    from openff.toolkit import Molecule, Topology
    from rdkit import Chem

    from polyzymd.config.validation import small_molecule_atom_names

    rdkit_molecule = Chem.AddHs(Chem.MolFromSmiles("OCc1ccccc1Cl"))
    molecule = Molecule.from_rdkit(rdkit_molecule)
    for atom in molecule.atoms:
        atom.metadata["residue_name"] = "LIG"
    built = [atom.name for atom in Topology.from_molecules([molecule]).to_openmm().atoms()]
    assert small_molecule_atom_names(rdkit_molecule) == built
    assert built[:3] == ["O1x", "C1x", "C2x"]


@pytest.mark.parametrize(
    ("selection", "names"),
    [
        ("name C1", ["C1"]),
        ("resname LIG and name C1x", ["C1x"]),
        ("(name C1 or name C2) and element C", ["C1", "C2"]),
        ("resname LIG and not name H1", None),
        ("name C1 C2", ["C1", "C2"]),
        ("name == 'C1'", ["C1"]),
        ("resid 1", []),
        ("name =~ 'C.*'", None),
    ],
)
def test_selection_atom_names(selection: str, names) -> None:
    from polyzymd.config.validation import selection_atom_names

    assert selection_atom_names(selection) == names


def test_a_restraint_on_a_ligand_atom_name_the_build_will_not_have_is_refused(
    tmp_path: Path,
) -> None:
    """name C1 matches nothing once OpenFF names benzene's atoms C1x ... H6x."""
    small_molecule = {"name": "benzene", "sdf_path": _benzene_sdf(tmp_path), "charges": NAGL}
    config = _buildable(
        tmp_path,
        small_molecule=small_molecule,
        restraints=[_ligand_restraint("name C1")],
    )
    with pytest.raises(ValueError) as error:
        require_inputs(config)
    message = str(error.value)
    assert "restraint 'in_cavity': atom2 selection 'name C1'" in message
    assert "C1x, C2x, C3x, C4x, C5x, C6x, H1x" in message
    assert "fix: in restraints.0.atom2.selection" in message


@pytest.mark.parametrize(
    "restraint",
    [
        _ligand_restraint("name C1x"),
        _ligand_restraint("resname LIG and name C3x"),
        _ligand_restraint("index 3"),
        _ligand_restraint("name C1", enabled=False),
    ],
)
def test_restraints_on_names_the_ligand_will_have_pass(tmp_path: Path, restraint) -> None:
    small_molecule = {"name": "benzene", "sdf_path": _benzene_sdf(tmp_path), "charges": NAGL}
    require_inputs(_buildable(tmp_path, small_molecule=small_molecule, restraints=[restraint]))


def test_an_sdf_that_names_every_atom_keeps_its_names(tmp_path: Path) -> None:
    """With unique names in an atom.prop._Name field, OpenFF keeps them, and so does the check."""
    from openff.toolkit import Molecule, Topology
    from rdkit import Chem

    from polyzymd.config.validation import small_molecule_atom_names

    path = _benzene_sdf(tmp_path)
    names = [f"C{i}" for i in range(1, 7)] + [f"H{i}" for i in range(1, 7)]
    text = path.read_text().replace("$$$$\n", f"> <atom.prop._Name>\n{' '.join(names)}\n\n$$$$\n")
    path.write_text(text)
    rdkit_molecule = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)))
    molecule = Molecule.from_file(str(path))
    built = [atom.name for atom in Topology.from_molecules([molecule]).to_openmm().atoms()]
    assert small_molecule_atom_names(rdkit_molecule) == built == names

    small_molecule = {"name": "benzene", "sdf_path": path, "charges": NAGL}
    restraint = _ligand_restraint("name C1")
    require_inputs(_buildable(tmp_path, small_molecule=small_molecule, restraints=[restraint]))


def test_collect_protocol_warnings_compares_stages_with_production(tmp_path: Path) -> None:
    """Constant stages and a final ramp's end are compared with thermodynamics.temperature."""
    data = _minimal_config_data(tmp_path)
    assert collect_protocol_warnings(SimulationConfig.model_validate(data)) == []

    data["thermodynamics"] = {"temperature": 330.0}
    data["simulation_phases"]["equilibration_stages"] = [
        {"name": "hold", "duration": 0.1, "temperature": 300.0, "ensemble": "NVT"},
        {
            "name": "heating",
            "ensemble": "NVT",
            "temperature_start": 100.0,
            "temperature_end": 300.0,
            "temperature_increment": 10.0,
            "temperature_interval_steps": 100,
        },
    ]
    warnings = collect_protocol_warnings(SimulationConfig.model_validate(data))

    assert len(warnings) == 2
    assert "stage 'hold' runs at 300 K" in warnings[0]
    assert "equilibration_stages.0.temperature" in warnings[0]
    assert "stage 'heating' ends at 300 K" in warnings[1]
    assert "equilibration_stages.1.temperature_end" in warnings[1]

    data["simulation_phases"]["equilibration_stages"][0]["temperature"] = 330.0
    data["simulation_phases"]["equilibration_stages"][1]["temperature_end"] = 330.0
    assert collect_protocol_warnings(SimulationConfig.model_validate(data)) == []
