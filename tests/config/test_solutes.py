"""The rc5 config layout: solutes and solvent, each with its own force field and charges."""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from polyzymd.config.schema import REWRITE_GUIDE, SimulationConfig

QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart"
SAGE = ["openff-2.3.0.offxml"]


@pytest.fixture
def data(tmp_path) -> dict:
    """The quickstart config as data, with a substrate SDF beside it."""
    config = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    config["solutes"][0]["pdb_path"] = str(QUICKSTART / "trpcage.pdb")
    (tmp_path / "ligand.sdf").write_text("")
    config["_ligand"] = str(tmp_path / "ligand.sdf")
    return config


def _ligand(data: dict, **extra) -> dict:
    entry = {
        "type": "small_molecule",
        "sdf_path": data["_ligand"],
        "force_field": SAGE,
        "charges": {"method": "nagl"},
    }
    return {**entry, **extra}


def _polymer(**extra) -> dict:
    entry = {
        "type": "polymer",
        "type_prefix": "TP",
        "monomers": [{"label": "A", "probability": 1.0}],
        "length": 3,
        "count": 2,
        "sdf_directory": "polymers",
        "force_field": SAGE,
        "charges": {"method": "preset"},
    }
    return {**entry, **extra}


def load(data: dict) -> SimulationConfig:
    data = copy.deepcopy(data)
    data.pop("_ligand", None)
    return SimulationConfig.model_validate(data)


def test_the_quickstart_reads_through_the_builder_views(data) -> None:
    """enzyme, substrate, polymers and force_field are views the builders still read."""
    data["solutes"].append(_ligand(data, name="substrate"))
    config = load(data)
    assert config.enzyme.name == "trpcage"
    assert config.substrate.name == "substrate"
    assert config.substrate.charge_method.value == "nagl"
    assert config.polymers is None
    assert config.force_field.protein == "ff14sb_off_impropers_0.0.4.offxml"
    assert config.force_field.small_molecule == "openff-2.3.0.offxml"
    assert config.force_field.nagl_model == "openff-gnn-am1bcc-1.0.0.pt"


@pytest.mark.parametrize("section", ["enzyme", "substrate", "polymers", "force_field"])
def test_an_older_config_is_told_where_its_sections_went(data, section) -> None:
    data[section] = {}
    with pytest.raises(ValidationError, match=f"{section} -> .*{REWRITE_GUIDE}"):
        load(data)


def test_solvent_primary_is_now_water(data) -> None:
    data["solvent"]["primary"] = {"type": "water", "model": "tip3p"}
    with pytest.raises(ValidationError, match="solvent.water"):
        load(data)


def test_name_and_segid_default_from_the_type(data) -> None:
    del data["solutes"][0]["name"]
    data["solutes"].append(_polymer())
    config = load(data)
    assert [(s.name, s.segid) for s in config.solutes] == [
        ("protein", "PROT"),
        ("polymer", "POLY"),
    ]


@pytest.mark.parametrize(
    "segid, message",
    [("TOOLONG", "segment ID"), ("SOL", "water and co-solvents"), ("ION", "ions")],
)
def test_a_segid_must_fit_a_psf_and_not_take_the_solvents(data, segid, message) -> None:
    data["solutes"][0]["segid"] = segid
    with pytest.raises(ValidationError, match=message):
        load(data)


def test_two_solutes_cannot_share_a_name_or_segid(data) -> None:
    data["solutes"].append(_ligand(data, name="trpcage"))
    with pytest.raises(ValidationError, match="used twice"):
        load(data)
    data["solutes"][-1] = _ligand(data, segid="PROT")
    with pytest.raises(ValidationError, match="segid PROT"):
        load(data)


def test_an_unknown_type_is_refused_not_guessed(data) -> None:
    data["solutes"].append({**_ligand(data), "type": "ligand"})
    with pytest.raises(ValidationError, match="type must be protein, small_molecule or polymer"):
        load(data)


def test_this_release_builds_one_solute_per_type_around_a_protein(data) -> None:
    data["solutes"] += [_ligand(data, name="a"), _ligand(data, name="b")]
    with pytest.raises(ValidationError, match="2 of type small_molecule"):
        load(data)
    data["solutes"] = [_ligand(data)]
    with pytest.raises(ValidationError, match="type: protein"):
        load(data)


def test_a_protein_takes_library_charges(data) -> None:
    data["solutes"][0]["charges"] = {"method": "nagl"}
    with pytest.raises(ValidationError, match="library charges"):
        load(data)


def test_nagl_names_its_model_and_other_methods_take_none(data) -> None:
    data["solutes"].append(_ligand(data))
    assert load(data).substrate.charges.model == "openff-gnn-am1bcc-1.0.0.pt"
    data["solutes"][-1] = _ligand(data, charges={"method": "am1bcc", "model": "x.pt"})
    with pytest.raises(ValidationError, match="does not use one"):
        load(data)


def test_polymer_charges_follow_how_the_chains_are_made(data) -> None:
    """Cached chains carry their SDF charges; dynamic chains are charged as they are built."""
    data["solutes"].append(_polymer(charges={"method": "nagl"}))
    with pytest.raises(ValidationError, match="method: preset"):
        load(data)
    data["solutes"][-1] = _polymer(
        generation_mode="dynamic",
        monomers=[{"label": "A", "probability": 1.0, "smiles": "C=C(C)C(=O)OC"}],
        reactions={"initiation": "default", "polymerization": "default", "termination": "default"},
    )
    with pytest.raises(ValidationError, match="nagl, am1bcc or espaloma"):
        load(data)


def test_every_molecule_but_the_protein_shares_the_solvent_force_field(data) -> None:
    data["solutes"].append(_ligand(data, force_field=["openff-2.2.0.offxml"]))
    with pytest.raises(ValidationError, match="must be the same"):
        load(data)


def test_one_nagl_model_charges_every_component(data) -> None:
    data["solutes"].append(_ligand(data, charges={"method": "nagl", "model": "a.pt"}))
    data["solvent"]["co_solvents"] = [
        {"name": "x", "smiles": "CCO", "count": 1, "charges": {"method": "nagl", "model": "b.pt"}}
    ]
    with pytest.raises(ValidationError, match="different NAGL models"):
        load(data)


def test_co_solvent_charges_say_what_is_built(data) -> None:
    """Shipped charges by default for a bundled molecule; any other must name its method."""
    data["solvent"]["co_solvents"] = [{"name": "dmso", "concentration": 0.2}]
    assert load(data).solvent.co_solvents[0].charges.method == "preset"
    data["solvent"]["co_solvents"] = [{"name": "x", "smiles": "CCO", "count": 1}]
    with pytest.raises(ValidationError, match="no charges shipped"):
        load(data)
    data["solvent"]["co_solvents"][0]["charges"] = {"method": "preset"}
    with pytest.raises(ValidationError, match="ships none"):
        load(data)
    data["solvent"]["co_solvents"] = [{"name": "dmso", "count": 1, "charge_method": "nagl"}]
    with pytest.raises(ValidationError, match="now 'charges'"):
        load(data)


def test_a_restraint_atom_names_a_solute(data) -> None:
    data["restraints"] = [
        {
            "type": "flat_bottom",
            "name": "r",
            "atom1": {"component": "trpcage", "selection": "resid 1 and name CA"},
            "atom2": {"component": "ligand", "selection": "name C1"},
            "distance": 5.0,
        }
    ]
    with pytest.raises(ValidationError, match="component 'ligand', which is not a solute"):
        load(data)
    del data["restraints"][0]["atom2"]["component"]
    with pytest.raises(ValidationError, match="component"):
        load(data)


@pytest.mark.parametrize("name", [" ", "x/../y"])
def test_a_solute_name_is_a_single_folder_name(data, name) -> None:
    """Run folders are named after solutes, so a name must be one folder name."""
    data["solutes"][0]["name"] = name
    with pytest.raises(ValidationError, match="non-empty, without '/'"):
        load(data)


@pytest.mark.parametrize("name", ["solvent", "ions", "SOL", "my enzyme"])
def test_a_solute_cannot_take_a_built_in_name_or_hold_spaces(data, name) -> None:
    """Built-in components and GROMACS index groups need names no solute shares."""
    data["solutes"][0]["name"] = name
    with pytest.raises(ValidationError, match="taken by a component|spaces"):
        load(data)


def test_the_old_co_solvent_key_says_bundled_molecules_kept_their_charges(data) -> None:
    """Rewriting charge_method: nagl as charges nagl on dmso would change its charges."""
    data["solvent"]["co_solvents"] = [{"name": "dmso", "count": 1, "charge_method": "nagl"}]
    with pytest.raises(ValidationError, match="shipped charges .* delete the key"):
        load(data)


def test_a_new_charge_method_on_a_bundled_co_solvent_is_a_new_analysis_identity(data) -> None:
    """Shipped charges keep the identity runs built before rc5 had; NAGL on dmso does not."""
    from polyzymd.analyses.identity import compute_config_hash

    data["solvent"]["co_solvents"] = [{"name": "dmso", "count": 10}]
    shipped = compute_config_hash(load(data))
    data["solvent"]["co_solvents"][0]["charges"] = {"method": "preset"}
    assert compute_config_hash(load(data)) == shipped
    data["solvent"]["co_solvents"][0]["charges"] = {"method": "nagl"}
    assert compute_config_hash(load(data)) != shipped


def test_preset_charges_need_the_shipped_molecule() -> None:
    """dmso with another SMILES has no shipped charges; validate says so before a build."""
    pytest.importorskip("rdkit")
    from polyzymd.config.schema import CoSolventSpec
    from polyzymd.config.validation import _require_shipped_molecule

    _require_shipped_molecule(CoSolventSpec(name="dmso", count=1))
    with pytest.raises(ValueError, match="not the dmso PolyzyMD ships"):
        _require_shipped_molecule(CoSolventSpec(name="dmso", smiles="CCO", count=1))
