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
PRESET = {"method": "preset"}
NAGL = {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"}


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
        "charges": NAGL,
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
    data["solutes"][0]["charges"] = NAGL
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
    data["solutes"].append(_polymer(charges=NAGL))
    with pytest.raises(ValidationError, match="method: preset"):
        load(data)
    data["solutes"][-1] = _polymer(
        generation_mode="dynamic",
        monomers=[{"label": "A", "probability": 1.0, "smiles": "[*]CC([*])(C)C(=O)OC"}],
        end_groups=_caps("[*][H]", "[*][H]"),
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
    """Every co-solvent names its charges; preset needs a molecule PolyzyMD ships charges for."""
    data["solvent"]["co_solvents"] = [{"name": "dmso", "concentration": 0.2}]
    with pytest.raises(
        ValidationError, match=r"needs charges, for example charges: \{method: preset\}"
    ):
        load(data)
    data["solvent"]["co_solvents"][0]["charges"] = PRESET
    assert load(data).solvent.co_solvents[0].charges.method == "preset"
    data["solvent"]["co_solvents"] = [{"name": "x", "smiles": "CCO", "count": 1}]
    with pytest.raises(ValidationError, match="needs charges"):
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
    with pytest.raises(
        ValidationError, match=r"shipped charges .* write charges: \{method: preset\}"
    ):
        load(data)


def test_a_new_charge_method_on_a_bundled_co_solvent_is_a_new_analysis_identity(data) -> None:
    """Shipped charges keep the identity runs built before rc5 had; NAGL on dmso does not."""
    from polyzymd.analyses.identity import compute_config_hash

    data["solvent"]["co_solvents"] = [{"name": "dmso", "count": 10, "charges": PRESET}]
    shipped = compute_config_hash(load(data))
    data["solvent"]["co_solvents"][0]["charges"] = {"method": "preset"}
    assert compute_config_hash(load(data)) == shipped
    data["solvent"]["co_solvents"][0]["charges"] = NAGL
    assert compute_config_hash(load(data)) != shipped


def test_preset_charges_need_the_shipped_molecule() -> None:
    """dmso with another SMILES has no shipped charges; validate says so before a build."""
    pytest.importorskip("rdkit")
    from polyzymd.config.schema import CoSolventSpec
    from polyzymd.config.validation import _require_shipped_molecule

    _require_shipped_molecule(CoSolventSpec(name="dmso", count=1, charges=PRESET))
    with pytest.raises(ValueError, match="not the dmso PolyzyMD ships"):
        _require_shipped_molecule(CoSolventSpec(name="dmso", smiles="CCO", count=1, charges=PRESET))


#: Star-sited repeat units of methyl methacrylate and HEMA.
MMA = "[*]CC([*])(C)C(=O)OC"
HEMA = "[*]CC([*])(C)C(=O)OCCO"


def _caps(head: str, tail: str) -> dict:
    """End groups with one fixed cap on each end (head label H, tail label T)."""
    return {
        "head": [{"label": "H", "smiles": head, "probability": 1.0}],
        "tail": [{"label": "T", "smiles": tail, "probability": 1.0}],
    }


def _dynamic(**extra) -> dict:
    """A dynamic MMA/HEMA polymer solute with hydrogen caps."""
    return _polymer(
        **{
            "generation_mode": "dynamic",
            "charges": NAGL,
            "sdf_directory": None,
            "monomers": [
                {"label": "A", "probability": 0.5, "name": "MMA", "smiles": MMA},
                {"label": "B", "probability": 0.5, "name": "HEMA", "smiles": HEMA},
            ],
            "end_groups": _caps("[*][H]", "[*]C"),
            **extra,
        }
    )


def test_a_dynamic_polymer_loads_with_star_sited_smiles(data) -> None:
    data["solutes"].append(_dynamic())
    polymers = load(data).polymers
    [head], [tail] = polymers.end_groups.head, polymers.end_groups.tail
    assert (head.label, head.smiles, head.probability) == ("H", "[*][H]", 1.0)
    assert (tail.label, tail.smiles, tail.residue_name) == ("T", "[*]C", None)
    assert [m.residue_name for m in polymers.monomers] == ["MMA", "HEM"]


def test_reactions_are_refused_with_the_rewrite_guide(data) -> None:
    """Configs written for the ATRP reaction templates are told how to rewrite them."""
    reactions = {"initiation": "default", "polymerization": "default", "termination": "default"}
    data["solutes"].append(_dynamic(reactions=reactions))
    with pytest.raises(ValidationError, match="rewrite_config_rc5") as error:
        load(data)
    assert "1.3.0rc5" in str(error.value)
    assert "end_groups" in str(error.value)


def test_a_dynamic_polymer_needs_its_end_groups(data) -> None:
    data["solutes"].append(_dynamic(end_groups=None))
    with pytest.raises(ValidationError, match="requires 'end_groups'"):
        load(data)


@pytest.mark.parametrize(
    ("monomer", "head", "message"),
    [
        ("CC(C)C(=O)OC", "[*][H]", "has 0 '\\*' sites"),
        ("[*]CC(C)C(=O)OC", "[*][H]", "has 1 '\\*' sites"),
        ("[*]CC([*])(C[*])C(=O)OC", "[*][H]", "has 3 '\\*' sites"),
        (MMA, "[H]", "end_groups.head 'H': '\\[H\\]' has 0"),
        (MMA, "[*]C[*]", "end_groups.head 'H': .* has 2"),
    ],
)
def test_star_counts_are_checked_at_load(data, monomer, head, message) -> None:
    """A repeat unit needs exactly two '*' and a cap exactly one."""
    data["solutes"].append(
        _dynamic(
            monomers=[{"label": "A", "probability": 1.0, "name": "X", "smiles": monomer}],
            end_groups=_caps(head, "[*][H]"),
        )
    )
    with pytest.raises(ValidationError, match=message):
        load(data)


@pytest.mark.parametrize(
    "names",
    [
        {"residue_name": "PMA"},
        {},
    ],
    ids=["explicit", "from_name"],
)
def test_monomers_cannot_share_a_residue_name(data, names) -> None:
    """Analyses tell monomers apart by residue name, given or taken from the name."""
    monomers = [
        {"label": "A", "probability": 0.5, "name": "HEMA", "smiles": HEMA, **names},
        {"label": "B", "probability": 0.5, "name": "HEMA2", "smiles": MMA, **names},
    ]
    data["solutes"].append(_dynamic(monomers=monomers))
    with pytest.raises(ValidationError, match="share residue_name"):
        load(data)


@pytest.mark.parametrize(
    ("second", "message"),
    [
        ({"label": "A", "name": "HEMA", "residue_name": "HEM"}, "share label \\['A'\\]"),
        ({"label": "B", "name": "MMA", "residue_name": "HEM"}, "share name \\['MMA'\\]"),
    ],
    ids=["label", "name"],
)
def test_monomers_cannot_share_a_label_or_name(data, second, message) -> None:
    """Each monomer is looked up by label and by name, so neither may repeat."""
    monomers = [
        {"label": "A", "probability": 0.5, "name": "MMA", "smiles": MMA, "residue_name": "MMA"},
        {"probability": 0.5, "smiles": HEMA, **second},
    ]
    data["solutes"].append(_dynamic(monomers=monomers))
    with pytest.raises(ValidationError, match=message):
        load(data)


def test_a_dynamic_polymer_config_checks_without_the_build_packages(data, tmp_path) -> None:
    """require_buildable(inputs=False) runs where MDAnalysis, RDKit and mBuild are not installed.

    The sim-cuda environments check a config before running it and have none of
    them; the star-site and residue-name checks that need them run in validate.
    """
    import json
    import subprocess
    import sys

    data["solutes"].append(
        _dynamic(
            monomers=[
                {
                    "label": "A",
                    "probability": 0.5,
                    "name": "MMA",
                    "smiles": MMA,
                    "residue_name": "MMA",
                },
                {"label": "B", "probability": 0.5, "name": "HEMA", "smiles": HEMA},
            ]
        )
    )
    data.pop("_ligand")
    config = tmp_path / "config.json"
    config.write_text(json.dumps(data))
    script = (
        "import json, sys\n"
        "for name in ('MDAnalysis', 'rdkit', 'mbuild'):\n"
        "    sys.modules[name] = None\n"
        "from polyzymd.config.schema import SimulationConfig\n"
        f"data = json.loads(open({str(config)!r}).read())\n"
        "SimulationConfig.model_validate(data).require_buildable(inputs=False)\n"
    )
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("length", [1, 2])
def test_dynamic_chains_without_a_repeat_unit_are_refused(data, length) -> None:
    """A dynamic chain's length counts its two caps, so it is at least 3."""
    data["solutes"].append(_dynamic(length=length))
    with pytest.raises(ValidationError, match="counts its two ends") as error:
        load(data)
    assert "at least 3" in str(error.value)


def test_a_dynamic_chain_of_one_repeat_unit_loads(data) -> None:
    data["solutes"].append(_dynamic(length=3))
    assert load(data).polymers.length == 3


def _option(label: str, smiles: str, probability: float, **extra) -> dict:
    return {"label": label, "smiles": smiles, "probability": probability, **extra}


def test_each_end_draws_from_options_with_probabilities(data) -> None:
    """Two head options and one own-residue tail cap load as given."""
    end_groups = {
        "head": [
            _option("H", "[*][H]", 0.25),
            _option("V", "[*]C=C(C)C(=O)OC", 0.75, name="MMA vinyl", residue_name="MMV"),
        ],
        "tail": [_option("H", "[*][H]", 1.0)],
    }
    data["solutes"].append(_dynamic(end_groups=end_groups, length=5))
    polymers = load(data).polymers
    assert [(o.label, o.probability) for o in polymers.end_groups.head] == [
        ("H", 0.25),
        ("V", 0.75),
    ]
    assert polymers.end_groups.head[1].residue_name == "MMV"
    assert polymers.end_groups.head[1].name == "MMA vinyl"


@pytest.mark.parametrize(
    ("head", "message"),
    [
        (
            [_option("H", "[*][H]", 0.5), _option("M", "[*]C", 0.4)],
            "end_groups.head probabilities must sum to 1.0, got 0.9",
        ),
        ([_option("H", "[*][H]", 0.5)], "a single cap needs probability 1.0"),
        (
            [_option("H", "[*][H]", 0.5), _option("H", "[*]C", 0.5)],
            "end_groups.head options share label \\['H'\\]",
        ),
        ([], "at least 1 item"),
        ([_option("HH", "[*][H]", 1.0)], "at most 1 character"),
    ],
    ids=["sum_0.9", "single_0.5", "duplicate_label", "empty", "long_label"],
)
def test_cap_options_are_checked(data, head, message) -> None:
    end_groups = {"head": head, "tail": [_option("H", "[*][H]", 1.0)]}
    data["solutes"].append(_dynamic(end_groups=end_groups))
    with pytest.raises(ValidationError, match=message):
        load(data)


def test_a_bare_smiles_cap_is_refused_with_the_list_form(data) -> None:
    data["solutes"].append(_dynamic(end_groups={"head": "[*][H]", "tail": "[*]C"}))
    with pytest.raises(ValidationError, match="end_groups.head is a list of caps") as error:
        load(data)
    assert "head: [{label: H, smiles: '[*][H]', probability: 1.0}]" in str(error.value)


#: Methyl methacrylate with the handedness of its quaternary carbon set.
MMA_R = "[*]C[C@@]([*])(C)C(=O)OC"
MMA_S = "[*]C[C@]([*])(C)C(=O)OC"


def test_stereo_letters_of_one_monomer_may_share_a_residue_name(data) -> None:
    """Two stereoisomers of MMA (the letters that set tacticity) are both residue MMA."""
    monomers = [
        {"label": "A", "probability": 0.5, "name": "MMA_R", "smiles": MMA_R, "residue_name": "MMA"},
        {"label": "B", "probability": 0.5, "name": "MMA_S", "smiles": MMA_S, "residue_name": "MMA"},
    ]
    data["solutes"].append(_dynamic(monomers=monomers))
    assert [m.residue_name for m in load(data).polymers.monomers] == ["MMA", "MMA"]


def test_monomers_that_are_not_stereoisomers_cannot_share_a_residue_name(data) -> None:
    monomers = [
        {"label": "A", "probability": 0.5, "name": "MMA_R", "smiles": MMA_R, "residue_name": "MMA"},
        {
            "label": "B",
            "probability": 0.5,
            "name": "HEMA_S",
            "smiles": "[*]C[C@]([*])(C)C(=O)OCCO",
            "residue_name": "MMA",
        },
    ]
    data["solutes"].append(_dynamic(monomers=monomers))
    with pytest.raises(ValidationError, match="only stereoisomers of one monomer"):
        load(data)


def test_two_monomers_with_the_same_smiles_cannot_share_a_residue_name(data) -> None:
    """Only stereoisomers share a residue name; the same SMILES twice is refused."""
    monomers = [
        {"label": "A", "probability": 0.5, "name": "MMA_1", "smiles": MMA_R, "residue_name": "MMA"},
        {"label": "B", "probability": 0.5, "name": "MMA_2", "smiles": MMA_R, "residue_name": "MMA"},
    ]
    data["solutes"].append(_dynamic(monomers=monomers))
    with pytest.raises(ValidationError, match="share residue_name"):
        load(data)
