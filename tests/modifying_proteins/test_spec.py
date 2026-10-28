"""Expanding and applying mutations: and modifications: blocks, and writing them back."""

import shutil
from pathlib import Path

import pytest
import yaml

from polyzymd import modifying_proteins as mp

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
CHITOBIOSE = ROOT / "tests" / "data" / "modifying_proteins" / "glycam_G42666HT.pdb"

GLYCAN_LINK = {
    "type": "glycan",
    "atom_name": "OG",
    "reaction": "O-glycosylation",
    "fragment_atom_name": "C1",
    "fragment_resnum": 2,
}


@pytest.fixture()
def library(tmp_path: Path) -> Path:
    folder = tmp_path / "glycans"
    folder.mkdir()
    shutil.copy(CHITOBIOSE, folder / "glycam_G42666HT.pdb")
    return folder


def test_group_mutation_expands_over_chains() -> None:
    mutations, _ = mp.expand(
        {"mutations": [{"chains": ["A", "C"], "residues": ["T301"], "to": "ASN"}]}
    )
    assert mutations == [{"to": "ASN", "site": "A:T301"}, {"to": "ASN", "site": "C:T301"}]


def test_glycan_group_expands_to_one_entry_per_site(library: Path) -> None:
    spec = {
        "modifications": [
            {
                "attach": {
                    **GLYCAN_LINK,
                    "chains": ["A", "B"],
                    "glycan_library": library.name,
                    "glycans": {"G42666HT": ["S13", "S20"]},
                }
            }
        ]
    }
    _, modifications = mp.expand(spec, base_dir=library.parent)
    sites = [entry["attach"]["site"] for entry in modifications]
    assert sites == ["A:S13", "A:S20", "B:S13", "B:S20"]
    first = modifications[0]["attach"]
    assert first["fragment"] == {"pdb_path": str(library / "glycam_G42666HT.pdb")}
    assert first["reaction"] == "O-glycosylation" and "glycans" not in first


def test_sites_file_names_each_site_and_its_glycan(library: Path) -> None:
    (library.parent / "sites.csv").write_text("site,glycan\nA:S13,G42666HT\nA:S20,G42666HT\n")
    spec = {
        "modifications": [
            {"attach": {**GLYCAN_LINK, "sites_file": "sites.csv", "glycan_library": "glycans"}}
        ]
    }
    _, modifications = mp.expand(spec, base_dir=library.parent)
    assert [entry["attach"]["site"] for entry in modifications] == ["A:S13", "A:S20"]
    assert all(entry["attach"]["fragment"]["pdb_path"].endswith(".pdb") for entry in modifications)


@pytest.mark.parametrize(
    ("spec", "message"),
    [
        ({"modifications": [{"mutate": {"site": "A:L2", "to": "ALA"}}]}, "under `mutations:`"),
        ({"modifications": [{"glycosylate": {"site": "A:N1"}}]}, "Unknown modification"),
        (
            {
                "modifications": [
                    {"attach": {"type": "generic", "sites": ["A:K8", "A:S13"], "name": "x"}}
                ]
            },
            "cannot set `name`",
        ),
        ({"mutations": [{"to": "ALA"}]}, "Give `site`"),
    ],
)
def test_expand_rejects_malformed_entries(spec, message) -> None:
    with pytest.raises(ValueError, match=message):
        mp.expand(spec)


def test_missing_glycan_file_names_what_is_expected(library: Path) -> None:
    spec = {
        "modifications": [
            {
                "attach": {
                    **GLYCAN_LINK,
                    "sites": ["A:S13"],
                    "glycans": {"G00000XX": ["A:S13"]},
                    "glycan_library": "glycans",
                }
            }
        ]
    }
    with pytest.raises(ValueError, match="glycam_G00000XX.pdb"):
        mp.expand(spec, base_dir=library.parent)


def test_spec_round_trip_rebuilds_the_same_protein(library: Path) -> None:
    spec = {
        "mutations": [{"site": "A:L2", "to": "ALA"}],
        "modifications": [
            {"deprotonate": {"site": "A:K8", "atom_name": "NZ"}},
            {
                "attach": {
                    "type": "generic",
                    "site": "A:K8",
                    "atom_name": "NZ",
                    "fragment": {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
                }
            },
            {
                "attach": {
                    **GLYCAN_LINK,
                    "glycan_library": "glycans",
                    "glycans": {"G42666HT": ["A:S13"]},
                }
            },
        ],
    }
    first = mp.load_protein(TRPCAGE)
    mp.apply_spec(first, spec, base_dir=library.parent, relax="none")

    written = yaml.safe_load(mp.to_yaml(first))
    assert [entry["attach"]["name"] for entry in written["modifications"][1:]] == ["M001", "G001"]
    second = mp.load_protein(TRPCAGE)
    mp.apply_spec(second, written, relax="none")

    assert second.atom_keys() == first.atom_keys()
    assert mp.modification_sites(second) == mp.modification_sites(first)
    assert mp.segments(second) == mp.segments(first)


def test_to_spec_refuses_a_fragment_it_cannot_write() -> None:
    from mbuild.biopolymers import prepare_fragment

    protein = mp.load_protein(TRPCAGE)
    acetyl = prepare_fragment("*C(=O)C", "ACY")  # built with mBuild directly
    mp.attach(protein, acetyl, type="generic", site="A:S13", atom_name="OG", relax=False)
    with pytest.raises(ValueError, match="built outside"):
        mp.to_spec(protein)
