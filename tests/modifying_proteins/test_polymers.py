"""Polymer chains grafted onto a protein: drawing, building and attaching them."""

from pathlib import Path

import pytest

from polyzymd import modifying_proteins as mp
from polyzymd.modifying_proteins.polymers import anchor_residue, draw_chain

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
#: NHS ester aminolysis: the amine bonds the ester carbon, N-hydroxysuccinimide leaves.
NHS_AMINOLYSIS = "[N:1]([H:2])[H].[C:3](=[O:4])[O][N]>>[N:1]([H:2])[C:3]=[O:4]"
CAP = [{"label": "H", "smiles": "[*][H]", "probability": 1.0}]


def chain(**extra) -> dict:
    """A 7-position SBMA chain (caps at 1 and 7) with the NHS linker pinned in the middle."""
    polymer = {
        "monomers": [
            {
                "label": "A",
                "name": "SBMA",
                "residue_name": "SBM",
                "smiles": "[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]",
                "probability": 1.0,
            },
            {
                "label": "N",
                "name": "NHS",
                "residue_name": "NHS",
                "smiles": "[*]CC([*])(C)C(=O)ON1C(=O)CCC1=O",
                "probability": 0.0,
            },
        ],
        "end_groups": {"head": CAP, "tail": CAP},
        "length": 7,
        "sequence": {"fixed": {4: "NHS"}},
        "random_seed": 1,
        "force_field": ["openff-2.3.0.offxml"],
        **extra,
    }
    return {key: value for key, value in polymer.items() if value is not None or key != "sequence"}


def test_the_drawn_chain_keeps_its_pinned_unit_and_caps() -> None:
    assert draw_chain(chain()) == "HAANAAH"
    assert draw_chain({**chain(), "drawn": "HNAAAAH"}) == "HNAAAAH"


def test_drawn_positions_need_a_replicate_or_a_seed() -> None:
    with pytest.raises(ValueError, match="needs a seed: the replicate number"):
        draw_chain(chain(random_seed=None))
    assert draw_chain(chain(random_seed=None, replicate=2)) == "HAANAAH"
    assert draw_chain(
        chain(random_seed=None, sequence={"mode": "explicit", "units": list("HAANAAH")})
    )


def test_attach_at_counts_the_whole_chain() -> None:
    assert anchor_residue(chain(), 4) == 3
    capped = chain()
    capped["end_groups"] = {"head": [{**CAP[0], "residue_name": "HCP"}], "tail": CAP}
    assert anchor_residue(capped, 4) == 4
    with pytest.raises(ValueError, match="outside the chain's positions 1..7"):
        anchor_residue(chain(), 8)
    with pytest.raises(ValueError, match="attach_at 1 is the head cap 'H', not a monomer"):
        anchor_residue(chain(), 1)


def test_a_cap_that_is_a_monomer_end_unit_can_be_attached_to() -> None:
    polymer = chain()
    polymer["end_groups"] = {
        "head": [{"label": "A", "smiles": "[*]C", "probability": 1.0}],
        "tail": CAP,
    }
    assert anchor_residue(polymer, 1) == 1


def test_a_grafted_chain_is_one_polymer_site_and_replays_to_the_same_atoms() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.deprotonate(protein, "A:K8", "NZ")
    record = mp.attach(
        protein,
        {"polymer": chain()},
        type="polymer",
        site="A:K8",
        atom_name="NZ",
        reaction=NHS_AMINOLYSIS,
        attach_at=4,
        relax=False,
    )
    assert record.name == "P001"
    assert record.arguments["attach_at"] == 4
    assert record.arguments["fragment"]["polymer"]["drawn"] == draw_chain(chain(), site="A:K8")
    assert [r[3] for r in record.new_residues] == ["SBM", "SBM", "NHS", "SBM", "SBM"]
    # NZ loses its second proton; the ester oxygen and the succinimide (12 atoms) leave.
    assert [key.atom_name for key in record.removed] == ["HZ2"]
    (site,) = mp.modification_sites(protein)
    assert site.types == frozenset({"polymer"})
    segments = mp.segments(protein)
    assert {segments[i] for i in site.atom_indices} >= {"P001"}

    replay = mp.load_protein(TRPCAGE)
    mp.apply_spec(replay, mp.to_spec(protein), base_dir=ROOT, relax="none")
    assert replay.atom_keys() == protein.atom_keys()


def _mixed(**extra) -> dict:
    """A 12-unit chain drawn half SBMA, half NHS-methacrylate, nothing pinned."""
    polymer = chain(length=14, sequence=None, **extra)
    polymer["monomers"][0]["probability"] = 0.5
    polymer["monomers"][1]["probability"] = 0.5
    return polymer


def test_each_site_of_an_entry_draws_its_own_chain() -> None:
    sites = ["A:K8", "A:K10", "A:K12"]
    drawn = [draw_chain(_mixed(), site=site) for site in sites]
    assert len(set(drawn)) == 3
    assert drawn == [draw_chain(_mixed(), site=site) for site in sites]
    assert draw_chain({**_mixed(), "drawn": drawn[0]}, site="A:K10") == drawn[0]


def test_a_graft_and_a_soluble_chain_draw_alike_from_one_seed() -> None:
    import random

    from polyzymd.builders.polymer import PolymerBuilder

    builder = PolymerBuilder(
        characters=["A", "N"], probabilities=[0.5, 0.5], length=14, type_prefix="TP"
    )
    builder._end_groups = {"head": CAP, "tail": CAP}
    random.seed(7)
    assert draw_chain(_mixed(random_seed=7)) == builder._draw_sequence()


def test_grafted_chains_are_drawn_anew_for_each_replicate() -> None:
    first, second = (
        draw_chain(_mixed(random_seed=None, replicate=n), site="A:K8") for n in (1, 2)
    )
    assert first != second
    assert first == draw_chain(_mixed(random_seed=None, replicate=1), site="A:K8")


def test_a_random_seed_gives_one_chain_for_every_replicate() -> None:
    first, second = (
        draw_chain(_mixed(random_seed=5, replicate=n), site="A:K8") for n in (1, 2)
    )
    assert first == second
