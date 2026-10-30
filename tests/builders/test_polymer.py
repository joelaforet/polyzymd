"""Sequences of polymer chains: statistical draws, blocks, explicit units and pinned positions."""

from __future__ import annotations

import random

from polyzymd.builders.polymer import PolymerBuilder, draw_units, generate_random_sequence

CAPS = {
    "head": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
    "tail": [{"label": "T", "smiles": "[*][H]", "probability": 1.0}],
}


def _builder(pattern=None, length=9) -> PolymerBuilder:
    builder = PolymerBuilder(
        characters=["A", "B", "N"],
        probabilities=[0.7, 0.3, 0.0],
        length=length,
        type_prefix="TP",
        sequence_pattern=pattern,
    )
    builder._end_groups = CAPS
    return builder


def test_an_open_pattern_draws_what_v1_3_drew() -> None:
    random.seed(4)
    expected = generate_random_sequence(7, ["A", "B"], [0.7, 0.3])
    random.seed(4)
    assert draw_units([None] * 7, ["A", "B"], [0.7, 0.3]) == expected


def test_without_a_pattern_the_chain_sequence_is_unchanged() -> None:
    random.seed(11)
    head = random.choices(["H"], weights=[1.0], k=1)[0]
    middle = generate_random_sequence(7, ["A", "B", "N"], [0.7, 0.3, 0.0])
    tail = random.choices(["T"], weights=[1.0], k=1)[0]
    random.seed(11)
    assert _builder()._draw_sequence() == head + middle + tail


def test_a_pinned_position_keeps_its_monomer_and_the_rest_are_drawn() -> None:
    # Positions count the whole chain: 1 and 9 are the caps, 5 the middle unit.
    pattern = [None, None, None, None, "N", None, None, None, None]
    random.seed(3)
    chains = [_builder(pattern)._draw_sequence() for _ in range(20)]
    assert all(chain[0] == "H" and chain[-1] == "T" for chain in chains)
    assert all(chain[1:-1][3] == "N" for chain in chains)
    assert all(set(chain[1:-1][:3] + chain[1:-1][4:]) <= {"A", "B"} for chain in chains)
    assert len(set(chains)) > 1


def test_a_full_pattern_is_the_same_for_every_chain() -> None:
    pattern = list("HAAANAAAT")
    random.seed(0)
    assert {_builder(pattern)._draw_sequence() for _ in range(5)} == {"HAAANAAAT"}


def test_build_from_config_follows_the_sequence_block(monkeypatch) -> None:
    from polyzymd.config.schema import PolymerSolute

    config = PolymerSolute(
        type="polymer",
        type_prefix="TP",
        generation_mode="dynamic",
        monomers=[
            {"label": "A", "name": "SBMA", "probability": 1.0, "smiles": "[*]CC([*])C"},
            {"label": "N", "name": "NHS", "probability": 0.0, "smiles": "[*]CC([*])CO"},
        ],
        end_groups=CAPS,
        length=9,
        count=3,
        force_field=["openff-2.3.0.offxml"],
        charges={"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"},
        sequence={
            "mode": "block",
            "blocks": [{"H": 1}, {"SBMA": 3}, {"NHS": 1}, {"SBMA": 3}, {"T": 1}],
        },
    )
    builder = PolymerBuilder(characters=["A"], probabilities=[1.0], length=9, type_prefix="TP")
    monkeypatch.setattr(builder, "_get_or_create_molecule", lambda sequence: sequence)
    molecules, counts = builder.build_from_config(config, seed=1)
    assert molecules == ["HAAANAAAT"] and counts == [3]
