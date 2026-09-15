"""Parity of the contract hydrogen-bond plugin with the package it replaced.

``hydrogen_bonds_reference.json`` was frozen from the pre-port
``polyzymd.analyses.hydrogen_bonds._mda`` kernel, on frames 500 to 600 of the
50:50 SBMA-EGMA LipA trajectory, before that package was deleted. It records
the topology and trajectory paths, the frame window, the campaign settings and
the commit it came from.

The old kernel on the base branch already restricted donors and acceptors to
nitrogen and oxygen, the correction that PR 102 made, so this is parity against
corrected behaviour and not against the carbon-donor counts the campaign
artifacts hold. Both implementations call the same
``MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis`` with the same
selections and cutoffs, so counts are integers that must agree exactly and
occupancies are ratios of integers that agree to floating-point noise.

Pairs that share an occupancy are ordered differently by the two
implementations, because the old kernel left ties in the order the pairs were
first seen and the port breaks them on the pair label. The occupancy vector is
therefore compared rank by rank, and each label the port reports is checked
against the occupancy the old kernel recorded for that same pair.

The tests skip when the trajectory is not on this machine, so a checkout
without the data still passes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.hydrogen_bonds import HydrogenBonds, HydrogenBondSettings

REFERENCE = Path(__file__).with_name("hydrogen_bonds_reference.json")
COUNT_TOLERANCE = 1e-12
OCCUPANCY_TOLERANCE = 1e-10


def _reference() -> dict[str, Any]:
    """Frozen values and the provenance block that describes them."""
    return json.loads(REFERENCE.read_text())


def _computed(reference: dict[str, Any]) -> tuple[dict[str, Any], np.ndarray]:
    """Run the port on the frozen window, keyed by observable name."""
    import MDAnalysis as mda

    from polyzymd.analyses.mda.frame_selection import FrameSelection

    topology, trajectory = Path(reference["topology"]), Path(reference["trajectory"])
    if not (topology.exists() and trajectory.exists()):
        pytest.skip(f"real trajectory not on this machine: {trajectory}")
    digest = hashlib.md5(topology.read_bytes()).hexdigest()
    assert digest == reference["topology_md5"], (
        f"{topology} is not the file the reference was frozen from "
        f"(md5 {digest}, expected {reference['topology_md5']})"
    )
    universe = mda.Universe(str(topology), str(trajectory))
    assert universe.atoms.n_atoms == reference["n_atoms"]
    observables, sidecars = HydrogenBonds().compute(
        universe,
        FrameSelection(**reference["frames"]),
        HydrogenBondSettings.model_validate(reference["settings"]),
    )
    return {observable.name: observable for observable in observables}, sidecars[
        "hydrogen_bond_events"
    ]


def test_counts_per_frame_match_the_frozen_reference() -> None:
    """Every partition reproduces the old kernel's bond count in every frame."""
    reference = _reference()
    observables, _ = _computed(reference)

    for summary, counts in reference["counts_per_frame"].items():
        observable = observables[f"hbonds_{summary}"]
        assert observable.kind == "mean_of_timeseries"
        assert observable.unit == "count"
        np.testing.assert_allclose(observable.values, counts, rtol=0, atol=COUNT_TOLERANCE)


def test_pair_occupancy_profiles_match_the_frozen_reference() -> None:
    """Ranked occupancies and the pairs behind them reproduce the old kernel."""
    reference = _reference()
    observables, _ = _computed(reference)
    top_n = reference["settings"]["top_n_pairs"]

    for summary, top_pairs in reference["top_pairs"].items():
        observable = observables[f"pair_occupancy_{summary}"]
        assert observable.kind == "profile"
        assert observable.index == [float(rank) for rank in range(top_n)]
        expected = [pair["occupancy"] for pair in top_pairs] + [0.0] * (top_n - len(top_pairs))
        np.testing.assert_allclose(observable.values, expected, rtol=0, atol=OCCUPANCY_TOLERANCE)
        assert observable.index_label == "occupancy rank"
        occupancy_of = reference["pair_occupancy"][summary]
        for label, value in zip(observable.metadata["pair_labels"], observable.values):
            if label:
                assert occupancy_of[label] == pytest.approx(value, abs=OCCUPANCY_TOLERANCE)


def test_event_sidecar_matches_the_frozen_event_count() -> None:
    """The raw event table the sidecar carries is the one MDAnalysis produced."""
    reference = _reference()
    _, events = _computed(reference)

    assert events.shape == (reference["n_events"], 6)
    assert np.all(events[:, 4] <= reference["settings"]["distance_cutoff"])
    assert np.all(events[:, 5] >= reference["settings"]["angle_cutoff"])


def test_reference_records_its_provenance() -> None:
    """The frozen file names the code, the window and the selections it came from."""
    reference = _reference()

    assert reference["polyzymd_commit"] == "6156c87e21e8b76f40c81378de9a5ccd9a142dcd"
    assert len(reference["topology_md5"]) == 32
    assert reference["frames"] == {"start": 500, "stop": 600, "step": 1}
    assert reference["donors_selection"].endswith("and element N O")
    assert reference["settings"]["donor_acceptor_elements"] == ["N", "O"]
