"""Tests for the distances contract plugin."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.distances import Distances, DistancesAnalysis, DistancesSettings
from polyzymd.analyses.mda import FrameSelection

mda = pytest.importorskip("MDAnalysis")

BOX = [40.0, 40.0, 40.0, 90.0, 90.0, 90.0]


def _universe(separation: float, n_frames: int = 4) -> Any:
    """Two atoms held ``separation`` angstrom apart for every frame."""

    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(
        2,
        n_residues=2,
        atom_resindex=[0, 1],
        residue_segindex=[0, 1],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["OG", "C13x"])
    universe.add_TopologyAttr("resname", ["SER", "RBY"])
    universe.add_TopologyAttr("resid", [1, 2])
    universe.add_TopologyAttr("masses", [16.0, 12.0])
    positions = np.asarray([[[0.0, 0.0, 0.0], [separation, 0.0, 0.0]]] * n_frames, dtype=np.float32)
    universe.load_new(positions, format=MemoryReader)
    for timestep in universe.trajectory:
        timestep.dimensions = BOX
    return universe


def _settings(**overrides: Any) -> DistancesSettings:
    """Settings for the one pair the fixtures measure."""

    payload: dict[str, Any] = {
        "pairs": [{"label": "Ser-Substrate", "selection_a": "name OG", "selection_b": "name C13x"}],
        "threshold": 3.5,
    }
    payload.update(overrides)
    return DistancesSettings(**payload)


def _frames() -> FrameSelection:
    """Frame selection covering the whole trajectory."""

    return FrameSelection(start=0, stop=None, step=1, timestep_ps=1.0)


def test_each_pair_reports_a_distance_and_a_contact_fraction() -> None:
    """A pair with a threshold yields a distance in angstrom and a fraction."""

    observables = Distances().compute(_universe(3.0), _frames(), _settings())

    assert [(obs.name, obs.kind, obs.unit) for obs in observables] == [
        ("Ser-Substrate", "mean_of_timeseries", "A"),
        ("Ser-Substrate below 3.5 A", "fraction", "fraction"),
    ]
    np.testing.assert_allclose(observables[0].values, [3.0] * 4, atol=1e-6)
    np.testing.assert_allclose(observables[1].values, [1.0] * 4)


def test_a_pair_without_any_threshold_reports_only_the_distance() -> None:
    """Dropping the global threshold drops the contact fraction with it."""

    observables = Distances().compute(_universe(3.0), _frames(), _settings(threshold=None))

    assert [obs.name for obs in observables] == ["Ser-Substrate"]


def test_the_pair_threshold_and_its_state_label_name_the_fraction() -> None:
    """A per-pair threshold overrides the global one and below_label names the state."""

    settings = _settings(
        pairs=[
            {
                "label": "Ser-Substrate",
                "selection_a": "name OG",
                "selection_b": "name C13x",
                "threshold": 10.0,
                "below_label": "Within 10 Angstrom",
            }
        ]
    )

    observables = Distances().compute(_universe(5.0), _frames(), settings)

    assert observables[1].name == "Ser-Substrate Within 10 Angstrom"
    np.testing.assert_allclose(observables[1].values, [1.0] * 4)


def test_deprecated_alignment_settings_warn_and_change_nothing() -> None:
    """Alignment settings are accepted for one release and ignored."""

    with pytest.warns(DeprecationWarning, match="align_trajectory"):
        settings = _settings(align_trajectory=True, alignment_mode="frame", alignment_frame=0)

    observables = Distances().compute(_universe(3.0), _frames(), settings)

    np.testing.assert_allclose(observables[0].values, [3.0] * 4, atol=1e-6)


def test_the_condition_aggregate_reports_replicate_level_uncertainty(
    run_contract_analysis: Any,
) -> None:
    """Three replicates give a mean, a SEM and an interval over replicates."""

    separations = {1: 3.0, 2: 4.0, 3: 5.0}
    artifact = run_contract_analysis(
        DistancesAnalysis,
        _settings(),
        lambda replicate: _universe(separations[replicate]),
    )
    aggregates = {
        payload["name"]: ObservableAggregate.model_validate(payload)
        for payload in artifact.payload["observables"]
    }

    distance = aggregates["Ser-Substrate"]
    assert distance.replicate_values == pytest.approx([3.0, 4.0, 5.0], abs=1e-6)
    assert distance.mean == pytest.approx(4.0, abs=1e-6)
    assert distance.sem == pytest.approx(1.0 / np.sqrt(3.0), abs=1e-6)
    assert distance.unit == "A" and distance.n_replicates == 3
    assert aggregates["Ser-Substrate below 3.5 A"].mean == pytest.approx(1.0 / 3.0)
