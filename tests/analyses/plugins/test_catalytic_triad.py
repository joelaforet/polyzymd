"""Tests for the catalytic triad contract plugin."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.catalytic_triad import (
    SIMULTANEOUS_CONTACT,
    CatalyticTriad,
    CatalyticTriadAnalysis,
    CatalyticTriadSettings,
)
from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.mda import FrameSelection

mda = pytest.importorskip("MDAnalysis")

BOX = [40.0, 40.0, 40.0, 90.0, 90.0, 90.0]

SETTINGS = CatalyticTriadSettings(
    name="LipA Catalytic Triad",
    threshold=3.5,
    pairs=[
        {"label": "Ser-His", "selection_a": "name OG", "selection_b": "name NE2"},
        {"label": "His-Asp", "selection_a": "name ND1", "selection_b": "name OD2"},
    ],
)


def _universe(ser_his: list[float], his_asp: list[float]) -> Any:
    """Four atoms whose two pair distances follow the given per-frame series."""

    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(
        4,
        n_residues=3,
        atom_resindex=[0, 1, 1, 2],
        residue_segindex=[0, 0, 0],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["OG", "NE2", "ND1", "OD2"])
    universe.add_TopologyAttr("resname", ["SER", "HIS", "ASP"])
    universe.add_TopologyAttr("resid", [77, 156, 133])
    universe.add_TopologyAttr("masses", [16.0, 14.0, 14.0, 16.0])
    frames = [
        [[0.0, 0.0, 0.0], [first, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 10.0 + second, 0.0]]
        for first, second in zip(ser_his, his_asp, strict=True)
    ]
    universe.load_new(np.asarray(frames, dtype=np.float32), format=MemoryReader)
    for timestep in universe.trajectory:
        timestep.dimensions = BOX
    return universe


def _frames() -> FrameSelection:
    """Frame selection covering the whole trajectory."""

    return FrameSelection(start=0, stop=None, step=1, timestep_ps=1.0)


def test_every_pair_reports_a_distance_and_a_contact_fraction() -> None:
    """Each pair gives a distance in angstrom and a fraction within the cutoff."""

    observables = CatalyticTriad().compute(
        _universe([3.0, 3.0, 4.0, 4.0], [3.0, 4.0, 3.0, 4.0]), _frames(), SETTINGS
    )

    assert [(obs.name, obs.kind, obs.unit) for obs in observables] == [
        ("Ser-His", "mean_of_timeseries", "A"),
        ("Ser-His within 3.5 A", "fraction", "fraction"),
        ("His-Asp", "mean_of_timeseries", "A"),
        ("His-Asp within 3.5 A", "fraction", "fraction"),
        (SIMULTANEOUS_CONTACT, "fraction", "fraction"),
    ]
    np.testing.assert_allclose(observables[0].values, [3.0, 3.0, 4.0, 4.0], atol=1e-6)
    np.testing.assert_allclose(observables[1].values, [1.0, 1.0, 0.0, 0.0])
    np.testing.assert_allclose(observables[3].values, [1.0, 0.0, 1.0, 0.0])


def test_simultaneous_contact_counts_only_frames_with_every_pair_inside() -> None:
    """The composite fraction is the AND of the per-pair indicators, not their mean."""

    observables = CatalyticTriad().compute(
        _universe([3.0, 3.0, 4.0, 4.0], [3.0, 4.0, 3.0, 4.0]), _frames(), SETTINGS
    )
    simultaneous = next(obs for obs in observables if obs.name == SIMULTANEOUS_CONTACT)

    np.testing.assert_allclose(simultaneous.values, [1.0, 0.0, 0.0, 0.0])


def test_fractions_are_stored_as_fractions_over_replicates(run_contract_analysis: Any) -> None:
    """The aggregate reports a fraction in [0, 1] with unit 'fraction', not a percent."""

    universes = {
        1: _universe([3.0, 3.0], [3.0, 3.0]),
        2: _universe([3.0, 4.0], [3.0, 3.0]),
        3: _universe([4.0, 4.0], [3.0, 3.0]),
    }
    artifact = run_contract_analysis(
        CatalyticTriadAnalysis, SETTINGS, lambda replicate: universes[replicate]
    )
    aggregates = {
        payload["name"]: ObservableAggregate.model_validate(payload)
        for payload in artifact.payload["observables"]
    }

    simultaneous = aggregates[SIMULTANEOUS_CONTACT]
    assert simultaneous.unit == "fraction"
    assert simultaneous.replicate_values == pytest.approx([1.0, 0.5, 0.0])
    assert simultaneous.mean == pytest.approx(0.5)
    assert aggregates["Ser-His"].unit == "A"
