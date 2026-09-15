"""Tests for the catalytic triad contract plugin."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.catalytic_triad import (
    SIMULTANEOUS_CONTACT,
    THRESHOLD_OPERATOR,
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
        n_segments=1,
        atom_resindex=[0, 1, 1, 2],
        residue_segindex=[0, 0, 0],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["OG", "NE2", "ND1", "OD2"])
    universe.add_TopologyAttr("resname", ["SER", "HIS", "ASP"])
    universe.add_TopologyAttr("resid", [77, 156, 133])
    universe.add_TopologyAttr("segid", ["A"])
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


def test_a_pair_exactly_at_the_cutoff_is_not_in_contact() -> None:
    """The cutoff comparison is strictly less than, in both places it is used."""

    universe = _universe([3.4, 3.5, 3.5], [3.4, 3.4, 3.5])

    observables = CatalyticTriad().compute(universe, _frames(), SETTINGS)
    by_name = {observable.name: observable for observable in observables}

    np.testing.assert_allclose(by_name["Ser-His within 3.5 A"].values, [1.0, 0.0, 0.0])
    np.testing.assert_allclose(by_name["His-Asp within 3.5 A"].values, [1.0, 1.0, 0.0])
    np.testing.assert_allclose(by_name[SIMULTANEOUS_CONTACT].values, [1.0, 0.0, 0.0])
    assert by_name[SIMULTANEOUS_CONTACT].metadata["threshold_operator"] == THRESHOLD_OPERATOR


def test_the_composite_fraction_is_tested_and_the_per_pair_ones_are_not() -> None:
    """Each pair fraction repeats its own distance; the composite one does not."""

    observables = CatalyticTriad().compute(_universe([3.0, 3.0], [3.0, 3.0]), _frames(), SETTINGS)
    tested = {observable.name: observable.tested for observable in observables}

    assert tested == {
        "Ser-His": True,
        "Ser-His within 3.5 A": False,
        "His-Asp": True,
        "His-Asp within 3.5 A": False,
        SIMULTANEOUS_CONTACT: True,
    }


def test_observables_carry_the_active_site_identity() -> None:
    """The name and description of the site travel with every observable."""

    settings = SETTINGS.model_copy(update={"description": "Ser-His-Asp relay"})

    observables = CatalyticTriad().compute(_universe([3.0], [3.0]), _frames(), settings)

    assert observables[0].metadata["active_site"] == "LipA Catalytic Triad"
    assert observables[0].metadata["description"] == "Ser-His-Asp relay"
    assert observables[-1].metadata["pairs"] == ["Ser-His", "His-Asp"]


def test_a_misspelled_setting_is_named_rather_than_absorbed() -> None:
    """An unknown key such as 'cutoff' warns instead of silently doing nothing."""

    with pytest.warns(UserWarning, match="cutoff"):
        CatalyticTriadSettings(
            pairs=[{"label": "a", "selection_a": "name OG", "selection_b": "name NE2"}],
            cutoff=3.5,
        )
