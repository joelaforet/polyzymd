"""Parity of the ported distance plugins against the packages they replaced.

The reference files were frozen from the old ``distances`` and
``catalytic_triad`` packages on real trajectory data before those packages were
deleted, using the campaign settings in the LipA_363K_REDO comparison file. Both
plugins now reach the same numbers through
:func:`polyzymd.analyses.mda.pair_distance.pair_distance_matrix`, which issues
the same ``calc_bonds`` call on the same coordinates, so the tolerance is
floating-point identity rather than a physical one. Alignment was already
removed from the base branch, so no alignment-related difference is expected.

The data lives on a local copy of the cluster trajectories. Without it the
tests skip, so a checkout with no mount still passes.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.catalytic_triad import (
    SIMULTANEOUS_CONTACT,
    CatalyticTriad,
    CatalyticTriadSettings,
)
from polyzymd.analyses.distances import DEPRECATED_KEYS, Distances, DistancesSettings
from polyzymd.analyses.mda import FrameSelection

mda = pytest.importorskip("MDAnalysis")

REFERENCES = Path(__file__).parent
ATOL = 1e-8


def _reference(name: str) -> dict[str, Any]:
    """Load one frozen reference file."""

    return json.loads((REFERENCES / f"{name}_reference.json").read_text(encoding="utf-8"))


def _settings_payload(reference: dict[str, Any]) -> dict[str, Any]:
    """Frozen settings without the keys this port dropped."""

    def strip(payload: dict[str, Any]) -> dict[str, Any]:
        kept = {key: value for key, value in payload.items() if key not in DEPRECATED_KEYS}
        if "pairs" in kept:
            kept["pairs"] = [strip(pair) for pair in kept["pairs"]]
        return kept

    return strip(reference["settings"])


def _universe(reference: dict[str, Any]) -> Any:
    """Rebuild the universe the reference was frozen from, or skip."""

    inputs = reference["inputs"]
    topology, trajectory = Path(inputs["topology"]), Path(inputs["trajectory"])
    if not topology.exists() or not trajectory.exists():
        pytest.skip(f"real trajectory data is not present at {trajectory}")
    return mda.Universe(str(topology), str(trajectory))


def _frames(reference: dict[str, Any]) -> FrameSelection:
    """Frame window the reference was frozen over."""

    window = reference["inputs"]["frames"]
    return FrameSelection(
        start=window["start"], stop=window["stop"], step=window["step"], timestep_ps=40.0
    )


def _series(observables: Any, name: str) -> np.ndarray:
    """Per-frame values of one observable."""

    return np.asarray(next(obs for obs in observables if obs.name == name).values, dtype=np.float64)


def test_distances_match_the_replaced_package() -> None:
    """Every configured pair reproduces the frozen per-frame distances."""

    reference = _reference("distances")
    settings = DistancesSettings(**_settings_payload(reference))
    observables = Distances().compute(_universe(reference), _frames(reference), settings)

    for pair in reference["pairs"]:
        expected = np.asarray(pair["distances"], dtype=np.float64)
        measured = _series(observables, pair["label"])
        np.testing.assert_allclose(measured, expected, rtol=0, atol=ATOL)
        assert float(np.mean(measured)) == pytest.approx(pair["mean_distance"], abs=ATOL)


def test_distance_contact_fractions_match_the_replaced_package() -> None:
    """The fraction below threshold is the same, now stored as a fraction."""

    reference = _reference("distances")
    settings = DistancesSettings(**_settings_payload(reference))
    observables = Distances().compute(_universe(reference), _frames(reference), settings)

    for pair, configured in zip(reference["pairs"], settings.pairs, strict=True):
        threshold = pair["threshold"]
        state = configured.below_label or f"below {float(threshold):g} A"
        measured = _series(observables, f"{pair['label']} {state}")
        assert float(np.mean(measured)) == pytest.approx(pair["fraction_below_threshold"], abs=ATOL)


def test_catalytic_triad_matches_the_replaced_package() -> None:
    """Triad distances and the simultaneous contact indicator are unchanged."""

    reference = _reference("catalytic_triad")
    settings = CatalyticTriadSettings(**_settings_payload(reference))
    observables = CatalyticTriad().compute(_universe(reference), _frames(reference), settings)

    for pair in reference["pairs"]:
        expected = np.asarray(pair["distances"], dtype=np.float64)
        np.testing.assert_allclose(_series(observables, pair["label"]), expected, rtol=0, atol=ATOL)
        within = _series(observables, f"{pair['label']} within {settings.threshold:g} A")
        assert float(np.mean(within)) == pytest.approx(pair["fraction_below_threshold"], abs=ATOL)

    simultaneous = _series(observables, SIMULTANEOUS_CONTACT)
    np.testing.assert_array_equal(
        simultaneous, np.asarray(reference["simultaneous_contact"], dtype=np.float64)
    )
    assert float(np.mean(simultaneous)) == pytest.approx(
        reference["simultaneous_contact_fraction"], abs=ATOL
    )
