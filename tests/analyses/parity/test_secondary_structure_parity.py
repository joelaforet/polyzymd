"""Parity of the ported secondary_structure plugin against the old implementation.

The reference file was frozen from the pre-port package
``src/polyzymd/analyses/secondary_structure/_mda.py`` on 100 frames of a real
trajectory, before that package was deleted. Both implementations call
``mdtraj.compute_dssp(simplified=True)`` on a topology built the same way from
the same MDAnalysis selection, so helix and strand must agree to floating-point
noise and the tolerance is 1e-12.

Coil is the one deliberate difference. The old encoder started from an
all-coil matrix and overwrote only H and E, so every residue mdtraj returned as
``NA`` was counted as coil. The port gives those residues their own
``ss_unassigned`` fraction, which makes the parity statement
``old_coil == new_coil + new_unassigned``. On this window mdtraj assigns every
residue, so the unassigned fraction is zero and coil matches directly; the test
asserts the sum anyway so it keeps holding on a window where it is not.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.secondary_structure import (
    SecondaryStructure,
    SecondaryStructureSettings,
)

REFERENCE = Path(__file__).with_name("secondary_structure_reference.json")
ATOL = 1e-12
DSSP_NAMES = ("helix", "strand", "coil", "unassigned")


@pytest.fixture(scope="module")
def reference() -> dict:
    """Frozen output of the pre-port implementation."""
    return json.loads(REFERENCE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def observables(reference: dict) -> dict:
    """Run the ported plugin over the frozen window, keyed by observable name."""
    mda = pytest.importorskip("MDAnalysis")
    pytest.importorskip("mdtraj")
    topology = Path(reference["topology"])
    trajectory = Path(reference["trajectory"])
    if not (topology.exists() and trajectory.exists()):
        pytest.skip(
            f"real trajectory data is not on this machine: {trajectory} is missing. "
            "Copy the run from the cluster to run the parity check."
        )
    digest = hashlib.md5(topology.read_bytes()).hexdigest()  # noqa: S324
    assert digest == reference["topology_md5"], (
        f"the topology at {topology} is not the one the reference was frozen from "
        f"({digest} instead of {reference['topology_md5']}); the parity numbers do not apply"
    )
    universe = mda.Universe(str(topology), str(trajectory))
    window = reference["frames"]
    frames = FrameSelection(start=window["start"], stop=window["stop"], step=window["step"])
    computed = SecondaryStructure().compute(
        universe, frames, SecondaryStructureSettings(**reference["settings"])
    )
    return {observable.name: observable for observable in computed}


def test_per_frame_helix_and_strand_match(observables: dict, reference: dict) -> None:
    """Helix and strand fractions per frame reproduce the old values exactly."""
    for label in ("helix", "strand"):
        np.testing.assert_allclose(
            observables[f"ss_{label}"].values,
            reference["old_per_frame_fractions"][label],
            rtol=0,
            atol=ATOL,
            err_msg=f"{label} fraction per frame drifted from the pre-port implementation",
        )


def test_coil_matches_once_unassigned_is_added_back(observables: dict, reference: dict) -> None:
    """Old coil equals new coil plus the residues the old code miscounted as coil."""
    recombined = np.asarray(observables["ss_coil"].values) + np.asarray(
        observables["ss_unassigned"].values
    )
    np.testing.assert_allclose(
        recombined, reference["old_per_frame_fractions"]["coil"], rtol=0, atol=ATOL
    )
    np.testing.assert_allclose(
        observables["ss_unassigned"].values,
        reference["unassigned_per_frame_fraction"],
        rtol=0,
        atol=ATOL,
    )


def test_window_means_match_old_overall_fractions(observables: dict, reference: dict) -> None:
    """The mean over frames reproduces the old whole-window fractions."""
    old = reference["old_overall_fractions"]
    for label in ("helix", "strand"):
        assert float(np.mean(observables[f"ss_{label}"].values)) == pytest.approx(
            old[label], abs=ATOL
        )
    coil = float(np.mean(observables["ss_coil"].values)) + float(
        np.mean(observables["ss_unassigned"].values)
    )
    assert coil == pytest.approx(old["coil"], abs=ATOL)
    assert sum(float(np.mean(observables[f"ss_{name}"].values)) for name in DSSP_NAMES) == (
        pytest.approx(1.0, abs=ATOL)
    )


def test_occupancy_profiles_match_old_persistence(observables: dict, reference: dict) -> None:
    """Per-residue occupancy reproduces the old per-residue persistence arrays."""
    for label in ("helix", "strand"):
        profile = observables[f"{label}_occupancy"]
        np.testing.assert_allclose(
            profile.values,
            reference[f"old_persistence_{label}"],
            rtol=0,
            atol=ATOL,
            err_msg=f"{label} occupancy profile drifted from the pre-port implementation",
        )
        assert profile.index == [float(value) for value in reference["residue_ids"]]
