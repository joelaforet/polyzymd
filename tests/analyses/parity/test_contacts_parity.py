"""Parity of the contract contacts plugin with the implementation it replaced.

``contacts_reference.json`` was frozen from the pre-port
``polyzymd.analyses.contacts._mda`` event detector, on frames 500 to 600 of the
50:50 SBMA-EGMA LipA trajectory, before that code was deleted. The file records
the topology and trajectory paths, the frame window, the settings, the commit
and the MDAnalysis version it came from.

The per-frame contact count, the per-frame coverage and the per-residue contact
fraction come from the same ``capped_distance`` call on the same frames, so they
must agree bit for bit and the tolerance is 1e-12. The window is frames 500 to
600 of the 50:50 run with the settings the LipA_363K_REDO campaign uses, a 4.0 A
cutoff between ``protein`` and ``resname SBM EGM``. The per-residue mean
residence time is the mean of the same event durations in a different summation
order, so it is checked to 1e-12 as well, which is far tighter than one frame of
the 40 ps time axis. The residence-time distribution is a histogram the old code
did not build, so the reference bins the old event durations with the same edges
the plugin uses.

The old plugin reported two replicate metrics. Its ``mean_contact_fraction`` is
the mean of the contact fraction profile, which is also the mean of
``coverage_per_frame``, and the test checks it both ways. Its ``coverage`` was
the share of residues touched at any point in the window, which the port
reports as ``coverage_any_frame``.

The tests skip when the trajectory is not on this machine, so a checkout without
the data still passes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.contacts import Contacts, ContactsSettings

REFERENCE = Path(__file__).with_name("contacts_reference.json")
TOLERANCE = 1e-12


def _reference() -> dict[str, Any]:
    """Frozen values and the provenance that describes them."""
    return json.loads(REFERENCE.read_text())


def _observables(reference: dict[str, Any]) -> tuple[dict[str, Any], np.ndarray]:
    """Run the port on the frozen window, or skip when the data is absent."""
    import MDAnalysis as mda

    from polyzymd.analyses.mda.frame_selection import FrameSelection

    topology, trajectory = Path(reference["topology"]), Path(reference["trajectory"])
    if not (topology.exists() and trajectory.exists()):
        pytest.skip(f"real trajectory not on this machine: {trajectory}")
    window = reference["frames"]
    observables, extras = Contacts().compute(
        mda.Universe(str(topology), str(trajectory)),
        FrameSelection(start=window["start"], stop=window["stop"], step=window["step"]),
        ContactsSettings(
            protein_selection=reference["settings"]["protein_selection"],
            polymer_selection=reference["settings"]["polymer_selection"],
            cutoff=reference["settings"]["cutoff"],
            residence_time_edges_ns=reference["residence_time_edges_ns"],
        ),
    )
    return {observable.name: observable for observable in observables}, extras["contact_events"]


def test_per_frame_series_match_the_old_detector() -> None:
    """Contact count and coverage reproduce the old per-frame contact sets."""
    reference = _reference()
    observables, _ = _observables(reference)

    counts = observables["contact_count"]
    coverage = observables["coverage_per_frame"]
    assert (counts.unit, coverage.unit) == ("count", "fraction")
    assert len(counts.values) == reference["frames"]["n_frames_used"]
    np.testing.assert_allclose(
        counts.values, reference["per_frame_contact_count"], rtol=0, atol=TOLERANCE
    )
    np.testing.assert_allclose(
        coverage.values, reference["per_frame_coverage"], rtol=0, atol=TOLERANCE
    )


def test_contact_fraction_profile_matches_the_old_replicate_metrics() -> None:
    """The profile, its mean and the old window coverage all reproduce."""
    reference = _reference()
    observables, _ = _observables(reference)

    profile = observables["contact_fraction"]
    assert profile.index == [float(resid) for resid in reference["protein_resids"]]
    np.testing.assert_allclose(
        profile.values, reference["contact_fraction"], rtol=0, atol=TOLERANCE
    )
    assert float(np.mean(profile.values)) == pytest.approx(
        reference["mean_contact_fraction"], abs=TOLERANCE
    )
    assert float(np.mean(observables["coverage_per_frame"].values)) == pytest.approx(
        reference["mean_contact_fraction"], abs=TOLERANCE
    )
    any_frame = observables["coverage_any_frame"]
    assert any_frame.values == pytest.approx([reference["coverage_any_frame"]], abs=TOLERANCE)
    assert any_frame.tested is False


def test_residence_times_match_the_old_event_durations() -> None:
    """The distribution and the per-residue means come from the same events."""
    reference = _reference()
    observables, events = _observables(reference)

    distribution = observables["residence_time_distribution"]
    assert distribution.index == reference["residence_time_edges_ns"][:-1]
    np.testing.assert_allclose(
        distribution.values, reference["residence_time_distribution"], rtol=0, atol=TOLERANCE
    )
    np.testing.assert_allclose(
        observables["mean_residence_time"].values,
        reference["mean_residence_time_ns"],
        rtol=0,
        atol=TOLERANCE,
    )
    assert len(events) == reference["n_contact_events"]


def test_event_sidecar_names_residue_chain_and_frames() -> None:
    """Every event row is a protein residue, a polymer chain and a frame span."""
    reference = _reference()
    _, events = _observables(reference)
    window = reference["frames"]

    assert events.shape == (reference["n_contact_events"], 4)
    assert set(events[:, 0].tolist()) <= set(reference["protein_resids"])
    assert events[:, 1].min() >= 0
    assert events[:, 2].min() >= window["start"]
    assert events[:, 3].max() < window["stop"]
    assert np.all(events[:, 3] >= events[:, 2])


def test_reference_records_its_provenance() -> None:
    """The frozen file names the code, the window and the settings it came from."""
    reference = _reference()

    assert "pre-port" in reference["source"]
    assert reference["polyzymd_commit"]
    assert reference["mdanalysis_version"]
    assert reference["settings"]["cutoff"] == 4.0
    assert reference["frames"] == {
        "start": 500,
        "stop": 600,
        "step": 1,
        "n_frames_used": 100,
    }
    topology = Path(reference["topology"])
    if not topology.exists():
        pytest.skip(f"real topology not on this machine: {topology}")
    assert hashlib.md5(topology.read_bytes()).hexdigest() == reference["topology_md5"]
