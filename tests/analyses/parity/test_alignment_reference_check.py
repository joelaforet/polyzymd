"""Real-trajectory check that the alignment reference modes give distinct answers.

The numbers in ``alignment_reference_check.json`` were recorded on one control
replicate of the LipA campaign after the 2026-09-12 fix to the "centroid" and
"frame" reference modes.

Parity tests read trajectories that are not in the repository, so they are off
unless you ask for them. Set ``POLYZYMD_REAL_DATA=1`` to run them, and
``POLYZYMD_REALDATA_ROOT`` to point at the replicate directory if it is not
where the JSON file records it. Both have to be satisfied or the test skips,
which is the normal case in continuous integration.

This test loads one replicate at its full length because ``align_trajectory``
calls ``AlignTraj(in_memory=True)``, which moves the whole trajectory into
memory rather than just the analyzed window. On the recorded dataset that peaks
near 2.5 GB. The behaviour predates this module and is not changed here.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

CHECK_FILE = Path(__file__).with_name("alignment_reference_check.json")
CHECK = json.loads(CHECK_FILE.read_text(encoding="utf-8"))
ROOT = Path(os.environ.get("POLYZYMD_REALDATA_ROOT", CHECK["dataset"]["root"]))
TOPOLOGY = ROOT / CHECK["dataset"]["topology"]
TRAJECTORY = ROOT / CHECK["dataset"]["trajectory"]

REAL_DATA_ENABLED = os.environ.get("POLYZYMD_REAL_DATA") == "1"

pytestmark = pytest.mark.skipif(
    not (REAL_DATA_ENABLED and TOPOLOGY.exists() and TRAJECTORY.exists()),
    reason=(
        "parity tests need POLYZYMD_REAL_DATA=1 and a real trajectory under "
        f"{ROOT} (override with POLYZYMD_REALDATA_ROOT)"
    ),
)


def _mean_rmsf(reference_mode: str, reference_frame: int | None) -> tuple[int | None, float]:
    """Align the production window one way and return the reference and mean RMSF."""

    import MDAnalysis as mda
    from MDAnalysis.analysis import rms

    from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory

    settings = CHECK["settings"]
    universe = mda.Universe(str(TOPOLOGY), str(TRAJECTORY))
    reference_index = align_trajectory(
        universe,
        AlignmentConfig(
            enabled=True,
            reference_mode=reference_mode,
            reference_frame=reference_frame,
            selection=settings["alignment_selection"],
            centroid_selection=settings["centroid_selection"],
        ),
        start_frame=settings["start_frame"],
        stop_frame=settings["stop_frame"],
        step_frame=settings["step_frame"],
    )
    atoms = universe.select_atoms(settings["selection"])
    assert len(atoms) == settings["n_atoms"]
    profile = rms.RMSF(atoms).run(
        start=settings["start_frame"],
        stop=settings["stop_frame"],
        step=settings["step_frame"],
    )
    return reference_index, float(np.asarray(profile.results.rmsf, dtype=float).mean())


def test_recorded_rmsf_means_are_reproduced() -> None:
    """Both reference modes reproduce their recorded mean RMSF and differ from each other."""

    expected = CHECK["after_fix"]
    tolerance = CHECK["tolerance"]

    frame_index, frame_mean = _mean_rmsf("frame", expected["frame"]["reference_frame_1indexed"])
    centroid_index, centroid_mean = _mean_rmsf("centroid", None)

    assert frame_index == expected["frame"]["resolved_reference_frame_0indexed"]
    assert centroid_index == expected["centroid"]["resolved_reference_frame_0indexed"]
    assert frame_mean == pytest.approx(expected["frame"]["mean_rmsf"], abs=tolerance)
    assert centroid_mean == pytest.approx(expected["centroid"]["mean_rmsf"], abs=tolerance)
    assert frame_mean != pytest.approx(centroid_mean, abs=1e-6)
