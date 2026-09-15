"""RMSF parity against the pre-port implementation on real trajectory data.

``rmsf_reference.json`` holds the per-residue profile the deleted
``rmsf/_mda.py`` produced for frames 500 to 600 of the control run of the LipA
363 K campaign, for three alignment reference modes. The file records the
topology and its md5, the trajectory, the window, the settings and the commit
the old code was exported from. The test skips when that trajectory is not on
this machine, so a checkout without the data still passes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from polyzymd.analyses.mda import FrameSelection
from polyzymd.analyses.rmsf import RMSF, RMSFSettings

REFERENCE_PATH = Path(__file__).parent / "rmsf_reference.json"
REFERENCE = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
MODES = sorted(REFERENCE["modes"])


def _universe():
    """Load the frozen window into memory, or skip when the data is absent."""
    mda = pytest.importorskip("MDAnalysis")
    topology = Path(REFERENCE["topology"])
    trajectory = Path(REFERENCE["trajectory"])
    if not (topology.exists() and trajectory.exists()):
        pytest.skip(f"real trajectory not on this machine: {trajectory}")
    digest = hashlib.md5()
    with topology.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    assert (
        digest.hexdigest() == REFERENCE["topology_md5"]
    ), f"topology at {topology} is not the one the reference was frozen from"
    window = REFERENCE["frame_window"]
    universe = mda.Universe(str(topology), str(trajectory))
    universe.transfer_to_memory(start=window["start"], stop=window["stop"])
    return universe


@pytest.mark.parametrize("mode", MODES)
def test_profile_matches_frozen_reference(mode: str) -> None:
    """The ported plugin reproduces the frozen per-residue profile.

    Both implementations make the same MDAnalysis calls in the same order on
    the same arrays, so the difference is exactly zero and the tolerance is
    zero. A non-zero difference here means the calculation changed, not that
    floating point drifted.
    """
    frozen = REFERENCE["modes"][mode]
    window = REFERENCE["frame_window"]
    n_frames = window["stop"] - window["start"]
    observables = RMSF().compute(
        _universe(),
        FrameSelection(start=0, stop=n_frames, step=1, n_frames_total=n_frames),
        RMSFSettings.model_validate(frozen["settings"]),
    )
    profile = next(observable for observable in observables if observable.name == "rmsf")
    assert profile.kind == "profile" and profile.unit == "A"
    assert [int(value) for value in profile.index] == frozen["residue_ids"]
    np.testing.assert_allclose(profile.values, frozen["rmsf_per_residue"], rtol=0.0, atol=0.0)
    assert float(np.mean(profile.values)) == frozen["mean_rmsf"]
