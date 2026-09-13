"""Parity of the ported rmsd plugin against the pre-port implementation.

The reference file holds the per-frame RMSD the old package produced on one real
replicate of the 363 K control condition, for four reference modes. It was
frozen from the old source exported with ``git archive`` at the commit it names,
so the numbers cannot drift with the branch.

Two tolerances appear below. ``average`` mode still superimposes the trajectory
in place, because its reference is a mean structure, and it reproduces the old
series exactly. The other modes no longer run that pass: MDAnalysis minimises
over rigid motions inside the measurement, so aligning first changed no value
except through the float32 rounding of the rewritten coordinates, which is worth
about 3e-7 A on a 1.3 A deviation. That rounding is what the loose tolerance
covers.

The test skips when the trajectory is not on this machine, so a checkout without
the data still passes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.rmsd import RMSD, RMSDRunSettings, RMSDSettings, observable_name

REFERENCE_PATH = Path(__file__).with_name("rmsd_reference.json")
REFERENCE = json.loads(REFERENCE_PATH.read_text())

#: Reproduced bit for bit, because the reference construction is unchanged.
ATOL_EXACT = 1e-12
#: float32 rounding the old in-place alignment pass introduced, now dropped.
ATOL_NO_ALIGNMENT_PASS = 1e-6
TOLERANCES = {"average": ATOL_EXACT}


def _skip_unless_available(paths: list[str]) -> None:
    """Skip the test when a recorded input file is not on this machine."""
    missing = [path for path in paths if not Path(path).exists()]
    if missing:
        pytest.skip(f"real-data parity inputs are not on this machine: {missing}")


@pytest.fixture(scope="module")
def universe_factory():
    """Return a factory that builds a fresh universe over the frozen inputs."""
    pytest.importorskip("MDAnalysis")
    _skip_unless_available([REFERENCE["topology"], REFERENCE["trajectory"]])
    import MDAnalysis as mda

    def build():
        return mda.Universe(REFERENCE["topology"], REFERENCE["trajectory"])

    return build


def test_the_frozen_inputs_are_the_files_that_were_measured() -> None:
    """A different topology would make every number below meaningless."""
    _skip_unless_available([REFERENCE["topology"]])
    digest = hashlib.md5(Path(REFERENCE["topology"]).read_bytes()).hexdigest()

    assert digest == REFERENCE["topology_md5"]


@pytest.mark.parametrize("case_name", sorted(REFERENCE["cases"]))
def test_per_frame_rmsd_matches_the_old_plugin(universe_factory, case_name: str) -> None:
    """Every per-frame value of every reference mode survives the port."""
    case = REFERENCE["cases"][case_name]
    stored = dict(case["settings"])
    if stored.get("reference_file"):
        _skip_unless_available([stored["reference_file"]])
    atol = TOLERANCES.get(case_name, ATOL_NO_ALIGNMENT_PASS)

    with pytest.warns(DeprecationWarning, match="convergence"):
        run = RMSDRunSettings.model_validate(stored)
    settings = RMSDSettings(runs=[run])

    universe = universe_factory()
    window = REFERENCE["frame_window"]
    frames = FrameSelection(
        start=window["start"],
        stop=window["stop"],
        step=window["step"],
        n_frames_total=len(universe.trajectory),
        timestep_ps=40.0,
    )
    observables = RMSD().compute(universe, frames, settings)

    assert [observable.name for observable in observables] == [observable_name(run)]
    observable = observables[0]
    assert observable.kind == "mean_of_timeseries"
    assert observable.unit == "A"
    values = np.asarray(observable.values, dtype=np.float64)
    expected = np.asarray(case["rmsd"], dtype=np.float64)
    np.testing.assert_allclose(values, expected, rtol=0, atol=atol)
    assert float(np.mean(values)) == pytest.approx(case["mean"], abs=atol)
    assert float(np.std(values)) == pytest.approx(case["std"], abs=atol)
    assert float(np.median(values)) == pytest.approx(case["median"], abs=atol)


def test_the_observable_name_states_the_reference_mode() -> None:
    """Two runs on one selection stay distinguishable by reference mode."""
    centroid = RMSDRunSettings(label="Whole Protein CA", reference_mode="centroid")
    frame = RMSDRunSettings(label="Whole Protein CA", reference_mode="frame")
    assert observable_name(centroid) == "rmsd_whole_protein_ca_ref_centroid"
    assert observable_name(frame) == "rmsd_whole_protein_ca_ref_frame"
