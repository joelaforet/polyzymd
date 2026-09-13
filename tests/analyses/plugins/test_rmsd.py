"""Known-answer and settings tests for the contract rmsd plugin.

The test universe is four unit-mass atoms on a cross whose scale grows frame by
frame. Scaling a centred, symmetric shape by ``s`` leaves no rotation to find,
so the minimised RMSD of frame ``i`` against the unit frame is exactly
``s_i - 1`` and the mean over the window is known by hand.
"""

from __future__ import annotations

from typing import Any

import pytest

from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.exceptions import ReplicateError, SelectionError
from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.rmsd import (
    RMSD,
    RMSDAnalysis,
    RMSDRunSettings,
    RMSDSettings,
    observable_name,
)
from tests.analyses.conftest import CROSS

SCALES = (1.0, 1.2, 1.4, 1.6, 1.8)
EXPECTED_MEAN = 0.4  # mean of |s - 1| over SCALES
TAIL_OFFSETS = (0.0, 0.1, 0.2, 0.3, 0.4)


def _growing_cross() -> Any:
    """Universe whose frames are the unit cross scaled by ``SCALES``."""
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("masses", [1.0] * 4)
    cross = np.asarray(CROSS, dtype=np.float32)
    universe.load_new(np.stack([cross * scale for scale in SCALES]), format=MemoryReader)
    return universe


def _core_and_tail() -> Any:
    """Rigid four-atom core plus a rigid three-atom tail that slides along z.

    Superposing on the core leaves the transform at the identity, because the
    core is the same in every frame and in the reference, so the deviation of
    the tail is exactly how far it slid. Superposing on the tail instead carries
    the tail onto itself and reports nothing.
    """
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    names = ["CA"] * 4 + ["TL"] * 3
    universe = mda.Universe.empty(
        len(names), n_residues=1, atom_resindex=[0] * len(names), trajectory=True
    )
    universe.add_TopologyAttr("masses", [1.0] * len(names))
    universe.add_TopologyAttr("names", names)
    cross = np.asarray(CROSS, dtype=np.float32)
    tail = np.asarray([[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]], dtype=np.float32)
    frames = [
        np.vstack([cross, tail + np.asarray([0.0, 0.0, offset], dtype=np.float32)])
        for offset in TAIL_OFFSETS
    ]
    universe.load_new(np.stack(frames), format=MemoryReader)
    return universe


def _run(**overrides: Any) -> RMSDRunSettings:
    """Run settings measuring every atom against the first frame."""
    return RMSDRunSettings.model_validate(
        {
            "label": "core",
            "selection": "all",
            "alignment_selection": "all",
            "reference_mode": "frame",
            "reference_frame": 0,
            **overrides,
        }
    )


def _frames(universe: Any) -> FrameSelection:
    """Whole-trajectory frame selection for a direct ``compute`` call."""
    return FrameSelection(start=0, stop=len(universe.trajectory), step=1)


def test_per_frame_rmsd_is_the_change_in_scale() -> None:
    """Each frame deviates from the reference by exactly the scale it grew by."""
    universe = _growing_cross()

    observables = RMSD().compute(universe, _frames(universe), RMSDSettings(runs=[_run()]))

    assert len(observables) == 1
    observable = observables[0]
    assert observable.name == "rmsd_core_ref_frame"
    assert observable.unit == "A"
    assert observable.kind == "mean_of_timeseries"
    assert observable.higher_is_better is False
    assert observable.values == pytest.approx([scale - 1.0 for scale in SCALES], abs=1e-6)


def test_the_average_reference_sits_at_the_mean_scale() -> None:
    """Against the average structure the deviation is measured from mean scale."""
    universe = _growing_cross()
    settings = RMSDSettings(runs=[_run(reference_mode="average")])

    observable = RMSD().compute(universe, _frames(universe), settings)[0]

    mean_scale = sum(SCALES) / len(SCALES)
    assert observable.name == "rmsd_core_ref_average"
    assert observable.values == pytest.approx([abs(s - mean_scale) for s in SCALES], abs=1e-6)


def test_the_condition_aggregate_is_a_mean_over_replicates(run_contract_analysis) -> None:
    """Three identical replicates give the known mean with no spread."""
    artifact = run_contract_analysis(
        RMSDAnalysis, RMSDSettings(runs=[_run()]), lambda replicate: _growing_cross()
    )

    aggregate = ObservableAggregate.model_validate(artifact.payload["observables"][0])

    assert aggregate.name == "rmsd_core_ref_frame"
    assert aggregate.unit == "A"
    assert aggregate.n_replicates == 3
    assert aggregate.replicate_values == pytest.approx([EXPECTED_MEAN] * 3, abs=1e-6)
    assert aggregate.mean == pytest.approx(EXPECTED_MEAN, abs=1e-6)
    assert aggregate.sem == pytest.approx(0.0, abs=1e-9)


def test_two_reference_modes_report_two_observables() -> None:
    """The reference mode is part of the name, so a table cannot hide it."""
    runs = [_run(label="core frame"), _run(label="core average", reference_mode="average")]
    universe = _growing_cross()

    observables = RMSD().compute(universe, _frames(universe), RMSDSettings(runs=runs))

    assert [observable.name for observable in observables] == [
        "rmsd_core_frame_ref_frame",
        "rmsd_core_average_ref_average",
    ]
    assert observable_name(runs[0]) == "rmsd_core_frame_ref_frame"


def test_an_empty_selection_raises_a_typed_error() -> None:
    """An empty selection is an error, not an RMSD of zero."""
    universe = _growing_cross()
    settings = RMSDSettings(runs=[_run(selection="index 99")])

    with pytest.raises(SelectionError, match="matched no atoms"):
        RMSD().compute(universe, _frames(universe), settings)


def test_the_dropped_convergence_settings_still_parse() -> None:
    """A comparison file from the old plugin loads, with a deprecation warning."""
    with pytest.warns(UserWarning, match="convergence_window_size_ns"):
        run = _run(convergence_window_size_ns=15.0, convergence_slope_threshold=0.0005)

    assert not hasattr(run, "convergence_window_size_ns")


def test_duplicate_run_labels_are_rejected() -> None:
    """Two runs cannot write one observable name."""
    with pytest.raises(ValueError, match="unique"):
        RMSDSettings(runs=[_run(), _run()])


def test_a_missing_external_reference_warns_then_fails_at_compute_time() -> None:
    """A comparison file naming a cluster path validates off the cluster."""
    with pytest.warns(UserWarning, match="not on this machine"):
        run = _run(reference_mode="external", reference_file="/nowhere/missing.pdb")
    universe = _growing_cross()

    with pytest.raises(ReplicateError, match="does not exist"):
        RMSD().compute(universe, _frames(universe), RMSDSettings(runs=[run]))


def test_external_mode_still_needs_a_reference_file() -> None:
    """The setting itself is required; only its presence on disk is deferred."""
    with pytest.raises(ValueError, match="needs reference_file"):
        _run(reference_mode="external")


def test_the_alignment_selection_chooses_what_superposition_minimises() -> None:
    """Superposing on the core measures the tail; superposing on the tail hides it."""
    universe = _core_and_tail()
    on_tail = _run(label="tail", selection="name TL", alignment_selection="name TL")
    on_core = _run(label="tail", selection="name TL", alignment_selection="name CA")

    settings = RMSDSettings(runs=[on_tail])
    hidden = RMSD().compute(universe, _frames(universe), settings)[0]
    settings = RMSDSettings(runs=[on_core])
    measured = RMSD().compute(universe, _frames(universe), settings)[0]

    assert hidden.values == pytest.approx([0.0] * len(TAIL_OFFSETS), abs=1e-6)
    assert measured.values == pytest.approx(list(TAIL_OFFSETS), abs=1e-6)


def test_two_labels_that_slug_the_same_are_rejected() -> None:
    """Two runs whose names collide would overwrite each other in the sidecar."""
    with pytest.raises(ValueError, match="unique observable names"):
        RMSDSettings(runs=[_run(label="Core Frame"), _run(label="core_frame")])
