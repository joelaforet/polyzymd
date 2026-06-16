"""Tests for MDAnalysis frame selection mapping."""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses.mda import FrameSelection
from polyzymd.analyses.shared.window import resolve_trajectory_window


def test_slice_selection_returns_non_none_run_kwargs() -> None:
    """Start, stop, and step should map directly to MDAnalysis kwargs."""

    selection = FrameSelection(start=2, stop=11, step=3, n_frames_total=12)

    assert selection.run_kwargs() == {"start": 2, "stop": 11, "step": 3}
    assert selection.n_frames_selected == 3


def test_slice_selection_omits_none_run_kwargs() -> None:
    """Unset slice arguments should not be forwarded to MDAnalysis."""

    selection = FrameSelection(stop=5)

    assert selection.run_kwargs() == {"stop": 5}
    assert selection.n_frames_selected == 5


def test_numpy_slice_scalars_return_python_run_kwargs() -> None:
    """NumPy slice scalars should normalize before MDAnalysis receives them."""

    selection = FrameSelection(
        start=np.int64(1),
        stop=np.int64(6),
        step=np.int64(2),
        equilibration_start=np.int64(1),
        equilibration_ps=np.float32(100.0),
        timestep_ps=np.float64(2.0),
        first_frame_time_ps=np.float64(10.0),
        selected_start_time_ps=np.float32(12.0),
        n_frames_total=np.int64(8),
    )

    run_kwargs = selection.run_kwargs()

    assert run_kwargs == {"start": 1, "stop": 6, "step": 2}
    assert all(
        isinstance(value, int) and not hasattr(value, "item") for value in run_kwargs.values()
    )
    assert isinstance(selection.equilibration_start, int)
    assert isinstance(selection.equilibration_ps, float)
    assert isinstance(selection.timestep_ps, float)
    assert isinstance(selection.first_frame_time_ps, float)
    assert isinstance(selection.selected_start_time_ps, float)
    assert isinstance(selection.n_frames_total, int)
    assert isinstance(selection.n_frames_selected, int)
    assert not hasattr(selection.equilibration_start, "item")
    assert not hasattr(selection.equilibration_ps, "item")
    assert not hasattr(selection.timestep_ps, "item")
    assert not hasattr(selection.first_frame_time_ps, "item")
    assert not hasattr(selection.selected_start_time_ps, "item")
    assert not hasattr(selection.n_frames_total, "item")
    assert not hasattr(selection.n_frames_selected, "item")


def test_frames_selection_returns_only_frames() -> None:
    """Explicit frames should be forwarded without slice arguments."""

    selection = FrameSelection(frames=[0, 3, 8], n_frames_total=10)

    assert selection.run_kwargs() == {"frames": [0, 3, 8]}
    assert selection.frames == (0, 3, 8)
    assert selection.n_frames_selected == 3


@pytest.mark.parametrize("field", ["start", "stop", "step"])
def test_frames_rejects_combined_slice_arguments(field: str) -> None:
    """MDAnalysis does not allow frames with start, stop, or step."""

    kwargs = {field: 1}

    with pytest.raises(ValueError, match="frames cannot be combined"):
        FrameSelection(frames=[0, 1], **kwargs)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"step": 0}, "step must be >= 1"),
        ({"start": -1}, "start must be >= 0"),
        ({"start": 5, "stop": 5}, "must be greater than start"),
        ({"start": 10, "n_frames_total": 10}, "outside the trajectory range"),
        ({"stop": 11, "n_frames_total": 10}, "exceeds trajectory length"),
    ],
)
def test_slice_selection_rejects_invalid_arguments(kwargs: dict[str, int], message: str) -> None:
    """Invalid slice arguments should fail before MDAnalysis is called."""

    with pytest.raises(ValueError, match=message):
        FrameSelection(**kwargs)


@pytest.mark.parametrize("frames", [[], "1:5", 3, 2.5])
def test_frames_rejects_invalid_selectors(frames: object) -> None:
    """Scalar, string-like, and empty frame selectors are ambiguous."""

    with pytest.raises(ValueError, match="frames"):
        FrameSelection(frames=frames)


@pytest.mark.parametrize(
    "frames",
    [
        [1.5],
        [None],
        ["0"],
        [0, 1.5],
        [True, 1],
        np.array([1.5]),
        np.array([None], dtype=object),
        np.array(["0"]),
    ],
)
def test_frames_rejects_non_integer_non_boolean_elements(frames: object) -> None:
    """Explicit frames should contain only integer indices or boolean masks."""

    with pytest.raises(ValueError, match="integer frame indices or boolean mask"):
        FrameSelection(frames=frames)


def test_numpy_integer_frame_indices_are_valid() -> None:
    """Array-like integer frame selectors should remain valid."""

    selection = FrameSelection(frames=np.array([0, 3, 8]), n_frames_total=10)

    run_frames = selection.run_kwargs()["frames"]

    assert run_frames == [0, 3, 8]
    assert all(isinstance(frame, int) and not isinstance(frame, bool) for frame in run_frames)
    assert selection.n_frames_selected == 3


def test_boolean_frame_mask_counts_selected_frames() -> None:
    """Boolean masks should record the number of selected frames."""

    selection = FrameSelection(frames=[True, False, True, False], n_frames_total=4)

    assert selection.run_kwargs() == {"frames": [True, False, True, False]}
    assert selection.n_frames_selected == 2


def test_numpy_boolean_frame_mask_counts_selected_frames() -> None:
    """Array-like boolean masks should be validated as masks."""

    selection = FrameSelection(frames=np.array([True, False, True, False]), n_frames_total=4)

    run_frames = selection.run_kwargs()["frames"]

    assert run_frames == [True, False, True, False]
    assert all(isinstance(frame, bool) for frame in run_frames)
    assert selection.n_frames_selected == 2


def test_boolean_frame_mask_rejects_length_mismatch() -> None:
    """Boolean masks should match known trajectory length."""

    with pytest.raises(ValueError, match="does not match trajectory length"):
        FrameSelection(frames=[True, False], n_frames_total=3)


def test_numpy_boolean_frame_mask_rejects_length_mismatch() -> None:
    """Array-like boolean masks should match known trajectory length."""

    with pytest.raises(ValueError, match="does not match trajectory length"):
        FrameSelection(frames=np.array([True, False]), n_frames_total=3)


def test_boolean_frame_mask_rejects_empty_selection() -> None:
    """Boolean masks should select at least one frame."""

    with pytest.raises(ValueError, match="select at least one frame"):
        FrameSelection(frames=[False, False], n_frames_total=2)


def test_numpy_boolean_frame_mask_rejects_empty_selection() -> None:
    """Array-like boolean masks should select at least one frame."""

    with pytest.raises(ValueError, match="select at least one frame"):
        FrameSelection(frames=np.array([False, False]), n_frames_total=2)


def test_from_trajectory_window_preserves_window_values() -> None:
    """Existing shared windows should bridge to frame selections."""

    window = resolve_trajectory_window(
        equilibration="500ps",
        n_frames_total=10,
        timestep_ps=100.0,
        stop=9,
        step=2,
    )

    selection = FrameSelection.from_trajectory_window(window)

    assert selection.run_kwargs() == {"start": 5, "stop": 9, "step": 2}
    assert selection.equilibration_start == 5
    assert selection.equilibration_ps == 500.0
    assert selection.timestep_ps == 100.0
    assert selection.first_frame_time_ps is None
    assert selection.selected_start_time_ps == pytest.approx(500.0)
    assert selection.equilibration_time_reference == "loaded_frame_zero"
    assert selection.n_frames_total == 10
    assert selection.n_frames_selected == 2


def test_from_equilibration_uses_shared_window_semantics() -> None:
    """Equilibration conversion should match shared trajectory windows."""

    selection = FrameSelection.from_equilibration(
        equilibration="850ps",
        n_frames_total=10,
        timestep_ps=100.0,
    )

    assert selection.run_kwargs() == {"start": 9, "stop": 10, "step": 1}
    assert selection.equilibration == "850ps"
    assert selection.equilibration_start == 9
    assert selection.equilibration_ps == 850.0
    assert selection.n_frames_selected == 1


def test_from_equilibration_preserves_warning_message() -> None:
    """Large but valid equilibration windows should keep warning provenance."""

    selection = FrameSelection.from_equilibration(
        equilibration="6ns",
        n_frames_total=10,
        timestep_ps=1000.0,
    )

    assert selection.warning_message is not None
    assert "Skipping 60.0% of trajectory" in selection.warning_message


def test_from_equilibration_preserves_timestamp_fields() -> None:
    """FrameSelection should expose timestamp-aware window provenance."""

    selection = FrameSelection.from_equilibration(
        equilibration="200ns",
        n_frames_total=14,
        timestep_ps=400.0,
        first_frame_time_ps=198_400.0,
    )

    assert selection.run_kwargs() == {"start": 4, "stop": 14, "step": 1}
    assert selection.equilibration == "200ns"
    assert selection.equilibration_start == 4
    assert selection.equilibration_ps == 200_000.0
    assert selection.first_frame_time_ps == 198_400.0
    assert selection.selected_start_time_ps == pytest.approx(200_000.0)
    assert selection.equilibration_time_reference == "trajectory_timestamp"
