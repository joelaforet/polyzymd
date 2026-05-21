"""Tests for centralized trajectory window helpers."""

from __future__ import annotations

import pytest

from polyzymd.analyses.shared.window import (
    resolve_replicate_trajectory_window,
    resolve_trajectory_window,
)


class _FakeLoader:
    """Minimal loader stub for trajectory window tests."""

    def __init__(self, timestep_ps: float = 100.0) -> None:
        self.timestep_ps = timestep_ps
        self.first_frame_time_ps: float | None = None
        self.calls: list[tuple[int, str]] = []

    def get_timestep(self, replicate: int, unit: str = "ps") -> float:
        """Return a fixed timestep and record the query."""

        self.calls.append((replicate, unit))
        if unit != "ps":
            raise ValueError(f"Unexpected unit: {unit}")
        return self.timestep_ps

    def get_first_frame_time(self, replicate: int, unit: str = "ps") -> float | None:
        """Return a fixed first-frame timestamp and record the query."""

        self.calls.append((replicate, f"first:{unit}"))
        if unit != "ps":
            raise ValueError(f"Unexpected unit: {unit}")
        return self.first_frame_time_ps


def test_resolve_trajectory_window_applies_equilibration_defaults() -> None:
    """Equilibration should define the default start frame."""

    window = resolve_trajectory_window(
        equilibration="1ns",
        n_frames_total=20,
        timestep_ps=100.0,
    )

    assert window.start == 10
    assert window.stop == 20
    assert window.step == 1
    assert window.equilibration_start == 10
    assert window.n_frames_selected == 10
    assert window.run_kwargs() == {"start": 10, "stop": 20, "step": 1}


def test_resolve_trajectory_window_accepts_custom_stop_and_step() -> None:
    """Custom stop and step should produce a concrete frame count."""

    window = resolve_trajectory_window(
        equilibration="0ns",
        n_frames_total=12,
        timestep_ps=50.0,
        start=2,
        stop=11,
        step=3,
    )

    assert window.start == 2
    assert window.stop == 11
    assert window.step == 3
    assert window.n_frames_selected == 3


def test_resolve_trajectory_window_rounds_equilibration_up_to_next_frame() -> None:
    """Non-integral equilibration should start at the first later frame."""

    window = resolve_trajectory_window(
        equilibration="850ps",
        n_frames_total=10,
        timestep_ps=100.0,
    )

    assert window.start == 9
    assert window.stop == 10
    assert window.equilibration_start == 9
    assert window.n_frames_selected == 1


def test_resolve_trajectory_window_rejects_start_before_equilibration() -> None:
    """Explicit starts may not move the window before equilibration."""

    with pytest.raises(ValueError, match="precedes equilibration-resolved start frame"):
        resolve_trajectory_window(
            equilibration="500ps",
            n_frames_total=10,
            timestep_ps=100.0,
            start=2,
        )


def test_resolve_trajectory_window_rejects_too_few_selected_frames() -> None:
    """Window validation should enforce the minimum selected frame count."""

    with pytest.raises(ValueError, match="need at least 4"):
        resolve_trajectory_window(
            equilibration="0ns",
            n_frames_total=10,
            timestep_ps=100.0,
            step=4,
            min_frames=4,
        )


def test_resolve_replicate_trajectory_window_uses_loader_timestep() -> None:
    """Loader-backed helper should fetch the timestep once in picoseconds."""

    loader = _FakeLoader(timestep_ps=200.0)

    window = resolve_replicate_trajectory_window(
        loader=loader,
        replicate=3,
        equilibration="400ps",
        n_frames_total=8,
        step=2,
    )

    assert loader.calls == [(3, "ps"), (3, "first:ps")]
    assert window.start == 2
    assert window.stop == 8
    assert window.step == 2
    assert window.n_frames_selected == 3


def test_resolve_trajectory_window_surfaces_large_equilibration_warning() -> None:
    """Large but valid equilibration windows should preserve the warning text."""

    window = resolve_trajectory_window(
        equilibration="6ns",
        n_frames_total=10,
        timestep_ps=1000.0,
    )

    assert window.warning_message is not None
    assert "Skipping 60.0% of trajectory" in window.warning_message


def test_resolve_trajectory_window_rejects_nonintegral_equilibration_near_end() -> None:
    """Equilibration should fail when no later discrete frame exists."""

    with pytest.raises(ValueError, match="leaves no frame at or after equilibration"):
        resolve_trajectory_window(
            equilibration="9.1ns",
            n_frames_total=10,
            timestep_ps=1000.0,
        )


def test_resolve_trajectory_window_uses_absolute_first_frame_timestamp() -> None:
    """Finite timestamps should make equilibration absolute trajectory time."""

    window = resolve_trajectory_window(
        equilibration="200ns",
        n_frames_total=14,
        timestep_ps=400.0,
        first_frame_time_ps=198_400.0,
    )

    assert window.start == 4
    assert window.equilibration_start == 4
    assert window.first_frame_time_ps == 198_400.0
    assert window.selected_start_time_ps == pytest.approx(200_000.0)
    assert window.equilibration_time_reference == "trajectory_timestamp"


@pytest.mark.parametrize("first_frame_time_ps", [None, float("nan"), float("inf"), "invalid"])
def test_resolve_trajectory_window_falls_back_without_valid_timestamp(
    first_frame_time_ps: object,
) -> None:
    """Missing or non-finite timestamps should preserve legacy behavior."""

    window = resolve_trajectory_window(
        equilibration="200ns",
        n_frames_total=600,
        timestep_ps=400.0,
        first_frame_time_ps=first_frame_time_ps,
    )

    assert window.start == 500
    assert window.selected_start_time_ps == pytest.approx(200_000.0)
    assert window.equilibration_time_reference == "loaded_frame_zero"


def test_resolve_trajectory_window_clamps_when_equilibration_precedes_first_frame() -> None:
    """Absolute equilibration before the loaded segment should start at frame zero."""

    window = resolve_trajectory_window(
        equilibration="100ns",
        n_frames_total=14,
        timestep_ps=400.0,
        first_frame_time_ps=198_400.0,
    )

    assert window.start == 0
    assert window.selected_start_time_ps == pytest.approx(198_400.0)
    assert window.equilibration_time_reference == "trajectory_timestamp"


def test_resolve_trajectory_window_rejects_absolute_equilibration_after_end() -> None:
    """Absolute equilibration beyond the last loaded timestamp should fail clearly."""

    with pytest.raises(ValueError, match="loaded trajectory timestamps"):
        resolve_trajectory_window(
            equilibration="205ns",
            n_frames_total=14,
            timestep_ps=400.0,
            first_frame_time_ps=198_400.0,
        )


def test_resolve_replicate_trajectory_window_accepts_timestep_override() -> None:
    """Loader-backed helper should preserve caller timestep override behavior."""

    loader = _FakeLoader(timestep_ps=999.0)
    loader.first_frame_time_ps = 198_400.0

    window = resolve_replicate_trajectory_window(
        loader=loader,
        replicate=3,
        equilibration="200ns",
        n_frames_total=14,
        timestep_ps=400.0,
    )

    assert loader.calls == [(3, "first:ps")]
    assert window.timestep_ps == 400.0
    assert window.start == 4
