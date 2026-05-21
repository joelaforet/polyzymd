"""Trajectory window helpers for trajectory-backed analyses.

This module centralizes the frame-window logic shared by analysis plugins that
need to combine equilibration skipping with MDAnalysis ``run()`` slice
arguments. The helpers return a validated window that can be passed directly to
trajectory-native runners without PolyzyMD re-owning the frame loop.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

from polyzymd.analyses.shared.diagnostics import validate_equilibration_time
from polyzymd.analyses.shared.loader import convert_time, parse_time_string

if TYPE_CHECKING:
    from polyzymd.analyses.shared.loader import TrajectoryLoader


def _equilibration_start_frame(equilibration_ps: float, timestep_ps: float) -> int:
    """Return the first frame at or after equilibration.

    Parameters
    ----------
    equilibration_ps : float
        Equilibration time in picoseconds.
    timestep_ps : float
        Time between consecutive frames in picoseconds.

    Returns
    -------
    int
        First frame index whose timestamp is greater than or equal to the
        equilibration time.
    """

    frame_position = equilibration_ps / timestep_ps
    rounded_position = round(frame_position)
    if math.isclose(frame_position, rounded_position, rel_tol=1e-12, abs_tol=1e-12):
        return int(rounded_position)
    return int(math.floor(frame_position)) + 1


@dataclass(frozen=True)
class TrajectoryWindow:
    """Validated frame window for a trajectory-backed analysis.

    Parameters
    ----------
    start : int
        Inclusive start frame for the analysis run.
    stop : int
        Exclusive stop frame for the analysis run.
    step : int
        Frame stride passed to ``runner.run(step=...)``.
    equilibration_start : int
        Start frame implied by the equilibration time alone.
    n_frames_total : int
        Total number of frames in the trajectory.
    n_frames_selected : int
        Number of frames selected by ``start``, ``stop``, and ``step``.
    timestep_ps : float
        Trajectory timestep in picoseconds.
    equilibration_ps : float
        Equilibration time converted to picoseconds.
    equilibration : str | None, optional
        Original equilibration time string used to resolve the window.
    first_frame_time_ps : float | None, optional
        Absolute MDAnalysis timestamp of the first loaded frame in
        picoseconds, when available.
    selected_start_time_ps : float | None, optional
        Timestamp of the selected start frame in the active time reference.
    equilibration_time_reference : str, optional
        Time reference used to interpret ``equilibration``. ``"trajectory_timestamp"``
        means absolute MDAnalysis timestamps were available; ``"loaded_frame_zero"``
        means the legacy loaded-frame-relative origin was used.
    warning_message : str | None
        Non-fatal equilibration warning generated during validation.
    """

    start: int
    stop: int
    step: int
    equilibration_start: int
    n_frames_total: int
    n_frames_selected: int
    timestep_ps: float
    equilibration_ps: float
    equilibration: str | None = None
    first_frame_time_ps: float | None = None
    selected_start_time_ps: float | None = None
    equilibration_time_reference: str = "loaded_frame_zero"
    warning_message: str | None = None

    def run_kwargs(self) -> dict[str, int]:
        """Return keyword arguments for ``MDAnalysis`` runner ``run()``.

        Returns
        -------
        dict[str, int]
            ``start``, ``stop``, and ``step`` values for ``run()``.
        """

        return {"start": self.start, "stop": self.stop, "step": self.step}


def resolve_replicate_trajectory_window(
    *,
    loader: "TrajectoryLoader",
    replicate: int,
    equilibration: str,
    n_frames_total: int,
    start: int | None = None,
    stop: int | None = None,
    step: int = 1,
    min_frames: int = 1,
    timestep_ps: float | None = None,
) -> TrajectoryWindow:
    """Resolve a validated window using loader trajectory timing metadata.

    Parameters
    ----------
    loader : TrajectoryLoader
        Loader for the replicate being analyzed.
    replicate : int
        Replicate number.
    equilibration : str
        Equilibration time string such as ``"10ns"``.
    n_frames_total : int
        Total number of frames in the trajectory.
    start : int | None, optional
        Absolute start frame for the analysis window. When ``None``, the
        equilibration-resolved start frame is used.
    stop : int | None, optional
        Absolute exclusive stop frame. When ``None``, the full remaining
        trajectory is used.
    step : int, optional
        Frame stride, by default 1.
    min_frames : int, optional
        Minimum required number of selected frames, by default 1.
    timestep_ps : float | None, optional
        Explicit timestep override in picoseconds. When ``None``, the loader
        timestep is used.

    Returns
    -------
    TrajectoryWindow
        Validated trajectory window with materialized ``run()`` arguments.
    """

    resolved_timestep_ps = (
        float(timestep_ps) if timestep_ps is not None else float(loader.get_timestep(replicate, unit="ps"))
    )
    try:
        first_frame_time_ps = loader.get_first_frame_time(replicate, unit="ps")
    except (AttributeError, TypeError, ValueError):
        first_frame_time_ps = None
    return resolve_trajectory_window(
        equilibration=equilibration,
        n_frames_total=n_frames_total,
        timestep_ps=resolved_timestep_ps,
        first_frame_time_ps=first_frame_time_ps,
        start=start,
        stop=stop,
        step=step,
        min_frames=min_frames,
    )


def resolve_trajectory_window(
    *,
    equilibration: str,
    n_frames_total: int,
    timestep_ps: float,
    start: int | None = None,
    stop: int | None = None,
    step: int = 1,
    min_frames: int = 1,
    first_frame_time_ps: float | None = None,
) -> TrajectoryWindow:
    """Resolve and validate a trajectory frame window.

    When the first loaded frame has a finite MDAnalysis timestamp,
    ``equilibration`` is interpreted as an absolute trajectory time. The start
    frame is the first loaded frame whose timestamp is greater than or equal to
    the equilibration time. When timestamp metadata is unavailable, the legacy
    loaded-frame-relative origin is used.

    Parameters
    ----------
    equilibration : str
        Equilibration time string such as ``"10ns"``.
    n_frames_total : int
        Total number of frames in the trajectory.
    timestep_ps : float
        Time between consecutive frames in picoseconds.
    start : int | None, optional
        Absolute start frame for the analysis window. When ``None``, the
        equilibration-resolved start frame is used.
    stop : int | None, optional
        Absolute exclusive stop frame. When ``None``, the trajectory end is
        used.
    step : int, optional
        Frame stride, by default 1.
    min_frames : int, optional
        Minimum required number of selected frames, by default 1.
    first_frame_time_ps : float | None, optional
        Absolute MDAnalysis timestamp of loaded frame 0 in picoseconds. Non-finite
        values are ignored and use the legacy loaded-frame-relative behavior.

    Returns
    -------
    TrajectoryWindow
        Validated frame window.

    Raises
    ------
    ValueError
        Raised when the timestep, equilibration, or window arguments are
        inconsistent with the trajectory.
    """

    if n_frames_total < 1:
        raise ValueError("Trajectory must contain at least one frame")
    if timestep_ps <= 0:
        raise ValueError(f"Trajectory timestep must be positive, got {timestep_ps}")
    if step < 1:
        raise ValueError(f"step must be >= 1, got {step}")
    if min_frames < 1:
        raise ValueError(f"min_frames must be >= 1, got {min_frames}")

    eq_value, eq_unit = parse_time_string(equilibration)
    equilibration_ps = convert_time(eq_value, eq_unit, "ps")
    finite_first_frame_time_ps = (
        float(first_frame_time_ps)
        if first_frame_time_ps is not None and math.isfinite(float(first_frame_time_ps))
        else None
    )
    if finite_first_frame_time_ps is None:
        equilibration_offset_ps = equilibration_ps
        time_reference = "loaded_frame_zero"
    else:
        equilibration_offset_ps = max(0.0, equilibration_ps - finite_first_frame_time_ps)
        time_reference = "trajectory_timestamp"
    equilibration_ns = equilibration_offset_ps / 1000.0
    trajectory_ns = (n_frames_total * timestep_ps) / 1000.0

    is_valid, warning_message = validate_equilibration_time(equilibration_ns, trajectory_ns)
    if not is_valid:
        raise ValueError(warning_message or "Invalid equilibration window")

    equilibration_start = _equilibration_start_frame(equilibration_offset_ps, timestep_ps)
    if equilibration_start >= n_frames_total:
        if finite_first_frame_time_ps is None:
            raise ValueError(
                "Equilibration time "
                f"({equilibration_ps / 1000.0:.3f} ns) leaves no frame at or after "
                f"equilibration in a {n_frames_total}-frame trajectory with timestep "
                f"{timestep_ps:.3f} ps"
            )
        last_frame_time_ps = finite_first_frame_time_ps + (n_frames_total - 1) * timestep_ps
        raise ValueError(
            "Equilibration time "
            f"({equilibration_ps:.3f} ps) leaves no frame at or after equilibration "
            f"in loaded trajectory timestamps {finite_first_frame_time_ps:.3f} ps to "
            f"{last_frame_time_ps:.3f} ps with timestep {timestep_ps:.3f} ps"
        )

    resolved_start = equilibration_start if start is None else start
    resolved_stop = n_frames_total if stop is None else stop

    if resolved_start < 0:
        raise ValueError(f"start must be >= 0, got {resolved_start}")
    if resolved_start < equilibration_start:
        raise ValueError(
            f"start={resolved_start} precedes equilibration-resolved start frame "
            f"{equilibration_start}"
        )
    if resolved_start >= n_frames_total:
        raise ValueError(
            f"start={resolved_start} is outside the trajectory range [0, {n_frames_total})"
        )
    if resolved_stop <= resolved_start:
        raise ValueError(f"stop={resolved_stop} must be greater than start={resolved_start}")
    if resolved_stop > n_frames_total:
        raise ValueError(f"stop={resolved_stop} exceeds trajectory length {n_frames_total}")

    n_frames_selected = len(range(resolved_start, resolved_stop, step))
    if n_frames_selected < min_frames:
        raise ValueError(
            "Trajectory window "
            f"[{resolved_start}:{resolved_stop}:{step}] selects {n_frames_selected} frame(s), "
            f"need at least {min_frames}"
        )

    selected_start_time_ps = (
        finite_first_frame_time_ps + resolved_start * timestep_ps
        if finite_first_frame_time_ps is not None
        else resolved_start * timestep_ps
    )

    return TrajectoryWindow(
        start=resolved_start,
        stop=resolved_stop,
        step=step,
        equilibration_start=equilibration_start,
        n_frames_total=n_frames_total,
        n_frames_selected=n_frames_selected,
        timestep_ps=float(timestep_ps),
        equilibration_ps=float(equilibration_ps),
        equilibration=equilibration,
        first_frame_time_ps=finite_first_frame_time_ps,
        selected_start_time_ps=float(selected_start_time_ps),
        equilibration_time_reference=time_reference,
        warning_message=warning_message,
    )
