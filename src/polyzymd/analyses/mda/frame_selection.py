"""Frame-selection helpers for MDAnalysis extension-layer jobs."""

from __future__ import annotations

from collections.abc import Iterable, Sized
from dataclasses import dataclass, field
from operator import index
from typing import TYPE_CHECKING, Any

from polyzymd.analyses.mda.base import MDARunKwargs
from polyzymd.analyses.shared.window import resolve_trajectory_window

if TYPE_CHECKING:
    from polyzymd.analyses.shared.window import TrajectoryWindow


def _is_scalar_frame_selector(frames: Any) -> bool:
    """Return whether ``frames`` is a scalar selector value.

    Parameters
    ----------
    frames : Any
        Candidate explicit frame selector.

    Returns
    -------
    bool
        ``True`` when the selector is a scalar that MDAnalysis should not
        receive through the ``frames`` keyword.
    """

    return isinstance(frames, (bool, int, float, complex))


def _is_boolean_frame_value(frame: Any) -> bool:
    """Return whether a frame selector value is a boolean mask entry.

    Parameters
    ----------
    frame : Any
        Candidate element from an explicit ``frames`` selector.

    Returns
    -------
    bool
        ``True`` for Python booleans and NumPy-style boolean scalar values.
    """

    if isinstance(frame, bool):
        return True
    frame_type = type(frame)
    return (
        frame_type.__name__ == "bool_"
        and frame_type.__module__.split(".", maxsplit=1)[0] == "numpy"
    )


def _is_integer_frame_value(frame: Any) -> bool:
    """Return whether a frame selector value is an integer index.

    Parameters
    ----------
    frame : Any
        Candidate element from an explicit ``frames`` selector.

    Returns
    -------
    bool
        ``True`` when ``frame`` can be used as an integer index and is not a
        boolean mask value.
    """

    if _is_boolean_frame_value(frame):
        return False
    try:
        index(frame)
    except TypeError:
        return False
    return True


def _coerce_frames(frames: Any) -> tuple[Any, ...]:
    """Validate and freeze an explicit frame selector.

    Parameters
    ----------
    frames : Any
        Candidate frame sequence or boolean mask.

    Returns
    -------
    tuple[Any, ...]
        Immutable frame selector suitable for forwarding to MDAnalysis.

    Raises
    ------
    ValueError
        Raised when ``frames`` is empty, scalar, string-like, or unsized.
    """

    if _is_scalar_frame_selector(frames):
        raise ValueError("frames must be a non-empty sequence or boolean mask, not a scalar")
    if isinstance(frames, (str, bytes)):
        raise ValueError("frames must be a non-empty sequence or boolean mask, not a string")
    if not isinstance(frames, Sized) or not isinstance(frames, Iterable):
        raise ValueError("frames must be a sized iterable sequence or boolean mask")

    try:
        frozen_frames = tuple(frames)
    except TypeError as exc:
        raise ValueError("frames must be a sized iterable sequence or boolean mask") from exc
    if len(frozen_frames) == 0:
        raise ValueError("frames must contain at least one entry")
    return frozen_frames


def _selected_count_from_frames(frames: tuple[Any, ...], n_frames_total: int | None) -> int:
    """Return the selected frame count for an explicit frame selector.

    Parameters
    ----------
    frames : tuple[Any, ...]
        Frozen frame sequence or boolean mask.
    n_frames_total : int | None
        Total trajectory length when known.

    Returns
    -------
    int
        Number of frames selected by ``frames``.

    Raises
    ------
    ValueError
        Raised when a boolean mask has the wrong length or selects no frames.
    """

    is_bool_mask = all(_is_boolean_frame_value(frame) for frame in frames)
    is_integer_indices = all(_is_integer_frame_value(frame) for frame in frames)
    if not is_bool_mask and not is_integer_indices:
        raise ValueError("frames must contain only integer frame indices or boolean mask values")
    if is_integer_indices:
        return len(frames)

    if n_frames_total is not None and len(frames) != n_frames_total:
        raise ValueError(
            "Boolean frame mask length "
            f"{len(frames)} does not match trajectory length {n_frames_total}"
        )

    n_frames_selected = sum(1 for frame in frames if frame)
    if n_frames_selected == 0:
        raise ValueError("frames must select at least one frame")
    return n_frames_selected


def _selected_count_from_slice(
    *,
    start: int | None,
    stop: int | None,
    step: int | None,
    n_frames_total: int | None,
) -> int | None:
    """Return the selected frame count for start/stop/step selectors.

    Parameters
    ----------
    start : int | None
        Inclusive start frame.
    stop : int | None
        Exclusive stop frame.
    step : int | None
        Frame stride.
    n_frames_total : int | None
        Total trajectory length when known.

    Returns
    -------
    int | None
        Number of selected frames when it can be resolved, otherwise ``None``.
    """

    resolved_start = 0 if start is None else start
    resolved_stop = n_frames_total if stop is None else stop
    if resolved_stop is None:
        return None
    resolved_step = 1 if step is None else step
    return len(range(resolved_start, resolved_stop, resolved_step))


@dataclass(frozen=True)
class FrameSelection:
    """Validated MDAnalysis ``run()`` frame-selection arguments.

    ``FrameSelection`` is the import-light bridge between PolyzyMD's
    equilibration/window semantics and MDAnalysis ``AnalysisBase.run`` keyword
    arguments. Explicit ``frames`` selectors are mutually exclusive with
    ``start``, ``stop``, and ``step`` because MDAnalysis treats them as separate
    selection modes.

    Parameters
    ----------
    start : int | None, optional
        Inclusive start frame passed to ``run(start=...)``.
    stop : int | None, optional
        Exclusive stop frame passed to ``run(stop=...)``.
    step : int | None, optional
        Frame stride passed to ``run(step=...)``.
    frames : Any, optional
        Explicit frame indices or boolean mask passed to ``run(frames=...)``.
    equilibration : str | None, optional
        Original equilibration setting used to derive this selection.
    equilibration_start : int | None, optional
        Start frame implied by equilibration alone.
    equilibration_ps : float | None, optional
        Equilibration time converted to picoseconds.
    timestep_ps : float | None, optional
        Trajectory timestep in picoseconds.
    first_frame_time_ps : float | None, optional
        Absolute MDAnalysis timestamp of loaded frame 0 in picoseconds, when
        available.
    selected_start_time_ps : float | None, optional
        Timestamp of the selected start frame in the active time reference.
    equilibration_time_reference : str | None, optional
        Time reference used to interpret ``equilibration``.
    n_frames_total : int | None, optional
        Total trajectory frame count when known.
    warning_message : str | None, optional
        Non-fatal warning from equilibration/window validation.
    """

    start: int | None = None
    stop: int | None = None
    step: int | None = None
    frames: Any = None
    equilibration: str | None = None
    equilibration_start: int | None = None
    equilibration_ps: float | None = None
    timestep_ps: float | None = None
    first_frame_time_ps: float | None = None
    selected_start_time_ps: float | None = None
    equilibration_time_reference: str | None = None
    n_frames_total: int | None = None
    warning_message: str | None = None
    n_frames_selected: int | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        """Validate mutually exclusive MDAnalysis frame-selection modes.

        Raises
        ------
        ValueError
            Raised when the selector would produce invalid MDAnalysis ``run``
            keyword arguments or an empty known frame selection.
        """

        self._validate_total_frame_count()

        if self.frames is not None:
            self._validate_explicit_frames_mode()
            return

        self._validate_slice_mode()

    def run_kwargs(self) -> MDARunKwargs:
        """Return keyword arguments for ``MDAnalysis`` ``AnalysisBase.run``.

        Returns
        -------
        MDARunKwargs
            ``frames`` alone when explicit frames are set, otherwise the
            non-``None`` ``start``, ``stop``, and ``step`` values.
        """

        if self.frames is not None:
            return MDARunKwargs(frames=list(self.frames))

        kwargs = MDARunKwargs()
        if self.start is not None:
            kwargs["start"] = self.start
        if self.stop is not None:
            kwargs["stop"] = self.stop
        if self.step is not None:
            kwargs["step"] = self.step
        return kwargs

    @classmethod
    def from_trajectory_window(cls, window: TrajectoryWindow) -> FrameSelection:
        """Build a frame selection from a resolved PolyzyMD trajectory window.

        Parameters
        ----------
        window : TrajectoryWindow
            Existing validated shared trajectory window.

        Returns
        -------
        FrameSelection
            Selection that forwards the same ``start``, ``stop``, and ``step``
            values to MDAnalysis while preserving window provenance.
        """

        return cls(
            start=window.start,
            stop=window.stop,
            step=window.step,
            equilibration=getattr(window, "equilibration", None),
            equilibration_start=window.equilibration_start,
            equilibration_ps=window.equilibration_ps,
            timestep_ps=window.timestep_ps,
            first_frame_time_ps=window.first_frame_time_ps,
            selected_start_time_ps=window.selected_start_time_ps,
            equilibration_time_reference=window.equilibration_time_reference,
            n_frames_total=window.n_frames_total,
            warning_message=window.warning_message,
        )

    @classmethod
    def from_equilibration(
        cls,
        *,
        equilibration: str,
        n_frames_total: int,
        timestep_ps: float,
        start: int | None = None,
        stop: int | None = None,
        step: int = 1,
        min_frames: int = 1,
        first_frame_time_ps: float | None = None,
    ) -> FrameSelection:
        """Resolve PolyzyMD equilibration/window settings to a frame selection.

        Finite first-frame timestamps make equilibration absolute in MDAnalysis
        trajectory time; missing timestamps keep the stale loaded-frame-relative
        origin.

        Parameters
        ----------
        equilibration : str
            Equilibration time string such as ``"10ns"``.
        n_frames_total : int
            Total number of trajectory frames.
        timestep_ps : float
            Trajectory timestep in picoseconds.
        start : int | None, optional
            Absolute start frame. When ``None``, equilibration determines the
            start frame.
        stop : int | None, optional
            Absolute exclusive stop frame. When ``None``, the trajectory end is
            used.
        step : int, optional
            Frame stride, by default 1.
        min_frames : int, optional
            Minimum required number of selected frames, by default 1.
        first_frame_time_ps : float | None, optional
            Absolute MDAnalysis timestamp of loaded frame 0 in picoseconds.

        Returns
        -------
        FrameSelection
            Validated selection suitable for MDAnalysis ``run()``.
        """

        window = resolve_trajectory_window(
            equilibration=equilibration,
            n_frames_total=n_frames_total,
            timestep_ps=timestep_ps,
            start=start,
            stop=stop,
            step=step,
            min_frames=min_frames,
            first_frame_time_ps=first_frame_time_ps,
        )
        selection = cls.from_trajectory_window(window)
        object.__setattr__(selection, "equilibration", equilibration)
        return selection

    def _validate_total_frame_count(self) -> None:
        """Validate the optional total frame count.

        Raises
        ------
        ValueError
            Raised when the known trajectory length is invalid.
        """

        if self.n_frames_total is not None and self.n_frames_total < 1:
            raise ValueError("n_frames_total must be >= 1 when provided")

    def _validate_explicit_frames_mode(self) -> None:
        """Validate explicit ``frames`` selectors.

        Raises
        ------
        ValueError
            Raised when ``frames`` is mixed with slice arguments or selects no
            frames.
        """

        if self.start is not None or self.stop is not None or self.step is not None:
            raise ValueError("frames cannot be combined with start, stop, or step")

        frozen_frames = _coerce_frames(self.frames)
        n_frames_selected = _selected_count_from_frames(frozen_frames, self.n_frames_total)
        object.__setattr__(self, "frames", frozen_frames)
        object.__setattr__(self, "n_frames_selected", n_frames_selected)

    def _validate_slice_mode(self) -> None:
        """Validate ``start``/``stop``/``step`` selectors.

        Raises
        ------
        ValueError
            Raised when the slice arguments are invalid or select no known
            frames.
        """

        if self.step is not None and self.step < 1:
            raise ValueError(f"step must be >= 1, got {self.step}")
        if self.start is not None and self.start < 0:
            raise ValueError(f"start must be >= 0, got {self.start}")
        if self.stop is not None and self.start is not None and self.stop <= self.start:
            raise ValueError(f"stop={self.stop} must be greater than start={self.start}")
        if self.n_frames_total is not None:
            if self.start is not None and self.start >= self.n_frames_total:
                raise ValueError(
                    f"start={self.start} is outside the trajectory range [0, {self.n_frames_total})"
                )
            if self.stop is not None and self.stop > self.n_frames_total:
                raise ValueError(
                    f"stop={self.stop} exceeds trajectory length {self.n_frames_total}"
                )

        n_frames_selected = _selected_count_from_slice(
            start=self.start,
            stop=self.stop,
            step=self.step,
            n_frames_total=self.n_frames_total,
        )
        if n_frames_selected == 0:
            raise ValueError("Frame selection must select at least one frame")
        object.__setattr__(self, "n_frames_selected", n_frames_selected)
