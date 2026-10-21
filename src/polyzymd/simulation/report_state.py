"""Keep the portable restart state in step with the trajectory frames.

A production segment writes one trajectory frame every ``report_interval``
steps.  When the job dies, the next segment resumes from a portable state
XML.  Two failures seen on Blanca in September 2026 came from resuming at a
step that did not match the last written frame:

- **Re-simulated boundary step.**  ``restart_state.xml`` used to be written
  only on a wall-clock timer.  After a hard kill the next segment resumed
  from a state up to one timer period older than the last frame, integrated
  forward again and wrote that frame's step a second time from a new
  stochastic branch.
- **Lost frame.**  When a trajectory write raised (``OSError: [Errno 116]
  Stale file handle``), the crash handler saved ``interrupted_state.xml`` at
  the step whose frame had just failed.  The next segment started from that
  step, so its first frame came one interval later and the failed step never
  appeared in any trajectory.

This module provides:

- :class:`ReportedStateTracker`, an OpenMM reporter placed after the DCD,
  state-data and checkpoint reporters.  At every report it atomically
  rewrites ``restart_state.xml`` and records the step of the last frame that
  every reporter wrote.
- :func:`save_state_after_crash`, used by the crash handlers.  It saves
  ``interrupted_state.xml`` only when no report step has passed since the
  last complete report, and otherwise writes the ``INTERRUPTED`` marker
  alone, so the next segment resumes from the report-aligned
  ``restart_state.xml``.
- :func:`read_state_step` and :func:`last_reported_frame`, which read the
  step of a state XML and the step and interval of the last trajectory frame
  of a segment from the DCD header and state-data CSV, without loading
  coordinates.
"""

from __future__ import annotations

import logging
import os
import re
import struct
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)

#: File name of the portable state rewritten at every trajectory report.
RESTART_STATE_NAME = "restart_state.xml"

_STEP_COUNT_RE = re.compile(r'<State\b[^>]*\bstepCount="(\d+)"')


def write_text_atomic(path: Path, text: str) -> None:
    """Write *text* to *path* through a temporary file and ``os.replace``.

    A reader, or a process killed mid-write, sees either the previous
    complete file or the new complete file, never a truncated one.  The
    temporary file sits in the same directory so the rename stays on one
    filesystem.

    Parameters
    ----------
    path : Path
        Destination file.
    text : str
        Content to write.
    """
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def read_state_step(path: Path) -> int | None:
    """Return the ``stepCount`` attribute of an OpenMM state XML.

    Only the first 4 KiB are read; OpenMM writes the attribute on the root
    ``<State>`` element at the top of the file.

    Parameters
    ----------
    path : Path
        Serialized ``openmm.State`` file.

    Returns
    -------
    int or None
        The integrator step count stored in the file, or ``None`` when the
        file cannot be read or carries no ``stepCount`` (OpenMM releases
        before 8.0 do not write it).
    """
    try:
        with Path(path).open("rb") as stream:
            head = stream.read(4096).decode("utf-8", errors="replace")
    except OSError:
        return None
    match = _STEP_COUNT_RE.search(head)
    return int(match.group(1)) if match else None


@dataclass(frozen=True)
class ReportedFrame:
    """Step and spacing of the last trajectory frame of a segment.

    Attributes
    ----------
    segment_index : int
        Segment whose trajectory holds the frame.
    last_step : int
        Integrator step of the last frame.
    report_interval : int or None
        Steps between frames, or ``None`` when the segment holds a single
        frame and no DCD header was readable.
    frames : int
        Number of frames in the segment.
    """

    segment_index: int
    last_step: int
    report_interval: int | None
    frames: int


def dcd_frame_info(dcd_path: Path) -> tuple[int, int, int] | None:
    """Read ``(frames, first_step, interval)`` from an OpenMM DCD header.

    OpenMM's ``DCDFile`` stores the frame count at byte 8, the first step at
    byte 12 and the interval at byte 16, and rewrites the count at every
    frame.  Files with fewer than 100 bytes hold no header.

    Returns
    -------
    tuple of int or None
        ``(frames, first_step, interval)``, or ``None`` when the header is
        missing, unreadable or does not start with the ``CORD`` signature.
    """
    try:
        with dcd_path.open("rb") as stream:
            header = stream.read(100)
    except OSError:
        return None
    if len(header) < 100 or header[4:8] != b"CORD":
        return None
    frames, first_step, interval = struct.unpack("<3i", header[8:20])
    return frames, first_step, interval


def _csv_steps(csv_path: Path) -> list[int]:
    """Return the ``Step`` column of an OpenMM state-data CSV, in file order."""
    steps: list[int] = []
    try:
        with csv_path.open() as stream:
            for line in stream:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or stripped.startswith('"#'):
                    continue
                try:
                    steps.append(int(float(stripped.split(",")[0].strip('"'))))
                except ValueError:
                    continue
    except OSError:
        return []
    return steps


def segment_frames(working_dir: Path, segment_index: int) -> ReportedFrame | None:
    """Return the last trajectory frame written by one production segment.

    The step comes from the ``Step`` column of the state-data CSV.  The DCD
    header gives only the frame count and interval: OpenMM's
    ``DCDReporter`` always writes ``reportInterval`` as the first step, so
    header steps count from the segment start, not from step 0.  OpenMM
    writes the trajectory frame before the state-data row of the same step,
    so when the DCD holds more frames than the CSV has rows, the missing
    rows are added at the frame interval.  When the header is missing or
    was rescaled for steps beyond 2**31 (OpenMM then stores the interval as
    1), the CSV alone is used.

    Parameters
    ----------
    working_dir : Path
        Replicate working directory.
    segment_index : int
        Production segment index.

    Returns
    -------
    ReportedFrame or None
        The last frame, or ``None`` when the state-data CSV holds no row.
    """
    seg_dir = Path(working_dir) / f"production_{segment_index}"
    dcd = dcd_frame_info(seg_dir / f"production_{segment_index}_trajectory.dcd")
    steps = _csv_steps(seg_dir / f"production_{segment_index}_state_data.csv")
    if not steps:
        return None
    if dcd is not None and dcd[0] > 0 and dcd[2] > 1:
        frames, _, interval = dcd
        return ReportedFrame(
            segment_index, steps[-1] + (frames - len(steps)) * interval, interval, frames
        )
    interval = steps[-1] - steps[-2] if len(steps) > 1 else None
    return ReportedFrame(segment_index, steps[-1], interval, len(steps))


def last_reported_frame(working_dir: Path, before_segment: int) -> ReportedFrame | None:
    """Return the last frame written by any segment before *before_segment*.

    Segments whose trajectory holds no frame (interrupted at start-up) are
    passed over, so the result is the frame the next segment must continue
    from.

    Parameters
    ----------
    working_dir : Path
        Replicate working directory.
    before_segment : int
        Segment about to run; segments ``before_segment - 1`` down to 0 are
        searched.

    Returns
    -------
    ReportedFrame or None
        The most recent frame, or ``None`` when no earlier segment wrote one.
    """
    for index in range(before_segment - 1, -1, -1):
        frame = segment_frames(working_dir, index)
        if frame is not None:
            return frame
    return None


def first_report_step(step: int, report_interval: int) -> int:
    """Return the first step above *step* at which OpenMM reporters fire.

    OpenMM reporters fire when the step count is a multiple of the interval,
    and a simulation that starts exactly on a multiple does not report that
    step.
    """
    return (step // report_interval + 1) * report_interval


class SkipReportsThrough:
    """Wrap an OpenMM reporter so it writes nothing at or before *last_step*.

    A segment that resumes from a state older than the previous segment's
    last frame integrates those report steps again; wrapping its reporters
    keeps the trajectory at one frame per step.

    Parameters
    ----------
    reporter : object
        OpenMM reporter to wrap.
    last_step : int
        Step of the last frame the trajectory already holds.
    """

    def __init__(self, reporter: Any, last_step: int) -> None:
        self._reporter = reporter
        self._last_step = int(last_step)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._reporter, name)

    def describeNextReport(self, simulation: Any) -> Any:
        return self._reporter.describeNextReport(simulation)

    def report(self, simulation: Any, state: Any) -> None:
        if simulation.currentStep > self._last_step:
            self._reporter.report(simulation, state)


class ReportedStateTracker:
    """OpenMM reporter that writes ``restart_state.xml`` at every frame.

    Append it to ``simulation.reporters`` after the DCD, state-data and
    checkpoint reporters.  OpenMM calls reporters in list order and stops at
    the first exception, so :meth:`report` runs only after every earlier
    reporter wrote the frame of that step.  It then serializes positions,
    velocities and context parameters to ``restart_state.xml`` through a
    temporary file and ``os.replace``, and records the step.

    The wall-clock restart checkpoint goes through :meth:`save_restart` so
    the tracker also knows the step of that file.

    Parameters
    ----------
    output_dir : Path
        Segment directory (``production_N/``).
    report_interval : int
        Steps between trajectory frames; must match the DCD reporter.
    start_step : int
        Integrator step count when the segment began.
    """

    def __init__(self, output_dir: Path, report_interval: int, start_step: int) -> None:
        if report_interval <= 0:
            raise ValueError("report_interval must be a positive integer")
        self._output_dir = Path(output_dir)
        self._report_interval = int(report_interval)
        self.start_step = int(start_step)
        #: Step of the last frame every reporter wrote, or None before the first.
        self.last_reported_step: int | None = None
        #: Frames written by this segment.
        self.frames_written = 0
        #: Step of the state currently in ``restart_state.xml``, or None.
        self.restart_state_step: int | None = None

    @property
    def report_interval(self) -> int:
        """Steps between trajectory frames."""
        return self._report_interval

    @property
    def restart_state_path(self) -> Path:
        """Path of the portable state kept in step with the frames."""
        return self._output_dir / RESTART_STATE_NAME

    def describeNextReport(self, simulation: Any) -> tuple[int, bool, bool, bool, bool, bool]:
        """Request a report at the next multiple of the report interval.

        Uses the six-element tuple form, which OpenMM 8.1 through 8.4 accept.
        No state data is requested because :meth:`report` reads the context
        directly.
        """
        steps = self._report_interval - simulation.currentStep % self._report_interval
        return (steps, False, False, False, False, False)

    def report(self, simulation: Any, state: Any) -> None:
        """Write ``restart_state.xml`` for the frame just written.

        Parameters
        ----------
        simulation : openmm.app.Simulation
            Simulation being reported.
        state : openmm.State
            State passed by OpenMM (unused; positions and velocities are read
            from the context with parameters included).
        """
        from openmm import XmlSerializer

        step = int(simulation.currentStep)
        # The earlier reporters already wrote this frame, so record it before
        # the state write: if that write fails, the crash handler still knows
        # the frame exists and saves the interrupted state at this step.
        self.last_reported_step = step
        self.frames_written += 1
        full_state = simulation.context.getState(
            getPositions=True, getVelocities=True, getParameters=True
        )
        write_text_atomic(self.restart_state_path, XmlSerializer.serialize(full_state))
        self.restart_state_step = step

    def save_restart(self, simulation: Any) -> None:
        """Write the wall-clock restart checkpoint and record its step.

        Parameters
        ----------
        simulation : openmm.app.Simulation
            Simulation to checkpoint.
        """
        from polyzymd.simulation.signals import save_restart_checkpoint

        step = int(simulation.currentStep)
        save_restart_checkpoint(simulation=simulation, output_dir=self._output_dir)
        self.restart_state_step = step

    def next_unreported_step(self) -> int:
        """Return the first report step whose frame this segment has not written."""
        base = self.last_reported_step if self.last_reported_step is not None else self.start_step
        return first_report_step(base, self._report_interval)

    def frame_fields(self) -> dict[str, int | None]:
        """Return the frame bookkeeping recorded on the segment in ``progress.json``."""
        return {
            "report_interval": self._report_interval,
            "start_step": self.start_step,
            "last_reported_step": self.last_reported_step,
            "samples_written": self.frames_written,
        }


def save_state_after_crash(
    simulation: Any,
    tracker: ReportedStateTracker,
    output_dir: Path,
    segment_index: int,
    total_steps: int,
) -> int:
    """Save what the next segment needs after an unexpected exception.

    The integrator may have passed a report step whose frame was not
    written, for example when the DCD write raised ``OSError``.  Saving the
    context there would make the next segment skip that frame.  So:

    - If the context step is below the next unreported report step, every
      frame up to the current step exists and the full interrupted state is
      saved with :func:`~polyzymd.simulation.signals.save_interrupted_state`.
    - Otherwise no ``interrupted_state.xml`` is written.  The ``INTERRUPTED``
      marker records the steps up to the state in ``restart_state.xml``,
      which the tracker wrote at the last complete frame (or the wall-clock
      checkpoint, if that is later but still before the lost frame), and the
      next segment resumes from that file and re-simulates the lost frame.

    Parameters
    ----------
    simulation : openmm.app.Simulation
        Simulation that raised.
    tracker : ReportedStateTracker
        Tracker attached to *simulation*.
    output_dir : Path
        Segment directory.
    segment_index : int
        Segment index, written to the marker.
    total_steps : int
        Steps planned for the segment, written to the marker.

    Returns
    -------
    int
        Steps from the segment start to the state the next segment will load.
    """
    from polyzymd.simulation.signals import save_interrupted_state, write_interrupted_marker

    current_step = int(simulation.context.getStepCount())
    lost_step = tracker.next_unreported_step()
    if current_step < lost_step:
        steps_completed = current_step - tracker.start_step
        save_interrupted_state(
            simulation=simulation,
            output_dir=output_dir,
            segment_index=segment_index,
            steps_completed=steps_completed,
            total_steps=total_steps,
        )
        tracker.restart_state_step = current_step
        return steps_completed

    resume_step = tracker.restart_state_step
    if resume_step is None or resume_step >= lost_step:
        resume_step = tracker.start_step
    steps_completed = resume_step - tracker.start_step
    LOGGER.error(
        f"Segment {segment_index}: the frame at step {lost_step} was not written "
        f"(integrator at step {current_step}); not saving interrupted_state.xml so the "
        f"next segment resumes from {RESTART_STATE_NAME} at step {resume_step} and "
        f"re-simulates that frame"
    )
    write_interrupted_marker(
        output_dir=output_dir,
        segment_index=segment_index,
        steps_completed=steps_completed,
        total_steps=total_steps,
    )
    return steps_completed
