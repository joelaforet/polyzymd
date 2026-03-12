"""
Progress tracking for self-resubmitting MD simulation jobs.

Tracks completed simulation segments across multiple SLURM job submissions,
enabling idempotent restart: each job scans the filesystem to determine what
work remains, runs until wall-time or completion, and resubmits itself if
the total requested simulation time has not been reached.

The progress file (``progress.json``) is the primary source of truth, with
filesystem scanning (``production_N/`` directories) used for validation on
startup. Writes are atomic (write to ``.tmp``, then rename) to prevent
corruption if the process is killed mid-write.
"""

from __future__ import annotations

import json
import logging
import os
import re
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List

from pydantic import BaseModel, Field

LOGGER = logging.getLogger(__name__)


class SegmentStatus(str, Enum):
    """Status of a simulation segment."""

    COMPLETED = "completed"
    INTERRUPTED = "interrupted"
    FAILED = "failed"
    RUNNING = "running"


class SimulationStatus(str, Enum):
    """Overall status of the simulation."""

    NOT_STARTED = "not_started"
    RUNNING = "running"
    COMPLETED = "completed"
    INTERRUPTED = "interrupted"
    FAILED = "failed"


class EquilibrationStageRecord(BaseModel):
    """Record of a single equilibration stage.

    Used for observability in ``progress.json`` and for detecting
    completed stages when resuming after a wall-time interruption.

    Parameters
    ----------
    index : int
        Stage index (0-based).
    name : str
        Stage name (e.g. ``"heating"``, ``"polymer_relaxation"``).
    status : SegmentStatus
        Current status of this stage.
    duration_ns : float
        Planned simulation time for this stage in nanoseconds.
    ensemble : str
        Thermodynamic ensemble (``"NVT"`` or ``"NPT"``).
    started_at : str
        ISO-format timestamp when the stage started.
    finished_at : str | None
        ISO-format timestamp when the stage finished.
    """

    index: int
    name: str = ""
    status: SegmentStatus = SegmentStatus.COMPLETED
    duration_ns: float = 0.0
    ensemble: str = "NVT"
    started_at: str = Field(default_factory=lambda: _now_iso())
    finished_at: str | None = None


class SegmentRecord(BaseModel):
    """Record of a single production segment.

    Parameters
    ----------
    index : int
        Segment index (0-based).
    steps_completed : int
        Number of MD steps completed in this segment.
    steps_requested : int
        Number of MD steps originally requested for this segment.
    samples_written : int
        Number of trajectory frames written.
    started_at : str
        ISO-format timestamp when the segment started.
    finished_at : str | None
        ISO-format timestamp when the segment finished (None if interrupted/running).
    status : SegmentStatus
        Current status of this segment.
    duration_ns : float
        Simulation time covered by this segment in nanoseconds.
    """

    index: int
    steps_completed: int = 0
    steps_requested: int = 0
    samples_written: int = 0
    started_at: str = Field(default_factory=lambda: _now_iso())
    finished_at: str | None = None
    status: SegmentStatus = SegmentStatus.RUNNING
    duration_ns: float = 0.0


class SimulationProgress(BaseModel):
    """Progress state for a self-resubmitting simulation.

    This is the primary data model written to ``progress.json`` in the
    simulation working directory. It tracks all completed segments and
    provides the information needed for each job to determine what work
    remains.

    Parameters
    ----------
    config_path : str
        Path to the YAML configuration file.
    total_steps_requested : int
        Total number of MD steps for the entire production run.
    total_samples_requested : int
        Total number of trajectory frames for the entire production run.
    timestep_fs : float
        Integration timestep in femtoseconds.
    segments : list of SegmentRecord
        Records for each completed/interrupted segment.
    status : SimulationStatus
        Overall simulation status.
    last_updated : str
        ISO-format timestamp of the last update.
    replicate : int
        Replicate index (1-based).
    """

    config_path: str = ""
    total_steps_requested: int = 0
    total_samples_requested: int = 0
    timestep_fs: float = 2.0
    equilibration_stages: List[EquilibrationStageRecord] = Field(default_factory=list)
    segments: List[SegmentRecord] = Field(default_factory=list)
    status: SimulationStatus = SimulationStatus.NOT_STARTED
    last_updated: str = Field(default_factory=lambda: _now_iso())
    replicate: int = 1

    @property
    def equilibration_complete(self) -> bool:
        """Whether all recorded equilibration stages completed successfully."""
        return len(self.equilibration_stages) > 0 and all(
            s.status == SegmentStatus.COMPLETED for s in self.equilibration_stages
        )

    @property
    def num_eq_stages_completed(self) -> int:
        """Number of equilibration stages that completed successfully."""
        return sum(1 for s in self.equilibration_stages if s.status == SegmentStatus.COMPLETED)

    @property
    def total_steps_completed(self) -> int:
        """Total MD steps completed across all segments (including interrupted)."""
        return sum(
            seg.steps_completed
            for seg in self.segments
            if seg.status in (SegmentStatus.COMPLETED, SegmentStatus.INTERRUPTED)
        )

    @property
    def total_samples_completed(self) -> int:
        """Total trajectory frames written across all completed segments."""
        return sum(
            seg.samples_written for seg in self.segments if seg.status == SegmentStatus.COMPLETED
        )

    @property
    def next_segment_index(self) -> int:
        """Index for the next segment to run.

        If the last segment was interrupted, the next segment gets a fresh
        index (the interrupted segment's output stays in its directory;
        the new segment resumes from the last completed state).
        """
        if not self.segments:
            return 0
        return max(seg.index for seg in self.segments) + 1

    @property
    def steps_remaining(self) -> int:
        """Number of MD steps still needed to reach the total."""
        return max(0, self.total_steps_requested - self.total_steps_completed)

    @property
    def is_complete(self) -> bool:
        """Whether the simulation has reached its total requested steps."""
        return self.total_steps_completed >= self.total_steps_requested

    @property
    def last_completed_segment(self) -> SegmentRecord | None:
        """The most recent segment that completed successfully."""
        completed = [s for s in self.segments if s.status == SegmentStatus.COMPLETED]
        if not completed:
            return None
        return max(completed, key=lambda s: s.index)

    @property
    def last_segment(self) -> SegmentRecord | None:
        """The most recent segment regardless of status."""
        if not self.segments:
            return None
        return max(self.segments, key=lambda s: s.index)

    def fraction_complete(self) -> float:
        """Fraction of the total simulation completed (0.0 to 1.0)."""
        if self.total_steps_requested == 0:
            return 0.0
        return min(1.0, self.total_steps_completed / self.total_steps_requested)

    def time_completed_ns(self) -> float:
        """Total simulation time completed in nanoseconds."""
        return sum(
            seg.duration_ns for seg in self.segments if seg.status == SegmentStatus.COMPLETED
        )


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------

PROGRESS_FILENAME = "progress.json"


def _now_iso() -> str:
    """Return current UTC time as ISO-format string."""
    return datetime.now(timezone.utc).isoformat()


def _progress_path(working_dir: Path) -> Path:
    """Return the path to the progress file."""
    return working_dir / PROGRESS_FILENAME


def load_progress(working_dir: str | Path) -> SimulationProgress | None:
    """Load progress from ``progress.json`` in *working_dir*.

    Parameters
    ----------
    working_dir : str or Path
        Simulation working directory.

    Returns
    -------
    SimulationProgress or None
        The progress state, or ``None`` if no progress file exists.
    """
    path = _progress_path(Path(working_dir))
    if not path.exists():
        return None

    try:
        with open(path, "r") as f:
            data = json.load(f)
        progress = SimulationProgress.model_validate(data)
        LOGGER.info(
            f"Loaded progress: {progress.total_steps_completed}/{progress.total_steps_requested} "
            f"steps ({progress.fraction_complete():.1%}), "
            f"{len(progress.segments)} segment(s)"
        )
        return progress
    except (json.JSONDecodeError, Exception) as exc:
        LOGGER.warning(f"Failed to load progress from {path}: {exc}")
        return None


def save_progress(working_dir: str | Path, progress: SimulationProgress) -> Path:
    """Save progress to ``progress.json`` atomically.

    Writes to a temporary file first, then renames (atomic on POSIX).

    Parameters
    ----------
    working_dir : str or Path
        Simulation working directory.
    progress : SimulationProgress
        Progress state to save.

    Returns
    -------
    Path
        Path to the saved progress file.
    """
    working_dir = Path(working_dir)
    working_dir.mkdir(parents=True, exist_ok=True)

    progress.last_updated = _now_iso()
    target = _progress_path(working_dir)
    tmp = target.with_suffix(".json.tmp")

    with open(tmp, "w") as f:
        json.dump(progress.model_dump(mode="json"), f, indent=2)

    os.replace(str(tmp), str(target))
    LOGGER.debug(f"Saved progress to {target}")
    return target


# ---------------------------------------------------------------------------
# Filesystem scanning
# ---------------------------------------------------------------------------

_PRODUCTION_DIR_RE = re.compile(r"^production_(\d+)$")
_EQ_STAGE_DIR_RE = re.compile(r"^equilibration_(\d+)_(.+)$")


def scan_equilibration_stages(working_dir: str | Path) -> List[EquilibrationStageRecord]:
    """Scan filesystem for completed equilibration stage directories.

    A stage is considered completed if its checkpoint file
    (``equilibration_N_name_checkpoint.chk``) exists. Directories
    without a checkpoint are recorded with ``FAILED`` status.

    Parameters
    ----------
    working_dir : str or Path
        Simulation working directory.

    Returns
    -------
    list of EquilibrationStageRecord
        Records sorted by stage index.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        return []

    records: List[EquilibrationStageRecord] = []
    for entry in sorted(working_dir.iterdir()):
        if not entry.is_dir():
            continue
        match = _EQ_STAGE_DIR_RE.match(entry.name)
        if not match:
            continue
        stage_idx = int(match.group(1))
        stage_name = match.group(2)
        chk = entry / f"equilibration_{stage_idx}_{stage_name}_checkpoint.chk"
        status = SegmentStatus.COMPLETED if chk.exists() else SegmentStatus.FAILED
        records.append(
            EquilibrationStageRecord(
                index=stage_idx,
                name=stage_name,
                status=status,
            )
        )

    records.sort(key=lambda r: r.index)
    LOGGER.debug(
        f"Equilibration scan: {sum(1 for r in records if r.status == SegmentStatus.COMPLETED)}"
        f"/{len(records)} stage(s) completed"
    )
    return records


def scan_filesystem(
    working_dir: str | Path,
    timestep_fs: float = 2.0,
) -> SimulationProgress:
    """Reconstruct progress by scanning ``production_N/`` directories.

    This is the validation / recovery path. It examines the filesystem to
    determine which segments completed, which were interrupted, and how
    many total steps have been completed.

    Parameters
    ----------
    working_dir : str or Path
        Simulation working directory.
    timestep_fs : float
        Integration timestep in femtoseconds (needed for step→ns conversion).

    Returns
    -------
    SimulationProgress
        Reconstructed progress state.
    """
    working_dir = Path(working_dir)
    segments: List[SegmentRecord] = []

    # Find all production_N/ directories
    prod_dirs: List[tuple[int, Path]] = []
    for entry in sorted(working_dir.iterdir()):
        if entry.is_dir():
            match = _PRODUCTION_DIR_RE.match(entry.name)
            if match:
                prod_dirs.append((int(match.group(1)), entry))

    for seg_idx, seg_dir in sorted(prod_dirs):
        record = _scan_segment_dir(seg_idx, seg_dir, timestep_fs)
        if record is not None:
            segments.append(record)

    progress = SimulationProgress(
        timestep_fs=timestep_fs,
        equilibration_stages=scan_equilibration_stages(working_dir),
        segments=segments,
    )

    # Determine overall status
    if not segments:
        progress.status = SimulationStatus.NOT_STARTED
    elif any(s.status == SegmentStatus.INTERRUPTED for s in segments):
        progress.status = SimulationStatus.INTERRUPTED
    elif any(s.status == SegmentStatus.FAILED for s in segments):
        progress.status = SimulationStatus.FAILED
    else:
        progress.status = SimulationStatus.RUNNING

    LOGGER.info(
        f"Filesystem scan: {progress.num_eq_stages_completed} eq stage(s), "
        f"{len(segments)} production segment(s), "
        f"{progress.total_steps_completed} total steps completed"
    )
    return progress


def _scan_segment_dir(
    seg_idx: int,
    seg_dir: Path,
    timestep_fs: float,
) -> SegmentRecord | None:
    """Scan a single ``production_N/`` directory to build a SegmentRecord.

    Parameters
    ----------
    seg_idx : int
        Segment index.
    seg_dir : Path
        Path to the segment directory.
    timestep_fs : float
        Integration timestep in femtoseconds.

    Returns
    -------
    SegmentRecord or None
        Record for this segment, or None if the directory is empty/invalid.
    """
    interrupted_marker = seg_dir / "INTERRUPTED"
    state_xml = seg_dir / f"production_{seg_idx}_state.xml"
    params_json = seg_dir / f"production_{seg_idx}_parameters.json"

    # Determine status
    if interrupted_marker.exists():
        status = SegmentStatus.INTERRUPTED
        # Read metadata from the INTERRUPTED marker
        steps_completed, total_steps = _parse_interrupted_marker(interrupted_marker)
        duration_ns = (steps_completed * timestep_fs) / 1e6
        return SegmentRecord(
            index=seg_idx,
            steps_completed=steps_completed,
            steps_requested=total_steps,
            samples_written=0,  # Unknown for interrupted segments
            status=status,
            duration_ns=duration_ns,
        )
    elif state_xml.exists():
        # Completed segment — read parameters to get step count
        status = SegmentStatus.COMPLETED
        steps_completed = 0
        samples_written = 0
        duration_ns = 0.0

        if params_json.exists():
            try:
                with open(params_json, "r") as f:
                    params = json.load(f)
                integ = params["__values__"]["integ_params"]["__values__"]
                total_time_val = integ["total_time"]["__values__"]["value"]
                time_unit = integ["total_time"]["__values__"]["unit"]
                ts_val = integ["time_step"]["__values__"]["value"]
                ts_unit = integ["time_step"]["__values__"]["unit"]
                num_samples = integ.get("num_samples", 0)

                # Convert to consistent units (ns and fs)
                duration_ns = _convert_to_ns(total_time_val, time_unit)
                ts_fs = _convert_to_fs(ts_val, ts_unit)
                steps_completed = int(duration_ns * 1e6 / ts_fs)
                samples_written = num_samples
            except (KeyError, ValueError, json.JSONDecodeError) as exc:
                LOGGER.warning(f"Could not parse parameters from {params_json}: {exc}")

        return SegmentRecord(
            index=seg_idx,
            steps_completed=steps_completed,
            steps_requested=steps_completed,
            samples_written=samples_written,
            status=status,
            duration_ns=duration_ns,
            finished_at=_now_iso(),  # Approximate
        )
    else:
        # No state.xml and no INTERRUPTED marker — check for checkpoint
        # from periodic CheckpointReporter (hard-kill recovery).
        checkpoint_chk = seg_dir / f"production_{seg_idx}_checkpoint.chk"
        state_data_csv = seg_dir / f"production_{seg_idx}_state_data.csv"

        if checkpoint_chk.exists():
            # Hard-killed segment: periodic checkpoint survives but the
            # signal handler never fired.  Treat as interrupted so the
            # continuation manager can recover from the checkpoint.
            steps_completed = (
                _estimate_steps_from_csv(state_data_csv) if state_data_csv.exists() else 0
            )
            duration_ns = (steps_completed * timestep_fs) / 1e6
            LOGGER.warning(
                f"Segment {seg_idx} directory has checkpoint but no state.xml "
                f"or INTERRUPTED marker — treating as interrupted "
                f"(hard-kill recovery, ~{steps_completed} steps)"
            )
            return SegmentRecord(
                index=seg_idx,
                steps_completed=steps_completed,
                steps_requested=0,  # Unknown — no INTERRUPTED metadata
                samples_written=0,
                status=SegmentStatus.INTERRUPTED,
                duration_ns=duration_ns,
            )
        else:
            # Truly failed: no recoverable files at all
            LOGGER.warning(
                f"Segment {seg_idx} directory exists but has no state.xml, "
                f"INTERRUPTED marker, or checkpoint — treating as failed"
            )
            return SegmentRecord(
                index=seg_idx,
                steps_completed=0,
                steps_requested=0,
                status=SegmentStatus.FAILED,
            )


def _parse_interrupted_marker(marker_path: Path) -> tuple[int, int]:
    """Parse an INTERRUPTED marker file for step counts.

    Parameters
    ----------
    marker_path : Path
        Path to the INTERRUPTED marker file.

    Returns
    -------
    tuple of (steps_completed, total_steps)
    """
    steps_completed = 0
    total_steps = 0

    try:
        text = marker_path.read_text()
        for line in text.strip().splitlines():
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            if key == "steps_completed":
                steps_completed = int(value)
            elif key == "total_steps":
                total_steps = int(value)
    except (ValueError, OSError) as exc:
        LOGGER.warning(f"Could not parse INTERRUPTED marker {marker_path}: {exc}")

    return steps_completed, total_steps


def _estimate_steps_from_csv(csv_path: Path) -> int:
    """Estimate the number of completed steps from a ``state_data.csv`` file.

    Reads the last non-empty data line and extracts the step count from
    the first column.  This is used to estimate progress for hard-killed
    segments that have no ``INTERRUPTED`` marker or ``state.xml`` but do
    have reporter output on disk.

    Parameters
    ----------
    csv_path : Path
        Path to the ``production_N_state_data.csv`` file.

    Returns
    -------
    int
        Estimated steps completed, or 0 if the file cannot be parsed.
    """
    try:
        with open(csv_path, "r") as f:
            lines = f.readlines()
        # Need at least a header + one data line
        if len(lines) < 2:
            return 0
        # Walk backwards to find the last non-empty line
        for line in reversed(lines):
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                # First column is the step count (may be quoted)
                step_str = stripped.split(",")[0].strip('"').strip()
                return int(float(step_str))
        return 0
    except (ValueError, IndexError, OSError) as exc:
        LOGGER.warning(f"Could not estimate steps from {csv_path}: {exc}")
        return 0


def _convert_to_ns(value: float, unit: str) -> float:
    """Convert a time value to nanoseconds."""
    unit = unit.lower().rstrip("s")
    if unit == "nanosecond":
        return value
    elif unit == "picosecond":
        return value / 1000.0
    elif unit == "femtosecond":
        return value / 1e6
    elif unit == "microsecond":
        return value * 1000.0
    else:
        LOGGER.warning(f"Unknown time unit '{unit}', assuming nanoseconds")
        return value


def _convert_to_fs(value: float, unit: str) -> float:
    """Convert a time value to femtoseconds."""
    unit = unit.lower().rstrip("s")
    if unit == "femtosecond":
        return value
    elif unit == "picosecond":
        return value * 1000.0
    elif unit == "nanosecond":
        return value * 1e6
    else:
        LOGGER.warning(f"Unknown time unit '{unit}', assuming femtoseconds")
        return value


# ---------------------------------------------------------------------------
# Progress helpers
# ---------------------------------------------------------------------------


def get_next_segment_info(
    progress: SimulationProgress,
    total_steps: int,
    total_samples: int,
) -> Dict[str, Any] | None:
    """Determine what the next segment should run.

    Uses a fixed report interval derived from the *global* config
    (``total_steps // total_samples``) so that every segment writes
    frames at the same cadence. This guarantees uniform time spacing
    across segments, enabling seamless trajectory concatenation.

    Parameters
    ----------
    progress : SimulationProgress
        Current progress state.
    total_steps : int
        Total steps requested for the full simulation.
    total_samples : int
        Total trajectory frames requested for the full simulation.

    Returns
    -------
    dict or None
        Dictionary with ``segment_index``, ``steps_to_run``,
        ``report_interval``, and ``samples_to_write`` for the next
        segment. Returns ``None`` if the simulation is already complete.
    """
    completed = progress.total_steps_completed
    remaining = total_steps - completed

    if remaining <= 0:
        return None

    # Fixed interval from the global config — same for every segment
    report_interval = max(1, total_steps // total_samples)
    samples = remaining // report_interval

    return {
        "segment_index": progress.next_segment_index,
        "steps_to_run": remaining,
        "report_interval": report_interval,
        "samples_to_write": samples,
    }


def validate_progress(
    working_dir: str | Path,
    progress: SimulationProgress,
    timestep_fs: float = 2.0,
) -> SimulationProgress:
    """Validate progress file against filesystem and reconcile.

    Loads the progress file and cross-checks it against the actual
    ``production_N/`` directories on disk. If discrepancies are found,
    the filesystem is treated as authoritative for segment status
    (completed vs interrupted), but the progress file's metadata
    (config_path, total_steps_requested, etc.) is preserved.

    Parameters
    ----------
    working_dir : str or Path
        Simulation working directory.
    progress : SimulationProgress
        Progress loaded from the progress file.
    timestep_fs : float
        Integration timestep in femtoseconds.

    Returns
    -------
    SimulationProgress
        Validated (and possibly corrected) progress state.
    """
    working_dir = Path(working_dir)
    fs_progress = scan_filesystem(working_dir, timestep_fs=timestep_fs)

    # Build a lookup of filesystem segment records by index
    fs_segments: Dict[int, SegmentRecord] = {s.index: s for s in fs_progress.segments}
    file_segments: Dict[int, SegmentRecord] = {s.index: s for s in progress.segments}

    reconciled: List[SegmentRecord] = []
    all_indices = sorted(set(fs_segments.keys()) | set(file_segments.keys()))

    for idx in all_indices:
        fs_rec = fs_segments.get(idx)
        file_rec = file_segments.get(idx)

        if fs_rec is not None and file_rec is not None:
            # Both exist — filesystem status is authoritative
            if fs_rec.status != file_rec.status:
                LOGGER.warning(
                    f"Segment {idx}: progress file says {file_rec.status.value}, "
                    f"filesystem says {fs_rec.status.value} — using filesystem"
                )
            # Use filesystem status but keep file_rec metadata if richer
            merged = SegmentRecord(
                index=idx,
                steps_completed=max(fs_rec.steps_completed, file_rec.steps_completed),
                steps_requested=max(fs_rec.steps_requested, file_rec.steps_requested),
                samples_written=max(fs_rec.samples_written, file_rec.samples_written),
                started_at=file_rec.started_at,
                finished_at=file_rec.finished_at or fs_rec.finished_at,
                status=fs_rec.status,
                duration_ns=max(fs_rec.duration_ns, file_rec.duration_ns),
            )
            reconciled.append(merged)
        elif fs_rec is not None:
            # Only on filesystem — progress file is missing this segment
            LOGGER.warning(f"Segment {idx}: found on filesystem but not in progress file — adding")
            reconciled.append(fs_rec)
        else:
            # Only in progress file — filesystem doesn't have it
            # This could happen if the directory was cleaned up
            assert file_rec is not None
            LOGGER.warning(
                f"Segment {idx}: in progress file but not on filesystem — keeping record"
            )
            reconciled.append(file_rec)

    progress.segments = reconciled

    # Recompute overall status
    if progress.is_complete:
        progress.status = SimulationStatus.COMPLETED
    elif any(s.status == SegmentStatus.INTERRUPTED for s in reconciled):
        progress.status = SimulationStatus.INTERRUPTED
    elif any(s.status == SegmentStatus.FAILED for s in reconciled):
        progress.status = SimulationStatus.FAILED
    elif reconciled:
        progress.status = SimulationStatus.RUNNING
    else:
        progress.status = SimulationStatus.NOT_STARTED

    LOGGER.info(
        f"Validated progress: {progress.total_steps_completed}/{progress.total_steps_requested} "
        f"steps ({progress.fraction_complete():.1%})"
    )
    return progress


def load_or_scan_progress(
    working_dir: str | Path,
    config_path: str = "",
    total_steps: int = 0,
    total_samples: int = 0,
    timestep_fs: float = 2.0,
    replicate: int = 1,
) -> SimulationProgress:
    """Load progress from file, falling back to filesystem scan.

    This is the primary entry point for jobs to determine current state.
    It loads ``progress.json`` if it exists, validates it against the
    filesystem, and returns the reconciled progress. If no progress file
    exists, it scans the filesystem and creates a new progress record.

    Parameters
    ----------
    working_dir : str or Path
        Simulation working directory.
    config_path : str
        Path to the YAML configuration file.
    total_steps : int
        Total steps requested (used when creating new progress).
    total_samples : int
        Total samples requested (used when creating new progress).
    timestep_fs : float
        Integration timestep in femtoseconds.
    replicate : int
        Replicate index (1-based).

    Returns
    -------
    SimulationProgress
        Current progress state, validated against filesystem.
    """
    working_dir = Path(working_dir)
    progress = load_progress(working_dir)

    if progress is not None:
        # Validate against filesystem
        progress = validate_progress(working_dir, progress, timestep_fs=timestep_fs)
        # Ensure metadata is set (may be missing from old progress files)
        if not progress.config_path and config_path:
            progress.config_path = config_path
        # Always sync total_steps/samples from config — the user may have
        # changed the production duration since the progress file was written.
        if total_steps > 0:
            progress.total_steps_requested = total_steps
        if total_samples > 0:
            progress.total_samples_requested = total_samples
        return progress

    # No progress file — scan filesystem
    LOGGER.info("No progress file found — scanning filesystem")
    progress = scan_filesystem(working_dir, timestep_fs=timestep_fs)
    progress.config_path = config_path
    progress.total_steps_requested = total_steps
    progress.total_samples_requested = total_samples
    progress.timestep_fs = timestep_fs
    progress.replicate = replicate

    # Check completion
    if progress.is_complete:
        progress.status = SimulationStatus.COMPLETED

    return progress
