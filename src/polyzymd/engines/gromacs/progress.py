"""GROMACS-specific progress tracking.

Writes and reads ``progress.json`` using the same model as OpenMM
(``SimulationProgress`` from ``polyzymd.simulation.progress``), adapted
for GROMACS flat directory layout.
"""

from __future__ import annotations

import logging
import re
from datetime import datetime, timezone
from pathlib import Path

from polyzymd.simulation.progress import (
    EquilibrationStageRecord,
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    SimulationStatus,
    load_progress,
    save_progress,
)

LOGGER = logging.getLogger(__name__)


def scan_gromacs_progress(
    working_dir: Path,
    config_path: str = "",
    replicate: int = 1,
    total_steps: int = 0,
    total_samples: int = 0,
    timestep_fs: float = 2.0,
) -> SimulationProgress:
    """Scan a GROMACS working directory for simulation progress.

    GROMACS uses a flat layout (no ``production_N/`` subdirectories).
    Progress is inferred from:

    - Equilibration: presence of ``eq_NN.gro`` files
    - Production: parsing ``prod.log`` for completed step counts
    - Completion: ``"Finished mdrun"`` marker in ``prod.log``

    Parameters
    ----------
    working_dir : Path
        GROMACS working directory.
    config_path : str, optional
        Source config path for metadata.
    replicate : int, optional
        Replicate index for metadata.
    total_steps : int, optional
        Total production steps from config.
    total_samples : int, optional
        Total trajectory samples from config.
    timestep_fs : float, optional
        Integration time step in femtoseconds.

    Returns
    -------
    SimulationProgress
        Reconstructed progress model.
    """
    working_dir = Path(working_dir)
    equilibration_stages = _scan_equilibration_gromacs(working_dir)

    log_info = _parse_gromacs_log(working_dir / "prod.log")
    steps_completed = int(log_info["steps_completed"])
    time_completed_ps = float(log_info["time_completed_ps"])
    nsteps_requested = int(log_info["nsteps_requested"])
    is_finished = bool(log_info["is_finished"])

    requested_steps = total_steps if total_steps > 0 else nsteps_requested
    if requested_steps <= 0:
        requested_steps = steps_completed

    segments: list[SegmentRecord] = []
    if steps_completed > 0 or is_finished:
        segment_status = SegmentStatus.COMPLETED if is_finished else SegmentStatus.INTERRUPTED
        duration_ns = max(time_completed_ps / 1000.0, (steps_completed * timestep_fs) / 1e6)
        segments.append(
            SegmentRecord(
                index=0,
                steps_completed=steps_completed,
                steps_requested=requested_steps,
                samples_written=0,
                status=segment_status,
                duration_ns=duration_ns,
                finished_at=datetime.now(timezone.utc).isoformat() if is_finished else None,
            )
        )

    progress = SimulationProgress(
        config_path=config_path,
        total_steps_requested=requested_steps,
        total_samples_requested=total_samples,
        timestep_fs=timestep_fs,
        equilibration_stages=equilibration_stages,
        segments=segments,
        replicate=replicate,
    )

    if is_finished or (requested_steps > 0 and progress.total_steps_completed >= requested_steps):
        progress.status = SimulationStatus.COMPLETED
    elif steps_completed > 0:
        progress.status = SimulationStatus.INTERRUPTED
    elif equilibration_stages:
        progress.status = SimulationStatus.RUNNING
    else:
        progress.status = SimulationStatus.NOT_STARTED

    return progress


def update_gromacs_progress(
    working_dir: Path,
    config_path: str = "",
    replicate: int = 1,
    mark_complete: bool = False,
) -> SimulationProgress:
    """Update ``progress.json`` for a GROMACS simulation.

    This function is designed to be called by the GROMACS SLURM wrapper after
    each ``mdrun`` invocation. It merges the latest ``prod.log`` scan into
    existing progress and records each restart as a new segment entry.

    Parameters
    ----------
    working_dir : Path
        GROMACS working directory.
    config_path : str, optional
        Simulation config path.
    replicate : int, optional
        Replicate index.
    mark_complete : bool, optional
        Force completion status after post-processing.

    Returns
    -------
    SimulationProgress
        Updated and saved progress state.
    """
    working_dir = Path(working_dir)
    existing = load_progress(working_dir)

    scanned = scan_gromacs_progress(
        working_dir=working_dir,
        config_path=config_path,
        replicate=replicate,
        total_steps=existing.total_steps_requested if existing else 0,
        total_samples=existing.total_samples_requested if existing else 0,
        timestep_fs=existing.timestep_fs if existing else 2.0,
    )

    if existing is None:
        progress = scanned
    else:
        progress = existing
        progress.equilibration_stages = scanned.equilibration_stages
        if config_path:
            progress.config_path = config_path
        progress.replicate = replicate

        if progress.total_steps_requested <= 0 and scanned.total_steps_requested > 0:
            progress.total_steps_requested = scanned.total_steps_requested
        if progress.total_samples_requested <= 0 and scanned.total_samples_requested > 0:
            progress.total_samples_requested = scanned.total_samples_requested

        old_steps = progress.total_steps_completed
        new_steps = scanned.total_steps_completed
        delta_steps = max(0, new_steps - old_steps)

        if delta_steps > 0:
            status = (
                SegmentStatus.COMPLETED
                if scanned.status == SimulationStatus.COMPLETED
                else SegmentStatus.INTERRUPTED
            )
            progress.segments.append(
                SegmentRecord(
                    index=progress.next_segment_index,
                    steps_completed=delta_steps,
                    steps_requested=max(delta_steps, scanned.total_steps_requested),
                    samples_written=0,
                    status=status,
                    duration_ns=(delta_steps * progress.timestep_fs) / 1e6,
                    finished_at=datetime.now(timezone.utc).isoformat(),
                )
            )
        elif not progress.segments and new_steps > 0 and scanned.segments:
            progress.segments.extend(scanned.segments)

        progress.status = scanned.status

    if mark_complete:
        progress.status = SimulationStatus.COMPLETED
        if progress.segments and progress.segments[-1].status != SegmentStatus.COMPLETED:
            progress.segments[-1].status = SegmentStatus.COMPLETED
            progress.segments[-1].finished_at = datetime.now(timezone.utc).isoformat()

    save_progress(working_dir, progress)
    return progress


def _parse_gromacs_log(log_path: Path) -> dict:
    """Parse ``prod.log`` for step counts and completion status.

    Parameters
    ----------
    log_path : Path
        Path to the GROMACS production log.

    Returns
    -------
    dict
        Dictionary containing:

        - ``steps_completed``: int
        - ``time_completed_ps``: float
        - ``is_finished``: bool
        - ``nsteps_requested``: int
    """
    if not log_path.exists():
        return {
            "steps_completed": 0,
            "time_completed_ps": 0.0,
            "is_finished": False,
            "nsteps_requested": 0,
        }

    text = log_path.read_text(errors="ignore")
    is_finished = "Finished mdrun" in text

    nsteps_match = re.search(r"\bnsteps\s*=\s*(\d+)", text)
    nsteps_requested = int(nsteps_match.group(1)) if nsteps_match else 0

    # Parse the step/time table lines and take the maximum observed step
    step_rows = re.findall(r"^\s*(\d+)\s+([0-9]+(?:\.[0-9]+)?)\s*$", text, flags=re.MULTILINE)
    steps_completed = 0
    time_completed_ps = 0.0
    for step_str, time_str in step_rows:
        step_val = int(step_str)
        if step_val >= steps_completed:
            steps_completed = step_val
            time_completed_ps = float(time_str)

    if steps_completed == 0:
        fallback_match = re.findall(r"Statistics over\s+(\d+)\s+steps", text)
        if fallback_match:
            steps_completed = max(int(item) for item in fallback_match)

    return {
        "steps_completed": steps_completed,
        "time_completed_ps": time_completed_ps,
        "is_finished": is_finished,
        "nsteps_requested": nsteps_requested,
    }


def _scan_equilibration_gromacs(working_dir: Path) -> list[EquilibrationStageRecord]:
    """Scan for completed GROMACS equilibration stages.

    Parameters
    ----------
    working_dir : Path
        GROMACS working directory.

    Returns
    -------
    list[EquilibrationStageRecord]
        Completed stage records sorted by index.
    """
    pattern = re.compile(r"^eq_(\d+)\.gro$")
    records: list[EquilibrationStageRecord] = []

    for path in sorted(working_dir.glob("eq_*.gro")):
        match = pattern.match(path.name)
        if match is None:
            continue
        idx = int(match.group(1))
        finished_at = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()
        records.append(
            EquilibrationStageRecord(
                index=idx - 1,
                name=f"eq_{idx:02d}",
                status=SegmentStatus.COMPLETED,
                finished_at=finished_at,
            )
        )

    records.sort(key=lambda item: item.index)
    return records


def load_or_scan_gromacs_progress(
    working_dir: Path,
    config_path: str = "",
    replicate: int = 1,
    total_steps: int = 0,
    total_samples: int = 0,
    timestep_fs: float = 2.0,
) -> SimulationProgress:
    """Load ``progress.json`` or scan filesystem as fallback.

    Parameters
    ----------
    working_dir : Path
        GROMACS working directory.
    config_path : str, optional
        Source config path.
    replicate : int, optional
        Replicate index.
    total_steps : int, optional
        Total production step count from config.
    total_samples : int, optional
        Total sample count from config.
    timestep_fs : float, optional
        Integration time step in femtoseconds.

    Returns
    -------
    SimulationProgress
        Current GROMACS progress model.
    """
    working_dir = Path(working_dir)
    progress = load_progress(working_dir)

    if progress is not None:
        scanned = scan_gromacs_progress(
            working_dir=working_dir,
            config_path=config_path or progress.config_path,
            replicate=replicate,
            total_steps=total_steps or progress.total_steps_requested,
            total_samples=total_samples or progress.total_samples_requested,
            timestep_fs=timestep_fs or progress.timestep_fs,
        )

        progress.equilibration_stages = scanned.equilibration_stages
        if scanned.total_steps_completed > progress.total_steps_completed:
            delta_steps = scanned.total_steps_completed - progress.total_steps_completed
            segment_status = (
                SegmentStatus.COMPLETED
                if scanned.status == SimulationStatus.COMPLETED
                else SegmentStatus.INTERRUPTED
            )
            progress.segments.append(
                SegmentRecord(
                    index=progress.next_segment_index,
                    steps_completed=delta_steps,
                    steps_requested=max(delta_steps, scanned.total_steps_requested),
                    samples_written=0,
                    status=segment_status,
                    duration_ns=(delta_steps * progress.timestep_fs) / 1e6,
                    finished_at=datetime.now(timezone.utc).isoformat(),
                )
            )

        if config_path:
            progress.config_path = config_path
        progress.replicate = replicate
        if total_steps > 0:
            progress.total_steps_requested = total_steps
        elif progress.total_steps_requested <= 0:
            progress.total_steps_requested = scanned.total_steps_requested
        if total_samples > 0:
            progress.total_samples_requested = total_samples
        elif progress.total_samples_requested <= 0:
            progress.total_samples_requested = scanned.total_samples_requested
        progress.timestep_fs = timestep_fs if timestep_fs > 0 else progress.timestep_fs
        progress.status = scanned.status
        return progress

    return scan_gromacs_progress(
        working_dir=working_dir,
        config_path=config_path,
        replicate=replicate,
        total_steps=total_steps,
        total_samples=total_samples,
        timestep_fs=timestep_fs,
    )
