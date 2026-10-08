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
    _mtime_iso,
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
            _segment_record(
                working_dir,
                index=0,
                steps_completed=steps_completed,
                steps_requested=requested_steps,
                status=segment_status,
                duration_ns=duration_ns,
                first_start=True,
                start_step=0,
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
    since: str = "",
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
    since : str, optional
        ISO time when this job started. Records without a version that
        started at or after it ran in this job, also when ``polyzymd status``
        saved them first.

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
        ran = [*scanned.equilibration_stages, *scanned.segments]
    else:
        progress = existing
        earlier = {record.index: record for record in existing.equilibration_stages}
        ran = [
            record
            for record in scanned.equilibration_stages
            if record.index not in earlier
            or earlier[record.index].finished_at != record.finished_at
        ]
        progress.equilibration_stages = _keep_provenance(
            scanned.equilibration_stages, existing.equilibration_stages
        )
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
        old_segments = len(progress.segments)

        if delta_steps > 0:
            status = (
                SegmentStatus.COMPLETED
                if scanned.status == SimulationStatus.COMPLETED
                else SegmentStatus.INTERRUPTED
            )
            progress.segments.append(
                _segment_record(
                    working_dir,
                    index=progress.next_segment_index,
                    steps_completed=delta_steps,
                    steps_requested=max(delta_steps, scanned.total_steps_requested),
                    status=status,
                    duration_ns=(delta_steps * progress.timestep_fs) / 1e6,
                    start_step=old_steps,
                )
            )
        elif not progress.segments and new_steps > 0 and scanned.segments:
            progress.segments.extend(scanned.segments)
        ran.extend(progress.segments[old_segments:])

        progress.status = scanned.status

    if mark_complete:
        progress.status = SimulationStatus.COMPLETED
        if progress.segments and progress.segments[-1].status != SegmentStatus.COMPLETED:
            progress.segments[-1].status = SegmentStatus.COMPLETED
            progress.segments[-1].finished_at = (
                progress.segments[-1].finished_at or datetime.now(timezone.utc).isoformat()
            )

    if since:
        ran.extend(
            record
            for record in started_since(
                [*progress.equilibration_stages, *progress.segments],
                datetime.fromisoformat(since),
            )
            if record.polyzymd_version is None
        )
    record_run_provenance(ran)
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
    try:
        text = log_path.read_text(errors="ignore")
    except OSError:
        return {
            "steps_completed": 0,
            "time_completed_ps": 0.0,
            "is_finished": False,
            "nsteps_requested": 0,
        }
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
        log = working_dir / f"eq_{idx:02d}.log"
        mdp = next(iter(sorted(working_dir.glob(f"eq_{idx:02d}_*.mdp"))), None)
        started_at, finished_at = _mdrun_times(log, first_start=True)
        mdp_values = _mdp_values(mdp)
        steps = int(_parse_gromacs_log(log)["steps_completed"])
        record = EquilibrationStageRecord(
            index=idx - 1,
            name=f"eq_{idx:02d}",
            status=SegmentStatus.COMPLETED,
            duration_ns=steps * _dt_ps(mdp_values) / 1000.0,
            ensemble="NVT" if mdp_values.get("pcoupl", "no").lower() == "no" else "NPT",
            finished_at=finished_at or _mtime_iso(path),
            seeds=_seeds(mdp, log),
        )
        if started_at is not None:
            record.started_at = started_at
        records.append(record)

    records.sort(key=lambda item: item.index)
    return records


#: ``Started mdrun on rank 0 Wed Oct  7 00:16:48 2026`` and its ``Finished`` twin.
_MDRUN_TIME = re.compile(r"^(Started|Finished) mdrun on rank 0 (.+?)\s*$", re.MULTILINE)


def _mdrun_times(log_path: Path, first_start: bool = False) -> tuple[str | None, str | None]:
    """Return when mdrun started and finished, from the lines it writes to ``log_path``.

    A restarted run appends to its log, so the start is the last one, or the
    first with ``first_start``. The finish is the one after the last start,
    so a run that is still going or was killed has none.
    Times are the local time of the reading machine, written as UTC ISO
    timestamps; ``None`` when the log has no such line.
    """
    try:
        text = log_path.read_text(errors="ignore")
    except OSError:
        return None, None
    started: str | None = None
    finished: str | None = None
    for kind, stamp in _MDRUN_TIME.findall(text):
        try:
            when = datetime.strptime(" ".join(stamp.split()), "%a %b %d %H:%M:%S %Y")
        except ValueError:
            continue
        iso = when.astimezone(timezone.utc).isoformat()
        if kind == "Started":
            if started is None or not first_start:
                started = iso
            finished = None
        else:
            finished = iso
    return started, finished


def _mdp_values(mdp_path: Path | None) -> dict[str, str]:
    """Return the ``key = value`` settings of an MDP file, keys with ``-`` written as ``_``."""
    values: dict[str, str] = {}
    try:
        text = mdp_path.read_text(errors="ignore") if mdp_path is not None else ""
    except OSError:
        return values
    for line in text.splitlines():
        key, sep, value = line.split(";", 1)[0].partition("=")
        if sep:
            values[key.strip().lower().replace("-", "_")] = value.strip()
    return values


def _dt_ps(mdp: dict[str, str]) -> float:
    """Return the MDP time step in ps; GROMACS uses 0.001 when ``dt`` is not set."""
    try:
        return float(mdp.get("dt", "0.001"))
    except ValueError:
        return 0.001


def _keep_provenance(
    scanned: list[EquilibrationStageRecord], earlier: list[EquilibrationStageRecord]
) -> list[EquilibrationStageRecord]:
    """Give the rescanned stage records the provenance that ``earlier`` recorded for them."""
    by_index = {record.index: record for record in earlier}
    for record in scanned:
        if record.index in by_index:
            record.polyzymd_version = by_index[record.index].polyzymd_version
            record.pixi_environment = by_index[record.index].pixi_environment
    return scanned


def started_since(
    records: list[EquilibrationStageRecord | SegmentRecord], since: datetime
) -> list[EquilibrationStageRecord | SegmentRecord]:
    """Return the records that started at or after ``since``."""
    return [record for record in records if datetime.fromisoformat(record.started_at) >= since]


def record_run_provenance(records: list[EquilibrationStageRecord | SegmentRecord]) -> None:
    """Record this process's PolyzyMD version and pixi environment on ``records``.

    Pass only the records of stages and segments that this process ran, so a
    scan never records a version that did not run them. ``openmm_version``
    stays None.
    """
    from polyzymd.utils.version import record_provenance

    found = record_provenance()
    for record in records:
        record.polyzymd_version = found["polyzymd_version"]
        record.pixi_environment = found["pixi_environment"]


def _seeds(mdp_path: Path | None, log_path: Path) -> dict[str, int | None] | None:
    """Return the ``ld_seed`` mdrun used and the ``gen_seed`` of a stage that drew velocities.

    ``ld_seed`` comes from the parameters mdrun writes to the log, which hold
    the seed GROMACS chose when the MDP said -1; ``gen_seed`` comes from the
    MDP file when it sets ``gen_vel = yes``. A seed of -1 (GROMACS chose one
    and did not log it) is recorded as None.
    """
    mdp = _mdp_values(mdp_path)
    seeds: dict[str, int] = {}
    try:
        logged = re.findall(r"^\s*ld-seed\s*=\s*(-?\d+)", log_path.read_text(errors="ignore"), re.M)
    except OSError:
        logged = []
    if logged:
        seeds["ld_seed"] = int(logged[-1])
    elif mdp.get("ld_seed", "").lstrip("-").isdigit():
        seeds["ld_seed"] = int(mdp["ld_seed"])
    if mdp.get("gen_vel", "no").lower() == "yes" and mdp.get("gen_seed", "").lstrip("-").isdigit():
        seeds["gen_seed"] = int(mdp["gen_seed"])
    return {name: (None if seed == -1 else seed) for name, seed in seeds.items()} or None


def _segment_record(
    working_dir: Path, start_step: int, first_start: bool = False, **fields
) -> SegmentRecord:
    """Return a production ``SegmentRecord`` with the times and seeds of ``prod.log``.

    The segment runs ``fields["steps_completed"]`` steps after ``start_step``.
    ``samples_written`` counts the ``prod.xtc`` frames of those steps, from
    ``nstxout-compressed`` in ``prod.mdp``; mdrun also writes step 0.
    ``first_start`` takes the start of the first mdrun of the log, for a
    record that covers every run of it.
    """
    log = working_dir / "prod.log"
    mdp = working_dir / "prod.mdp"
    started_at, finished_at = _mdrun_times(log, first_start=first_start)
    interval = _mdp_values(mdp).get("nstxout_compressed", "0")
    interval = int(interval) if interval.isdigit() else 0
    end_step = start_step + fields["steps_completed"]
    samples = 0
    if interval > 0:
        samples = end_step // interval - start_step // interval + int(start_step == 0)
    # A run that was killed, or is still going, ends at its last log write.
    record = SegmentRecord(
        samples_written=samples,
        finished_at=finished_at or _mtime_iso(log),
        seeds=_seeds(mdp, log),
        **fields,
    )
    if started_at is not None:
        record.started_at = started_at
    return record


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

        progress.equilibration_stages = _keep_provenance(
            scanned.equilibration_stages, progress.equilibration_stages
        )
        if scanned.total_steps_completed > progress.total_steps_completed:
            delta_steps = scanned.total_steps_completed - progress.total_steps_completed
            segment_status = (
                SegmentStatus.COMPLETED
                if scanned.status == SimulationStatus.COMPLETED
                else SegmentStatus.INTERRUPTED
            )
            progress.segments.append(
                _segment_record(
                    working_dir,
                    index=progress.next_segment_index,
                    steps_completed=delta_steps,
                    steps_requested=max(delta_steps, scanned.total_steps_requested),
                    status=segment_status,
                    duration_ns=(delta_steps * progress.timestep_fs) / 1e6,
                    start_step=progress.total_steps_completed,
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
