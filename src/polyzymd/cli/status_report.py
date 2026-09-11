"""Agent-oriented status reporting for ``polyzymd status``.

The original ``polyzymd status`` printed a progress bar per replicate. That
answers "how far along is it?" but not the two questions an operator (or an
LLM agent acting for one) actually asks about a preemptable restart chain:

* **Is anything still driving this replicate?** A ``progress.json`` status of
  ``interrupted`` is ambiguous: the chain may be mid-resubmit, preempted and
  requeued, or dead after a ``FATAL`` in the wrapper. Only SLURM knows.
* **When will it finish?** The segment records already carry wall-clock
  timestamps and step counts, so a throughput and an ETA cost nothing.

This module joins the three sources that answer those questions —
``progress.json`` (via the engine's ``load_or_scan_progress``), a single
``squeue`` call for the user's jobs, and the newest SLURM log per replicate —
into one :class:`SystemReport` per config, and renders it either as a compact
plain-text block designed to be cheap in tokens or as JSON.

Everything that touches the filesystem or SLURM is injectable so the logic is
unit-testable without a scheduler.
"""

from __future__ import annotations

import getpass
import json
import logging
import os
import re
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

LOGGER = logging.getLogger(__name__)

# Minimum wall time a segment must have run before its throughput is trusted.
_MIN_RATE_WINDOW = timedelta(minutes=10)
# Minimum steps a segment must have advanced before its throughput is trusted.
_MIN_RATE_STEPS = 1_000

# Lines in a SLURM log that explain why a chain stopped. Ordered so that the
# most specific wrapper verdicts win over generic Python tracebacks.
_ERROR_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^FATAL:.*"),
    re.compile(r"^CONCURRENT:.*"),
    re.compile(r"^STOP:.*"),
    re.compile(r"^Segment \d+ failed:.*"),
    re.compile(r"^Validation error:.*"),
    re.compile(r"^slurmstepd: error:.*"),
    re.compile(r"^\w*Error: .*"),
)
_LOG_TAIL_LINES = 80
_ERROR_MAX_CHARS = 160


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class SlurmJob:
    """One queued or running SLURM job, as reported by ``squeue``."""

    job_id: str
    name: str
    state: str
    elapsed: str = ""
    node: str = ""
    reason: str = ""


@dataclass
class ReplicateReport:
    """Everything the report knows about one replicate."""

    replicate: int
    directory: str | None
    progress_status: str
    completed_ns: float
    total_ns: float
    fraction: float
    verdict: str
    jobs: list[SlurmJob] = field(default_factory=list)
    rate_ns_per_day: float | None = None
    eta_days: float | None = None
    last_error: str | None = None
    last_log: str | None = None
    note: str | None = None

    def to_dict(self) -> dict:
        data = asdict(self)
        return data


@dataclass
class SystemReport:
    """Report for one config file (one system, N replicates)."""

    name: str
    config_path: str
    scratch_directory: str
    replicates: list[ReplicateReport]
    error: str | None = None

    def counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for rep in self.replicates:
            out[rep.verdict] = out.get(rep.verdict, 0) + 1
        return out

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "config_path": self.config_path,
            "scratch_directory": self.scratch_directory,
            "error": self.error,
            "counts": self.counts(),
            "replicates": [r.to_dict() for r in self.replicates],
        }


# Verdict vocabulary. Kept deliberately small so an agent can branch on it.
VERDICT_COMPLETED = "completed"
VERDICT_RUNNING = "running"
VERDICT_QUEUED = "queued"
VERDICT_DEAD = "dead"
VERDICT_NOT_STARTED = "not_started"
VERDICT_NOT_FOUND = "not_found"
VERDICT_ORDER = (
    VERDICT_COMPLETED,
    VERDICT_RUNNING,
    VERDICT_QUEUED,
    VERDICT_DEAD,
    VERDICT_NOT_STARTED,
    VERDICT_NOT_FOUND,
)


# ---------------------------------------------------------------------------
# SLURM
# ---------------------------------------------------------------------------


def query_user_jobs(user: str | None = None, timeout: int = 20) -> list[SlurmJob] | None:
    """Return every queued or running job for *user* in one ``squeue`` call.

    Returns ``None`` (not an empty list) when ``squeue`` is unavailable or
    fails, so callers can distinguish "no jobs" from "could not ask".
    """
    if user is None:
        user = os.environ.get("USER") or getpass.getuser()
    cmd = [
        "squeue",
        "--noheader",
        "-u",
        user,
        "--states",
        "RUNNING,PENDING,COMPLETING,CONFIGURING,SUSPENDED,REQUEUED",
        "--format",
        "%i|%j|%t|%M|%N|%R",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError) as exc:
        LOGGER.warning(f"squeue unavailable ({exc}); SLURM state will be omitted")
        return None
    if result.returncode != 0:
        LOGGER.warning(
            f"squeue exited {result.returncode}; SLURM state will be omitted: "
            f"{result.stderr.strip()[:200]}"
        )
        return None
    return parse_squeue_output(result.stdout)


def parse_squeue_output(text: str) -> list[SlurmJob]:
    """Parse ``squeue --format "%i|%j|%t|%M|%N|%R"`` output."""
    jobs: list[SlurmJob] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        parts = line.split("|")
        if len(parts) < 3:
            continue
        parts += [""] * (6 - len(parts))
        job_id, name, state, elapsed, node, reason = parts[:6]
        jobs.append(
            SlurmJob(
                job_id=job_id.strip(),
                name=name.strip(),
                state=state.strip(),
                elapsed=elapsed.strip(),
                node=node.strip(),
                reason=reason.strip(),
            )
        )
    return jobs


def jobs_by_name(jobs: Iterable[SlurmJob]) -> dict[str, list[SlurmJob]]:
    out: dict[str, list[SlurmJob]] = {}
    for job in jobs:
        out.setdefault(job.name, []).append(job)
    return out


# ---------------------------------------------------------------------------
# Throughput / ETA
# ---------------------------------------------------------------------------


def _parse_iso(value: str | None) -> datetime | None:
    if not value:
        return None
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed


def estimate_rate_ns_per_day(progress, now: datetime, *, live: bool) -> float | None:
    """Estimate throughput from the most recent informative segment.

    Preference order:

    1. The newest segment with a ``finished_at`` timestamp and enough steps
       (a completed or cleanly interrupted segment has a known wall window).
    2. If *live* is true (SLURM reports a running job), the newest segment
       without ``finished_at``, timed from ``started_at`` to *now*.

    Returns ``None`` when no segment gives a trustworthy window.
    """
    segments = sorted(getattr(progress, "segments", []) or [], key=lambda s: s.index)
    timestep_fs = float(getattr(progress, "timestep_fs", 0.0) or 0.0)
    if timestep_fs <= 0:
        return None

    def _rate(steps: int, start: datetime, end: datetime) -> float | None:
        window = end - start
        if steps < _MIN_RATE_STEPS or window < _MIN_RATE_WINDOW:
            return None
        ns = steps * timestep_fs / 1e6
        return ns / (window.total_seconds() / 86400.0)

    for seg in reversed(segments):
        start = _parse_iso(seg.started_at)
        end = _parse_iso(seg.finished_at)
        if start is None or end is None:
            continue
        rate = _rate(seg.steps_completed, start, end)
        if rate is not None:
            return rate

    if live:
        for seg in reversed(segments):
            if seg.finished_at:
                continue
            start = _parse_iso(seg.started_at)
            if start is None:
                continue
            rate = _rate(seg.steps_completed, start, now)
            if rate is not None:
                return rate
    return None


def estimate_eta_days(remaining_ns: float, rate_ns_per_day: float | None) -> float | None:
    if rate_ns_per_day is None or rate_ns_per_day <= 0:
        return None
    if remaining_ns <= 0:
        return 0.0
    return remaining_ns / rate_ns_per_day


# ---------------------------------------------------------------------------
# SLURM log inspection
# ---------------------------------------------------------------------------


def find_latest_log(logs_dir: Path | None, job_name: str) -> Path | None:
    """Newest ``<job_name>.<jobid>.out`` in *logs_dir* by modification time."""
    if logs_dir is None or not logs_dir.is_dir():
        return None
    candidates = list(logs_dir.glob(f"{job_name}.*.out"))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def _tail_lines(path: Path, n: int) -> list[str]:
    try:
        with path.open("rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            block = min(size, 64 * 1024)
            fh.seek(size - block)
            data = fh.read().decode("utf-8", errors="replace")
    except OSError:
        return []
    return data.splitlines()[-n:]


def last_error_line(log_path: Path | None) -> str | None:
    """Return the most informative failure line near the end of a SLURM log."""
    if log_path is None:
        return None
    lines = _tail_lines(log_path, _LOG_TAIL_LINES)
    best: tuple[int, str] | None = None
    for line in lines:
        stripped = line.strip()
        for rank, pattern in enumerate(_ERROR_PATTERNS):
            if pattern.match(stripped):
                # Later lines override earlier ones at the same or better rank.
                if best is None or rank <= best[0]:
                    best = (rank, stripped)
                break
    if best is None:
        return None
    text = best[1]
    if len(text) > _ERROR_MAX_CHARS:
        text = text[: _ERROR_MAX_CHARS - 1] + "…"
    return text


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------


def classify(progress_status: str, jobs: Sequence[SlurmJob], has_dir: bool) -> str:
    """Map (progress.json status, live SLURM jobs) to a verdict."""
    if not has_dir:
        return VERDICT_NOT_FOUND
    if progress_status == "completed":
        return VERDICT_COMPLETED
    if any(j.state in ("R", "CG", "CF", "S") for j in jobs):
        return VERDICT_RUNNING
    if jobs:
        return VERDICT_QUEUED
    if progress_status == "not_started":
        return VERDICT_NOT_STARTED
    return VERDICT_DEAD


def build_system_report(
    sim_config,
    config_path: str | Path,
    *,
    engine_inst,
    jobs_index: Mapping[str, list[SlurmJob]] | None,
    now: datetime | None = None,
    job_name_fn: Callable[[object, int], str] | None = None,
    save_progress_fn: Callable[[Path, object], object] | None = None,
) -> SystemReport:
    """Assemble a :class:`SystemReport` for one loaded config.

    Parameters
    ----------
    sim_config
        A loaded ``SimulationConfig``.
    config_path
        Path the config was loaded from (echoed in the report so an agent can
        copy it into a follow-up command).
    engine_inst
        Engine created for *sim_config*; supplies working-directory resolution
        and ``load_or_scan_progress``.
    jobs_index
        ``{job_name: [SlurmJob, ...]}`` for the user's live jobs, or ``None``
        when SLURM could not be queried.
    """
    from polyzymd.simulation.progress import SimulationStatus

    if now is None:
        now = datetime.now(timezone.utc)
    if job_name_fn is None:
        from polyzymd.workflow.daisy_chain import create_job_name

        job_name_fn = create_job_name

    prod = sim_config.simulation_phases.production
    total_ns = float(prod.duration)
    dir_name = sim_config._format_run_directory_name(1)
    system_name = dir_name.rsplit("_run", 1)[0] if "_run" in dir_name else dir_name
    scratch = Path(sim_config.output.effective_scratch_directory)

    try:
        logs_dir = Path(sim_config.output.get_slurm_logs_directory())
    except Exception:  # pragma: no cover - defensive against mocked configs
        logs_dir = None

    replicates: list[ReplicateReport] = []
    for rep_num, rep_path in sorted(dict(sim_config.discover_replicate_dirs()).items()):
        job_name = job_name_fn(sim_config, rep_num)
        jobs = list(jobs_index.get(job_name, [])) if jobs_index is not None else []

        if rep_path is None:
            replicates.append(
                ReplicateReport(
                    replicate=rep_num,
                    directory=None,
                    progress_status="not_found",
                    completed_ns=0.0,
                    total_ns=total_ns,
                    fraction=0.0,
                    verdict=VERDICT_NOT_FOUND,
                    jobs=jobs,
                )
            )
            continue

        engine_dir = engine_inst.resolve_engine_working_directory(rep_path)
        progress = engine_inst.load_or_scan_progress(engine_dir, rep_num)
        if save_progress_fn is not None:
            save_progress_fn(engine_dir, progress)

        status_val = progress.status
        status_str = status_val.value if hasattr(status_val, "value") else str(status_val)
        steps_for_display = progress.total_steps_completed
        if status_val == SimulationStatus.FAILED:
            steps_for_display = max((seg.steps_completed for seg in progress.segments), default=0)
        fraction = (
            min(1.0, steps_for_display / progress.total_steps_requested)
            if progress.total_steps_requested
            else 0.0
        )
        completed_ns = steps_for_display * progress.timestep_fs / 1e6

        verdict = classify(status_str, jobs, has_dir=True)
        live = verdict == VERDICT_RUNNING
        rate = None
        eta = None
        if verdict != VERDICT_COMPLETED:
            rate = estimate_rate_ns_per_day(progress, now, live=live)
            if verdict in (VERDICT_RUNNING, VERDICT_QUEUED):
                eta = estimate_eta_days(total_ns - completed_ns, rate)

        last_err = None
        last_log = None
        if verdict in (VERDICT_DEAD, VERDICT_NOT_STARTED):
            log_path = find_latest_log(logs_dir, job_name)
            if log_path is None and logs_dir is not None and logs_dir.is_dir():
                # Build logs are named build_r<N>_<jobid>.out
                builds = list(logs_dir.glob(f"build_r{rep_num}_*.out"))
                if builds:
                    log_path = max(builds, key=lambda p: p.stat().st_mtime)
            if log_path is not None:
                last_log = log_path.name
                last_err = last_error_line(log_path)

        note = None
        if jobs_index is None:
            note = "slurm unavailable"

        replicates.append(
            ReplicateReport(
                replicate=rep_num,
                directory=str(rep_path),
                progress_status=status_str,
                completed_ns=completed_ns,
                total_ns=total_ns,
                fraction=fraction,
                verdict=verdict,
                jobs=jobs,
                rate_ns_per_day=rate,
                eta_days=eta,
                last_error=last_err,
                last_log=last_log,
                note=note,
            )
        )

    return SystemReport(
        name=system_name,
        config_path=str(config_path),
        scratch_directory=str(scratch),
        replicates=replicates,
    )


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _fmt_days(days: float | None) -> str:
    if days is None:
        return "eta ?"
    if days < 1.0:
        return f"eta {days * 24:.0f}h"
    return f"eta {days:.1f}d"


def _fmt_job(job: SlurmJob) -> str:
    parts = [f"job {job.job_id}", job.state]
    if job.state == "R":
        if job.elapsed:
            parts.append(job.elapsed)
        if job.node:
            parts.append(job.node)
    elif job.reason:
        parts.append(f"({job.reason})")
    return " ".join(parts)


def render_replicate_line(rep: ReplicateReport, label_width: int) -> str:
    label = f"run{rep.replicate}"
    pct = f"{rep.fraction * 100:3.0f}%"
    ns = f"{rep.completed_ns:6.1f}/{rep.total_ns:.0f}ns"
    verdict = rep.verdict.upper()
    fields = [f"{label:<{label_width}}", ns, pct, f"{verdict:<11}"]

    if rep.verdict == VERDICT_COMPLETED:
        pass
    elif rep.verdict in (VERDICT_RUNNING, VERDICT_QUEUED):
        fields.extend(_fmt_job(j) for j in rep.jobs)
        if rep.rate_ns_per_day is not None:
            fields.append(f"{rep.rate_ns_per_day:.0f}ns/d")
        fields.append(_fmt_days(rep.eta_days))
    elif rep.verdict == VERDICT_NOT_FOUND:
        fields.append("no directory in scratch")
    else:
        if rep.note:
            fields.append(rep.note)
        else:
            fields.append("no job")
        if rep.last_error:
            src = f" [{rep.last_log}]" if rep.last_log else ""
            fields.append(f"last: {rep.last_error}{src}")
        elif rep.last_log:
            fields.append(f"last log {rep.last_log} (no error line found)")
        else:
            fields.append("no slurm log found")
    return "  ".join(fields)


def render_agent(
    reports: Sequence[SystemReport],
    *,
    now: datetime | None = None,
    slurm_available: bool = True,
    preset_hint: str | None = None,
) -> str:
    """Compact, colour-free, fixed-vocabulary text for agents and terminals."""
    if now is None:
        now = datetime.now(timezone.utc)

    totals: dict[str, int] = {}
    n_reps = 0
    for rep_report in reports:
        for verdict, count in rep_report.counts().items():
            totals[verdict] = totals.get(verdict, 0) + count
            n_reps += count

    lines: list[str] = []
    summary = ", ".join(f"{totals[v]} {v}" for v in VERDICT_ORDER if totals.get(v))
    header = (
        f"# polyzymd status  {now.strftime('%Y-%m-%d %H:%M UTC')}  "
        f"{len(reports)} system(s)  {n_reps} replicate(s): {summary or 'none'}"
    )
    lines.append(header)
    if not slurm_available:
        lines.append(
            "# WARNING: squeue unavailable — DEAD/QUEUED cannot be distinguished from RUNNING; "
            "verdicts below fall back to progress.json only"
        )

    for report in reports:
        lines.append("")
        lines.append(f"## {report.name}  ({report.config_path})")
        if report.error:
            lines.append(f"ERROR {report.error}")
            continue
        if not report.replicates:
            lines.append(f"no replicate directories in {report.scratch_directory}")
            continue
        width = max(len(f"run{r.replicate}") for r in report.replicates)
        for rep in report.replicates:
            lines.append(render_replicate_line(rep, width))

    dead = [(r, rep) for r in reports for rep in r.replicates if rep.verdict == VERDICT_DEAD]
    if dead:
        lines.append("")
        lines.append("# dead chains — resume from checkpoint with:")
        by_cfg: dict[str, list[int]] = {}
        for report, rep in dead:
            by_cfg.setdefault(report.config_path, []).append(rep.replicate)
        preset = f" --preset {preset_hint}" if preset_hint else " --preset <preset>"
        for cfg, reps in by_cfg.items():
            rep_arg = ",".join(str(r) for r in sorted(reps))
            lines.append(f"polyzymd submit -c {cfg} -r {rep_arg}{preset}")
    return "\n".join(lines) + "\n"


def render_json(
    reports: Sequence[SystemReport], *, now: datetime | None = None, slurm_available: bool = True
) -> str:
    if now is None:
        now = datetime.now(timezone.utc)
    payload = {
        "generated_at": now.isoformat(),
        "slurm_available": slurm_available,
        "systems": [r.to_dict() for r in reports],
    }
    return json.dumps(payload, indent=1, default=str) + "\n"


# ---------------------------------------------------------------------------
# Config discovery
# ---------------------------------------------------------------------------


def discover_config_files(roots: Iterable[str | Path], *, max_depth: int = 3) -> list[Path]:
    """Find ``config.yaml`` files under *roots*, at most *max_depth* levels deep.

    Depth is bounded so that pointing ``--all`` at a large project tree does
    not walk into trajectory directories on slow network filesystems.
    """
    found: list[Path] = []
    for root in roots:
        root_path = Path(root)
        if root_path.is_file():
            found.append(root_path)
            continue
        if not root_path.is_dir():
            continue
        base_depth = len(root_path.resolve().parts)
        for dirpath, dirnames, filenames in os.walk(root_path):
            depth = len(Path(dirpath).resolve().parts) - base_depth
            if depth >= max_depth:
                dirnames[:] = []
            dirnames[:] = [d for d in dirnames if not d.startswith(".")]
            if "config.yaml" in filenames:
                found.append(Path(dirpath) / "config.yaml")
    return sorted(set(found))
