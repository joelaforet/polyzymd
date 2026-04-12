"""Replicate-level SLURM orchestration for analysis comparisons.

This module provides a shared DAG submission layer for analysis plugins:

- one replicate worker per (condition, replicate)
- one aggregate worker per condition
- one finalizer worker per analysis comparison

The DAG parallelizes at the `compute_replicate()` boundary, with one SLURM job
per (condition, replicate) pair. This boundary is the plugin API's atomic unit.
Sub-replicate parallelism (for example, per-run work inside SASA-style
calculations) is intentionally handled inside each plugin's
`compute_replicate()` implementation. Plugins can use internal
threading/multiprocessing for that finer-grained work when needed.
"""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import stat
import subprocess
import time
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Literal, cast

from pydantic import BaseModel, Field, model_validator

from polyzymd.analyses.base import Analysis
from polyzymd.analyses.orchestrator import prepare_comparison_run
from polyzymd.analyses.shared.paths import sanitize_label
from polyzymd.config.comparison import ComparisonConfig

LOGGER = logging.getLogger(__name__)

_MEM_PATTERN = re.compile(r"^\d+(?:[KMGT]B?)$", flags=re.IGNORECASE)
_TIME_PATTERN = re.compile(r"^(?:\d+-)?\d{1,3}:\d{2}:\d{2}$")
_SLURM_JOB_ID_PATTERN = re.compile(r"^\d+(?:_\d+)?$")
_SBATCH_PARSE_FAILURE_MARKER = "SBATCH_PARSE_FAILURE"


class AnalysisSlurmResources(BaseModel):
    """SLURM resource settings for analysis workers."""

    pixi_path: str = Field(default="pixi")
    partition: str = "aa100"
    qos: str | None = None
    account: str | None = None
    ntasks: int = Field(default=1, ge=1, le=256)
    cpus_per_task: int = Field(default=1, ge=1, le=256)
    mem: str = "4G"
    time: str = "01:00:00"
    max_retries: int = Field(default=3, ge=1)
    mail_user: str | None = None
    mail_type: str = "FAIL"

    @model_validator(mode="after")
    def _validate_slurm_values(self) -> AnalysisSlurmResources:
        """Validate user-provided values used in SLURM scripts.

        Returns
        -------
        AnalysisSlurmResources
            Validated resource object.

        Raises
        ------
        ValueError
            If a field contains unsafe shell characters or invalid format.
        """
        self.pixi_path = _sanitize_slurm_value(self.pixi_path, "pixi_path")
        self.partition = _sanitize_slurm_value(self.partition, "partition")
        self.mem = _sanitize_slurm_value(self.mem, "mem")
        self.time = _sanitize_slurm_value(self.time, "time")
        self.mail_type = _sanitize_slurm_value(self.mail_type, "mail_type")
        if self.qos is not None:
            self.qos = _sanitize_slurm_value(self.qos, "qos")
        if self.account is not None:
            self.account = _sanitize_slurm_value(self.account, "account")
        if self.mail_user is not None:
            self.mail_user = _sanitize_slurm_value(self.mail_user, "mail_user")

        if not _MEM_PATTERN.match(self.mem):
            raise ValueError("Invalid mem format. Expected values like '4G', '8000M', or '16GB'.")
        if not _TIME_PATTERN.match(self.time):
            raise ValueError("Invalid time format. Expected 'HH:MM:SS' or 'D-HH:MM:SS'.")
        return self


class ReplicateTaskSpec(BaseModel):
    """Task spec for one replicate job."""

    condition_index: int
    replicate: int
    condition_label: str
    condition_slug: str


class ConditionTaskSpec(BaseModel):
    """Task spec for one condition aggregate job."""

    condition_index: int
    condition_label: str
    condition_slug: str
    replicate_specs: list[ReplicateTaskSpec]


class AnalysisJobManifest(BaseModel):
    """Snapshot of inputs needed to run analysis workers."""

    analysis_name: str
    comparison_yaml: str
    condition_specs: list[ConditionTaskSpec]
    settings_snapshot: dict[str, Any]
    snapshot_hash: str = ""
    partial_policy: Literal["strict", "allow_partial"] = "strict"
    equilibration: str
    recompute: bool = False
    resources: AnalysisSlurmResources
    created_at: str

    @model_validator(mode="after")
    def _populate_snapshot_hash(self) -> AnalysisJobManifest:
        """Backfill snapshot hash for backward compatibility.

        Returns
        -------
        AnalysisJobManifest
            Manifest with a snapshot hash.
        """
        if not self.snapshot_hash:
            self.snapshot_hash = compute_manifest_snapshot_hash(
                analysis_name=self.analysis_name,
                settings_snapshot=self.settings_snapshot,
                condition_specs=self.condition_specs,
                equilibration=self.equilibration,
            )
        return self

    def save(self, path: Path) -> Path:
        """Save manifest as JSON."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path) -> AnalysisJobManifest:
        """Load manifest from JSON."""
        return cls.model_validate_json(path.read_text())


class SubmittedJobGraph(BaseModel):
    """Submitted SLURM job IDs for analysis DAG nodes."""

    replicate_jobs: dict[tuple[int, int], str]
    array_jobs: dict[str, str] | None = None
    aggregator_jobs: dict[int, str]
    finalizer_job_id: str

    def save(self, path: Path) -> Path:
        """Save graph as JSON with portable keys."""
        payload = {
            "replicate_jobs": {
                f"{cond_idx}:{rep}": job_id
                for (cond_idx, rep), job_id in self.replicate_jobs.items()
            },
            "array_jobs": self.array_jobs,
            "aggregator_jobs": {
                str(cond_idx): job_id for cond_idx, job_id in self.aggregator_jobs.items()
            },
            "finalizer_job_id": self.finalizer_job_id,
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2))
        return path

    @classmethod
    def load(cls, path: Path) -> SubmittedJobGraph:
        """Load graph from JSON with tuple/int key reconstruction."""
        payload = json.loads(path.read_text())
        raw_replicate_jobs = {
            tuple(map(int, key.split(":"))): value
            for key, value in payload.get("replicate_jobs", {}).items()
        }
        replicate_jobs = cast(dict[tuple[int, int], str], raw_replicate_jobs)
        aggregator_jobs = {
            int(key): value for key, value in payload.get("aggregator_jobs", {}).items()
        }
        return cls(
            replicate_jobs=replicate_jobs,
            array_jobs=payload.get("array_jobs"),
            aggregator_jobs=aggregator_jobs,
            finalizer_job_id=payload["finalizer_job_id"],
        )


class TaskStatus(BaseModel):
    """Task status persisted by worker wrappers."""

    state: Literal["pending", "running", "succeeded", "failed", "retrying"]
    attempt_count: int = Field(default=0, ge=0)
    error_message: str | None = None
    last_updated: str
    slurm_job_id: str | None = None


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sanitize_slurm_value(value: str, field_name: str) -> str:
    """Validate a SLURM value before writing it into scripts.

    Parameters
    ----------
    value : str
        User-provided value.
    field_name : str
        Name of the field being validated.

    Returns
    -------
    str
        Sanitized value.

    Raises
    ------
    ValueError
        If value contains unsafe characters or is too long.
    """
    if len(value) > 128:
        raise ValueError(f"Unsafe SLURM value for '{field_name}': exceeds 128 characters")
    if any(ch.isspace() for ch in value):
        raise ValueError(f"Unsafe SLURM value for '{field_name}': whitespace is not allowed")
    if any(token in value for token in ("\n", "\r", ";", "`", "$(", "${", "|", "<", ">")):
        raise ValueError(
            f"Unsafe SLURM value for '{field_name}': contains disallowed shell characters"
        )
    return value


def _sanitize_path_for_script(path: Path) -> str:
    """Validate a filesystem path before script interpolation.

    Parameters
    ----------
    path : Path
        Filesystem path that will be interpolated into generated script content.

    Returns
    -------
    str
        Resolved path string when no unsafe tokens are present.

    Raises
    ------
    ValueError
        If the resolved path includes shell-unsafe or quote-breaking characters.
    """
    raw_path_str = str(path)
    if any(0 <= ord(char) <= 0x1F or ord(char) == 0x7F for char in raw_path_str):
        raise ValueError(
            "Unsafe path for generated scripts: found ASCII control character in "
            f"{raw_path_str!r}. Please rename the path to remove control characters"
        )

    path_str = str(path.resolve())
    if any(0 <= ord(char) <= 0x1F or ord(char) == 0x7F for char in path_str):
        raise ValueError(
            "Unsafe path for generated scripts: found ASCII control character in "
            f"{path_str!r}. Please rename the path to remove control characters"
        )

    disallowed_tokens = ("'", '"', "`", "$(", "${", "|", "<", ">", ";")
    for token in disallowed_tokens:
        if token in path_str:
            display_token = token.encode("unicode_escape").decode("ascii")
            raise ValueError(
                "Unsafe path for generated scripts: "
                f"found disallowed token '{display_token}' in {path_str!r}. "
                "Please rename the path to remove this character sequence"
            )
    return path_str


def _pixi_run_prefix(resources: AnalysisSlurmResources) -> str:
    """Return the pixi command prefix used in generated scripts."""
    pixi = resources.pixi_path
    if pixi == "pixi":
        detected = shutil.which("pixi")
        if detected is not None:
            pixi = detected
    return f'"{pixi}" run -e build'


def _condition_specs_from_conditions(conditions: list[Any]) -> list[ConditionTaskSpec]:
    """Build manifest condition task specs from prepared conditions."""
    condition_specs: list[ConditionTaskSpec] = []
    for cond_idx, condition in enumerate(conditions):
        slug = sanitize_label(condition.label)
        reps = [
            ReplicateTaskSpec(
                condition_index=cond_idx,
                replicate=rep,
                condition_label=condition.label,
                condition_slug=slug,
            )
            for rep in condition.replicates
        ]
        condition_specs.append(
            ConditionTaskSpec(
                condition_index=cond_idx,
                condition_label=condition.label,
                condition_slug=slug,
                replicate_specs=reps,
            )
        )
    return condition_specs


def compute_manifest_snapshot_hash(
    analysis_name: str,
    settings_snapshot: dict[str, Any],
    condition_specs: list[ConditionTaskSpec],
    equilibration: str,
) -> str:
    """Compute deterministic hash for manifest-sensitive comparison inputs."""
    payload = {
        "analysis_name": analysis_name,
        "settings_snapshot": settings_snapshot,
        "condition_specs": [spec.model_dump(mode="json") for spec in condition_specs],
        "equilibration": equilibration,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256(encoded).hexdigest()


def validate_manifest_snapshot(
    manifest: AnalysisJobManifest,
    analysis: Analysis,
    config: ComparisonConfig,
) -> tuple[list[Any], str, Path]:
    """Validate that live comparison inputs match the manifest snapshot.

    Returns
    -------
    tuple[list[Any], str, Path]
        Prepared conditions, resolved equilibration, and analysis root.

    Raises
    ------
    RuntimeError
        If current config/plugin settings drift from the submitted manifest.
    """
    prepared = prepare_comparison_run(
        analysis,
        config,
        manifest.equilibration,
    )
    valid_conditions = prepared["valid_conditions"]
    live_settings = prepared["settings"]
    resolved_equilibration = prepared["equilibration"]
    analysis_root = prepared["analysis_root"]
    live_settings_snapshot = (
        live_settings.model_dump(mode="json") if hasattr(live_settings, "model_dump") else {}
    )
    live_condition_specs = _condition_specs_from_conditions(valid_conditions)
    live_hash = compute_manifest_snapshot_hash(
        analysis_name=manifest.analysis_name,
        settings_snapshot=live_settings_snapshot,
        condition_specs=live_condition_specs,
        equilibration=resolved_equilibration,
    )
    if live_hash != manifest.snapshot_hash:
        raise RuntimeError(
            "Manifest/config drift detected. The current comparison.yaml or plugin settings no "
            "longer match the submitted manifest snapshot. Re-submit the analysis jobs."
        )
    return valid_conditions, resolved_equilibration, analysis_root


def _manifest_path(hpc_dir: Path) -> Path:
    return hpc_dir / "manifest.json"


def _graph_path(hpc_dir: Path) -> Path:
    return hpc_dir / "job_graph.json"


def _submission_error_path(hpc_dir: Path) -> Path:
    return hpc_dir / "submission_error.json"


def _scripts_dir(hpc_dir: Path) -> Path:
    return hpc_dir / "scripts"


def _logs_dir(hpc_dir: Path) -> Path:
    return hpc_dir / "logs"


def _rep_status_path(hpc_dir: Path, condition_slug: str, replicate: int) -> Path:
    return hpc_dir / "status" / "replicates" / condition_slug / f"rep_{replicate}.json"


def _cond_status_path(hpc_dir: Path, condition_slug: str) -> Path:
    return hpc_dir / "status" / "conditions" / f"{condition_slug}.json"


def _final_status_path(hpc_dir: Path) -> Path:
    return hpc_dir / "status" / "finalize.json"


def update_task_status(
    status_path: Path,
    state: Literal["pending", "running", "succeeded", "failed", "retrying"],
    attempt_count: int,
    error_message: str | None = None,
) -> None:
    """Atomically write a task status JSON file."""
    status_path.parent.mkdir(parents=True, exist_ok=True)
    payload = TaskStatus(
        state=state,
        attempt_count=attempt_count,
        error_message=error_message,
        last_updated=_utc_now(),
        slurm_job_id=os.getenv("SLURM_JOB_ID"),
    )
    with NamedTemporaryFile("w", encoding="utf-8", dir=status_path.parent, delete=False) as tmp:
        tmp.write(payload.model_dump_json(indent=2))
        tmp_path = Path(tmp.name)
    os.replace(tmp_path, status_path)


def _status_update_python(
    status_path: Path,
    state: str,
    resources: AnalysisSlurmResources,
    error_expr: str = "None",
) -> str:
    status_path_str = _sanitize_path_for_script(status_path)
    return (
        f'{_pixi_run_prefix(resources)} python -c "from pathlib import Path; '
        "from polyzymd.workflow.analysis_slurm import update_task_status; "
        f"update_task_status(Path(r'{status_path_str}'), '{state}', int($ATTEMPT), {error_expr})\""
    )


def _status_attempt_python(status_path: Path, resources: AnalysisSlurmResources) -> str:
    status_path_str = _sanitize_path_for_script(status_path)
    return (
        f'{_pixi_run_prefix(resources)} python -c "import json; '
        f"d=json.load(open(r'{status_path_str}')); "
        "print(d.get('attempt_count', 0))\""
    )


def _slurm_header(resources: AnalysisSlurmResources, job_name: str, log_path: Path) -> str:
    log_path_str = _sanitize_path_for_script(log_path)
    lines = [
        "#!/bin/bash",
        "#SBATCH --requeue",
        f"#SBATCH --partition={resources.partition}",
        f"#SBATCH --job-name={job_name}",
        f'#SBATCH --output="{log_path_str}"',
        f"#SBATCH --ntasks={resources.ntasks}",
        f"#SBATCH --cpus-per-task={resources.cpus_per_task}",
        f"#SBATCH --mem={resources.mem}",
        f"#SBATCH --time={resources.time}",
    ]
    if resources.qos:
        lines.append(f"#SBATCH --qos={resources.qos}")
    if resources.account:
        lines.append(f"#SBATCH --account={resources.account}")
    if resources.mail_user:
        lines.append(f"#SBATCH --mail-type={resources.mail_type}")
        lines.append(f"#SBATCH --mail-user={resources.mail_user}")
    return "\n".join(lines)


def _ensure_layout(hpc_dir: Path, manifest: AnalysisJobManifest) -> None:
    _scripts_dir(hpc_dir).mkdir(parents=True, exist_ok=True)
    _logs_dir(hpc_dir).mkdir(parents=True, exist_ok=True)
    (hpc_dir / "status" / "replicates").mkdir(parents=True, exist_ok=True)
    (hpc_dir / "status" / "conditions").mkdir(parents=True, exist_ok=True)

    for cond_spec in manifest.condition_specs:
        for rep_spec in cond_spec.replicate_specs:
            status_path = _rep_status_path(hpc_dir, cond_spec.condition_slug, rep_spec.replicate)
            if not status_path.exists():
                update_task_status(status_path, "pending", 0)
        cond_status = _cond_status_path(hpc_dir, cond_spec.condition_slug)
        if not cond_status.exists():
            update_task_status(cond_status, "pending", 0)

    fin = _final_status_path(hpc_dir)
    if not fin.exists():
        update_task_status(fin, "pending", 0)


def build_manifest(
    analysis: Analysis,
    config: ComparisonConfig,
    resources: AnalysisSlurmResources,
    recompute: bool,
    equilibration: str | None,
    allow_partial: bool = False,
) -> AnalysisJobManifest:
    """Build submission manifest from comparison config and plugin settings."""
    prepared = prepare_comparison_run(
        analysis,
        config,
        equilibration,
    )
    valid_conditions = prepared["valid_conditions"]
    settings = prepared["settings"]
    resolved_equilibration = prepared["equilibration"]
    condition_specs = _condition_specs_from_conditions(valid_conditions)
    settings_snapshot = settings.model_dump(mode="json") if hasattr(settings, "model_dump") else {}
    snapshot_hash = compute_manifest_snapshot_hash(
        analysis_name=analysis.name,
        settings_snapshot=settings_snapshot,
        condition_specs=condition_specs,
        equilibration=resolved_equilibration,
    )

    return AnalysisJobManifest(
        analysis_name=analysis.name,
        comparison_yaml=str(Path(config.source_path).resolve()) if config.source_path else "",
        condition_specs=condition_specs,
        settings_snapshot=settings_snapshot,
        snapshot_hash=snapshot_hash,
        partial_policy="allow_partial" if allow_partial else "strict",
        equilibration=resolved_equilibration,
        recompute=recompute,
        resources=resources,
        created_at=_utc_now(),
    )


def generate_replicate_script(
    manifest: AnalysisJobManifest,
    task_spec: ReplicateTaskSpec,
    resources: AnalysisSlurmResources,
    hpc_dir: Path,
) -> Path:
    """Generate a replicate worker script with automatic retries."""
    script_path = (
        _scripts_dir(hpc_dir) / f"replicate__{task_spec.condition_slug}__r{task_spec.replicate}.sh"
    )
    log_path = (
        _logs_dir(hpc_dir) / f"replicate__{task_spec.condition_slug}__r{task_spec.replicate}.%j.out"
    )
    status_path = _rep_status_path(hpc_dir, task_spec.condition_slug, task_spec.replicate)
    manifest_path = _manifest_path(hpc_dir)
    status_path_str = _sanitize_path_for_script(status_path)
    manifest_path_str = _sanitize_path_for_script(manifest_path)
    _sanitize_path_for_script(log_path)

    header = _slurm_header(
        resources,
        f"pzmd_r_{task_spec.condition_slug}_{task_spec.replicate}",
        log_path,
    )
    worker_cmd = (
        f"{_pixi_run_prefix(resources)} polyzymd compare worker-replicate "
        f'--manifest "{manifest_path_str}" '
        f"--condition-index {task_spec.condition_index} "
        f"--replicate {task_spec.replicate}"
    )
    error_expr = "'worker exit code ' + str($EXIT_CODE)"

    script = (
        f"{header}\n\n"
        "set -u\n\n"
        f'STATUS_FILE="{status_path_str}"\n'
        f'MANIFEST="{manifest_path_str}"\n'
        f"MAX_RETRIES={resources.max_retries}\n"
        f'ATTEMPT=$({_status_attempt_python(status_path, resources)} 2>/dev/null || echo "0")\n'
        "ATTEMPT=$((ATTEMPT + 1))\n"
        f"{_status_update_python(status_path, 'running', resources)}\n"
        f"{worker_cmd}\n"
        "EXIT_CODE=$?\n"
        "if [ $EXIT_CODE -ne 0 ]; then\n"
        "    if [ $ATTEMPT -lt $MAX_RETRIES ]; then\n"
        f"        {_status_update_python(status_path, 'retrying', resources, error_expr)}\n"
        "        scontrol requeue $SLURM_JOB_ID\n"
        "        exit 0\n"
        "    fi\n"
        f"    {_status_update_python(status_path, 'failed', resources, error_expr)}\n"
        "    exit $EXIT_CODE\n"
        "fi\n"
        f"{_status_update_python(status_path, 'succeeded', resources)}\n"
        "exit 0\n"
    )

    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text(script)
    script_path.chmod(script_path.stat().st_mode | stat.S_IXUSR)
    return script_path


def generate_aggregate_script(
    manifest: AnalysisJobManifest,
    cond_spec: ConditionTaskSpec,
    resources: AnalysisSlurmResources,
    hpc_dir: Path,
) -> Path:
    """Generate an aggregate worker script with retries."""
    script_path = _scripts_dir(hpc_dir) / f"aggregate__{cond_spec.condition_slug}.sh"
    log_path = _logs_dir(hpc_dir) / f"aggregate__{cond_spec.condition_slug}.%j.out"
    status_path = _cond_status_path(hpc_dir, cond_spec.condition_slug)
    manifest_path = _manifest_path(hpc_dir)
    status_path_str = _sanitize_path_for_script(status_path)
    manifest_path_str = _sanitize_path_for_script(manifest_path)
    _sanitize_path_for_script(log_path)

    header = _slurm_header(resources, f"pzmd_a_{cond_spec.condition_slug}", log_path)
    worker_cmd = (
        f"{_pixi_run_prefix(resources)} polyzymd compare worker-aggregate "
        f'--manifest "{manifest_path_str}" --condition-index {cond_spec.condition_index}'
    )
    error_expr = "'worker exit code ' + str($EXIT_CODE)"

    script = (
        f"{header}\n\n"
        "set -u\n\n"
        f'STATUS_FILE="{status_path_str}"\n'
        f"MAX_RETRIES={resources.max_retries}\n"
        f'ATTEMPT=$({_status_attempt_python(status_path, resources)} 2>/dev/null || echo "0")\n'
        "ATTEMPT=$((ATTEMPT + 1))\n"
        f"{_status_update_python(status_path, 'running', resources)}\n"
        f"{worker_cmd}\n"
        "EXIT_CODE=$?\n"
        "if [ $EXIT_CODE -ne 0 ]; then\n"
        "    if [ $ATTEMPT -lt $MAX_RETRIES ]; then\n"
        f"        {_status_update_python(status_path, 'retrying', resources, error_expr)}\n"
        "        scontrol requeue $SLURM_JOB_ID\n"
        "        exit 0\n"
        "    fi\n"
        f"    {_status_update_python(status_path, 'failed', resources, error_expr)}\n"
        "    exit $EXIT_CODE\n"
        "fi\n"
        f"{_status_update_python(status_path, 'succeeded', resources)}\n"
        "exit 0\n"
    )

    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text(script)
    script_path.chmod(script_path.stat().st_mode | stat.S_IXUSR)
    return script_path


def _array_spec(replicates: list[int]) -> str:
    """Build a SLURM array selector specification.

    Parameters
    ----------
    replicates : list[int]
        Replicate IDs to include in one array job.

    Returns
    -------
    str
        Comma-separated array spec such as ``"1,2,3"``.

    Raises
    ------
    ValueError
        If no replicate IDs are provided.
    """
    ordered = sorted(set(replicates))
    if not ordered:
        raise ValueError("At least one replicate is required for an array submission")
    return ",".join(str(rep) for rep in ordered)


def generate_array_script(
    cond_spec: ConditionTaskSpec,
    manifest: AnalysisJobManifest,
    resources: AnalysisSlurmResources,
    replicates: list[int],
    hpc_dir: Path,
) -> Path:
    """Generate one array worker script for all replicates of a condition.

    Parameters
    ----------
    cond_spec : ConditionTaskSpec
        Condition task specification from the manifest.
    manifest : AnalysisJobManifest
        Submission manifest used by worker commands.
    resources : AnalysisSlurmResources
        SLURM resource settings used in script header.
    replicates : list[int]
        Replicate IDs included in this array job.
    hpc_dir : Path
        Root directory where scripts and logs are written.

    Returns
    -------
    Path
        Generated executable script path.
    """
    script_path = _scripts_dir(hpc_dir) / f"replicate_array__{cond_spec.condition_slug}.sh"
    log_path = _logs_dir(hpc_dir) / f"replicate_array__{cond_spec.condition_slug}.%a.log"
    manifest_path = _manifest_path(hpc_dir)
    manifest_path_str = _sanitize_path_for_script(manifest_path)
    _sanitize_path_for_script(log_path)

    header = _slurm_header(resources, f"pzmd_ra_{cond_spec.condition_slug}", log_path)
    array_spec = _array_spec(replicates)

    case_branches: list[str] = []
    for replicate in sorted(set(replicates)):
        status_path = _rep_status_path(hpc_dir, cond_spec.condition_slug, replicate)
        error_expr = "'worker exit code ' + str($EXIT_CODE)"
        worker_cmd = (
            f"{_pixi_run_prefix(resources)} polyzymd compare worker-replicate "
            f'--manifest "{manifest_path_str}" '
            f"--condition-index {cond_spec.condition_index} "
            "--replicate $REP"
        )
        branch = (
            f"{replicate})\n"
            f'    ATTEMPT=$({_status_attempt_python(status_path, resources)} 2>/dev/null || echo "0")\n'
            "    ATTEMPT=$((ATTEMPT + 1))\n"
            f"    {_status_update_python(status_path, 'running', resources)}\n"
            f"    {worker_cmd}\n"
            "    EXIT_CODE=$?\n"
            "    if [ $EXIT_CODE -ne 0 ]; then\n"
            "        if [ $ATTEMPT -lt $MAX_RETRIES ]; then\n"
            f"            {_status_update_python(status_path, 'retrying', resources, error_expr)}\n"
            "            scontrol requeue ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}\n"
            "            exit 0\n"
            "        fi\n"
            f"        {_status_update_python(status_path, 'failed', resources, error_expr)}\n"
            "        exit $EXIT_CODE\n"
            "    fi\n"
            f"    {_status_update_python(status_path, 'succeeded', resources)}\n"
            "    ;;\n"
        )
        case_branches.append(branch)

    script = (
        f"{header}\n"
        f"#SBATCH --array={array_spec}\n\n"
        "set -u\n\n"
        f"MAX_RETRIES={resources.max_retries}\n"
        "REP=$SLURM_ARRAY_TASK_ID\n\n"
        'case "$REP" in\n'
        + "".join(case_branches)
        + "*)\n"
        + '    echo "Unknown array task id: $REP" >&2\n'
        + "    exit 2\n"
        + "    ;;\n"
        + "esac\n\n"
        + "exit 0\n"
    )

    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text(script)
    script_path.chmod(script_path.stat().st_mode | stat.S_IXUSR)
    return script_path


def generate_finalize_script(
    manifest: AnalysisJobManifest,
    resources: AnalysisSlurmResources,
    hpc_dir: Path,
) -> Path:
    """Generate the final compare+plot worker script with retries."""
    script_path = _scripts_dir(hpc_dir) / "finalize.sh"
    log_path = _logs_dir(hpc_dir) / "finalize.%j.out"
    status_path = _final_status_path(hpc_dir)
    manifest_path = _manifest_path(hpc_dir)
    status_path_str = _sanitize_path_for_script(status_path)
    manifest_path_str = _sanitize_path_for_script(manifest_path)
    _sanitize_path_for_script(log_path)

    header = _slurm_header(resources, f"pzmd_f_{manifest.analysis_name}", log_path)
    worker_cmd = (
        f"{_pixi_run_prefix(resources)} polyzymd compare worker-finalize "
        f'--manifest "{manifest_path_str}"'
    )
    error_expr = "'worker exit code ' + str($EXIT_CODE)"

    script = (
        f"{header}\n\n"
        "set -u\n\n"
        f'STATUS_FILE="{status_path_str}"\n'
        f"MAX_RETRIES={resources.max_retries}\n"
        f'ATTEMPT=$({_status_attempt_python(status_path, resources)} 2>/dev/null || echo "0")\n'
        "ATTEMPT=$((ATTEMPT + 1))\n"
        f"{_status_update_python(status_path, 'running', resources)}\n"
        f"{worker_cmd}\n"
        "EXIT_CODE=$?\n"
        "if [ $EXIT_CODE -ne 0 ]; then\n"
        "    if [ $ATTEMPT -lt $MAX_RETRIES ]; then\n"
        f"        {_status_update_python(status_path, 'retrying', resources, error_expr)}\n"
        "        scontrol requeue $SLURM_JOB_ID\n"
        "        exit 0\n"
        "    fi\n"
        f"    {_status_update_python(status_path, 'failed', resources, error_expr)}\n"
        "    exit $EXIT_CODE\n"
        "fi\n"
        f"{_status_update_python(status_path, 'succeeded', resources)}\n"
        "exit 0\n"
    )

    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text(script)
    script_path.chmod(script_path.stat().st_mode | stat.S_IXUSR)
    return script_path


def _submit_sbatch(script_path: Path, dependency: str | None = None) -> str:
    cmd = ["sbatch"]
    if dependency:
        cmd.extend(["--dependency", dependency])
    cmd.append(str(script_path))
    try:
        completed = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "SLURM is not available: 'sbatch' not found on PATH. The HPC submission "
            "commands require a SLURM cluster. Run analysis locally with "
            "'polyzymd compare run' instead."
        ) from exc
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        details = stderr if stderr else "No stderr output from sbatch"
        raise RuntimeError(f"SLURM submission failed for {script_path}: {details}") from exc
    output = completed.stdout.strip()
    match = re.search(r"Submitted batch job\s+(\d+)", output)
    if match is None:
        raise RuntimeError(
            f"{_SBATCH_PARSE_FAILURE_MARKER}: Could not parse job id from sbatch output. "
            f"Raw stdout: {output!r}"
        )
    return match.group(1)


def _cancel_jobs(job_ids: list[str]) -> dict[str, dict[str, Any]]:
    """Attempt to cancel submitted SLURM jobs.

    Parameters
    ----------
    job_ids : list[str]
        Job IDs to cancel.

    Returns
    -------
    dict[str, dict[str, Any]]
        Per-job cancellation outcomes keyed by job ID with retry metadata.

    Notes
    -----
    Cancellation is best-effort. Failures are logged and suppressed.
    """
    results: dict[str, dict[str, Any]] = {}
    if not job_ids:
        return results

    for job_id in job_ids:
        if not _SLURM_JOB_ID_PATTERN.fullmatch(job_id):
            LOGGER.warning("Skipping invalid SLURM job id for cancellation: %r", job_id)
            results[job_id] = {
                "attempted": False,
                "cancelled": False,
                "attempts": 0,
                "error": "invalid_job_id",
            }
            continue

        attempts = 0
        cancelled = False
        error: str | None = None

        for attempt in (1, 2):
            attempts = attempt
            LOGGER.info("Cancelling SLURM job %s (attempt %d/2)", job_id, attempt)
            try:
                subprocess.run(["scancel", job_id], capture_output=True, text=True, check=True)
            except FileNotFoundError as exc:
                error = f"scancel_not_found: {exc}"
                LOGGER.warning("Could not cancel submitted jobs: 'scancel' not found (%s)", exc)
                break
            except subprocess.CalledProcessError as exc:
                stderr = (exc.stderr or "").strip()
                error = stderr if stderr else "No stderr output from scancel"
                LOGGER.warning(
                    "Failed to cancel submitted job %s on attempt %d/2: %s",
                    job_id,
                    attempt,
                    error,
                )
                if attempt == 1:
                    time.sleep(2)
                    continue
            else:
                cancelled = True
                error = None
                break

            break

        results[job_id] = {
            "attempted": True,
            "cancelled": cancelled,
            "attempts": attempts,
            "error": error,
        }

    return results


def _query_sacct(job_ids: list[str]) -> dict[str, str]:
    """Query sacct for base SLURM job states.

    Parameters
    ----------
    job_ids : list[str]
        SLURM job IDs to query.

    Returns
    -------
    dict[str, str]
        Mapping of base job ID to SLURM state string.

    Notes
    -----
    This parser ignores sub-step rows such as ``12345.batch`` and only keeps
    top-level job records.

    If ``sacct`` emits duplicate rows for the same top-level job ID, conflict
    resolution is deterministic:

    - terminal states are preferred over non-terminal states
    - among terminal states, the most recent terminal row is preferred
    """
    if not job_ids:
        return {}

    try:
        completed = subprocess.run(
            [
                "sacct",
                "-j",
                ",".join(job_ids),
                "--format=JobIDRaw,State",
                "--noheader",
                "--parsable2",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError as exc:
        LOGGER.warning("Could not reconcile status with sacct: command not found (%s)", exc)
        return {}
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        details = stderr if stderr else "No stderr output from sacct"
        LOGGER.warning("Could not reconcile status with sacct: %s", details)
        return {}

    states: dict[str, str] = {}
    for raw_line in completed.stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        parts = line.split("|")
        if len(parts) != 2:
            continue

        job_id, state = parts[0].strip(), parts[1].strip()
        if not job_id or not state:
            continue
        if "." in job_id:
            continue

        existing = states.get(job_id)
        if existing is None:
            states[job_id] = state
            continue

        existing_terminal = _is_terminal_slurm_state(existing)
        current_terminal = _is_terminal_slurm_state(state)

        if current_terminal and not existing_terminal:
            states[job_id] = state
            continue
        if current_terminal and existing_terminal:
            states[job_id] = state
            continue
        if not current_terminal and not existing_terminal:
            states[job_id] = state

    return states


def _is_terminal_slurm_state(slurm_state: str) -> bool:
    """Return whether a SLURM state string is terminal.

    Parameters
    ----------
    slurm_state : str
        Raw SLURM state string from ``sacct``.

    Returns
    -------
    bool
        ``True`` when the state is terminal, otherwise ``False``.
    """
    normalized = slurm_state.strip().upper().split(maxsplit=1)[0]
    if normalized.startswith("CANCELLED"):
        return True
    return normalized in {
        "COMPLETED",
        "FAILED",
        "TIMEOUT",
        "OUT_OF_MEMORY",
        "NODE_FAIL",
        "BOOT_FAIL",
        "DEADLINE",
        "REVOKED",
        "PREEMPTED",
    }


def _map_slurm_state(slurm_state: str) -> str | None:
    """Map SLURM states to local task states.

    Parameters
    ----------
    slurm_state : str
        Raw SLURM state as returned by ``sacct``.

    Returns
    -------
    str | None
        New local state if a terminal transition is needed, otherwise ``None``.
    """
    normalized = slurm_state.strip().upper()
    normalized = normalized.split(maxsplit=1)[0]

    if normalized == "COMPLETED":
        return "succeeded"
    if normalized in {
        "FAILED",
        "OUT_OF_MEMORY",
        "TIMEOUT",
        "NODE_FAIL",
        "BOOT_FAIL",
        "DEADLINE",
        "REVOKED",
    }:
        return "failed"
    if normalized.startswith("CANCELLED") or normalized == "PREEMPTED":
        return "failed"
    if normalized in {"RUNNING", "PENDING", "SUSPENDED", "REQUEUED"}:
        return None

    LOGGER.debug("Unknown SLURM state %r during reconciliation", slurm_state)
    return None


def reconcile_status_with_slurm(hpc_dir: Path) -> dict[str, Any]:
    """Reconcile local status files with live SLURM accounting state.

    Parameters
    ----------
    hpc_dir : Path
        Root HPC artifact directory for one analysis submission.

    Returns
    -------
    dict[str, Any]
        Summary with checked file count, update count, and per-file changes.
    """
    rep_root = hpc_dir / "status" / "replicates"
    cond_root = hpc_dir / "status" / "conditions"

    status_paths: list[Path] = []
    if rep_root.exists():
        for cond_dir in sorted(p for p in rep_root.iterdir() if p.is_dir()):
            status_paths.extend(sorted(cond_dir.glob("rep_*.json")))
    if cond_root.exists():
        status_paths.extend(sorted(cond_root.glob("*.json")))
    finalize_path = _final_status_path(hpc_dir)
    if finalize_path.exists():
        status_paths.append(finalize_path)

    actionable: list[tuple[Path, dict[str, Any], str]] = []
    for status_path in status_paths:
        try:
            payload = json.loads(status_path.read_text())
            if not isinstance(payload, dict):
                raise ValueError("Status payload must be a JSON object")

            local_state_raw = payload["state"]
            job_id_raw = payload["slurm_job_id"]
            if not isinstance(local_state_raw, str):
                raise ValueError("Status payload field 'state' must be a string")
            if job_id_raw is not None and not isinstance(job_id_raw, str):
                raise ValueError("Status payload field 'slurm_job_id' must be a string or null")
        except (
            json.JSONDecodeError,
            KeyError,
            ValueError,
            FileNotFoundError,
            PermissionError,
        ) as exc:
            LOGGER.warning(
                "Skipping unreadable status file during reconciliation: %s (%s)", status_path, exc
            )
            continue

        local_state = local_state_raw.strip().lower()
        job_id = job_id_raw
        if (
            local_state in {"running", "pending", "retrying"}
            and isinstance(job_id, str)
            and job_id.strip()
        ):
            actionable.append((status_path, payload, job_id.strip()))

    if not actionable:
        return {"checked": 0, "updated": 0, "changes": []}

    unique_job_ids = sorted({job_id for _, _, job_id in actionable})
    sacct_states = _query_sacct(unique_job_ids)

    changes: list[dict[str, str]] = []
    for status_path, payload, job_id in actionable:
        slurm_state = sacct_states.get(job_id)
        if slurm_state is None:
            continue

        mapped_state = _map_slurm_state(slurm_state)
        if mapped_state is None:
            continue

        old_state = str(payload.get("state", "")).strip().lower()
        if mapped_state == old_state:
            continue

        try:
            latest_payload = json.loads(status_path.read_text())
            if not isinstance(latest_payload, dict):
                raise ValueError("Status payload must be a JSON object")
            latest_state_raw = latest_payload["state"]
            if not isinstance(latest_state_raw, str):
                raise ValueError("Status payload field 'state' must be a string")
            latest_state = latest_state_raw.strip().lower()
        except (
            json.JSONDecodeError,
            KeyError,
            ValueError,
            FileNotFoundError,
            PermissionError,
        ) as exc:
            LOGGER.warning(
                "Skipping unreadable status file during reconciliation: %s (%s)", status_path, exc
            )
            continue

        if latest_state != old_state:
            LOGGER.debug(
                "Skipping status reconciliation for %s because state changed from %s to %s",
                status_path,
                old_state,
                latest_state,
            )
            continue

        latest_job_id = latest_payload.get("slurm_job_id")
        if latest_job_id != job_id:
            LOGGER.debug(
                "Skipping status reconciliation for %s because slurm_job_id changed from %s to %s",
                status_path,
                job_id,
                latest_job_id,
            )
            continue

        latest_payload["state"] = mapped_state
        latest_payload["reconciled_from"] = slurm_state
        latest_payload["reconciled_at"] = _utc_now()

        if mapped_state == "failed" and (
            slurm_state.upper().startswith("CANCELLED")
            or slurm_state.upper().startswith("PREEMPTED")
        ):
            latest_payload["error_message"] = f"Reconciled from SLURM state {slurm_state}"

        try:
            with NamedTemporaryFile(
                "w", encoding="utf-8", dir=status_path.parent, delete=False
            ) as tmp:
                tmp.write(json.dumps(latest_payload, indent=2))
                tmp_path = Path(tmp.name)
            os.replace(tmp_path, status_path)
        except OSError as write_exc:
            LOGGER.warning("Failed to write reconciled status for %s: %s", status_path, write_exc)
            if tmp_path.exists():
                tmp_path.unlink(missing_ok=True)
            continue

        changes.append(
            {
                "job_id": job_id,
                "path": str(status_path),
                "old_state": old_state,
                "new_state": mapped_state,
                "slurm_state": slurm_state,
            }
        )

    return {
        "checked": len(actionable),
        "updated": len(changes),
        "changes": changes,
    }


def submit_analysis_graph(
    manifest: AnalysisJobManifest,
    resources: AnalysisSlurmResources,
    hpc_dir: Path,
) -> SubmittedJobGraph:
    """Submit replicate, aggregate, and finalizer jobs with dependencies.

    Parameters
    ----------
    manifest : AnalysisJobManifest
        Submission manifest describing all condition and replicate tasks.
    resources : AnalysisSlurmResources
        SLURM resource settings used for generated scripts.
    hpc_dir : Path
        Root directory where scripts, logs, and submission metadata are stored.

    Returns
    -------
    SubmittedJobGraph
        Graph of submitted replicate, aggregate, and finalizer job IDs.

    Raises
    ------
    RuntimeError
        Propagated if any ``sbatch`` submission fails.
    """
    _ensure_layout(hpc_dir, manifest)
    submission_error_path = _submission_error_path(hpc_dir)
    if submission_error_path.exists():
        submission_error_path.unlink()
    manifest.save(_manifest_path(hpc_dir))

    replicate_jobs: dict[tuple[int, int], str] = {}
    aggregator_jobs: dict[int, str] = {}
    submitted_job_ids: list[str] = []

    try:
        for cond_spec in manifest.condition_specs:
            for rep_spec in cond_spec.replicate_specs:
                script = generate_replicate_script(manifest, rep_spec, resources, hpc_dir)
                job_id = _submit_sbatch(script)
                replicate_jobs[(rep_spec.condition_index, rep_spec.replicate)] = job_id
                submitted_job_ids.append(job_id)

        for cond_spec in manifest.condition_specs:
            replicate_ids = [
                replicate_jobs[(cond_spec.condition_index, rep_spec.replicate)]
                for rep_spec in cond_spec.replicate_specs
            ]
            dependency = "afterany:" + ":".join(replicate_ids)
            script = generate_aggregate_script(manifest, cond_spec, resources, hpc_dir)
            job_id = _submit_sbatch(script, dependency=dependency)
            aggregator_jobs[cond_spec.condition_index] = job_id
            submitted_job_ids.append(job_id)

        aggregate_dependency = "afterany:" + ":".join(
            aggregator_jobs[idx] for idx in sorted(aggregator_jobs.keys())
        )
        finalize_script = generate_finalize_script(manifest, resources, hpc_dir)
        finalizer_job_id = _submit_sbatch(finalize_script, dependency=aggregate_dependency)
        submitted_job_ids.append(finalizer_job_id)
    except (RuntimeError, subprocess.SubprocessError, OSError) as exc:
        cancel_results: dict[str, dict[str, Any]] = {}
        if submitted_job_ids:
            cancel_results = _cancel_jobs(submitted_job_ids)
        cancelled_count = sum(
            1 for result in cancel_results.values() if result.get("cancelled") is True
        )
        LOGGER.warning(
            "Analysis DAG submission failed and cancellation succeeded for %d/%d tracked jobs",
            cancelled_count,
            len(submitted_job_ids),
        )
        if _SBATCH_PARSE_FAILURE_MARKER in str(exc):
            LOGGER.warning(
                "Submission failed while parsing sbatch output; one or more jobs may not have "
                "been tracked for rollback"
            )

        raw_sbatch_stdout: str | None = None
        if _SBATCH_PARSE_FAILURE_MARKER in str(exc):
            marker = "Raw stdout:"
            if marker in str(exc):
                raw_sbatch_stdout = str(exc).split(marker, maxsplit=1)[1].strip()

        error_payload = {
            "error": str(exc),
            "cancelled_job_ids": submitted_job_ids,
            "cancellation_results": cancel_results,
            "raw_sbatch_stdout": raw_sbatch_stdout,
            "timestamp": _utc_now(),
        }
        try:
            submission_error_path.write_text(json.dumps(error_payload, indent=2))
        except (OSError, TypeError, ValueError) as sidecar_exc:
            LOGGER.warning(
                "Failed to write submission error sidecar at %s: %s",
                submission_error_path,
                sidecar_exc,
            )
        raise

    graph = SubmittedJobGraph(
        replicate_jobs=replicate_jobs,
        aggregator_jobs=aggregator_jobs,
        finalizer_job_id=finalizer_job_id,
    )
    graph.save(_graph_path(hpc_dir))
    return graph


def submit_analysis_graph_with_arrays(
    manifest: AnalysisJobManifest,
    resources: AnalysisSlurmResources,
    hpc_dir: Path,
) -> SubmittedJobGraph:
    """Submit one array per condition plus aggregate and finalizer jobs.

    Parameters
    ----------
    manifest : AnalysisJobManifest
        Submission manifest describing condition and replicate tasks.
    resources : AnalysisSlurmResources
        SLURM resource settings used for generated scripts.
    hpc_dir : Path
        Root directory where scripts, logs, and submission metadata are stored.

    Returns
    -------
    SubmittedJobGraph
        Graph of submitted array, aggregate, and finalizer job IDs.

    Raises
    ------
    RuntimeError
        Propagated if any ``sbatch`` submission fails.
    """
    _ensure_layout(hpc_dir, manifest)
    submission_error_path = _submission_error_path(hpc_dir)
    if submission_error_path.exists():
        submission_error_path.unlink()
    manifest.save(_manifest_path(hpc_dir))

    array_jobs: dict[str, str] = {}
    aggregator_jobs: dict[int, str] = {}
    submitted_job_ids: list[str] = []

    try:
        for cond_spec in manifest.condition_specs:
            replicates = [rep_spec.replicate for rep_spec in cond_spec.replicate_specs]
            script = generate_array_script(cond_spec, manifest, resources, replicates, hpc_dir)
            job_id = _submit_sbatch(script)
            array_jobs[cond_spec.condition_slug] = job_id
            submitted_job_ids.append(job_id)

        for cond_spec in manifest.condition_specs:
            parent_array_job_id = array_jobs[cond_spec.condition_slug]
            dependency = f"afterany:{parent_array_job_id}"
            script = generate_aggregate_script(manifest, cond_spec, resources, hpc_dir)
            job_id = _submit_sbatch(script, dependency=dependency)
            aggregator_jobs[cond_spec.condition_index] = job_id
            submitted_job_ids.append(job_id)

        aggregate_dependency = "afterany:" + ":".join(
            aggregator_jobs[idx] for idx in sorted(aggregator_jobs.keys())
        )
        finalize_script = generate_finalize_script(manifest, resources, hpc_dir)
        finalizer_job_id = _submit_sbatch(finalize_script, dependency=aggregate_dependency)
        submitted_job_ids.append(finalizer_job_id)
    except (RuntimeError, subprocess.SubprocessError, OSError) as exc:
        cancel_results: dict[str, dict[str, Any]] = {}
        if submitted_job_ids:
            cancel_results = _cancel_jobs(submitted_job_ids)
        cancelled_count = sum(
            1 for result in cancel_results.values() if result.get("cancelled") is True
        )
        LOGGER.warning(
            "Analysis DAG array submission failed and cancellation succeeded for %d/%d tracked jobs",
            cancelled_count,
            len(submitted_job_ids),
        )
        if _SBATCH_PARSE_FAILURE_MARKER in str(exc):
            LOGGER.warning(
                "Submission failed while parsing sbatch output; one or more jobs may not have "
                "been tracked for rollback"
            )

        raw_sbatch_stdout: str | None = None
        if _SBATCH_PARSE_FAILURE_MARKER in str(exc):
            marker = "Raw stdout:"
            if marker in str(exc):
                raw_sbatch_stdout = str(exc).split(marker, maxsplit=1)[1].strip()

        error_payload = {
            "error": str(exc),
            "cancelled_job_ids": submitted_job_ids,
            "cancellation_results": cancel_results,
            "raw_sbatch_stdout": raw_sbatch_stdout,
            "timestamp": _utc_now(),
        }
        try:
            submission_error_path.write_text(json.dumps(error_payload, indent=2))
        except (OSError, TypeError, ValueError) as sidecar_exc:
            LOGGER.warning(
                "Failed to write submission error sidecar at %s: %s",
                submission_error_path,
                sidecar_exc,
            )
        raise

    graph = SubmittedJobGraph(
        replicate_jobs={},
        array_jobs=array_jobs,
        aggregator_jobs=aggregator_jobs,
        finalizer_job_id=finalizer_job_id,
    )
    graph.save(_graph_path(hpc_dir))
    return graph


def read_analysis_status(hpc_dir: Path) -> dict[str, Any]:
    """Read all status files for one analysis HPC run."""
    status_root = hpc_dir / "status"
    rep_root = status_root / "replicates"
    cond_root = status_root / "conditions"

    replicate_summary: dict[str, dict[str, dict[str, Any]]] = {}
    warnings: list[str] = []
    if rep_root.exists():
        for cond_dir in sorted(p for p in rep_root.iterdir() if p.is_dir()):
            per_rep: dict[str, dict[str, Any]] = {}
            for status_file in sorted(cond_dir.glob("rep_*.json")):
                try:
                    status = TaskStatus.model_validate_json(status_file.read_text())
                    per_rep[status_file.stem] = status.model_dump()
                except (
                    json.JSONDecodeError,
                    KeyError,
                    ValueError,
                    FileNotFoundError,
                    PermissionError,
                ) as exc:
                    warning = f"Corrupted status file: {status_file} ({type(exc).__name__}: {exc})"
                    LOGGER.warning(warning)
                    warnings.append(warning)
                    per_rep[status_file.stem] = {
                        "state": "unknown",
                        "attempt_count": 0,
                        "error_message": warning,
                        "last_updated": None,
                        "slurm_job_id": None,
                    }
            replicate_summary[cond_dir.name] = per_rep

    condition_summary: dict[str, dict[str, Any]] = {}
    if cond_root.exists():
        for status_file in sorted(cond_root.glob("*.json")):
            try:
                status = TaskStatus.model_validate_json(status_file.read_text())
                condition_summary[status_file.stem] = status.model_dump()
            except (
                json.JSONDecodeError,
                KeyError,
                ValueError,
                FileNotFoundError,
                PermissionError,
            ) as exc:
                warning = f"Corrupted status file: {status_file} ({type(exc).__name__}: {exc})"
                LOGGER.warning(warning)
                warnings.append(warning)
                condition_summary[status_file.stem] = {
                    "state": "unknown",
                    "attempt_count": 0,
                    "error_message": warning,
                    "last_updated": None,
                    "slurm_job_id": None,
                }

    finalize = None
    finalize_path = _final_status_path(hpc_dir)
    if finalize_path.exists():
        try:
            finalize = TaskStatus.model_validate_json(finalize_path.read_text()).model_dump()
        except (
            json.JSONDecodeError,
            KeyError,
            ValueError,
            FileNotFoundError,
            PermissionError,
        ) as exc:
            warning = f"Corrupted status file: {finalize_path} ({type(exc).__name__}: {exc})"
            LOGGER.warning(warning)
            warnings.append(warning)
            finalize = {
                "state": "unknown",
                "attempt_count": 0,
                "error_message": warning,
                "last_updated": None,
                "slurm_job_id": None,
            }

    states: list[str] = []
    for cond_map in replicate_summary.values():
        states.extend(item["state"] for item in cond_map.values())
    states.extend(item["state"] for item in condition_summary.values())
    if finalize is not None:
        states.append(finalize["state"])

    counts = {
        state: states.count(state)
        for state in [
            "pending",
            "running",
            "retrying",
            "succeeded",
            "failed",
            "unknown",
        ]
    }

    return {
        "replicates": replicate_summary,
        "conditions": condition_summary,
        "finalize": finalize,
        "counts": counts,
        "warnings": warnings,
    }
