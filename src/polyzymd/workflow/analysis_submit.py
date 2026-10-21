"""Submit ``polyzymd analyze`` to SLURM as one job per replicate and a report job.

:func:`write_submission` writes, for one ``polyzymd analyze`` command, a
SLURM array with one task per condition and replicate, each measuring and
storing its replicate's result, and a report job that runs the full command
once they have all ended, reusing every stored result for the statistics
and figures. :func:`submit` submits them with ``sbatch``. A task that fails
leaves its replicate unmeasured, and the report job measures it itself.
"""

from __future__ import annotations

import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Sequence

#: SLURM settings of each analysis preset. Analysis runs on CPUs, so these are
#: CPU partitions; a value of ``None`` leaves the directive out. Add a preset
#: here for another cluster.
ANALYSIS_PRESETS: dict[str, dict[str, str | None]] = {
    "alpine-cpu": {"partition": "acpu", "qos": "cpu-normal", "account": None},
    "blanca-shirts": {
        "partition": "blanca-shirts",
        "qos": "blanca-shirts",
        "account": "blanca-shirts",
    },
    "blanca-chbe-rdi": {
        "partition": "blanca-chbe-rdi",
        "qos": "blanca-chbe-rdi",
        "account": "blanca-chbe-rdi",
    },
    "bridges2-rm": {"partition": "RM-shared", "qos": None, "account": None},
}

_SAFE = re.compile(r"^[A-Za-z0-9_.:,/+=@-]+$")


@dataclass
class Resources:
    """SLURM settings of the submitted jobs."""

    partition: str | None = None
    qos: str | None = None
    account: str | None = None
    time: str = "12:00:00"
    mem: str = "16G"
    cpus: int = 2

    @classmethod
    def from_preset(cls, preset: str | None, **overrides: object) -> Resources:
        """Return the settings of ``preset``, with every override that is not ``None``."""
        from polyzymd.analyses.exceptions import ProtocolError

        values: dict[str, object] = {}
        if preset is not None:
            if preset not in ANALYSIS_PRESETS:
                raise ProtocolError(
                    f"No analysis preset named {preset!r}.",
                    hint=f"Use one of {', '.join(sorted(ANALYSIS_PRESETS))}, or give "
                    "--partition, --account and --qos.",
                )
            values.update(ANALYSIS_PRESETS[preset])
        values.update({key: value for key, value in overrides.items() if value is not None})
        resources = cls(**values)  # type: ignore[arg-type]
        for name in ("partition", "qos", "account", "time", "mem"):
            value = getattr(resources, name)
            if value is not None and not _SAFE.match(str(value)):
                raise ProtocolError(
                    f"The SLURM setting {name}={value!r} has characters a batch script cannot "
                    "carry safely.",
                    hint="Use letters, digits and . : , / + = @ - _ only.",
                )
        return resources

    def directives(self) -> list[str]:
        """Return the ``#SBATCH`` lines of these settings."""
        lines = [
            f"#SBATCH --{name}={value}"
            for name, value in (
                ("partition", self.partition),
                ("qos", self.qos),
                ("account", self.account),
            )
            if value
        ]
        return [
            *lines,
            "#SBATCH --ntasks=1",
            f"#SBATCH --cpus-per-task={int(self.cpus)}",
            f"#SBATCH --mem={self.mem}",
            f"#SBATCH --time={self.time}",
        ]


@dataclass
class Submission:
    """The files of one submission, and the job IDs once submitted."""

    folder: Path
    tasks: Path
    array: Path
    report: Path
    n_tasks: int
    array_id: str | None = None
    report_id: str | None = None
    report_output: Path | None = None
    notes: list[str] = field(default_factory=list)


def polyzymd_command() -> tuple[str, list[str]]:
    """Return the command that runs this PolyzyMD, and the environment lines the jobs need.

    The jobs run the Python interpreter and ``PYTHONPATH`` of the process
    that submits them, so they measure with the same PolyzyMD and packages,
    whether that is a pixi environment, a development install or a copy of
    the source.
    """
    import os
    import sys

    lines = ["export MPLBACKEND=Agg"]
    if os.environ.get("PYTHONPATH"):
        lines.insert(0, f"export PYTHONPATH={shlex.quote(os.environ['PYTHONPATH'])}")
    return f"{shlex.quote(sys.executable)} -m polyzymd.cli.main", lines


def write_submission(
    name: str,
    tasks: Sequence[tuple[Path, str, int]],
    task_options: Sequence[str],
    report_arguments: Sequence[str],
    resources: Resources,
    output_dir: Path,
    command: str,
    working_dir: Path,
    report_output: Path | None = None,
    json_report: bool = True,
    environment: Sequence[str] = (),
    study_file: Path | None = None,
) -> Submission:
    """Write the task list, the array script and the report script of one submission.

    ``tasks`` holds the config, label and replicate of every array task.
    Each task runs ``<command> analyze <name> -c <config> --label <label>
    --replicates <replicate> <task_options>``, and the report job runs
    ``<command> analyze <name> <report_arguments>``; both run in
    ``working_dir``, so relative paths resolve as they did when submitting.
    The files go to ``<output_dir>/slurm/<name>_<time>/``, with the job logs
    under its ``logs/``; the report job writes its report to
    ``report_output``, by default ``report.json`` (or ``report.txt``) there.
    With ``study_file``, each task runs ``<command> analyze <name> --study
    <study_file> --label <label> --replicates <replicate> <task_options>``
    instead, so ``name`` is a run of that study file.
    """
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    folder = Path(output_dir).resolve() / "slurm" / f"{name}_{stamp}"
    (folder / "logs").mkdir(parents=True, exist_ok=True)
    task_file = folder / "tasks.tsv"
    task_file.write_text(
        "".join(f"{config}\t{label}\t{replicate}\n" for config, label, replicate in tasks)
    )
    header = ["#!/bin/bash", *resources.directives(), f'#SBATCH --chdir="{working_dir}"']
    setup = ["set -euo pipefail", *environment]
    options = " ".join(shlex.quote(option) for option in task_options)
    array = folder / "replicates.sbatch"
    array.write_text(
        "\n".join(
            [
                *header,
                f"#SBATCH --job-name={name}-replicates",
                f"#SBATCH --array=0-{len(tasks) - 1}",
                f'#SBATCH --output="{folder}/logs/replicate.%A_%a.out"',
                "# One condition and replicate per task: measure it and store its result.",
                *setup,
                f"IFS=$'\\t' read -r config label replicate < <(sed -n \"$((SLURM_ARRAY_TASK_ID + 1))p\" {shlex.quote(str(task_file))})",
                f"{command} analyze {shlex.quote(name)} "
                + (
                    f"--study {shlex.quote(str(study_file))} "
                    if study_file is not None
                    else '-c "$config" '
                )
                + f'--label "$label" --replicates "$replicate" {options}',
                "",
            ]
        )
    )
    output = (
        Path(report_output).resolve()
        if report_output
        else folder / ("report.json" if json_report else "report.txt")
    )
    report_arguments = [*report_arguments, "-o", str(output)]
    report = folder / "report.sbatch"
    report.write_text(
        "\n".join(
            [
                *header,
                f"#SBATCH --job-name={name}-report",
                f'#SBATCH --output="{folder}/logs/report.%j.out"',
                "# Every condition: reads the stored results, measures any replicate a",
                "# task left unmeasured, and writes the report and figures.",
                *setup,
                f"{command} analyze {shlex.quote(name)} "
                + " ".join(shlex.quote(argument) for argument in report_arguments),
                "",
            ]
        )
    )
    return Submission(folder, task_file, array, report, len(tasks), report_output=output)


def submit(submission: Submission) -> Submission:
    """Submit the array, then the report job to start once every task has ended."""
    from polyzymd.analyses.exceptions import ProtocolError

    if shutil.which("sbatch") is None:
        raise ProtocolError(
            "--submit needs sbatch, which is not on PATH.",
            hint="Load your cluster's SLURM module first, such as 'module load slurm/blanca', "
            f"or submit the scripts in {submission.folder} yourself; --dry-run only writes them.",
        )

    def sbatch(*arguments: str) -> str:
        result = subprocess.run(
            ["sbatch", "--parsable", *arguments], capture_output=True, text=True
        )
        if result.returncode != 0:
            raise ProtocolError(
                f"sbatch refused the job: {result.stderr.strip() or result.stdout.strip()}",
                hint=f"Check the SLURM settings, or submit the scripts in {submission.folder} "
                "yourself.",
            )
        return result.stdout.strip().split(";")[0]

    submission.array_id = sbatch(str(submission.array))
    submission.report_id = sbatch(
        f"--dependency=afterany:{submission.array_id}", str(submission.report)
    )
    return submission
