"""SLURM submission helpers."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path


def make_log_folder(script_path: Path | str) -> None:
    """Create the folder of the script's absolute ``#SBATCH --output`` log, which sbatch does not."""
    try:
        text = Path(script_path).read_text()
    except OSError:
        return  # sbatch says that the script is missing
    for line in text.splitlines():
        match = re.match(r"^#SBATCH\s+--output[=\s]\s*\"?([^\"]+?)\"?\s*$", line)
        if match and Path(match.group(1)).is_absolute():
            Path(match.group(1)).parent.mkdir(parents=True, exist_ok=True)


def require_sbatch(script_path: Path | str) -> None:
    """Raise when ``sbatch`` is not on ``PATH``, naming the script that was written."""
    if shutil.which("sbatch") is None:
        raise RuntimeError(
            "sbatch not found. Load your scheduler module first, for example "
            f"`ml slurm/blanca` on CU Boulder Blanca. The job script was written to {script_path}; "
            f"submit it with `sbatch --export=NONE {script_path}` or run submit again."
        )


def run_sbatch(script_path: Path | str) -> subprocess.CompletedProcess[str]:
    """Submit a script with ``sbatch --export=NONE``.

    The job starts from a clean login environment, as OpenMM jobs do, and
    runs its own ``module load`` lines. Nothing is loaded on the submitting
    host. The folder of the job's log is created first (:func:`make_log_folder`).

    Parameters
    ----------
    script_path : Path or str
        Path to the SLURM batch script.

    Returns
    -------
    subprocess.CompletedProcess[str]
        Result of the ``sbatch`` invocation.

    Raises
    ------
    RuntimeError
        If ``sbatch`` is not on ``PATH``.
    """
    require_sbatch(script_path)
    make_log_folder(script_path)
    return subprocess.run(
        ["sbatch", "--export=NONE", str(script_path)], capture_output=True, text=True, check=False
    )
