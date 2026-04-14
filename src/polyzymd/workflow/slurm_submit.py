"""SLURM submission helpers with optional module preloading."""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path


def run_sbatch(
    script_path: Path | str,
    *,
    module_load: str | None = None,
    extra_args: tuple[str, ...] = (),
) -> subprocess.CompletedProcess[str]:
    """Submit a script via ``sbatch``, optionally preloading HPC modules.

    On clusters where SLURM itself is loaded via ``module load`` (e.g.
    ``module load slurm/blanca``), the feature database is only available
    after the module is sourced. When *module_load* is provided the
    command is wrapped in ``bash -lc`` so that ``sbatch`` sees the full
    feature set.

    Parameters
    ----------
    script_path : Path or str
        Path to the SLURM batch script.
    module_load : str or None
        Shell command(s) to source before ``sbatch`` (e.g.
        ``"module load slurm/blanca gcc/11.2.0 gromacs/2024.2"``).
    extra_args : tuple of str
        Extra CLI arguments passed to ``sbatch`` (before the script path).

    Returns
    -------
    subprocess.CompletedProcess[str]
        Result of the ``sbatch`` invocation.
    """
    sbatch_parts = ["sbatch", *extra_args, str(script_path)]

    if module_load:
        sbatch_cmd = " ".join(shlex.quote(part) for part in sbatch_parts)
        shell_cmd = f"{module_load}; {sbatch_cmd}"
        return subprocess.run(
            ["bash", "-lc", shell_cmd],
            capture_output=True,
            text=True,
            check=False,
        )

    return subprocess.run(sbatch_parts, capture_output=True, text=True, check=False)
