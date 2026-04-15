"""GROMACS binary discovery and classification helpers."""

from __future__ import annotations

import shutil
from os import getenv
from pathlib import Path


def is_mpi_binary(binary: str) -> bool:
    """Check whether a GROMACS binary is a real-MPI build.

    Real-MPI builds (``gmx_mpi``, ``gmx_mpi_d``) do not support the
    ``-ntmpi`` flag and require an external MPI launcher (``mpirun``).
    Thread-MPI builds (``gmx``, ``gmx_d``) use ``-ntmpi`` natively.

    Parameters
    ----------
    binary : str
        Binary name or full path.

    Returns
    -------
    bool
        ``True`` for real-MPI builds, ``False`` for thread-MPI builds.
    """
    return "_mpi" in Path(binary).name


def resolve_gromacs_binary(config: object = None, *, gpu: bool = False) -> str:
    """Resolve the GROMACS binary path.

    Precedence:

    - ``config.gromacs.gmx_binary``
    - ``$GMX_BIN`` environment variable
    - ``gmx`` from ``PATH``
    - ``gmx_mpi`` from ``PATH`` (CPU mode only)

    In GPU mode the function **rejects** real-MPI binaries (``gmx_mpi``)
    because CUDA support on most HPC clusters is compiled only into the
    thread-MPI ``gmx`` binary.  PATH discovery also skips ``gmx_mpi``
    when ``gpu=True``.

    Parameters
    ----------
    config : object, optional
        Simulation config that may include a gromacs section.
    gpu : bool, optional
        Whether GPU mode is enabled.  When ``True``, real-MPI binaries
        are rejected with a clear error.

    Returns
    -------
    str
        Resolved binary path or executable name.

    Raises
    ------
    FileNotFoundError
        If no candidate binary can be resolved.
    ValueError
        If GPU mode is enabled but the resolved binary is a real-MPI
        build (which lacks CUDA support).
    """
    if config is not None:
        gromacs_cfg = getattr(config, "gromacs", None)
        if gromacs_cfg is not None:
            cfg_binary = getattr(gromacs_cfg, "gmx_binary", None)
            if cfg_binary:
                binary = str(cfg_binary)
                if gpu and is_mpi_binary(binary):
                    raise ValueError(
                        f"GPU mode requires a thread-MPI GROMACS binary (gmx), "
                        f"but config specifies '{binary}'. Real-MPI builds "
                        f"(gmx_mpi) typically lack CUDA support. Set "
                        f"gmx_binary to 'gmx' or remove it to auto-resolve."
                    )
                return binary

    env_binary = getenv("GMX_BIN")
    if env_binary:
        if gpu and is_mpi_binary(env_binary):
            raise ValueError(
                f"GPU mode requires a thread-MPI GROMACS binary (gmx), "
                f"but $GMX_BIN is set to '{env_binary}'. Real-MPI builds "
                f"(gmx_mpi) typically lack CUDA support. Unset GMX_BIN or "
                f"point it to a thread-MPI binary."
            )
        return env_binary

    gmx_path = shutil.which("gmx")
    if gmx_path:
        return gmx_path

    # In GPU mode, skip gmx_mpi fallback — it won't have CUDA.
    if not gpu:
        gmx_mpi_path = shutil.which("gmx_mpi")
        if gmx_mpi_path:
            return gmx_mpi_path

    raise FileNotFoundError(
        "Could not resolve GROMACS binary. Set config.gromacs.gmx_binary, "
        "set GMX_BIN, or install gmx on PATH"
        + (" (GPU mode requires the thread-MPI gmx binary)" if gpu else "")
    )
