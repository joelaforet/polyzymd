"""GROMACS binary discovery helpers."""

from __future__ import annotations

import shutil
from os import getenv


def resolve_gromacs_binary(config: object = None) -> str:
    """Resolve the GROMACS binary path.

    Precedence:

    - ``config.gromacs.gmx_binary``
    - ``$GMX_BIN`` environment variable
    - ``gmx`` from ``PATH``
    - ``gmx_mpi`` from ``PATH``

    Parameters
    ----------
    config : object, optional
        Simulation config that may include a gromacs section.

    Returns
    -------
    str
        Resolved binary path or executable name.

    Raises
    ------
    FileNotFoundError
        If no candidate binary can be resolved.
    """
    if config is not None:
        gromacs_cfg = getattr(config, "gromacs", None)
        if gromacs_cfg is not None:
            cfg_binary = getattr(gromacs_cfg, "gmx_binary", None)
            if cfg_binary:
                return str(cfg_binary)

    env_binary = getenv("GMX_BIN")
    if env_binary:
        return env_binary

    gmx_path = shutil.which("gmx")
    if gmx_path:
        return gmx_path

    gmx_mpi_path = shutil.which("gmx_mpi")
    if gmx_mpi_path:
        return gmx_mpi_path

    raise FileNotFoundError(
        "Could not resolve GROMACS binary. Set config.gromacs.gmx_binary, "
        "set GMX_BIN, or install gmx/gmx_mpi on PATH"
    )
