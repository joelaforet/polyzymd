"""Version helpers shared by build, simulation and analysis provenance records."""

from __future__ import annotations

import os
import socket
from importlib.metadata import PackageNotFoundError
from typing import Any


def get_polyzymd_version() -> str:
    """Return the installed PolyzyMD version.

    Returns
    -------
    str
        Installed package version, or ``"unknown"`` when package metadata is
        unavailable in an editable or source-tree execution context.
    """
    try:
        from importlib.metadata import version

        return version("polyzymd")
    except (ImportError, PackageNotFoundError):
        return "unknown"


def get_openmm_version() -> str | None:
    """Return the OpenMM version string, or ``None`` when OpenMM is not importable."""
    try:
        from openmm import version

        return str(version.full_version)
    except ImportError:
        return None


def runtime_provenance(context: Any = None) -> dict[str, Any]:
    """Describe the software and host that a simulation phase runs under.

    Parameters
    ----------
    context : openmm.Context, optional
        The Context that runs the phase. When given, ``openmm_platform``
        records its platform and property values
        (:func:`polyzymd.simulation.platform.platform_record`).

    Returns
    -------
    dict
        ``polyzymd_version``, ``openmm_version``, ``pixi_environment``
        (``$PIXI_ENVIRONMENT_NAME``), ``hostname`` and ``slurm_job_id``
        (``$SLURM_JOB_ID``).  Unavailable values are ``None``.  Recorded in
        ``progress.json`` segment records and ``production_N_parameters.json``
        so that a restart chain that silently switched environment or OpenMM
        build can be detected after the fact.
    """
    try:
        hostname: str | None = socket.gethostname()
    except OSError:
        hostname = None
    found = {
        "polyzymd_version": get_polyzymd_version(),
        "openmm_version": get_openmm_version(),
        "pixi_environment": os.environ.get("PIXI_ENVIRONMENT_NAME"),
        "hostname": hostname,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    }
    if context is not None:
        from polyzymd.simulation.platform import platform_record

        found["openmm_platform"] = platform_record(context)
    return found


RECORD_PROVENANCE_KEYS = ("polyzymd_version", "openmm_version", "pixi_environment")


def record_provenance(context: Any = None) -> dict[str, Any]:
    """Subset of :func:`runtime_provenance` stored on progress records.

    With ``context``, it also holds ``openmm_platform``.
    """
    full = runtime_provenance(context)
    return {key: full[key] for key in (*RECORD_PROVENANCE_KEYS, "openmm_platform") if key in full}
