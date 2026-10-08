"""Version helpers shared by build, simulation and analysis provenance records."""

from __future__ import annotations

import os
import socket
import sys
from importlib.metadata import PackageNotFoundError
from pathlib import Path
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


def package_version(module: str) -> str | None:
    """Return the version of the package that ``import module`` loads, or ``None``.

    A package that reports version ``0.0.0`` (a conda build without its
    version in the metadata, such as OpenFF Interchange) has the version of
    its conda package record (``conda-meta/<name>-<version>-<build>.json`` of
    the environment), or ``None``.
    """
    import json

    try:
        imported = __import__(module, fromlist=["__version__"])
        found = str(getattr(imported, "__version__", None) or getattr(imported, "version", None))
        version: str | None = None if found == "0.0.0" else found
    except Exception:  # noqa: BLE001 - an absent or broken package is recorded as absent
        version = None
    if version is None:
        name = module.lower().replace(".", "-")
        for record in Path(sys.prefix, "conda-meta").glob(f"{name}-*.json"):
            try:
                found_record = json.loads(record.read_text())
            except (OSError, ValueError):
                continue
            if isinstance(found_record, dict) and found_record.get("name") == name:
                version = found_record.get("version")
    return version


def build_versions() -> dict[str, str | None]:
    """Return the versions of the OpenFF packages that parameterize a system, for either engine.

    ``openff_toolkit_version`` and ``openff_interchange_version``
    (:func:`package_version`), as ``build_manifest.json`` records them.
    """
    return {
        "openff_toolkit_version": package_version("openff.toolkit"),
        "openff_interchange_version": package_version("openff.interchange"),
    }


def pixi_workspace() -> Path | None:
    """Return the pixi workspace whose environment runs Python, or ``None``.

    Pixi installs an environment in ``<workspace>/.pixi/envs/<name>``, which is
    ``sys.prefix`` when Python runs from it, with or without ``pixi run``.
    """
    prefix = Path(sys.prefix)
    if prefix.parent.name == "envs" and prefix.parent.parent.name == ".pixi":
        return prefix.parents[2]
    return None


def pixi_environment() -> str | None:
    """Return the name of the pixi environment that runs Python, or ``None``.

    ``$PIXI_ENVIRONMENT_NAME`` when ``pixi run`` set it, otherwise the name of
    the environment folder when Python runs from a pixi workspace
    (:func:`pixi_workspace`).
    """
    return os.environ.get("PIXI_ENVIRONMENT_NAME") or (
        Path(sys.prefix).name if pixi_workspace() else None
    )


def runtime_provenance(simulation: Any = None) -> dict[str, Any]:
    """Describe the software and host that a simulation phase runs under.

    Parameters
    ----------
    simulation : openmm.app.Simulation, optional
        The Simulation that runs the phase. When given, ``openmm_platform``
        records the platform and property values of its Context
        (:func:`polyzymd.simulation.platform.platform_record`).

    Returns
    -------
    dict
        ``polyzymd_version``, ``openmm_version``, ``pixi_environment``
        (:func:`pixi_environment`), ``hostname`` and ``slurm_job_id``
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
        "pixi_environment": pixi_environment(),
        "hostname": hostname,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    }
    if simulation is not None:
        from polyzymd.simulation.platform import platform_record

        found["openmm_platform"] = platform_record(simulation.context)
    return found


RECORD_PROVENANCE_KEYS = ("polyzymd_version", "openmm_version", "pixi_environment")


def record_provenance(simulation: Any = None) -> dict[str, Any]:
    """Subset of :func:`runtime_provenance` stored on progress records.

    With ``simulation``, it also holds ``openmm_platform``.
    """
    full = runtime_provenance(simulation)
    return {key: full[key] for key in (*RECORD_PROVENANCE_KEYS, "openmm_platform") if key in full}
