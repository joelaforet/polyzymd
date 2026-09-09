"""Version helpers shared by build, simulation and analysis provenance records."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError


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
