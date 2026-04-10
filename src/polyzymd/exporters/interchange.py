"""Engine-agnostic export dispatch for PolyzyMD systems.

Routes export requests to engine-specific exporters based on format string.
Supported formats: gromacs (stable), lammps (planned), amber (planned).
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from openff.interchange import Interchange

    from polyzymd.config.schema import SimulationConfig
    from polyzymd.core.atom_groups import SystemComponentInfo

# Supported export formats
SUPPORTED_FORMATS = ("gromacs", "lammps", "amber")


class ExportFormat(str, Enum):
    """Supported MD engine export formats."""

    GROMACS = "gromacs"
    LAMMPS = "lammps"
    AMBER = "amber"


def get_supported_formats() -> tuple[str, ...]:
    """Return the tuple of supported export format names.

    Returns
    -------
    tuple[str, ...]
        Supported export format names in canonical lowercase form
    """

    return SUPPORTED_FORMATS


def export_system(
    interchange: "Interchange",
    config: "SimulationConfig",
    output_dir: str | Path,
    fmt: str | ExportFormat,
    component_info: "SystemComponentInfo | None" = None,
    prefix: str | None = None,
    gmx_command: str = "gmx",
) -> dict[str, Any]:
    """Export a parameterized system to a requested MD engine format.

    Parameters
    ----------
    interchange : Interchange
        OpenFF Interchange object with the parameterized system
    config : SimulationConfig
        PolyzyMD simulation configuration
    output_dir : str or Path
        Directory where exported files will be written
    fmt : str or ExportFormat
        Target MD engine format (``"gromacs"``, ``"lammps"``, ``"amber"``)
    component_info : SystemComponentInfo or None, optional
        Component metadata used for position restraints (GROMACS-only for now)
    prefix : str or None, optional
        Filename prefix for exported files. If None, exporter derives a default
    gmx_command : str, optional
        GROMACS executable command used by generated run scripts, by default ``"gmx"``

    Returns
    -------
    dict[str, Any]
        Mapping of export artifact keys to generated file paths

    Raises
    ------
    ValueError
        If ``fmt`` is not one of the supported export formats
    NotImplementedError
        If the requested format is recognized but not yet implemented
    """

    fmt_str = fmt.value if isinstance(fmt, ExportFormat) else fmt.lower().strip()

    if fmt_str not in SUPPORTED_FORMATS:
        raise ValueError(
            f"Unsupported export format: {fmt_str!r}. "
            f"Supported formats: {', '.join(SUPPORTED_FORMATS)}"
        )

    if fmt_str == "gromacs":
        return _export_gromacs(interchange, config, output_dir, component_info, prefix, gmx_command)
    if fmt_str == "lammps":
        return _export_lammps(interchange, config, output_dir, prefix)
    if fmt_str == "amber":
        return _export_amber(interchange, config, output_dir, prefix)

    # Defensive fallback for future format handling changes
    raise ValueError(f"Unhandled export format: {fmt_str!r}")


def _export_gromacs(
    interchange: "Interchange",
    config: "SimulationConfig",
    output_dir: str | Path,
    component_info: "SystemComponentInfo | None",
    prefix: str | None,
    gmx_command: str,
) -> dict[str, Any]:
    """Export to GROMACS format via :class:`GromacsExporter`.

    Parameters
    ----------
    interchange : Interchange
        Parameterized system to export
    config : SimulationConfig
        Simulation configuration used by the exporter
    output_dir : str or Path
        Destination directory for output files
    component_info : SystemComponentInfo or None
        Component metadata used to generate position restraints
    prefix : str or None
        Optional filename prefix
    gmx_command : str
        GROMACS executable command for run-script generation

    Returns
    -------
    dict[str, Any]
        GROMACS export artifacts keyed by type
    """

    from polyzymd.exporters.gromacs import GromacsExporter

    exporter = GromacsExporter(
        interchange=interchange,
        config=config,
        component_info=component_info,
    )
    return exporter.export(
        output_dir=output_dir,
        prefix=prefix,
        gmx_command=gmx_command,
    )


def _export_lammps(
    interchange: "Interchange",
    config: "SimulationConfig",
    output_dir: str | Path,
    prefix: str | None,
) -> dict[str, Any]:
    """Export to LAMMPS format.

    Parameters
    ----------
    interchange : Interchange
        Parameterized system to export
    config : SimulationConfig
        Simulation configuration
    output_dir : str or Path
        Destination directory for output files
    prefix : str or None
        Optional filename prefix

    Returns
    -------
    dict[str, Any]
        Placeholder return type for future implementation

    Raises
    ------
    NotImplementedError
        LAMMPS export is planned for v1.4.0
    """

    _ = (interchange, config, output_dir, prefix)
    raise NotImplementedError(
        "LAMMPS export is planned for v1.4.0. "
        "See https://github.com/Shirts-Lab-Linux/polyzymd/issues for updates."
    )


def _export_amber(
    interchange: "Interchange",
    config: "SimulationConfig",
    output_dir: str | Path,
    prefix: str | None,
) -> dict[str, Any]:
    """Export to AMBER format.

    Parameters
    ----------
    interchange : Interchange
        Parameterized system to export
    config : SimulationConfig
        Simulation configuration
    output_dir : str or Path
        Destination directory for output files
    prefix : str or None
        Optional filename prefix

    Returns
    -------
    dict[str, Any]
        Placeholder return type for future implementation

    Raises
    ------
    NotImplementedError
        AMBER export is planned for v1.4.0
    """

    _ = (interchange, config, output_dir, prefix)
    raise NotImplementedError(
        "AMBER export is planned for v1.4.0. "
        "See https://github.com/Shirts-Lab-Linux/polyzymd/issues for updates."
    )
