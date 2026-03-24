"""Shared utility functions for comparator modules.

Extracted from individual comparators to eliminate duplicated file-location
logic across contacts.py, exposure.py, binding_free_energy.py, etc.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def find_analysis_dir(
    sim_config: Any,
    analysis_subdir: str = "analysis",
    cond_config_path: Path | None = None,
) -> Path:
    """Get analysis directory with fallback to condition config parent.

    Checks multiple locations in order:
    1. sim_config.output.projects_directory / analysis_subdir
    2. cond_config_path.parent / analysis_subdir (if provided and exists)

    This allows cached results to be found even when the original
    projects_directory points to a remote/unavailable location.

    Parameters
    ----------
    sim_config : SimulationConfig
        Simulation configuration object.
    analysis_subdir : str, optional
        Subdirectory path relative to the project root, by default ``"analysis"``.
        Examples: ``"analysis/contacts"``, ``"analysis"``.
    cond_config_path : Path, optional
        Path to the condition's config.yaml file. Used as fallback
        location for finding cached results.

    Returns
    -------
    Path
        Analysis directory path (primary location, or fallback if it exists).
    """
    # Primary location: projects_directory
    primary_dir = sim_config.output.projects_directory / analysis_subdir
    if primary_dir.exists():
        return primary_dir

    # Fallback: config file's parent directory
    if cond_config_path is not None:
        fallback_dir = cond_config_path.parent / analysis_subdir
        if fallback_dir.exists():
            return fallback_dir

    # Return primary path even if doesn't exist (for error messages)
    return primary_dir


def find_replicate_result(
    sim_config: Any,
    replicate: int,
    result_filename: str,
    analysis_subdir: str = "analysis/contacts",
    cond_config_path: Path | None = None,
) -> Path:
    """Find path to an existing replicate result file.

    Checks multiple locations in order:
    1. sim_config.output.projects_directory / analysis_subdir / result_filename
    2. cond_config_path.parent / analysis_subdir / result_filename (if provided)

    This allows cached results to be found even when the original
    projects_directory points to a remote/unavailable location.

    Parameters
    ----------
    sim_config : SimulationConfig
        Simulation configuration object.
    replicate : int
        Replicate number (for documentation; not used directly since
        *result_filename* should already encode the replicate).
    result_filename : str
        Filename to look for, e.g. ``f"contacts_rep{replicate}.json"``.
    analysis_subdir : str, optional
        Subdirectory path relative to the project root,
        by default ``"analysis/contacts"``.
    cond_config_path : Path, optional
        Path to the condition's config.yaml file. Used as fallback
        location for finding cached results.

    Returns
    -------
    Path
        Path to the result file (primary location, or fallback if it exists).
    """
    # Primary location: projects_directory
    primary = sim_config.output.projects_directory / analysis_subdir / result_filename
    if primary.exists():
        return primary

    # Fallback: config file's parent directory
    if cond_config_path is not None:
        fallback = cond_config_path.parent / analysis_subdir / result_filename
        if fallback.exists():
            return fallback

    # Return primary path even if doesn't exist (for error messages)
    return primary


def parse_equilibration_time(eq_string: str) -> tuple[float, str]:
    """Parse an equilibration time string like '10ns' or '500ps'.

    Returns the numeric value and unit as-is, without converting between
    units. If no unit suffix is found, defaults to ``"ns"``.

    Parameters
    ----------
    eq_string : str
        Equilibration time string, e.g. ``"10ns"``, ``"500ps"``, ``"10"``.

    Returns
    -------
    tuple[float, str]
        ``(value, unit)`` tuple, e.g. ``(10.0, "ns")`` or ``(500.0, "ps")``.
    """
    eq_str = eq_string.lower()
    if eq_str.endswith("ns"):
        return float(eq_str[:-2]), "ns"
    elif eq_str.endswith("ps"):
        return float(eq_str[:-2]), "ps"
    else:
        return float(eq_str), "ns"


def format_replicate_range(replicates: list[int]) -> str:
    """Format a list of replicate numbers into a compact string.

    Consecutive ranges are collapsed (e.g. ``[1, 2, 3]`` → ``"reps1-3"``),
    while non-consecutive lists are joined with underscores
    (e.g. ``[1, 3, 5]`` → ``"reps1_3_5"``).

    Parameters
    ----------
    replicates : list[int]
        Replicate numbers (need not be sorted).

    Returns
    -------
    str
        Formatted replicate string, e.g. ``"reps1-3"`` or ``"reps1_3_5"``.
    """
    reps = sorted(replicates)
    if reps == list(range(reps[0], reps[-1] + 1)):
        return f"reps{reps[0]}-{reps[-1]}"
    return "reps" + "_".join(map(str, reps))


def sanitize_label(label: str) -> str:
    """Convert a condition label to a filesystem-safe directory name.

    Replaces ``%`` with ``pct``, spaces with underscores, and strips any
    remaining characters that are not alphanumeric, hyphens, underscores,
    or dots.  Consecutive underscores are collapsed.

    Parameters
    ----------
    label : str
        Condition label, e.g. ``"SBMA-EGMA 25%"`` or
        ``"No Polymer (Control)"``.

    Returns
    -------
    str
        Sanitized string safe for use as a directory name,
        e.g. ``"SBMA-EGMA_25pct"`` or ``"No_Polymer_Control"``.
    """
    import re

    s = label.strip()
    s = s.replace("%", "pct")
    s = s.replace(" ", "_")
    # Keep only word chars (alphanumeric + underscore), hyphens, and dots
    s = re.sub(r"[^\w\-.]", "_", s)
    # Collapse consecutive underscores
    s = re.sub(r"_+", "_", s)
    return s.strip("_")
