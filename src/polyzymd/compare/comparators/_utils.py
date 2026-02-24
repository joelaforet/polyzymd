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
