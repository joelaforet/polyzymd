from __future__ import annotations

import re
from pathlib import Path


def resolve_condition_output_dir(
    source_path: Path | None,
    label: str,
    analysis_subdir: str,
) -> Path | None:
    """Resolve a condition-specific output directory for analysis results.

    When running in comparison mode (source_path is set), returns a
    condition-specific path under the comparison project directory. This
    prevents cache collisions when multiple conditions share the same
    projects_directory.

    When source_path is None (standalone/programmatic usage), returns None.

    Parameters
    ----------
    source_path : Path | None
        Path to the comparison configuration file.
    label : str
        Condition label.
    analysis_subdir : str
        Analysis subdirectory name.

    Returns
    -------
    Path | None
        Condition-specific output directory, or None when source_path is not set.
    """
    if source_path is None:
        return None

    comparison_dir = source_path.parent
    return comparison_dir / "analysis" / sanitize_label(label) / analysis_subdir


def resolve_analysis_dir(
    projects_dir: Path,
    analysis_subdir: str,
    cond_config_path: Path | None = None,
    source_path: Path | None = None,
    label: str | None = None,
) -> Path:
    """Get analysis directory with multi-location fallback.

    Priority order:

    1. Comparison-mode per-condition dir (if source_path and label provided)
    2. projects_dir / analysis_subdir
    3. cond_config_path.parent / analysis_subdir (if provided and exists)

    Returns the first existing path, or the primary path for error messages.

    Parameters
    ----------
    projects_dir : Path
        Root project directory.
    analysis_subdir : str
        Analysis subdirectory name.
    cond_config_path : Path | None, optional
        Condition configuration path, by default None.
    source_path : Path | None, optional
        Comparison configuration path, by default None.
    label : str | None, optional
        Condition label for comparison-mode resolution, by default None.

    Returns
    -------
    Path
        The resolved analysis directory.
    """
    if source_path is not None and label is not None:
        condition_dir = resolve_condition_output_dir(source_path, label, analysis_subdir)
        if condition_dir is not None and condition_dir.exists():
            return condition_dir

    primary = projects_dir / analysis_subdir
    if primary.exists():
        return primary

    if cond_config_path is not None:
        fallback = cond_config_path.parent / analysis_subdir
        if fallback.exists():
            return fallback

    return primary


def resolve_aggregated_dir(analysis_dir: Path) -> Path:
    """Get the aggregated results subdirectory.

    Parameters
    ----------
    analysis_dir : Path
        Base analysis directory.

    Returns
    -------
    Path
        Aggregated results subdirectory.
    """
    return analysis_dir / "aggregated"


def sanitize_label(label: str) -> str:
    """Convert a condition label to a filesystem-safe directory name.

    Replaces ``%`` with ``pct``, spaces with underscores, strips remaining
    non-alphanumeric chars (except hyphens, underscores, dots), and collapses
    consecutive underscores.

    Parameters
    ----------
    label : str
        Original condition label.

    Returns
    -------
    str
        Filesystem-safe label.
    """
    sanitized = label.strip()
    sanitized = sanitized.replace("%", "pct")
    sanitized = sanitized.replace(" ", "_")
    sanitized = re.sub(r"[^\w\-.]", "_", sanitized)
    sanitized = re.sub(r"_+", "_", sanitized)
    return sanitized.strip("_")
