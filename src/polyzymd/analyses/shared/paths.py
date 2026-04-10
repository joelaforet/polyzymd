"""Path utilities for the analysis plugin system."""

from __future__ import annotations

import re


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
