"""Path utilities for the analysis plugin system."""

from __future__ import annotations

import re
from collections.abc import Sequence


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


def format_replicate_cache_token(replicates: Sequence[int]) -> str:
    """Format replicate IDs for cache filenames without range collisions.

    Contiguous replicate IDs are compacted as a range, while non-contiguous
    IDs are listed explicitly so ``(1, 3)`` cannot collide with ``(1, 2, 3)``.

    Parameters
    ----------
    replicates : Sequence[int]
        Iterable of replicate IDs.

    Returns
    -------
    str
        Cache-safe token such as ``"reps1-3"`` or ``"reps1_3"``.
    """

    try:
        reps = sorted({int(rep) for rep in replicates})
    except (TypeError, ValueError):
        return "no_replicates"
    if not reps:
        return "no_replicates"
    if reps == list(range(reps[0], reps[-1] + 1)):
        return f"reps{reps[0]}-{reps[-1]}"
    return "reps" + "_".join(str(rep) for rep in reps)
