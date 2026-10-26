"""Paths as PolyzyMD prints them and writes them into job scripts.

A config keeps its output folders as written (``../../../runs/...``), so
its hash does not change. Where a path is printed or written into a job
script, its ``..`` parts are resolved as the operating system resolves them,
so a folder reached through a symlink is named where it really is.
"""

from __future__ import annotations

import os
from pathlib import Path


def tidy_path(path: str | os.PathLike[str]) -> Path:
    """Return ``path`` absolute, with ``..`` parts resolved by the operating system.

    A path without ``..`` is only made absolute, so a symlinked folder such as
    a scratch link keeps the name it was given.
    """
    absolute = Path(path).absolute()
    if ".." not in absolute.parts:
        return absolute
    try:
        return absolute.resolve()
    except (OSError, RuntimeError):
        return Path(os.path.normpath(absolute))


def shown_path(path: object) -> str:
    """Return ``path`` without ``..`` parts, for printing (:func:`tidy_path`)."""
    return str(tidy_path(str(path)))
