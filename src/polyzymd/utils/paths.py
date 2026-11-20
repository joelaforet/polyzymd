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


def names_existing(path: str | os.PathLike[str], *, file: bool = False) -> bool:
    """Return whether ``path`` names an existing file or folder (only a file with ``file``).

    A setting such as a selection may be probed in case it names a file. Text
    that cannot be a path, such as a name part longer than the operating
    system allows (``File name too long``) or a null byte, names nothing, so
    the answer is False instead of an error.
    """
    try:
        candidate = Path(path)
        return candidate.is_file() if file else candidate.exists()
    except (OSError, ValueError):
        return False
