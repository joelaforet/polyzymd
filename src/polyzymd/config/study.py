"""The study folder: one directory holding everything a study needs.

A study is the directory that holds a ``study.yaml`` file. Inside it, each
condition is its own ``polyzymd init`` project under ``conditions/``, each
comparison is its own ``polyzymd compare init`` project under
``comparisons/``, and analyses written against the PolyzyMD API live in one
``analyses/`` folder that every comparison shares. Zipping the folder, or
keeping it in git, keeps every input and every analysis of the study together.

A comparison that is not inside a study still finds analyses in an
``analyses/`` folder next to its ``comparison.yaml``.
"""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict, Field

STUDY_FILE = "study.yaml"
"""Name of the file that marks the root of a study."""

ANALYSES_DIR = "analyses"
"""Default folder, relative to the study root, holding the study's analyses."""


class StudyConfig(BaseModel):
    """Contents of ``study.yaml``."""

    model_config = ConfigDict(extra="forbid")

    name: str = Field(min_length=1, description="Study name")
    description: str | None = Field(default=None, description="What the study asks")
    analyses: Path = Field(
        default=Path(ANALYSES_DIR),
        description="Folder holding the study's analyses, relative to the study root",
    )

    @classmethod
    def from_yaml(cls, path: Path | str) -> StudyConfig:
        """Read and validate a ``study.yaml`` file."""
        data = yaml.safe_load(Path(path).read_text()) or {}
        return cls.model_validate(data)


def find_study_root(start: Path | str) -> Path | None:
    """Return the nearest directory at or above ``start`` holding ``study.yaml``.

    Parameters
    ----------
    start : Path or str
        A file or directory inside the study.

    Returns
    -------
    Path or None
        The study root, or ``None`` when no ancestor holds ``study.yaml``.
    """
    path = Path(start).resolve()
    directory = path if path.is_dir() else path.parent
    for candidate in (directory, *directory.parents):
        if (candidate / STUDY_FILE).is_file():
            return candidate
    return None


def analyses_directory(start: Path | str) -> Path | None:
    """Return the analyses folder that applies to a file or directory.

    Inside a study this is the folder ``study.yaml`` names, ``analyses/`` by
    default. Outside a study it is an ``analyses/`` folder in the directory of
    ``start`` itself.

    Parameters
    ----------
    start : Path or str
        A ``comparison.yaml``, a simulation ``config.yaml``, or a directory.

    Returns
    -------
    Path or None
        The folder, or ``None`` when it does not exist.
    """
    root = find_study_root(start)
    if root is not None:
        directory = root / StudyConfig.from_yaml(root / STUDY_FILE).analyses
    else:
        path = Path(start).resolve()
        directory = (path if path.is_dir() else path.parent) / ANALYSES_DIR
    return directory if directory.is_dir() else None
