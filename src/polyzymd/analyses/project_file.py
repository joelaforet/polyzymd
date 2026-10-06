"""Read ``project.yaml``, the analyses a paper runs in each of its studies.

A project folder holds one paper: its ``project.yaml`` lists the study
folders, one per protein, and the analyses every study runs; see the
"Projects and studies" explanation page. A study finds its project in the
folder above it (:func:`find_project`), and :func:`project_analyses` gives
the project's analyses for that study, with the study's ``region`` and
``structure`` names resolved (:func:`resolve_names`).
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.statistics_plan import read_stats_plan

#: File name that marks a project folder.
PROJECT_FILE = "project.yaml"

_TOP_KEYS = ("polyzymd", "studies", "analyses", "stats", "metadata")
#: Key of a project analysis that limits it to some studies.
STUDIES_KEY = "studies"
# Names are letters, digits, _ and -, and may start with a digit (4TGL_open).
_REGION = re.compile(r"\bregion\s+([\w][\w-]*)")
_STRUCTURE = re.compile(r"^\s*structure\s+(\S+)\s*$")


@dataclass(frozen=True)
class ProjectFile:
    """The contents of one ``project.yaml``.

    ``studies`` maps each study's label to its folder. ``analyses`` holds
    each analysis entry as written, with its ``studies:`` list (``None`` for
    every study) under :data:`STUDIES_KEY` removed into ``runs_in``.
    """

    path: Path
    studies: dict[str, Path]
    analyses: dict[str, dict[str, Any]]
    runs_in: dict[str, list[str] | None]
    polyzymd: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    #: The paper's statistical plan (``stats: {plan: file.py:function}``).
    stats: Any = None

    @property
    def root(self) -> Path:
        """The project folder: the directory holding ``project.yaml``."""
        return self.path.parent

    def label_of(self, study_dir: Path) -> str | None:
        """Return the label of the study in ``study_dir``, or ``None`` when it is not listed."""
        study_dir = Path(study_dir).resolve()
        for label, folder in self.studies.items():
            if folder == study_dir:
                return label
        return None


def find_project_file(path: str | Path) -> Path:
    """Return the ``project.yaml`` that ``path`` names: the file itself, or the one in that folder."""
    path = Path(path).expanduser().resolve()
    file = path / PROJECT_FILE if path.is_dir() else path
    if not file.is_file():
        raise ProtocolError(
            f"No {PROJECT_FILE} at {path}.",
            hint="Pass the project folder or its project.yaml.",
        )
    return file


def load_project_file(path: str | Path) -> ProjectFile:
    """Read and check ``project.yaml``.

    Study folders are resolved against the project folder and must each hold
    a ``study.yaml``. Every ``studies:`` list of an analysis must name listed
    studies.

    Raises
    ------
    ProtocolError
        If the file is missing or not YAML, a key is unknown, a study folder
        has no ``study.yaml``, or an analysis lists an unknown study.
    """
    import yaml

    from polyzymd.analyses.study_file import STUDY_FILE, _unknown

    file = find_project_file(path)
    try:
        raw = yaml.safe_load(file.read_text()) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ProtocolError(
            f"Cannot read {file}: {exc}", hint="Check that it is valid YAML."
        ) from exc
    if not isinstance(raw, Mapping):
        raise ProtocolError(f"{file} must be a YAML mapping.", hint="See the documented example.")
    _unknown(raw, _TOP_KEYS, str(file))
    studies_raw = raw.get("studies")
    if not isinstance(studies_raw, Mapping) or not studies_raw:
        raise ProtocolError(
            f"{file}: studies must map each study's label to its folder.",
            hint="For example 'studies: {lipa363: lipa363, rml333: rml333}'.",
        )
    studies = {
        str(label): (file.parent / Path(str(folder)).expanduser()).resolve()
        for label, folder in studies_raw.items()
    }
    for label, folder in studies.items():
        if not (folder / STUDY_FILE).is_file():
            raise ProtocolError(
                f"{file}: study {label} has no {STUDY_FILE} in {folder}.",
                hint="Point each study at a folder holding its study.yaml.",
            )
    analyses_raw = raw.get("analyses") or {}
    if not isinstance(analyses_raw, Mapping):
        raise ProtocolError(
            f"{file}: analyses must map each run name to its settings.",
            hint="For example 'analyses: {rg: {}}'.",
        )
    analyses: dict[str, dict[str, Any]] = {}
    runs_in: dict[str, list[str] | None] = {}
    for run, entry in analyses_raw.items():
        entry = dict(entry or {})
        listed = entry.pop(STUDIES_KEY, None)
        if listed is not None:
            if not isinstance(listed, list) or not listed:
                raise ProtocolError(
                    f"{file}: analyses.{run}.studies must list study labels.",
                    hint=f"For example 'studies: [{next(iter(studies))}]'.",
                )
            unknown = [str(name) for name in listed if str(name) not in studies]
            if unknown:
                raise ProtocolError(
                    f"{file}: analyses.{run} lists studies {unknown} that the project does not.",
                    hint=f"Use labels from studies: {', '.join(studies)}.",
                )
            listed = [str(name) for name in listed]
        analyses[str(run)] = entry
        runs_in[str(run)] = listed
    metadata = raw.get("metadata") or {}
    if not isinstance(metadata, Mapping):
        raise ProtocolError(f"{file}: metadata must be a mapping.", hint="Leave it out for now.")
    version = raw.get("polyzymd")
    return ProjectFile(
        path=file,
        studies=studies,
        analyses=analyses,
        runs_in=runs_in,
        polyzymd=None if version is None else str(version),
        metadata=dict(metadata),
        stats=read_stats_plan(raw.get("stats"), f"{file}: stats", file.parent),
    )


def find_project(study_dir: Path) -> tuple[ProjectFile, str] | None:
    """Return the project whose ``project.yaml`` in the folder above lists ``study_dir``, with its label."""
    candidate = Path(study_dir).resolve().parent / PROJECT_FILE
    if not candidate.is_file():
        return None
    project = load_project_file(candidate)
    label = project.label_of(study_dir)
    return None if label is None else (project, label)


def resolve_names(
    value: Any, regions: Mapping[str, str], structures: Mapping[str, Path], where: str
) -> Any:
    """Replace ``region <name>`` and ``structure <name>`` with a study's selection and file.

    ``region <name>`` anywhere in a string becomes the region's selection in
    parentheses; a string that is only ``structure <name>`` becomes the path
    of that structure. Mappings and lists are resolved item by item; other
    values are returned as they are.

    Raises
    ------
    ProtocolError
        If a name is not one of the study's regions or structures.
    """
    if isinstance(value, Mapping):
        return {key: resolve_names(item, regions, structures, where) for key, item in value.items()}
    if isinstance(value, list):
        return [resolve_names(item, regions, structures, where) for item in value]
    if not isinstance(value, str):
        return value
    if value.strip().startswith("structure ") and not _STRUCTURE.match(value):
        raise ProtocolError(
            f"{where}: {value!r} must be 'structure <name>' with one name.",
            hint="Write the value as structure followed by one name from structures:.",
        )
    structure = _STRUCTURE.match(value)
    if structure:
        name = structure.group(1)
        if name not in structures:
            raise ProtocolError(
                f"{where} uses structure {name}, which this study does not define.",
                hint="Add it under structures: in the study.yaml"
                + (f" (it has {', '.join(structures)})" if structures else "")
                + ", or limit the analysis with studies: in project.yaml.",
            )
        return str(structures[name])

    def region(match: re.Match) -> str:
        name = match.group(1)
        if name not in regions:
            raise ProtocolError(
                f"{where} uses region {name}, which this study does not define.",
                hint="Add it under regions: in the study.yaml"
                + (f" (it has {', '.join(regions)})" if regions else "")
                + ", or limit the analysis with studies: in project.yaml.",
            )
        return f"({regions[name]})"

    return _REGION.sub(region, value)
