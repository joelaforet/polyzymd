"""Read ``project.yaml``, the analyses a paper runs in each of its studies.

A project folder holds one paper: its ``project.yaml`` lists the study
folders and the analyses every study runs; see the
"Projects and studies" explanation page. This module contains:

- :func:`load_project_file`, which reads and checks ``project.yaml`` into a
  :class:`ProjectFile`.
- :func:`find_project`, which looks for the project of a study in the
  folder above the study.
- :func:`resolve_names`, which replaces ``region <name>`` and
  ``structure <name>`` in a project analysis's settings with one study's
  selection and file.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError

#: File name that marks a project folder.
PROJECT_FILE = "project.yaml"

_TOP_KEYS = ("studies", "analyses", "metadata")
#: Key of a project analysis that limits it to some studies.
STUDIES_KEY = "studies"
# Names are letters, digits, _ and -, and may start with a digit (4TGL_open).
_REGION = re.compile(r"\bregion\s+([\w][\w-]*)")
_STRUCTURE = re.compile(r"^\s*structure\s+(\S+)\s*$")
#: Setting keys whose values are atom selections, besides keys that contain "selection".
_SELECTION_KEYS = {
    "groups",
    "regions",
    "core",
    "target",
    "contexts",
    "donors",
    "hydrogens",
    "acceptors",
}


@dataclass(frozen=True)
class ProjectFile:
    """The contents of one ``project.yaml``, as read by :func:`load_project_file`.

    Attributes
    ----------
    path : Path
        Absolute path of the ``project.yaml``.
    studies : dict of str to Path
        Each study's label mapped to its absolute folder, in file order.
    analyses : dict of str to dict
        Each analysis entry as written, by run name, with its ``studies:``
        key (:data:`STUDIES_KEY`) removed.
    runs_in : dict of str to list of str or None
        For each run name, the labels its ``studies:`` key listed, or
        ``None`` when the analysis runs in every study.
    metadata : dict
        The ``metadata:`` mapping, used for citation and deposit files.
    """

    path: Path
    studies: dict[str, Path]
    analyses: dict[str, dict[str, Any]]
    runs_in: dict[str, list[str] | None]
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def root(self) -> Path:
        """The project folder: the directory holding ``project.yaml``."""
        return self.path.parent

    def label_of(self, study_dir: Path) -> str | None:
        """Return the label of the study in ``study_dir``, or ``None`` when it is not listed.

        Parameters
        ----------
        study_dir : Path
            A study folder; it is resolved before comparison with the
            project's study folders.

        Returns
        -------
        str or None
            The study's label, or ``None`` when no listed study has that folder.
        """
        study_dir = Path(study_dir).resolve()
        for label, folder in self.studies.items():
            if folder == study_dir:
                return label
        return None


def find_project_file(path: str | Path) -> Path:
    """Return the ``project.yaml`` that ``path`` names: the file itself, or the one in that folder.

    Parameters
    ----------
    path : str or Path
        A project folder or a ``project.yaml``; ``~`` is expanded.

    Returns
    -------
    Path
        The absolute path of the file.

    Raises
    ------
    ProtocolError
        If no such file exists.
    """
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

    The top-level keys may be ``studies``, ``analyses`` and
    ``metadata``; ``studies`` is required and must not be empty. Study
    folders are resolved against the project folder and must each hold a
    ``study.yaml``. Every ``studies:`` list of an analysis must be a
    non-empty list of listed study labels.

    Parameters
    ----------
    path : str or Path
        The project folder or its ``project.yaml``.

    Returns
    -------
    ProjectFile
        The checked contents of the file.

    Raises
    ------
    ProtocolError
        If the file is missing or not a YAML mapping, a top-level key is
        unknown, ``studies`` is missing or empty, a study folder has no
        ``study.yaml``, ``analyses`` or ``metadata`` is not a mapping, an
        analysis's ``studies:`` is not a non-empty list or names an unknown
        study.
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
        # A study finds its project in the folder above it, and the project
        # publishes only what is under it, so each study is directly inside.
        if folder.parent != file.parent.resolve():
            raise ProtocolError(
                f"{file}: study {label} is {folder}, not a folder directly inside the project.",
                hint=f"Move it to {file.parent / folder.name} and list it as "
                f"'{label}: {folder.name}'.",
            )
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
        if entry is not None and not isinstance(entry, Mapping):
            raise ProtocolError(
                f"{file}: analyses.{run} must be a mapping of settings, got {entry!r}.",
                hint=f"For example '{run}: {{selection: protein}}', or '{run}: {{}}' for the "
                "defaults.",
            )
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
    return ProjectFile(
        path=file,
        studies=studies,
        analyses=analyses,
        runs_in=runs_in,
        metadata=dict(metadata),
    )


def find_project(study_dir: Path) -> tuple[ProjectFile, str] | None:
    """Return the project that lists ``study_dir``, with the study's label.

    Only the folder directly above ``study_dir`` is searched for a
    ``project.yaml``.

    Parameters
    ----------
    study_dir : Path
        A study folder.

    Returns
    -------
    tuple of (ProjectFile, str) or None
        The project file and the study's label in it, or ``None`` when the
        folder above holds no ``project.yaml`` or that file does not list
        ``study_dir``.

    Raises
    ------
    ProtocolError
        If the ``project.yaml`` found cannot be read or fails its checks
        (see :func:`load_project_file`).
    """
    candidate = Path(study_dir).resolve().parent / PROJECT_FILE
    if not candidate.is_file():
        return None
    project = load_project_file(candidate)
    label = project.label_of(study_dir)
    return None if label is None else (project, label)


def resolve_names(
    value: Any,
    regions: Mapping[str, str],
    structures: Mapping[str, Path],
    where: str,
    selection: bool = False,
) -> Any:
    """Replace ``region <name>`` and ``structure <name>`` with a study's selection and file.

    ``region <name>`` in a selection becomes the region's selection in
    parentheses. A selection is a value under a key that contains
    ``selection`` (``selection_a``, ``selections``) or under one of
    :data:`_SELECTION_KEYS`; other strings, such as labels, keep the words
    as written. A string that is only ``structure <name>`` becomes the path
    of that structure, as a string. Mappings and lists are resolved item by
    item; other values are returned as they are.

    Parameters
    ----------
    value : Any
        A project analysis setting, or a mapping or list of them.
    regions : Mapping of str to str
        The study's region names mapped to their selections.
    structures : Mapping of str to Path
        The study's structure names mapped to their files.
    where : str
        Location used in error messages.
    selection : bool, optional
        Whether ``value`` is a selection itself.

    Returns
    -------
    Any
        ``value`` with every name replaced.

    Raises
    ------
    ProtocolError
        If a name is not one of the study's regions or structures, or a
        string starting with ``structure`` is not ``structure <name>`` with
        exactly one name.
    """
    if isinstance(value, Mapping):
        return {
            key: resolve_names(
                item,
                regions,
                structures,
                where,
                selection or "selection" in str(key) or key in _SELECTION_KEYS,
            )
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [resolve_names(item, regions, structures, where, selection) for item in value]
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
        """Return the selection of the matched region name, in parentheses."""
        name = match.group(1)
        if name not in regions:
            raise ProtocolError(
                f"{where} uses region {name}, which this study does not define.",
                hint="Add it under regions: in the study.yaml"
                + (f" (it has {', '.join(regions)})" if regions else "")
                + ", or limit the analysis with studies: in project.yaml.",
            )
        return f"({regions[name]})"

    return _REGION.sub(region, value) if selection else value
