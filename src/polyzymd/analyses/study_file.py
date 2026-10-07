"""Read ``study.yaml``, the analysis protocol of a study folder.

A study folder holds one MD study; see the "Study folders" explanation page.
``study.yaml`` names the conditions (each a simulation ``config.yaml``,
control first), the equilibration window applied to every replicate (an
analysis may set its own), and the settings of every analysis. :func:`load_study_file` reads and
checks it: every key must be one this module knows, and a misspelt key is
refused with the nearest known spelling, so a typo never falls back to a
default silently.
"""

from __future__ import annotations

import difflib
import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError

#: File name that marks a study folder.
STUDY_FILE = "study.yaml"
#: Folder, beside ``study.yaml``, that holds each analysis's results.
RESULTS_FOLDER = "results"
#: Where this machine keeps each condition's runs; never committed or published.
DATA_FILE = "data.local.yaml"

_TOP_KEYS = (
    "description",
    "equilibration",
    "stride",
    "until",
    "replicates",
    "conditions",
    "structures",
    "regions",
    "analyses",
    "metadata",
)
#: Keys of a condition written as a mapping.
_CONDITION_KEYS = ("config", "factors")
_ENTRY_KEYS = ("analysis",)
#: Keys of any ``analyses:`` entry that set its own analysis window.
WINDOW_KEYS = ("equilibration", "until", "stride")
#: Keys of an ``analyses:`` entry that runs your own function.
USER_KEYS = (
    "function",
    "kind",
    "unit",
    "selections",
    "universe",
    "settings",
    "labels",
    "reduce",
    "allow_empty",
    "missing",
    "parts",
)
#: How a user function is run: once per replicate, or once per frame.
USER_KINDS = ("per_replicate", "timeseries")


@dataclass(frozen=True)
class UserFunction:
    """An ``analyses:`` entry that runs a function from a Python file of the study.

    Attributes
    ----------
    file : Path
        The Python file, ``function:`` before the colon, relative to the study.
    qualname : str
        The function's name in that file, after the colon.
    kind : str
        ``"per_replicate"``: called once per replicate with ``frames=`` the
        production frame indices, returning one value (or one per label with
        ``labels: returned``). ``"timeseries"``: called once per production
        frame, returning one number, and each replicate's series reduced to
        one value by ``reduce``.
    unit : str or None
        Unit of the values, recorded and printed.
    selections : dict of str to str
        Keyword argument to MDAnalysis selection; each is passed as the
        replicate's ``AtomGroup``.
    universe : str or None
        Keyword argument that receives the replicate's ``Universe``.
    settings : dict
        Keyword arguments passed as they are.
    labels : str or None
        ``"returned"`` when the function returns ``(labels, values)``.
    reduce : str
        How a timeseries becomes one value per replicate, ``"mean"`` by default.
    allow_empty : bool
        Pass a selection that matches no atoms, such as a polymer selection
        in a no-polymer control, to the function as an empty AtomGroup, as
        the shipped analyses do; ``False`` refuses such a replicate.
    parts : list of str or None
        Names of several quantities the function measures in one pass: a
        timeseries function returns, per frame, a dict with these keys (or
        a sequence in this order), a per_replicate function one row per
        part. Each part is stored and reported as its own result.
    missing : float or None
        For ``labels: returned``, the value a replicate gets for a label that
        other replicates returned and it did not, such as ``nan``; ``None``
        refuses such a replicate.
    """

    file: Path
    qualname: str
    kind: str
    unit: str | None = None
    selections: dict[str, str] = field(default_factory=dict)
    universe: str | None = None
    settings: dict[str, Any] = field(default_factory=dict)
    labels: str | None = None
    reduce: str = "mean"
    allow_empty: bool = False
    missing: float | None = None
    parts: list[str] | None = None


@dataclass(frozen=True)
class AnalysisEntry:
    """One entry of ``analyses:``: an analysis and the settings it runs with.

    ``run`` is the entry's key, which names its results folder. ``analysis``
    is the shipped analysis it runs, the key itself unless given, or
    ``None`` when ``function`` names your own function instead.
    ``equilibration``, ``until`` and ``stride`` are the entry's own analysis
    window and frame stride, or ``None`` to use the study's.
    """

    run: str
    analysis: str | None
    settings: dict[str, Any] = field(default_factory=dict)
    function: UserFunction | None = None
    equilibration: str | None = None
    until: str | None = None
    stride: int | None = None


@dataclass(frozen=True)
class StudyFile:
    """The contents of one ``study.yaml``, with condition paths made absolute."""

    path: Path
    equilibration: str
    conditions: dict[str, Path]
    analyses: dict[str, AnalysisEntry]
    stride: int = 1
    replicates: list[int] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    data: dict[str, Path] = field(default_factory=dict)
    until: str | None = None
    #: What the study simulates, such as the protein and temperature.
    description: str | None = None
    #: Named structure files, used as ``structure <name>``.
    structures: dict[str, Path] = field(default_factory=dict)
    #: Named selections, used as ``region <name>``.
    regions: dict[str, str] = field(default_factory=dict)
    #: Each condition's factors, such as ``{"sbma_fraction": 0.5}``.
    factors: dict[str, dict[str, Any]] = field(default_factory=dict)
    #: The project this study belongs to, and its label there.
    project: Any = None
    project_label: str | None = None

    @property
    def root(self) -> Path:
        """The study folder: the directory holding ``study.yaml``."""
        return self.path.parent

    def results_dir(self, run: str) -> Path:
        """The folder of one analysis run's stored results, report and figures."""
        return self.root / RESULTS_FOLDER / run

    def window(self, run: str) -> tuple[str, str | None]:
        """Return the equilibration and ``until`` of ``run``: its entry's, else the study's.

        A run the file does not list runs with the study's window.
        """
        entry = self.analyses.get(run)
        if entry is None:
            return self.equilibration, self.until
        return (
            entry.equilibration or self.equilibration,
            entry.until if entry.until is not None else self.until,
        )

    def stride_of(self, run: str) -> int:
        """Return the frame stride of ``run``: its entry's, else the study's."""
        entry = self.analyses.get(run)
        return entry.stride if entry is not None and entry.stride is not None else self.stride


def find_study_file(path: str | Path) -> Path:
    """Return the ``study.yaml`` that ``path`` names: the file itself, or the one in that folder."""
    candidate = Path(path).expanduser().resolve()
    if candidate.is_dir():
        candidate = candidate / STUDY_FILE
    if not candidate.is_file():
        raise ProtocolError(
            f"No study file at {candidate}.",
            hint=f"Give the path of a {STUDY_FILE}, or of the study folder that holds one.",
        )
    return candidate


def check_folder_names(labels: Any, where: str, study_folders: bool = True) -> None:
    """Raise ProtocolError when two condition labels give one folder name.

    Results are stored in one folder per label, so ``SBMA 50`` and
    ``SBMA 50%`` would share one. With ``study_folders`` (a study.yaml), the
    lower-case folder names of ``conditions/`` and of the deposit are
    checked too, so ``WT`` and ``wt`` are refused there; conditions made in
    code (``Study.from_configs``) have no such folders.
    """
    from polyzymd.analyses.study_scaffold import condition_folder
    from polyzymd.analyses.timeseries import _safe

    for name in (_safe, condition_folder) if study_folders else (_safe,):
        seen: dict[str, str] = {}
        for label in map(str, labels):
            other = seen.setdefault(name(label), label)
            if other != label:
                raise ProtocolError(
                    f"{where}: the conditions {other!r} and {label!r} give one folder name, "
                    f"{name(label)!r}, so their results would overwrite each other.",
                    hint="Rename one of them.",
                )


def _unknown(keys: Any, known: tuple[str, ...], where: str) -> None:
    """Refuse every key of ``keys`` not in ``known``, naming the closest known key."""
    for key in keys:
        if key in known:
            continue
        close = difflib.get_close_matches(str(key), known, n=1)
        raise ProtocolError(
            f"{where} has an unknown key {key!r}.",
            hint=(f"Did you mean {close[0]!r}? " if close else "")
            + f"The keys it takes are {', '.join(known)}.",
        )


def _analysis_entry(run: str, raw: Any, path: Path) -> AnalysisEntry:
    """Check one ``analyses:`` entry against the settings of the analysis it names."""
    from polyzymd.analyses.protocols import ANALYSES

    if raw is None:
        raw = {}
    if not isinstance(raw, Mapping):
        raise ProtocolError(
            f"{path}: analyses.{run} must be a mapping of settings, got {raw!r}.",
            hint=f"Write it as '{run}: {{setting: value}}', or '{run}: {{}}' for the defaults.",
        )
    raw = dict(raw)
    window = {key: raw.pop(key) for key in WINDOW_KEYS if key in raw}
    equilibration = (
        None
        if window.get("equilibration") is None
        else _equilibration(window["equilibration"], f"{path}: analyses.{run}")
    )
    until = _until(window.get("until"), f"{path}: analyses.{run}")
    stride = window.get("stride")
    if stride is not None and (
        isinstance(stride, bool) or not isinstance(stride, int) or stride < 1
    ):
        raise ProtocolError(
            f"{path}: analyses.{run}: stride must be a whole number of at least 1, got {stride!r}.",
            hint="Leave it out to use the study's stride.",
        )
    own = {"equilibration": equilibration, "until": until, "stride": stride}
    if "function" in raw:
        return AnalysisEntry(run, None, {}, _user_function(run, raw, path), **own)
    settings = dict(raw)
    analysis = str(settings.pop("analysis", run))
    if analysis not in ANALYSES:
        close = difflib.get_close_matches(analysis, list(ANALYSES), n=1)
        raise ProtocolError(
            f"{path}: analyses.{run} names the analysis {analysis!r}, which PolyzyMD does not ship.",
            hint=(f"Did you mean {close[0]!r}? " if close else "")
            + f"Use one of {', '.join(ANALYSES)}"
            + ("" if "analysis" in raw else f", or add 'analysis: NAME' to the {run!r} entry")
            + ".",
        )
    _unknown(
        settings,
        (*_ENTRY_KEYS, *WINDOW_KEYS, *ANALYSES[analysis].defaults),
        f"{path}: analyses.{run}",
    )
    # A relative file, such as reference_file: structures/ref.pdb, is relative
    # to the file that lists the analysis, whatever the shell's folder.
    return AnalysisEntry(run, analysis, _resolve_files(settings, path.parent), **own)


def _user_function(run: str, raw: Mapping, path: Path) -> UserFunction:
    """Check an ``analyses:`` entry that names a function in a Python file of the study."""
    where = f"{path}: analyses.{run}"
    _unknown(raw, USER_KEYS, where)
    spec = str(raw["function"])
    file, colon, qualname = spec.rpartition(":")
    if not colon or not file or not qualname:
        raise ProtocolError(
            f"{where}: function must be 'path/to/file.py:function_name', got {spec!r}.",
            hint="For example 'function: analyses/lid.py:lid_distance', relative to study.yaml.",
        )
    location = (path.parent / file).resolve()
    if not location.is_file():
        raise ProtocolError(
            f"{where}: no Python file at {location}.",
            hint="Give the file relative to study.yaml, for example analyses/lid.py.",
        )
    kind = raw.get("kind")
    if kind not in USER_KINDS:
        raise ProtocolError(
            f"{where}: kind must be one of {', '.join(USER_KINDS)}, got {kind!r}.",
            hint="per_replicate calls the function once per replicate with frames=; "
            "timeseries calls it once per frame and averages each replicate's series.",
        )
    for key in ("selections", "settings"):
        if not isinstance(raw.get(key, {}), Mapping):
            raise ProtocolError(
                f"{where}: {key} must be a mapping of keyword argument to value.",
                hint="For example 'selections: {lid: \"resid 140-150 and name CA\"}'.",
            )
    labels = raw.get("labels")
    if labels not in (None, "returned") or (labels and kind != "per_replicate"):
        raise ProtocolError(
            f"{where}: labels can only be 'returned', for a per_replicate function.",
            hint="Return (labels, values) from the function and write 'labels: returned'.",
        )
    reduce = str(raw.get("reduce", "mean"))
    if reduce != "mean" and kind != "timeseries":
        raise ProtocolError(
            f"{where}: reduce applies only to a timeseries function.",
            hint="Leave reduce out for a per_replicate function.",
        )
    clash = set(raw.get("selections") or {}) & set(raw.get("settings") or {})
    if clash or (
        raw.get("universe") and raw["universe"] in clash | set(raw.get("selections") or {})
    ):
        raise ProtocolError(
            f"{where}: a keyword argument is given twice: {sorted(clash) or raw['universe']}.",
            hint="Give each keyword argument in only one of selections, settings and universe.",
        )
    missing = raw.get("missing")
    if missing is not None:
        if labels != "returned":
            raise ProtocolError(
                f"{where}: missing applies only to a per_replicate function with labels: returned.",
                hint="Leave missing out, or return (labels, values) and write labels: returned.",
            )
        if isinstance(missing, bool) or not isinstance(missing, (int, float)):
            raise ProtocolError(
                f"{where}: missing must be a number, not {missing!r}.",
                hint="Write a number, or .nan for a label without a value.",
            )
    parts = raw.get("parts")
    if parts is not None:
        if not isinstance(parts, list) or not parts or len(set(map(str, parts))) != len(parts):
            raise ProtocolError(
                f"{where}: parts must be a list of distinct names, not {parts!r}.",
                hint="For example 'parts: [area, contacts, gyration]'.",
            )
        parts = [str(part) for part in parts]
    unit = raw.get("unit")
    return UserFunction(
        file=location,
        qualname=qualname,
        kind=kind,
        unit=None if unit is None else str(unit),
        selections={str(k): str(v) for k, v in (raw.get("selections") or {}).items()},
        universe=None if raw.get("universe") is None else str(raw["universe"]),
        settings=_resolve_files(dict(raw.get("settings") or {}), path.parent),
        labels=labels,
        reduce=reduce,
        allow_empty=bool(raw.get("allow_empty", False)),
        missing=None if missing is None else float(missing),
        parts=parts,
    )


def read_data_file(path: Path, labels: Any) -> dict[str, Path]:
    """Read ``data.local.yaml``: condition label to the directory holding its runs.

    Relative directories are taken from the study folder. A missing file
    gives an empty mapping, so each condition's config says where its runs
    are, as it does for someone running the simulations.

    Raises
    ------
    ProtocolError
        If the file is not a mapping of known condition labels to paths.
    """
    import yaml

    if not path.is_file():
        return {}
    try:
        raw = yaml.safe_load(path.read_text()) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ProtocolError(
            f"Cannot read {path}: {exc}", hint="Check that it is valid YAML."
        ) from exc
    if not isinstance(raw, Mapping):
        raise ProtocolError(
            f"{path} must map each condition label to the directory holding its runs.",
            hint="Write it with polyzymd study locate DIR, or see data.example.yaml.",
        )
    known = tuple(labels)
    _unknown(raw, known, str(path))
    return {
        str(label): (path.parent / Path(str(folder)).expanduser()).resolve()
        for label, folder in raw.items()
    }


#: Suffixes that mark a setting as naming a file, for ``study check``.
FILE_SUFFIXES = (
    ".pdb",
    ".gro",
    ".tpr",
    ".top",
    ".itp",
    ".xtc",
    ".dcd",
    ".sdf",
    ".mol2",
    ".json",
    ".yaml",
    ".yml",
    ".csv",
    ".tsv",
    ".txt",
    ".npy",
    ".npz",
    ".dat",
)


def _resolve_files(value: Any, folder: Path) -> Any:
    """Make a relative path that names an existing file or folder under ``folder`` absolute.

    A study's own function receives its ``settings`` as written, and runs
    with the shell's working directory, so ``structures/ref.pdb`` is resolved
    here against the folder of the file that lists it, as the shipped
    analyses resolve theirs; other values are returned as they are.
    """
    if isinstance(value, Mapping):
        return {key: _resolve_files(item, folder) for key, item in value.items()}
    if isinstance(value, list):
        return [_resolve_files(item, folder) for item in value]
    if (
        isinstance(value, str)
        and value
        and "\n" not in value
        and len(value) < 1024
        and not Path(value).is_absolute()
        and (folder / value).exists()
    ):
        return str((folder / value).resolve())
    return value


def missing_files(value: Any) -> list[str]:
    """Return the settings values that look like file names but name no existing file.

    A string ending in one of :data:`FILE_SUFFIXES` that is not an existing
    path after :func:`_resolve_files`; ``study check`` reports them before
    any analysis runs.
    """
    if isinstance(value, Mapping):
        return [m for item in value.values() for m in missing_files(item)]
    if isinstance(value, list):
        return [m for item in value for m in missing_files(item)]
    if (
        isinstance(value, str)
        and value.lower().endswith(FILE_SUFFIXES)
        and not Path(value).exists()
    ):
        return [value]
    return []


def entry_record(protocol: StudyFile, run: str) -> dict[str, Any]:
    """Return everything an ``analyses:`` entry says, resolved for the study, fit to publish.

    The shipped analysis or the function (file relative to the study, name,
    kind, unit, selections, settings, labels, reduce, allow_empty, missing,
    parts), the window and the stride. A report records it, so ``freeze``
    can tell when any of it changed since.
    """
    entry = protocol.analyses[run]
    project_root = protocol.project.root if protocol.project is not None else None
    equilibration, until = protocol.window(run)
    record: dict[str, Any] = {
        "analysis": entry.analysis,
        "settings": portable(entry.settings, protocol.root, project_root),
        "equilibration": equilibration,
        "until": until,
        "stride": protocol.stride_of(run),
    }
    function = entry.function
    if function is not None:
        record["function"] = {
            "file": portable(str(function.file), protocol.root, project_root),
            "qualname": function.qualname,
            "kind": function.kind,
            "unit": function.unit,
            "selections": dict(function.selections),
            "universe": function.universe,
            "settings": portable(function.settings, protocol.root, project_root),
            "labels": function.labels,
            "reduce": function.reduce,
            "allow_empty": function.allow_empty,
            "missing": function.missing,
            "parts": function.parts,
        }
    return json.loads(json.dumps(record, default=str))


def portable(value: Any, root: Path, project_root: Path | None = None) -> Any:
    """Return ``value`` with absolute paths made fit to publish, for reports and manifests.

    A path inside the study folder ``root`` or its project folder becomes a
    path relative to ``root`` (``structures/ref.pdb``, ``../analyses/f.py``);
    any other absolute path, which names a place on one machine, becomes its
    file name. Mappings and lists are converted item by item.
    """
    import os

    if isinstance(value, Mapping):
        return {key: portable(item, root, project_root) for key, item in value.items()}
    if isinstance(value, list):
        return [portable(item, root, project_root) for item in value]
    if isinstance(value, str) and Path(value).is_absolute():
        path = Path(value)
        for base in (root, project_root):
            if base is not None and path.is_relative_to(Path(base).resolve()):
                return Path(os.path.relpath(path, Path(root).resolve())).as_posix()
        return path.name
    return value


def outside_configs(protocol: StudyFile) -> dict[str, Path]:
    """Return the condition configs outside the study folder and its project folder, by label.

    Freeze deposits only those folders, so it refuses these configs.
    """
    inside = [protocol.root.resolve()]
    if protocol.project is not None:
        inside.append(protocol.project.root.resolve())
    return {
        label: config
        for label, config in protocol.conditions.items()
        if not any(config.is_relative_to(folder) for folder in inside)
    }


def _condition(label: str, value: Any, file: Path) -> tuple[Path, dict[str, Any]]:
    """Read one condition: a config path (or its folder), or ``{config: ..., factors: {...}}``."""
    factors: dict[str, Any] = {}
    if isinstance(value, Mapping):
        _unknown(value, _CONDITION_KEYS, f"{file}: conditions.{label}")
        if "config" not in value:
            raise ProtocolError(
                f"{file}: conditions.{label} has no config.",
                hint="Write it as '{config: conditions/<name>, factors: {name: value}}'.",
            )
        factors = value.get("factors") or {}
        if not isinstance(factors, Mapping) or any(
            isinstance(v, (Mapping, list)) for v in factors.values()
        ):
            raise ProtocolError(
                f"{file}: conditions.{label}.factors must map names to single values.",
                hint="For example 'factors: {sbma_fraction: 0.5}'.",
            )
        factors = {str(k): v for k, v in factors.items()}
        value = value["config"]
    path = (file.parent / Path(str(value)).expanduser()).resolve()
    if path.is_dir():
        path = path / "config.yaml"
    return path, factors


def _named(value: Any, key: str, file: Path) -> dict[str, Any]:
    """Read ``structures:`` or ``regions:``: a mapping of names to values."""
    if value is None:
        return {}
    if not isinstance(value, Mapping) or any(
        not isinstance(name, str) or not name.replace("_", "a").replace("-", "a").isalnum()
        for name in value
    ):
        raise ProtocolError(
            f"{file}: {key} must map names (letters, digits, _ and -) to values.",
            hint=f"For example '{key}: {{core: resid 5-120}}'."
            if key == "regions"
            else f"For example '{key}: {{reference: structures/ref.pdb}}'.",
        )
    return dict(value)


def _equilibration(value: Any, where: Any) -> str:
    """Check an equilibration window, such as ``100ns``."""
    from polyzymd.analyses.shared.loader import parse_time_string

    text = str(value)
    try:
        parse_time_string(text)
    except ValueError as exc:
        raise ProtocolError(
            f"{where}: cannot read equilibration {text!r}: {exc}",
            hint="Write it as a time such as '100ns', '500ps' or '0ns'.",
        ) from exc
    return text


def _until(value: Any, file: Any) -> str | None:
    """Check ``until:``, the end of a common analysis window, or ``common``."""
    if value is None:
        return None
    if str(value) == "common":
        return "common"
    from polyzymd.analyses.shared.loader import parse_time_string

    try:
        parse_time_string(str(value))
    except ValueError as exc:
        raise ProtocolError(
            f"{file}: cannot read until {value!r}: {exc}",
            hint="Write it as a time such as '38ns', 'common' to end every replicate at the "
            "shortest one's last time, or leave it out to use every production frame.",
        ) from exc
    return str(value)


def load_study_file(path: str | Path) -> StudyFile:
    """Read and check ``study.yaml``.

    Condition paths are resolved against the folder holding the file, so
    the folder can be moved as a whole. ``replicates`` takes a list or a
    range string such as ``"1-5"``. ``data`` holds ``data.local.yaml`` beside
    the file, when there is one (:func:`read_data_file`).

    Raises
    ------
    ProtocolError
        If the file is missing or not YAML, a key is unknown, a required key
        (``equilibration``, ``conditions``) is missing, a value has the
        wrong type, or two condition labels give one folder name.
    """
    import yaml

    file = find_study_file(path)
    try:
        raw = yaml.safe_load(file.read_text()) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ProtocolError(
            f"Cannot read {file}: {exc}", hint="Check that the file is valid YAML."
        ) from exc
    if not isinstance(raw, Mapping):
        raise ProtocolError(
            f"{file} must be a YAML mapping.", hint="Start from the documented example."
        )
    _unknown(raw, _TOP_KEYS, str(file))
    for required in ("equilibration", "conditions"):
        if required not in raw:
            raise ProtocolError(
                f"{file} has no {required!r}.",
                hint="A study.yaml names its conditions (control first) and one equilibration "
                "window, for example 'equilibration: 100ns'.",
            )

    equilibration = _equilibration(raw["equilibration"], file)

    conditions_raw = raw["conditions"]
    if not isinstance(conditions_raw, Mapping) or not conditions_raw:
        raise ProtocolError(
            f"{file}: conditions must map each condition label to its config.yaml.",
            hint="For example 'conditions: {No polymer: conditions/no_polymer/config.yaml}'.",
        )
    conditions: dict[str, Path] = {}
    factors: dict[str, dict[str, Any]] = {}
    for label, value in conditions_raw.items():
        conditions[str(label)], factors[str(label)] = _condition(str(label), value, file)
    check_folder_names(conditions, str(file))

    structures = _named(raw.get("structures"), "structures", file)
    structures = {
        name: (file.parent / Path(str(where)).expanduser()).resolve()
        for name, where in structures.items()
    }
    for name, where in structures.items():
        if not where.is_file():
            raise ProtocolError(
                f"{file}: structure {name} is not a file: {where}.",
                hint="Give each structure's path relative to study.yaml, such as structures/ref.pdb.",
            )
    regions = {
        name: str(text) for name, text in _named(raw.get("regions"), "regions", file).items()
    }

    stride = raw.get("stride", 1)
    if isinstance(stride, bool) or not isinstance(stride, int) or stride < 1:
        raise ProtocolError(
            f"{file}: stride must be a whole number of at least 1, got {stride!r}.",
            hint="Leave it out to measure every production frame.",
        )

    replicates = raw.get("replicates")
    if replicates is not None:
        from polyzymd.utils.replicates import parse_replicate_range

        try:
            replicates = (
                parse_replicate_range(replicates)
                if isinstance(replicates, str)
                else [int(index) for index in replicates]
            )
        except (TypeError, ValueError) as exc:
            raise ProtocolError(
                f"{file}: cannot read replicates {raw['replicates']!r}: {exc}",
                hint="Write a list such as [1, 2, 3] or a range such as '1-5'.",
            ) from exc

    analyses_raw = raw.get("analyses") or {}
    if not isinstance(analyses_raw, Mapping):
        raise ProtocolError(
            f"{file}: analyses must map each run name to its settings.",
            hint="For example 'analyses: {contacts: {method: occlusion}}'.",
        )
    from polyzymd.analyses.project_file import find_project, resolve_names

    found = find_project(file.parent)
    project, project_label = found if found is not None else (None, None)
    analyses: dict[str, AnalysisEntry] = {}
    if project is not None:
        # The project's analyses, for this study: its regions and structures resolved.
        for run, entry in project.analyses.items():
            listed = project.runs_in[run]
            if listed is not None and project_label not in listed:
                continue
            where = f"{project.path}: analyses.{run} (in study {project_label})"
            resolved = resolve_names(entry, regions, structures, where)
            analyses[run] = _analysis_entry(run, resolved, project.path)
    for run, entry in analyses_raw.items():
        if str(run) in analyses:
            raise ProtocolError(
                f"{file}: analyses.{run} is also an analysis of {project.path}.",
                hint="Give the study's own analysis another name, or change it in project.yaml.",
            )
        resolved = resolve_names(entry, regions, structures, f"{file}: analyses.{run}")
        analyses[str(run)] = _analysis_entry(str(run), resolved, file)

    metadata = raw.get("metadata") or {}
    if not isinstance(metadata, Mapping):
        raise ProtocolError(f"{file}: metadata must be a mapping.", hint="Leave it out for now.")

    data = read_data_file(file.parent / DATA_FILE, conditions)
    return StudyFile(
        path=file,
        equilibration=equilibration,
        conditions=conditions,
        analyses=analyses,
        stride=stride,
        replicates=replicates,
        # A study of a project publishes with the project's metadata unless it has its own.
        metadata=dict(metadata) or (dict(project.metadata) if project is not None else {}),
        data=data,
        until=_until(raw.get("until"), file),
        description=None if raw.get("description") is None else str(raw["description"]),
        structures=structures,
        regions=regions,
        factors={label: value for label, value in factors.items() if value},
        project=project,
        project_label=project_label,
    )
