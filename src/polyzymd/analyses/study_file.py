"""Read ``study.yaml``, the analysis protocol of a study folder.

A study folder holds one MD study; see the "Study folders" explanation page.
``study.yaml`` names the conditions (each a simulation ``config.yaml``,
control first), the one equilibration window applied to every replicate,
and the settings of every analysis. :func:`load_study_file` reads and
checks it: every key must be one this module knows, and a misspelt key is
refused with the nearest known spelling, so a typo never falls back to a
default silently.
"""

from __future__ import annotations

import difflib
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
    "polyzymd",
    "equilibration",
    "stride",
    "until",
    "replicates",
    "conditions",
    "analyses",
    "metadata",
)
_ENTRY_KEYS = ("analysis",)
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


@dataclass(frozen=True)
class AnalysisEntry:
    """One entry of ``analyses:``: an analysis and the settings it runs with.

    ``run`` is the entry's key, which names its results folder. ``analysis``
    is the shipped analysis it runs, the key itself unless given, or
    ``None`` when ``function`` names your own function instead.
    """

    run: str
    analysis: str | None
    settings: dict[str, Any] = field(default_factory=dict)
    function: UserFunction | None = None


@dataclass(frozen=True)
class StudyFile:
    """The contents of one ``study.yaml``, with condition paths made absolute."""

    path: Path
    equilibration: str
    conditions: dict[str, Path]
    analyses: dict[str, AnalysisEntry]
    polyzymd: str | None = None
    stride: int = 1
    replicates: list[int] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    data: dict[str, Path] = field(default_factory=dict)
    until: str | None = None

    @property
    def root(self) -> Path:
        """The study folder: the directory holding ``study.yaml``."""
        return self.path.parent

    def results_dir(self, run: str) -> Path:
        """The folder of one analysis run's stored results, report and figures."""
        return self.root / RESULTS_FOLDER / run


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
    from polyzymd.analyses.protocols import FUNCTION_ANALYSES

    if raw is None:
        raw = {}
    if not isinstance(raw, Mapping):
        raise ProtocolError(
            f"{path}: analyses.{run} must be a mapping of settings, got {raw!r}.",
            hint=f"Write it as '{run}: {{setting: value}}', or '{run}: {{}}' for the defaults.",
        )
    if "function" in raw:
        return AnalysisEntry(run, None, {}, _user_function(run, raw, path))
    settings = dict(raw)
    analysis = str(settings.pop("analysis", run))
    if analysis not in FUNCTION_ANALYSES:
        close = difflib.get_close_matches(analysis, list(FUNCTION_ANALYSES), n=1)
        raise ProtocolError(
            f"{path}: analyses.{run} names the analysis {analysis!r}, which PolyzyMD does not ship.",
            hint=(f"Did you mean {close[0]!r}? " if close else "")
            + f"Use one of {', '.join(FUNCTION_ANALYSES)}"
            + ("" if "analysis" in raw else f", or add 'analysis: NAME' to the {run!r} entry")
            + ".",
        )
    _unknown(settings, (*_ENTRY_KEYS, *FUNCTION_ANALYSES[analysis]), f"{path}: analyses.{run}")
    return AnalysisEntry(run, analysis, settings)


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
    unit = raw.get("unit")
    return UserFunction(
        file=location,
        qualname=qualname,
        kind=kind,
        unit=None if unit is None else str(unit),
        selections={str(k): str(v) for k, v in (raw.get("selections") or {}).items()},
        universe=None if raw.get("universe") is None else str(raw["universe"]),
        settings=dict(raw.get("settings") or {}),
        labels=labels,
        reduce=reduce,
        allow_empty=bool(raw.get("allow_empty", False)),
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


def _until(value: Any, file: Path) -> str | None:
    """Check ``until:``, the end of a common analysis window."""
    if value is None:
        return None
    from polyzymd.analyses.shared.loader import parse_time_string

    try:
        parse_time_string(str(value))
    except ValueError as exc:
        raise ProtocolError(
            f"{file}: cannot read until {value!r}: {exc}",
            hint="Write it as a time such as '38ns', or leave it out to use every production frame.",
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
        (``equilibration``, ``conditions``) is missing, or a value has the
        wrong type.
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

    from polyzymd.analyses.shared.loader import parse_time_string

    equilibration = str(raw["equilibration"])
    try:
        parse_time_string(equilibration)
    except ValueError as exc:
        raise ProtocolError(
            f"{file}: cannot read equilibration {equilibration!r}: {exc}",
            hint="Write it as a time such as '100ns', '500ps' or '0ns'.",
        ) from exc

    conditions_raw = raw["conditions"]
    if not isinstance(conditions_raw, Mapping) or not conditions_raw:
        raise ProtocolError(
            f"{file}: conditions must map each condition label to its config.yaml.",
            hint="For example 'conditions: {No polymer: conditions/no_polymer/config.yaml}'.",
        )
    conditions = {
        str(label): (file.parent / Path(str(config)).expanduser()).resolve()
        for label, config in conditions_raw.items()
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
    analyses = {
        str(run): _analysis_entry(str(run), entry, file) for run, entry in analyses_raw.items()
    }

    metadata = raw.get("metadata") or {}
    if not isinstance(metadata, Mapping):
        raise ProtocolError(f"{file}: metadata must be a mapping.", hint="Leave it out for now.")

    version = raw.get("polyzymd")
    data = read_data_file(file.parent / DATA_FILE, conditions)
    return StudyFile(
        path=file,
        equilibration=equilibration,
        conditions=conditions,
        analyses=analyses,
        polyzymd=None if version is None else str(version),
        stride=stride,
        replicates=replicates,
        metadata=dict(metadata),
        data=data,
        until=_until(raw.get("until"), file),
    )
