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

_TOP_KEYS = (
    "polyzymd",
    "equilibration",
    "stride",
    "replicates",
    "conditions",
    "analyses",
    "metadata",
)
_ENTRY_KEYS = ("analysis",)


@dataclass(frozen=True)
class AnalysisEntry:
    """One entry of ``analyses:``: a shipped analysis and the settings it runs with.

    ``run`` is the entry's key, which names its results folder; ``analysis``
    is the shipped analysis it runs, the key itself unless given.
    """

    run: str
    analysis: str
    settings: dict[str, Any] = field(default_factory=dict)


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


def load_study_file(path: str | Path) -> StudyFile:
    """Read and check ``study.yaml``.

    Condition paths are resolved against the folder holding the file, so
    the folder can be moved as a whole. ``replicates`` takes a list or a
    range string such as ``"1-5"``.

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
    return StudyFile(
        path=file,
        equilibration=equilibration,
        conditions=conditions,
        analyses=analyses,
        polyzymd=None if version is None else str(version),
        stride=stride,
        replicates=replicates,
        metadata=dict(metadata),
    )
