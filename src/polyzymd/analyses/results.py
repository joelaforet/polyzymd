"""Every number a study's comparisons report, as three tidy tables.

:func:`load_results` reads the comparison results the runner wrote and returns
one row per reported quantity, so a figure script or an agent works from a
table instead of from artifact JSON::

    from polyzymd.analyses import load_results

    results = load_results("paper1_lipase_thermal")        # a study root
    frame = results.to_dataframe("conditions")               # needs pandas
    results.to_csv("paper1_lipase_thermal/results")

``conditions`` holds one row per condition and scalar observable: mean, SEM,
95 percent interval over replicates, the replicate values and their count.
``comparisons`` holds one row per pairwise test, with the adjusted p-value.
``profiles`` holds one row per index of a per-residue or per-bin observable.
Only the ``comparison/<analysis>/result.json`` files are read, so a published
study works without its trajectories.
"""

from __future__ import annotations

import csv
import glob
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import yaml

TABLES = ("conditions", "comparisons", "profiles")

_CONDITION_COLUMNS = (
    "study",
    "comparison",
    "analysis",
    "condition",
    "is_control",
    "observable",
    "kind",
    "unit",
    "n_replicates",
    "mean",
    "sem",
    "ci95_low",
    "ci95_high",
    "coverage",
    "n_eff_min",
    "replicate_values",
    "equilibration",
)
_COMPARISON_COLUMNS = (
    "study",
    "comparison",
    "analysis",
    "observable",
    "kind",
    "unit",
    "control",
    "condition",
    "n_control",
    "n_condition",
    "delta",
    "percent_change",
    "test",
    "p_value",
    "p_adjusted",
    "correction",
    "cohens_d",
    "significant",
    "testable",
    "note",
)
_PROFILE_COLUMNS = (
    "study",
    "comparison",
    "analysis",
    "condition",
    "observable",
    "unit",
    "index_label",
    "index",
    "mean",
    "sem",
    "n_replicates",
)


@dataclass
class Results:
    """The rows of every comparison read, and what could not be read."""

    conditions: list[dict[str, Any]] = field(default_factory=list)
    comparisons: list[dict[str, Any]] = field(default_factory=list)
    profiles: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dataframe(self, table: str = "conditions") -> Any:
        """One table as a pandas DataFrame. Requires pandas."""
        import pandas as pd

        return pd.DataFrame(self._rows(table), columns=_columns(table))

    def to_csv(self, directory: Path | str) -> list[Path]:
        """Write ``conditions.csv``, ``comparisons.csv`` and ``profiles.csv``."""
        target = Path(directory)
        target.mkdir(parents=True, exist_ok=True)
        written = []
        for table in TABLES:
            path = target / f"{table}.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=_columns(table))
                writer.writeheader()
                for row in self._rows(table):
                    writer.writerow({key: _csv_value(value) for key, value in row.items()})
            written.append(path)
        return written

    def _rows(self, table: str) -> list[dict[str, Any]]:
        if table not in TABLES:
            raise ValueError(f"table must be one of {TABLES}, got {table!r}")
        return getattr(self, table)


def load_results(
    sources: Path | str | Sequence[Path | str],
    *,
    analyses: Iterable[str] | None = None,
) -> Results:
    """Read the results of one or more comparisons into tidy rows.

    Parameters
    ----------
    sources : path, glob or sequence of them
        A study root (every ``comparisons/*/comparison.yaml`` below it), a
        comparison folder, a ``comparison.yaml``, or a glob matching any of
        these.
    analyses : iterable of str, optional
        Only these analyses. Defaults to every analysis with results.

    Returns
    -------
    Results
        The three tables. A comparison with no results on disk is named in
        ``warnings`` rather than raising.
    """
    from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError

    wanted = set(analyses) if analyses is not None else None
    results = Results()
    for directory in _comparison_dirs(sources):
        study, name = _study_name(directory), _comparison_name(directory)
        found = sorted((directory / "comparison").glob("*/result.json"))
        if wanted is not None:
            found = [path for path in found if path.parent.name in wanted]
        if not found:
            results.warnings.append(f"{directory}: no comparison results on disk")
            continue
        for path in found:
            try:
                artifact = ArtifactStore(path.parent).read_comparison_result("result.json")
            except ArtifactStoreError as exc:
                results.warnings.append(f"{path}: unreadable comparison result: {exc}")
                continue
            base = {"study": study, "comparison": name, "analysis": artifact.analysis_name}
            _add_rows(results, base, artifact)
    return results


def _add_rows(results: Results, base: dict[str, Any], artifact: Any) -> None:
    """Append one comparison artifact's rows to the three tables."""
    control = artifact.control_label
    equilibration = artifact.metadata.get("equilibration")
    for label, aggregates in (artifact.payload.get("conditions") or {}).items():
        for item in aggregates:
            if item.get("kind") == "profile":
                for position, index in enumerate(item.get("index") or []):
                    results.profiles.append(
                        {
                            **base,
                            "condition": label,
                            "observable": item["name"],
                            "unit": item.get("unit"),
                            "index_label": item.get("index_label"),
                            "index": index,
                            "mean": _at(item.get("profile_mean"), position),
                            "sem": _at(item.get("profile_sem"), position),
                            "n_replicates": item.get("n_replicates"),
                        }
                    )
                continue
            results.conditions.append(
                {
                    **base,
                    "condition": label,
                    "is_control": label == control,
                    "observable": item["name"],
                    "kind": item.get("kind"),
                    "unit": item.get("unit"),
                    "n_replicates": item.get("n_replicates"),
                    "mean": item.get("mean"),
                    "sem": item.get("sem"),
                    "ci95_low": item.get("ci95_low"),
                    "ci95_high": item.get("ci95_high"),
                    "coverage": item.get("coverage"),
                    "n_eff_min": item.get("n_eff_min"),
                    "replicate_values": list(item.get("replicate_values") or []),
                    "equilibration": equilibration,
                }
            )
    for item in artifact.payload.get("comparisons") or []:
        row = {**base, **{key: item.get(key) for key in _COMPARISON_COLUMNS if key not in base}}
        row["observable"] = item.get("name")
        results.comparisons.append(row)


def _comparison_dirs(sources: Path | str | Sequence[Path | str]) -> list[Path]:
    """Comparison folders named by the sources, in a stable order."""
    items = [sources] if isinstance(sources, (str, Path)) else list(sources)
    directories: list[Path] = []
    for item in items:
        matches = [Path(match) for match in glob.glob(str(item))] or [Path(item)]
        for path in matches:
            if path.is_file():
                directories.append(path.parent)
            elif (path / "comparisons").is_dir():
                directories += [
                    yaml_path.parent for yaml_path in path.glob("comparisons/*/comparison.yaml")
                ]
            elif path.is_dir():
                directories.append(path)
    unique = sorted({directory.resolve() for directory in directories})
    return unique


def _comparison_name(directory: Path) -> str:
    """The ``name`` in a comparison's YAML, or its folder name."""
    path = directory / "comparison.yaml"
    if path.is_file():
        try:
            data = yaml.safe_load(path.read_text()) or {}
        except yaml.YAMLError:
            data = {}
        if isinstance(data, dict) and data.get("name"):
            return str(data["name"])
    return directory.name


def _study_name(directory: Path) -> str | None:
    """Name of the study a comparison belongs to, if it is inside one."""
    from polyzymd.config.study import STUDY_FILE, StudyConfig, find_study_root

    root = find_study_root(directory)
    if root is None:
        return None
    try:
        return StudyConfig.from_yaml(root / STUDY_FILE).name
    except (OSError, ValueError):
        return root.name


def _columns(table: str) -> tuple[str, ...]:
    return {
        "conditions": _CONDITION_COLUMNS,
        "comparisons": _COMPARISON_COLUMNS,
        "profiles": _PROFILE_COLUMNS,
    }[table]


def _at(values: Sequence[Any] | None, position: int) -> Any:
    return values[position] if values is not None and position < len(values) else None


def _csv_value(value: Any) -> Any:
    """Lists become JSON so a CSV cell round-trips; ``None`` becomes empty."""
    if isinstance(value, (list, tuple)):
        return json.dumps(list(value))
    return "" if value is None else value
