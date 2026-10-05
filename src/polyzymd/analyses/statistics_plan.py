"""Statistics on stored results: replicate tables, trend tests and a study's own plan.

Slice P2 of the "Projects and studies" design:

- :func:`replicate_table` gives one row per replicate of an analysis run
  (and per label or part), the sampling unit of every test, from stored
  results only.
- :func:`trend_tests` fits, for each numeric factor the conditions declare,
  the slope of the replicate values against it, with the replicate as the
  unit, and corrects the family of factors with Benjamini-Hochberg.
- :func:`run_stats_plan` runs the function a ``stats:`` entry names on a
  :class:`~polyzymd.analyses.study.Study` or
  :class:`~polyzymd.analyses.project.Project` and stores what it returns with
  a record of its code and of the reports it read.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError

#: Folder, under a study's or project's ``results/``, of each stats plan's output.
STATS_FOLDER = "stats"
#: File holding a stats plan's record: its code and the reports it read.
STATS_RECORD = "record.json"
_REDUCE = {"mean": "mean", "fraction": "mean", "std": "std"}


@dataclass(frozen=True)
class StatsPlan:
    """A ``stats:`` entry: the function ``qualname`` in the Python file ``file``."""

    file: Path
    qualname: str


def read_stats_plan(raw: Any, where: str, folder: Path) -> StatsPlan | None:
    """Read ``stats: {plan: path/to/file.py:function}``; ``None`` when absent."""
    if raw is None:
        return None
    spec = raw.get("plan") if isinstance(raw, Mapping) else None
    if not isinstance(spec, str) or set(raw) != {"plan"}:
        raise ProtocolError(
            f"{where}: stats must be {{plan: path/to/file.py:function}}.",
            hint="For example 'stats: {plan: stats/plan.py:plan}', relative to the file.",
        )
    file, colon, qualname = spec.rpartition(":")
    location = (folder / file).resolve()
    if not colon or not file or not qualname or not location.is_file():
        raise ProtocolError(
            f"{where}: stats plan {spec!r} does not name a function in an existing Python file.",
            hint="Write it as 'path/to/file.py:function_name', relative to the file.",
        )
    return StatsPlan(location, qualname)


def replicate_table(study: Any, run: str) -> Any:
    """Return one row per replicate of ``run`` in ``study``, from its stored results.

    Per-frame values are reduced over each replicate's frames the way the
    analysis reduces them: the entry's ``reduce`` for the study's own
    function (``mean``, ``fraction`` as the mean, or ``std`` with ``ddof=1``),
    otherwise the mean. Labelled and multi-part results keep one row per
    label and part.

    Returns
    -------
    pandas.DataFrame
        Columns ``study``, ``condition``, ``replicate``, ``part``, ``label``,
        ``value``, ``unit``, then one per factor the conditions declare.
    """
    stored = study.results(run)
    table = stored.table
    entry = study.protocol.analyses.get(run) if study.protocol is not None else None
    reduce = entry.function.reduce if entry is not None and entry.function is not None else "mean"
    keys = ["condition", "replicate", "part", "label"]
    per_frame = table["frame"].notna()
    framed = table[per_frame]
    if len(framed):
        how = _REDUCE.get(reduce, "mean")
        grouped = framed.groupby(keys, dropna=False, sort=False)
        values = grouped["value"].std(ddof=1) if how == "std" else grouped["value"].mean()
        framed = values.reset_index().assign(unit=grouped["unit"].first().to_numpy())
    rows = table[~per_frame][[*keys, "value", "unit"]]
    import pandas as pd

    result = pd.concat([framed[[*keys, "value", "unit"]], rows], ignore_index=True)
    # The study's condition order, control first, then replicate.
    order = {label: i for i, label in enumerate(study.labels)}
    result = result.sort_values(
        ["condition", "replicate"],
        key=lambda column: column.map(order) if column.name == "condition" else column,
        kind="stable",
    ).reset_index(drop=True)
    label = getattr(study.protocol, "project_label", None) or (
        study.root.name if study.root is not None else None
    )
    result.insert(0, "study", label)
    factors = getattr(study.protocol, "factors", None) or {}
    for name in dict.fromkeys(n for values in factors.values() for n in values):
        result[name] = result["condition"].map(lambda c, name=name: factors.get(c, {}).get(name))
    return result


def trend_tests(report: Any, factors: Mapping[str, Mapping[str, Any]]) -> list[Any]:
    """Fit the slope of the replicate values against each numeric factor of the conditions.

    Uses the per-replicate values of the report's conditions (the values its
    comparisons used), with the replicate as the unit, over the conditions
    that declare the factor. One slope per factor, by ordinary least squares,
    with a two-sided t test of zero slope and a 95 percent interval; the
    factors form one Benjamini-Hochberg family. Labelled results (one value
    per residue, say) get none.
    """
    from scipy import stats

    from polyzymd.analyses.protocols import TrendReport
    from polyzymd.analyses.shared.inferential_statistics import benjamini_hochberg

    if any(item.entry is not None for item in report.conditions):
        return []
    names = [
        name
        for name in dict.fromkeys(n for values in factors.values() for n in values)
        if all(
            isinstance(values.get(name), (int, float)) and not isinstance(values.get(name), bool)
            for values in factors.values()
            if name in values
        )
    ]
    trends = []
    for name in names:
        x, y, used = [], [], []
        for item in report.conditions:
            level = factors.get(item.label, {}).get(name)
            if level is None:
                continue
            used.append(item.label)
            for value in item.replicate_values:
                x.append(float(level))
                y.append(float(value))
        trend = TrendReport(factor=name, conditions=used, n_replicates=len(y))
        if len(set(x)) >= 2 and len(y) >= 3 and len(set(y)) > 1:
            fit = stats.linregress(x, y)
            half = stats.t.ppf(0.975, len(y) - 2) * fit.stderr
            trend = trend.model_copy(
                update={
                    "slope": float(fit.slope),
                    "slope_ci95": (float(fit.slope - half), float(fit.slope + half)),
                    "p": float(fit.pvalue),
                    "r_squared": float(fit.rvalue**2),
                    "testable": True,
                }
            )
        trends.append(trend)
    corrected = benjamini_hochberg([t.p for t in trends])
    return [
        t.model_copy(
            update={
                "p_adjusted": c.adjusted_p_value,
                "significant": c.significant,
                "family_size": sum(1 for x in trends if x.p is not None),
            }
        )
        for t, c in zip(trends, corrected, strict=True)
    ]


def trend_sentence(metric: str, unit: str | None, trend: Any) -> str:
    """Write the verdict sentence of one trend test."""
    from polyzymd.analyses.protocols import VERDICT_NOT_TESTABLE, _interval, _num

    if not trend.testable:
        return (
            f"{VERDICT_NOT_TESTABLE}: trend of {metric} with {trend.factor} needs at least two "
            f"levels, three replicates and values that vary (n {trend.n_replicates})"
        )
    per = f" {unit} per unit {trend.factor}" if unit else f" per unit {trend.factor}"
    evidence = (
        f"slope {_num(trend.slope)}{per}, 95% CI {_interval(trend.slope_ci95)}, "
        f"p_adj {_num(trend.p_adjusted)}, n {trend.n_replicates} replicates over "
        f"{len(trend.conditions)} conditions"
    )
    if trend.significant:
        direction = "rises" if trend.slope > 0 else "falls"
        return f"{metric} {direction} with {trend.factor} ({evidence})"
    return f"no linear trend of {metric} with {trend.factor} detected ({evidence})"


def _file_hash(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _report_hashes(target: Any) -> dict[str, str]:
    """Return the SHA-256 of every stored report the target's studies hold, by study/run."""
    from polyzymd.analyses.results import REPORT_FILE

    studies = (
        [(label, target[label]) for label in target.labels]
        if hasattr(target, "runs_in")
        else [(None, target)]
    )
    hashes = {}
    for label, study in studies:
        for run in study.protocol.analyses:
            report = study.protocol.results_dir(run) / REPORT_FILE
            if report.is_file():
                hashes[f"{label}/{run}" if label else run] = _file_hash(report)
    return hashes


def stats_folder(target: Any, plan: StatsPlan) -> Path:
    """Return the folder of a stats plan's output: ``<root>/results/stats/<function>``."""
    return Path(target.root) / "results" / STATS_FOLDER / plan.qualname


def run_stats_plan(target: Any, plan: StatsPlan) -> Path:
    """Run a ``stats:`` plan on a Study or Project and store what it returns.

    The function receives the target and returns a mapping of names to
    tables (``pandas.DataFrame``, written as ``<name>.csv``) or JSON values
    (written together to ``values.json``). ``record.json`` holds the SHA-256
    of the plan's file, the function name, the SHA-256 of every stored report
    the target held when it ran, and the PolyzyMD version, so
    :func:`stats_status` can tell when it is out of date.

    Returns
    -------
    Path
        The folder holding the output.
    """
    import pandas as pd

    import polyzymd
    from polyzymd.analyses.user_functions import load_function

    function = load_function(plan.file, plan.qualname)
    inputs = _report_hashes(target)
    returned = function(target)
    if not isinstance(returned, Mapping):
        raise ProtocolError(
            f"stats plan {plan.qualname} returned {type(returned).__name__}, not a mapping.",
            hint="Return a dict of names to tables (DataFrames) or numbers and strings.",
        )
    folder = stats_folder(target, plan)
    folder.mkdir(parents=True, exist_ok=True)
    for old in folder.glob("*.csv"):
        old.unlink()
    values: dict[str, Any] = {}
    for name, value in returned.items():
        if isinstance(value, pd.DataFrame):
            value.to_csv(folder / f"{name}.csv", index=False)
        else:
            values[str(name)] = value
    (folder / "values.json").write_text(json.dumps(values, indent=1, default=_json) + "\n")
    root = Path(target.root)
    (folder / STATS_RECORD).write_text(
        json.dumps(
            {
                "plan": {
                    "file": str(plan.file.relative_to(root))
                    if plan.file.is_relative_to(root)
                    else str(plan.file),
                    "function": plan.qualname,
                    "sha256": _file_hash(plan.file),
                },
                "inputs": inputs,
                "polyzymd_version": polyzymd.__version__,
            },
            indent=1,
        )
        + "\n"
    )
    return folder


def stats_status(target: Any, plan: StatsPlan) -> str:
    """Say whether a stats plan's stored output matches its code and the current reports."""
    record_path = stats_folder(target, plan) / STATS_RECORD
    if not record_path.is_file():
        return "not run"
    record = json.loads(record_path.read_text())
    reasons = []
    if record.get("plan", {}).get("sha256") != _file_hash(plan.file):
        reasons.append("its code changed")
    if record.get("inputs") != _report_hashes(target):
        reasons.append("the analysis reports changed")
    return "up to date" if not reasons else "stale: " + " and ".join(reasons)


def _json(value: Any) -> Any:
    """Turn NumPy and pandas scalars into plain JSON values."""
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "tolist"):
        return value.tolist()
    return str(value)
