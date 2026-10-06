"""Statistics on stored results: replicate tables, trend tests and user statistics plans.

This module reads stored analysis results only; it loads no trajectory. Its
main functions are:

- :func:`replicate_table`, which returns one row per replicate of an
  analysis run (and per label or part), the sampling unit of every test.
- :func:`trend_tests`, which fits, for each numeric factor the conditions
  declare, the slope of the replicate values against the factor, with the
  replicate as the unit, and corrects the family of factors with
  Benjamini-Hochberg.
- :func:`read_stats_plan`, which reads the ``stats:`` entry of a
  ``study.yaml`` or ``project.yaml``.
- :func:`run_stats_plan`, which runs the function a ``stats:`` entry names on
  a :class:`~polyzymd.analyses.study.Study` or
  :class:`~polyzymd.analyses.project.Project` and writes what it returns,
  with a record of the function's code and of the reports it read.
- :func:`stats_status`, which compares that record with the current code and
  reports.
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
    """A ``stats:`` entry: the function ``qualname`` in the Python file ``file``.

    :func:`read_stats_plan` builds it from ``stats: {plan: file.py:function}``
    and :func:`run_stats_plan` runs it.

    Attributes
    ----------
    file : Path
        Absolute path of the Python file that defines the function.
    qualname : str
        Name of the function in that file. It also names the output folder
        ``results/stats/<qualname>``.
    """

    file: Path
    qualname: str


def read_stats_plan(raw: Any, where: str, folder: Path) -> StatsPlan | None:
    """Read a ``stats: {plan: path/to/file.py:function}`` entry.

    The value of ``plan`` is split at its last colon into a file path and a
    function name. The path is resolved against ``folder`` and must name an
    existing file; whether the file defines the function is checked only when
    the plan runs.

    Parameters
    ----------
    raw : Any
        The value of the ``stats`` key as parsed from YAML, or ``None`` when
        the key is absent.
    where : str
        Location used in error messages, such as ``"<file>: stats"``.
    folder : Path
        Folder that relative plan paths are resolved against: the folder
        of the ``study.yaml`` or ``project.yaml``.

    Returns
    -------
    StatsPlan or None
        The plan, or ``None`` when ``raw`` is ``None``.

    Raises
    ------
    ProtocolError
        If ``raw`` is not a mapping with the single key ``plan`` holding a
        string, if the string is not of the form ``file:function``, or if
        the file does not exist.
    """
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

    Reads ``study.results(run)``; no trajectory is loaded. Per-frame values
    are reduced over each replicate's frames the way the analysis reduces
    them: with the ``reduce`` of the analysis entry's function when the entry
    has one (``mean``; ``fraction`` as the mean; ``std`` with ``ddof=1``),
    otherwise with the mean. Values stored once per replicate are kept as
    they are. Labelled and multi-part results keep one row per label and
    part. Rows are sorted by the study's condition order, then by replicate.

    Parameters
    ----------
    study : Study
        The :class:`~polyzymd.analyses.study.Study` whose results are read.
    run : str
        Name of the analysis run.

    Returns
    -------
    pandas.DataFrame
        Columns ``study``, ``condition``, ``replicate``, ``part``, ``label``,
        ``value``, ``unit``, then one per factor the conditions declare
        (``None`` where a condition does not declare it). ``study`` holds the
        study's label in its project, or else the name of its folder.

    Raises
    ------
    ProtocolError
        If ``run`` has no stored results (raised by ``study.results``).
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
        # A frame that is not finite makes its replicate's value NaN, as in the
        # report, rather than being skipped as pandas does by default.
        values = grouped["value"].agg(
            (lambda v: v.std(ddof=1, skipna=False))
            if how == "std"
            else (lambda v: v.mean(skipna=False))
        )
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
    """Fit the slope of the condition means against each numeric factor of the conditions.

    A factor is numeric when every condition that declares it gives an
    ``int`` or ``float`` (not a ``bool``). For each numeric factor, each
    condition that declares it contributes one point: the mean of its
    replicate values (the values its comparisons used), at the factor's
    level. The condition, not the replicate, is the unit, because the
    factor varies only between conditions; replicates of one condition
    measure that condition's run-to-run scatter and say nothing more about
    the line between conditions (Hurlbert, 1984; Lazic, 2010). The fit is
    ordinary least squares (:func:`scipy.stats.linregress`) on the means,
    with a two-sided t test of zero slope and a 95 percent t interval, both
    on ``k - 2`` degrees of freedom for ``k`` conditions. The p-values of all
    factors form one Benjamini-Hochberg family.

    Parameters
    ----------
    report : ProtocolReport
        A study's report of one analysis; its ``conditions`` give each
        condition's label, ``entry`` and ``replicate_values``.
    factors : Mapping of str to Mapping of str to Any
        Factor levels by condition label, then by factor name.

    Returns
    -------
    list of TrendReport
        One :class:`~polyzymd.analyses.protocols.TrendReport` per numeric
        factor. A factor with a replicate value that is not finite, fewer
        than three levels (two make it a pairwise comparison), or condition
        means that are all equal has ``testable=False``, its ``reason``, and
        no slope. The list is empty when the report holds labelled results
        (any condition with an ``entry``, such as one value per residue).

    Notes
    -----
    With three to five conditions the test has one to three degrees of
    freedom, so only a clear, steady change across conditions reaches
    significance. The Benjamini-Hochberg step-up procedure (Benjamini and
    Hochberg, 1995) is applied with
    :func:`~polyzymd.analyses.shared.inferential_statistics.benjamini_hochberg`
    at its default ``alpha`` of 0.05.

    References
    ----------
    Hurlbert, S. H. (1984). Pseudoreplication and the design of ecological
    field experiments. Ecological Monographs 54, 187-211.
    Lazic, S. E. (2010). The problem of pseudoreplication in neuroscientific
    studies: is it affecting your analysis? BMC Neuroscience 11, 5.
    doi:10.1186/1471-2202-11-5
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
        levels, means, used, n_values, bad = [], [], [], 0, 0
        for item in report.conditions:
            level = factors.get(item.label, {}).get(name)
            if level is None:
                continue
            values = [float(v) for v in item.replicate_values]
            used.append(item.label)
            n_values += len(values)
            bad += sum(1 for v in values if not math.isfinite(v))
            if values:
                levels.append(float(level))
                means.append(sum(values) / len(values))
        trend = TrendReport(factor=name, conditions=used, n_replicates=n_values)
        if bad:
            reason = f"{bad} replicate value{'s are' if bad != 1 else ' is'} not finite"
        elif len(set(levels)) < 3:
            n = len(set(levels))
            reason = f"{n} level{'s' if n != 1 else ''}; a trend needs at least three"
        elif len(set(means)) == 1:
            reason = "the condition means do not vary"
        else:
            reason = None
        if reason is not None:
            trend = trend.model_copy(update={"reason": reason})
        else:
            fit = stats.linregress(levels, means)
            dof = len(means) - 2
            half = stats.t.ppf(0.975, dof) * fit.stderr
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
    """Write the verdict sentence of one trend test.

    A testable trend gives the slope, its 95 percent interval, the adjusted
    p-value and the number of replicates and conditions, and says whether
    ``metric`` rises or falls with the factor (when significant) or that no
    linear trend was detected. An untestable trend gives the reason.

    Parameters
    ----------
    metric : str
        Name of the measured quantity, used in the sentence.
    unit : str or None
        Unit of the metric; the slope is given in ``<unit> per unit <factor>``.
    trend : TrendReport
        One result of :func:`trend_tests`.

    Returns
    -------
    str
        The sentence, without a final period.
    """
    from polyzymd.analyses.protocols import VERDICT_NOT_TESTABLE, _interval, _num

    if not trend.testable:
        return (
            f"{VERDICT_NOT_TESTABLE}: trend of {metric} with {trend.factor}: {trend.reason} "
            f"({len(trend.conditions)} conditions, {trend.n_replicates} replicates)"
        )
    per = f" {unit} per unit {trend.factor}" if unit else f" per unit {trend.factor}"
    evidence = (
        f"slope {_num(trend.slope)}{per}, 95% CI {_interval(trend.slope_ci95)}, "
        f"p_adj {_num(trend.p_adjusted)}, fitted on {len(trend.conditions)} condition means "
        f"of {trend.n_replicates} replicates"
    )
    if trend.significant:
        direction = "rises" if trend.slope > 0 else "falls"
        return f"{metric} {direction} with {trend.factor} ({evidence})"
    return f"no linear trend of {metric} with {trend.factor} detected ({evidence})"


def _file_hash(path: Path) -> str:
    """Return the SHA-256 hex digest of the bytes of ``path``."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _report_hashes(target: Any) -> dict[str, str]:
    """Return the SHA-256 of every stored report of the target, keyed ``<study>/<run>``.

    For a Project the key is ``<study label>/<run>``; for a Study it is the
    run name alone. Runs without a stored ``report.json`` are left out.
    """
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


def _protocol_hashes(target: Any) -> dict[str, str]:
    """Return a hash of each file that defines the target's factors and analyses.

    Every study's ``study.yaml`` (keyed ``<study label>/study.yaml`` for a
    Project) and, for a Project, ``project.yaml``, each hashed without its
    ``metadata:`` block: editing a factor, a region or an analysis changes
    what :func:`replicate_table` gives a plan, while filling in authors or a
    DOI does not.
    """
    if hasattr(target, "runs_in"):
        hashes = {"project.yaml": _science_hash(target.protocol.path)}
        for label in target.labels:
            hashes[f"{label}/study.yaml"] = _science_hash(target[label].protocol.path)
        return hashes
    return {"study.yaml": _science_hash(target.protocol.path)}


def _science_hash(path: Path) -> str:
    """Return the SHA-256 of a study.yaml or project.yaml without its ``metadata:`` block."""
    import yaml

    data = yaml.safe_load(Path(path).read_text()) or {}
    data.pop("metadata", None)
    return hashlib.sha256(json.dumps(data, sort_keys=True, default=str).encode()).hexdigest()


def _plan_hash(plan: StatsPlan) -> str:
    """Return the SHA-256 of the plan's folder: its file and what it may use there.

    The plan's folder is on ``sys.path`` while it runs, so a helper module or
    data file it reads from there is part of the plan
    (:func:`~polyzymd.analyses.timeseries.code_files`).
    """
    from polyzymd.analyses.timeseries import folder_hash

    return folder_hash(plan.file)


def stats_folder(target: Any, plan: StatsPlan) -> Path:
    """Return the folder of a stats plan's output: ``<root>/results/stats/<function>``.

    The folder is not created.

    Parameters
    ----------
    target : Study or Project
        The study or project the plan runs on; its ``root`` is used.
    plan : StatsPlan
        The plan; its ``qualname`` names the folder.

    Returns
    -------
    Path
        The output folder.
    """
    return Path(target.root) / "results" / STATS_FOLDER / plan.qualname


def run_stats_plan(target: Any, plan: StatsPlan) -> Path:
    """Run a ``stats:`` plan on a Study or Project and write what it returns.

    The plan's function is loaded with
    :func:`~polyzymd.analyses.user_functions.load_function`, called with the
    target, and must return a mapping of names to values. Each
    ``pandas.DataFrame`` value is written to ``<name>.csv`` without its
    index; all other values are written together to ``values.json``
    (NumPy and pandas scalars become plain JSON values). CSV files left from
    an earlier run are deleted first. ``record.json`` holds the plan's file
    (relative to the target's root when inside it), the function name, the
    SHA-256 of the file, the SHA-256 of every stored report the target held
    before the function ran, and the PolyzyMD version; :func:`stats_status`
    compares it with the current state.

    Parameters
    ----------
    target : Study or Project
        The :class:`~polyzymd.analyses.study.Study` or
        :class:`~polyzymd.analyses.project.Project` passed to the function.
    plan : StatsPlan
        The plan to run.

    Returns
    -------
    Path
        The output folder, ``<root>/results/stats/<function>``.

    Raises
    ------
    ProtocolError
        If the function cannot be loaded or does not return a mapping.
    """
    import pandas as pd

    import polyzymd
    from polyzymd.analyses.study_file import portable
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
                    "file": portable(str(plan.file), root),
                    "function": plan.qualname,
                    "sha256": _plan_hash(plan),
                    "hash_of": "module_folder",
                },
                "inputs": inputs,
                "protocol": _protocol_hashes(target),
                "polyzymd_version": polyzymd.__version__,
            },
            indent=1,
        )
        + "\n"
    )
    return folder


def stats_status(target: Any, plan: StatsPlan) -> str:
    """Return whether a stats plan's stored output matches its code and the current reports.

    Compares the SHA-256 of the plan's code (its file and what it may use in
    its folder), of the target's stored reports, and of its ``study.yaml`` and
    ``project.yaml`` files (factors, regions, analyses) with those in the
    plan's ``record.json``.

    Parameters
    ----------
    target : Study or Project
        The study or project the plan runs on.
    plan : StatsPlan
        The plan to check.

    Returns
    -------
    str
        ``"not run"`` when there is no record, ``"up to date"`` when both
        match, otherwise ``"stale: "`` followed by what changed (its code,
        the analysis reports, or both).
    """
    record_path = stats_folder(target, plan) / STATS_RECORD
    if not record_path.is_file():
        return "not run"
    record = json.loads(record_path.read_text())
    reasons = []
    if record.get("plan", {}).get("sha256") != _plan_hash(plan):
        reasons.append("its code changed")
    if record.get("inputs") != _report_hashes(target):
        reasons.append("the analysis reports changed")
    if record.get("protocol") != _protocol_hashes(target):
        reasons.append("the study or project files changed")
    return "up to date" if not reasons else "stale: " + " and ".join(reasons)


def _json(value: Any) -> Any:
    """Convert a value :func:`json.dumps` cannot serialise, as its ``default`` hook.

    A NumPy or pandas scalar is unwrapped to its Python number, with a
    non-finite float as ``None``; an array becomes a list; anything else
    becomes its string.
    """
    if hasattr(value, "tolist") and getattr(value, "ndim", 0) > 0:
        return value.tolist()
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)
