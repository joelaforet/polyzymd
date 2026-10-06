"""Statistics on stored results: replicate tables and trend tests.

This module reads stored analysis results only; it loads no trajectory.

- :func:`replicate_table` returns one row per replicate of an analysis run
  (and per label or part), the sampling unit of every test.
- :func:`trend_tests` fits, for each numeric factor the conditions declare,
  the slope of the condition means against the factor, and corrects the
  family of factors with Benjamini-Hochberg.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

_REDUCE = {"mean": "mean", "fraction": "mean", "std": "std"}


def replicate_table(study: Any, run: str) -> Any:
    """Return one row per replicate of ``run`` in ``study``, from its stored results.

    Reads ``study.results(run)``; no trajectory is loaded. Per-frame values
    are reduced over each replicate's frames the way the analysis reduces
    them: with the ``reduce`` of the analysis entry's function when the entry
    has one (``mean``; ``fraction`` as the mean; ``std`` with ``ddof=1``),
    otherwise with the mean. Values stored once per replicate are kept as
    they are. Labelled and multi-part results keep one row per label and
    part. Only the replicates the run's report lists are kept, when it lists
    them. Rows are sorted by the study's condition order, then by replicate.

    Parameters
    ----------
    study : Study
        The :class:`~polyzymd.analyses.study.Study` whose results are read.
    run : str
        Name of the analysis run.

    Returns
    -------
    pandas.DataFrame
        Columns ``study``, ``condition``, ``replicate``, ``name`` (the stored
        quantity), ``part``, ``label``, ``value``, ``unit``, then one per factor the conditions declare
        (``None`` where a condition does not declare it). ``study`` holds the
        study's label in its project, or else the name of its folder.

    Raises
    ------
    ProtocolError
        If ``run`` has no stored results (raised by ``study.results``).
    """
    stored = study.results(run)
    table = stored.table
    # Only the replicates the report used, so the table and the report agree.
    used = {
        (c.label, int(r))
        for c in (stored.report.conditions if stored.report is not None else [])
        for r in c.replicates
    }
    if used:
        pairs = zip(table["condition"], table["replicate"])
        table = table[[(c, int(r)) in used for c, r in pairs]]
    entry = study.protocol.analyses.get(run) if study.protocol is not None else None
    reduce = entry.function.reduce if entry is not None and entry.function is not None else "mean"
    # A run can store several quantities (two hydrogen-bond summaries), so the
    # quantity's name is part of every row's identity.
    keys = ["condition", "replicate", "name", "part", "label"]
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
        One :class:`~polyzymd.analyses.protocols.TrendReport` per factor. A
        factor that is not numeric, or has a replicate value that is not
        finite, fewer than three levels (two make it a pairwise comparison),
        or condition means that are all equal, has ``testable=False``, its
        ``reason``, and no slope. The list is empty when the report holds labelled results
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
    trends = []
    for name in dict.fromkeys(n for values in factors.values() for n in values):
        numeric = all(
            isinstance(values.get(name), (int, float)) and not isinstance(values.get(name), bool)
            for values in factors.values()
            if name in values
        )
        levels, means, used, n_values, bad = [], [], [], 0, 0
        for item in report.conditions:
            level = factors.get(item.label, {}).get(name)
            if level is None:
                continue
            values = [float(v) for v in item.replicate_values]
            used.append(item.label)
            n_values += len(values)
            bad += sum(1 for v in values if not math.isfinite(v))
            if values and numeric:
                levels.append(float(level))
                means.append(sum(values) / len(values))
        trend = TrendReport(factor=name, conditions=used, n_replicates=n_values)
        if not numeric:
            reason = "its levels are not all numbers (YAML reads 1e-3 as text; write 1.0e-3)"
        elif bad:
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
