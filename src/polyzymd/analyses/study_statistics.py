"""Statistics on stored results: replicate tables, trend tests and comparison pairs.

This module reads stored analysis results only; it loads no trajectory.

- :func:`replicate_table` returns one row per replicate of an analysis run
  (and per label or part), the sampling unit of every test.
- :func:`comparison_pairs` gives each condition the control it is compared
  with: the study's first condition, or the control of its own stratum.
- :func:`trend_tests` fits, for each numeric factor the conditions declare,
  the slope of the condition means against the factor, and corrects the
  family of factors with Benjamini-Hochberg.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
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


def stratum_control(
    reference: str, factors: Mapping[str, Mapping[str, Any]], within: Sequence[str]
) -> dict[str, Any]:
    """Return the factor values of the control of each stratum, taken from condition ``reference``.

    The values are those of every factor that is not a ``within`` factor,
    with ``None`` for a factor that ``reference`` does not declare. So the
    control of a stratum is its condition whose other factors equal those of
    ``reference``.
    """
    names = dict.fromkeys(name for values in factors.values() for name in values)
    return {name: factors.get(reference, {}).get(name) for name in names if name not in within}


def comparison_pairs(
    chosen: Sequence[str],
    labels: Sequence[str],
    factors: Mapping[str, Mapping[str, Any]],
    within: Sequence[str] = (),
    control: str | Mapping[str, Any] | None = None,
) -> list[tuple[str, str, dict[str, Any] | None]]:
    """Return ``(control, condition, stratum)`` for each condition of ``chosen`` to compare.

    Without ``within`` every condition of ``chosen`` is compared with one
    control: ``control``, or the first of ``labels``; ``stratum`` is ``None``.

    With ``within``, a list of factor names, the conditions of ``labels``
    that share their values of those factors form one stratum, and each
    condition is compared with the control of its stratum; ``stratum`` maps
    each ``within`` factor to its value there. The control of a stratum is
    its one condition whose factors have the values of ``control``, a
    mapping of factor names to values. A condition label as ``control``
    (by default the first of ``labels``) stands for its other factors
    (:func:`stratum_control`). A control is never compared with itself.

    Raises
    ------
    ProtocolError
        If ``control`` maps factor values without ``within``, names a
        ``within`` factor, a condition does not declare a ``within`` factor,
        or a stratum of ``chosen`` has no control or more than one.
    """
    from polyzymd.analyses.exceptions import ProtocolError

    within = list(within)
    if not within:
        if isinstance(control, Mapping):
            raise ProtocolError(
                f"The control {dict(control)} is given as factor values without within.",
                hint="Name the factors that form each stratum with within, or give the control "
                "condition's label.",
            )
        first = control or labels[0]
        return [(first, label, None) for label in chosen if label != first]
    for label in labels:
        missing = [name for name in within if name not in factors.get(label, {})]
        if missing:
            raise ProtocolError(
                f"Condition {label} has no factor {', '.join(missing)}, which within needs.",
                hint="Give every condition each within factor under factors:.",
            )
    if isinstance(control, Mapping):
        if set(control) & set(within):
            raise ProtocolError(
                f"The control {dict(control)} names a within factor.",
                hint="The control's values are of the factors that vary inside a stratum; "
                "give control= with those factors only.",
            )
        wanted = dict(control)
    else:
        wanted = stratum_control(control or labels[0], factors, within)

    def stratum(label: str) -> tuple:
        return tuple(factors[label][name] for name in within)

    controls: dict[tuple, str] = {}
    problems = []
    for level in dict.fromkeys(stratum(label) for label in chosen):
        found = [
            label
            for label in labels
            if stratum(label) == level
            and all(factors[label].get(key) == value for key, value in wanted.items())
        ]
        if len(found) == 1:
            controls[level] = found[0]
            continue
        where = ", ".join(f"{name} {value}" for name, value in zip(within, level, strict=True))
        count = "no control" if not found else f"{len(found)} controls ({', '.join(found)})"
        problems.append(f"the stratum {where} has {count}")
    if problems:
        raise ProtocolError(
            f"Each stratum needs one control, the condition whose factors are {wanted}, but "
            f"{'; '.join(problems)}.",
            hint="Give each stratum one control condition, or set the control's factor values "
            "with control.",
        )
    return [
        (controls[stratum(label)], label, dict(zip(within, stratum(label), strict=True)))
        for label in chosen
        if label != controls[stratum(label)]
    ]


def _number_or_number_text(level: Any) -> bool:
    """Return whether ``level`` is a number, or text such as ``"1e-3"`` that reads as one."""
    if isinstance(level, str):
        try:
            float(level)
        except ValueError:
            return False
        return True
    return isinstance(level, (int, float)) and not isinstance(level, bool)


def trend_tests(report: Any, factors: Mapping[str, Mapping[str, Any]]) -> list[Any]:
    """Fit the slope of the condition means against each numeric factor of the conditions.

    A factor is numeric when every condition that declares it gives an
    ``int`` or ``float`` (not a ``bool``). A factor with any other level,
    such as ``polymer: PEG`` or ``true``, gets no trend. For each numeric factor, each
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
        One :class:`~polyzymd.analyses.protocols.TrendReport` per factor whose
        levels are all numbers or text that reads as a number. A factor with
        a level in text (YAML reads ``1e-3`` as text), a replicate value that
        is not finite, fewer than three levels (two make it a pairwise
        comparison), or condition means that are all equal, has ``testable=False``, its
        ``reason``, and no slope. A condition without partner
        (``no_partner``) is left out of the fit, and ``reason`` says so. The list is empty when the report holds labelled results
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

    from polyzymd.analyses.protocols import FLAT_TREND, TrendReport
    from polyzymd.analyses.shared.inferential_statistics import benjamini_hochberg

    if any(item.entry is not None for item in report.conditions):
        return []
    trends = []
    for name in dict.fromkeys(n for values in factors.values() for n in values):
        declared = [values[name] for values in factors.values() if name in values]
        if not all(_number_or_number_text(level) for level in declared):
            continue
        numeric = not any(isinstance(level, str) for level in declared)
        levels, means, used, n_values, bad, left = [], [], [], 0, 0, []
        for item in report.conditions:
            level = factors.get(item.label, {}).get(name)
            if level is None:
                continue
            if getattr(item, "no_partner", None):
                # Its 0 is not a measurement, so it is no point of the fit.
                left.append(f"{item.label} left out: no partner ({item.no_partner})")
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
            reason = "some levels are text (YAML reads 1e-3 as text; write 1.0e-3)"
        elif bad:
            reason = f"{bad} replicate value{'s are' if bad != 1 else ' is'} not finite"
        elif len(set(levels)) < 3:
            n = len(set(levels))
            reason = f"{n} level{'s' if n != 1 else ''}; a trend needs at least three"
        elif len(set(means)) == 1:
            reason = FLAT_TREND
        else:
            reason = None
        note = "; ".join(left) or None
        if reason is not None:
            trend = trend.model_copy(update={"reason": "; ".join(filter(None, [reason, note]))})
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
                    "reason": note,
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
    from polyzymd.analyses.protocols import FLAT_TREND, VERDICT_NOT_TESTABLE, _interval, _num

    if not trend.testable:
        flat = trend.reason.startswith(FLAT_TREND)
        head = "no trend" if flat else f"{VERDICT_NOT_TESTABLE}: trend"
        return (
            f"{head} of {metric} with {trend.factor}: {trend.reason} "
            f"({len(trend.conditions)} conditions, {trend.n_replicates} replicates)"
        )
    per = f" {unit} per unit {trend.factor}" if unit else f" per unit {trend.factor}"
    evidence = (
        f"slope {_num(trend.slope)}{per}, 95% CI {_interval(trend.slope_ci95)}, "
        f"p_adj {_num(trend.p_adjusted)}, fitted on {len(trend.conditions)} condition means "
        f"of {trend.n_replicates} replicates"
        + (f"; {trend.reason}" if trend.reason else "")
    )
    if trend.significant:
        direction = "rises" if trend.slope > 0 else "falls"
        return f"{metric} {direction} with {trend.factor} ({evidence})"
    return f"no linear trend of {metric} with {trend.factor} detected ({evidence})"
