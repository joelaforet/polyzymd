"""Shared statistical utilities for analysis comparisons.

This module provides:

- :func:`pairwise_comparisons` — t-test + Cohen's d + percent-change for
  pairs of conditions (control-vs-all or all-pairs).
- :func:`anova_test` — one-way ANOVA across multiple conditions.
- :func:`rank_conditions` — rank conditions by a metric value.
- :func:`default_scalar_comparison` — the full "simple comparison" pipeline
  used by the default :meth:`Analysis.compare`.
- :func:`interpret_direction` — pick a direction label from percent-change.

These are **utility functions**, not methods on a base class.  Complex
analyses override ``compare()`` entirely and call whichever functions they
need.

See Also
--------
polyzymd.analyses.shared.inferential_statistics : Lower-level statistical
    primitives (t-test, Cohen's d, ANOVA) that this module wraps.
"""

from __future__ import annotations

import json
import logging
import math
from datetime import datetime
from typing import Any

from polyzymd.analyses.base import (
    ANOVAResult,
    ComparisonResult,
    ConditionSummary,
    MetricValue,
    PairwiseResult,
)
from polyzymd.analyses.shared.multi_run_formatting import (
    SINGLE_REPLICATE_SEM_NOTE,
    format_sem_value,
    is_sem_estimable,
)

logger = logging.getLogger("polyzymd.analyses")

NOT_TESTABLE_SINGLETON_NOTE = (
    "Inferential statistics require at least two replicates per condition."
)


# ---------------------------------------------------------------------------
# Direction interpretation
# ---------------------------------------------------------------------------


def interpret_direction(
    pct_change: float,
    direction_labels: tuple[str, str, str] = ("decreased", "unchanged", "increased"),
    threshold: float = 1.0,
) -> str:
    """Interpret percent-change as a direction label.

    Parameters
    ----------
    pct_change : float
        Percent change (negative = decrease, positive = increase).
    direction_labels : tuple[str, str, str]
        ``(negative_label, unchanged_label, positive_label)``.
    threshold : float
        Absolute percent-change below which the result is "unchanged".

    Returns
    -------
    str
        One of the three direction labels.
    """
    if math.isnan(pct_change):
        return direction_labels[1]
    if math.isinf(pct_change):
        return direction_labels[2] if pct_change > 0 else direction_labels[0]
    if abs(pct_change) < threshold:
        return direction_labels[1]
    return direction_labels[0] if pct_change < 0 else direction_labels[2]


def format_pct(pct: float) -> str:
    """Format percent change values for human-readable output.

    Parameters
    ----------
    pct : float
        Percent change value.

    Returns
    -------
    str
        Formatted percent value.

        - ``+inf`` as ``"new (baseline=0)"``
        - ``-inf`` as ``"gone (current=0)"``
        - ``nan`` as ``"undefined"``
        - finite values as signed one-decimal percentages
    """
    pct = pct + 0.0

    if math.isnan(pct):
        return "undefined"
    if math.isinf(pct):
        return "new (baseline=0)" if pct > 0 else "gone (current=0)"
    # Canonical percent format for finite values
    return f"{pct:+.1f}%"


def _format_pct(pct: float) -> str:
    """Backward-compatible alias for ``format_pct``.

    Parameters
    ----------
    pct : float
        Percent change value.

    Returns
    -------
    str
        Formatted percent value.
    """
    return format_pct(pct)


# ---------------------------------------------------------------------------
# Pairwise comparisons
# ---------------------------------------------------------------------------


def pairwise_comparisons(
    metrics_by_condition: dict[str, MetricValue],
    control_label: str | None = None,
    ttest_method: str = "student",
    posthoc_method: str = "ttest_bh",
    fdr_alpha: float = 0.05,
) -> list[PairwiseResult]:
    """Compute pairwise statistical comparisons for a single metric.

    Parameters
    ----------
    metrics_by_condition : dict[str, MetricValue]
        Mapping ``condition_label -> MetricValue``.
    control_label : str | None
        If provided, compare all conditions against this control.
        Otherwise, compare all unique pairs.
    ttest_method : str, optional
        T-test method to use, ``"student"`` or ``"welch"``, by default
        ``"student"``. Only used when ``posthoc_method`` is ``"ttest_bh"``.
    posthoc_method : str, optional
        Pairwise post-hoc method: ``"ttest_bh"`` or ``"tukey_hsd"``, by
        default ``"ttest_bh"``.
    fdr_alpha : float, optional
        Significance threshold for Tukey HSD pairwise comparisons,
        by default 0.05.  Must be in ``(0, 1]``.

    Returns
    -------
    list[PairwiseResult]
        Pairwise comparison results.
    """
    from polyzymd.analyses.shared.inferential_statistics import (
        cohens_d,
        independent_ttest,
        percent_change,
        tukey_hsd,
    )

    # Validate fdr_alpha for all post-hoc methods (BH validates internally,
    # but Tukey uses it directly for significance — guard here uniformly).
    if not (0.0 < fdr_alpha <= 1.0) or math.isnan(fdr_alpha):
        raise ValueError(f"fdr_alpha must be in (0, 1], got {fdr_alpha!r}")

    labels = list(metrics_by_condition.keys())
    results: list[PairwiseResult] = []

    if posthoc_method == "ttest_bh":
        pairs: list[tuple[str, str]]
        if control_label and control_label in metrics_by_condition:
            pairs = [(control_label, lb) for lb in labels if lb != control_label]
        else:
            pairs = [
                (labels[i], labels[j])
                for i in range(len(labels))
                for j in range(i + 1, len(labels))
            ]

        for label_a, label_b in pairs:
            mv_a = metrics_by_condition[label_a]
            mv_b = metrics_by_condition[label_b]
            testable = len(mv_a.replicate_values) >= 2 and len(mv_b.replicate_values) >= 2

            ttest = independent_ttest(
                mv_a.replicate_values,
                mv_b.replicate_values,
                method=ttest_method,
            )
            effect = cohens_d(mv_a.replicate_values, mv_b.replicate_values)
            pct = percent_change(mv_a.mean, mv_b.mean)
            direction = interpret_direction(pct, mv_a.direction_labels)

            results.append(
                PairwiseResult(
                    condition_a=label_a,
                    condition_b=label_b,
                    metric=mv_a.name,
                    t_statistic=ttest.t_statistic,
                    p_value=ttest.p_value,
                    p_value_adjusted=None,
                    posthoc_method=posthoc_method,
                    cohens_d=effect.cohens_d,
                    effect_size_interpretation=effect.interpretation,
                    direction=direction,
                    significant=ttest.significant if testable else False,
                    percent_change=pct,
                    testable=testable,
                    note=None if testable else NOT_TESTABLE_SINGLETON_NOTE,
                )
            )
    elif posthoc_method == "tukey_hsd":
        if len(labels) < 2:
            return results  # Not enough conditions for pairwise Tukey
        group_arrays = [metrics_by_condition[label].replicate_values for label in labels]
        if any(len(group) < 2 for group in group_arrays):
            tukey_results = []
            for i in range(len(labels)):
                for j in range(i + 1, len(labels)):
                    tukey_results.append((i, j, float("nan"), False))
        else:
            tukey_results = [
                (item.group_i, item.group_j, item.p_value, True)
                for item in tukey_hsd(*group_arrays)
            ]

        for group_i, group_j, p_value, testable in tukey_results:
            label_a = labels[group_i]
            label_b = labels[group_j]

            if (
                control_label
                and control_label in metrics_by_condition
                and control_label not in (label_a, label_b)
            ):
                continue

            if control_label and control_label in metrics_by_condition and label_b == control_label:
                label_a, label_b = label_b, label_a

            mv_a = metrics_by_condition[label_a]
            mv_b = metrics_by_condition[label_b]
            effect = cohens_d(mv_a.replicate_values, mv_b.replicate_values)
            pct = percent_change(mv_a.mean, mv_b.mean)
            direction = interpret_direction(pct, mv_a.direction_labels)

            results.append(
                PairwiseResult(
                    condition_a=label_a,
                    condition_b=label_b,
                    metric=mv_a.name,
                    t_statistic=float("nan"),
                    p_value=p_value,
                    p_value_adjusted=p_value if testable else None,
                    posthoc_method=posthoc_method,
                    cohens_d=effect.cohens_d,
                    effect_size_interpretation=effect.interpretation,
                    direction=direction,
                    significant=p_value <= fdr_alpha if testable else False,
                    percent_change=pct,
                    testable=testable,
                    note=None if testable else NOT_TESTABLE_SINGLETON_NOTE,
                )
            )
    else:
        raise ValueError(
            f"Unknown posthoc method {posthoc_method!r}; expected 'ttest_bh' or 'tukey_hsd'"
        )

    return results


# ---------------------------------------------------------------------------
# ANOVA
# ---------------------------------------------------------------------------


def anova_test(
    metrics_by_condition: dict[str, MetricValue],
    metric_name: str = "default",
    alpha: float = 0.05,
) -> ANOVAResult | None:
    """Run one-way ANOVA across conditions for a single metric.

    Parameters
    ----------
    metrics_by_condition : dict[str, MetricValue]
        Mapping ``condition_label -> MetricValue``.
    metric_name : str
        Label for the metric in the result.
    alpha : float, optional
        Inclusive significance threshold for the ANOVA test (``p <= alpha``),
        by default 0.05.

    Returns
    -------
    ANOVAResult | None
        ANOVA result, or ``None`` if fewer than 3 conditions.
    """
    if len(metrics_by_condition) < 3:
        return None

    from polyzymd.analyses.shared.inferential_statistics import one_way_anova

    groups = [mv.replicate_values for mv in metrics_by_condition.values()]
    testable = all(len(group) >= 2 for group in groups)
    result = one_way_anova(*groups)

    return ANOVAResult(
        metric=metric_name,
        f_statistic=result.f_statistic,
        p_value=result.p_value,
        significant=result.p_value <= alpha if testable else False,
        testable=testable,
        note=None if testable else NOT_TESTABLE_SINGLETON_NOTE,
    )


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------


def rank_conditions(
    metrics_by_condition: dict[str, MetricValue],
) -> list[str]:
    """Rank condition labels by metric value.

    Respects ``MetricValue.higher_is_better``: if ``True``, highest mean comes
    first; if ``False``, lowest mean comes first. If ``None``, conditions are
    sorted by metric value descending for neutral display.

    Parameters
    ----------
    metrics_by_condition : dict[str, MetricValue]
        Mapping ``condition_label -> MetricValue``.

    Returns
    -------
    list[str]
        Condition labels in display order.
    """
    if not metrics_by_condition:
        return []
    # All MetricValues should agree on higher_is_better; use the first
    first_label, first_mv = next(iter(metrics_by_condition.items()))
    mismatches = [
        (label, mv.higher_is_better)
        for label, mv in metrics_by_condition.items()
        if mv.higher_is_better != first_mv.higher_is_better
    ]
    if mismatches:
        mismatch_text = ", ".join(
            f"{label}={higher_is_better!r}" for label, higher_is_better in mismatches
        )
        raise ValueError(
            "Inconsistent MetricValue.higher_is_better values across conditions: "
            f"{first_label}={first_mv.higher_is_better!r}; mismatches: {mismatch_text}"
        )
    if first_mv.higher_is_better is None:
        return sorted(
            metrics_by_condition.keys(),
            key=lambda lb: metrics_by_condition[lb].mean,
            reverse=True,
        )
    reverse = first_mv.higher_is_better
    return sorted(
        metrics_by_condition.keys(),
        key=lambda lb: metrics_by_condition[lb].mean,
        reverse=reverse,
    )


# ---------------------------------------------------------------------------
# Scalar comparison utility (full statistical pipeline)
# ---------------------------------------------------------------------------


def default_scalar_comparison(
    analysis_name: str,
    project_name: str,
    metrics_by_condition: dict[str, dict[str, MetricValue]],
    control_label: str | None = None,
    equilibration: str = "0ns",
    fdr_alpha: float = 0.05,
    ttest_method: str = "student",
    posthoc_method: str = "ttest_bh",
) -> ComparisonResult:
    """Run the framework scalar comparison statistics pipeline.

    This utility converts condition-level scalar metrics into the legacy
    ``ComparisonResult`` statistical model used internally by the framework and
    by tests of the statistical helpers. Built-in MDAnalysis plugins should
    expose canonical comparison artifacts; the artifact comparison layer calls
    this function as an implementation detail before wrapping the statistics in
    a ``ComparisonArtifact``.

    For each metric:
    1. Pairwise t-tests + Cohen's d + percent-change.
       Raw pairwise p-values are computed for all pair tests.
    2. ANOVA (if 3+ conditions).
    3. Ranking.
    4. Benjamini-Hochberg (BH) correction is applied once across the
       full family of all pairwise tests across all metrics.

    Parameters
    ----------
    analysis_name : str
        Analysis identifier (e.g. ``"rmsf"``).
    project_name : str
        Comparison project name.
    metrics_by_condition : dict[str, dict[str, MetricValue]]
        Outer key: condition label.  Inner key: metric name.
        Example: ``{"No Polymer": {"mean_rmsf": MetricValue(...)}, ...}``.
    control_label : str | None
        Control condition label.
    equilibration : str
        Equilibration time string.
    fdr_alpha : float, optional
        Significance threshold for pairwise tests and ANOVA.  Used as
        the BH false-discovery-rate threshold (``"ttest_bh"``) or the
        Tukey family-wise threshold (``"tukey_hsd"``), by default 0.05.
    ttest_method : str, optional
        T-test method to use, ``"student"`` or ``"welch"``, by default
        ``"student"``. Only used when ``posthoc_method`` is ``"ttest_bh"``.
    posthoc_method : str, optional
        Pairwise post-hoc method: ``"ttest_bh"`` (pairwise t-tests + BH)
        or ``"tukey_hsd"`` (Tukey HSD), by default ``"ttest_bh"``.

    Returns
    -------
    ComparisonResult
        Structured, serializable comparison result with ``.save()``
        and ``.load()`` methods.
    """
    from polyzymd import __version__

    if not metrics_by_condition:
        raise ValueError("metrics_by_condition is empty — need at least 1 condition.")

    baseline_label, baseline_metrics = next(iter(metrics_by_condition.items()))
    baseline_keys = set(baseline_metrics.keys())
    key_mismatches: list[str] = []
    for label, condition_metrics in metrics_by_condition.items():
        metric_keys = set(condition_metrics.keys())
        missing_keys = sorted(baseline_keys - metric_keys)
        extra_keys = sorted(metric_keys - baseline_keys)
        if missing_keys or extra_keys:
            key_mismatches.append(f"{label}: missing={missing_keys}, extra={extra_keys}")
    if key_mismatches:
        mismatch_details = "; ".join(key_mismatches)
        raise ValueError(
            "Inconsistent metric keys across conditions. "
            f"Baseline '{baseline_label}' keys={sorted(baseline_keys)}; "
            f"differences: {mismatch_details}"
        )

    # Discover metric names from ALL conditions (union), preserving first-seen order
    seen: set[str] = set()
    metric_names: list[str] = []
    for cond_metrics in metrics_by_condition.values():
        for mn in cond_metrics:
            if mn not in seen:
                seen.add(mn)
                metric_names.append(mn)

    if not metric_names:
        raise ValueError("No metric names found in metrics_by_condition.")

    # Warn if control label is specified but not present
    if control_label and control_label not in metrics_by_condition:
        logger.warning(
            "Control label %r not found in metrics_by_condition — "
            "falling back to all-pairs comparison.",
            control_label,
        )
        control_label = None

    all_pairwise: list[PairwiseResult] = []
    all_anova: list[ANOVAResult] = []
    all_rankings: dict[str, list[str]] = {}

    for metric_name in metric_names:
        # Build per-condition MetricValue for this metric
        per_cond: dict[str, MetricValue] = {}
        for label, metrics in metrics_by_condition.items():
            if metric_name in metrics:
                per_cond[label] = metrics[metric_name]

        if not per_cond:
            continue

        # Pairwise comparisons require at least 2 conditions
        if len(per_cond) >= 2:
            pw = pairwise_comparisons(
                per_cond,
                control_label,
                ttest_method=ttest_method,
                posthoc_method=posthoc_method,
                fdr_alpha=fdr_alpha,
            )
            all_pairwise.extend(pw)

        # ANOVA requires at least 3 conditions
        if len(per_cond) >= 3:
            anova = anova_test(per_cond, metric_name, alpha=fdr_alpha)
            if anova is not None:
                all_anova.append(anova)

        # Ranking always works when at least 1 condition exists
        all_rankings[metric_name] = rank_conditions(per_cond)

    # Apply BH correction across ALL pairwise results (full family)
    if all_pairwise and posthoc_method == "ttest_bh":
        from polyzymd.analyses.shared.inferential_statistics import benjamini_hochberg

        raw_p_values = [r.p_value if r.testable else None for r in all_pairwise]
        bh_results = benjamini_hochberg(raw_p_values, alpha=fdr_alpha)
        for result, bh in zip(all_pairwise, bh_results, strict=False):
            result.p_value_adjusted = bh.adjusted_p_value
            if not result.testable:
                result.significant = False
                continue
            p_for_significance = (
                bh.adjusted_p_value if bh.adjusted_p_value is not None else result.p_value
            )
            result.significant = p_for_significance <= fdr_alpha

    # Build condition summaries
    condition_summaries: list[ConditionSummary] = []
    for label, metrics in metrics_by_condition.items():
        extra: dict[str, Any] = {}
        for metric_name, mv in metrics.items():
            extra[f"{metric_name}_mean"] = mv.mean
            extra[f"{metric_name}_sem"] = mv.sem
            extra[f"{metric_name}_replicate_values"] = mv.replicate_values
        if metrics:
            rep_counts = [len(mv.replicate_values) for mv in metrics.values()]
            n_reps = min(rep_counts)
            if len(set(rep_counts)) > 1:
                logger.warning(
                    "Condition %r has inconsistent replicate counts across metrics: %s; "
                    "using minimum (%d)",
                    label,
                    dict(zip(metrics.keys(), rep_counts)),
                    n_reps,
                )
        else:
            n_reps = 0
        condition_summaries.append(ConditionSummary(label=label, n_replicates=n_reps, **extra))

    # Use the first metric's ranking as the primary ranking
    primary_ranking = all_rankings.get(metric_names[0], []) if metric_names else []

    return ComparisonResult(
        analysis_type=analysis_name,
        name=project_name,
        control_label=control_label,
        fdr_alpha=fdr_alpha,
        ttest_method=ttest_method,
        posthoc_method=posthoc_method,
        conditions=condition_summaries,
        pairwise_comparisons=all_pairwise,
        anova=all_anova if all_anova else None,
        ranking=primary_ranking,
        rankings_by_metric=all_rankings if len(metric_names) > 1 else None,
        equilibration_time=equilibration,
        created_at=datetime.now().isoformat(),
        polyzymd_version=__version__,
    )


# ---------------------------------------------------------------------------
# Generic scalar comparison formatter
# ---------------------------------------------------------------------------


def format_scalar_comparison(
    result: ComparisonResult,
    *,
    title: str = "Comparison",
    metric_label: str = "Value",
    metric_unit: str = "",
    metric_key: str | None = None,
    output_format: str = "text",
    higher_is_better: bool | None = True,
) -> str:
    """Format a framework ``ComparisonResult`` for CLI display.

    This formatter is retained for framework/statistics callers that still work
    directly with ``ComparisonResult`` instances. Built-in plugin formatters
    should format canonical ``ComparisonArtifact`` payloads with
    :func:`format_scalar_comparison_artifact_payload` instead.

    Parameters
    ----------
    result : ComparisonResult
        The result to format.
    title : str
        Display title (e.g. ``"RMSF Comparison"``).
    metric_label : str
        Human-readable metric name for table headers (e.g. ``"Mean RMSF"``).
    metric_unit : str
        Unit string appended to values (e.g. ``"A"`` for Angstroms).
    metric_key : str | None
        Key prefix in ``ConditionSummary`` extra fields.  If ``None``,
        auto-detected from the first pairwise comparison's metric name.
    output_format : str
        ``"text"``, ``"markdown"``, or ``"json"``.
    higher_is_better : bool | None
        Affects ranking wording. If ``None``, ranking direction is reported as
        neutral.

    Returns
    -------
    str
        Formatted output string.
    """
    if output_format == "json":
        return result.model_dump_json(indent=2)

    # Auto-detect metric key from pairwise comparisons
    if metric_key is None and result.pairwise_comparisons:
        metric_key = result.pairwise_comparisons[0].metric
    if metric_key is None:
        metric_key = "value"

    # Helper to get metric value from a condition summary
    def _get_mean(cond: ConditionSummary) -> float:
        return getattr(cond, f"{metric_key}_mean", 0.0)

    def _get_sem(cond: ConditionSummary) -> float:
        return getattr(cond, f"{metric_key}_sem", 0.0)

    def _get_cond(label: str) -> ConditionSummary:
        for c in result.conditions:
            if c.label == label:
                return c
        raise KeyError(f"Condition '{label}' not found")

    unit_str = f" {metric_unit}" if metric_unit else ""

    if output_format == "markdown":
        return _format_scalar_markdown(
            result,
            title,
            metric_label,
            unit_str,
            metric_key,
            _get_mean,
            _get_sem,
            _get_cond,
            higher_is_better,
        )

    return _format_scalar_text(
        result,
        title,
        metric_label,
        unit_str,
        metric_key,
        _get_mean,
        _get_sem,
        _get_cond,
        higher_is_better,
    )


def format_scalar_comparison_artifact_payload(
    payload: dict[str, Any],
    *,
    title: str = "Comparison",
    metric_label: str = "Value",
    metric_unit: str = "",
    metric_key: str | None = None,
    output_format: str = "text",
    higher_is_better: bool | None = True,
) -> str:
    """Format a canonical comparison artifact payload for CLI display.

    Parameters
    ----------
    payload : dict[str, Any]
        Canonical comparison artifact payload produced by the default scalar
        comparison pipeline.
    title : str, optional
        Display title, by default ``"Comparison"``.
    metric_label : str, optional
        Human-readable metric label, by default ``"Value"``.
    metric_unit : str, optional
        Unit suffix for metric values, by default ``""``.
    metric_key : str | None, optional
        Metric key to format. If omitted, the first ranking or pairwise metric
        in the payload is used.
    output_format : str, optional
        ``"text"``, ``"markdown"``, or ``"json"``, by default ``"text"``.
    higher_is_better : bool | None, optional
        Ranking direction used for interpretation text, by default ``True``.

    Returns
    -------
    str
        Formatted comparison output.
    """

    if output_format == "json":
        return json.dumps(payload, indent=2)

    metric_key = _select_payload_metric_key(payload, metric_key)
    unit_str = f" {metric_unit}" if metric_unit else ""
    conditions = payload.get("condition_summaries", [])
    condition_by_label = {str(item.get("label", "")): item for item in conditions}
    ranking = _payload_ranking(payload, metric_key)
    pairwise = [
        item
        for item in payload.get("pairwise_comparisons", [])
        if not metric_key or item.get("metric") == metric_key
    ]
    statistical_parameters = payload.get("statistical_parameters", {})
    comparison_name = str(statistical_parameters.get("project_name", payload.get("name", "")))
    equilibration = str(statistical_parameters.get("equilibration", "0ns"))
    control_label = payload.get("effective_control") or payload.get("control_label")

    if output_format == "markdown":
        lines = [f"# {title}: {comparison_name}", ""]
        lines.append(f"**Equilibration:** {equilibration}")
        if control_label:
            lines.append(f"**Control:** {control_label}")
        lines.extend(["", f"## Condition ranking ({metric_label})", ""])
        lines.append(f"| Rank | Condition | {metric_label} | SEM |")
        lines.append("|---:|---|---:|---:|")
        for item in ranking:
            label = str(item.get("label", ""))
            condition = condition_by_label.get(label, {})
            lines.append(
                f"| {item.get('rank', '')} | {label} | "
                f"{_payload_metric(condition, metric_key, 'mean'):.4f}{unit_str} | "
                f"{_payload_metric(condition, metric_key, 'sem'):.4f}{unit_str} |"
            )
        return "\n".join(lines)

    lines = ["", f"{title}: {comparison_name}", "=" * 60, f"Equilibration: {equilibration}"]
    if control_label:
        lines.append(f"Control: {control_label}")
    lines.extend(["", f"Condition ranking ({metric_label}):"])
    direction = _payload_direction_label(higher_is_better)
    if direction:
        lines.append(direction)
    for item in ranking:
        label = str(item.get("label", ""))
        condition = condition_by_label.get(label, {})
        lines.append(
            f"  {item.get('rank', '')}. {label}: "
            f"{_payload_metric(condition, metric_key, 'mean'):.4f}{unit_str} ± "
            f"{_payload_metric(condition, metric_key, 'sem'):.4f}{unit_str}"
        )
    if pairwise:
        lines.extend(["", "Pairwise comparisons:"])
        for item in pairwise:
            lines.append(
                f"  {item.get('condition_a')} vs {item.get('condition_b')}: "
                f"Δ={float(item.get('effect_size', item.get('mean_difference', 0.0))):.4f}, "
                f"p={float(item.get('p_value', 1.0)):.4g}"
            )
    return "\n".join(lines)


def _select_payload_metric_key(payload: dict[str, Any], metric_key: str | None) -> str:
    """Return the metric key to format from a comparison payload."""

    if metric_key is not None:
        return metric_key
    rankings_by_metric = payload.get("rankings_by_metric") or {}
    if rankings_by_metric:
        return str(next(iter(rankings_by_metric)))
    pairwise = payload.get("pairwise_comparisons", [])
    if pairwise:
        return str(pairwise[0].get("metric", "value"))
    return "value"


def _payload_ranking(payload: dict[str, Any], metric_key: str) -> list[dict[str, Any]]:
    """Return the selected ranking list from a comparison payload."""

    rankings_by_metric = payload.get("rankings_by_metric") or {}
    ranking = rankings_by_metric.get(metric_key, payload.get("ranking", []))
    normalized: list[dict[str, Any]] = []
    for index, item in enumerate(ranking, start=1):
        if isinstance(item, dict):
            normalized.append(item)
        elif isinstance(item, str):
            normalized.append({"rank": index, "label": item})
    return normalized


def _payload_metric(condition: dict[str, Any], metric_key: str, suffix: str) -> float:
    """Return a scalar metric from a payload condition summary."""

    value = condition.get(f"{metric_key}_{suffix}", condition.get(suffix, 0.0))
    return float(value if value is not None else 0.0)


def _payload_direction_label(higher_is_better: bool | None) -> str:
    """Return a ranking direction label for payload formatting."""

    if higher_is_better is True:
        return "Higher values rank first."
    if higher_is_better is False:
        return "Lower values rank first."
    return ""


def _format_scalar_text(
    result: ComparisonResult,
    title: str,
    metric_label: str,
    unit_str: str,
    metric_key: str,
    _get_mean,
    _get_sem,
    _get_cond,
    higher_is_better: bool | None,
) -> str:
    """Format as ASCII table (internal)."""
    alpha = result.fdr_alpha if result.fdr_alpha is not None else 0.05

    selected_ranking = result.ranking
    if result.rankings_by_metric and metric_key in result.rankings_by_metric:
        selected_ranking = result.rankings_by_metric[metric_key]

    selected_pairwise = result.pairwise_comparisons
    if metric_key and any(comp.metric == metric_key for comp in result.pairwise_comparisons):
        selected_pairwise = [
            comp for comp in result.pairwise_comparisons if comp.metric == metric_key
        ]

    selected_anova = result.anova
    if result.anova and any(a.metric == metric_key for a in result.anova):
        selected_anova = [a for a in result.anova if a.metric == metric_key]

    lines: list[str] = []

    # Header
    lines.append("")
    lines.append(f"{title}: {result.name}")
    lines.append("=" * 60)
    lines.append(f"Equilibration: {result.equilibration_time}")
    if result.control_label:
        lines.append(f"Control: {result.control_label}")
    lines.append("")

    # Condition ranking table
    if higher_is_better is None:
        rank_desc = "no direction preference"
    else:
        rank_desc = "highest first" if higher_is_better else "lowest first"
    lines.append(f"Condition Summary (ranked by {metric_label}, {rank_desc})")
    lines.append("-" * 60)
    header = f"{'Rank':<5} {'Condition':<20} {metric_label:<14} {'SEM':<10} {'N':<4}"
    lines.append(header)
    lines.append("-" * 60)

    for rank, label in enumerate(selected_ranking, 1):
        cond = _get_cond(label)
        marker = "*" if label == result.control_label else " "
        mean_val = _get_mean(cond)
        sem_val = _get_sem(cond)
        sem_str = format_sem_value(sem_val, cond.n_replicates, precision=4)
        lines.append(
            f"{rank:<5} {cond.label:<20} {mean_val:>10.4f}{unit_str}  "
            f"{sem_str:>8}  {cond.n_replicates:<4}{marker}"
        )

    lines.append("-" * 60)
    if result.control_label:
        lines.append("* = control condition")
    if any(not is_sem_estimable(_get_cond(label).n_replicates) for label in selected_ranking):
        lines.append(SINGLE_REPLICATE_SEM_NOTE)
    lines.append("")

    # Pairwise comparisons
    if selected_pairwise:
        lines.append("Pairwise Comparisons")
        lines.append("-" * 80)
        cohens_label = "Cohen's d"
        posthoc_method = selected_pairwise[0].posthoc_method if selected_pairwise else "ttest_bh"
        has_adjusted = any(c.p_value_adjusted is not None for c in selected_pairwise)
        if posthoc_method == "ttest_bh":
            if has_adjusted:
                header = (
                    f"{'Comparison':<30} {'% Change':<10} {'t-stat':<10} {'p-value':<12} "
                    f"{'p (adj)':<12} {cohens_label:<10} {'Effect':<12}"
                )
            else:
                header = (
                    f"{'Comparison':<30} {'% Change':<10} {'t-stat':<10} {'p-value':<12} "
                    f"{cohens_label:<10} {'Effect':<12}"
                )
        elif posthoc_method == "tukey_hsd":
            header = (
                f"{'Comparison':<30} {'% Change':<10} {'Tukey p':<12} "
                f"{cohens_label:<10} {'Effect':<12}"
            )
        else:
            raise ValueError(f"Unknown posthoc method {posthoc_method!r}")
        lines.append(header)
        lines.append("-" * 80)

        for comp in selected_pairwise:
            name = f"{comp.condition_b} vs {comp.condition_a}"
            sig_marker = "*" if comp.significant else ""
            pct_str = format_pct(comp.percent_change)
            d_str = f"{comp.cohens_d:.2f}" if comp.testable else "n/a"
            if comp.posthoc_method == "ttest_bh":
                p_str = f"{comp.p_value:.4f}{sig_marker}" if comp.testable else "not testable"
                p_adj = comp.p_value_adjusted
                p_adj_str = "n/a" if p_adj is None else f"{p_adj:.4f}{sig_marker}"
                t_str = (
                    f"{comp.t_statistic:.3f}"
                    if comp.testable and not math.isnan(comp.t_statistic)
                    else "n/a"
                )
                if has_adjusted:
                    lines.append(
                        f"{name:<30} {pct_str:<10} {t_str:<10} {p_str:<12} {p_adj_str:<12} "
                        f"{d_str:<10} {comp.effect_size_interpretation:<12}"
                    )
                else:
                    lines.append(
                        f"{name:<30} {pct_str:<10} {t_str:<10} {p_str:<12} "
                        f"{d_str:<10} {comp.effect_size_interpretation:<12}"
                    )
            elif comp.posthoc_method == "tukey_hsd":
                p_str = f"{comp.p_value:.4f}{sig_marker}" if comp.testable else "not testable"
                lines.append(
                    f"{name:<30} {pct_str:<10} {p_str:<12} "
                    f"{d_str:<10} {comp.effect_size_interpretation:<12}"
                )
            else:
                raise ValueError(f"Unknown posthoc method {comp.posthoc_method!r}")

        lines.append("-" * 80)
        if posthoc_method == "ttest_bh" and has_adjusted and result.fdr_alpha is not None:
            lines.append(f"* p_adj <= {result.fdr_alpha} (BH-corrected)")
        elif posthoc_method == "ttest_bh" and has_adjusted:
            lines.append("* BH-corrected p_adj significant")
        elif posthoc_method == "ttest_bh":
            lines.append(f"* p <= {alpha:g}")
        elif posthoc_method == "tukey_hsd":
            lines.append(f"* Tukey HSD p <= {alpha:g}")
        else:
            raise ValueError(f"Unknown posthoc method {posthoc_method!r}")
        if any(not comp.testable for comp in selected_pairwise):
            lines.append(
                "Not testable: inferential statistics require at least two replicates per condition"
            )
        lines.append("")

    # ANOVA
    if selected_anova:
        lines.append("One-way ANOVA")
        lines.append("-" * 40)
        for a in selected_anova:
            sig = "Yes" if a.significant else "No"
            if a.metric != "default":
                lines.append(f"Metric: {a.metric}")
            if a.testable:
                lines.append(f"F-statistic: {a.f_statistic:.3f}")
                lines.append(f"p-value:     {a.p_value:.4f}")
                lines.append(f"Significant: {sig} (alpha={alpha:g})")
            else:
                lines.append("F-statistic: n/a")
                lines.append("p-value:     not testable")
                lines.append("Significant: No")
                lines.append(f"Note: {a.note or NOT_TESTABLE_SINGLETON_NOTE}")
        lines.append("")

    # Interpretation
    if selected_ranking:
        lines.append("Interpretation")
        lines.append("-" * 60)
        top_label = selected_ranking[0]
        top_cond = _get_cond(top_label)
        top_val = _get_mean(top_cond)
        if higher_is_better is None:
            lines.append(
                f"Highest value: {top_label} (listed by value; no better/worse direction implied)"
            )
        else:
            lines.append(f"Best: {top_label} ({metric_label} = {top_val:.4f}{unit_str})")

        if result.control_label and top_label != result.control_label:
            from polyzymd.analyses.shared.inferential_statistics import percent_change

            ctrl = _get_cond(result.control_label)
            ctrl_val = _get_mean(ctrl)
            pct = percent_change(ctrl_val, top_val)
            if not math.isnan(pct):
                direction = interpret_direction(pct, ("lower", "unchanged", "higher"))
                if direction == "unchanged":
                    lines.append(f"  -> unchanged relative to control ({result.control_label})")
                else:
                    magnitude = format_pct(pct).lstrip("+-")
                    lines.append(
                        f"  -> {magnitude} {direction} than control ({result.control_label})"
                    )
        lines.append("")

    lines.append(f"Analysis completed: {result.created_at}")
    lines.append(f"PolyzyMD version: {result.polyzymd_version}")
    lines.append("")

    return "\n".join(lines)


def _format_scalar_markdown(
    result: ComparisonResult,
    title: str,
    metric_label: str,
    unit_str: str,
    metric_key: str,
    _get_mean,
    _get_sem,
    _get_cond,
    higher_is_better: bool | None,
) -> str:
    """Format as Markdown (internal)."""
    alpha = result.fdr_alpha if result.fdr_alpha is not None else 0.05

    selected_ranking = result.ranking
    if result.rankings_by_metric and metric_key in result.rankings_by_metric:
        selected_ranking = result.rankings_by_metric[metric_key]

    selected_pairwise = result.pairwise_comparisons
    if metric_key and any(comp.metric == metric_key for comp in result.pairwise_comparisons):
        selected_pairwise = [
            comp for comp in result.pairwise_comparisons if comp.metric == metric_key
        ]

    selected_anova = result.anova
    if result.anova and any(a.metric == metric_key for a in result.anova):
        selected_anova = [a for a in result.anova if a.metric == metric_key]

    lines: list[str] = []

    lines.append(f"# {title}: {result.name}")
    lines.append("")
    lines.append("## Parameters")
    lines.append("")
    lines.append(f"- **Equilibration:** {result.equilibration_time}")
    if result.control_label:
        lines.append(f"- **Control:** {result.control_label}")
    lines.append(f"- **Date:** {result.created_at}")
    lines.append(f"- **Version:** {result.polyzymd_version}")
    lines.append("")

    # Condition table
    lines.append("## Condition Summary")
    lines.append("")
    lines.append(f"| Rank | Condition | {metric_label}{unit_str} | SEM | N |")
    lines.append(
        "|------|-----------|" + "-" * (len(metric_label) + len(unit_str) + 2) + "|-----|---|"
    )

    for rank, label in enumerate(selected_ranking, 1):
        cond = _get_cond(label)
        marker = " (control)" if label == result.control_label else ""
        mean_val = _get_mean(cond)
        sem_val = _get_sem(cond)
        sem_str = format_sem_value(sem_val, cond.n_replicates, precision=4)
        lines.append(
            f"| {rank} | **{cond.label}**{marker} | "
            f"{mean_val:.4f} | {sem_str} | {cond.n_replicates} |"
        )

    if any(not is_sem_estimable(_get_cond(label).n_replicates) for label in selected_ranking):
        lines.append("")
        lines.append(f"*{SINGLE_REPLICATE_SEM_NOTE}.*")

    lines.append("")

    # Pairwise
    if selected_pairwise:
        lines.append("## Statistical Comparisons")
        lines.append("")
        posthoc_method = selected_pairwise[0].posthoc_method if selected_pairwise else "ttest_bh"
        has_adjusted = any(c.p_value_adjusted is not None for c in selected_pairwise)
        if posthoc_method == "ttest_bh":
            if has_adjusted:
                lines.append(
                    "| Comparison | % Change | t-stat | p-value | p (adj) | Cohen's d | Effect | Sig |"
                )
                lines.append(
                    "|------------|----------|--------|---------|---------|-----------|--------|-----|"
                )
            else:
                lines.append(
                    "| Comparison | % Change | t-stat | p-value | Cohen's d | Effect | Sig |"
                )
                lines.append(
                    "|------------|----------|--------|---------|-----------|--------|-----|"
                )
        elif posthoc_method == "tukey_hsd":
            lines.append("| Comparison | % Change | Tukey p | Cohen's d | Effect | Sig |")
            lines.append("|------------|----------|---------|-----------|--------|-----|")
        else:
            raise ValueError(f"Unknown posthoc method {posthoc_method!r}")
        for comp in selected_pairwise:
            name = f"{comp.condition_b} vs {comp.condition_a}"
            sig = "Yes" if comp.significant else "No"
            d_str = f"{comp.cohens_d:.2f}" if comp.testable else "n/a"
            if comp.posthoc_method == "ttest_bh":
                p_adj = "n/a" if comp.p_value_adjusted is None else f"{comp.p_value_adjusted:.4f}"
                p_value = f"{comp.p_value:.4f}" if comp.testable else "not testable"
                t_str = (
                    f"{comp.t_statistic:.3f}"
                    if comp.testable and not math.isnan(comp.t_statistic)
                    else "n/a"
                )
                if has_adjusted:
                    lines.append(
                        f"| {name} | {format_pct(comp.percent_change)} | "
                        f"{t_str} | {p_value} | {p_adj} | {d_str} | "
                        f"{comp.effect_size_interpretation} | {sig} |"
                    )
                else:
                    lines.append(
                        f"| {name} | {format_pct(comp.percent_change)} | "
                        f"{t_str} | {p_value} | {d_str} | "
                        f"{comp.effect_size_interpretation} | {sig} |"
                    )
            elif comp.posthoc_method == "tukey_hsd":
                p_value = f"{comp.p_value:.4f}" if comp.testable else "not testable"
                lines.append(
                    f"| {name} | {format_pct(comp.percent_change)} | "
                    f"{p_value} | {d_str} | "
                    f"{comp.effect_size_interpretation} | {sig} |"
                )
            else:
                raise ValueError(f"Unknown posthoc method {comp.posthoc_method!r}")
        if posthoc_method == "ttest_bh" and has_adjusted and result.fdr_alpha is not None:
            lines.append("")
            lines.append(f"*Significance uses BH-adjusted p-values: p_adj <= {result.fdr_alpha}.*")
        elif posthoc_method == "ttest_bh" and not has_adjusted:
            lines.append("")
            lines.append(f"*Significance uses raw p <= {alpha:g}.*")
        elif posthoc_method == "tukey_hsd":
            lines.append("")
            lines.append(f"*Significance uses Tukey HSD p <= {alpha:g}.*")
        if any(not comp.testable for comp in selected_pairwise):
            lines.append("")
            lines.append(
                "*Not testable: inferential statistics require at least two replicates per condition.*"
            )
        lines.append("")

    # ANOVA
    if selected_anova:
        lines.append("## ANOVA")
        lines.append("")
        for a in selected_anova:
            sig = "Yes" if a.significant else "No"
            if a.testable:
                lines.append(f"- **F-statistic:** {a.f_statistic:.3f}")
                lines.append(f"- **p-value:** {a.p_value:.4f}")
            else:
                lines.append("- **F-statistic:** n/a")
                lines.append("- **p-value:** not testable")
                lines.append(f"- **Note:** {a.note or NOT_TESTABLE_SINGLETON_NOTE}")
            lines.append(f"- **Significant:** {sig}")
        lines.append("")

    # Key findings
    if selected_ranking:
        lines.append("## Key Findings")
        lines.append("")
        top_label = selected_ranking[0]
        top_cond = _get_cond(top_label)
        top_val = _get_mean(top_cond)
        if higher_is_better is None:
            lines.append(
                "1. **Highest value:** "
                f"{top_label} (listed by value; no better/worse direction implied)"
            )
        else:
            lines.append(
                f"1. **Best condition:** {top_label} ({metric_label} = {top_val:.4f}{unit_str})"
            )
        if result.control_label and top_label != result.control_label:
            from polyzymd.analyses.shared.inferential_statistics import percent_change

            ctrl = _get_cond(result.control_label)
            ctrl_val = _get_mean(ctrl)
            pct = percent_change(ctrl_val, top_val)
            if not math.isnan(pct):
                direction = interpret_direction(pct, ("lower", "unchanged", "higher"))
                if direction == "unchanged":
                    lines.append("2. unchanged relative to control")
                else:
                    magnitude = format_pct(pct).lstrip("+-")
                    lines.append(f"2. {magnitude} {direction} than control")
        lines.append("")

    return "\n".join(lines)
