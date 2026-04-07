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
polyzymd.compare.statistics : Lower-level statistical primitives (t-test,
    Cohen's d, ANOVA) that this module wraps.
"""

from __future__ import annotations

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

logger = logging.getLogger("polyzymd.analyses")


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

        - ``+inf`` as ``"+∞%"``
        - ``-inf`` as ``"-∞%"``
        - ``nan`` as ``"n/a"``
        - finite values as signed one-decimal percentages
    """
    pct = pct + 0.0

    if math.isnan(pct):
        return "n/a"
    if math.isinf(pct):
        return "+∞%" if pct > 0 else "-∞%"
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
) -> list[PairwiseResult]:
    """Compute pairwise statistical comparisons for a single metric.

    .. note::

       When 3+ conditions are compared, the pairwise p-values are **not**
       corrected for multiple comparisons (no Bonferroni or
       Benjamini-Hochberg adjustment).  The accompanying ANOVA provides an
       omnibus test, but individual pairwise p-values should be interpreted
       cautiously.  Multiple-comparison correction is a planned enhancement.

    Parameters
    ----------
    metrics_by_condition : dict[str, MetricValue]
        Mapping ``condition_label -> MetricValue``.
    control_label : str | None
        If provided, compare all conditions against this control.
        Otherwise, compare all unique pairs.

    Returns
    -------
    list[PairwiseResult]
        Pairwise comparison results.
    """
    from polyzymd.compare.statistics import (
        cohens_d,
        independent_ttest,
        percent_change,
    )

    labels = list(metrics_by_condition.keys())
    results: list[PairwiseResult] = []

    pairs: list[tuple[str, str]]
    if control_label and control_label in metrics_by_condition:
        pairs = [(control_label, lb) for lb in labels if lb != control_label]
    else:
        pairs = [
            (labels[i], labels[j]) for i in range(len(labels)) for j in range(i + 1, len(labels))
        ]

    for label_a, label_b in pairs:
        mv_a = metrics_by_condition[label_a]
        mv_b = metrics_by_condition[label_b]

        ttest = independent_ttest(mv_a.replicate_values, mv_b.replicate_values)
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
                cohens_d=effect.cohens_d,
                effect_size_interpretation=effect.interpretation,
                direction=direction,
                significant=ttest.significant,
                percent_change=pct,
            )
        )

    return results


# ---------------------------------------------------------------------------
# ANOVA
# ---------------------------------------------------------------------------


def anova_test(
    metrics_by_condition: dict[str, MetricValue],
    metric_name: str = "default",
) -> ANOVAResult | None:
    """Run one-way ANOVA across conditions for a single metric.

    Parameters
    ----------
    metrics_by_condition : dict[str, MetricValue]
        Mapping ``condition_label -> MetricValue``.
    metric_name : str
        Label for the metric in the result.

    Returns
    -------
    ANOVAResult | None
        ANOVA result, or ``None`` if fewer than 3 conditions.
    """
    if len(metrics_by_condition) < 3:
        return None

    from polyzymd.compare.statistics import one_way_anova

    groups = [mv.replicate_values for mv in metrics_by_condition.values()]
    result = one_way_anova(*groups)

    return ANOVAResult(
        metric=metric_name,
        f_statistic=result.f_statistic,
        p_value=result.p_value,
        significant=result.significant,
    )


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------


def rank_conditions(
    metrics_by_condition: dict[str, MetricValue],
) -> list[str]:
    """Rank condition labels by metric value.

    Respects ``MetricValue.higher_is_better``: if ``True``, highest
    mean comes first; if ``False``, lowest mean comes first.

    Parameters
    ----------
    metrics_by_condition : dict[str, MetricValue]
        Mapping ``condition_label -> MetricValue``.

    Returns
    -------
    list[str]
        Condition labels in ranked order (best first).
    """
    if not metrics_by_condition:
        return []
    # All MetricValues should agree on higher_is_better; use the first
    first_mv = next(iter(metrics_by_condition.values()))
    reverse = first_mv.higher_is_better
    return sorted(
        metrics_by_condition.keys(),
        key=lambda lb: metrics_by_condition[lb].mean,
        reverse=reverse,
    )


# ---------------------------------------------------------------------------
# Default scalar comparison (full pipeline)
# ---------------------------------------------------------------------------


def default_scalar_comparison(
    analysis_name: str,
    project_name: str,
    metrics_by_condition: dict[str, dict[str, MetricValue]],
    control_label: str | None = None,
    equilibration: str = "0ns",
) -> ComparisonResult:
    """Run the standard scalar comparison pipeline.

    This is the default implementation used by :meth:`Analysis.compare`
    for analyses with one or more simple scalar metrics (RMSF, triad,
    secondary structure).

    For each metric:
    1. Pairwise t-tests + Cohen's d + percent-change.
    2. ANOVA (if 3+ conditions).
    3. Ranking.

    .. note::

       Pairwise t-tests use Welch's test (unequal variance).  When 3+
       conditions are compared, pairwise p-values are **not** corrected
       for multiple comparisons.  ANOVA provides an omnibus test; interpret
       individual pairwise p-values with caution.  Multiple-comparison
       correction (Bonferroni / Benjamini-Hochberg) is a planned
       enhancement.

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

    Returns
    -------
    ComparisonResult
        Structured, serializable comparison result with ``.save()``
        and ``.load()`` methods.
    """
    from polyzymd import __version__

    if not metrics_by_condition:
        raise ValueError("metrics_by_condition is empty — need at least 2 conditions.")

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

        if len(per_cond) < 2:
            continue

        # Pairwise
        pw = pairwise_comparisons(per_cond, control_label)
        all_pairwise.extend(pw)

        # ANOVA
        anova = anova_test(per_cond, metric_name)
        if anova is not None:
            all_anova.append(anova)

        # Ranking
        all_rankings[metric_name] = rank_conditions(per_cond)

    # Build condition summaries
    condition_summaries: list[ConditionSummary] = []
    for label, metrics in metrics_by_condition.items():
        extra: dict[str, Any] = {}
        for metric_name, mv in metrics.items():
            extra[f"{metric_name}_mean"] = mv.mean
            extra[f"{metric_name}_sem"] = mv.sem
            extra[f"{metric_name}_replicate_values"] = mv.replicate_values
        n_reps = len(next(iter(metrics.values())).replicate_values) if metrics else 0
        condition_summaries.append(ConditionSummary(label=label, n_replicates=n_reps, **extra))

    # Use the first metric's ranking as the primary ranking
    primary_ranking = all_rankings.get(metric_names[0], []) if metric_names else []

    return ComparisonResult(
        analysis_type=analysis_name,
        name=project_name,
        control_label=control_label,
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
    higher_is_better: bool = True,
) -> str:
    """Format a :class:`ComparisonResult` for CLI display.

    This is a generic formatter for analyses that use the default scalar
    comparison pipeline.  It handles text and markdown output with
    condition rankings, pairwise statistical tests, ANOVA, and
    interpretation.

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
    higher_is_better : bool
        Affects interpretation wording (lower RMSF = more stable, etc.).

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


def _format_scalar_text(
    result: ComparisonResult,
    title: str,
    metric_label: str,
    unit_str: str,
    metric_key: str,
    _get_mean,
    _get_sem,
    _get_cond,
    higher_is_better: bool,
) -> str:
    """Format as ASCII table (internal)."""
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
    rank_desc = "lowest first" if not higher_is_better else "highest first"
    lines.append(f"Condition Summary (ranked by {metric_label}, {rank_desc})")
    lines.append("-" * 60)
    header = f"{'Rank':<5} {'Condition':<20} {metric_label:<14} {'SEM':<10} {'N':<4}"
    lines.append(header)
    lines.append("-" * 60)

    for rank, label in enumerate(result.ranking, 1):
        cond = _get_cond(label)
        marker = "*" if label == result.control_label else " "
        mean_val = _get_mean(cond)
        sem_val = _get_sem(cond)
        lines.append(
            f"{rank:<5} {cond.label:<20} {mean_val:>10.4f}{unit_str}  "
            f"{sem_val:>8.4f}  {cond.n_replicates:<4}{marker}"
        )

    lines.append("-" * 60)
    if result.control_label:
        lines.append("* = control condition")
    lines.append("")

    # Pairwise comparisons
    if result.pairwise_comparisons:
        lines.append("Pairwise Comparisons")
        lines.append("-" * 80)
        cohens_label = "Cohen's d"
        header = (
            f"{'Comparison':<30} {'% Change':<10} {'p-value':<12} {cohens_label:<10} {'Effect':<12}"
        )
        lines.append(header)
        lines.append("-" * 80)

        for comp in result.pairwise_comparisons:
            name = f"{comp.condition_b} vs {comp.condition_a}"
            sig_marker = "*" if comp.significant else ""
            p_str = f"{comp.p_value:.4f}{sig_marker}"
            pct_str = format_pct(comp.percent_change)
            d_str = f"{comp.cohens_d:.2f}"
            lines.append(
                f"{name:<30} {pct_str:<10} {p_str:<12} "
                f"{d_str:<10} {comp.effect_size_interpretation:<12}"
            )

        lines.append("-" * 80)
        lines.append("* p < 0.05")
        lines.append("")

    # ANOVA
    if result.anova:
        lines.append("One-way ANOVA")
        lines.append("-" * 40)
        for a in result.anova:
            sig = "Yes" if a.significant else "No"
            if a.metric != "default":
                lines.append(f"Metric: {a.metric}")
            lines.append(f"F-statistic: {a.f_statistic:.3f}")
            lines.append(f"p-value:     {a.p_value:.4f}")
            lines.append(f"Significant: {sig} (alpha=0.05)")
        lines.append("")

    # Interpretation
    if result.ranking:
        lines.append("Interpretation")
        lines.append("-" * 60)
        best = result.ranking[0]
        best_cond = _get_cond(best)
        best_val = _get_mean(best_cond)
        lines.append(f"Best: {best} ({metric_label} = {best_val:.4f}{unit_str})")

        if result.control_label and best != result.control_label:
            from polyzymd.compare.statistics import percent_change

            ctrl = _get_cond(result.control_label)
            ctrl_val = _get_mean(ctrl)
            pct = percent_change(ctrl_val, best_val)
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
    higher_is_better: bool,
) -> str:
    """Format as Markdown (internal)."""
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

    for rank, label in enumerate(result.ranking, 1):
        cond = _get_cond(label)
        marker = " (control)" if label == result.control_label else ""
        mean_val = _get_mean(cond)
        sem_val = _get_sem(cond)
        lines.append(
            f"| {rank} | **{cond.label}**{marker} | "
            f"{mean_val:.4f} | {sem_val:.4f} | {cond.n_replicates} |"
        )

    lines.append("")

    # Pairwise
    if result.pairwise_comparisons:
        lines.append("## Statistical Comparisons")
        lines.append("")
        lines.append("| Comparison | % Change | p-value | Cohen's d | Effect | Sig |")
        lines.append("|------------|----------|---------|-----------|--------|-----|")
        for comp in result.pairwise_comparisons:
            name = f"{comp.condition_b} vs {comp.condition_a}"
            sig = "Yes" if comp.significant else "No"
            lines.append(
                f"| {name} | {format_pct(comp.percent_change)} | "
                f"{comp.p_value:.4f} | {comp.cohens_d:.2f} | "
                f"{comp.effect_size_interpretation} | {sig} |"
            )
        lines.append("")

    # ANOVA
    if result.anova:
        lines.append("## ANOVA")
        lines.append("")
        for a in result.anova:
            sig = "Yes" if a.significant else "No"
            lines.append(f"- **F-statistic:** {a.f_statistic:.3f}")
            lines.append(f"- **p-value:** {a.p_value:.4f}")
            lines.append(f"- **Significant:** {sig}")
        lines.append("")

    # Key findings
    if result.ranking:
        lines.append("## Key Findings")
        lines.append("")
        best = result.ranking[0]
        best_cond = _get_cond(best)
        best_val = _get_mean(best_cond)
        lines.append(f"1. **Best condition:** {best} ({metric_label} = {best_val:.4f}{unit_str})")
        if result.control_label and best != result.control_label:
            from polyzymd.compare.statistics import percent_change

            ctrl = _get_cond(result.control_label)
            ctrl_val = _get_mean(ctrl)
            pct = percent_change(ctrl_val, best_val)
            if not math.isnan(pct):
                direction = interpret_direction(pct, ("lower", "unchanged", "higher"))
                if direction == "unchanged":
                    lines.append("2. unchanged relative to control")
                else:
                    magnitude = format_pct(pct).lstrip("+-")
                    lines.append(f"2. {magnitude} {direction} than control")
        lines.append("")

    return "\n".join(lines)
