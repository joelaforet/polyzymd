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
from datetime import datetime
from typing import Any

from polyzymd.analyses.base import MetricValue

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
    if abs(pct_change) < threshold:
        return direction_labels[1]
    return direction_labels[0] if pct_change < 0 else direction_labels[2]


# ---------------------------------------------------------------------------
# Pairwise comparisons
# ---------------------------------------------------------------------------


def pairwise_comparisons(
    metrics_by_condition: dict[str, MetricValue],
    control_label: str | None = None,
) -> list[dict[str, Any]]:
    """Compute pairwise statistical comparisons for a single metric.

    Parameters
    ----------
    metrics_by_condition : dict[str, MetricValue]
        Mapping ``condition_label -> MetricValue``.
    control_label : str | None
        If provided, compare all conditions against this control.
        Otherwise, compare all unique pairs.

    Returns
    -------
    list[dict[str, Any]]
        Each dict contains: ``condition_a``, ``condition_b``, ``metric``,
        ``t_statistic``, ``p_value``, ``cohens_d``,
        ``effect_size_interpretation``, ``direction``, ``significant``,
        ``percent_change``.
    """
    from polyzymd.compare.statistics import (
        cohens_d,
        independent_ttest,
        percent_change,
    )

    labels = list(metrics_by_condition.keys())
    results: list[dict[str, Any]] = []

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
        effect = cohens_d(mv_a.replicate_values, mv_b.replicate_values, rmsf_mode=False)
        pct = percent_change(mv_a.mean, mv_b.mean)
        direction = interpret_direction(pct, mv_a.direction_labels)

        results.append(
            {
                "condition_a": label_a,
                "condition_b": label_b,
                "metric": mv_a.name,
                "t_statistic": ttest.t_statistic,
                "p_value": ttest.p_value,
                "cohens_d": effect.cohens_d,
                "effect_size_interpretation": effect.interpretation,
                "direction": direction,
                "significant": ttest.significant,
                "percent_change": pct,
            }
        )

    return results


# ---------------------------------------------------------------------------
# ANOVA
# ---------------------------------------------------------------------------


def anova_test(
    metrics_by_condition: dict[str, MetricValue],
    metric_name: str = "default",
) -> dict[str, Any] | None:
    """Run one-way ANOVA across conditions for a single metric.

    Parameters
    ----------
    metrics_by_condition : dict[str, MetricValue]
        Mapping ``condition_label -> MetricValue``.
    metric_name : str
        Label for the metric in the result dict.

    Returns
    -------
    dict[str, Any] | None
        Dict with ``metric``, ``f_statistic``, ``p_value``, ``significant``.
        ``None`` if fewer than 3 conditions.
    """
    if len(metrics_by_condition) < 3:
        return None

    from polyzymd.compare.statistics import one_way_anova

    groups = [mv.replicate_values for mv in metrics_by_condition.values()]
    result = one_way_anova(*groups)

    return {
        "metric": metric_name,
        "f_statistic": result.f_statistic,
        "p_value": result.p_value,
        "significant": result.significant,
    }


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
) -> dict[str, Any]:
    """Run the standard scalar comparison pipeline.

    This is the default implementation used by :meth:`Analysis.compare`
    for analyses with one or more simple scalar metrics (RMSF, triad,
    secondary structure).

    For each metric:
    1. Pairwise t-tests + Cohen's d + percent-change.
    2. ANOVA (if 3+ conditions).
    3. Ranking.

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
    dict[str, Any]
        Serializable comparison result with keys: ``analysis_type``,
        ``name``, ``control_label``, ``conditions``, ``pairwise_comparisons``,
        ``anova``, ``ranking``, ``equilibration_time``, ``created_at``,
        ``polyzymd_version``.  This can be persisted directly or wrapped
        in an analysis-specific result model.
    """
    from polyzymd import __version__

    # Discover metric names from the first condition
    first_cond = next(iter(metrics_by_condition.values()))
    metric_names = list(first_cond.keys())

    all_pairwise: list[dict[str, Any]] = []
    all_anova: list[dict[str, Any]] = []
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
    condition_summaries = []
    for label, metrics in metrics_by_condition.items():
        summary: dict[str, Any] = {"label": label}
        for metric_name, mv in metrics.items():
            summary[f"{metric_name}_mean"] = mv.mean
            summary[f"{metric_name}_sem"] = mv.sem
            summary[f"{metric_name}_replicate_values"] = mv.replicate_values
        summary["n_replicates"] = (
            len(next(iter(metrics.values())).replicate_values) if metrics else 0
        )
        condition_summaries.append(summary)

    # Use the first metric's ranking as the primary ranking
    primary_ranking = all_rankings.get(metric_names[0], []) if metric_names else []

    return {
        "analysis_type": analysis_name,
        "name": project_name,
        "control_label": control_label,
        "conditions": condition_summaries,
        "pairwise_comparisons": all_pairwise,
        "anova": all_anova if all_anova else None,
        "ranking": primary_ranking,
        "rankings_by_metric": all_rankings if len(metric_names) > 1 else None,
        "equilibration_time": equilibration,
        "created_at": datetime.now().isoformat(),
        "polyzymd_version": __version__,
    }
