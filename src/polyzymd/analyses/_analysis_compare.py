"""Default comparison mechanics for analysis plugins."""

from __future__ import annotations

import logging
from typing import Any

from pydantic import BaseModel

from polyzymd.analyses._comparison_models import MetricValue
from polyzymd.analyses.exceptions import PluginContractError

logger = logging.getLogger("polyzymd.analyses")


def default_compare(analysis: Any, ctx: Any) -> BaseModel | None:
    """Compare aggregated results across conditions with the scalar pipeline.

    Parameters
    ----------
    analysis : Any
        Analysis instance.
    ctx : Any
        Comparison context.

    Returns
    -------
    BaseModel | None
        Comparison result, or ``None`` if no metrics are available.
    """
    from polyzymd.analyses.stats import default_scalar_comparison

    metrics_by_condition: dict[str, dict[str, MetricValue]] = {}
    for condition in ctx.conditions:
        summary = ctx.aggregated_results.get(condition.label)
        if summary is None:
            agg_dir_parent = ctx.analysis_dirs.get(condition.label)
            if agg_dir_parent is None:
                logger.warning(
                    "%s: no analysis directory for condition %r — skipping.",
                    analysis.name,
                    condition.label,
                )
                continue
            agg_dir = agg_dir_parent / "aggregated"
            summary = analysis._load_aggregated_result(agg_dir)
            if summary is None:
                logger.warning(
                    "%s: missing aggregated result for condition %r — skipping.",
                    analysis.name,
                    condition.label,
                )
                continue

        extracted = analysis.extract_metrics(summary)
        if not isinstance(extracted, dict):
            raise PluginContractError(
                f"Plugin '{analysis.name}' extract_metrics() must return dict[str, MetricValue] "
                f"for condition '{condition.label}', got {type(extracted).__name__}"
            )
        for metric_key, metric_value in extracted.items():
            if not isinstance(metric_value, MetricValue):
                raise PluginContractError(
                    f"Plugin '{analysis.name}' extract_metrics() returned invalid value for "
                    f"key '{metric_key}' in condition '{condition.label}': expected MetricValue, "
                    f"got {type(metric_value).__name__}"
                )
        if not extracted:
            raise PluginContractError(
                f"Plugin '{analysis.name}' extract_metrics() returned empty dict for "
                f"condition '{condition.label}' — implement extract_metrics() or override compare()"
            )
        metrics_by_condition[condition.label] = extracted

    if not metrics_by_condition:
        logger.warning("%s: no conditions have metrics — skipping comparison.", analysis.name)
        return None

    return default_scalar_comparison(
        analysis_name=analysis.name,
        project_name=ctx.name,
        metrics_by_condition=metrics_by_condition,
        control_label=ctx.effective_control,
        equilibration=ctx.equilibration,
        fdr_alpha=ctx.fdr_alpha,
        ttest_method=ctx.ttest_method,
        posthoc_method=ctx.posthoc_method,
    )
