"""Default comparison mechanics for analysis plugins."""

from __future__ import annotations

import logging
from typing import Any

from pydantic import BaseModel

from polyzymd.analyses._framework.comparison_models import MetricValue
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.mda.artifacts import ConditionArtifact
from polyzymd.analyses.mda.comparison import (
    MDAComparisonContext,
    MDAComparisonError,
    compare_condition_artifacts,
)

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
    summaries_by_condition: dict[str, tuple[Any, Any, str | None]] = {}
    for condition in ctx.conditions:
        summary = ctx.aggregated_results.get(condition.label)
        source: str | None = f"comparison context for {condition.label}"
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
            source = str(analysis.aggregate_result_path(agg_dir))
            summary = analysis._load_aggregated_result(agg_dir)
            if summary is None:
                logger.warning(
                    "%s: missing aggregated result for condition %r — skipping.",
                    analysis.name,
                    condition.label,
                )
                continue

        summaries_by_condition[condition.label] = (condition, summary, source)

    if not summaries_by_condition:
        logger.warning("%s: no conditions have metrics — skipping comparison.", analysis.name)
        return None

    condition_artifacts = [
        summary
        for _, summary, _ in summaries_by_condition.values()
        if isinstance(summary, ConditionArtifact)
    ]
    if condition_artifacts:
        if len(condition_artifacts) != len(summaries_by_condition):
            artifact_labels = [artifact.condition_label for artifact in condition_artifacts]
            noncanonical_labels = [
                label
                for label, (_, summary, _) in summaries_by_condition.items()
                if not isinstance(summary, ConditionArtifact)
            ]
            raise MDAComparisonError(
                f"{analysis.name}: mixed MDA condition artifacts and condition aggregate results "
                f"cannot be compared together; MDA labels={artifact_labels}, "
                f"stale labels={noncanonical_labels}. Recompute this analysis or clear stale "
                "aggregate result.json files before comparing."
            )
        return compare_condition_artifacts(
            condition_artifacts,
            MDAComparisonContext(
                analysis_name=analysis.name,
                project_name=ctx.name,
                expected_condition_labels=[condition.label for condition in ctx.conditions],
                expected_replicates_by_condition={
                    condition.label: condition.replicates for condition in ctx.conditions
                },
                control_label=ctx.control_label,
                effective_control=ctx.effective_control,
                equilibration=ctx.equilibration,
                settings_fingerprint=analysis.aggregate_settings_fingerprint(ctx.settings),
                min_replicates=analysis.min_replicates,
                fdr_alpha=getattr(ctx, "fdr_alpha", 0.05),
                ttest_method=getattr(ctx, "ttest_method", "student"),
                posthoc_method=getattr(ctx, "posthoc_method", "ttest_bh"),
            ),
        )

    from polyzymd.analyses.stats import default_scalar_comparison

    metrics_by_condition: dict[str, dict[str, MetricValue]] = {}
    for condition, summary, source in summaries_by_condition.values():
        summary = analysis.validate_aggregated_result(
            summary,
            condition=condition,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            source=source,
            expected_replicates=condition.replicates,
            allow_replicate_subset=True,
        )
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
