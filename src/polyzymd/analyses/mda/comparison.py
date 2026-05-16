"""Comparison engine for MDAnalysis condition artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from polyzymd.analyses._comparison_models import ComparisonResult, MetricValue
from polyzymd.analyses.mda.aggregation import AggregatedMetric
from polyzymd.analyses.mda.artifacts import ComparisonArtifact, ConditionArtifact
from polyzymd.analyses.mda.base import MDAnalysisExtensionError


class MDAComparisonError(MDAnalysisExtensionError):
    """Error raised when condition artifacts cannot be compared."""


@dataclass(frozen=True)
class MDAComparisonContext:
    """Identity and statistical controls for condition-artifact comparison."""

    analysis_name: str
    project_name: str
    expected_condition_labels: Sequence[str] | None = None
    control_label: str | None = None
    effective_control: str | None = None
    equilibration: str = "0ns"
    settings_fingerprint: str | None = None
    fdr_alpha: float = 0.05
    ttest_method: str = "student"
    posthoc_method: str = "ttest_bh"

    def __post_init__(self) -> None:
        """Normalize expected condition labels and reject ambiguous inputs."""

        if self.expected_condition_labels is None:
            return
        labels = tuple(self.expected_condition_labels)
        duplicates = sorted({label for label in labels if labels.count(label) > 1})
        if duplicates:
            raise MDAComparisonError(
                "Expected condition labels must be unique; " f"duplicates: {', '.join(duplicates)}"
            )
        object.__setattr__(self, "expected_condition_labels", labels)


@dataclass(frozen=True)
class _DerivedMetricStatistics:
    """Statistics derived directly from replicate-level metric values."""

    values: list[float]
    mean: float
    std: float
    sem: float


def compare_condition_artifacts(
    artifacts: Sequence[ConditionArtifact],
    ctx: MDAComparisonContext,
) -> ComparisonArtifact:
    """Compare aggregate condition artifacts with replicate-level statistics.

    Parameters
    ----------
    artifacts : sequence of ConditionArtifact
        Condition artifacts produced by MDAnalysis extension-layer aggregation.
    ctx : MDAComparisonContext
        Comparison identity, expected condition labels, and statistical controls.

    Returns
    -------
    ComparisonArtifact
        Stable comparison artifact containing scalar statistics and provenance.
    """

    normalized = _validate_condition_artifacts(artifacts, ctx)
    metrics_by_condition, metric_metadata = _build_metric_values(normalized, ctx)

    from polyzymd.analyses.stats import default_scalar_comparison

    comparison = default_scalar_comparison(
        analysis_name=ctx.analysis_name,
        project_name=ctx.project_name,
        metrics_by_condition=metrics_by_condition,
        control_label=ctx.effective_control,
        equilibration=ctx.equilibration,
        fdr_alpha=ctx.fdr_alpha,
        ttest_method=ctx.ttest_method,
        posthoc_method=ctx.posthoc_method,
    )
    return _comparison_result_to_artifact(comparison, normalized, ctx, metric_metadata)


def _validate_condition_artifacts(
    artifacts: Sequence[ConditionArtifact],
    ctx: MDAComparisonContext,
) -> list[ConditionArtifact]:
    """Validate condition-artifact identity before statistical comparison.

    Parameters
    ----------
    artifacts : sequence of ConditionArtifact
        Candidate artifacts.
    ctx : MDAComparisonContext
        Expected analysis and condition identity.

    Returns
    -------
    list of ConditionArtifact
        Artifacts in input order when validation passes.
    """

    if not artifacts:
        raise MDAComparisonError(f"{ctx.analysis_name}: no condition artifacts to compare")
    normalized = list(artifacts)
    labels = [artifact.condition_label for artifact in normalized]
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise MDAComparisonError(
            f"{ctx.analysis_name}: duplicate condition artifact labels: {', '.join(duplicates)}"
        )

    if ctx.expected_condition_labels is not None:
        expected = list(ctx.expected_condition_labels)
        missing = sorted(set(expected) - set(labels))
        unexpected = sorted(set(labels) - set(expected))
        if missing or unexpected:
            raise MDAComparisonError(
                f"{ctx.analysis_name}: condition artifact labels do not match expected labels; "
                f"missing={missing}, unexpected={unexpected}. Recompute the analysis or clear "
                "stale aggregate result.json files."
            )

    for artifact in normalized:
        _validate_artifact_identity(artifact, ctx)
    return normalized


def _validate_artifact_identity(artifact: ConditionArtifact, ctx: MDAComparisonContext) -> None:
    """Validate one condition artifact against the comparison context.

    Parameters
    ----------
    artifact : ConditionArtifact
        Candidate condition artifact.
    ctx : MDAComparisonContext
        Expected analysis identity.
    """

    if artifact.analysis_name != ctx.analysis_name:
        raise MDAComparisonError(
            f"{ctx.analysis_name}: condition {artifact.condition_label!r} analysis mismatch; "
            f"artifact has {artifact.analysis_name!r}"
        )
    if ctx.settings_fingerprint is not None:
        stored = artifact.metadata.get("settings_fingerprint")
        if stored != ctx.settings_fingerprint:
            raise MDAComparisonError(
                f"{ctx.analysis_name}: condition {artifact.condition_label!r} settings "
                f"fingerprint mismatch; expected {ctx.settings_fingerprint!r}, got {stored!r}. "
                "Recompute the analysis or clear stale aggregate result.json files."
            )
    if len(set(artifact.replicates)) != len(artifact.replicates):
        raise MDAComparisonError(
            f"{ctx.analysis_name}: condition {artifact.condition_label!r} has duplicate replicates"
        )


def _build_metric_values(
    artifacts: Sequence[ConditionArtifact],
    ctx: MDAComparisonContext,
) -> tuple[dict[str, dict[str, MetricValue]], dict[str, dict[str, Any]]]:
    """Build ``MetricValue`` inputs from condition artifacts.

    Parameters
    ----------
    artifacts : sequence of ConditionArtifact
        Validated condition artifacts.
    ctx : MDAComparisonContext
        Comparison context used for diagnostics.

    Returns
    -------
    tuple of dict
        Metrics keyed by condition label and metric metadata keyed by metric name.
    """

    metrics_by_condition: dict[str, dict[str, MetricValue]] = {}
    metric_metadata: dict[str, dict[str, Any]] = {}
    expected_metric_keys: set[str] | None = None
    baseline_label: str | None = None
    key_mismatches: list[str] = []

    for artifact in artifacts:
        metric_payload = artifact.payload.get("metrics")
        if not isinstance(metric_payload, Mapping):
            raise MDAComparisonError(
                f"{ctx.analysis_name}: condition {artifact.condition_label!r} is missing "
                "payload['metrics'] mapping"
            )
        condition_metrics: dict[str, MetricValue] = {}
        metric_keys: set[str] = set()
        for metric_key, raw_metric in metric_payload.items():
            if not isinstance(metric_key, str):
                raise MDAComparisonError(
                    f"{ctx.analysis_name}: condition {artifact.condition_label!r} has non-string "
                    f"metric key {metric_key!r}"
                )
            metric_keys.add(metric_key)
            summary = _validated_aggregated_metric(raw_metric, artifact, ctx, metric_key)
            metadata = _metric_metadata_for(artifact, metric_key, raw_metric)
            _merge_metric_metadata(metric_metadata, metric_key, metadata, artifact, ctx)
            condition_metrics[metric_key] = MetricValue(
                name=metric_key,
                mean=summary.mean,
                sem=summary.sem,
                replicate_values=summary.values,
                higher_is_better=metadata.get("higher_is_better", True),
                direction_labels=metadata.get(
                    "direction_labels", ("decreased", "unchanged", "increased")
                ),
            )
        if not condition_metrics:
            raise MDAComparisonError(
                f"{ctx.analysis_name}: condition {artifact.condition_label!r} declares no metrics"
            )
        if expected_metric_keys is None:
            expected_metric_keys = metric_keys
            baseline_label = artifact.condition_label
        elif metric_keys != expected_metric_keys:
            missing = sorted(expected_metric_keys - metric_keys)
            extra = sorted(metric_keys - expected_metric_keys)
            key_mismatches.append(f"{artifact.condition_label}: missing={missing}, extra={extra}")
        metrics_by_condition[artifact.condition_label] = condition_metrics

    if key_mismatches:
        raise MDAComparisonError(
            f"{ctx.analysis_name}: inconsistent metric keys across condition artifacts; "
            f"baseline {baseline_label!r} keys={sorted(expected_metric_keys or [])}; "
            f"differences: {'; '.join(key_mismatches)}"
        )
    return metrics_by_condition, metric_metadata


def _validated_aggregated_metric(
    raw_metric: Any,
    artifact: ConditionArtifact,
    ctx: MDAComparisonContext,
    metric_key: str,
) -> _DerivedMetricStatistics:
    """Validate one raw metric payload with the aggregation schema.

    Parameters
    ----------
    raw_metric : Any
        Candidate ``AggregatedMetric`` payload.
    artifact : ConditionArtifact
        Source condition artifact.
    ctx : MDAComparisonContext
        Comparison context used for diagnostics.
    metric_key : str
        Metric key within ``payload['metrics']``.

    Returns
    -------
    _DerivedMetricStatistics
        Metric statistics derived from validated replicate values.
    """

    if not isinstance(raw_metric, Mapping):
        raise MDAComparisonError(
            f"{ctx.analysis_name}: condition {artifact.condition_label!r} metric {metric_key!r} "
            f"must be a mapping, got {type(raw_metric).__name__}"
        )
    raw_values = raw_metric.get("values")
    if not isinstance(raw_values, Sequence) or isinstance(raw_values, (str, bytes, bytearray)):
        raise MDAComparisonError(
            f"{ctx.analysis_name}: condition {artifact.condition_label!r} metric {metric_key!r} "
            "must provide a sequence of replicate values"
        )
    for index, value in enumerate(raw_values):
        _validate_finite_scalar(
            value,
            ctx=ctx,
            condition_label=artifact.condition_label,
            metric_key=metric_key,
            value_name=f"values[{index}]",
        )
    try:
        summary = AggregatedMetric.model_validate(raw_metric)
    except ValueError as exc:
        raise MDAComparisonError(
            f"{ctx.analysis_name}: condition {artifact.condition_label!r} metric {metric_key!r} "
            f"does not match AggregatedMetric schema: {exc}"
        ) from exc

    for value_name in ("mean", "sem", "std"):
        _validate_finite_scalar(
            getattr(summary, value_name),
            ctx=ctx,
            condition_label=artifact.condition_label,
            metric_key=metric_key,
            value_name=value_name,
        )
    if summary.name != metric_key:
        raise MDAComparisonError(
            f"{ctx.analysis_name}: condition {artifact.condition_label!r} metric key {metric_key!r} "
            f"does not match AggregatedMetric.name {summary.name!r}"
        )
    if summary.n != len(summary.values) or summary.n != len(artifact.replicates):
        raise MDAComparisonError(
            f"{ctx.analysis_name}: condition {artifact.condition_label!r} metric {metric_key!r} "
            "has inconsistent replicate counts; expected n == len(values) == len(replicates), "
            f"got n={summary.n}, len(values)={len(summary.values)}, "
            f"len(replicates)={len(artifact.replicates)}"
        )
    derived = _derive_metric_statistics(summary.values)
    _validate_stored_statistic_matches(
        summary.mean,
        derived.mean,
        statistic_name="mean",
        artifact=artifact,
        ctx=ctx,
        metric_key=metric_key,
    )
    _validate_stored_statistic_matches(
        summary.std,
        derived.std,
        statistic_name="std",
        artifact=artifact,
        ctx=ctx,
        metric_key=metric_key,
    )
    _validate_stored_statistic_matches(
        summary.sem,
        derived.sem,
        statistic_name="sem",
        artifact=artifact,
        ctx=ctx,
        metric_key=metric_key,
    )
    return derived


def _derive_metric_statistics(values: Sequence[float]) -> _DerivedMetricStatistics:
    """Derive comparison statistics from replicate-level values.

    Parameters
    ----------
    values : sequence of float
        Replicate-level metric values from an ``AggregatedMetric`` payload.

    Returns
    -------
    _DerivedMetricStatistics
        Mean, sample standard deviation, and SEM derived from ``values``.
    """

    derived_values = [float(value) for value in values]
    mean = sum(derived_values) / len(derived_values)
    if len(derived_values) == 1:
        std = 0.0
    else:
        std = math.sqrt(
            sum((value - mean) ** 2 for value in derived_values) / (len(derived_values) - 1)
        )
    sem = std / math.sqrt(len(derived_values))
    return _DerivedMetricStatistics(values=derived_values, mean=mean, std=std, sem=sem)


def _validate_stored_statistic_matches(
    stored: float,
    derived: float,
    *,
    statistic_name: str,
    artifact: ConditionArtifact,
    ctx: MDAComparisonContext,
    metric_key: str,
) -> None:
    """Reject stale aggregate statistics that disagree with replicate values.

    Parameters
    ----------
    stored : float
        Statistic stored in the condition artifact.
    derived : float
        Statistic recalculated from ``AggregatedMetric.values``.
    statistic_name : str
        Name used in diagnostics.
    artifact : ConditionArtifact
        Source condition artifact.
    ctx : MDAComparisonContext
        Comparison context used for diagnostics.
    metric_key : str
        Metric key.
    """

    if math.isclose(stored, derived, rel_tol=1e-9, abs_tol=1e-12):
        return
    raise MDAComparisonError(
        f"{ctx.analysis_name}: condition {artifact.condition_label!r} metric {metric_key!r} "
        f"stored {statistic_name}={stored!r} does not match value-derived "
        f"{statistic_name}={derived!r}; recompute the analysis or clear stale aggregate "
        "result.json files"
    )


def _validate_finite_scalar(
    value: Any,
    *,
    ctx: MDAComparisonContext,
    condition_label: str,
    metric_key: str,
    value_name: str,
) -> None:
    """Validate that a metric payload value is a finite scalar.

    Parameters
    ----------
    value : Any
        Candidate scalar value.
    ctx : MDAComparisonContext
        Comparison context used for diagnostics.
    condition_label : str
        Source condition label.
    metric_key : str
        Metric key.
    value_name : str
        Field name used in diagnostics.
    """

    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise MDAComparisonError(
            f"{ctx.analysis_name}: condition {condition_label!r} metric {metric_key!r} "
            f"has non-finite scalar {value_name}={value!r}"
        )


def _metric_metadata_for(
    artifact: ConditionArtifact,
    metric_key: str,
    raw_metric: Any,
) -> dict[str, Any]:
    """Extract optional comparison metadata for one metric.

    Parameters
    ----------
    artifact : ConditionArtifact
        Source condition artifact.
    metric_key : str
        Metric key.
    raw_metric : Any
        Raw metric payload mapping.

    Returns
    -------
    dict[str, Any]
        Normalized metadata for ``MetricValue`` construction and artifact output.
    """

    metadata: dict[str, Any] = {}
    payload_metadata = artifact.payload.get("metric_metadata")
    if isinstance(payload_metadata, Mapping):
        metric_metadata = payload_metadata.get(metric_key)
        if isinstance(metric_metadata, Mapping):
            metadata.update(metric_metadata)
    if isinstance(raw_metric, Mapping):
        for key in ("higher_is_better", "direction_labels", "label", "unit"):
            if key in raw_metric:
                metadata[key] = raw_metric[key]
    if "higher_is_better" in metadata and not isinstance(
        metadata["higher_is_better"], (bool, type(None))
    ):
        raise MDAComparisonError(
            f"{artifact.analysis_name}: metric {metric_key!r} higher_is_better must be "
            "bool or None"
        )
    if "direction_labels" in metadata:
        labels = metadata["direction_labels"]
        if (
            not isinstance(labels, Sequence)
            or isinstance(labels, (str, bytes, bytearray))
            or len(labels) != 3
            or any(not isinstance(label, str) for label in labels)
        ):
            raise MDAComparisonError(
                f"{artifact.analysis_name}: metric {metric_key!r} direction_labels must be "
                "three strings"
            )
        metadata["direction_labels"] = tuple(labels)
    metadata.setdefault("higher_is_better", True)
    metadata.setdefault("direction_labels", ("decreased", "unchanged", "increased"))
    return metadata


def _merge_metric_metadata(
    metric_metadata: dict[str, dict[str, Any]],
    metric_key: str,
    metadata: Mapping[str, Any],
    artifact: ConditionArtifact,
    ctx: MDAComparisonContext,
) -> None:
    """Record metric metadata and reject cross-condition disagreements.

    Parameters
    ----------
    metric_metadata : dict[str, dict[str, Any]]
        Accumulated metadata by metric key.
    metric_key : str
        Metric key being merged.
    metadata : mapping
        Metadata for the current condition.
    artifact : ConditionArtifact
        Source condition artifact.
    ctx : MDAComparisonContext
        Comparison context used for diagnostics.
    """

    normalized = dict(metadata)
    if metric_key not in metric_metadata:
        metric_metadata[metric_key] = normalized
        return
    previous = metric_metadata[metric_key]
    for key in ("higher_is_better", "direction_labels"):
        if previous.get(key) != normalized.get(key):
            raise MDAComparisonError(
                f"{ctx.analysis_name}: condition {artifact.condition_label!r} metric {metric_key!r} "
                f"metadata mismatch for {key}; expected {previous.get(key)!r}, "
                f"got {normalized.get(key)!r}"
            )


def _comparison_result_to_artifact(
    comparison: ComparisonResult,
    artifacts: Sequence[ConditionArtifact],
    ctx: MDAComparisonContext,
    metric_metadata: Mapping[str, Mapping[str, Any]],
) -> ComparisonArtifact:
    """Convert the scalar comparison model into an MDA comparison artifact.

    Parameters
    ----------
    comparison : ComparisonResult
        Result from the shared scalar statistics pipeline.
    artifacts : sequence of ConditionArtifact
        Source condition artifacts.
    ctx : MDAComparisonContext
        Comparison context.
    metric_metadata : mapping
        Metadata keyed by metric name.

    Returns
    -------
    ComparisonArtifact
        Serialized artifact envelope for cross-condition comparison.
    """

    condition_labels = [artifact.condition_label for artifact in artifacts]
    payload = {
        "condition_summaries": [condition.model_dump() for condition in comparison.conditions],
        "pairwise_comparisons": [
            pairwise.model_dump() for pairwise in comparison.pairwise_comparisons
        ],
        "anova": [result.model_dump() for result in comparison.anova or []],
        "ranking": comparison.ranking,
        "rankings_by_metric": comparison.rankings_by_metric,
        "metric_metadata": _json_metric_metadata(metric_metadata),
        "statistical_parameters": {
            "fdr_alpha": ctx.fdr_alpha,
            "ttest_method": ctx.ttest_method,
            "posthoc_method": ctx.posthoc_method,
            "control_label": ctx.control_label,
            "effective_control": ctx.effective_control,
            "equilibration": ctx.equilibration,
        },
    }
    metadata: dict[str, Any] = {
        "n_conditions": len(condition_labels),
        "metrics": list(metric_metadata.keys()),
        "comparison_result_model": type(comparison).__name__,
    }
    if ctx.settings_fingerprint is not None:
        metadata["settings_fingerprint"] = ctx.settings_fingerprint

    return ComparisonArtifact(
        analysis_name=ctx.analysis_name,
        conditions=condition_labels,
        control_label=ctx.control_label,
        effective_control=ctx.effective_control,
        payload=payload,
        provenance={
            "source": "mda_condition_artifact_comparison",
            "source_condition_labels": condition_labels,
            "source_replicates": {
                artifact.condition_label: list(artifact.replicates) for artifact in artifacts
            },
            "source_condition_artifacts": [
                {
                    "condition_label": artifact.condition_label,
                    "replicates": list(artifact.replicates),
                    "source_replicates": list(artifact.source_replicates),
                    "skipped_replicates": list(artifact.skipped_replicates),
                }
                for artifact in artifacts
            ],
            "skipped_replicates": {
                artifact.condition_label: list(artifact.skipped_replicates)
                for artifact in artifacts
            },
        },
        metadata=metadata,
        warnings=_combined_warnings(artifacts),
    )


def _json_metric_metadata(
    metric_metadata: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Convert metric metadata to JSON-friendly dictionaries.

    Parameters
    ----------
    metric_metadata : mapping
        Metadata keyed by metric name.

    Returns
    -------
    dict[str, dict[str, Any]]
        JSON-friendly metadata.
    """

    converted: dict[str, dict[str, Any]] = {}
    for metric_key, metadata in metric_metadata.items():
        converted[metric_key] = {
            key: list(value) if isinstance(value, tuple) else value
            for key, value in metadata.items()
        }
    return converted


def _combined_warnings(artifacts: Sequence[ConditionArtifact]) -> list[str]:
    """Return de-duplicated warnings from source condition artifacts.

    Parameters
    ----------
    artifacts : sequence of ConditionArtifact
        Source artifacts.

    Returns
    -------
    list[str]
        Warnings in first-seen order.
    """

    warnings: list[str] = []
    seen: set[str] = set()
    for artifact in artifacts:
        for warning in artifact.warnings:
            if warning in seen:
                continue
            warnings.append(warning)
            seen.add(warning)
    return warnings
