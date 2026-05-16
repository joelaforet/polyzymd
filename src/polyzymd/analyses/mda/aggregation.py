"""Condition aggregation for MDAnalysis replicate artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Protocol

from pydantic import BaseModel, Field

from polyzymd.analyses.mda.artifacts import ConditionArtifact, ReplicateArtifact
from polyzymd.analyses.mda.base import MDAnalysisExtensionError
from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError


class MDAAggregationError(MDAnalysisExtensionError):
    """Error raised when MDAnalysis replicate artifacts cannot be aggregated."""


class AggregatedMetric(BaseModel):
    """Summary statistics for one metric across biological replicates."""

    name: str
    values: list[float] = Field(default_factory=list)
    mean: float
    sem: float
    std: float
    n: int = Field(ge=1)


class ReplicateMetricPolicy(Protocol):
    """Protocol for reducing one replicate artifact to scalar metrics."""

    def extract_metrics(self, artifact: ReplicateArtifact) -> Mapping[str, float]:
        """Extract one scalar value per metric from a replicate artifact.

        Parameters
        ----------
        artifact : ReplicateArtifact
            Artifact produced for one replicate.

        Returns
        -------
        Mapping[str, float]
            Metric names mapped to one replicate-level scalar each.
        """


class ExplicitReplicateMetricPolicy:
    """Extract explicitly declared replicate-level scalar metrics.

    The default policy deliberately reads only ``payload["metrics"]`` or
    ``payload["replicate_metrics"]``. It does not reduce arrays, events, job
    tables, or frame-level values because those reductions are analysis-specific
    scientific choices.
    """

    def extract_metrics(self, artifact: ReplicateArtifact) -> Mapping[str, float]:
        """Return validated scalar metrics from a replicate artifact.

        Parameters
        ----------
        artifact : ReplicateArtifact
            Artifact produced for one replicate.

        Returns
        -------
        Mapping[str, float]
            Validated replicate-level scalar metrics.
        """

        payload = artifact.payload
        metric_payload = payload.get("metrics")
        if metric_payload is None:
            metric_payload = payload.get("replicate_metrics")
        if metric_payload is None:
            raise MDAAggregationError(
                f"{artifact.analysis_name}: replicate {artifact.replicate} is missing explicit "
                "payload['metrics'] or payload['replicate_metrics'] scalar metrics"
            )
        if not isinstance(metric_payload, Mapping):
            raise MDAAggregationError(
                f"{artifact.analysis_name}: replicate {artifact.replicate} metrics must be a mapping, "
                f"got {type(metric_payload).__name__}"
            )
        metrics: dict[str, float] = {}
        for name, value in metric_payload.items():
            if not isinstance(name, str):
                raise MDAAggregationError(
                    f"{artifact.analysis_name}: replicate {artifact.replicate} has non-string "
                    f"metric name {name!r}"
                )
            metrics[name] = _validate_metric_scalar(
                value,
                analysis_name=artifact.analysis_name,
                replicate=artifact.replicate,
                metric_name=name,
            )
        if not metrics:
            raise MDAAggregationError(
                f"{artifact.analysis_name}: replicate {artifact.replicate} declares no metrics"
            )
        return metrics


@dataclass(frozen=True)
class MDAAggregationContext:
    """Identity and provenance expected during condition aggregation."""

    analysis_name: str
    condition_label: str
    expected_replicates: tuple[int, ...]
    settings_fingerprint: str | None = None
    min_replicates: int = 1
    allow_partial: bool = False
    require_compatible_frame_selection: bool = True
    expected_frame_selection: Mapping[str, Any] | None = None
    validate_sidecars: bool = True
    artifact_stores: Mapping[int, ArtifactStore] = field(default_factory=dict)
    source_replicates: Sequence[Mapping[str, Any]] = ()
    skipped_replicates: Sequence[Mapping[str, Any]] = ()

    def __post_init__(self) -> None:
        """Normalize replicate identity and validate minimum count."""

        expected = tuple(int(rep) for rep in self.expected_replicates)
        if len(set(expected)) != len(expected):
            raise MDAAggregationError("Expected replicate list contains duplicates")
        if self.min_replicates < 1:
            raise MDAAggregationError("min_replicates must be at least 1")
        object.__setattr__(self, "expected_replicates", expected)
        object.__setattr__(
            self,
            "source_replicates",
            tuple(dict(entry) for entry in self.source_replicates),
        )
        object.__setattr__(
            self,
            "skipped_replicates",
            tuple(dict(entry) for entry in self.skipped_replicates),
        )


def aggregate_replicate_artifacts(
    artifacts: Sequence[ReplicateArtifact],
    ctx: MDAAggregationContext,
    policy: ReplicateMetricPolicy | None = None,
) -> ConditionArtifact:
    """Aggregate replicate artifacts into a condition artifact.

    Parameters
    ----------
    artifacts : sequence of ReplicateArtifact
        Replicate artifacts to aggregate.
    ctx : MDAAggregationContext
        Expected condition identity and provenance.
    policy : ReplicateMetricPolicy or None, optional
        Metric extraction policy, by default ``ExplicitReplicateMetricPolicy``.

    Returns
    -------
    ConditionArtifact
        Aggregated condition artifact containing replicate-level statistics.
    """

    metric_policy = policy or ExplicitReplicateMetricPolicy()
    normalized = _validate_artifact_set(artifacts, ctx)
    frame_selection = _validate_artifact_provenance(normalized, ctx)
    replicate_metrics = _extract_replicate_metrics(normalized, metric_policy)
    metric_summaries = _summarize_metrics(replicate_metrics, ctx)
    replicate_ids = [artifact.replicate for artifact in normalized]
    source_replicates = list(ctx.source_replicates) or _source_replicates_from_artifacts(normalized)
    skipped_replicates = _skipped_replicates(normalized, ctx)
    metadata: dict[str, Any] = {
        "n_replicates": len(replicate_ids),
        "metric_policy": type(metric_policy).__name__,
    }
    if ctx.settings_fingerprint is not None:
        metadata["settings_fingerprint"] = ctx.settings_fingerprint
    return ConditionArtifact(
        analysis_name=ctx.analysis_name,
        condition_label=ctx.condition_label,
        replicates=replicate_ids,
        payload={
            "metrics": {name: summary.model_dump() for name, summary in metric_summaries.items()},
            "replicate_metrics": {
                str(replicate): dict(metrics) for replicate, metrics in replicate_metrics.items()
            },
            "n_replicates": len(replicate_ids),
        },
        provenance={
            "source": "mda_replicate_artifact_aggregation",
            "frame_selection": frame_selection,
        },
        metadata=metadata,
        source_replicates=source_replicates,
        skipped_replicates=skipped_replicates,
        warnings=_combined_warnings(normalized),
    )


def aggregate_replicate_artifacts_from_disk(
    analysis_dir: Path,
    ctx: MDAAggregationContext,
    policy: ReplicateMetricPolicy | None = None,
    *,
    artifact_path: str | Path = "result.json",
) -> ConditionArtifact:
    """Load replicate artifacts from disk and aggregate them.

    Parameters
    ----------
    analysis_dir : Path
        Condition analysis directory containing ``run_N`` subdirectories.
    ctx : MDAAggregationContext
        Expected condition identity and aggregation policy controls.
    policy : ReplicateMetricPolicy or None, optional
        Optional custom metric extraction policy.
    artifact_path : str or Path, optional
        Store-relative replicate artifact filename, by default ``"result.json"``.

    Returns
    -------
    ConditionArtifact
        Aggregated condition artifact.
    """

    artifacts: list[ReplicateArtifact] = []
    stores: dict[int, ArtifactStore] = {}
    sources: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = [dict(entry) for entry in ctx.skipped_replicates]
    for replicate in ctx.expected_replicates:
        run_dir = analysis_dir / f"run_{replicate}"
        store = ArtifactStore(run_dir)
        result_path = run_dir / artifact_path
        if not result_path.exists():
            _record_or_raise_skip(
                ctx,
                skipped,
                replicate=replicate,
                reason="missing artifact",
                path=result_path,
            )
            continue
        try:
            artifact = store.read_replicate_result(artifact_path)
            source = store.source_artifact_ref(artifact_path)
        except ArtifactStoreError as exc:
            _record_or_raise_skip(
                ctx,
                skipped,
                replicate=replicate,
                reason=f"stale or invalid artifact: {exc}",
                path=result_path,
            )
            continue
        artifacts.append(artifact)
        stores[replicate] = store
        sources.append({"replicate": replicate, "artifact": source})
    disk_ctx = replace(
        ctx, artifact_stores=stores, source_replicates=sources, skipped_replicates=skipped
    )
    return aggregate_replicate_artifacts(artifacts, disk_ctx, policy)


def _validate_metric_scalar(
    value: Any,
    *,
    analysis_name: str,
    replicate: int,
    metric_name: str,
) -> float:
    """Validate one replicate-level metric scalar."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise MDAAggregationError(
            f"{analysis_name}: replicate {replicate} metric {metric_name!r} must be one "
            f"finite scalar per replicate, got {type(value).__name__}"
        )
    scalar = float(value)
    if not math.isfinite(scalar):
        raise MDAAggregationError(
            f"{analysis_name}: replicate {replicate} metric {metric_name!r} is non-finite"
        )
    return scalar


def _validate_artifact_set(
    artifacts: Sequence[ReplicateArtifact],
    ctx: MDAAggregationContext,
) -> list[ReplicateArtifact]:
    """Validate artifact identity, duplicate, missing, and sidecar provenance."""

    expected = set(ctx.expected_replicates)
    seen: set[int] = set()
    normalized: list[ReplicateArtifact] = []
    for artifact in artifacts:
        if not isinstance(artifact, ReplicateArtifact):
            raise MDAAggregationError(
                f"{ctx.analysis_name}: expected ReplicateArtifact, got {type(artifact).__name__}"
            )
        if artifact.analysis_name != ctx.analysis_name:
            raise MDAAggregationError(
                f"Artifact analysis mismatch for replicate {artifact.replicate}: "
                f"stored {artifact.analysis_name!r}, expected {ctx.analysis_name!r}"
            )
        if artifact.condition_label != ctx.condition_label:
            raise MDAAggregationError(
                f"Artifact condition mismatch for replicate {artifact.replicate}: "
                f"stored {artifact.condition_label!r}, expected {ctx.condition_label!r}"
            )
        if artifact.replicate not in expected:
            raise MDAAggregationError(
                f"Unexpected replicate artifact {artifact.replicate}; expected "
                f"{list(ctx.expected_replicates)}"
            )
        if artifact.replicate in seen:
            raise MDAAggregationError(f"Duplicate replicate artifact {artifact.replicate}")
        seen.add(artifact.replicate)
        _validate_settings_fingerprint(artifact, ctx)
        _validate_artifact_sidecars(artifact, ctx)
        normalized.append(artifact)
    missing = [replicate for replicate in ctx.expected_replicates if replicate not in seen]
    if missing and not ctx.allow_partial:
        raise MDAAggregationError(
            f"{ctx.analysis_name}: missing replicate artifact(s) for {missing}; "
            "set allow_partial=True to aggregate available replicates"
        )
    if len(normalized) < ctx.min_replicates:
        raise MDAAggregationError(
            f"{ctx.analysis_name}: only {len(normalized)} replicate artifact(s) available, "
            f"need at least {ctx.min_replicates}"
        )
    return sorted(normalized, key=lambda artifact: artifact.replicate)


def _validate_settings_fingerprint(artifact: ReplicateArtifact, ctx: MDAAggregationContext) -> None:
    """Validate artifact settings identity when expected."""

    if ctx.settings_fingerprint is None:
        return
    stored = artifact.metadata.get("settings_fingerprint")
    if stored is None:
        stored = artifact.metadata.get("settings_fp")
    if stored is None:
        raise MDAAggregationError(
            f"{ctx.analysis_name}: replicate {artifact.replicate} is missing settings fingerprint"
        )
    if str(stored) != str(ctx.settings_fingerprint):
        raise MDAAggregationError(
            f"{ctx.analysis_name}: settings fingerprint mismatch for replicate {artifact.replicate}: "
            f"stored {stored}, expected {ctx.settings_fingerprint}"
        )


def _validate_artifact_sidecars(artifact: ReplicateArtifact, ctx: MDAAggregationContext) -> None:
    """Validate sidecars using the per-replicate artifact store."""

    if not ctx.validate_sidecars:
        return
    if not artifact.sidecars:
        return
    store = ctx.artifact_stores.get(artifact.replicate)
    if store is None:
        raise MDAAggregationError(
            f"{ctx.analysis_name}: replicate {artifact.replicate} has sidecars but no "
            "ArtifactStore was provided for validation"
        )
    for sidecar in artifact.sidecars:
        try:
            store.validate_sidecar(sidecar)
        except ArtifactStoreError as exc:
            raise MDAAggregationError(
                f"{ctx.analysis_name}: stale sidecar for replicate {artifact.replicate}: {exc}"
            ) from exc


def _validate_artifact_provenance(
    artifacts: Sequence[ReplicateArtifact],
    ctx: MDAAggregationContext,
) -> dict[str, Any] | None:
    """Validate compatible frame-selection provenance across replicates."""

    if not ctx.require_compatible_frame_selection:
        return (
            dict(ctx.expected_frame_selection) if ctx.expected_frame_selection is not None else None
        )
    expected = (
        dict(ctx.expected_frame_selection) if ctx.expected_frame_selection is not None else None
    )
    for artifact in artifacts:
        frame_selection = artifact.provenance.get("frame_selection")
        if not isinstance(frame_selection, Mapping):
            raise MDAAggregationError(
                f"{ctx.analysis_name}: replicate {artifact.replicate} is missing "
                "frame-selection provenance"
            )
        frame_payload = dict(frame_selection)
        if expected is None:
            expected = frame_payload
            continue
        if frame_payload != expected:
            raise MDAAggregationError(
                f"{ctx.analysis_name}: incompatible frame-selection provenance for "
                f"replicate {artifact.replicate}"
            )
    return expected


def _extract_replicate_metrics(
    artifacts: Sequence[ReplicateArtifact],
    policy: ReplicateMetricPolicy,
) -> dict[int, dict[str, float]]:
    """Extract and validate identical metric sets from all artifacts."""

    extracted: dict[int, dict[str, float]] = {}
    expected_names: tuple[str, ...] | None = None
    for artifact in artifacts:
        metrics = dict(policy.extract_metrics(artifact))
        for name, value in metrics.items():
            metrics[name] = _validate_metric_scalar(
                value,
                analysis_name=artifact.analysis_name,
                replicate=artifact.replicate,
                metric_name=name,
            )
        names = tuple(sorted(metrics))
        if expected_names is None:
            expected_names = names
        elif names != expected_names:
            raise MDAAggregationError(
                f"{artifact.analysis_name}: metric set mismatch for replicate {artifact.replicate}: "
                f"stored {list(names)}, expected {list(expected_names)}"
            )
        extracted[artifact.replicate] = metrics
    return extracted


def _summarize_metrics(
    replicate_metrics: Mapping[int, Mapping[str, float]],
    ctx: MDAAggregationContext,
) -> dict[str, AggregatedMetric]:
    """Compute mean, sample standard deviation, and SEM per metric."""

    metric_names = sorted(next(iter(replicate_metrics.values())).keys())
    summaries: dict[str, AggregatedMetric] = {}
    for name in metric_names:
        values = [
            float(replicate_metrics[replicate][name]) for replicate in sorted(replicate_metrics)
        ]
        mean = sum(values) / len(values)
        if len(values) == 1:
            std = 0.0
        else:
            std = math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))
        sem = std / math.sqrt(len(values))
        summaries[name] = AggregatedMetric(
            name=name,
            values=values,
            mean=mean,
            sem=sem,
            std=std,
            n=len(values),
        )
    if not summaries:
        raise MDAAggregationError(f"{ctx.analysis_name}: no metrics available to aggregate")
    return summaries


def _source_replicates_from_artifacts(
    artifacts: Sequence[ReplicateArtifact],
) -> list[dict[str, Any]]:
    """Build source-replicate records when no file hashes are available."""

    return [{"replicate": artifact.replicate} for artifact in artifacts]


def _skipped_replicates(
    artifacts: Sequence[ReplicateArtifact],
    ctx: MDAAggregationContext,
) -> list[dict[str, Any]]:
    """Return explicit skipped-replicate provenance for partial aggregation."""

    skipped = [dict(entry) for entry in ctx.skipped_replicates]
    if not ctx.allow_partial:
        return skipped
    existing = {int(entry["replicate"]) for entry in skipped if "replicate" in entry}
    present = {artifact.replicate for artifact in artifacts}
    for replicate in ctx.expected_replicates:
        if replicate not in present and replicate not in existing:
            skipped.append({"replicate": replicate, "reason": "missing artifact"})
    return skipped


def _combined_warnings(artifacts: Sequence[ReplicateArtifact]) -> list[str]:
    """Combine replicate warnings while preserving first-seen order."""

    warnings: list[str] = []
    seen: set[str] = set()
    for artifact in artifacts:
        for warning in artifact.warnings:
            text = str(warning)
            if text not in seen:
                seen.add(text)
                warnings.append(text)
    return warnings


def _record_or_raise_skip(
    ctx: MDAAggregationContext,
    skipped: list[dict[str, Any]],
    *,
    replicate: int,
    reason: str,
    path: Path,
) -> None:
    """Record a skipped replicate or raise for default fail-fast behavior."""

    if not ctx.allow_partial:
        raise MDAAggregationError(
            f"{ctx.analysis_name}: cannot aggregate replicate {replicate} from {path}: {reason}"
        )
    skipped.append({"replicate": replicate, "reason": reason, "path": str(path)})
