"""MDAnalysis-native catalytic-triad jobs and artifact adapters."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses._results_base import get_polyzymd_version
from polyzymd.analyses.catalytic_triad._results import TriadAggregatedResult, TriadResult
from polyzymd.analyses.distances._mda import (
    DistancesRunnerPayload,
    _pair_payload_to_json,
    mdanalysis_version,
    payload_from_pair_distance_analysis,
    resolve_distance_pairs,
)
from polyzymd.analyses.distances._results import (
    DistancePairAggregatedResult,
    DistancePairResult,
    DistanceResultMetadata,
)
from polyzymd.analyses.mda import (
    ArtifactStore,
    ConditionArtifact,
    MDAAnalysisJob,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    PairDistanceSpec,
    ReplicateArtifact,
    build_pair_distance_analysis,
    pair_distance_version,
)
from polyzymd.analyses.mda.plugin import frame_selection_payload, strict_json_payload
from polyzymd.analyses.shared.config_hash import compute_config_hash
from polyzymd.analyses.shared.loader import parse_time_string

if TYPE_CHECKING:
    from polyzymd.analyses.catalytic_triad import CatalyticTriadSettings
    from polyzymd.analyses.mda import ArtifactSidecarRef, MDAReplicateJobContext

LOGGER = logging.getLogger(__name__)
SIMULTANEOUS_CONTACT_METRIC = "simultaneous_contact_fraction"
SIMULTANEOUS_CONTACT_METADATA = {
    "label": "Simultaneous Contact",
    "unit": "%",
    "higher_is_better": True,
    "direction_labels": ("worsening", "unchanged", "improving"),
}


def build_triad_jobs(
    ctx: MDAReplicateJobContext,
    settings: CatalyticTriadSettings,
) -> list[MDAAnalysisJob]:
    """Build the pair-distance job that feeds catalytic-triad reduction.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        Framework-provided MDAnalysis replicate context.
    settings : CatalyticTriadSettings
        Resolved catalytic-triad settings.

    Returns
    -------
    list of MDAAnalysisJob
        Single pair-distance job whose pair x frame matrix is reduced by the
        triad collector.
    """

    from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory

    alignment = AlignmentConfig()
    if alignment.enabled:
        align_trajectory(
            ctx.universe,
            alignment,
            start_frame=ctx.frame_selection.start,
            stop_frame=ctx.frame_selection.stop,
            step_frame=ctx.frame_selection.step,
        )
    resolved_pairs = resolve_distance_pairs(
        universe=ctx.universe,
        pairs=settings.get_pair_selections(),
        thresholds=[float(settings.threshold)] * settings.n_pairs,
        pair_labels=settings.get_pair_labels(),
    )
    metadata = {
        **ctx.universe_policy.metadata,
        "triad_settings": settings.model_dump(mode="json"),
        "pair_distance_version": pair_distance_version(),
        "pbc_policy": {"use_pbc": True, "box_source": "timestep.dimensions"},
        "composite_reducer": {
            "metric": SIMULTANEOUS_CONTACT_METRIC,
            "threshold_operator": "strict_less_than",
        },
    }
    return [
        MDAAnalysisJob(
            name="triad_pair_distances",
            analysis=build_pair_distance_analysis(
                universe=ctx.universe,
                pairs=resolved_pairs,
                use_pbc=True,
            ),
            frame_selection=ctx.frame_selection,
            universe_policy=MDAUniversePolicy(
                condition_label=ctx.replicate_context.condition.label,
                replicate=ctx.replicate,
                provenance=ctx.universe_policy.provenance,
                metadata=metadata,
            ),
        )
    ]


class TriadArtifactCollector:
    """Collect pair-distance jobs into a catalytic-triad replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Map completed pair-distance results to triad JSON and NPZ sidecars.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context.
        completed_jobs : sequence of MDAJobResult
            Completed pair-distance job results.

        Returns
        -------
        ReplicateArtifact
            Canonical replicate artifact with summary JSON and one matrix
            sidecar.
        """

        if len(completed_jobs) != 1:
            raise ValueError(
                f"Catalytic triad expected one pair-distance job, got {len(completed_jobs)}"
            )
        settings = ctx.settings
        job = completed_jobs[0]
        resolved_pairs = _pair_specs_from_settings(settings)
        payload = payload_from_pair_distance_analysis(
            analysis=job.analysis,
            resolved_pairs=resolved_pairs,
            timestep_ps=float(ctx.frame_selection.timestep_ps or 1.0),
            step=int(ctx.frame_selection.step or 1),
            n_frames_total=_n_frames_total(ctx, job),
        )
        distance_matrix = _distance_matrix(payload)
        simultaneous = _simultaneous_contact_vector(distance_matrix, float(settings.threshold))
        simultaneous_fraction = float(np.mean(simultaneous)) if simultaneous.size else 0.0
        n_frames_simultaneous = int(np.count_nonzero(simultaneous))
        sidecar = _write_replicate_sidecar(ctx, job, payload, distance_matrix, simultaneous)
        pair_payloads = [
            _pair_payload_to_json(pair_payload) for pair_payload in payload.pair_payloads
        ]
        sim_stats = _simultaneous_contact_statistics(
            simultaneous=simultaneous,
            fraction=simultaneous_fraction,
            timestep_ps=_effective_timestep_ps(payload, ctx),
        )
        eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
        config_hash = compute_config_hash(ctx.replicate_context.sim_config)
        warnings = _unique_warnings([*ctx.warnings, *payload.warnings])
        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "triad_name": settings.name,
                "triad_description": settings.description,
                "pairs": pair_payloads,
                "pair_results": pair_payloads,
                "threshold": float(settings.threshold),
                "simultaneous_contact_fraction": simultaneous_fraction,
                "n_frames_simultaneous": n_frames_simultaneous,
                "sim_contact_sem": sim_stats["sem"],
                "sim_contact_correlation_time": sim_stats["correlation_time"],
                "sim_contact_correlation_time_unit": sim_stats["correlation_time_unit"],
                "sim_contact_n_independent": sim_stats["n_independent"],
                "sim_contact_warning": sim_stats["warning"],
                "metrics": {SIMULTANEOUS_CONTACT_METRIC: simultaneous_fraction * 100.0},
                "replicate_metrics": {SIMULTANEOUS_CONTACT_METRIC: simultaneous_fraction * 100.0},
                "n_pairs": len(pair_payloads),
                "n_frames_total": payload.n_frames_total,
                "n_frames_used": payload.n_frames_used,
            },
            sidecars=[sidecar],
            provenance={
                "source": "mda_pair_distance_composite_reducer",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
                "pbc_policy": {"use_pbc": True, "box_source": "timestep.dimensions"},
                "composite_reducer": {
                    "metric": SIMULTANEOUS_CONTACT_METRIC,
                    "threshold": float(settings.threshold),
                    "threshold_operator": "strict_less_than",
                },
            },
            metadata={
                "result_kind": "catalytic_triad_mda_replicate",
                "settings_fingerprint": ctx.settings_fingerprint,
                "config_hash": config_hash,
                "polyzymd_version": get_polyzymd_version(),
                "mdanalysis_version": mdanalysis_version(),
                "pair_distance_version": pair_distance_version(),
                "equilibration_time": eq_value,
                "equilibration_unit": eq_unit,
                "selection_string": _combined_selection(settings.get_pair_selections()),
            },
            warnings=warnings,
        )


def aggregate_triad_artifacts(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: CatalyticTriadSettings,
    equilibration: str,
    output_dir: Path,
    result_path: Path,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> ConditionArtifact:
    """Aggregate catalytic-triad replicate artifacts into a condition artifact.

    Parameters
    ----------
    condition_label : str
        Condition label being aggregated.
    replicates : sequence of int
        Expected replicate IDs.
    settings : CatalyticTriadSettings
        Active triad settings.
    equilibration : str
        Equilibration string from the framework context.
    output_dir : Path
        Aggregated output directory.
    result_path : Path
        Canonical condition artifact path.
    artifacts : sequence of ReplicateArtifact
        Per-replicate triad artifacts.
    settings_fingerprint : str
        Active settings fingerprint.

    Returns
    -------
    ConditionArtifact
        Canonical condition artifact with legacy-compatible summaries.
    """

    ordered_artifacts = _validate_and_order_artifacts(
        condition_label=condition_label,
        expected_replicates=replicates,
        settings=settings,
        settings_fingerprint=settings_fingerprint,
        artifacts=artifacts,
        analysis_dir=output_dir.parent,
    )
    legacy_result = _aggregate_legacy_result(
        condition_label=condition_label,
        replicates=replicates,
        settings=settings,
        equilibration=equilibration,
        artifacts=ordered_artifacts,
        settings_fingerprint=settings_fingerprint,
    )
    metrics, replicate_metrics = _condition_metrics(legacy_result)
    artifact = ConditionArtifact(
        analysis_name="catalytic_triad",
        condition_label=condition_label,
        replicates=[int(rep) for rep in replicates],
        payload={
            "triad_name": legacy_result.triad_name,
            "triad_description": legacy_result.triad_description,
            "pairs": [pair.model_dump(mode="json") for pair in legacy_result.pair_results],
            "pair_results": [pair.model_dump(mode="json") for pair in legacy_result.pair_results],
            "threshold": legacy_result.threshold,
            "overall_simultaneous_contact": legacy_result.overall_simultaneous_contact,
            "sem_simultaneous_contact": legacy_result.sem_simultaneous_contact,
            "per_replicate_simultaneous": legacy_result.per_replicate_simultaneous,
            "metrics": metrics,
            "replicate_metrics": replicate_metrics,
            "metric_metadata": {SIMULTANEOUS_CONTACT_METRIC: SIMULTANEOUS_CONTACT_METADATA},
            "n_replicates": legacy_result.n_replicates,
        },
        provenance={
            "source": "triad_replicate_artifact_aggregation",
            "frame_selection": ordered_artifacts[0].provenance.get("frame_selection"),
            "composite_reducer": {
                "metric": SIMULTANEOUS_CONTACT_METRIC,
                "threshold": float(settings.threshold),
                "threshold_operator": "strict_less_than",
            },
        },
        metadata={
            "result_kind": "catalytic_triad_mda_condition",
            "settings_fingerprint": settings_fingerprint,
            "config_hash": legacy_result.config_hash,
            "polyzymd_version": legacy_result.polyzymd_version,
            "mdanalysis_version": mdanalysis_version(),
            "pair_distance_version": pair_distance_version(),
            "equilibration_time": legacy_result.equilibration_time,
            "equilibration_unit": legacy_result.equilibration_unit,
            "selection_string": legacy_result.selection_string,
            "source_result_files": [
                str(output_dir.parent / f"run_{replicate}" / "result.json")
                for replicate in replicates
            ],
            "n_replicates": legacy_result.n_replicates,
        },
        source_replicates=[
            {
                "replicate": int(replicate),
                "path": str(output_dir.parent / f"run_{replicate}" / "result.json"),
            }
            for replicate in replicates
        ],
        warnings=_combined_warnings(ordered_artifacts),
    )
    ArtifactStore(result_path.parent).write_condition_result(artifact, result_path.name)
    return artifact


def artifact_to_triad_result(artifact: ReplicateArtifact) -> TriadResult:
    """Convert a replicate artifact to the established triad result model."""

    if artifact.analysis_name != "catalytic_triad":
        raise ValueError(f"Expected catalytic_triad artifact, got {artifact.analysis_name!r}")
    metadata = artifact.metadata
    payload = artifact.payload
    return TriadResult(
        config_hash=str(metadata.get("config_hash", "unknown")),
        polyzymd_version=str(metadata.get("polyzymd_version", get_polyzymd_version())),
        replicate=artifact.replicate,
        equilibration_time=float(metadata.get("equilibration_time", 0.0)),
        equilibration_unit=str(metadata.get("equilibration_unit", "ns")),
        selection_string=str(metadata.get("selection_string", "")),
        triad_name=str(payload.get("triad_name", "catalytic_triad")),
        triad_description=payload.get("triad_description"),
        pair_results=[DistancePairResult.model_validate(pair) for pair in payload.get("pairs", [])],
        threshold=float(payload.get("threshold", 3.5)),
        simultaneous_contact_fraction=float(payload.get("simultaneous_contact_fraction", 0.0)),
        n_frames_simultaneous=int(payload.get("n_frames_simultaneous", 0) or 0),
        simultaneous_contact_timeseries=None,
        sim_contact_sem=payload.get("sim_contact_sem"),
        sim_contact_correlation_time=payload.get("sim_contact_correlation_time"),
        sim_contact_correlation_time_unit=payload.get("sim_contact_correlation_time_unit"),
        sim_contact_n_independent=payload.get("sim_contact_n_independent"),
        sim_contact_warning=payload.get("sim_contact_warning"),
        n_frames_total=int(payload.get("n_frames_total", 0) or 0),
        n_frames_used=int(payload.get("n_frames_used", 0) or 0),
        settings_fingerprint=metadata.get("settings_fingerprint"),
    )


def condition_artifact_to_legacy_result(artifact: Any) -> TriadAggregatedResult:
    """Convert a condition artifact to the established aggregate model."""

    if isinstance(artifact, TriadAggregatedResult):
        return artifact
    if not isinstance(artifact, ConditionArtifact):
        raise TypeError(
            "Catalytic-triad aggregate adapters require canonical ConditionArtifact inputs. "
            f"Got {type(artifact).__name__}; recompute the condition or clear stale caches."
        )
    metadata = artifact.metadata
    payload = artifact.payload
    return TriadAggregatedResult(
        config_hash=str(metadata.get("config_hash", "unknown")),
        polyzymd_version=str(metadata.get("polyzymd_version", get_polyzymd_version())),
        replicate=None,
        equilibration_time=float(metadata.get("equilibration_time", 0.0)),
        equilibration_unit=str(metadata.get("equilibration_unit", "ns")),
        selection_string=str(metadata.get("selection_string", "")),
        replicates=[int(rep) for rep in artifact.replicates],
        n_replicates=len(artifact.replicates),
        triad_name=str(payload.get("triad_name", "catalytic_triad")),
        triad_description=payload.get("triad_description"),
        pair_results=[
            DistancePairAggregatedResult.model_validate(pair) for pair in payload.get("pairs", [])
        ],
        threshold=float(payload.get("threshold", 3.5)),
        overall_simultaneous_contact=float(payload.get("overall_simultaneous_contact", 0.0)),
        sem_simultaneous_contact=float(payload.get("sem_simultaneous_contact", 0.0)),
        per_replicate_simultaneous=[
            float(value) for value in payload.get("per_replicate_simultaneous", [])
        ],
        source_result_files=[str(path) for path in metadata.get("source_result_files", [])],
        settings_fingerprint=metadata.get("settings_fingerprint"),
    )


def _aggregate_legacy_result(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: CatalyticTriadSettings,
    equilibration: str,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> TriadAggregatedResult:
    """Build the established aggregate model from replicate artifacts."""

    from polyzymd.analyses.shared.aggregation import aggregate_distance_pair_stats
    from polyzymd.analyses.shared.statistics import compute_sem

    del condition_label
    eq_value, eq_unit = parse_time_string(equilibration)
    legacy_results = [artifact_to_triad_result(artifact) for artifact in artifacts]
    first = legacy_results[0]
    metadata = DistanceResultMetadata(
        config_hash=first.config_hash,
        polyzymd_version=get_polyzymd_version(),
        replicate=None,
        equilibration_time=eq_value,
        equilibration_unit=eq_unit,
    )
    aggregated_pairs: list[DistancePairAggregatedResult] = []
    for pair_idx, pair in enumerate(first.pair_results):
        stats = aggregate_distance_pair_stats(legacy_results, pair_idx)
        aggregated_pairs.append(
            DistancePairAggregatedResult.from_aggregated_stats(
                metadata,
                pair,
                stats,
                replicates=replicates,
                threshold=settings.threshold,
                pair_label=pair.pair_label,
            )
        )
    per_replicate_simultaneous = [result.simultaneous_contact_fraction for result in legacy_results]
    sim_stats = compute_sem(per_replicate_simultaneous)
    return TriadAggregatedResult(
        config_hash=first.config_hash,
        polyzymd_version=metadata.polyzymd_version,
        replicate=None,
        equilibration_time=eq_value,
        equilibration_unit=eq_unit,
        selection_string=_combined_selection(settings.get_pair_selections()),
        replicates=[int(rep) for rep in replicates],
        n_replicates=len(replicates),
        triad_name=settings.name,
        triad_description=settings.description,
        pair_results=aggregated_pairs,
        threshold=settings.threshold,
        overall_simultaneous_contact=sim_stats.mean,
        sem_simultaneous_contact=sim_stats.sem,
        per_replicate_simultaneous=per_replicate_simultaneous,
        source_result_files=[],
        settings_fingerprint=settings_fingerprint,
    )


def _validate_and_order_artifacts(
    *,
    condition_label: str,
    expected_replicates: Sequence[int],
    settings: CatalyticTriadSettings,
    settings_fingerprint: str,
    artifacts: Sequence[ReplicateArtifact],
    analysis_dir: Path,
) -> list[ReplicateArtifact]:
    """Validate triad artifact identity, settings, pairs, and sidecars."""

    expected = [int(rep) for rep in expected_replicates]
    by_replicate: dict[int, ReplicateArtifact] = {}
    for artifact in artifacts:
        if artifact.analysis_name != "catalytic_triad":
            raise ValueError(f"Expected catalytic_triad artifact, got {artifact.analysis_name!r}")
        if artifact.condition_label != condition_label:
            raise ValueError(
                f"Catalytic-triad artifact condition mismatch: expected {condition_label!r}, "
                f"got {artifact.condition_label!r}"
            )
        if artifact.replicate in by_replicate:
            raise ValueError(
                f"Duplicate catalytic-triad artifact for replicate {artifact.replicate}"
            )
        if artifact.metadata.get("settings_fingerprint") != settings_fingerprint:
            raise ValueError(
                f"Catalytic-triad artifact replicate {artifact.replicate} has settings "
                f"fingerprint {artifact.metadata.get('settings_fingerprint')}, expected "
                f"{settings_fingerprint}. Recompute the condition or clear stale caches."
            )
        _validate_pair_payloads(artifact, settings)
        store = ArtifactStore(analysis_dir / f"run_{artifact.replicate}")
        store.validate_sidecar(_triad_distance_sidecar(artifact))
        by_replicate[int(artifact.replicate)] = artifact
    observed = sorted(by_replicate)
    if observed != sorted(expected):
        raise ValueError(
            f"Catalytic-triad artifacts for condition '{condition_label}' do not match expected "
            f"replicates {expected}: found {observed}. Recompute missing replicates or clear "
            "stale caches before aggregating."
        )
    return [by_replicate[rep] for rep in expected]


def _validate_pair_payloads(artifact: ReplicateArtifact, settings: CatalyticTriadSettings) -> None:
    """Validate pair payloads against active triad settings."""

    pairs = artifact.payload.get("pairs", [])
    if len(pairs) != len(settings.pairs):
        raise ValueError(
            f"Catalytic-triad replicate {artifact.replicate} has {len(pairs)} pairs, "
            f"expected {len(settings.pairs)}"
        )
    for idx, (payload, pair_setting) in enumerate(zip(pairs, settings.pairs, strict=True), start=1):
        if payload.get("pair_label") != pair_setting.label:
            raise ValueError(
                f"Catalytic-triad replicate {artifact.replicate} pair {idx} label mismatch"
            )
        if (
            payload.get("selection1") != pair_setting.selection_a
            or payload.get("selection2") != pair_setting.selection_b
        ):
            raise ValueError(
                f"Catalytic-triad replicate {artifact.replicate} pair {idx} selection mismatch"
            )
        if payload.get("threshold") != settings.threshold:
            raise ValueError(
                f"Catalytic-triad replicate {artifact.replicate} pair {idx} threshold mismatch"
            )
    if artifact.payload.get("threshold") != settings.threshold:
        raise ValueError(
            f"Catalytic-triad replicate {artifact.replicate} threshold mismatch: "
            f"{artifact.payload.get('threshold')!r} != {settings.threshold!r}"
        )


def _condition_metrics(
    legacy_result: TriadAggregatedResult,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, float]]]:
    """Build percent-scaled condition metrics from a triad aggregate."""

    values = [float(value) * 100.0 for value in legacy_result.per_replicate_simultaneous]
    metric = _metric_summary(
        SIMULTANEOUS_CONTACT_METRIC,
        values,
        legacy_result.overall_simultaneous_contact * 100.0,
        legacy_result.sem_simultaneous_contact * 100.0,
    )
    metric.update(SIMULTANEOUS_CONTACT_METADATA)
    replicate_metrics = {
        str(replicate): {SIMULTANEOUS_CONTACT_METRIC: value}
        for replicate, value in zip(legacy_result.replicates, values, strict=True)
    }
    return {SIMULTANEOUS_CONTACT_METRIC: metric}, replicate_metrics


def _metric_summary(name: str, values: Sequence[float], mean: float, sem: float) -> dict[str, Any]:
    """Return an artifact metric summary."""

    array = np.asarray(values, dtype=np.float64)
    return {
        "name": name,
        "values": [float(value) for value in values],
        "mean": float(mean),
        "sem": float(sem),
        "std": float(np.std(array, ddof=1)) if array.size > 1 else 0.0,
        "n": int(array.size),
    }


def _write_replicate_sidecar(
    ctx: MDACollectorContext,
    job: MDAJobResult,
    payload: DistancesRunnerPayload,
    distance_matrix: NDArray[np.float64],
    simultaneous: NDArray[np.bool_],
) -> ArtifactSidecarRef:
    """Write the triad pair x frame matrix and simultaneous vector sidecar."""

    thresholds = np.asarray(
        [np.nan if pair.threshold is None else pair.threshold for pair in payload.pair_payloads],
        dtype=np.float64,
    )
    return ctx.artifact_store.write_npz_sidecar(
        "sidecars/triad_pair_distances.npz",
        distance_matrix=distance_matrix,
        frames=np.asarray(payload.frames if payload.frames is not None else [], dtype=np.int64),
        time_ns=np.asarray(
            payload.time_ns if payload.time_ns is not None else [], dtype=np.float64
        ),
        thresholds=thresholds,
        pair_labels=np.asarray([pair.pair_label for pair in payload.pair_payloads]),
        selections_a=np.asarray([pair.selection1 for pair in payload.pair_payloads]),
        selections_b=np.asarray([pair.selection2 for pair in payload.pair_payloads]),
        simultaneous_contact=simultaneous.astype(np.bool_),
        metadata={
            "job_name": job.name,
            "kind": "triad_distance_matrix",
            "layout": "pair_x_frame",
            "threshold_operator": "strict_less_than",
        },
    )


def _pair_specs_from_settings(settings: CatalyticTriadSettings) -> list[PairDistanceSpec]:
    """Build metadata-only pair specs from triad settings."""

    return [
        PairDistanceSpec(
            label=pair.label,
            selection_a=pair.selection_a,
            selection_b=pair.selection_b,
            atoms_a=None,
            atoms_b=None,
            mode_a=None,
            mode_b=None,
            threshold=float(settings.threshold),
        )
        for pair in settings.pairs
    ]


def _distance_matrix(payload: DistancesRunnerPayload) -> NDArray[np.float64]:
    """Return a finite pair x frame distance matrix from distance payloads."""

    if not payload.pair_payloads:
        raise ValueError("Catalytic-triad pair-distance job produced no pair payloads")
    matrix = np.vstack([pair.distances for pair in payload.pair_payloads]).astype(np.float64)
    if matrix.ndim != 2 or matrix.shape[1] != payload.n_frames_used:
        raise ValueError(
            "Catalytic-triad distance matrix shape does not match the analyzed frame count"
        )
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Catalytic-triad distance matrix contains non-finite values")
    return matrix


def _simultaneous_contact_vector(
    distance_matrix: NDArray[np.float64], threshold: float
) -> NDArray[np.bool_]:
    """Return frames where every pair is strictly below the contact threshold."""

    return np.all(distance_matrix < float(threshold), axis=0).astype(np.bool_)


def _simultaneous_contact_statistics(
    *,
    simultaneous: NDArray[np.bool_],
    fraction: float,
    timestep_ps: float,
) -> dict[str, float | int | str | None]:
    """Return autocorrelation-corrected statistics for simultaneous contact."""

    if simultaneous.size < 20:
        return {
            "sem": None,
            "correlation_time": None,
            "correlation_time_unit": None,
            "n_independent": None,
            "warning": None,
        }
    from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

    try:
        tau_result = estimate_correlation_time(
            simultaneous.astype(np.float64),
            timestep=timestep_ps,
            timestep_unit="ps",
            method="integration",
            n_frames=int(simultaneous.size),
        )
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
        LOGGER.warning("Autocorrelation analysis for contact timeseries failed: %s", exc)
        sem = float(np.sqrt(fraction * (1.0 - fraction) / float(simultaneous.size)))
        return {
            "sem": sem,
            "correlation_time": None,
            "correlation_time_unit": None,
            "n_independent": None,
            "warning": str(exc),
        }
    n_independent = int(tau_result.n_independent)
    sem = (
        float(np.sqrt(fraction * (1.0 - fraction) / float(n_independent)))
        if n_independent > 0
        else None
    )
    return {
        "sem": sem,
        "correlation_time": tau_result.tau,
        "correlation_time_unit": tau_result.tau_unit,
        "n_independent": n_independent,
        "warning": tau_result.warning,
    }


def _effective_timestep_ps(payload: DistancesRunnerPayload, ctx: MDACollectorContext) -> float:
    """Return the effective analyzed-frame timestep in picoseconds."""

    time_ns = payload.time_ns
    if time_ns is not None and time_ns.size >= 2:
        diffs_ps = np.diff(time_ns.astype(np.float64)) * 1000.0
        finite = diffs_ps[np.isfinite(diffs_ps) & (diffs_ps > 0)]
        if finite.size:
            return float(np.median(finite))
    return float(ctx.frame_selection.timestep_ps or 1.0) * float(ctx.frame_selection.step or 1)


def _n_frames_total(ctx: MDACollectorContext, job: MDAJobResult) -> int:
    """Return total trajectory frame count with a result-matrix fallback."""

    if ctx.frame_selection.n_frames_total is not None:
        return int(ctx.frame_selection.n_frames_total)
    matrix = np.asarray(getattr(job.results, "distance_matrix", []))
    if matrix.ndim == 2:
        return int(matrix.shape[1])
    return 0


def _combined_selection(pairs: Sequence[tuple[str, str]]) -> str:
    """Return the established combined selection string."""

    return "; ".join(f"({selection1} : {selection2})" for selection1, selection2 in pairs)


def _unique_warnings(warnings: Sequence[str]) -> list[str]:
    """Return unique warnings while preserving first-seen order."""

    unique: list[str] = []
    seen: set[str] = set()
    for warning in warnings:
        text = str(warning)
        if text not in seen:
            seen.add(text)
            unique.append(text)
    return unique


def _combined_warnings(artifacts: Sequence[ReplicateArtifact]) -> list[str]:
    """Combine warnings from replicate artifacts."""

    return _unique_warnings([warning for artifact in artifacts for warning in artifact.warnings])


def _triad_distance_sidecar(artifact: ReplicateArtifact) -> ArtifactSidecarRef:
    """Return the triad distance-matrix sidecar from a replicate artifact."""

    matches = [
        sidecar
        for sidecar in artifact.sidecars
        if sidecar.metadata.get("kind") == "triad_distance_matrix"
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Catalytic-triad replicate {artifact.replicate} requires exactly one triad "
            f"distance sidecar, found {len(matches)}"
        )
    return matches[0]


def load_replicate_distance_sidecar(
    analysis_dir: Path,
    replicate: int,
) -> tuple[ReplicateArtifact, Path]:
    """Load and validate a replicate triad artifact and distance sidecar path."""

    run_dir = analysis_dir / f"run_{replicate}"
    store = ArtifactStore(run_dir)
    artifact = store.read_replicate_result("result.json")
    sidecar = _triad_distance_sidecar(artifact)
    return artifact, store.validate_sidecar(sidecar)


def load_condition_artifact(aggregated_dir: Path) -> ConditionArtifact | None:
    """Load the canonical condition artifact if it exists."""

    result_path = aggregated_dir / "result.json"
    if not result_path.exists():
        return None
    return ArtifactStore(aggregated_dir).read_condition_result("result.json")
