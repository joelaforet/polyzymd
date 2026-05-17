"""MDAnalysis-native distances jobs and artifact adapters."""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses._results_base import get_polyzymd_version
from polyzymd.analyses.distances._results import (
    DistanceAggregatedResult,
    DistancePairAggregatedResult,
    DistancePairResult,
    DistanceResult,
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
    from polyzymd.analyses.distances import DistancesSettings
    from polyzymd.analyses.mda import ArtifactSidecarRef, MDAReplicateJobContext

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class DistancePairPayload:
    """Trajectory-native payload for one distance pair."""

    pair_label: str
    selection1: str
    selection2: str
    distances: NDArray[np.float64]
    mean_distance: float
    std_distance: float
    median_distance: float
    min_distance: float
    max_distance: float
    sem_distance: float | None
    correlation_time: float | None
    correlation_time_unit: str | None
    n_independent_frames: int | None
    statistical_inefficiency: float | None
    autocorrelation_warning: str | None
    threshold: float | None
    fraction_below_threshold: float | None
    histogram_edges: NDArray[np.float64]
    histogram_counts: NDArray[np.int64]
    kde_x: NDArray[np.float64] | None
    kde_y: NDArray[np.float64] | None
    kde_peak: float | None
    kde_bandwidth: float | None
    n_frames_total: int
    n_frames_used: int


@dataclass(frozen=True)
class DistanceReplicatePayload:
    """Collection of distance-pair payloads for one replicate."""

    n_frames_total: int
    n_frames_used: int
    pair_payloads: list[DistancePairPayload]
    frames: NDArray[np.int64] | None = None
    time_ns: NDArray[np.float64] | None = None
    warnings: tuple[str, ...] = ()


def build_distance_jobs(
    ctx: MDAReplicateJobContext, settings: DistancesSettings
) -> list[MDAAnalysisJob]:
    """Build the custom pair-distance MDAnalysis job for one replicate.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        Framework-provided MDAnalysis replicate context.
    settings : DistancesSettings
        Resolved distances settings.

    Returns
    -------
    list of MDAAnalysisJob
        Single job measuring all configured pairs as a pair x frame matrix.
    """

    from polyzymd.analyses.shared.alignment import align_trajectory

    alignment = settings.get_alignment_config()
    if getattr(alignment, "enabled", False):
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
        thresholds=settings.get_pair_thresholds(),
        pair_labels=settings.get_pair_labels(),
    )
    metadata = {
        **ctx.universe_policy.metadata,
        "distance_settings": settings.model_dump(mode="json"),
        "pair_distance_version": pair_distance_version(),
        "pbc_policy": {"use_pbc": settings.use_pbc, "box_source": "timestep.dimensions"},
    }
    return [
        MDAAnalysisJob(
            name="pair_distances",
            analysis=build_pair_distance_analysis(
                universe=ctx.universe,
                pairs=resolved_pairs,
                use_pbc=settings.use_pbc,
            ),
            frame_selection=ctx.frame_selection,
            backend_policy=ctx.backend_policy,
            universe_policy=MDAUniversePolicy(
                condition_label=ctx.replicate_context.condition.label,
                replicate=ctx.replicate,
                provenance=ctx.universe_policy.provenance,
                metadata=metadata,
            ),
        )
    ]


def compute_distance_payloads(
    *,
    universe: Any,
    pairs: Sequence[tuple[str, str]],
    thresholds: Sequence[float | None],
    start: int,
    stop: int,
    step: int,
    timestep_ps: float,
    use_pbc: bool,
    alignment: Any,
    pair_label_func: Callable[[str, str], str],
) -> DistanceReplicatePayload:
    """Compute distance payloads through the shared pair-distance primitive.

    This compatibility helper keeps direct callers and the catalytic-triad
    runner on the same MDAnalysis ``AnalysisBase`` implementation used by the
    migrated distances plugin.
    """

    from polyzymd.analyses.shared.alignment import align_trajectory

    if getattr(alignment, "enabled", False):
        align_trajectory(universe, alignment, start_frame=start, stop_frame=stop, step_frame=step)
    resolved_pairs = resolve_distance_pairs(
        universe=universe,
        pairs=pairs,
        thresholds=thresholds,
        pair_label_func=pair_label_func,
    )
    analysis = build_pair_distance_analysis(
        universe=universe,
        pairs=resolved_pairs,
        use_pbc=use_pbc,
    ).run(start=start, stop=stop, step=step)
    return payload_from_pair_distance_analysis(
        analysis=analysis,
        resolved_pairs=resolved_pairs,
        timestep_ps=timestep_ps,
        step=step,
        n_frames_total=len(universe.trajectory),
    )


def resolve_distance_pairs(
    *,
    universe: Any,
    pairs: Sequence[tuple[str, str]],
    thresholds: Sequence[float | None],
    pair_label_func: Callable[[str, str], str] | None = None,
    pair_labels: Sequence[str] | None = None,
) -> list[PairDistanceSpec]:
    """Resolve and validate configured distance-pair selections."""

    from polyzymd.analyses.shared.diagnostics import (
        get_selection_diagnostics,
        warn_if_multi_chain_selection,
    )
    from polyzymd.analyses.shared.selections import parse_selection_string

    if pair_labels is not None and len(pair_labels) != len(pairs):
        raise ValueError(
            f"Distances pair label count {len(pair_labels)} does not match pair count {len(pairs)}"
        )
    if pair_labels is None and pair_label_func is None:
        raise ValueError("Either pair_labels or pair_label_func must be provided")

    resolved_pairs: list[PairDistanceSpec] = []
    for pair_index, ((selection1, selection2), threshold) in enumerate(
        zip(pairs, thresholds, strict=True)
    ):
        parsed1 = parse_selection_string(selection1)
        parsed2 = parse_selection_string(selection2)
        atoms1 = universe.select_atoms(parsed1.selection)
        atoms2 = universe.select_atoms(parsed2.selection)
        if len(atoms1) == 0:
            diagnostics = get_selection_diagnostics(universe, selection1)
            raise ValueError(f"Selection '{selection1}' matched no atoms.\n\n{diagnostics}")
        if len(atoms2) == 0:
            diagnostics = get_selection_diagnostics(universe, selection2)
            raise ValueError(f"Selection '{selection2}' matched no atoms.\n\n{diagnostics}")
        pair_label = (
            str(pair_labels[pair_index])
            if pair_labels is not None
            else str(pair_label_func(selection1, selection2))
        )
        warn_if_multi_chain_selection(atoms1, selection1, f"for distance pair '{pair_label}'")
        warn_if_multi_chain_selection(atoms2, selection2, f"for distance pair '{pair_label}'")
        resolved_pairs.append(
            PairDistanceSpec(
                label=pair_label,
                selection_a=selection1,
                selection_b=selection2,
                atoms_a=atoms1,
                atoms_b=atoms2,
                mode_a=parsed1.mode,
                mode_b=parsed2.mode,
                threshold=threshold,
            )
        )
    return resolved_pairs


def payload_from_pair_distance_analysis(
    *,
    analysis: Any,
    resolved_pairs: Sequence[PairDistanceSpec],
    timestep_ps: float,
    step: int,
    n_frames_total: int,
) -> DistanceReplicatePayload:
    """Build summary payloads from a completed pair-distance analysis."""

    distance_matrix = np.asarray(analysis.results.distance_matrix, dtype=np.float64)
    frames = np.asarray(getattr(analysis.results, "frames", []), dtype=np.int64)
    time_ns = _analysis_time_ns(analysis, frames, timestep_ps)
    effective_timestep_ps = _effective_timestep_ps(time_ns, timestep_ps, step)
    pair_payloads = [
        _summarize_distance_series(
            distances=np.asarray(distance_matrix[pair_index], dtype=np.float64),
            pair_label=pair.label,
            selection1=pair.selection_a,
            selection2=pair.selection_b,
            threshold=pair.threshold,
            timestep_ps=effective_timestep_ps,
            n_frames_total=n_frames_total,
        )
        for pair_index, pair in enumerate(resolved_pairs)
    ]
    return DistanceReplicatePayload(
        n_frames_total=n_frames_total,
        n_frames_used=int(distance_matrix.shape[1]) if distance_matrix.ndim == 2 else 0,
        pair_payloads=pair_payloads,
        frames=frames,
        time_ns=time_ns,
        warnings=tuple(str(w) for w in getattr(analysis.results, "warnings", [])),
    )


class DistanceArtifactCollector:
    """Collect completed pair-distance jobs into a replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Map completed distance jobs to summary JSON and NPZ sidecars."""

        if len(completed_jobs) != 1:
            raise ValueError(f"Distances expected one pair-distance job, got {len(completed_jobs)}")
        job = completed_jobs[0]
        settings = ctx.settings
        resolved_pairs = _pair_specs_from_settings(settings)
        payload = payload_from_pair_distance_analysis(
            analysis=job.analysis,
            resolved_pairs=resolved_pairs,
            timestep_ps=float(ctx.frame_selection.timestep_ps or 1.0),
            step=int(ctx.frame_selection.step or 1),
            n_frames_total=_n_frames_total(ctx, job),
        )
        sidecar = _write_replicate_sidecar(ctx, job, payload)
        pair_payloads = [
            _pair_payload_to_json(pair_payload) for pair_payload in payload.pair_payloads
        ]
        metrics = _replicate_metrics(pair_payloads)
        eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
        config_hash = compute_config_hash(ctx.replicate_context.sim_config)
        warnings = _unique_warnings([*ctx.warnings, *payload.warnings])
        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "pairs": pair_payloads,
                "pair_results": pair_payloads,
                "metrics": metrics,
                "replicate_metrics": metrics,
                "n_pairs": len(pair_payloads),
                "n_frames_total": payload.n_frames_total,
                "n_frames_used": payload.n_frames_used,
            },
            sidecars=[sidecar],
            provenance={
                "source": "mda_pair_distance_jobs",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
                "pbc_policy": {"use_pbc": settings.use_pbc, "box_source": "timestep.dimensions"},
            },
            metadata={
                "result_kind": "distances_mda_replicate",
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


def aggregate_distance_artifacts(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: DistancesSettings,
    equilibration: str,
    output_dir: Path,
    result_path: Path,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> ConditionArtifact:
    """Aggregate distance replicate artifacts into a condition artifact."""

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
    sidecar = _write_condition_sidecar(output_dir, ordered_artifacts, legacy_result)
    metrics, replicate_metrics = _condition_metrics(legacy_result)
    artifact = ConditionArtifact(
        analysis_name="distances",
        condition_label=condition_label,
        replicates=[int(rep) for rep in replicates],
        payload={
            "pairs": [pair.model_dump(mode="json") for pair in legacy_result.pair_results],
            "pair_results": [pair.model_dump(mode="json") for pair in legacy_result.pair_results],
            "metrics": metrics,
            "replicate_metrics": replicate_metrics,
            "n_replicates": legacy_result.n_replicates,
        },
        sidecars=[sidecar],
        provenance={
            "source": "distance_replicate_artifact_aggregation",
            "frame_selection": ordered_artifacts[0].provenance.get("frame_selection"),
        },
        metadata={
            "result_kind": "distances_mda_condition",
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


def artifact_to_distance_result(artifact: ReplicateArtifact) -> DistanceResult:
    """Convert a replicate artifact to the established distance result model."""

    if artifact.analysis_name != "distances":
        raise ValueError(f"Expected distances artifact, got {artifact.analysis_name!r}")
    metadata = artifact.metadata
    payload = artifact.payload
    return DistanceResult(
        config_hash=str(metadata.get("config_hash", "unknown")),
        polyzymd_version=str(metadata.get("polyzymd_version", get_polyzymd_version())),
        replicate=artifact.replicate,
        equilibration_time=float(metadata.get("equilibration_time", 0.0)),
        equilibration_unit=str(metadata.get("equilibration_unit", "ns")),
        selection_string=str(metadata.get("selection_string", "")),
        pair_results=[DistancePairResult.model_validate(pair) for pair in payload.get("pairs", [])],
        n_frames_total=int(payload.get("n_frames_total", 0) or 0),
        n_frames_used=int(payload.get("n_frames_used", 0) or 0),
        settings_fingerprint=metadata.get("settings_fingerprint"),
        trajectory_files=[],
    )


def condition_artifact_to_legacy_result(artifact: Any) -> DistanceAggregatedResult:
    """Convert a condition artifact to the established aggregate model."""

    if isinstance(artifact, DistanceAggregatedResult):
        return artifact
    if not isinstance(artifact, ConditionArtifact):
        return artifact
    metadata = artifact.metadata
    return DistanceAggregatedResult(
        config_hash=str(metadata.get("config_hash", "unknown")),
        polyzymd_version=str(metadata.get("polyzymd_version", get_polyzymd_version())),
        replicate=None,
        equilibration_time=float(metadata.get("equilibration_time", 0.0)),
        equilibration_unit=str(metadata.get("equilibration_unit", "ns")),
        selection_string=str(metadata.get("selection_string", "")),
        replicates=[int(rep) for rep in artifact.replicates],
        n_replicates=len(artifact.replicates),
        pair_results=[
            DistancePairAggregatedResult.model_validate(pair)
            for pair in artifact.payload.get("pairs", [])
        ],
        settings_fingerprint=metadata.get("settings_fingerprint"),
        source_result_files=[str(path) for path in metadata.get("source_result_files", [])],
    )


def mdanalysis_version() -> str:
    """Return the lazily imported MDAnalysis version string."""

    try:
        import MDAnalysis as mda
    except ImportError:
        return "unknown"
    return str(getattr(mda, "__version__", "unknown"))


def _summarize_distance_series(
    *,
    distances: NDArray[np.float64],
    pair_label: str,
    selection1: str,
    selection2: str,
    threshold: float | None,
    timestep_ps: float,
    n_frames_total: int,
) -> DistancePairPayload:
    """Summarize one distance series into the established payload schema."""

    from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

    n_frames_used = len(distances)
    mean_distance = float(np.mean(distances))
    std_distance = float(np.std(distances))
    median_distance = float(np.median(distances))
    min_distance = float(np.min(distances))
    max_distance = float(np.max(distances))
    fraction_below_threshold = (
        float(np.mean(distances < threshold)) if threshold is not None else None
    )
    histogram_counts, histogram_edges = np.histogram(distances, bins=50)
    kde_x: NDArray[np.float64] | None = None
    kde_y: NDArray[np.float64] | None = None
    kde_peak: float | None = None
    kde_bandwidth: float | None = None
    try:
        from scipy.stats import gaussian_kde

        if len(distances) > 10:
            kde = gaussian_kde(distances)
            kde_bandwidth = float(kde.factor) * float(np.std(distances))
            kde_x = np.linspace(max(0.0, min_distance - 0.5), max_distance + 0.5, 200)
            kde_y = np.asarray(kde(kde_x), dtype=np.float64)
            kde_peak = float(kde_x[int(np.argmax(kde_y))])
    except ImportError:
        pass
    except (ValueError, np.linalg.LinAlgError, RuntimeError) as exc:
        LOGGER.warning("KDE computation failed for %s: %s", pair_label, exc)

    sem_distance: float | None = None
    correlation_time: float | None = None
    correlation_time_unit: str | None = None
    n_independent_frames: int | None = None
    statistical_inefficiency: float | None = None
    autocorrelation_warning: str | None = None
    if len(distances) >= 20:
        try:
            tau_result = estimate_correlation_time(
                distances,
                timestep=timestep_ps,
                timestep_unit="ps",
                method="integration",
                n_frames=n_frames_used,
            )
            correlation_time = tau_result.tau
            correlation_time_unit = tau_result.tau_unit
            n_independent_frames = tau_result.n_independent
            statistical_inefficiency = tau_result.statistical_inefficiency
            autocorrelation_warning = tau_result.warning
            if n_independent_frames > 0:
                sem_distance = float(std_distance / np.sqrt(float(n_independent_frames)))
        except (ValueError, np.linalg.LinAlgError, RuntimeError) as exc:
            LOGGER.warning("Autocorrelation analysis failed for %s: %s", pair_label, exc)
            sem_distance = float(std_distance / np.sqrt(float(n_frames_used)))

    return DistancePairPayload(
        pair_label=pair_label,
        selection1=selection1,
        selection2=selection2,
        distances=np.asarray(distances, dtype=np.float64),
        mean_distance=mean_distance,
        std_distance=std_distance,
        median_distance=median_distance,
        min_distance=min_distance,
        max_distance=max_distance,
        sem_distance=sem_distance,
        correlation_time=correlation_time,
        correlation_time_unit=correlation_time_unit,
        n_independent_frames=n_independent_frames,
        statistical_inefficiency=statistical_inefficiency,
        autocorrelation_warning=autocorrelation_warning,
        threshold=threshold,
        fraction_below_threshold=fraction_below_threshold,
        histogram_edges=np.asarray(histogram_edges, dtype=np.float64),
        histogram_counts=np.asarray(histogram_counts, dtype=np.int64),
        kde_x=kde_x,
        kde_y=kde_y,
        kde_peak=kde_peak,
        kde_bandwidth=kde_bandwidth,
        n_frames_total=n_frames_total,
        n_frames_used=n_frames_used,
    )


def _write_replicate_sidecar(
    ctx: MDACollectorContext,
    job: MDAJobResult,
    payload: DistanceReplicatePayload,
) -> ArtifactSidecarRef:
    """Write the replicate pair x frame distance matrix sidecar."""

    distance_matrix = np.vstack([pair.distances for pair in payload.pair_payloads]).astype(
        np.float64
    )
    thresholds = np.asarray(
        [np.nan if pair.threshold is None else pair.threshold for pair in payload.pair_payloads],
        dtype=np.float64,
    )
    return ctx.artifact_store.write_npz_sidecar(
        Path("sidecars") / _sidecar_filename("pair_distances", 0),
        distance_matrix=distance_matrix,
        frames=np.asarray(payload.frames if payload.frames is not None else [], dtype=np.int64),
        time_ns=np.asarray(
            payload.time_ns if payload.time_ns is not None else [], dtype=np.float64
        ),
        thresholds=thresholds,
        pair_labels=np.asarray([pair.pair_label for pair in payload.pair_payloads]),
        selections_a=np.asarray([pair.selection1 for pair in payload.pair_payloads]),
        selections_b=np.asarray([pair.selection2 for pair in payload.pair_payloads]),
        metadata={"job_name": job.name, "kind": "distance_matrix", "layout": "pair_x_frame"},
    )


def _write_condition_sidecar(
    output_dir: Path,
    artifacts: Sequence[ReplicateArtifact],
    legacy_result: DistanceAggregatedResult,
) -> ArtifactSidecarRef:
    """Write pooled condition-level pair-distance sidecar arrays."""

    arrays_by_pair: dict[int, list[NDArray[np.float64]]] = {
        idx: [] for idx in range(len(legacy_result.pair_results))
    }
    replicate_ids: list[int] = []
    for artifact in artifacts:
        sidecar = _distance_matrix_sidecar(artifact)
        sidecar_path = ArtifactStore(
            output_dir.parent / f"run_{artifact.replicate}"
        ).validate_sidecar(sidecar)
        with np.load(sidecar_path) as npz_data:
            matrix = np.asarray(npz_data["distance_matrix"], dtype=np.float64)
            for pair_index in arrays_by_pair:
                arrays_by_pair[pair_index].append(matrix[pair_index])
        replicate_ids.append(int(artifact.replicate))

    sidecar_arrays: dict[str, Any] = {
        "replicates": np.asarray(replicate_ids, dtype=np.int64),
        "pair_labels": np.asarray([pair.pair_label for pair in legacy_result.pair_results]),
    }
    for pair_index, arrays in arrays_by_pair.items():
        sidecar_arrays[f"pair_{pair_index}_distances"] = np.concatenate(arrays).astype(np.float64)
    return ArtifactStore(output_dir).write_npz_sidecar(
        "sidecars/pooled_distances.npz",
        metadata={"kind": "pooled_distance_matrices", "layout": "pair_index_concatenated"},
        **sidecar_arrays,
    )


def _aggregate_legacy_result(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: DistancesSettings,
    equilibration: str,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> DistanceAggregatedResult:
    """Build the established aggregate model from replicate artifacts."""

    from polyzymd.analyses.shared.aggregation import aggregate_distance_pair_stats

    eq_value, eq_unit = parse_time_string(equilibration)
    legacy_results = [artifact_to_distance_result(artifact) for artifact in artifacts]
    first = legacy_results[0]
    metadata = DistanceResultMetadata(
        config_hash=first.config_hash,
        polyzymd_version=get_polyzymd_version(),
        replicate=None,
        equilibration_time=eq_value,
        equilibration_unit=eq_unit,
    )
    aggregated_pairs: list[DistancePairAggregatedResult] = []
    thresholds = settings.get_pair_thresholds()
    for pair_idx, pair in enumerate(first.pair_results):
        stats = aggregate_distance_pair_stats(legacy_results, pair_idx)
        aggregated_pairs.append(
            DistancePairAggregatedResult.from_aggregated_stats(
                metadata,
                pair,
                stats,
                replicates=replicates,
                threshold=thresholds[pair_idx] if pair_idx < len(thresholds) else None,
            )
        )
    return DistanceAggregatedResult.from_pair_results(
        metadata,
        aggregated_pairs,
        replicates=replicates,
        source_result_files=[],
        settings_fingerprint=settings_fingerprint,
        selection_string=_combined_selection(settings.get_pair_selections()),
    )


def _validate_and_order_artifacts(
    *,
    condition_label: str,
    expected_replicates: Sequence[int],
    settings: DistancesSettings,
    settings_fingerprint: str,
    artifacts: Sequence[ReplicateArtifact],
    analysis_dir: Path,
) -> list[ReplicateArtifact]:
    """Validate distance artifact identity, settings, pairs, and sidecars."""

    expected = [int(rep) for rep in expected_replicates]
    expected_pairs = settings.get_pair_selections()
    expected_thresholds = settings.get_pair_thresholds()
    by_replicate: dict[int, ReplicateArtifact] = {}
    for artifact in artifacts:
        if artifact.analysis_name != "distances":
            raise ValueError(f"Expected distances artifact, got {artifact.analysis_name!r}")
        if artifact.condition_label != condition_label:
            raise ValueError(
                f"Distances artifact condition mismatch: expected {condition_label!r}, "
                f"got {artifact.condition_label!r}"
            )
        if artifact.replicate in by_replicate:
            raise ValueError(f"Duplicate distances artifact for replicate {artifact.replicate}")
        if artifact.metadata.get("settings_fingerprint") != settings_fingerprint:
            raise ValueError(
                f"Distances artifact replicate {artifact.replicate} has settings fingerprint "
                f"{artifact.metadata.get('settings_fingerprint')}, expected {settings_fingerprint}"
            )
        _validate_pair_payloads(
            artifact,
            expected_pairs,
            expected_thresholds,
            expected_labels=settings.get_pair_labels(),
        )
        for sidecar in artifact.sidecars:
            ArtifactStore(analysis_dir / f"run_{artifact.replicate}").validate_sidecar(sidecar)
        by_replicate[int(artifact.replicate)] = artifact
    observed = sorted(by_replicate)
    if observed != sorted(expected):
        raise ValueError(
            f"Distance artifacts for condition '{condition_label}' do not match expected "
            f"replicates {expected}: found {observed}. Recompute missing replicates or clear "
            "stale caches before aggregating."
        )
    return [by_replicate[rep] for rep in expected]


def _validate_pair_payloads(
    artifact: ReplicateArtifact,
    expected_pairs: Sequence[tuple[str, str]],
    expected_thresholds: Sequence[float | None],
    expected_labels: Sequence[str] | None = None,
) -> None:
    """Validate pair payloads against active settings."""

    pairs = artifact.payload.get("pairs", [])
    if len(pairs) != len(expected_pairs):
        raise ValueError(
            f"Distances replicate {artifact.replicate} has {len(pairs)} pairs, "
            f"expected {len(expected_pairs)}"
        )
    labels = expected_labels if expected_labels is not None else [None] * len(expected_pairs)
    for idx, (payload, (selection_a, selection_b), threshold, label) in enumerate(
        zip(pairs, expected_pairs, expected_thresholds, labels, strict=True), start=1
    ):
        if payload.get("selection1") != selection_a or payload.get("selection2") != selection_b:
            raise ValueError(
                f"Distances replicate {artifact.replicate} pair {idx} selection mismatch"
            )
        if label is not None and payload.get("pair_label") != label:
            raise ValueError(f"Distances replicate {artifact.replicate} pair {idx} label mismatch")
        if payload.get("threshold") != threshold:
            raise ValueError(
                f"Distances replicate {artifact.replicate} pair {idx} threshold mismatch"
            )


def _pair_payload_to_json(payload: DistancePairPayload) -> dict[str, Any]:
    """Return JSON summary for one pair without raw distance arrays."""

    result = DistancePairResult.from_runner_payload(
        DistanceResultMetadata(
            config_hash="unknown",
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=0.0,
            equilibration_unit="ns",
        ),
        payload,
        store_distributions=False,
    ).model_dump(mode="json")
    for key in (
        "config_hash",
        "polyzymd_version",
        "replicate",
        "equilibration_time",
        "equilibration_unit",
        "selection_string",
        "analysis_type",
        "created_at",
        "distances",
    ):
        result.pop(key, None)
    return result


def _replicate_metrics(pair_payloads: Sequence[Mapping[str, Any]]) -> dict[str, float]:
    """Build replicate-level scalar metrics from pair summaries."""

    metrics: dict[str, float] = {}
    for pair in pair_payloads:
        label = str(pair.get("pair_label", "pair"))
        metrics[f"{label}.mean_distance"] = float(pair["mean_distance"])
        if pair.get("fraction_below_threshold") is not None:
            metrics[f"{label}.fraction_below_threshold"] = float(pair["fraction_below_threshold"])
    return metrics


def _condition_metrics(
    legacy_result: DistanceAggregatedResult,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, float]]]:
    """Build condition metric dictionaries from an aggregate model."""

    metrics: dict[str, dict[str, Any]] = {}
    replicate_metrics: dict[str, dict[str, float]] = {
        str(replicate): {} for replicate in legacy_result.replicates
    }
    for pair in legacy_result.pair_results:
        mean_name = f"{pair.pair_label}.mean_distance"
        values = [float(value) for value in pair.per_replicate_means]
        metrics[mean_name] = _metric_summary(mean_name, values, pair.overall_mean, pair.overall_sem)
        for replicate, value in zip(pair.replicates, values, strict=True):
            replicate_metrics.setdefault(str(replicate), {})[mean_name] = float(value)
        if pair.per_replicate_fractions_below is not None:
            frac_name = f"{pair.pair_label}.fraction_below_threshold"
            frac_values = [float(value) for value in pair.per_replicate_fractions_below]
            metrics[frac_name] = _metric_summary(
                frac_name,
                frac_values,
                float(pair.overall_fraction_below or 0.0),
                float(pair.sem_fraction_below or 0.0),
            )
            for replicate, value in zip(pair.replicates, frac_values, strict=True):
                replicate_metrics.setdefault(str(replicate), {})[frac_name] = float(value)
    return metrics, replicate_metrics


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


def _analysis_time_ns(
    analysis: Any,
    frames: NDArray[np.int64],
    timestep_ps: float,
) -> NDArray[np.float64]:
    """Return analysis times in nanoseconds with a frame-based fallback."""

    times_ps = np.asarray(getattr(analysis.results, "times_ps", []), dtype=np.float64)
    if times_ps.size == frames.size and times_ps.size > 0:
        return times_ps / 1000.0
    return frames.astype(np.float64) * float(timestep_ps) / 1000.0


def _n_frames_total(ctx: MDACollectorContext, job: MDAJobResult) -> int:
    """Return total trajectory frame count with a result-matrix fallback."""

    if ctx.frame_selection.n_frames_total is not None:
        return int(ctx.frame_selection.n_frames_total)
    matrix = np.asarray(getattr(job.results, "distance_matrix", []))
    if matrix.ndim == 2:
        return int(matrix.shape[1])
    return 0


def _effective_timestep_ps(time_ns: NDArray[np.float64], timestep_ps: float, step: int) -> float:
    """Return effective timestep between analyzed samples in picoseconds."""

    if time_ns.size >= 2:
        diffs_ps = np.diff(time_ns) * 1000.0
        finite = diffs_ps[np.isfinite(diffs_ps) & (diffs_ps > 0)]
        if finite.size:
            return float(np.median(finite))
    return float(timestep_ps) * float(step)


def _distance_matrix_sidecar(artifact: ReplicateArtifact) -> ArtifactSidecarRef:
    """Return the distance-matrix sidecar from a replicate artifact."""

    matches = [
        sidecar
        for sidecar in artifact.sidecars
        if sidecar.metadata.get("kind") == "distance_matrix"
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Distances replicate {artifact.replicate} requires exactly one distance matrix "
            f"sidecar, found {len(matches)}"
        )
    return matches[0]


def _combined_selection(pairs: Sequence[tuple[str, str]]) -> str:
    """Return the established combined selection string."""

    return "; ".join(f"({selection1} : {selection2})" for selection1, selection2 in pairs)


def _pair_specs_from_settings(settings: DistancesSettings) -> list[PairDistanceSpec]:
    """Build metadata-only pair specs from settings for result collection."""

    return [
        PairDistanceSpec(
            label=pair.label,
            selection_a=pair.selection_a,
            selection_b=pair.selection_b,
            atoms_a=None,
            atoms_b=None,
            mode_a=None,
            mode_b=None,
            threshold=threshold,
        )
        for pair, threshold in zip(settings.pairs, settings.get_pair_thresholds(), strict=True)
    ]


def _sidecar_filename(label: str, job_index: int) -> str:
    """Return a collision-resistant sidecar filename."""

    digest = sha256(label.encode("utf-8")).hexdigest()[:10]
    return f"{job_index:02d}_{digest}_distances.npz"


def _lazy_pair_label(selection_a: str, selection_b: str) -> str:
    """Import the public label helper lazily to avoid circular imports."""

    from polyzymd.analyses.distances import _make_pair_label

    return _make_pair_label(selection_a, selection_b)


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
