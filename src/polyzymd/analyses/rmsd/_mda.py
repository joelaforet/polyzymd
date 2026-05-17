"""MDAnalysis-native RMSD job helpers and artifact collectors."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence
from hashlib import sha256
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses._framework.cache_identity import compute_config_hash
from polyzymd.analyses._framework.results_base import get_polyzymd_version
from polyzymd.analyses.mda import (
    MDAAnalysisJob,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.plugin import frame_selection_payload, strict_json_payload
from polyzymd.analyses.rmsd._results import RMSDAggregatedResult, RMSDRunAggregatedResult
from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory
from polyzymd.analyses.shared.loader import parse_time_string
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.mda import ArtifactSidecarRef, FrameSelection, MDAReplicateJobContext
    from polyzymd.analyses.rmsd import RMSDRunSettings, RMSDSettings

LOGGER = logging.getLogger(__name__)


def build_rmsd_jobs(
    ctx: MDAReplicateJobContext, runs: Sequence[RMSDRunSettings]
) -> list[MDAAnalysisJob]:
    """Build one MDAnalysis RMSD job per configured run.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        Framework-provided MDA replicate context.
    runs : sequence of RMSDRunSettings
        RMSD run settings from the user-facing YAML configuration.

    Returns
    -------
    list of MDAAnalysisJob
        Jobs that construct a fresh universe during execution.
    """

    jobs: list[MDAAnalysisJob] = []
    for run in runs:
        run_payload = run.model_dump(mode="json")
        jobs.append(
            MDAAnalysisJob(
                name=run.label,
                analysis_factory=_make_rmsd_analysis_factory(ctx, run),
                frame_selection=ctx.frame_selection,
                backend_policy=ctx.backend_policy,
                universe_policy=MDAUniversePolicy(
                    condition_label=ctx.replicate_context.condition.label,
                    replicate=ctx.replicate,
                    provenance=ctx.universe_policy.provenance,
                    metadata={
                        **ctx.universe_policy.metadata,
                        "fresh_universe_per_run": True,
                        "rmsd_run": run_payload,
                    },
                ),
            )
        )
    return jobs


class RMSDArtifactCollector:
    """Collect completed MDAnalysis RMSD jobs into a replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Map completed RMSD jobs to an artifact with NPZ sidecars.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context for one replicate.
        completed_jobs : sequence of MDAJobResult
            Completed MDAnalysis RMSD jobs.

        Returns
        -------
        ReplicateArtifact
            PolyzyMD-owned replicate artifact with scalar metrics and sidecars.
        """

        config_hash = compute_config_hash(ctx.replicate_context.sim_config)
        eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
        run_payloads: list[dict[str, Any]] = []
        sidecars: list[ArtifactSidecarRef] = []
        metrics: dict[str, float] = {}

        for job_index, job in enumerate(completed_jobs):
            run_payload, run_sidecar = _collect_run(ctx, job, job_index=job_index)
            run_payloads.append(run_payload)
            sidecars.append(run_sidecar)
            metrics[f"{job.name}.mean_rmsd"] = run_payload["mean_rmsd"]

        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "runs": run_payloads,
                "metrics": metrics,
                "replicate_metrics": metrics,
                "n_runs": len(run_payloads),
                "n_frames_total": ctx.frame_selection.n_frames_total,
                "n_frames_used": ctx.frame_selection.n_frames_selected,
            },
            sidecars=sidecars,
            provenance={
                "source": "mda_rmsd_jobs",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
            },
            metadata={
                "result_kind": "rmsd_mda_replicate",
                "settings_fingerprint": ctx.settings_fingerprint,
                "config_hash": config_hash,
                "polyzymd_version": get_polyzymd_version(),
                "equilibration_time": eq_value,
                "equilibration_unit": eq_unit,
            },
            warnings=list(ctx.warnings),
        )


def aggregate_rmsd_artifacts(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: RMSDSettings,
    equilibration: str,
    output_dir: Path,
    result_path: Path,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> Any:
    """Aggregate RMSD replicate artifacts into a condition artifact.

    Parameters
    ----------
    condition_label : str
        Label for the condition being aggregated.
    replicates : sequence of int
        Replicate IDs represented in ``artifacts``.
    settings : RMSDSettings
        Active RMSD settings.
    equilibration : str
        Equilibration string from the framework context.
    output_dir : Path
        Aggregated output directory.
    result_path : Path
        Canonical condition artifact path.
    artifacts : sequence of ReplicateArtifact
        Per-replicate RMSD artifacts.
    settings_fingerprint : str
        Active settings fingerprint.

    Returns
    -------
    ConditionArtifact
        Aggregated RMSD condition artifact.
    """

    from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact

    if not artifacts:
        raise ValueError(
            f"RMSD aggregation for condition '{condition_label}' requires at least one "
            "replicate artifact. No replicate inputs were provided."
        )

    run_labels = [run.label for run in settings.runs]
    ordered_artifacts = _validate_and_order_artifacts(
        condition_label=condition_label,
        expected_replicates=replicates,
        run_labels=run_labels,
        settings_fingerprint=settings_fingerprint,
        artifacts=artifacts,
    )
    eq_value, eq_unit = parse_time_string(equilibration)
    first = ordered_artifacts[0]
    config_hash = str(first.metadata.get("config_hash", "unknown"))
    replicate_order = [int(rep) for rep in replicates]

    if len(replicate_order) == 1:
        LOGGER.warning(
            "Only one replicate available for RMSD aggregation in condition '%s'; "
            "replicate-level SEM is reported as 0.0",
            condition_label,
        )

    aggregated_runs: list[dict[str, Any]] = []
    metrics: dict[str, dict[str, Any]] = {}
    replicate_metrics: dict[str, dict[str, float]] = {str(rep): {} for rep in replicate_order}

    for run_label in run_labels:
        run_entries = [_run_payload_for(artifact, run_label) for artifact in ordered_artifacts]
        per_means = [_finite_float(entry["mean_rmsd"], field="mean_rmsd") for entry in run_entries]
        per_stds = [_finite_float(entry["std_rmsd"], field="std_rmsd") for entry in run_entries]
        per_medians = [
            _finite_float(entry["median_rmsd"], field="median_rmsd") for entry in run_entries
        ]
        per_convergence_times = [entry.get("convergence_time_ns") for entry in run_entries]
        per_assessable = [bool(entry.get("convergence_assessable", False)) for entry in run_entries]

        n_converged = sum(time is not None for time in per_convergence_times)
        n_assessable = sum(per_assessable)
        finite_convergence_times = [
            float(time)
            for time in per_convergence_times
            if time is not None and math.isfinite(time)
        ]
        mean_stats = compute_sem(per_means)
        overall_median = float(np.median(np.asarray(per_medians, dtype=np.float64)))
        metric_name = f"{run_label}.mean_rmsd"

        for replicate, mean_value in zip(replicate_order, per_means):
            replicate_metrics[str(replicate)][metric_name] = mean_value

        metrics[metric_name] = {
            "name": metric_name,
            "values": per_means,
            "mean": mean_stats.mean,
            "sem": mean_stats.sem,
            "std": (
                float(np.std(np.asarray(per_means, dtype=np.float64), ddof=1))
                if len(per_means) > 1
                else 0.0
            ),
            "n": len(per_means),
        }

        template = run_entries[0]
        aggregated_runs.append(
            {
                "config_hash": config_hash,
                "polyzymd_version": get_polyzymd_version(),
                "replicate": None,
                "equilibration_time": eq_value,
                "equilibration_unit": eq_unit,
                "selection_string": template["selection"],
                "replicates": replicate_order,
                "n_replicates": len(replicate_order),
                "run_label": run_label,
                "selection": template["selection"],
                "alignment_selection": template["alignment_selection"],
                "overall_mean": mean_stats.mean,
                "overall_sem": mean_stats.sem,
                "overall_median": overall_median,
                "per_replicate_means": per_means,
                "per_replicate_stds": per_stds,
                "per_replicate_medians": per_medians,
                "per_replicate_convergence_times_ns": per_convergence_times,
                "per_replicate_convergence_assessable": per_assessable,
                "n_converged_replicates": n_converged,
                "n_assessable_replicates": n_assessable,
                "convergence_fraction": (
                    float(n_converged) / float(n_assessable) if n_assessable > 0 else 0.0
                ),
                "all_converged": n_assessable > 0 and n_converged == n_assessable,
                "mean_convergence_time_ns": (
                    float(np.mean(np.asarray(finite_convergence_times, dtype=np.float64)))
                    if finite_convergence_times
                    else None
                ),
                "median_convergence_time_ns": (
                    float(np.median(np.asarray(finite_convergence_times, dtype=np.float64)))
                    if finite_convergence_times
                    else None
                ),
            }
        )

    source_result_files = _source_result_files(output_dir, replicate_order)
    artifact = ConditionArtifact(
        analysis_name="rmsd",
        condition_label=condition_label,
        replicates=replicate_order,
        payload={
            "runs": aggregated_runs,
            "metrics": metrics,
            "replicate_metrics": replicate_metrics,
            "n_replicates": len(replicate_order),
        },
        provenance={
            "source": "rmsd_replicate_artifact_aggregation",
            "frame_selection": first.provenance.get("frame_selection"),
        },
        metadata={
            "result_kind": "rmsd_mda_condition",
            "settings_fingerprint": settings_fingerprint,
            "config_hash": config_hash,
            "polyzymd_version": get_polyzymd_version(),
            "equilibration_time": eq_value,
            "equilibration_unit": eq_unit,
            "selection_string": "; ".join(run.selection for run in settings.runs),
            "source_result_files": source_result_files,
            "n_replicates": len(replicate_order),
        },
        source_replicates=[
            {"replicate": replicate, "path": path} for replicate, path in source_result_files
        ],
        warnings=_combined_warnings(ordered_artifacts),
    )
    ArtifactStore(result_path.parent).write_condition_result(artifact, result_path.name)
    return artifact


def condition_artifact_to_legacy_result(artifact: Any) -> RMSDAggregatedResult:
    """Convert an RMSD condition artifact to the established comparison model.

    Parameters
    ----------
    artifact : Any
        Condition artifact produced by RMSD aggregation.

    Returns
    -------
    RMSDAggregatedResult
        In-memory compatibility model used by existing comparison formatters.
    """

    if isinstance(artifact, RMSDAggregatedResult):
        return artifact
    run_payloads = artifact.payload.get("runs") if hasattr(artifact, "payload") else None
    if not isinstance(run_payloads, list):
        raise TypeError("RMSD condition artifact is missing payload['runs']")
    run_results = [RMSDRunAggregatedResult.model_validate(payload) for payload in run_payloads]
    metadata = dict(getattr(artifact, "metadata", {}) or {})
    return RMSDAggregatedResult(
        config_hash=str(metadata.get("config_hash", "unknown")),
        polyzymd_version=str(metadata.get("polyzymd_version", get_polyzymd_version())),
        replicate=None,
        equilibration_time=float(metadata.get("equilibration_time", 0.0)),
        equilibration_unit=str(metadata.get("equilibration_unit", "ns")),
        selection_string=str(metadata.get("selection_string", "")),
        replicates=[int(rep) for rep in artifact.replicates],
        n_replicates=len(artifact.replicates),
        run_results=run_results,
        settings_fingerprint=metadata.get("settings_fingerprint"),
        source_result_files=[str(entry) for entry in metadata.get("source_result_files", [])],
    )


def _make_rmsd_analysis_factory(ctx: MDAReplicateJobContext, run: RMSDRunSettings) -> Any:
    """Create a zero-argument factory for one fresh-universe RMSD analysis."""

    def factory() -> Any:
        from polyzymd.analyses.shared.loader import TrajectoryLoader

        loader = TrajectoryLoader(ctx.replicate_context.sim_config)
        universe = loader.load_universe(ctx.replicate, cache=False)
        return _build_rmsd_analysis(universe=universe, run=run, frame_selection=ctx.frame_selection)

    return factory


def _build_rmsd_analysis(
    *, universe: Any, run: RMSDRunSettings, frame_selection: FrameSelection
) -> Any:
    """Build an MDAnalysis RMSD analysis object after reference preparation."""

    from MDAnalysis.analysis.rms import RMSD

    start, stop, step = _slice_bounds(frame_selection)
    centroid_selection = run.centroid_selection
    if run.reference_mode == "centroid" and centroid_selection is None:
        centroid_selection = run.alignment_selection
        LOGGER.info(
            "Run '%s': centroid_selection not set, using alignment_selection='%s'",
            run.label,
            centroid_selection,
        )

    reference_frame_1indexed = run.reference_frame + 1 if run.reference_mode == "frame" else None
    alignment_config = AlignmentConfig(
        enabled=True,
        reference_mode=run.reference_mode,
        reference_frame=reference_frame_1indexed,
        selection=run.alignment_selection,
        centroid_selection=centroid_selection or run.alignment_selection,
        reference_file=(Path(run.reference_file) if run.reference_file is not None else None),
    )
    ref_frame_idx = align_trajectory(
        universe,
        alignment_config,
        start_frame=start,
        stop_frame=stop,
        step_frame=step,
    )

    atom_group = universe.select_atoms(run.selection)
    if len(atom_group) == 0:
        raise ValueError(f"Run '{run.label}' selection matched no atoms: {run.selection!r}")

    _reference_universe, reference_atom_group = _build_reference_structure(
        universe=universe,
        atom_group=atom_group,
        run=run,
        start=start,
        stop=stop,
        step=step,
        ref_frame_idx=ref_frame_idx,
    )
    rmsd_analysis = RMSD(atom_group, reference=reference_atom_group, select="all", ref_frame=0)
    rmsd_analysis._polyzymd_rmsd_metadata = {
        "reference_frame": ref_frame_idx + 1 if ref_frame_idx is not None else None,
        "reference_file": run.reference_file,
    }
    return rmsd_analysis


def _build_reference_structure(
    *,
    universe: Any,
    atom_group: Any,
    run: RMSDRunSettings,
    start: int,
    stop: int,
    step: int,
    ref_frame_idx: int | None,
) -> tuple[Any, Any]:
    """Materialize a reference structure for MDAnalysis RMSD."""

    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    if run.reference_mode in {"centroid", "frame"}:
        if ref_frame_idx is None:
            raise ValueError(
                f"Run '{run.label}' expected a reference frame for mode '{run.reference_mode}'"
            )
        universe.trajectory[ref_frame_idx]
        ref_positions = atom_group.positions.copy().astype(np.float64)
    elif run.reference_mode == "average":
        mean_analysis = _build_mean_position_analysis(atom_group).run(
            start=start,
            stop=stop,
            step=step,
        )
        ref_positions = np.asarray(mean_analysis.results.mean_positions, dtype=np.float64)
    elif run.reference_mode == "external":
        if run.reference_file is None:
            raise ValueError(
                f"Run '{run.label}' requires reference_file when reference_mode='external'"
            )
        ref_path = Path(run.reference_file)
        ref_universe = mda.Universe(str(ref_path))
        ref_atoms = ref_universe.select_atoms(run.selection)
        if len(ref_atoms) == 0:
            raise ValueError(
                f"Run '{run.label}' external PDB '{ref_path.name}' has no atoms matching "
                f"selection {run.selection!r}."
            )
        if len(ref_atoms) != len(atom_group):
            raise ValueError(
                f"Run '{run.label}' atom count mismatch between trajectory ({len(atom_group)}) "
                f"and external PDB ({len(ref_atoms)}) for selection {run.selection!r}."
            )
        ref_positions = ref_atoms.positions.copy().astype(np.float64)
    else:
        raise ValueError(f"Unsupported RMSD reference_mode: {run.reference_mode!r}")

    reference_universe = mda.Merge(atom_group)
    reference_universe.load_new(ref_positions[np.newaxis, :, :], format=MemoryReader)
    return reference_universe, reference_universe.atoms


def _build_mean_position_analysis(atoms: Any) -> Any:
    """Build an AnalysisBase object that accumulates mean coordinates."""

    from MDAnalysis.analysis.base import AnalysisBase

    class _MeanPositionAnalysis(AnalysisBase):
        """Accumulate average coordinates over the analyzed frames."""

        def __init__(self, atom_group: Any) -> None:
            super().__init__(atom_group.universe.trajectory)
            self._atom_group = atom_group

        def _prepare(self) -> None:
            """Initialize accumulation arrays before frame iteration."""

            self._position_sum = np.zeros_like(self._atom_group.positions, dtype=np.float64)
            self._n_frames = 0

        def _single_frame(self) -> None:
            """Accumulate coordinates for the current frame."""

            self._position_sum += self._atom_group.positions.astype(np.float64)
            self._n_frames += 1

        def _conclude(self) -> None:
            """Store mean positions after frame iteration."""

            if self._n_frames == 0:
                raise ValueError("Mean-position analysis requires at least one frame")
            self.results.mean_positions = self._position_sum / float(self._n_frames)

    return _MeanPositionAnalysis(atoms)


def _collect_run(
    ctx: MDACollectorContext, job: MDAJobResult, *, job_index: int
) -> tuple[dict[str, Any], ArtifactSidecarRef]:
    """Collect one completed RMSD job into summary and sidecar payloads."""

    from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time
    from polyzymd.analyses.shared.convergence import find_convergence_time

    run_settings = _run_settings_from_job(job)
    rmsd_table = np.asarray(job.results.rmsd, dtype=np.float64)
    if rmsd_table.ndim != 2 or rmsd_table.shape[1] < 3 or rmsd_table.shape[0] == 0:
        raise ValueError(f"RMSD job '{job.name}' produced an invalid results.rmsd array")
    rmsd_values = np.asarray(rmsd_table[:, 2], dtype=np.float64)
    if not np.all(np.isfinite(rmsd_values)):
        raise ValueError(f"RMSD job '{job.name}' produced non-finite RMSD values")

    raw_timestep_ps = float(ctx.frame_selection.timestep_ps or 1.0)
    frame_stride = int(ctx.frame_selection.step or 1)
    frames = np.asarray(rmsd_table[:, 0], dtype=np.int64)
    time_ps = _time_ps_from_rmsd_table(rmsd_table, job_name=job.name)
    time_ns = time_ps / 1000.0
    effective_timestep_ps = _effective_timestep_ps(time_ps, raw_timestep_ps, frame_stride)
    stats = _summarize_rmsd_values(rmsd_values)

    sem_rmsd: float | None = None
    correlation_time: float | None = None
    correlation_time_unit: str | None = None
    n_independent_frames: int | None = None
    statistical_inefficiency: float | None = None
    autocorrelation_warning: str | None = None
    if len(rmsd_values) >= 20:
        tau_result = estimate_correlation_time(
            rmsd_values,
            timestep=effective_timestep_ps,
            timestep_unit="ps",
            method="integration",
            n_frames=len(rmsd_values),
        )
        correlation_time = tau_result.tau
        correlation_time_unit = tau_result.tau_unit
        n_independent_frames = tau_result.n_independent
        statistical_inefficiency = tau_result.statistical_inefficiency
        autocorrelation_warning = tau_result.warning
        if n_independent_frames > 0:
            sem_rmsd = float(stats["std_rmsd"] / np.sqrt(float(n_independent_frames)))

    convergence_result = find_convergence_time(
        time_ns,
        rmsd_values,
        window_size_ns=float(run_settings["convergence_window_size_ns"]),
        step_size_ns=float(run_settings["convergence_step_size_ns"]),
        slope_threshold=float(run_settings["convergence_slope_threshold"]),
        sustained_for_ns=float(run_settings["convergence_sustained_for_ns"]),
    )
    sidecar = ctx.artifact_store.write_npz_sidecar(
        Path("sidecars") / _sidecar_filename(job.name, job_index),
        rmsd_values=rmsd_values,
        time_ns=time_ns,
        frames=frames,
        convergence_window_start_ns=np.asarray(
            convergence_result.window_start_times_ns, dtype=np.float64
        ),
        convergence_window_mean_rmsd=np.asarray(
            convergence_result.window_mean_values, dtype=np.float64
        ),
        convergence_slope_time_ns=np.asarray(convergence_result.slope_times_ns, dtype=np.float64),
        convergence_slope=np.asarray(convergence_result.slopes, dtype=np.float64),
        convergence_converged=np.asarray(convergence_result.converged, dtype=np.bool_),
        convergence_time_ns=np.asarray(
            (
                np.nan
                if convergence_result.convergence_time_ns is None
                else convergence_result.convergence_time_ns
            ),
            dtype=np.float64,
        ),
        raw_timestep_ps=np.asarray(raw_timestep_ps, dtype=np.float64),
        frame_stride=np.asarray(frame_stride, dtype=np.int64),
        effective_timestep_ps=np.asarray(effective_timestep_ps, dtype=np.float64),
        metadata={"run_label": job.name, "kind": "rmsd_timeseries"},
    )
    reference_metadata = getattr(job.analysis, "_polyzymd_rmsd_metadata", {})
    payload = {
        "run_label": job.name,
        "selection": str(run_settings["selection"]),
        "alignment_selection": str(run_settings["alignment_selection"]),
        "reference_mode": str(run_settings["reference_mode"]),
        "reference_frame": reference_metadata.get("reference_frame"),
        "reference_file": reference_metadata.get("reference_file"),
        **stats,
        "sem_rmsd": sem_rmsd,
        "correlation_time": correlation_time,
        "correlation_time_unit": correlation_time_unit,
        "n_independent_frames": n_independent_frames,
        "statistical_inefficiency": statistical_inefficiency,
        "autocorrelation_warning": autocorrelation_warning,
        "converged": bool(convergence_result.converged),
        "convergence_assessable": bool(convergence_result.assessable),
        "convergence_time_ns": convergence_result.convergence_time_ns,
        "convergence_message": convergence_result.message,
        "n_frames_total": ctx.frame_selection.n_frames_total,
        "n_frames_used": len(rmsd_values),
        "time_unit": "ns",
        "timestep_ps": effective_timestep_ps,
        "raw_timestep_ps": raw_timestep_ps,
        "frame_stride": frame_stride,
        "sidecar": sidecar.model_dump(mode="json"),
        "npz_path": sidecar.path,
    }
    return payload, sidecar


def _run_settings_from_job(job: MDAJobResult) -> Mapping[str, Any]:
    """Return RMSD run settings stored in job metadata."""

    metadata = dict(job.universe_policy.metadata)
    run_settings = metadata.get("rmsd_run")
    if not isinstance(run_settings, Mapping):
        raise ValueError(f"RMSD job '{job.name}' is missing RMSD run metadata")
    return run_settings


def _summarize_rmsd_values(rmsd_values: NDArray[np.float64]) -> dict[str, float]:
    """Compute finite scalar summaries from one RMSD timeseries."""

    return {
        "mean_rmsd": float(np.mean(rmsd_values)),
        "std_rmsd": float(np.std(rmsd_values, ddof=0)),
        "median_rmsd": float(np.median(rmsd_values)),
        "min_rmsd": float(np.min(rmsd_values)),
        "max_rmsd": float(np.max(rmsd_values)),
        "final_rmsd": float(rmsd_values[-1]),
    }


def _time_ps_from_rmsd_table(
    rmsd_table: NDArray[np.float64], *, job_name: str
) -> NDArray[np.float64]:
    """Return MDAnalysis-reported time values in picoseconds.

    Parameters
    ----------
    rmsd_table : NDArray[np.float64]
        ``RMSD.results.rmsd`` table with MDAnalysis frame, time, and RMSD columns.
    job_name : str
        Job label used in diagnostics.

    Returns
    -------
    NDArray[np.float64]
        Time values from the MDAnalysis results table in picoseconds.
    """

    time_ps = np.asarray(rmsd_table[:, 1], dtype=np.float64)
    if time_ps.ndim != 1 or len(time_ps) == 0 or not np.all(np.isfinite(time_ps)):
        raise ValueError(f"RMSD job '{job_name}' produced invalid results.rmsd time values")
    return time_ps


def _effective_timestep_ps(
    time_ps: NDArray[np.float64], raw_timestep_ps: float, frame_stride: int
) -> float:
    """Estimate effective sample spacing from MDAnalysis time values.

    Parameters
    ----------
    time_ps : NDArray[np.float64]
        MDAnalysis-reported time values in picoseconds.
    raw_timestep_ps : float
        Fallback raw trajectory timestep in picoseconds.
    frame_stride : int
        Frame stride used for the job.

    Returns
    -------
    float
        Median spacing in picoseconds when available, otherwise the legacy
        raw timestep multiplied by frame stride.
    """

    if len(time_ps) < 2:
        return float(raw_timestep_ps) * float(frame_stride)
    deltas = np.diff(time_ps)
    finite_positive = deltas[np.isfinite(deltas) & (deltas > 0.0)]
    if finite_positive.size == 0:
        return float(raw_timestep_ps) * float(frame_stride)
    return float(np.median(finite_positive))


def _slice_bounds(frame_selection: FrameSelection) -> tuple[int, int, int]:
    """Return slice bounds for alignment and reference construction."""

    if frame_selection.frames is not None:
        raise ValueError("RMSD jobs require start/stop/step frame selection")
    if frame_selection.n_frames_total is None:
        raise ValueError("RMSD jobs require known total frame count")
    return (
        0 if frame_selection.start is None else int(frame_selection.start),
        (
            frame_selection.n_frames_total
            if frame_selection.stop is None
            else int(frame_selection.stop)
        ),
        1 if frame_selection.step is None else int(frame_selection.step),
    )


def _validate_and_order_artifacts(
    *,
    condition_label: str,
    expected_replicates: Sequence[int],
    run_labels: Sequence[str],
    settings_fingerprint: str,
    artifacts: Sequence[ReplicateArtifact],
) -> list[ReplicateArtifact]:
    """Validate RMSD replicate artifact coverage and ordering."""

    expected = [int(rep) for rep in expected_replicates]
    by_rep: dict[int, ReplicateArtifact] = {}
    for artifact in artifacts:
        if artifact.analysis_name != "rmsd" or artifact.condition_label != condition_label:
            raise ValueError("RMSD aggregation received an artifact with mismatched identity")
        if artifact.metadata.get("settings_fingerprint") != settings_fingerprint:
            raise ValueError(
                f"RMSD replicate {artifact.replicate} settings fingerprint mismatch: "
                f"stored {artifact.metadata.get('settings_fingerprint')}, expected {settings_fingerprint}"
            )
        if artifact.replicate in by_rep:
            raise ValueError(f"Duplicate RMSD artifact for replicate {artifact.replicate}")
        observed_labels = {
            _run_payload_label(payload) for payload in artifact.payload.get("runs", [])
        }
        missing_labels = sorted(set(run_labels) - observed_labels)
        if missing_labels:
            raise ValueError(
                f"RMSD replicate {artifact.replicate} is missing configured run(s) {missing_labels}"
            )
        by_rep[artifact.replicate] = artifact
    observed = sorted(by_rep)
    if observed != sorted(expected):
        raise ValueError(
            f"RMSD aggregation for condition '{condition_label}' is incomplete. Expected "
            f"replicate artifacts for {sorted(expected)}, found {observed}."
        )
    return [by_rep[rep] for rep in expected]


def _run_payload_for(artifact: ReplicateArtifact, run_label: str) -> Mapping[str, Any]:
    """Return one run payload from a replicate artifact."""

    for payload in artifact.payload.get("runs", []):
        if _run_payload_label(payload) == run_label:
            return payload
    raise ValueError(f"RMSD replicate {artifact.replicate} missing run {run_label!r}")


def _run_payload_label(payload: Any) -> str:
    """Return a run label from a payload mapping."""

    if not isinstance(payload, Mapping):
        raise ValueError("RMSD run payload must be a mapping")
    label = payload.get("run_label")
    if not isinstance(label, str) or not label:
        raise ValueError("RMSD run payload is missing a non-empty run_label")
    return label


def _finite_float(value: Any, *, field: str) -> float:
    """Validate and return one finite floating-point field."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"RMSD field {field!r} must be a finite scalar")
    scalar = float(value)
    if not math.isfinite(scalar):
        raise ValueError(f"RMSD field {field!r} is non-finite")
    return scalar


def _source_result_files(output_dir: Path, replicates: Sequence[int]) -> list[tuple[int, str]]:
    """Return canonical source result paths for aggregate provenance."""

    analysis_dir = output_dir.parent
    return [
        (int(replicate), str(analysis_dir / f"run_{replicate}" / "result.json"))
        for replicate in replicates
    ]


def _combined_warnings(artifacts: Sequence[ReplicateArtifact]) -> list[str]:
    """Combine warnings from replicate artifacts in first-seen order."""

    warnings: list[str] = []
    seen: set[str] = set()
    for artifact in artifacts:
        for warning in artifact.warnings:
            text = str(warning)
            if text not in seen:
                seen.add(text)
                warnings.append(text)
    return warnings


def _safe_label(value: str) -> str:
    """Return a filesystem-safe label token."""

    return value.replace(" ", "_").replace("-", "_").replace("/", "_").lower()


def _sidecar_filename(run_label: str, job_index: int) -> str:
    """Return a collision-proof sidecar filename for one RMSD run.

    Parameters
    ----------
    run_label : str
        Original user-facing run label.
    job_index : int
        Zero-based job order from the completed job sequence.

    Returns
    -------
    str
        Stable sidecar filename including job order, sanitized label, and a hash
        of the original label.
    """

    label_hash = sha256(run_label.encode("utf-8")).hexdigest()[:12]
    return f"rmsd_{job_index:03d}_{_safe_label(run_label)}_{label_hash}_timeseries.npz"
