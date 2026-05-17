"""MDAnalysis-compatible SASA jobs and canonical artifacts."""

from __future__ import annotations

import hashlib
import logging
import math
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses._results_base import get_polyzymd_version
from polyzymd.analyses.mda import (
    ArtifactStore,
    ConditionArtifact,
    FrameSelection,
    MDAAnalysisJob,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.plugin import frame_selection_payload, strict_json_payload
from polyzymd.analyses.shared.config_hash import compute_config_hash
from polyzymd.analyses.shared.loader import parse_time_string
from polyzymd.analyses.shared.sasa import (
    NM2_TO_A2,
    SASAComputationResult,
    resolve_selection_indices,
    validate_target_subset,
)
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.mda import ArtifactSidecarRef, MDAReplicateJobContext
    from polyzymd.analyses.sasa import SASASettings

LOGGER = logging.getLogger(__name__)

MEAN_SASA_METRIC = "mean_sasa"
SASA_SIDECAR_PREFIX = "sidecars/sasa"
SASA_METRIC_METADATA: dict[str, Any] = {
    "higher_is_better": False,
    "direction_labels": ("shielding", "unchanged", "exposure"),
    "units": "A^2",
    "description": "Mean solvent-accessible surface area over selected frames",
    "statistical_policy": "mean_based",
}


def build_sasa_jobs(ctx: MDAReplicateJobContext, settings: SASASettings) -> list[MDAAnalysisJob]:
    """Build AnalysisBase-compatible SASA jobs for one replicate.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        Framework-provided MDAnalysis replicate context.
    settings : SASASettings
        User-facing SASA settings.

    Returns
    -------
    list of MDAAnalysisJob
        One job per configured SASA run with run-specific stride folded into
        the frame selection.
    """

    jobs: list[MDAAnalysisJob] = []
    backend_policy = getattr(ctx, "backend_policy", None)
    for order, run in enumerate(settings.runs):
        context_selection = run.context_selection or run.target_selection
        run_frame_selection = _frame_selection_with_stride(ctx.frame_selection, run.stride)
        job_name = _job_name(run.label, order)
        jobs.append(
            MDAAnalysisJob(
                name=job_name,
                analysis=build_sasa_surface_area_analysis(
                    universe=ctx.universe,
                    run_label=run.label,
                    target_selection=run.target_selection,
                    context_selection=context_selection,
                    probe_radius_nm=settings.probe_radius_nm,
                    n_sphere_points=settings.n_sphere_points,
                    chunk_size=settings.chunk_size,
                    timestep_ps=ctx.frame_selection.timestep_ps,
                ),
                frame_selection=run_frame_selection,
                **({"backend_policy": backend_policy} if backend_policy is not None else {}),
                universe_policy=MDAUniversePolicy(
                    condition_label=ctx.replicate_context.condition.label,
                    replicate=ctx.replicate,
                    provenance=ctx.universe_policy.provenance,
                    metadata={
                        **ctx.universe_policy.metadata,
                        "sasa_run_label": run.label,
                        "sasa_run_order": order,
                        "sasa_target_selection": run.target_selection,
                        "sasa_context_selection": context_selection,
                        "sasa_probe_radius_nm": settings.probe_radius_nm,
                        "sasa_n_sphere_points": settings.n_sphere_points,
                        "sasa_chunk_size": settings.chunk_size,
                        "sasa_stride": run.stride,
                    },
                ),
            )
        )
    return jobs


def build_sasa_surface_area_analysis(
    *,
    universe: Any,
    run_label: str,
    target_selection: str,
    context_selection: str,
    probe_radius_nm: float,
    n_sphere_points: int,
    chunk_size: int,
    timestep_ps: float | None,
) -> Any:
    """Build an ``AnalysisBase`` wrapper around the MDTraj SASA kernel.

    MDAnalysis owns frame iteration through ``AnalysisBase``. The wrapper keeps
    only one bounded coordinate chunk in memory and evaluates MDTraj
    Shrake-Rupley whenever that buffer reaches ``chunk_size``.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe for one replicate.
    run_label : str
        User-facing run label.
    target_selection : str
        Atom selection whose SASA is reported.
    context_selection : str
        Atom selection considered by the surface-area computation.
    probe_radius_nm : float
        Shrake-Rupley probe radius in nm.
    n_sphere_points : int
        Number of Shrake-Rupley sphere points.
    chunk_size : int
        Number of selected frames per MDTraj kernel call.
    timestep_ps : float or None
        Raw trajectory timestep in picoseconds when known.

    Returns
    -------
    Any
        AnalysisBase-compatible object with raw SASA arrays in ``results``.
    """

    from MDAnalysis.analysis.base import AnalysisBase

    class SASASurfaceAreaAnalysis(AnalysisBase):
        """Collect bounded coordinate chunks and compute target SASA arrays."""

        def __init__(self) -> None:
            target_atoms, target_indices = resolve_selection_indices(
                universe,
                target_selection,
                role="target",
                run_label=run_label,
            )
            context_atoms, context_indices = resolve_selection_indices(
                universe,
                context_selection,
                role="context",
                run_label=run_label,
            )
            if target_indices.size > 0 and context_indices.size > 0:
                validate_target_subset(
                    target_indices,
                    context_indices,
                    run_label=run_label,
                    target_selection=target_selection,
                    context_selection=context_selection,
                )
            super().__init__(universe.trajectory)
            self._target_atoms = target_atoms
            self._context_atoms = context_atoms
            self._target_indices = target_indices
            self._context_indices = context_indices
            self._raw_timestep_ps = float(timestep_ps) if timestep_ps is not None else None
            self._zero_atom_selection = target_indices.size == 0 or context_indices.size == 0
            self._chunk_size = int(chunk_size)
            if self._chunk_size <= 0:
                raise ValueError("chunk_size must be >= 1")
            if not self._zero_atom_selection:
                self._topology = _mdtraj_topology_from_atom_group(context_atoms)
                self._target_local_indices = _target_local_indices(target_indices, context_indices)
                self._residue_items = _residue_items(target_atoms, self._target_local_indices)

        def _prepare(self) -> None:
            """Initialize coordinate buffering for selected frames."""

            self._position_chunk: list[NDArray[np.float32]] = []
            self._atom_sasa_chunks: list[NDArray[np.float64]] = []
            self._residue_sasa_chunks: list[NDArray[np.float64]] = []
            self._total_sasa_chunks: list[NDArray[np.float64]] = []
            self.results.max_buffered_coordinate_frames = 0
            self.results.run_label = run_label
            self.results.target_selection = target_selection
            self.results.context_selection = context_selection
            self.results.probe_radius_nm = float(probe_radius_nm)
            self.results.n_sphere_points = int(n_sphere_points)
            self.results.chunk_size = self._chunk_size

        def _single_frame(self) -> None:
            """Append current context coordinates and flush full chunks."""

            if self._zero_atom_selection:
                return
            self._position_chunk.append(
                np.asarray(self._context_atoms.positions, dtype=np.float32).copy()
            )
            self.results.max_buffered_coordinate_frames = max(
                int(self.results.max_buffered_coordinate_frames), len(self._position_chunk)
            )
            if len(self._position_chunk) >= self._chunk_size:
                self._flush_coordinate_chunk()

        def _conclude(self) -> None:
            """Compute raw SASA arrays and expose identity metadata."""

            frame_indices = _analysis_frame_indices(self)
            time_ns = _analysis_time_ns(self, frame_indices, self._raw_timestep_ps)
            if self._zero_atom_selection:
                _record_zero_selection_results(
                    self.results,
                    frame_indices=frame_indices,
                    time_ns=time_ns,
                    target_indices=self._target_indices,
                    context_indices=self._context_indices,
                )
                return

            self._flush_coordinate_chunk()
            atom_sasa_a2 = _concatenate_sasa_chunks(
                self._atom_sasa_chunks,
                n_frames=int(frame_indices.size),
                width=int(self._target_local_indices.size),
            )
            residue_sasa_a2 = _concatenate_sasa_chunks(
                self._residue_sasa_chunks,
                n_frames=int(frame_indices.size),
                width=len(self._residue_items),
            )
            total_sasa_a2 = _concatenate_sasa_vectors(
                self._total_sasa_chunks, n_frames=int(frame_indices.size)
            )
            if atom_sasa_a2.shape[0] != frame_indices.size:
                raise ValueError("SASA chunk output does not match selected frame count")
            residue_identity = _residue_identity(self._residue_items)
            self.results.atom_sasa_a2 = atom_sasa_a2
            self.results.residue_sasa_a2 = residue_sasa_a2
            self.results.total_sasa_a2 = total_sasa_a2
            self.results.frames = frame_indices
            self.results.time_ns = time_ns
            self.results.target_atom_indices = self._target_indices
            self.results.context_atom_indices = self._context_indices
            self.results.residue_keys = residue_identity["residue_keys"]
            self.results.residue_chainids = residue_identity["residue_chainids"]
            self.results.residue_resids = residue_identity["residue_resids"]
            self.results.residue_resnames = residue_identity["residue_resnames"]
            self.results.n_frames = int(frame_indices.size)
            self.results.n_target_atoms = int(self._target_indices.size)
            self.results.n_context_atoms = int(self._context_indices.size)
            self.results.n_target_residues = len(residue_identity["residue_keys"])
            self.results.zero_atom_selection = False

        def _flush_coordinate_chunk(self) -> None:
            """Compute SASA for the current bounded coordinate chunk."""

            if not self._position_chunk:
                return
            positions = np.asarray(self._position_chunk, dtype=np.float32)
            raw = _compute_sasa_chunk_from_positions(
                positions_angstrom=positions,
                topology=self._topology,
                target_local_indices=self._target_local_indices,
                residue_items=self._residue_items,
                probe_radius_nm=probe_radius_nm,
                n_sphere_points=n_sphere_points,
            )
            self._atom_sasa_chunks.append(raw["atom_sasa_a2"])
            self._residue_sasa_chunks.append(raw["residue_sasa_a2"])
            self._total_sasa_chunks.append(raw["total_sasa_a2"])
            self._position_chunk.clear()

    return SASASurfaceAreaAnalysis()


class SASAArtifactCollector:
    """Collect completed SASA jobs into a sidecar-backed replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Map completed SASA jobs to JSON summaries and NPZ sidecars.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context for one replicate.
        completed_jobs : sequence of MDAJobResult
            Completed SASA jobs.

        Returns
        -------
        ReplicateArtifact
            Canonical replicate artifact with raw arrays in sidecars only.
        """

        if not completed_jobs:
            raise ValueError("SASA collection expects at least one completed job")

        run_payloads: list[dict[str, Any]] = []
        sidecars = []
        sampled_frame_indices: set[int] = set()
        n_frames_total = int(ctx.frame_selection.n_frames_total or 0)
        for order, job in enumerate(completed_jobs):
            raw = _raw_result_from_job(job)
            sampled_frame_indices.update(int(frame) for frame in raw.frames.tolist())
            run_summary = _summarize_raw_sasa(
                raw=raw,
                run_label=str(job.results.run_label),
                target_selection=str(job.results.target_selection),
                context_selection=str(job.results.context_selection),
                n_frames_total=n_frames_total,
                raw_timestep_ps=ctx.frame_selection.timestep_ps,
                frame_stride=_effective_frame_stride(ctx.frame_selection, job.frame_selection),
            )
            sidecar_path = _sidecar_path(str(job.results.run_label), order)
            sidecar = ctx.artifact_store.write_npz_sidecar(
                sidecar_path,
                atom_sasa_a2=raw.atom_sasa_a2,
                residue_sasa_a2=raw.residue_sasa_a2,
                total_sasa_a2=raw.total_sasa_a2,
                frames=raw.frames,
                time_ns=raw.time_ns,
                target_atom_indices=raw.target_atom_indices,
                context_atom_indices=raw.context_atom_indices,
                residue_keys=np.asarray(raw.residue_keys, dtype="U64"),
                residue_chainids=np.asarray(raw.residue_chainids, dtype="U16"),
                residue_resids=np.asarray(raw.residue_resids, dtype=np.int64),
                residue_resnames=np.asarray(raw.residue_resnames, dtype="U16"),
                metadata={
                    "kind": "sasa_surface_area_arrays",
                    "job_name": job.name,
                    "run_label": str(job.results.run_label),
                    "target_selection": str(job.results.target_selection),
                    "context_selection": str(job.results.context_selection),
                    "units": "A^2",
                    "probe_radius_nm": float(job.results.probe_radius_nm),
                    "n_sphere_points": int(job.results.n_sphere_points),
                    "frame_selection": frame_selection_payload(job.frame_selection),
                    "shape": {
                        "atom_sasa_a2": list(raw.atom_sasa_a2.shape),
                        "residue_sasa_a2": list(raw.residue_sasa_a2.shape),
                        "total_sasa_a2": list(raw.total_sasa_a2.shape),
                    },
                },
            )
            run_summary.update(
                {
                    "sidecar": sidecar.model_dump(mode="json"),
                    "sidecar_path": sidecar.path,
                    "raw_npz_path": sidecar.path,
                    "npz_path": sidecar.path,
                    "metadata_path": None,
                    "raw_metadata_path": None,
                    "probe_radius_nm": float(job.results.probe_radius_nm),
                    "n_sphere_points": int(job.results.n_sphere_points),
                }
            )
            run_payloads.append(run_summary)
            sidecars.append(sidecar)

        eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
        config_hash = compute_config_hash(ctx.replicate_context.sim_config)
        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "run_results": run_payloads,
                "n_runs": len(run_payloads),
                "n_frames_total": n_frames_total,
                "n_frames_used": len(sampled_frame_indices),
                "trajectory_files": _trajectory_files(ctx.universe_policy.as_dict()),
                "metrics": _replicate_metrics(run_payloads),
                "replicate_metrics": _replicate_metrics(run_payloads),
                "metric_metadata": _metric_metadata(run_payloads),
            },
            sidecars=sidecars,
            provenance={
                "source": "mda_sasa_surface_area_jobs",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
            },
            metadata={
                "result_kind": "sasa_mda_replicate",
                "settings_fingerprint": ctx.settings_fingerprint,
                "config_hash": config_hash,
                "polyzymd_version": get_polyzymd_version(),
                "mdanalysis_version": mdanalysis_version(),
                "mdtraj_version": mdtraj_version(),
                "equilibration_time": eq_value,
                "equilibration_unit": eq_unit,
                "selection_string": "; ".join(
                    f"{run['target_selection']}|{run['context_selection']}" for run in run_payloads
                ),
                "statistical_policy": SASA_METRIC_METADATA["statistical_policy"],
            },
            warnings=_unique_warnings(ctx.warnings),
        )


def aggregate_sasa_artifacts(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: SASASettings,
    equilibration: str,
    output_dir: Path,
    result_path: Path,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> ConditionArtifact:
    """Aggregate SASA replicate artifacts into one condition artifact.

    Parameters
    ----------
    condition_label : str
        Condition label being aggregated.
    replicates : sequence of int
        Expected replicate IDs.
    settings : SASASettings
        Active SASA settings.
    equilibration : str
        Equilibration string from the framework context.
    output_dir : Path
        Aggregated output directory.
    result_path : Path
        Canonical condition artifact path.
    artifacts : sequence of ReplicateArtifact
        Per-replicate SASA artifacts.
    settings_fingerprint : str
        Active settings fingerprint.

    Returns
    -------
    ConditionArtifact
        Canonical condition artifact with per-run summaries.
    """

    ordered_artifacts = _validate_and_order_artifacts(
        condition_label=condition_label,
        expected_replicates=replicates,
        settings=settings,
        settings_fingerprint=settings_fingerprint,
        artifacts=artifacts,
    )
    analysis_dir = output_dir.parent
    run_results: list[dict[str, Any]] = []
    for run in settings.runs:
        entries = [_run_payload_for(artifact, run.label) for artifact in ordered_artifacts]
        _validate_structural_consistency(condition_label, run.label, entries)
        sidecar_rows = [
            _load_sasa_sidecar(artifact, run.label, analysis_dir / f"run_{artifact.replicate}")
            for artifact in ordered_artifacts
        ]
        _validate_residue_identity(condition_label, run.label, ordered_artifacts, sidecar_rows)
        run_results.append(
            _aggregate_run_payload(
                condition_label=condition_label,
                run_label=run.label,
                entries=entries,
                sidecar_rows=sidecar_rows,
                replicates=replicates,
            )
        )

    eq_value, eq_unit = parse_time_string(equilibration)
    config_hash = str(ordered_artifacts[0].metadata.get("config_hash", "unknown"))
    artifact = ConditionArtifact(
        analysis_name="sasa",
        condition_label=condition_label,
        replicates=[int(replicate) for replicate in replicates],
        payload={
            "run_results": run_results,
            "n_runs": len(run_results),
            "n_replicates": len(ordered_artifacts),
            "settings_fingerprint": settings_fingerprint,
            "metrics": _condition_metrics(run_results),
            "replicate_metrics": _condition_replicate_metrics(run_results, replicates),
            "metric_metadata": _metric_metadata(run_results),
        },
        provenance={
            "source": "sasa_replicate_artifact_aggregation",
            "frame_selection": ordered_artifacts[0].provenance.get("frame_selection"),
            "statistical_policy": SASA_METRIC_METADATA["statistical_policy"],
        },
        metadata={
            "result_kind": "sasa_mda_condition",
            "settings_fingerprint": settings_fingerprint,
            "config_hash": config_hash,
            "polyzymd_version": get_polyzymd_version(),
            "mdanalysis_version": mdanalysis_version(),
            "mdtraj_version": mdtraj_version(),
            "equilibration_time": eq_value,
            "equilibration_unit": eq_unit,
            "selection_string": "; ".join(
                f"{run.target_selection}|{run.context_selection or run.target_selection}"
                for run in settings.runs
            ),
            "source_result_files": [
                str(output_dir.parent / f"run_{artifact.replicate}" / "result.json")
                for artifact in ordered_artifacts
            ],
            "n_replicates": len(ordered_artifacts),
            "statistical_policy": SASA_METRIC_METADATA["statistical_policy"],
        },
        source_replicates=[
            {
                "replicate": int(artifact.replicate),
                "path": str(output_dir.parent / f"run_{artifact.replicate}" / "result.json"),
            }
            for artifact in ordered_artifacts
        ],
        warnings=_combined_warnings(ordered_artifacts),
    )
    ArtifactStore(result_path.parent).write_condition_result(artifact, result_path.name)
    return artifact


def condition_artifact_to_legacy_result(artifact: ConditionArtifact) -> Any:
    """Adapt a canonical condition artifact to the legacy aggregate model.

    Parameters
    ----------
    artifact : ConditionArtifact
        Canonical SASA condition artifact.

    Returns
    -------
    SASAAggregatedResult
        Legacy-shaped aggregate model used by existing comparison and formatters.
    """

    from polyzymd.analyses.sasa._results import SASAAggregatedResult, SASARunAggregatedResult

    payload = artifact.payload
    return SASAAggregatedResult(
        config_hash=str(artifact.metadata.get("config_hash", "unknown")),
        polyzymd_version=str(artifact.metadata.get("polyzymd_version", get_polyzymd_version())),
        replicate=None,
        equilibration_time=float(artifact.metadata.get("equilibration_time", 0.0)),
        equilibration_unit=str(artifact.metadata.get("equilibration_unit", "ns")),
        selection_string=str(artifact.metadata.get("selection_string", "")),
        replicates=[int(replicate) for replicate in artifact.replicates],
        n_replicates=int(payload.get("n_replicates", len(artifact.replicates))),
        run_results=[
            SASARunAggregatedResult(
                config_hash=str(artifact.metadata.get("config_hash", "unknown")),
                polyzymd_version=str(
                    artifact.metadata.get("polyzymd_version", get_polyzymd_version())
                ),
                replicate=None,
                equilibration_time=float(artifact.metadata.get("equilibration_time", 0.0)),
                equilibration_unit=str(artifact.metadata.get("equilibration_unit", "ns")),
                selection_string=str(run_result.get("target_selection", "")),
                replicates=[int(replicate) for replicate in run_result.get("replicates", [])],
                n_replicates=int(run_result.get("n_replicates", 0)),
                run_label=str(run_result["run_label"]),
                target_selection=str(run_result["target_selection"]),
                context_selection=str(run_result["context_selection"]),
                overall_mean=_legacy_float(run_result.get("overall_mean")),
                overall_sem=_legacy_float(run_result.get("overall_sem")),
                overall_median=_legacy_float(run_result.get("overall_median")),
                overall_min=_legacy_float(run_result.get("overall_min")),
                overall_max=_legacy_float(run_result.get("overall_max")),
                overall_final=_legacy_float(run_result.get("overall_final")),
                per_replicate_means=_legacy_float_list(run_result["per_replicate_means"]),
                per_replicate_stds=_legacy_float_list(run_result["per_replicate_stds"]),
                per_replicate_medians=_legacy_float_list(run_result["per_replicate_medians"]),
                per_replicate_mins=_legacy_float_list(run_result["per_replicate_mins"]),
                per_replicate_maxs=_legacy_float_list(run_result["per_replicate_maxs"]),
                per_replicate_finals=_legacy_float_list(run_result["per_replicate_finals"]),
                n_target_atoms=int(run_result["n_target_atoms"]),
                n_context_atoms=run_result.get("n_context_atoms"),
                per_replicate_context_atom_counts=list(
                    run_result.get("per_replicate_context_atom_counts", [])
                ),
                n_context_atoms_variable=bool(run_result.get("n_context_atoms_variable", False)),
                n_target_residues=int(run_result["n_target_residues"]),
                zero_atom_selection=bool(run_result.get("zero_atom_selection", False)),
                residue_keys=list(run_result.get("residue_keys", [])),
                residue_chainids=list(run_result.get("residue_chainids", [])),
                residue_resids=[int(value) for value in run_result.get("residue_resids", [])],
                residue_resnames=list(run_result.get("residue_resnames", [])),
                per_residue_mean_sasa=list(run_result.get("per_residue_mean_sasa", [])),
                per_residue_sem_sasa=list(run_result.get("per_residue_sem_sasa", [])),
            )
            for run_result in payload.get("run_results", [])
        ],
        settings_fingerprint=artifact.metadata.get("settings_fingerprint"),
        source_result_files=list(artifact.metadata.get("source_result_files", [])),
    )


def load_condition_artifact(aggregated_dir: Path) -> ConditionArtifact | None:
    """Load a canonical SASA condition artifact if it exists."""

    result_path = aggregated_dir / "result.json"
    if not result_path.exists():
        return None
    return ArtifactStore(aggregated_dir).read_condition_result("result.json")


def load_replicate_sasa_sidecar(
    artifact: ReplicateArtifact,
    run_label: str,
    run_dir: Path,
) -> SASAComputationResult:
    """Load and validate the SASA NPZ sidecar for one run.

    Parameters
    ----------
    artifact : ReplicateArtifact
        Replicate artifact that references the sidecar.
    run_label : str
        SASA run label to load.
    run_dir : Path
        Replicate artifact-store root.

    Returns
    -------
    SASAComputationResult
        Raw SASA arrays and identity metadata.
    """

    sidecar = _sasa_sidecar(artifact, run_label)
    sidecar_path = ArtifactStore(run_dir).validate_sidecar(sidecar)
    try:
        with np.load(sidecar_path) as data:
            raw = SASAComputationResult(
                atom_sasa_a2=np.asarray(data["atom_sasa_a2"], dtype=np.float64),
                residue_sasa_a2=np.asarray(data["residue_sasa_a2"], dtype=np.float64),
                total_sasa_a2=np.asarray(data["total_sasa_a2"], dtype=np.float64),
                frames=np.asarray(data["frames"], dtype=np.int64),
                time_ns=np.asarray(data["time_ns"], dtype=np.float64),
                target_atom_indices=np.asarray(data["target_atom_indices"], dtype=np.int64),
                context_atom_indices=np.asarray(data["context_atom_indices"], dtype=np.int64),
                residue_keys=[str(value) for value in data["residue_keys"].tolist()],
                residue_chainids=[str(value) for value in data["residue_chainids"].tolist()],
                residue_resids=[int(value) for value in data["residue_resids"].tolist()],
                residue_resnames=[str(value) for value in data["residue_resnames"].tolist()],
            )
    except (OSError, KeyError, ValueError) as exc:
        raise ValueError(
            f"Failed to load SASA sidecar for replicate {artifact.replicate} run "
            f"'{run_label}': {exc}"
        ) from exc
    _validate_raw_sasa_shapes(raw, run_label=run_label, replicate=artifact.replicate)
    _validate_sidecar_matches_payload(artifact, run_label, raw)
    return raw


def _compute_sasa_from_positions(
    *,
    positions_angstrom: NDArray[np.float32],
    target_atoms: Any,
    context_atoms: Any,
    target_indices: NDArray[np.int64],
    context_indices: NDArray[np.int64],
    probe_radius_nm: float,
    n_sphere_points: int,
    chunk_size: int,
) -> SASAComputationResult:
    """Run MDTraj Shrake-Rupley on buffered AnalysisBase coordinates."""

    import mdtraj as md

    if chunk_size <= 0:
        raise ValueError("chunk_size must be >= 1")
    topology = _mdtraj_topology_from_atom_group(context_atoms)
    target_local_indices = _target_local_indices(target_indices, context_indices)
    residue_items = _residue_items(target_atoms, target_local_indices)
    n_frames = int(positions_angstrom.shape[0])
    atom_sasa_target_a2 = np.empty((n_frames, target_local_indices.size), dtype=np.float64)
    residue_sasa_a2 = np.empty((n_frames, len(residue_items)), dtype=np.float64)
    total_sasa_a2 = np.empty(n_frames, dtype=np.float64)
    positions_nm = positions_angstrom.astype(np.float32, copy=False) / 10.0

    for chunk_start in range(0, n_frames, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n_frames)
        mdtraj_traj = md.Trajectory(xyz=positions_nm[chunk_start:chunk_end], topology=topology)
        atom_sasa_nm2 = np.asarray(
            md.shrake_rupley(
                mdtraj_traj,
                mode="atom",
                probe_radius=probe_radius_nm,
                n_sphere_points=n_sphere_points,
            ),
            dtype=np.float64,
        )
        atom_sasa_target_chunk = atom_sasa_nm2[:, target_local_indices] * NM2_TO_A2
        atom_sasa_target_a2[chunk_start:chunk_end, :] = atom_sasa_target_chunk
        for idx, (_, atom_locals) in enumerate(residue_items):
            residue_sasa_a2[chunk_start:chunk_end, idx] = (
                np.sum(atom_sasa_nm2[:, atom_locals], axis=1) * NM2_TO_A2
            )
        total_sasa_a2[chunk_start:chunk_end] = np.sum(atom_sasa_target_chunk, axis=1)

    residue_identity = _residue_identity(residue_items)
    return SASAComputationResult(
        atom_sasa_a2=atom_sasa_target_a2,
        residue_sasa_a2=residue_sasa_a2,
        total_sasa_a2=total_sasa_a2,
        frames=np.empty((0,), dtype=np.int64),
        time_ns=np.empty((0,), dtype=np.float64),
        target_atom_indices=target_indices,
        context_atom_indices=context_indices,
        residue_keys=residue_identity["residue_keys"],
        residue_chainids=residue_identity["residue_chainids"],
        residue_resids=residue_identity["residue_resids"],
        residue_resnames=residue_identity["residue_resnames"],
    )


def _compute_sasa_chunk_from_positions(
    *,
    positions_angstrom: NDArray[np.float32],
    topology: Any,
    target_local_indices: NDArray[np.int64],
    residue_items: Sequence[Any],
    probe_radius_nm: float,
    n_sphere_points: int,
) -> dict[str, NDArray[np.float64]]:
    """Run MDTraj Shrake-Rupley for one bounded coordinate chunk.

    Parameters
    ----------
    positions_angstrom : ndarray of float32
        Context coordinates for at most ``chunk_size`` selected frames.
    topology : Any
        MDTraj topology for the context atom group.
    target_local_indices : ndarray of int64
        Context-local indices of target atoms.
    residue_items : sequence
        Residue identity keys paired with context-local atom indices.
    probe_radius_nm : float
        Shrake-Rupley probe radius in nm.
    n_sphere_points : int
        Number of Shrake-Rupley sphere points.

    Returns
    -------
    dict[str, ndarray]
        Atom, residue, and total SASA arrays for the input chunk.
    """

    import mdtraj as md

    n_frames = int(positions_angstrom.shape[0])
    atom_sasa_target_a2 = np.empty((n_frames, target_local_indices.size), dtype=np.float64)
    residue_sasa_a2 = np.empty((n_frames, len(residue_items)), dtype=np.float64)
    total_sasa_a2 = np.empty(n_frames, dtype=np.float64)
    positions_nm = positions_angstrom.astype(np.float32, copy=False) / 10.0
    mdtraj_traj = md.Trajectory(xyz=positions_nm, topology=topology)
    atom_sasa_nm2 = np.asarray(
        md.shrake_rupley(
            mdtraj_traj,
            mode="atom",
            probe_radius=probe_radius_nm,
            n_sphere_points=n_sphere_points,
        ),
        dtype=np.float64,
    )
    atom_sasa_target_chunk = atom_sasa_nm2[:, target_local_indices] * NM2_TO_A2
    atom_sasa_target_a2[:, :] = atom_sasa_target_chunk
    for idx, (_, atom_locals) in enumerate(residue_items):
        residue_sasa_a2[:, idx] = np.sum(atom_sasa_nm2[:, atom_locals], axis=1) * NM2_TO_A2
    total_sasa_a2[:] = np.sum(atom_sasa_target_chunk, axis=1)
    return {
        "atom_sasa_a2": atom_sasa_target_a2,
        "residue_sasa_a2": residue_sasa_a2,
        "total_sasa_a2": total_sasa_a2,
    }


def _concatenate_sasa_chunks(
    chunks: Sequence[NDArray[np.float64]], *, n_frames: int, width: int
) -> NDArray[np.float64]:
    """Concatenate two-dimensional SASA chunks with an empty-frame fallback."""

    if chunks:
        return np.concatenate(chunks, axis=0).astype(np.float64, copy=False)
    return np.empty((n_frames, width), dtype=np.float64)


def _concatenate_sasa_vectors(
    chunks: Sequence[NDArray[np.float64]], *, n_frames: int
) -> NDArray[np.float64]:
    """Concatenate one-dimensional SASA chunks with an empty-frame fallback."""

    if chunks:
        return np.concatenate(chunks, axis=0).astype(np.float64, copy=False)
    return np.empty(n_frames, dtype=np.float64)


def _mdtraj_topology_from_atom_group(atom_group: Any) -> Any:
    """Build an MDTraj topology from an MDAnalysis atom group."""

    import mdtraj as md

    if hasattr(atom_group, "write"):
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=True) as tmp_pdb:
            atom_group.write(tmp_pdb.name)
            return md.load(tmp_pdb.name).topology

    topology = md.Topology()
    chain = topology.add_chain()
    current_residue_key: tuple[int, int, str] | None = None
    md_residue = None
    for atom in atom_group:
        residue = getattr(atom, "residue", None)
        resid = int(getattr(atom, "resid", getattr(residue, "resid", 1)))
        resname = str(getattr(atom, "resname", getattr(residue, "resname", "UNK")))
        resindex = int(getattr(atom, "resindex", getattr(residue, "ix", resid)))
        residue_key = (resindex, resid, resname)
        if residue_key != current_residue_key:
            md_residue = topology.add_residue(resname, chain, resSeq=resid)
            current_residue_key = residue_key
        topology.add_atom(str(getattr(atom, "name", "C")), _mdtraj_element(md, atom), md_residue)
    return topology


def _mdtraj_element(md: Any, atom: Any) -> Any:
    """Return a best-effort MDTraj element for an atom-like object."""

    symbol = str(getattr(atom, "element", "") or "").strip()
    if not symbol:
        atom_name = str(getattr(atom, "name", "")).strip()
        symbol = "".join(char for char in atom_name if char.isalpha())[:1]
    symbol = symbol.capitalize() if symbol else "C"
    try:
        return md.element.get_by_symbol(symbol)
    except KeyError:
        return md.element.carbon


def _record_zero_selection_results(
    results: Any,
    *,
    frame_indices: NDArray[np.int64],
    time_ns: NDArray[np.float64],
    target_indices: NDArray[np.int64],
    context_indices: NDArray[np.int64],
) -> None:
    """Populate raw result arrays for a target or context zero-selection run."""

    if target_indices.size == 0:
        LOGGER.warning(
            "Run '%s' target selection matched zero atoms (%r); returning missing SASA metrics",
            results.run_label,
            results.target_selection,
        )
    if context_indices.size == 0:
        LOGGER.warning(
            "Run '%s' context selection matched zero atoms (%r); returning missing SASA metrics",
            results.run_label,
            results.context_selection,
        )
    n_frames = int(frame_indices.size)
    results.atom_sasa_a2 = np.empty((n_frames, target_indices.size), dtype=np.float64)
    results.residue_sasa_a2 = np.empty((n_frames, 0), dtype=np.float64)
    results.total_sasa_a2 = np.zeros(n_frames, dtype=np.float64)
    results.frames = frame_indices
    results.time_ns = time_ns
    results.target_atom_indices = target_indices
    results.context_atom_indices = context_indices
    results.residue_keys = []
    results.residue_chainids = []
    results.residue_resids = []
    results.residue_resnames = []
    results.n_frames = n_frames
    results.n_target_atoms = int(target_indices.size)
    results.n_context_atoms = int(context_indices.size)
    results.n_target_residues = 0
    results.zero_atom_selection = True


def _raw_result_from_job(job: MDAJobResult) -> SASAComputationResult:
    """Extract and validate raw SASA arrays from one completed job."""

    raw = SASAComputationResult(
        atom_sasa_a2=np.asarray(job.results.atom_sasa_a2, dtype=np.float64),
        residue_sasa_a2=np.asarray(job.results.residue_sasa_a2, dtype=np.float64),
        total_sasa_a2=np.asarray(job.results.total_sasa_a2, dtype=np.float64),
        frames=np.asarray(job.results.frames, dtype=np.int64),
        time_ns=np.asarray(job.results.time_ns, dtype=np.float64),
        target_atom_indices=np.asarray(job.results.target_atom_indices, dtype=np.int64),
        context_atom_indices=np.asarray(job.results.context_atom_indices, dtype=np.int64),
        residue_keys=[str(value) for value in job.results.residue_keys],
        residue_chainids=[str(value) for value in job.results.residue_chainids],
        residue_resids=[int(value) for value in job.results.residue_resids],
        residue_resnames=[str(value) for value in job.results.residue_resnames],
    )
    _validate_raw_sasa_shapes(raw, run_label=str(job.results.run_label), replicate=None)
    return raw


def _summarize_raw_sasa(
    *,
    raw: SASAComputationResult,
    run_label: str,
    target_selection: str,
    context_selection: str,
    n_frames_total: int,
    raw_timestep_ps: float | None,
    frame_stride: int | None,
) -> dict[str, Any]:
    """Summarize raw SASA arrays into JSON-safe replicate metadata."""

    from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

    zero_atom = raw.target_atom_indices.size == 0 or raw.context_atom_indices.size == 0
    total = np.asarray(raw.total_sasa_a2, dtype=np.float64)
    finite_total = total[np.isfinite(total)] if not zero_atom else np.asarray([], dtype=np.float64)
    missing_sasa_reason: str | None = None
    if zero_atom:
        missing_sasa_reason = "zero_atom_selection"
        mean_sasa = std_sasa = median_sasa = min_sasa = max_sasa = final_sasa = None
    elif finite_total.size:
        mean_sasa = float(np.mean(finite_total))
        std_sasa = float(np.std(finite_total, ddof=0))
        median_sasa = float(np.median(finite_total))
        min_sasa = float(np.min(finite_total))
        max_sasa = float(np.max(finite_total))
        final_sasa = _json_float_or_none(total[-1])
    else:
        missing_sasa_reason = "no_finite_sasa_samples"
        mean_sasa = std_sasa = median_sasa = min_sasa = max_sasa = final_sasa = None

    sem_sasa: float | None = None
    correlation_time: float | None = None
    correlation_time_unit: str | None = None
    n_independent_frames: int | None = None
    statistical_inefficiency: float | None = None
    autocorrelation_warning: str | None = None
    effective_timestep_ps = _effective_timestep_ps(raw_timestep_ps, frame_stride)
    if finite_total.size >= 20 and std_sasa is not None:
        try:
            tau = estimate_correlation_time(
                finite_total,
                timestep=effective_timestep_ps,
                timestep_unit="ps",
                method="integration",
                n_frames=len(finite_total),
            )
            correlation_time = tau.tau
            correlation_time_unit = tau.tau_unit
            n_independent_frames = tau.n_independent
            statistical_inefficiency = tau.statistical_inefficiency
            autocorrelation_warning = tau.warning
            if n_independent_frames and n_independent_frames > 0 and np.isfinite(std_sasa):
                sem_sasa = float(std_sasa / np.sqrt(float(n_independent_frames)))
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            LOGGER.warning("Autocorrelation analysis failed for SASA run '%s': %s", run_label, exc)
            if np.isfinite(std_sasa):
                sem_sasa = float(std_sasa / np.sqrt(float(finite_total.size)))

    return {
        "run_label": run_label,
        "target_selection": target_selection,
        "context_selection": context_selection,
        "mean_sasa": mean_sasa,
        "std_sasa": std_sasa,
        "median_sasa": median_sasa,
        "min_sasa": min_sasa,
        "max_sasa": max_sasa,
        "final_sasa": final_sasa,
        "sem_sasa": sem_sasa,
        "missing_sasa_reason": missing_sasa_reason,
        "correlation_time": correlation_time,
        "correlation_time_unit": correlation_time_unit,
        "n_independent_frames": n_independent_frames,
        "statistical_inefficiency": statistical_inefficiency,
        "autocorrelation_warning": autocorrelation_warning,
        "n_frames_total": n_frames_total,
        "n_frames_used": int(raw.frames.size),
        "n_target_atoms": int(raw.target_atom_indices.size),
        "n_context_atoms": int(raw.context_atom_indices.size),
        "n_target_residues": len(raw.residue_keys),
        "zero_atom_selection": zero_atom,
        "time_unit": "ns",
        "timestep_ps": effective_timestep_ps,
        "raw_timestep_ps": raw_timestep_ps,
        "frame_stride": frame_stride,
        "frames_first": int(raw.frames[0]) if raw.frames.size else None,
        "frames_last": int(raw.frames[-1]) if raw.frames.size else None,
    }


def _aggregate_run_payload(
    *,
    condition_label: str,
    run_label: str,
    entries: Sequence[Mapping[str, Any]],
    sidecar_rows: Sequence[SASAComputationResult],
    replicates: Sequence[int],
) -> dict[str, Any]:
    """Aggregate one SASA run across replicate entries and sidecars."""

    del condition_label
    per_means = [_json_float_or_none(entry.get("mean_sasa")) for entry in entries]
    per_stds = [_json_float_or_none(entry.get("std_sasa")) for entry in entries]
    per_medians = [_json_float_or_none(entry.get("median_sasa")) for entry in entries]
    per_mins = [_json_float_or_none(entry.get("min_sasa")) for entry in entries]
    per_maxs = [_json_float_or_none(entry.get("max_sasa")) for entry in entries]
    per_finals = [_json_float_or_none(entry.get("final_sasa")) for entry in entries]
    finite_means = np.asarray([value for value in per_means if value is not None], dtype=np.float64)
    if finite_means.size:
        mean_stats = compute_sem(finite_means)
        overall_mean = _json_float_or_none(mean_stats.mean)
        overall_sem = _json_float_or_none(mean_stats.sem)
    else:
        overall_mean = None
        overall_sem = None

    residue_matrix: list[NDArray[np.float64]] = []
    residue_keys: list[str] = []
    residue_chainids: list[str] = []
    residue_resids: list[int] = []
    residue_resnames: list[str] = []
    for raw in sidecar_rows:
        if raw.residue_sasa_a2.size == 0:
            continue
        residue_matrix.append(np.nanmean(raw.residue_sasa_a2, axis=0).astype(np.float64))
        if not residue_keys:
            residue_keys = list(raw.residue_keys)
            residue_chainids = list(raw.residue_chainids)
            residue_resids = list(raw.residue_resids)
            residue_resnames = list(raw.residue_resnames)

    if residue_matrix:
        stacked = np.stack(residue_matrix, axis=0)
        per_residue_mean = np.nanmean(stacked, axis=0).astype(np.float64)
        if stacked.shape[0] > 1:
            per_residue_sem = np.nanstd(stacked, axis=0, ddof=1) / math.sqrt(stacked.shape[0])
        else:
            per_residue_sem = np.zeros(stacked.shape[1], dtype=np.float64)
    else:
        per_residue_mean = np.asarray([], dtype=np.float64)
        per_residue_sem = np.asarray([], dtype=np.float64)

    template = entries[0]
    context_atom_counts = [int(entry["n_context_atoms"]) for entry in entries]
    context_count_is_variable = len(set(context_atom_counts)) > 1
    return {
        "run_label": run_label,
        "target_selection": str(template["target_selection"]),
        "context_selection": str(template["context_selection"]),
        "replicates": [int(replicate) for replicate in replicates],
        "n_replicates": len(replicates),
        "overall_mean": overall_mean,
        "overall_sem": overall_sem,
        "overall_median": _mean_optional_values(per_medians),
        "overall_min": _min_optional_values(per_mins),
        "overall_max": _max_optional_values(per_maxs),
        "overall_final": _mean_optional_values(per_finals),
        "per_replicate_means": per_means,
        "per_replicate_stds": per_stds,
        "per_replicate_medians": per_medians,
        "per_replicate_mins": per_mins,
        "per_replicate_maxs": per_maxs,
        "per_replicate_finals": per_finals,
        "n_target_atoms": int(template["n_target_atoms"]),
        "n_context_atoms": None if context_count_is_variable else context_atom_counts[0],
        "per_replicate_context_atom_counts": context_atom_counts,
        "n_context_atoms_variable": context_count_is_variable,
        "n_target_residues": int(template["n_target_residues"]),
        "zero_atom_selection": all(bool(entry["zero_atom_selection"]) for entry in entries),
        "residue_keys": residue_keys,
        "residue_chainids": residue_chainids,
        "residue_resids": residue_resids,
        "residue_resnames": residue_resnames,
        "per_residue_mean_sasa": per_residue_mean.tolist(),
        "per_residue_sem_sasa": per_residue_sem.tolist(),
    }


def _validate_and_order_artifacts(
    *,
    condition_label: str,
    expected_replicates: Sequence[int],
    settings: SASASettings,
    settings_fingerprint: str,
    artifacts: Sequence[ReplicateArtifact],
) -> list[ReplicateArtifact]:
    """Validate replicate artifacts and return them in requested order."""

    by_replicate: dict[int, ReplicateArtifact] = {}
    expected_run_labels = [run.label for run in settings.runs]
    for artifact in artifacts:
        if artifact.analysis_name != "sasa":
            raise ValueError(
                f"SASA aggregation received an artifact for {artifact.analysis_name!r}"
            )
        if artifact.condition_label != condition_label:
            raise ValueError(
                f"SASA artifact condition mismatch: expected {condition_label!r}, "
                f"got {artifact.condition_label!r}"
            )
        stored_fingerprint = artifact.metadata.get("settings_fingerprint")
        if stored_fingerprint != settings_fingerprint:
            raise ValueError(
                "SASA settings fingerprint mismatch for replicate "
                f"{artifact.replicate}: expected {settings_fingerprint!r}, got "
                f"{stored_fingerprint!r}. Recompute stale caches."
            )
        observed = [
            str(entry.get("run_label")) for entry in artifact.payload.get("run_results", [])
        ]
        missing = [label for label in expected_run_labels if label not in observed]
        duplicates = sorted({label for label in observed if observed.count(label) > 1})
        if missing or duplicates:
            raise ValueError(
                f"SASA replicate {artifact.replicate} has incomplete run coverage. "
                f"Missing={missing}, duplicates={duplicates}."
            )
        if artifact.replicate in by_replicate:
            raise ValueError(f"Duplicate SASA replicate artifact {artifact.replicate}")
        by_replicate[int(artifact.replicate)] = artifact

    ordered: list[ReplicateArtifact] = []
    missing_replicates: list[int] = []
    for replicate in expected_replicates:
        artifact = by_replicate.get(int(replicate))
        if artifact is None:
            missing_replicates.append(int(replicate))
            continue
        ordered.append(artifact)
    if missing_replicates:
        raise ValueError(f"Missing SASA replicate artifacts: {missing_replicates}")
    return ordered


def _load_sasa_sidecar(
    artifact: ReplicateArtifact,
    run_label: str,
    run_dir: Path,
) -> SASAComputationResult:
    """Load the sidecar for aggregation with contextual errors."""

    try:
        return load_replicate_sasa_sidecar(artifact, run_label, run_dir)
    except Exception as exc:
        raise ValueError(
            f"SASA aggregation failed to validate sidecar for replicate {artifact.replicate} "
            f"run '{run_label}': {exc}"
        ) from exc


def _run_payload_for(artifact: ReplicateArtifact, run_label: str) -> Mapping[str, Any]:
    """Return the run payload for ``run_label`` from an artifact."""

    for entry in artifact.payload.get("run_results", []):
        if isinstance(entry, Mapping) and entry.get("run_label") == run_label:
            return entry
    raise ValueError(f"SASA replicate {artifact.replicate} is missing run {run_label!r}")


def _validate_structural_consistency(
    condition_label: str,
    run_label: str,
    entries: Sequence[Mapping[str, Any]],
) -> None:
    """Validate target atom and residue counts across entries."""

    first = entries[0]
    expected = (int(first["n_target_atoms"]), int(first["n_target_residues"]))
    for entry in entries[1:]:
        current = (int(entry["n_target_atoms"]), int(entry["n_target_residues"]))
        if current != expected:
            raise ValueError(
                f"SASA aggregation for condition '{condition_label}' run '{run_label}' found "
                f"target metadata mismatch: {current} != {expected}."
            )


def _validate_residue_identity(
    condition_label: str,
    run_label: str,
    artifacts: Sequence[ReplicateArtifact],
    sidecar_rows: Sequence[SASAComputationResult],
) -> None:
    """Validate residue identity and order across non-empty sidecars."""

    reference: SASAComputationResult | None = None
    reference_replicate: int | None = None
    for artifact, raw in zip(artifacts, sidecar_rows, strict=True):
        if raw.residue_sasa_a2.size == 0:
            continue
        if reference is None:
            reference = raw
            reference_replicate = artifact.replicate
            continue
        if (
            raw.residue_keys != reference.residue_keys
            or raw.residue_chainids != reference.residue_chainids
            or raw.residue_resids != reference.residue_resids
            or raw.residue_resnames != reference.residue_resnames
        ):
            raise ValueError(
                f"SASA aggregation for condition '{condition_label}' run '{run_label}' found "
                "residue metadata mismatch between replicates "
                f"{reference_replicate} and {artifact.replicate}."
            )


def _validate_raw_sasa_shapes(
    raw: SASAComputationResult,
    *,
    run_label: str,
    replicate: int | None,
) -> None:
    """Validate raw SASA array dimensionality and identity lengths."""

    context = f"SASA run '{run_label}'"
    if replicate is not None:
        context += f" replicate {replicate}"
    n_frames = int(raw.total_sasa_a2.shape[0])
    if raw.total_sasa_a2.ndim != 1 or raw.frames.ndim != 1 or raw.time_ns.ndim != 1:
        raise ValueError(f"{context} has invalid one-dimensional arrays")
    if raw.frames.size != n_frames or raw.time_ns.size != n_frames:
        raise ValueError(f"{context} frame/time arrays do not match total SASA length")
    if raw.atom_sasa_a2.ndim != 2 or raw.atom_sasa_a2.shape[0] != n_frames:
        raise ValueError(f"{context} atom SASA array has invalid shape")
    if raw.residue_sasa_a2.ndim != 2 or raw.residue_sasa_a2.shape[0] != n_frames:
        raise ValueError(f"{context} residue SASA array has invalid shape")
    if raw.atom_sasa_a2.shape[1] != raw.target_atom_indices.size:
        raise ValueError(f"{context} target atom identity does not match atom SASA width")
    if not (
        len(raw.residue_keys)
        == len(raw.residue_chainids)
        == len(raw.residue_resids)
        == len(raw.residue_resnames)
        == raw.residue_sasa_a2.shape[1]
    ):
        raise ValueError(f"{context} residue identity does not match residue SASA width")


def _validate_sidecar_matches_payload(
    artifact: ReplicateArtifact,
    run_label: str,
    raw: SASAComputationResult,
) -> None:
    """Validate sidecar shapes against replicate JSON summary metadata."""

    entry = _run_payload_for(artifact, run_label)
    checks = {
        "n_target_atoms": int(raw.target_atom_indices.size),
        "n_context_atoms": int(raw.context_atom_indices.size),
        "n_target_residues": len(raw.residue_keys),
        "n_frames_used": int(raw.frames.size),
    }
    for key, value in checks.items():
        if int(entry.get(key, -1)) != value:
            raise ValueError(
                f"SASA sidecar identity mismatch for replicate {artifact.replicate} "
                f"run '{run_label}': {key} differs from result.json metadata."
            )


def _sasa_sidecar(artifact: ReplicateArtifact, run_label: str) -> ArtifactSidecarRef:
    """Return the SASA sidecar reference for ``run_label``."""

    entry = _run_payload_for(artifact, run_label)
    sidecar_path = entry.get("sidecar_path") or entry.get("raw_npz_path") or entry.get("npz_path")
    for sidecar in artifact.sidecars:
        if sidecar.metadata.get("run_label") == run_label or sidecar.path == sidecar_path:
            return sidecar
    raise ValueError(
        f"SASA replicate {artifact.replicate} is missing sidecar for run {run_label!r}"
    )


def _analysis_frame_indices(analysis: Any) -> NDArray[np.int64]:
    """Return frame indices selected by AnalysisBase."""

    frames = getattr(analysis, "frames", None)
    if frames is None:
        n_frames = len(getattr(analysis, "times", []))
        return np.arange(n_frames, dtype=np.int64)
    return np.asarray(frames, dtype=np.int64)


def _analysis_time_ns(
    analysis: Any,
    frame_indices: NDArray[np.int64],
    raw_timestep_ps: float | None,
) -> NDArray[np.float64]:
    """Return selected frame times in ns."""

    times = getattr(analysis, "times", None)
    if times is not None:
        time_ps = np.asarray(times, dtype=np.float64)
        if time_ps.shape[0] == frame_indices.shape[0]:
            return time_ps / 1000.0
    timestep = (
        raw_timestep_ps if raw_timestep_ps is not None and math.isfinite(raw_timestep_ps) else 1.0
    )
    return frame_indices.astype(np.float64) * float(timestep) / 1000.0


def _frame_selection_with_stride(frame_selection: FrameSelection, stride: int) -> FrameSelection:
    """Return a frame selection with SASA run-specific stride applied."""

    if stride <= 0:
        raise ValueError("SASA stride must be >= 1")
    if frame_selection.frames is not None:
        frames = tuple(frame_selection.frames)
        if all(_is_boolean_frame_value(frame) for frame in frames):
            selected = tuple(index for index, value in enumerate(frames) if value)[::stride]
        else:
            selected = tuple(int(frame) for frame in frames)[::stride]
        return FrameSelection(
            frames=selected,
            equilibration=frame_selection.equilibration,
            equilibration_start=frame_selection.equilibration_start,
            equilibration_ps=frame_selection.equilibration_ps,
            timestep_ps=frame_selection.timestep_ps,
            n_frames_total=frame_selection.n_frames_total,
            warning_message=frame_selection.warning_message,
        )

    return FrameSelection(
        start=frame_selection.start,
        stop=frame_selection.stop,
        step=(frame_selection.step or 1) * stride,
        equilibration=frame_selection.equilibration,
        equilibration_start=frame_selection.equilibration_start,
        equilibration_ps=frame_selection.equilibration_ps,
        timestep_ps=frame_selection.timestep_ps,
        n_frames_total=frame_selection.n_frames_total,
        warning_message=frame_selection.warning_message,
    )


def _is_boolean_frame_value(frame: Any) -> bool:
    """Return whether a frame selector value is boolean-like."""

    if isinstance(frame, bool):
        return True
    frame_type = type(frame)
    return frame_type.__name__ == "bool_" and frame_type.__module__.split(".", 1)[0] == "numpy"


def _target_local_indices(
    target_indices: NDArray[np.int64], context_indices: NDArray[np.int64]
) -> NDArray[np.int64]:
    """Map universe-global target atom indices to context-local indices."""

    context_index_to_local = {int(idx): i for i, idx in enumerate(context_indices.tolist())}
    return np.asarray(
        [context_index_to_local[int(idx)] for idx in target_indices.tolist()],
        dtype=np.int64,
    )


def _residue_items(target_atoms: Any, target_local_indices: NDArray[np.int64]) -> list[Any]:
    """Return residue-to-target-local atom index groups."""

    residue_to_indices: dict[tuple[str, int, str], list[int]] = {}
    for atom_local, atom in zip(target_local_indices.tolist(), target_atoms, strict=True):
        chain_id = str(getattr(atom, "chainID", getattr(atom, "chainid", "")))
        key = (chain_id, int(atom.resid), str(atom.resname))
        residue_to_indices.setdefault(key, []).append(int(atom_local))
    return list(residue_to_indices.items())


def _residue_identity(residue_items: Sequence[Any]) -> dict[str, list[Any]]:
    """Return residue identity arrays for grouped target atoms."""

    return {
        "residue_keys": [
            f"{chain}:{resid}:{resname}" for (chain, resid, resname), _ in residue_items
        ],
        "residue_chainids": [chain for (chain, _, _), _ in residue_items],
        "residue_resids": [int(resid) for (_, resid, _), _ in residue_items],
        "residue_resnames": [resname for (_, _, resname), _ in residue_items],
    }


def _effective_frame_stride(base: FrameSelection, actual: FrameSelection) -> int | None:
    """Return the factor between base and job frame-selection strides."""

    if base.frames is not None or actual.frames is not None:
        return None
    base_step = base.step or 1
    actual_step = actual.step or 1
    if base_step <= 0:
        return actual_step
    return max(1, actual_step // base_step)


def _effective_timestep_ps(raw_timestep_ps: float | None, frame_stride: int | None) -> float | None:
    """Return sample spacing after stride when known."""

    if raw_timestep_ps is None:
        return None
    return float(raw_timestep_ps) * float(frame_stride or 1)


def _sidecar_path(run_label: str, order: int) -> str:
    """Return a stable store-relative sidecar path for a run."""

    token_payload = f"{order}|{run_label}"
    digest = hashlib.sha256(token_payload.encode("utf-8")).hexdigest()[:12]
    return f"{SASA_SIDECAR_PREFIX}_{order:02d}_{_safe_label(run_label)}_{digest}.npz"


def _job_name(run_label: str, order: int) -> str:
    """Return a unique MDA job name for a configured SASA run."""

    return f"sasa:{order}:{_safe_label(run_label)}"


def _safe_label(run_label: str) -> str:
    """Return a filesystem-safe run label token."""

    safe = "".join(char.lower() if char.isalnum() else "_" for char in run_label.strip())
    return "_".join(part for part in safe.split("_") if part) or "run"


def _json_float_or_none(value: Any) -> float | None:
    """Return a finite float or ``None`` for artifact-safe missing values."""

    if value is None:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def _legacy_float(value: Any) -> float:
    """Return a float for legacy result models, preserving missing values as NaN."""

    numeric = _json_float_or_none(value)
    return float("nan") if numeric is None else numeric


def _legacy_float_list(values: Sequence[Any]) -> list[float]:
    """Return floats for legacy result models from artifact-safe values."""

    return [_legacy_float(value) for value in values]


def _mean_optional_values(values: Sequence[float | None]) -> float | None:
    """Return the mean of finite values or ``None`` when none exist."""

    finite = [float(value) for value in values if value is not None]
    if not finite:
        return None
    return _json_float_or_none(np.mean(np.asarray(finite, dtype=np.float64)))


def _min_optional_values(values: Sequence[float | None]) -> float | None:
    """Return the minimum finite value or ``None`` when none exist."""

    finite = [float(value) for value in values if value is not None]
    return min(finite) if finite else None


def _max_optional_values(values: Sequence[float | None]) -> float | None:
    """Return the maximum finite value or ``None`` when none exist."""

    finite = [float(value) for value in values if value is not None]
    return max(finite) if finite else None


def _replicate_metrics(run_payloads: Sequence[Mapping[str, Any]]) -> dict[str, float | None]:
    """Return scalar replicate metrics keyed by run label."""

    return {
        _metric_key(str(run["run_label"])): _json_float_or_none(run.get("mean_sasa"))
        for run in run_payloads
    }


def _condition_metrics(run_results: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    """Return condition-level metric summaries keyed by run label."""

    metrics: dict[str, dict[str, Any]] = {}
    for run in run_results:
        values = [_json_float_or_none(value) for value in run["per_replicate_means"]]
        finite = [value for value in values if value is not None]
        if len(finite) > 1:
            std = _json_float_or_none(np.std(finite, ddof=1))
        else:
            std = 0.0 if finite else None
        metrics[_metric_key(str(run["run_label"]))] = {
            "name": _metric_key(str(run["run_label"])),
            "values": values,
            "mean": _json_float_or_none(run.get("overall_mean")),
            "sem": _json_float_or_none(run.get("overall_sem")),
            "std": std,
            "n": len(finite),
            **SASA_METRIC_METADATA,
        }
    return metrics


def _condition_replicate_metrics(
    run_results: Sequence[Mapping[str, Any]], replicates: Sequence[int]
) -> dict[str, dict[str, float | None]]:
    """Return per-replicate metric payloads for condition artifacts."""

    by_replicate = {str(int(replicate)): {} for replicate in replicates}
    for run in run_results:
        metric_key = _metric_key(str(run["run_label"]))
        for replicate, value in zip(replicates, run["per_replicate_means"], strict=True):
            by_replicate[str(int(replicate))][metric_key] = _json_float_or_none(value)
    return by_replicate


def _metric_metadata(run_payloads: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    """Return metric metadata for all SASA runs."""

    return {_metric_key(str(run["run_label"])): dict(SASA_METRIC_METADATA) for run in run_payloads}


def _metric_key(run_label: str) -> str:
    """Return the scalar metric key for a run label."""

    return f"{MEAN_SASA_METRIC}:{run_label}"


def _trajectory_files(universe_policy: Mapping[str, Any]) -> list[str]:
    """Extract trajectory file paths from universe policy provenance."""

    provenance = universe_policy.get("provenance")
    if not isinstance(provenance, Mapping):
        return []
    trajectory = provenance.get("trajectory")
    if isinstance(trajectory, Mapping):
        path = trajectory.get("path")
        return [str(path)] if path is not None else []
    if isinstance(trajectory, Sequence) and not isinstance(trajectory, (str, bytes)):
        paths = []
        for item in trajectory:
            if isinstance(item, Mapping) and item.get("path") is not None:
                paths.append(str(item["path"]))
        return paths
    return []


def _unique_warnings(warnings: Sequence[str]) -> list[str]:
    """Return de-duplicated warning messages preserving first occurrence."""

    unique: list[str] = []
    seen: set[str] = set()
    for warning in warnings:
        if warning in seen:
            continue
        unique.append(str(warning))
        seen.add(str(warning))
    return unique


def _combined_warnings(artifacts: Sequence[ReplicateArtifact]) -> list[str]:
    """Return de-duplicated warnings from source replicate artifacts."""

    warnings: list[str] = []
    seen: set[str] = set()
    for artifact in artifacts:
        for warning in artifact.warnings:
            if warning in seen:
                continue
            warnings.append(warning)
            seen.add(warning)
    return warnings


def mdanalysis_version() -> str | None:
    """Return the imported MDAnalysis version when available."""

    try:
        import MDAnalysis as mda
    except ImportError:
        return None
    return str(getattr(mda, "__version__", "unknown"))


def mdtraj_version() -> str | None:
    """Return the imported MDTraj version when available."""

    try:
        import mdtraj as md
    except ImportError:
        return None
    return str(getattr(md, "__version__", "unknown"))
