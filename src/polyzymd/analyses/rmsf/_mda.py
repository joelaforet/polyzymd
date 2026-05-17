"""MDAnalysis-native RMSF profile jobs and artifact adapters."""

from __future__ import annotations

import hashlib
import logging
import math
from collections.abc import Mapping, Sequence
from operator import index as operator_index
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
from polyzymd.analyses.rmsf._results import RMSFAggregatedResult, RMSFResult
from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory
from polyzymd.analyses.shared.config_hash import compute_config_hash
from polyzymd.analyses.shared.diagnostics import get_selection_diagnostics
from polyzymd.analyses.shared.loader import parse_time_string
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.mda import ArtifactSidecarRef, MDAReplicateJobContext
    from polyzymd.analyses.rmsf import RMSFSettings

LOGGER = logging.getLogger(__name__)
AUTOCORRELATION_FRAME_THRESHOLD = 100
MEAN_RMSF_METRIC = "mean_rmsf"
RMSF_PROFILE_VERSION = "1"
RMSF_METRIC_METADATA = {
    "label": "Mean RMSF",
    "unit": "A",
    "higher_is_better": False,
    "direction_labels": ("stabilizing", "unchanged", "destabilizing"),
    "statistical_policy": {
        "metric_classification": "variance_based",
        "frame_strategy": "subsample_by_autocorrelation_when_estimated",
        "uncertainty_level": "biological_replicate_sem",
    },
}


def build_rmsf_jobs(ctx: MDAReplicateJobContext, settings: RMSFSettings) -> list[MDAAnalysisJob]:
    """Build one MDAnalysis-compatible RMSF profile job for a replicate.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        Framework-provided MDAnalysis replicate context.
    settings : RMSFSettings
        Resolved RMSF settings.

    Returns
    -------
    list of MDAAnalysisJob
        Single RMSF profile job with explicit independent-frame selection.
    """

    job_input = prepare_rmsf_profile_input(
        universe=ctx.universe,
        settings=settings,
        frame_selection=ctx.frame_selection,
        condition_label=ctx.replicate_context.condition.label,
        replicate=ctx.replicate,
    )
    selected_frame_selection = _selected_frame_selection(ctx.frame_selection, job_input.frames)
    backend_policy = getattr(ctx, "backend_policy", None)
    return [
        MDAAnalysisJob(
            name="rmsf_profile",
            analysis=RMSFProfileAnalysis(
                job_input.atoms,
                reference_positions=job_input.reference_positions,
            ),
            frame_selection=selected_frame_selection,
            **({"backend_policy": backend_policy} if backend_policy is not None else {}),
            universe_policy=MDAUniversePolicy(
                condition_label=ctx.replicate_context.condition.label,
                replicate=ctx.replicate,
                provenance=ctx.universe_policy.provenance,
                metadata={
                    **ctx.universe_policy.metadata,
                    "rmsf_settings": settings.model_dump(mode="json"),
                    "rmsf_profile_version": RMSF_PROFILE_VERSION,
                    "reference_metadata": job_input.reference_metadata,
                    "alignment_metadata": job_input.alignment_metadata,
                    "autocorrelation_metadata": job_input.autocorrelation_metadata,
                    "residue_profile": job_input.residue_profile,
                },
            ),
        )
    ]


class RMSFArtifactCollector:
    """Collect completed RMSF profile jobs into a replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Map one completed RMSF job to JSON metadata and an NPZ sidecar.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context for one replicate.
        completed_jobs : sequence of MDAJobResult
            Completed RMSF jobs.

        Returns
        -------
        ReplicateArtifact
            Canonical replicate artifact with profile arrays in a sidecar.
        """

        if len(completed_jobs) != 1:
            raise ValueError(f"RMSF expected one profile job, got {len(completed_jobs)}")

        job = completed_jobs[0]
        metadata = dict(job.universe_policy.metadata)
        residue_profile = _require_mapping(metadata.get("residue_profile"), "residue_profile")
        reference_metadata = _require_mapping(
            metadata.get("reference_metadata"), "reference_metadata"
        )
        alignment_metadata = _require_mapping(
            metadata.get("alignment_metadata"), "alignment_metadata"
        )
        autocorrelation_metadata = _require_mapping(
            metadata.get("autocorrelation_metadata"), "autocorrelation_metadata"
        )
        rmsf_values = np.asarray(job.results.rmsf_values, dtype=np.float64)
        if rmsf_values.ndim != 1 or rmsf_values.size == 0 or not np.all(np.isfinite(rmsf_values)):
            raise ValueError("RMSF profile job produced invalid RMSF values")

        residue_ids = _int_array(residue_profile.get("residue_ids"), "residue_ids")
        residue_indices = _int_array(residue_profile.get("residue_indices"), "residue_indices")
        residue_names = [str(value) for value in residue_profile.get("residue_names", [])]
        identity_keys = [str(value) for value in residue_profile.get("identity_keys", [])]
        if not (
            len(residue_ids)
            == len(residue_indices)
            == len(residue_names)
            == len(identity_keys)
            == len(rmsf_values)
        ):
            raise ValueError("RMSF residue identity arrays do not match RMSF profile length")

        selected_frames = np.asarray(job.frame_selection.frames, dtype=np.int64)
        sidecar = ctx.artifact_store.write_npz_sidecar(
            "sidecars/rmsf_profile.npz",
            residue_ids=residue_ids,
            residue_indices=residue_indices,
            residue_names=np.asarray(residue_names),
            identity_keys=np.asarray(identity_keys),
            rmsf_values=rmsf_values,
            selected_frames=selected_frames,
            atom_counts_per_residue=_int_array(
                residue_profile.get("atom_counts_per_residue"), "atom_counts_per_residue"
            ),
            metadata={"kind": "rmsf_residue_profile", "job_name": job.name},
        )
        eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
        config_hash = compute_config_hash(ctx.replicate_context.sim_config)
        stats = _profile_statistics(rmsf_values)
        profile_payload = {
            "residue_ids": residue_ids.tolist(),
            "residue_indices": residue_indices.tolist(),
            "residue_names": residue_names,
            "identity_keys": identity_keys,
            "rmsf_values": rmsf_values.tolist(),
            "n_residues": int(rmsf_values.size),
            "sidecar": sidecar.model_dump(mode="json"),
            "npz_path": sidecar.path,
        }
        frame_metadata = {
            "n_frames_total": ctx.frame_selection.n_frames_total,
            "n_frames_window": autocorrelation_metadata.get("n_frames_window"),
            "n_frames_used": int(selected_frames.size),
            "selected_frames": selected_frames.tolist(),
            "timestep_ps": ctx.frame_selection.timestep_ps,
            "frame_stride": ctx.frame_selection.step,
        }
        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "profile": profile_payload,
                "residue_profile": profile_payload,
                "metrics": {MEAN_RMSF_METRIC: stats["mean_rmsf"]},
                "replicate_metrics": {MEAN_RMSF_METRIC: stats["mean_rmsf"]},
                "metric_metadata": {MEAN_RMSF_METRIC: RMSF_METRIC_METADATA},
                "frame_metadata": frame_metadata,
                "autocorrelation": dict(autocorrelation_metadata),
                "reference": dict(reference_metadata),
                "alignment": dict(alignment_metadata),
                **stats,
            },
            sidecars=[sidecar],
            provenance={
                "source": "mda_rmsf_profile_job",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "analysis_frame_selection": frame_selection_payload(job.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
            },
            metadata={
                "result_kind": "rmsf_mda_replicate",
                "settings_fingerprint": ctx.settings_fingerprint,
                "config_hash": config_hash,
                "polyzymd_version": get_polyzymd_version(),
                "mdanalysis_version": mdanalysis_version(),
                "rmsf_profile_version": RMSF_PROFILE_VERSION,
                "equilibration_time": eq_value,
                "equilibration_unit": eq_unit,
                "selection_string": str(ctx.settings.selection),
                "statistical_policy": RMSF_METRIC_METADATA["statistical_policy"],
            },
            warnings=_unique_warnings(ctx.warnings),
        )


def aggregate_rmsf_artifacts(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: RMSFSettings,
    equilibration: str,
    output_dir: Path,
    result_path: Path,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> ConditionArtifact:
    """Aggregate RMSF replicate artifacts into one condition artifact.

    Parameters
    ----------
    condition_label : str
        Condition label being aggregated.
    replicates : sequence of int
        Expected replicate IDs.
    settings : RMSFSettings
        Active RMSF settings.
    equilibration : str
        Equilibration string from the framework context.
    output_dir : Path
        Aggregated output directory.
    result_path : Path
        Canonical condition artifact path.
    artifacts : sequence of ReplicateArtifact
        Per-replicate RMSF artifacts.
    settings_fingerprint : str
        Active settings fingerprint.

    Returns
    -------
    ConditionArtifact
        Canonical condition artifact with per-residue profile summaries.
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
        replicates=replicates,
        settings=settings,
        equilibration=equilibration,
        artifacts=ordered_artifacts,
        settings_fingerprint=settings_fingerprint,
    )
    values = [float(value) for value in legacy_result.per_replicate_mean_rmsf]
    metric = _metric_summary(
        MEAN_RMSF_METRIC,
        values,
        legacy_result.overall_mean_rmsf,
        legacy_result.overall_sem_rmsf,
    )
    metric.update(RMSF_METRIC_METADATA)
    replicate_metrics = {
        str(replicate): {MEAN_RMSF_METRIC: value}
        for replicate, value in zip(legacy_result.replicates, values, strict=True)
    }
    artifact = ConditionArtifact(
        analysis_name="rmsf",
        condition_label=condition_label,
        replicates=[int(rep) for rep in replicates],
        payload={
            "profile": {
                "residue_ids": legacy_result.residue_ids,
                "residue_names": legacy_result.residue_names,
                "residue_indices": ordered_artifacts[0].payload["profile"]["residue_indices"],
                "identity_keys": ordered_artifacts[0].payload["profile"]["identity_keys"],
                "mean_rmsf_per_residue": legacy_result.mean_rmsf_per_residue,
                "sem_rmsf_per_residue": legacy_result.sem_rmsf_per_residue,
            },
            "residue_ids": legacy_result.residue_ids,
            "residue_names": legacy_result.residue_names,
            "mean_rmsf_per_residue": legacy_result.mean_rmsf_per_residue,
            "sem_rmsf_per_residue": legacy_result.sem_rmsf_per_residue,
            "per_replicate_mean_rmsf": values,
            "overall_mean_rmsf": legacy_result.overall_mean_rmsf,
            "overall_sem_rmsf": legacy_result.overall_sem_rmsf,
            "overall_min_rmsf": legacy_result.overall_min_rmsf,
            "overall_max_rmsf": legacy_result.overall_max_rmsf,
            "metrics": {MEAN_RMSF_METRIC: metric},
            "replicate_metrics": replicate_metrics,
            "metric_metadata": {MEAN_RMSF_METRIC: RMSF_METRIC_METADATA},
            "n_replicates": legacy_result.n_replicates,
        },
        provenance={
            "source": "rmsf_replicate_artifact_aggregation",
            "frame_selection": ordered_artifacts[0].provenance.get("frame_selection"),
            "statistical_policy": RMSF_METRIC_METADATA["statistical_policy"],
        },
        metadata={
            "result_kind": "rmsf_mda_condition",
            "settings_fingerprint": settings_fingerprint,
            "config_hash": legacy_result.config_hash,
            "polyzymd_version": legacy_result.polyzymd_version,
            "mdanalysis_version": mdanalysis_version(),
            "rmsf_profile_version": RMSF_PROFILE_VERSION,
            "equilibration_time": legacy_result.equilibration_time,
            "equilibration_unit": legacy_result.equilibration_unit,
            "selection_string": legacy_result.selection_string,
            "source_result_files": [
                str(output_dir.parent / f"run_{replicate}" / "result.json")
                for replicate in replicates
            ],
            "n_replicates": legacy_result.n_replicates,
            "statistical_policy": RMSF_METRIC_METADATA["statistical_policy"],
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


class RMSFProfileInput:
    """Prepared inputs for one RMSF profile job."""

    def __init__(
        self,
        *,
        atoms: Any,
        frames: NDArray[np.int64],
        reference_positions: NDArray[np.float64] | None,
        residue_profile: dict[str, Any],
        reference_metadata: dict[str, Any],
        alignment_metadata: dict[str, Any],
        autocorrelation_metadata: dict[str, Any],
    ) -> None:
        """Store prepared RMSF job inputs."""

        self.atoms = atoms
        self.frames = frames
        self.reference_positions = reference_positions
        self.residue_profile = residue_profile
        self.reference_metadata = reference_metadata
        self.alignment_metadata = alignment_metadata
        self.autocorrelation_metadata = autocorrelation_metadata


class RMSFProfileAnalysis:
    """Factory wrapper for an MDAnalysis ``AnalysisBase`` RMSF profile."""

    def __new__(cls, atoms: Any, reference_positions: NDArray[np.float64] | None = None) -> Any:
        """Return an ``AnalysisBase`` instance without importing MDAnalysis eagerly."""

        from MDAnalysis.analysis.base import AnalysisBase

        class _RMSFProfileAnalysis(AnalysisBase):
            """Compute per-residue RMSF over frames owned by ``AnalysisBase``."""

            def __init__(
                self,
                atom_group: Any,
                ref_positions: NDArray[np.float64] | None,
            ) -> None:
                super().__init__(atom_group.universe.trajectory)
                self._atom_group = atom_group
                self._reference_positions = (
                    None if ref_positions is None else np.asarray(ref_positions, dtype=np.float64)
                )

            def _prepare(self) -> None:
                """Initialize streaming accumulators."""

                n_atoms = len(self._atom_group)
                self._n_frames = 0
                if self._reference_positions is None:
                    self._mean = np.zeros((n_atoms, 3), dtype=np.float64)
                    self._m2 = np.zeros((n_atoms, 3), dtype=np.float64)
                else:
                    self._sq_diff_sum = np.zeros(n_atoms, dtype=np.float64)

            def _single_frame(self) -> None:
                """Accumulate one trajectory frame."""

                positions = self._atom_group.positions.astype(np.float64)
                self._n_frames += 1
                if self._reference_positions is None:
                    delta = positions - self._mean
                    self._mean += delta / float(self._n_frames)
                    delta2 = positions - self._mean
                    self._m2 += delta * delta2
                    return
                diff = positions - self._reference_positions
                self._sq_diff_sum += np.sum(diff**2, axis=1)

            def _conclude(self) -> None:
                """Store per-residue RMSF profile values."""

                if self._n_frames == 0:
                    raise ValueError("RMSF analysis requires at least one frame")
                if self._reference_positions is None:
                    atom_rmsf = np.sqrt(np.sum(self._m2 / float(self._n_frames), axis=1))
                else:
                    atom_rmsf = np.sqrt(self._sq_diff_sum / float(self._n_frames))
                self.results.atom_rmsf_values = np.asarray(atom_rmsf, dtype=np.float64)
                self.results.rmsf_values = aggregate_atom_rmsf_per_residue(
                    self._atom_group,
                    atom_rmsf,
                )

        return _RMSFProfileAnalysis(atoms, reference_positions)


def prepare_rmsf_profile_input(
    *,
    universe: Any,
    settings: RMSFSettings,
    frame_selection: FrameSelection,
    condition_label: str,
    replicate: int,
) -> RMSFProfileInput:
    """Validate selections, align the trajectory, and choose RMSF frames.

    Parameters
    ----------
    universe : Any
        Loaded MDAnalysis universe.
    settings : RMSFSettings
        Active RMSF settings.
    frame_selection : FrameSelection
        Framework-resolved frame window.
    condition_label : str
        Condition label for diagnostics.
    replicate : int
        Replicate ID for diagnostics.

    Returns
    -------
    RMSFProfileInput
        Prepared atom group, selected frames, reference positions, and metadata.
    """

    _validate_explicit_frame_reference_mode(frame_selection, settings)
    alignment_start, alignment_stop, alignment_step = _alignment_bounds(frame_selection)
    frames = _selected_frames_from_frame_selection(frame_selection)
    _validate_reference_settings(settings)
    atoms = universe.select_atoms(settings.selection)
    if len(atoms) == 0:
        diagnostics = get_selection_diagnostics(universe, settings.selection)
        raise ValueError(
            f"RMSF selection {settings.selection!r} matched no atoms for condition "
            f"'{condition_label}' replicate {replicate}.\n\n{diagnostics}"
        )
    _validate_external_reference_selections_before_alignment(
        universe=universe,
        settings=settings,
        atoms=atoms,
        condition_label=condition_label,
        replicate=replicate,
    )

    alignment_config = AlignmentConfig(
        enabled=True,
        reference_mode=settings.reference_mode,
        reference_frame=settings.reference_frame,
        selection=settings.alignment_selection,
        centroid_selection=settings.centroid_selection,
        reference_file=(
            Path(settings.reference_file) if settings.reference_file is not None else None
        ),
    )
    ref_frame_idx = align_trajectory(
        universe,
        alignment_config,
        start_frame=alignment_start,
        stop_frame=alignment_stop,
        step_frame=alignment_step,
    )
    autocorrelation_metadata = _autocorrelation_metadata(
        universe=universe,
        atoms=atoms,
        frames=frames,
        timestep_ps=float(frame_selection.timestep_ps or 1.0),
        step=alignment_step,
        start=alignment_start,
        stop=alignment_stop,
        explicit_frame_selection=frame_selection.frames is not None,
    )
    selected_frames = np.asarray(autocorrelation_metadata["selected_frames"], dtype=np.int64)
    reference_positions = _external_reference_positions(universe, atoms, settings)
    reference_file_identity = (
        external_reference_file_identity(settings.reference_file)
        if settings.reference_mode == "external"
        else None
    )
    residue_profile = residue_profile_identity(atoms)
    return RMSFProfileInput(
        atoms=atoms,
        frames=selected_frames,
        reference_positions=reference_positions,
        residue_profile=residue_profile,
        reference_metadata={
            "reference_mode": settings.reference_mode,
            "reference_frame": ref_frame_idx + 1 if ref_frame_idx is not None else None,
            "reference_file": settings.reference_file,
            "reference_file_identity": reference_file_identity,
            "external_reference_used": reference_positions is not None,
            "trajectory_selection_identity": _selection_identity_payload(atoms),
        },
        alignment_metadata={
            "alignment_selection": settings.alignment_selection,
            "centroid_selection": settings.centroid_selection,
        },
        autocorrelation_metadata=autocorrelation_metadata,
    )


def aggregate_atom_rmsf_per_residue(
    atoms: Any, atom_rmsf: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Reduce per-atom RMSF values to residue means.

    Parameters
    ----------
    atoms : Any
        Selected atom group.
    atom_rmsf : ndarray of float
        Per-atom RMSF values in selection order.

    Returns
    -------
    ndarray of float
        Per-residue RMSF profile.
    """

    atom_indices = np.asarray(atoms.indices, dtype=np.int64)
    profile = np.zeros(len(atoms.residues), dtype=np.float64)
    for index, residue in enumerate(atoms.residues):
        residue_indices = np.asarray(residue.atoms.indices, dtype=np.int64)
        mask = np.isin(atom_indices, residue_indices)
        if not np.any(mask):
            raise ValueError(f"Residue {getattr(residue, 'resid', index)} has no selected atoms")
        profile[index] = float(np.mean(atom_rmsf[mask]))
    return profile


def residue_profile_identity(atoms: Any) -> dict[str, Any]:
    """Return residue identity arrays for the selected atom group.

    Parameters
    ----------
    atoms : Any
        Selected atom group.

    Returns
    -------
    dict[str, Any]
        JSON-compatible residue IDs, indices, names, identity keys, and atom counts.
    """

    atom_indices = np.asarray(atoms.indices, dtype=np.int64)
    residue_ids: list[int] = []
    residue_indices: list[int] = []
    residue_names: list[str] = []
    identity_keys: list[str] = []
    atom_counts: list[int] = []
    for fallback_index, residue in enumerate(atoms.residues):
        residue_index = int(getattr(residue, "ix", fallback_index))
        resid = int(residue.resid)
        resname = str(getattr(residue, "resname", ""))
        segid = str(getattr(getattr(residue, "segment", None), "segid", ""))
        residue_atom_indices = np.asarray(residue.atoms.indices, dtype=np.int64)
        atom_count = int(np.count_nonzero(np.isin(atom_indices, residue_atom_indices)))
        residue_ids.append(resid)
        residue_indices.append(residue_index)
        residue_names.append(resname)
        identity_keys.append(f"{segid}:{residue_index}:{resid}:{resname}")
        atom_counts.append(atom_count)
    if not residue_ids:
        raise ValueError("RMSF selection produced no residues")
    return {
        "residue_ids": residue_ids,
        "residue_indices": residue_indices,
        "residue_names": residue_names,
        "identity_keys": identity_keys,
        "atom_counts_per_residue": atom_counts,
    }


def external_reference_file_identity(reference_file: str | Path | None) -> dict[str, Any] | None:
    """Return stable identity metadata for an external reference file.

    Parameters
    ----------
    reference_file : str, Path, or None
        Configured external reference path.

    Returns
    -------
    dict or None
        JSON-compatible path, size, and SHA-256 metadata, or ``None`` when no
        reference file is configured.
    """

    if reference_file is None:
        return None
    path = Path(reference_file).expanduser()
    resolved = path.resolve(strict=False)
    identity: dict[str, Any] = {
        "path": str(reference_file),
        "resolved_path": str(resolved),
        "exists": path.exists(),
    }
    if not path.exists():
        return identity

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    stat = path.stat()
    identity.update(
        {
            "size_bytes": int(stat.st_size),
            "sha256": digest.hexdigest(),
        }
    )
    return identity


def artifact_to_rmsf_result(artifact: ReplicateArtifact) -> RMSFResult:
    """Convert a replicate artifact to the established RMSF result model."""

    if artifact.analysis_name != "rmsf":
        raise ValueError(f"Expected rmsf artifact, got {artifact.analysis_name!r}")
    metadata = artifact.metadata
    payload = artifact.payload
    profile = _require_mapping(payload.get("profile"), "profile")
    reference = _require_mapping(payload.get("reference"), "reference")
    alignment = _require_mapping(payload.get("alignment"), "alignment")
    autocorrelation = _require_mapping(payload.get("autocorrelation"), "autocorrelation")
    frame_metadata = _require_mapping(payload.get("frame_metadata"), "frame_metadata")
    return RMSFResult(
        config_hash=str(metadata.get("config_hash", "unknown")),
        polyzymd_version=str(metadata.get("polyzymd_version", get_polyzymd_version())),
        replicate=artifact.replicate,
        equilibration_time=float(metadata.get("equilibration_time", 0.0)),
        equilibration_unit=str(metadata.get("equilibration_unit", "ns")),
        selection_string=str(metadata.get("selection_string", "")),
        correlation_time=autocorrelation.get("correlation_time"),
        correlation_time_unit=autocorrelation.get("correlation_time_unit"),
        n_independent_frames=autocorrelation.get("n_independent_frames"),
        residue_ids=[int(value) for value in profile.get("residue_ids", [])],
        residue_names=[str(value) for value in profile.get("residue_names", [])],
        rmsf_values=[float(value) for value in profile.get("rmsf_values", [])],
        mean_rmsf=float(payload.get("mean_rmsf", 0.0)),
        std_rmsf=float(payload.get("std_rmsf", 0.0)),
        min_rmsf=float(payload.get("min_rmsf", 0.0)),
        max_rmsf=float(payload.get("max_rmsf", 0.0)),
        reference_mode=reference.get("reference_mode"),
        reference_frame=reference.get("reference_frame"),
        alignment_selection=alignment.get("alignment_selection"),
        reference_file=reference.get("reference_file"),
        n_frames_total=int(frame_metadata.get("n_frames_total", 0) or 0),
        n_frames_used=int(frame_metadata.get("n_frames_used", 0) or 0),
        settings_fingerprint=metadata.get("settings_fingerprint"),
        trajectory_files=[],
    )


def condition_artifact_to_legacy_result(artifact: Any) -> RMSFAggregatedResult:
    """Convert a condition artifact to the established aggregate model."""

    if isinstance(artifact, RMSFAggregatedResult):
        return artifact
    if not isinstance(artifact, ConditionArtifact):
        raise TypeError(
            "RMSF aggregate adapters require canonical ConditionArtifact inputs. "
            f"Got {type(artifact).__name__}; recompute the condition or clear stale caches."
        )
    metadata = artifact.metadata
    payload = artifact.payload
    return RMSFAggregatedResult(
        config_hash=str(metadata.get("config_hash", "unknown")),
        polyzymd_version=str(metadata.get("polyzymd_version", get_polyzymd_version())),
        replicate=None,
        equilibration_time=float(metadata.get("equilibration_time", 0.0)),
        equilibration_unit=str(metadata.get("equilibration_unit", "ns")),
        selection_string=str(metadata.get("selection_string", "")),
        replicates=[int(rep) for rep in artifact.replicates],
        n_replicates=len(artifact.replicates),
        residue_ids=[int(value) for value in payload.get("residue_ids", [])],
        residue_names=[str(value) for value in payload.get("residue_names", [])],
        mean_rmsf_per_residue=[float(value) for value in payload.get("mean_rmsf_per_residue", [])],
        sem_rmsf_per_residue=[float(value) for value in payload.get("sem_rmsf_per_residue", [])],
        per_replicate_mean_rmsf=[
            float(value) for value in payload.get("per_replicate_mean_rmsf", [])
        ],
        overall_mean_rmsf=float(payload.get("overall_mean_rmsf", 0.0)),
        overall_sem_rmsf=float(payload.get("overall_sem_rmsf", 0.0)),
        overall_min_rmsf=float(payload.get("overall_min_rmsf", 0.0)),
        overall_max_rmsf=float(payload.get("overall_max_rmsf", 0.0)),
        source_result_files=[str(path) for path in metadata.get("source_result_files", [])],
        settings_fingerprint=metadata.get("settings_fingerprint"),
    )


def load_condition_artifact(aggregated_dir: Path) -> ConditionArtifact | None:
    """Load the canonical RMSF condition artifact if present."""

    result_path = aggregated_dir / "result.json"
    if not result_path.exists():
        return None
    return ArtifactStore(aggregated_dir).read_condition_result("result.json")


def mdanalysis_version() -> str | None:
    """Return the MDAnalysis version without importing it at module import time."""

    try:
        import MDAnalysis as mda
    except ImportError:
        return None
    return str(getattr(mda, "__version__", "unknown"))


def _aggregate_legacy_result(
    *,
    replicates: Sequence[int],
    settings: RMSFSettings,
    equilibration: str,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> RMSFAggregatedResult:
    """Build the established aggregate model from replicate artifacts."""

    legacy_results = [artifact_to_rmsf_result(artifact) for artifact in artifacts]
    first = legacy_results[0]
    eq_value, eq_unit = parse_time_string(equilibration)
    per_replicate_rmsf = [
        np.asarray(result.rmsf_values, dtype=np.float64) for result in legacy_results
    ]
    matrix = np.vstack(per_replicate_rmsf)
    means = np.mean(matrix, axis=0)
    sems = np.zeros(matrix.shape[1], dtype=np.float64)
    if matrix.shape[0] > 1:
        sems = np.std(matrix, axis=0, ddof=1) / math.sqrt(matrix.shape[0])
    per_replicate_means = [float(result.mean_rmsf) for result in legacy_results]
    overall_stats = compute_sem(per_replicate_means)
    return RMSFAggregatedResult(
        config_hash=first.config_hash,
        polyzymd_version=get_polyzymd_version(),
        replicate=None,
        equilibration_time=eq_value,
        equilibration_unit=eq_unit,
        selection_string=settings.selection,
        replicates=[int(rep) for rep in replicates],
        n_replicates=len(replicates),
        residue_ids=first.residue_ids,
        residue_names=first.residue_names,
        mean_rmsf_per_residue=means.tolist(),
        sem_rmsf_per_residue=sems.tolist(),
        per_replicate_mean_rmsf=per_replicate_means,
        overall_mean_rmsf=overall_stats.mean,
        overall_sem_rmsf=overall_stats.sem,
        overall_min_rmsf=float(np.min(means)),
        overall_max_rmsf=float(np.max(means)),
        settings_fingerprint=settings_fingerprint,
        source_result_files=[],
    )


def _validate_and_order_artifacts(
    *,
    condition_label: str,
    expected_replicates: Sequence[int],
    settings: RMSFSettings,
    settings_fingerprint: str,
    artifacts: Sequence[ReplicateArtifact],
    analysis_dir: Path,
) -> list[ReplicateArtifact]:
    """Validate RMSF artifact identity, settings, residue identity, and sidecars."""

    expected = [int(rep) for rep in expected_replicates]
    by_replicate: dict[int, ReplicateArtifact] = {}
    expected_identity: tuple[str, ...] | None = None
    for artifact in artifacts:
        if artifact.analysis_name != "rmsf":
            raise ValueError(f"Expected rmsf artifact, got {artifact.analysis_name!r}")
        if artifact.condition_label != condition_label:
            raise ValueError(
                f"RMSF artifact condition mismatch: expected {condition_label!r}, "
                f"got {artifact.condition_label!r}"
            )
        if artifact.replicate in by_replicate:
            raise ValueError(f"Duplicate RMSF artifact for replicate {artifact.replicate}")
        if artifact.metadata.get("settings_fingerprint") != settings_fingerprint:
            raise ValueError(
                f"RMSF artifact replicate {artifact.replicate} has settings fingerprint "
                f"{artifact.metadata.get('settings_fingerprint')}, expected {settings_fingerprint}. "
                "Recompute the condition or clear stale caches."
            )
        if artifact.metadata.get("selection_string") != settings.selection:
            raise ValueError(f"RMSF replicate {artifact.replicate} selection mismatch")
        _validate_reference_file_identity(artifact, settings)
        identity = tuple(_profile_identity_keys(artifact))
        if expected_identity is None:
            expected_identity = identity
        elif identity != expected_identity:
            raise ValueError(
                f"RMSF residue identity mismatch for replicate {artifact.replicate}. "
                "Recompute with consistent topology/selection before aggregating."
            )
        store = ArtifactStore(analysis_dir / f"run_{artifact.replicate}")
        store.validate_sidecar(_rmsf_profile_sidecar(artifact))
        by_replicate[int(artifact.replicate)] = artifact
    observed = sorted(by_replicate)
    if observed != sorted(expected):
        raise ValueError(
            f"RMSF artifacts for condition '{condition_label}' do not match expected replicates "
            f"{expected}: found {observed}. Recompute missing replicates or clear stale caches "
            "before aggregating."
        )
    return [by_replicate[rep] for rep in expected]


def _validate_reference_settings(settings: RMSFSettings) -> None:
    """Validate reference settings that require files or frame IDs."""

    if settings.reference_mode == "frame" and settings.reference_frame is None:
        raise ValueError("reference_frame is required when reference_mode='frame'")
    if settings.reference_mode == "external":
        if settings.reference_file is None:
            raise ValueError(
                "reference_file is required when reference_mode='external'. Provide a valid PDB."
            )
        ref_path = Path(settings.reference_file)
        if not ref_path.exists():
            raise ValueError(
                f"reference_file does not exist: {ref_path}. Provide a valid external PDB."
            )


def _validate_reference_file_identity(artifact: ReplicateArtifact, settings: RMSFSettings) -> None:
    """Validate cached external reference file identity against current contents."""

    if settings.reference_mode != "external":
        return
    reference = _require_mapping(artifact.payload.get("reference"), "reference")
    expected_identity = external_reference_file_identity(settings.reference_file)
    stored_identity = reference.get("reference_file_identity")
    if stored_identity != expected_identity:
        raise ValueError(
            f"RMSF replicate {artifact.replicate} reference file identity mismatch. "
            "The configured external reference file changed after this artifact was cached; "
            "recompute the condition or clear stale caches."
        )


def _validate_external_reference_selections_before_alignment(
    *,
    universe: Any,
    settings: RMSFSettings,
    atoms: Any,
    condition_label: str,
    replicate: int,
) -> None:
    """Validate external-reference selections before alignment mutates coordinates."""

    if settings.reference_mode != "external" or settings.reference_file is None:
        return
    ref_path = Path(settings.reference_file)
    ref_universe = _load_external_reference_universe(ref_path)
    selections = [
        ("RMSF", settings.selection, atoms),
        ("alignment", settings.alignment_selection, None),
    ]
    for purpose, selection, selected_atoms in selections:
        trajectory_atoms = (
            selected_atoms if selected_atoms is not None else universe.select_atoms(selection)
        )
        if len(trajectory_atoms) == 0:
            diagnostics = get_selection_diagnostics(universe, selection)
            raise ValueError(
                f"RMSF external-reference {purpose.lower()} selection {selection!r} matched no "
                f"trajectory atoms for condition '{condition_label}' replicate {replicate}.\n\n"
                f"{diagnostics}"
            )
        ref_atoms = ref_universe.select_atoms(selection)
        _validate_reference_selection_identity(
            trajectory_atoms=trajectory_atoms,
            reference_atoms=ref_atoms,
            reference_path=ref_path,
            selection=selection,
            purpose=purpose,
        )


def _external_reference_positions(
    universe: Any,
    atoms: Any,
    settings: RMSFSettings,
) -> NDArray[np.float64] | None:
    """Load external RMSF reference positions when configured."""

    del universe
    if settings.reference_mode != "external" or settings.reference_file is None:
        return None

    ref_path = Path(settings.reference_file)
    ref_universe = _load_external_reference_universe(ref_path)
    ref_atoms = ref_universe.select_atoms(settings.selection)
    _validate_reference_selection_identity(
        trajectory_atoms=atoms,
        reference_atoms=ref_atoms,
        reference_path=ref_path,
        selection=settings.selection,
        purpose="RMSF",
    )
    return ref_atoms.positions.copy().astype(np.float64)


def _load_external_reference_universe(reference_path: Path) -> Any:
    """Load an external reference universe lazily."""

    import MDAnalysis as mda

    return mda.Universe(str(reference_path))


def _validate_reference_selection_identity(
    *,
    trajectory_atoms: Any,
    reference_atoms: Any,
    reference_path: Path,
    selection: str,
    purpose: str,
) -> None:
    """Validate external reference atom and residue identity/order."""

    if len(reference_atoms) == 0:
        raise ValueError(
            f"External RMSF reference {reference_path} has no atoms matching {purpose.lower()} "
            f"selection {selection!r}."
        )
    if len(reference_atoms) != len(trajectory_atoms):
        raise ValueError(
            f"External RMSF reference atom count ({len(reference_atoms)}) does not match "
            f"trajectory {purpose.lower()} selection ({len(trajectory_atoms)}) for {selection!r}."
        )

    trajectory_residues = residue_profile_identity(trajectory_atoms)
    reference_residues = residue_profile_identity(reference_atoms)
    residue_fields = ("identity_keys", "atom_counts_per_residue")
    for field_name in residue_fields:
        if trajectory_residues[field_name] != reference_residues[field_name]:
            raise ValueError(
                f"External RMSF reference {reference_path} residue identity/order does not match "
                f"the trajectory {purpose.lower()} selection {selection!r}. Use a reference "
                "with the same selected residues in the same order."
            )

    trajectory_atom_keys = _atom_identity_keys(trajectory_atoms)
    reference_atom_keys = _atom_identity_keys(reference_atoms)
    if trajectory_atom_keys and reference_atom_keys and trajectory_atom_keys != reference_atom_keys:
        raise ValueError(
            f"External RMSF reference {reference_path} atom identity/order does not match the "
            f"trajectory {purpose.lower()} selection {selection!r}. Use a reference with the "
            "same selected atoms in the same order."
        )


def _selection_identity_payload(atoms: Any) -> dict[str, Any]:
    """Return residue and atom identity metadata for a selected atom group."""

    residue_identity = residue_profile_identity(atoms)
    return {
        "residue_identity_keys": residue_identity["identity_keys"],
        "atom_counts_per_residue": residue_identity["atom_counts_per_residue"],
        "atom_identity_keys": _atom_identity_keys(atoms),
    }


def _atom_identity_keys(atoms: Any) -> list[str]:
    """Return selected atom identity keys in selection order when available."""

    try:
        iterator = iter(atoms)
    except TypeError:
        return []

    keys: list[str] = []
    for selected_index, atom in enumerate(iterator):
        residue = getattr(atom, "residue", None)
        residue_index = getattr(residue, "ix", "")
        resid = getattr(residue, "resid", getattr(atom, "resid", ""))
        resname = getattr(residue, "resname", getattr(atom, "resname", ""))
        segid = getattr(getattr(residue, "segment", None), "segid", getattr(atom, "segid", ""))
        atom_name = getattr(atom, "name", "")
        keys.append(f"{selected_index}:{segid}:{residue_index}:{resid}:{resname}:{atom_name}")
    return keys


def _autocorrelation_metadata(
    *,
    universe: Any,
    atoms: Any,
    frames: NDArray[np.int64],
    timestep_ps: float,
    step: int,
    start: int,
    stop: int,
    explicit_frame_selection: bool = False,
) -> dict[str, Any]:
    """Return independent-frame metadata for variance-based RMSF."""

    metadata: dict[str, Any] = {
        "autocorrelation_analyzed": False,
        "correlation_time": None,
        "correlation_time_unit": None,
        "n_independent_frames": None,
        "statistical_inefficiency": None,
        "warning": None,
        "n_frames_window": int(frames.size),
        "selected_frames": frames.tolist(),
        "explicit_frame_selection": explicit_frame_selection,
    }
    if explicit_frame_selection:
        metadata["warning"] = (
            "RMSF autocorrelation subsampling is skipped for explicit frame selections; "
            "the provided frames are preserved exactly."
        )
        return metadata
    if frames.size <= AUTOCORRELATION_FRAME_THRESHOLD:
        return metadata

    from polyzymd.analyses.shared.autocorrelation import (
        compute_acf,
        estimate_correlation_time,
        get_independent_indices,
    )

    rmsd_values = compute_rmsd_timeseries(atoms, start_frame=start, stop_frame=stop, step=step)
    acf_result = compute_acf(rmsd_values, timestep=timestep_ps * step, timestep_unit="ps")
    tau_result = estimate_correlation_time(acf_result, n_frames=int(frames.size))
    independent_frames = get_independent_indices(
        n_frames=stop,
        correlation_time=tau_result.tau,
        timestep=timestep_ps,
        start_frame=start,
    )
    independent_frames = independent_frames[independent_frames < stop]
    if step > 1:
        independent_frames = independent_frames[np.isin(independent_frames, frames)]
    if independent_frames.size == 0:
        independent_frames = frames

    metadata.update(
        {
            "autocorrelation_analyzed": True,
            "correlation_time": tau_result.tau,
            "correlation_time_unit": tau_result.tau_unit,
            "n_independent_frames": tau_result.n_independent,
            "statistical_inefficiency": tau_result.statistical_inefficiency,
            "warning": tau_result.warning,
            "selected_frames": independent_frames.astype(np.int64).tolist(),
        }
    )
    return metadata


def compute_rmsd_timeseries(
    atoms: Any,
    *,
    start_frame: int,
    stop_frame: int,
    step: int,
) -> NDArray[np.float64]:
    """Compute an RMSD series through an ``AnalysisBase`` helper."""

    atoms.universe.trajectory[start_frame]
    reference_positions = atoms.positions.copy().astype(np.float64)
    analysis = _build_rmsd_timeseries_analysis(atoms, reference_positions).run(
        start=start_frame,
        stop=stop_frame,
        step=step,
    )
    return np.asarray(analysis.results.rmsd_values, dtype=np.float64)


def _build_rmsd_timeseries_analysis(atoms: Any, reference_positions: NDArray[np.float64]) -> Any:
    """Build an AnalysisBase object for RMSD autocorrelation input."""

    from MDAnalysis.analysis.base import AnalysisBase

    class _RMSDTimeseriesAnalysis(AnalysisBase):
        """Collect RMSD values against fixed reference coordinates."""

        def __init__(self, atom_group: Any, ref_positions: NDArray[np.float64]) -> None:
            super().__init__(atom_group.universe.trajectory)
            self._atom_group = atom_group
            self._ref_positions = ref_positions

        def _prepare(self) -> None:
            """Initialize the RMSD value list."""

            self.results.rmsd_values = []

        def _single_frame(self) -> None:
            """Append the current-frame RMSD."""

            diff = self._atom_group.positions - self._ref_positions
            self.results.rmsd_values.append(float(np.sqrt(np.mean(np.sum(diff**2, axis=1)))))

        def _conclude(self) -> None:
            """Store collected RMSD values as an array."""

            self.results.rmsd_values = np.asarray(self.results.rmsd_values, dtype=np.float64)

    return _RMSDTimeseriesAnalysis(atoms, reference_positions)


def _selected_frame_selection(base: FrameSelection, frames: NDArray[np.int64]) -> FrameSelection:
    """Return an explicit-frame selection preserving window provenance."""

    return FrameSelection(
        frames=[int(frame) for frame in frames],
        equilibration=base.equilibration,
        equilibration_start=base.equilibration_start,
        equilibration_ps=base.equilibration_ps,
        timestep_ps=base.timestep_ps,
        n_frames_total=base.n_frames_total,
        warning_message=base.warning_message,
    )


def _validate_explicit_frame_reference_mode(
    frame_selection: FrameSelection,
    settings: RMSFSettings,
) -> None:
    """Reject explicit frames for reference modes that build references from slices."""

    if frame_selection.frames is None:
        return
    if settings.reference_mode not in {"average", "centroid"}:
        return
    raise ValueError(
        "RMSF explicit FrameSelection.frames is not supported with "
        f"reference_mode={settings.reference_mode!r} because that reference mode builds the "
        "alignment/reference structure from a contiguous trajectory slice. Use reference_mode "
        "'frame' or 'external', or provide a start/stop/step slice selection."
    )


def _selected_frames_from_frame_selection(frame_selection: FrameSelection) -> NDArray[np.int64]:
    """Return exact frame indices selected for RMSF analysis."""

    if frame_selection.frames is not None:
        return _explicit_frame_indices(frame_selection)
    if frame_selection.n_frames_total is None:
        raise ValueError("RMSF jobs require known total frame count")
    start, stop, step = _alignment_bounds(frame_selection)
    return np.arange(start, stop, step, dtype=np.int64)


def _alignment_bounds(frame_selection: FrameSelection) -> tuple[int, int, int]:
    """Return slice bounds used only for pre-RMSF alignment/reference helpers."""

    if frame_selection.frames is not None:
        frame_array = _explicit_frame_indices(frame_selection)
        return int(np.min(frame_array)), int(np.max(frame_array)) + 1, 1
    if frame_selection.n_frames_total is None:
        raise ValueError("RMSF jobs require known total frame count")
    return (
        0 if frame_selection.start is None else int(frame_selection.start),
        (
            frame_selection.n_frames_total
            if frame_selection.stop is None
            else int(frame_selection.stop)
        ),
        1 if frame_selection.step is None else int(frame_selection.step),
    )


def _explicit_frame_indices(frame_selection: FrameSelection) -> NDArray[np.int64]:
    """Normalize an explicit frame selector to exact integer frame indices."""

    frames = frame_selection.frames
    if frames is None:
        raise ValueError("RMSF explicit frame indices require FrameSelection.frames")
    raw_frames = tuple(frames)
    if not raw_frames:
        raise ValueError("RMSF frame selection cannot be empty")
    if all(_is_boolean_frame_value(frame) for frame in raw_frames):
        if (
            frame_selection.n_frames_total is not None
            and len(raw_frames) != frame_selection.n_frames_total
        ):
            raise ValueError(
                "RMSF boolean frame mask length "
                f"{len(raw_frames)} does not match trajectory length "
                f"{frame_selection.n_frames_total}"
            )
        frame_array = np.asarray(raw_frames, dtype=bool)
        selected = np.flatnonzero(frame_array).astype(np.int64)
    else:
        if not all(_is_integer_frame_value(frame) for frame in raw_frames):
            raise ValueError("RMSF explicit frames must be integer indices or a boolean mask")
        selected = np.asarray([operator_index(frame) for frame in raw_frames], dtype=np.int64)

    if selected.size == 0:
        raise ValueError("RMSF frame selection cannot be empty")
    if np.any(selected < 0):
        raise ValueError("RMSF explicit frame indices must be non-negative")
    if frame_selection.n_frames_total is not None and np.any(
        selected >= frame_selection.n_frames_total
    ):
        raise ValueError(
            "RMSF explicit frame indices must be within trajectory range "
            f"[0, {frame_selection.n_frames_total})"
        )
    return selected


def _is_boolean_frame_value(frame: Any) -> bool:
    """Return whether a frame selector value is boolean-like."""

    if isinstance(frame, bool):
        return True
    frame_type = type(frame)
    return (
        frame_type.__name__ == "bool_"
        and frame_type.__module__.split(".", maxsplit=1)[0] == "numpy"
    )


def _is_integer_frame_value(frame: Any) -> bool:
    """Return whether a frame selector value is integer-like."""

    if _is_boolean_frame_value(frame):
        return False
    try:
        operator_index(frame)
    except TypeError:
        return False
    return True


def _profile_statistics(rmsf_values: NDArray[np.float64]) -> dict[str, float]:
    """Return finite scalar RMSF profile statistics."""

    return {
        "mean_rmsf": float(np.mean(rmsf_values)),
        "std_rmsf": float(np.std(rmsf_values, ddof=0)),
        "min_rmsf": float(np.min(rmsf_values)),
        "max_rmsf": float(np.max(rmsf_values)),
    }


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


def _rmsf_profile_sidecar(artifact: ReplicateArtifact) -> ArtifactSidecarRef:
    """Return the RMSF profile sidecar from a replicate artifact."""

    matches = [
        sidecar
        for sidecar in artifact.sidecars
        if sidecar.metadata.get("kind") == "rmsf_residue_profile"
    ]
    if len(matches) != 1:
        raise ValueError(
            f"RMSF replicate {artifact.replicate} requires exactly one profile sidecar, "
            f"found {len(matches)}"
        )
    return matches[0]


def _profile_identity_keys(artifact: ReplicateArtifact) -> list[str]:
    """Return identity keys from a replicate artifact profile."""

    profile = _require_mapping(artifact.payload.get("profile"), "profile")
    keys = profile.get("identity_keys")
    if not isinstance(keys, Sequence) or isinstance(keys, (str, bytes, bytearray)):
        raise ValueError(f"RMSF replicate {artifact.replicate} is missing residue identity keys")
    return [str(key) for key in keys]


def _int_array(value: Any, field_name: str) -> NDArray[np.int64]:
    """Validate an integer sequence and return an array."""

    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise ValueError(f"RMSF field {field_name!r} must be a sequence")
    return np.asarray([int(item) for item in value], dtype=np.int64)


def _require_mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    """Validate a mapping payload field."""

    if not isinstance(value, Mapping):
        raise ValueError(f"RMSF payload field {field_name!r} must be a mapping")
    return value


def _combined_warnings(artifacts: Sequence[ReplicateArtifact]) -> list[str]:
    """Combine warnings from replicate artifacts."""

    return _unique_warnings([warning for artifact in artifacts for warning in artifact.warnings])


def _unique_warnings(warnings: Sequence[str]) -> list[str]:
    """Return unique warnings while preserving first-seen order."""

    unique: list[str] = []
    seen: set[str] = set()
    for warning in warnings:
        text = str(warning)
        if text in seen:
            continue
        seen.add(text)
        unique.append(text)
    return unique
