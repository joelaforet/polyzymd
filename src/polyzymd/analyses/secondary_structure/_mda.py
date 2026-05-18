"""MDAnalysis-compatible DSSP jobs and secondary-structure artifacts."""

from __future__ import annotations

import logging
import math
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses._framework.cache_identity import compute_config_hash
from polyzymd.analyses._framework.results_base import get_polyzymd_version
from polyzymd.analyses.mda import (
    ArtifactStore,
    ConditionArtifact,
    MDAAnalysisJob,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.plugin import frame_selection_payload, strict_json_payload
from polyzymd.analyses.shared.loader import parse_time_string
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.mda import ArtifactSidecarRef, MDAReplicateJobContext
    from polyzymd.analyses.secondary_structure import SecondaryStructureSettings

LOGGER = logging.getLogger(__name__)

DSSP_STATE_ENCODING: dict[str, int] = {"C": 0, "H": 1, "E": 2}
DSSP_STATE_LABELS: dict[str, str] = {"C": "coil", "H": "helix", "E": "strand"}
DSSP_MATRIX_SIDECAR = "dssp_state_matrix.npz"
HELIX_FRACTION_METRIC = "helix_fraction"
SS_METRIC_METADATA: dict[str, Any] = {
    "higher_is_better": True,
    "direction_labels": ("destabilizing", "unchanged", "stabilizing"),
    "description": "Mean fraction of DSSP states assigned as helix",
}


def build_secondary_structure_jobs(
    ctx: MDAReplicateJobContext,
    settings: SecondaryStructureSettings,
) -> list[MDAAnalysisJob]:
    """Build the single DSSP categorical-state job for one replicate.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        Framework-provided MDAnalysis replicate context.
    settings : SecondaryStructureSettings
        User-facing secondary-structure settings.

    Returns
    -------
    list of MDAAnalysisJob
        One job backed by an AnalysisBase-compatible DSSP adapter.
    """

    selection = _protein_chain_selection(settings.chain_id)
    return [
        MDAAnalysisJob(
            name="dssp",
            analysis=build_dssp_analysis(
                universe=ctx.universe,
                chain_id=settings.chain_id,
                timestep_ps=ctx.frame_selection.timestep_ps,
            ),
            frame_selection=ctx.frame_selection,
            backend_policy=ctx.backend_policy,
            universe_policy=MDAUniversePolicy(
                condition_label=ctx.replicate_context.condition.label,
                replicate=ctx.replicate,
                provenance=ctx.universe_policy.provenance,
                metadata={
                    **ctx.universe_policy.metadata,
                    "secondary_structure_selection": selection,
                    "dssp_state_encoding": dict(DSSP_STATE_ENCODING),
                },
            ),
        )
    ]


def build_dssp_analysis(*, universe: Any, chain_id: str, timestep_ps: float | None) -> Any:
    """Build an AnalysisBase-compatible DSSP adapter.

    The adapter lets MDAnalysis own frame iteration, including explicit frame
    selectors. MDTraj is used only at conclusion time as the DSSP kernel.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe for one replicate.
    chain_id : str
        Protein chain ID following the PolyzyMD chain convention.
    timestep_ps : float or None
        Trajectory timestep in picoseconds when available.

    Returns
    -------
    Any
        AnalysisBase-compatible object with a ``results`` namespace.
    """

    from MDAnalysis.analysis.base import AnalysisBase

    selection = _protein_chain_selection(chain_id)

    class DSSPAnalysis(AnalysisBase):
        """Collect protein-chain coordinates and compute DSSP states."""

        def __init__(self) -> None:
            atom_group = universe.select_atoms(selection)
            if len(atom_group) == 0:
                raise ValueError(
                    f"Secondary-structure selection matched no atoms: {selection!r}. "
                    "Check the configured chain_id."
                )
            super().__init__(universe.trajectory)
            self._atom_group = atom_group
            self._raw_timestep_ps = float(timestep_ps) if timestep_ps is not None else None

        def _prepare(self) -> None:
            """Initialize per-frame coordinate storage."""

            self.results.positions_angstrom = []
            self.results.selection = selection

        def _single_frame(self) -> None:
            """Store selected atom coordinates for the current trajectory frame."""

            self.results.positions_angstrom.append(
                np.asarray(self._atom_group.positions, dtype=np.float32).copy()
            )

        def _conclude(self) -> None:
            """Run DSSP over selected frames and expose categorical arrays."""

            positions = np.asarray(self.results.positions_angstrom, dtype=np.float32)
            if positions.size == 0:
                raise ValueError("DSSP analysis selected no frames")
            topology, residue_identity = _build_mdtraj_topology(self._atom_group)
            state_matrix = _compute_dssp_state_matrix(positions, topology)
            frame_indices = _analysis_frame_indices(self, state_matrix.shape[0])
            time_ps = _analysis_time_ps(self, state_matrix.shape[0], self._raw_timestep_ps)

            self.results.state_matrix = state_matrix
            self.results.frame_indices = frame_indices
            self.results.time_ps = time_ps
            self.results.residue_ids = residue_identity["residue_ids"]
            self.results.residue_names = residue_identity["residue_names"]
            self.results.residue_indices = residue_identity["residue_indices"]
            self.results.identity_keys = residue_identity["identity_keys"]
            self.results.selection = selection
            self.results.n_frames = int(state_matrix.shape[0])
            self.results.n_residues = int(state_matrix.shape[1])

    return DSSPAnalysis()


class SecondaryStructureArtifactCollector:
    """Collect completed DSSP jobs into a canonical replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Map a completed DSSP job to a sidecar-backed artifact.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context for one replicate.
        completed_jobs : sequence of MDAJobResult
            Completed MDAnalysis-compatible DSSP jobs.

        Returns
        -------
        ReplicateArtifact
            Canonical replicate artifact with a DSSP matrix sidecar.
        """

        if len(completed_jobs) != 1:
            raise ValueError(
                "Secondary-structure collection expects exactly one DSSP job, "
                f"got {len(completed_jobs)}"
            )
        job = completed_jobs[0]
        state_matrix = _validated_state_matrix(job.results.state_matrix)
        residue_ids = _int_list(job.results.residue_ids, field="residue_ids")
        residue_indices = _int_list(job.results.residue_indices, field="residue_indices")
        residue_names = _str_list(job.results.residue_names, field="residue_names")
        identity_keys = _str_list(job.results.identity_keys, field="identity_keys")
        frame_indices = _int_list(job.results.frame_indices, field="frame_indices")
        time_ps = np.asarray(job.results.time_ps, dtype=np.float64)
        _validate_matrix_identity(
            state_matrix=state_matrix,
            residue_ids=residue_ids,
            residue_names=residue_names,
            residue_indices=residue_indices,
            identity_keys=identity_keys,
            frame_indices=frame_indices,
            time_ps=time_ps,
        )

        sidecar = ctx.artifact_store.write_npz_sidecar(
            DSSP_MATRIX_SIDECAR,
            state_matrix=state_matrix.astype(np.int8),
            residue_ids=np.asarray(residue_ids, dtype=np.int64),
            residue_names=np.asarray(residue_names, dtype="U8"),
            residue_indices=np.asarray(residue_indices, dtype=np.int64),
            identity_keys=np.asarray(identity_keys, dtype="U64"),
            frame_indices=np.asarray(frame_indices, dtype=np.int64),
            time_ps=time_ps.astype(np.float64),
            metadata={
                "state_encoding": dict(DSSP_STATE_ENCODING),
                "matrix_shape": list(state_matrix.shape),
            },
        )
        fractions = _state_fraction_summary(state_matrix)
        config_hash = compute_config_hash(ctx.replicate_context.sim_config)
        eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
        metric_value = fractions["overall_helix_fraction"]

        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "selection": str(job.results.selection),
                "state_encoding": dict(DSSP_STATE_ENCODING),
                "state_labels": dict(DSSP_STATE_LABELS),
                "residue_ids": residue_ids,
                "residue_names": residue_names,
                "residue_indices": residue_indices,
                "identity_keys": identity_keys,
                "n_frames": int(state_matrix.shape[0]),
                "n_residues": int(state_matrix.shape[1]),
                "sidecar_key": DSSP_MATRIX_SIDECAR,
                "metrics": {HELIX_FRACTION_METRIC: metric_value},
                "replicate_metrics": {HELIX_FRACTION_METRIC: metric_value},
                **fractions,
            },
            sidecars=[sidecar],
            provenance={
                "source": "mda_secondary_structure_dssp_job",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
            },
            metadata={
                "result_kind": "secondary_structure_mda_replicate",
                "settings_fingerprint": ctx.settings_fingerprint,
                "config_hash": config_hash,
                "polyzymd_version": get_polyzymd_version(),
                "mdanalysis_version": mdanalysis_version(),
                "mdtraj_version": mdtraj_version(),
                "equilibration_time": eq_value,
                "equilibration_unit": eq_unit,
                "selection_string": str(job.results.selection),
                "dssp_state_encoding": dict(DSSP_STATE_ENCODING),
            },
            warnings=list(ctx.warnings),
        )


def aggregate_secondary_structure_artifacts(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: SecondaryStructureSettings,
    equilibration: str,
    output_dir: Path,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> ConditionArtifact:
    """Aggregate secondary-structure replicate artifacts into a condition artifact.

    Parameters
    ----------
    condition_label : str
        Label for the condition being aggregated.
    replicates : sequence of int
        Replicate IDs represented in ``artifacts``.
    settings : SecondaryStructureSettings
        Active secondary-structure settings.
    equilibration : str
        Equilibration string from the framework context.
    output_dir : Path
        Aggregated output directory.
    artifacts : sequence of ReplicateArtifact
        Per-replicate secondary-structure artifacts.
    settings_fingerprint : str
        Active settings fingerprint.

    Returns
    -------
    ConditionArtifact
        Aggregated condition artifact.
    """

    del settings
    if not artifacts:
        raise ValueError(
            f"Secondary-structure aggregation for condition '{condition_label}' requires at "
            "least one replicate artifact. No replicate inputs were provided."
        )
    ordered_artifacts = _validate_and_order_artifacts(
        condition_label=condition_label,
        expected_replicates=replicates,
        settings_fingerprint=settings_fingerprint,
        artifacts=artifacts,
    )
    matrices = _load_validated_matrices(
        artifacts=ordered_artifacts,
        analysis_dir=output_dir.parent,
    )
    first_artifact, first_matrix = matrices[0]
    residue_ids = _int_list(first_artifact.payload["residue_ids"], field="residue_ids")
    residue_names = _str_list(first_artifact.payload["residue_names"], field="residue_names")
    residue_indices = _int_list(first_artifact.payload["residue_indices"], field="residue_indices")
    identity_keys = _str_list(first_artifact.payload["identity_keys"], field="identity_keys")

    per_replicate: dict[str, list[float]] = {"coil": [], "helix": [], "strand": []}
    persistence: dict[str, list[NDArray[np.float64]]] = {
        "coil": [],
        "helix": [],
        "strand": [],
    }
    for artifact, matrix in matrices:
        _validate_matching_residue_identity(first_artifact, artifact)
        summary = _state_fraction_summary(matrix)
        per_replicate["coil"].append(summary["overall_coil_fraction"])
        per_replicate["helix"].append(summary["overall_helix_fraction"])
        per_replicate["strand"].append(summary["overall_strand_fraction"])
        persistence["coil"].append((matrix == DSSP_STATE_ENCODING["C"]).mean(axis=0))
        persistence["helix"].append((matrix == DSSP_STATE_ENCODING["H"]).mean(axis=0))
        persistence["strand"].append((matrix == DSSP_STATE_ENCODING["E"]).mean(axis=0))

    persistence_payload = _aggregate_persistence(persistence)
    n_reps = len(ordered_artifacts)
    helix_stats = compute_sem(per_replicate["helix"])
    strand_stats = compute_sem(per_replicate["strand"])
    coil_stats = compute_sem(per_replicate["coil"])
    metric = _metric_summary(HELIX_FRACTION_METRIC, per_replicate["helix"], helix_stats)
    metric.update(SS_METRIC_METADATA)
    replicate_metrics = {
        str(artifact.replicate): {HELIX_FRACTION_METRIC: value}
        for artifact, value in zip(ordered_artifacts, per_replicate["helix"], strict=True)
    }
    eq_value, eq_unit = parse_time_string(equilibration)
    config_hash = str(first_artifact.metadata.get("config_hash", "unknown"))
    artifact = ConditionArtifact(
        analysis_name="secondary_structure",
        condition_label=condition_label,
        replicates=[int(replicate) for replicate in replicates],
        payload={
            "state_encoding": dict(DSSP_STATE_ENCODING),
            "state_labels": dict(DSSP_STATE_LABELS),
            "residue_ids": residue_ids,
            "residue_names": residue_names,
            "residue_indices": residue_indices,
            "identity_keys": identity_keys,
            "n_replicates": n_reps,
            "n_residues": int(first_matrix.shape[1]),
            "mean_overall_helix": helix_stats.mean,
            "sem_overall_helix": helix_stats.sem,
            "mean_overall_strand": strand_stats.mean,
            "sem_overall_strand": strand_stats.sem,
            "mean_overall_coil": coil_stats.mean,
            "sem_overall_coil": coil_stats.sem,
            "per_replicate_helix": per_replicate["helix"],
            "per_replicate_strand": per_replicate["strand"],
            "per_replicate_coil": per_replicate["coil"],
            "metrics": {HELIX_FRACTION_METRIC: metric},
            "replicate_metrics": replicate_metrics,
            "metric_metadata": {HELIX_FRACTION_METRIC: SS_METRIC_METADATA},
            **persistence_payload,
        },
        provenance={
            "source": "secondary_structure_replicate_artifact_aggregation",
            "frame_selection": first_artifact.provenance.get("frame_selection"),
            "state_encoding": dict(DSSP_STATE_ENCODING),
        },
        metadata={
            "result_kind": "secondary_structure_mda_condition",
            "settings_fingerprint": settings_fingerprint,
            "config_hash": config_hash,
            "polyzymd_version": get_polyzymd_version(),
            "mdanalysis_version": mdanalysis_version(),
            "mdtraj_version": mdtraj_version(),
            "equilibration_time": eq_value,
            "equilibration_unit": eq_unit,
            "selection_string": first_artifact.payload.get("selection"),
            "dssp_state_encoding": dict(DSSP_STATE_ENCODING),
            "n_replicates": n_reps,
            "source_result_files": [
                str(output_dir.parent / f"run_{artifact.replicate}" / "result.json")
                for artifact in ordered_artifacts
            ],
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
    return artifact


def load_replicate_matrix(artifact: ReplicateArtifact, run_dir: Path) -> NDArray[np.int8]:
    """Load and validate the DSSP matrix sidecar for one replicate.

    Parameters
    ----------
    artifact : ReplicateArtifact
        Replicate artifact that references the sidecar.
    run_dir : Path
        Replicate artifact-store root.

    Returns
    -------
    ndarray of int8
        Validated DSSP state matrix.
    """

    sidecar = _dssp_sidecar(artifact)
    store = ArtifactStore(run_dir)
    sidecar_path = store.validate_sidecar(sidecar)
    try:
        with np.load(sidecar_path) as data:
            matrix = _validated_state_matrix(data["state_matrix"])
            residue_ids = _int_list(data["residue_ids"].tolist(), field="residue_ids")
            residue_names = _str_list(data["residue_names"].tolist(), field="residue_names")
            residue_indices = _int_list(data["residue_indices"].tolist(), field="residue_indices")
            identity_keys = _str_list(data["identity_keys"].tolist(), field="identity_keys")
            frame_indices = _int_list(data["frame_indices"].tolist(), field="frame_indices")
            time_ps = np.asarray(data["time_ps"], dtype=np.float64)
    except (OSError, KeyError, ValueError) as exc:
        raise ValueError(
            f"Failed to load DSSP state-matrix sidecar for replicate {artifact.replicate}: {exc}"
        ) from exc
    _validate_matrix_identity(
        state_matrix=matrix,
        residue_ids=residue_ids,
        residue_names=residue_names,
        residue_indices=residue_indices,
        identity_keys=identity_keys,
        frame_indices=frame_indices,
        time_ps=time_ps,
    )
    _validate_artifact_sidecar_identity(artifact, residue_ids, residue_names, residue_indices)
    return matrix


def _protein_chain_selection(chain_id: str) -> str:
    """Return the MDAnalysis selection for a protein chain."""

    return f"protein and chainid {chain_id}"


def _build_mdtraj_topology(atom_group: Any) -> tuple[Any, dict[str, list[Any]]]:
    """Build an MDTraj topology from the selected MDAnalysis atoms."""

    import mdtraj as md

    topology = md.Topology()
    chain = topology.add_chain()
    residue_ids: list[int] = []
    residue_names: list[str] = []
    residue_indices: list[int] = []
    identity_keys: list[str] = []

    for residue in atom_group.residues:
        residue_id = int(
            getattr(residue, "resid", getattr(residue, "resnum", len(residue_ids) + 1))
        )
        residue_name = str(getattr(residue, "resname", getattr(residue, "name", "UNK"))).upper()
        residue_index = int(getattr(residue, "ix", getattr(residue, "resindex", len(residue_ids))))
        md_residue = topology.add_residue(residue_name, chain, resSeq=residue_id)
        residue_ids.append(residue_id)
        residue_names.append(residue_name)
        residue_indices.append(residue_index)
        identity_keys.append(f"{residue_index}:{residue_id}:{residue_name}")
        for atom in residue.atoms:
            topology.add_atom(str(atom.name), _mdtraj_element(md, atom), md_residue)

    return topology, {
        "residue_ids": residue_ids,
        "residue_names": residue_names,
        "residue_indices": residue_indices,
        "identity_keys": identity_keys,
    }


def _mdtraj_element(md: Any, atom: Any) -> Any:
    """Return a best-effort MDTraj element for an MDAnalysis atom."""

    symbol = str(getattr(atom, "element", "") or "").strip()
    if not symbol:
        atom_name = str(getattr(atom, "name", "")).strip()
        symbol = "".join(char for char in atom_name if char.isalpha())[:1]
    symbol = symbol.capitalize() if symbol else "C"
    try:
        return md.element.get_by_symbol(symbol)
    except KeyError:
        return md.element.carbon


def _compute_dssp_state_matrix(
    positions_angstrom: NDArray[np.float32], topology: Any
) -> NDArray[np.int8]:
    """Compute and encode simplified DSSP states with MDTraj."""

    import mdtraj as md

    trajectory = md.Trajectory(
        xyz=np.asarray(positions_angstrom, dtype=np.float32) / 10.0, top=topology
    )
    raw_dssp = md.compute_dssp(trajectory, simplified=True)
    return encode_dssp_matrix(raw_dssp)


def encode_dssp_matrix(dssp_raw: NDArray[Any]) -> NDArray[np.int8]:
    """Encode simplified DSSP characters as categorical integer states.

    Parameters
    ----------
    dssp_raw : ndarray
        Character matrix with simplified DSSP states ``C``, ``H``, and ``E``.

    Returns
    -------
    ndarray of int8
        Encoded state matrix with ``C=0``, ``H=1``, and ``E=2``.
    """

    raw = np.asarray(dssp_raw)
    if raw.ndim != 2:
        raise ValueError(f"DSSP state matrix must be 2-dimensional, got shape {raw.shape}")
    matrix = np.zeros(raw.shape, dtype=np.int8)
    for char, code in DSSP_STATE_ENCODING.items():
        matrix[raw == char] = code
    return matrix


def _analysis_frame_indices(analysis: Any, n_frames: int) -> NDArray[np.int64]:
    """Return frame indices selected by AnalysisBase."""

    frames = getattr(analysis, "frames", None)
    if frames is None:
        return np.arange(n_frames, dtype=np.int64)
    frame_array = np.asarray(frames, dtype=np.int64)
    if frame_array.shape[0] != n_frames:
        return np.arange(n_frames, dtype=np.int64)
    return frame_array


def _analysis_time_ps(
    analysis: Any,
    n_frames: int,
    raw_timestep_ps: float | None,
) -> NDArray[np.float64]:
    """Return selected frame times in picoseconds."""

    times = getattr(analysis, "times", None)
    if times is not None:
        time_array = np.asarray(times, dtype=np.float64)
        if time_array.shape[0] == n_frames:
            return time_array
    timestep = (
        raw_timestep_ps if raw_timestep_ps is not None and math.isfinite(raw_timestep_ps) else 1.0
    )
    return _analysis_frame_indices(analysis, n_frames).astype(np.float64) * float(timestep)


def _validated_state_matrix(value: Any) -> NDArray[np.int8]:
    """Validate a DSSP categorical matrix."""

    matrix = np.asarray(value)
    if matrix.ndim != 2:
        raise ValueError(f"DSSP state matrix must be 2-dimensional, got shape {matrix.shape}")
    if matrix.shape[0] < 1 or matrix.shape[1] < 1:
        raise ValueError("DSSP state matrix must contain at least one frame and one residue")
    matrix = matrix.astype(np.int8, copy=False)
    valid = set(DSSP_STATE_ENCODING.values())
    observed = set(np.unique(matrix).astype(int).tolist())
    invalid = sorted(observed - valid)
    if invalid:
        raise ValueError(f"DSSP state matrix contains invalid state codes: {invalid}")
    return matrix


def _validate_matrix_identity(
    *,
    state_matrix: NDArray[np.int8],
    residue_ids: Sequence[int],
    residue_names: Sequence[str],
    residue_indices: Sequence[int],
    identity_keys: Sequence[str],
    frame_indices: Sequence[int],
    time_ps: NDArray[np.float64],
) -> None:
    """Validate sidecar identity arrays against the matrix shape."""

    n_frames, n_residues = state_matrix.shape
    if not (
        len(residue_ids)
        == len(residue_names)
        == len(residue_indices)
        == len(identity_keys)
        == n_residues
    ):
        raise ValueError("DSSP residue identity arrays do not match matrix residue dimension")
    if len(frame_indices) != n_frames or len(time_ps) != n_frames:
        raise ValueError("DSSP frame identity arrays do not match matrix frame dimension")


def _state_fraction_summary(matrix: NDArray[np.int8]) -> dict[str, Any]:
    """Compute per-residue persistence and overall state fractions."""

    total_entries = int(matrix.size)
    return {
        "persistence_coil": (matrix == DSSP_STATE_ENCODING["C"]).mean(axis=0).tolist(),
        "persistence_helix": (matrix == DSSP_STATE_ENCODING["H"]).mean(axis=0).tolist(),
        "persistence_strand": (matrix == DSSP_STATE_ENCODING["E"]).mean(axis=0).tolist(),
        "overall_coil_fraction": float(np.sum(matrix == DSSP_STATE_ENCODING["C"])) / total_entries,
        "overall_helix_fraction": float(np.sum(matrix == DSSP_STATE_ENCODING["H"])) / total_entries,
        "overall_strand_fraction": float(np.sum(matrix == DSSP_STATE_ENCODING["E"]))
        / total_entries,
    }


def _aggregate_persistence(
    persistence: dict[str, list[NDArray[np.float64]]],
) -> dict[str, list[float]]:
    """Aggregate per-residue state persistence arrays across replicates."""

    payload: dict[str, list[float]] = {}
    for state, values in persistence.items():
        stack = np.vstack(values).astype(np.float64)
        payload[f"mean_persistence_{state}"] = stack.mean(axis=0).tolist()
        if stack.shape[0] == 1:
            sem = np.zeros(stack.shape[1], dtype=np.float64)
        else:
            sem = stack.std(axis=0, ddof=1) / math.sqrt(stack.shape[0])
        payload[f"sem_persistence_{state}"] = sem.tolist()
    return payload


def _metric_summary(name: str, values: Sequence[float], stats: Any) -> dict[str, Any]:
    """Build an aggregated metric payload."""

    value_array = np.asarray(values, dtype=np.float64)
    return {
        "name": name,
        "values": [float(value) for value in values],
        "mean": float(stats.mean),
        "sem": float(stats.sem),
        "std": float(np.std(value_array, ddof=1)) if len(value_array) > 1 else 0.0,
        "n": int(len(value_array)),
    }


def _validate_and_order_artifacts(
    *,
    condition_label: str,
    expected_replicates: Sequence[int],
    settings_fingerprint: str,
    artifacts: Sequence[ReplicateArtifact],
) -> list[ReplicateArtifact]:
    """Validate replicate artifacts and return them in requested order."""

    by_replicate: dict[int, ReplicateArtifact] = {}
    for artifact in artifacts:
        if artifact.analysis_name != "secondary_structure":
            raise ValueError(
                "Secondary-structure aggregation received an artifact for "
                f"{artifact.analysis_name!r}"
            )
        if artifact.condition_label != condition_label:
            raise ValueError(
                f"Secondary-structure artifact condition mismatch: expected {condition_label!r}, "
                f"got {artifact.condition_label!r}"
            )
        stored_fingerprint = artifact.metadata.get("settings_fingerprint")
        if stored_fingerprint != settings_fingerprint:
            raise ValueError(
                "Secondary-structure settings fingerprint mismatch for replicate "
                f"{artifact.replicate}: expected {settings_fingerprint!r}, got "
                f"{stored_fingerprint!r}. Recompute stale caches."
            )
        if artifact.replicate in by_replicate:
            raise ValueError(f"Duplicate secondary-structure replicate {artifact.replicate}")
        by_replicate[int(artifact.replicate)] = artifact

    ordered: list[ReplicateArtifact] = []
    missing: list[int] = []
    for replicate in expected_replicates:
        artifact = by_replicate.get(int(replicate))
        if artifact is None:
            missing.append(int(replicate))
            continue
        ordered.append(artifact)
    if missing:
        raise ValueError(f"Missing secondary-structure replicate artifacts: {missing}")
    return ordered


def _load_validated_matrices(
    *,
    artifacts: Sequence[ReplicateArtifact],
    analysis_dir: Path,
) -> list[tuple[ReplicateArtifact, NDArray[np.int8]]]:
    """Load sidecar matrices for all replicate artifacts."""

    matrices: list[tuple[ReplicateArtifact, NDArray[np.int8]]] = []
    for artifact in artifacts:
        run_dir = analysis_dir / f"run_{artifact.replicate}"
        matrix = load_replicate_matrix(artifact, run_dir)
        matrices.append((artifact, matrix))
    return matrices


def _validate_matching_residue_identity(
    reference: ReplicateArtifact,
    candidate: ReplicateArtifact,
) -> None:
    """Validate residue identity and order against a reference artifact."""

    for key in ("residue_ids", "residue_names", "residue_indices", "identity_keys"):
        if list(reference.payload.get(key, [])) != list(candidate.payload.get(key, [])):
            raise ValueError(
                "Secondary-structure residue identity mismatch between replicates: "
                f"replicate {reference.replicate} and {candidate.replicate} differ in {key}."
            )


def _validate_artifact_sidecar_identity(
    artifact: ReplicateArtifact,
    residue_ids: Sequence[int],
    residue_names: Sequence[str],
    residue_indices: Sequence[int],
) -> None:
    """Validate sidecar residue identity arrays against artifact JSON metadata."""

    checks = {
        "residue_ids": list(residue_ids),
        "residue_names": list(residue_names),
        "residue_indices": list(residue_indices),
    }
    for key, sidecar_value in checks.items():
        if list(artifact.payload.get(key, [])) != sidecar_value:
            raise ValueError(
                f"Secondary-structure sidecar identity mismatch for replicate {artifact.replicate}: "
                f"{key} differs from result.json metadata."
            )


def _dssp_sidecar(artifact: ReplicateArtifact) -> ArtifactSidecarRef:
    """Return the DSSP state-matrix sidecar reference from an artifact."""

    for sidecar in artifact.sidecars:
        if sidecar.path == DSSP_MATRIX_SIDECAR:
            return sidecar
    raise ValueError(
        f"Secondary-structure replicate {artifact.replicate} is missing {DSSP_MATRIX_SIDECAR}"
    )


def _int_list(values: Any, *, field: str) -> list[int]:
    """Coerce a JSON or NumPy sequence to a list of integers."""

    try:
        return [int(value) for value in values]
    except TypeError as exc:
        raise ValueError(f"{field} must be a sequence of integers") from exc


def _str_list(values: Any, *, field: str) -> list[str]:
    """Coerce a JSON or NumPy sequence to a list of strings."""

    try:
        return [str(value) for value in values]
    except TypeError as exc:
        raise ValueError(f"{field} must be a sequence of strings") from exc


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
