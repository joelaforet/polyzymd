"""Shared SASA computation utilities.

This module provides reusable helpers for SASA analyses that need
selection validation, atom-level and residue-level SASA traces, and
raw artifact persistence.
"""

from __future__ import annotations

import json
import logging
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses.shared.diagnostics import warn_if_multi_chain_selection

LOGGER = logging.getLogger(__name__)

NM2_TO_A2 = 100.0


@dataclass(frozen=True)
class SASAComputationResult:
    """Container for raw SASA computation outputs.

    Parameters
    ----------
    atom_sasa_a2 : NDArray[np.float64]
        Per-atom SASA trace in Å², shape ``(n_frames, n_target_atoms)``.
    residue_sasa_a2 : NDArray[np.float64]
        Per-residue SASA trace in Å², shape ``(n_frames, n_target_residues)``.
    total_sasa_a2 : NDArray[np.float64]
        Total target SASA in Å², shape ``(n_frames,)``.
    frames : NDArray[np.int64]
        0-indexed frame indices used for analysis.
    time_ns : NDArray[np.float64]
        Time axis in ns corresponding to ``frames``.
    target_atom_indices : NDArray[np.int64]
        Universe-global atom indices for target atoms.
    context_atom_indices : NDArray[np.int64]
        Universe-global atom indices for context atoms.
    residue_keys : list[str]
        Residue identity keys in ``chainID:resid:resname`` format.
    residue_chainids : list[str]
        Chain IDs for each residue key.
    residue_resids : list[int]
        Residue IDs for each residue key.
    residue_resnames : list[str]
        Residue names for each residue key.
    """

    atom_sasa_a2: NDArray[np.float64]
    residue_sasa_a2: NDArray[np.float64]
    total_sasa_a2: NDArray[np.float64]
    frames: NDArray[np.int64]
    time_ns: NDArray[np.float64]
    target_atom_indices: NDArray[np.int64]
    context_atom_indices: NDArray[np.int64]
    residue_keys: list[str]
    residue_chainids: list[str]
    residue_resids: list[int]
    residue_resnames: list[str]


def resolve_selection_indices(
    universe: Any,
    selection: str,
    *,
    role: str,
    run_label: str,
) -> tuple[Any, NDArray[np.int64]]:
    """Resolve an MDAnalysis selection to atom indices.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe.
    selection : str
        MDAnalysis selection string.
    role : str
        Selection role used in warning context.
    run_label : str
        Human-readable run label.

    Returns
    -------
    tuple[Any, NDArray[np.int64]]
        Selected atom group and its universe-global atom indices.
    """
    atoms = universe.select_atoms(selection)
    context = f"for SASA {role} selection in run '{run_label}'"
    if role == "context":
        context += " (multi-chain context selections are often intentional for SASA)"
    warn_if_multi_chain_selection(atoms, selection, context=context)
    return atoms, np.asarray(atoms.indices, dtype=np.int64)


def validate_target_subset(
    target_indices: NDArray[np.int64],
    context_indices: NDArray[np.int64],
    *,
    run_label: str,
    target_selection: str,
    context_selection: str,
) -> None:
    """Validate target-selection atoms are a subset of context-selection atoms.

    Raises
    ------
    ValueError
        If any target atom is absent from the context selection.
    """
    if target_indices.size == 0:
        return

    context_set = set(context_indices.tolist())
    missing = [idx for idx in target_indices.tolist() if idx not in context_set]
    if not missing:
        return

    raise ValueError(
        "SASA run '{label}' requires target_selection atoms to be a subset of "
        "context_selection atoms. Missing {count} target atoms in context. "
        "target_selection={target!r}, context_selection={context!r}".format(
            label=run_label,
            count=len(missing),
            target=target_selection,
            context=context_selection,
        )
    )


def compute_sasa(
    universe: Any,
    *,
    run_label: str,
    target_selection: str,
    context_selection: str,
    probe_radius_nm: float,
    n_sphere_points: int,
    start_frame: int,
    stop_frame: int,
    timestep_ps: float,
    chunk_size: int = 100,
    stride: int = 1,
) -> SASAComputationResult:
    """Compute target SASA over a selected context.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe.
    run_label : str
        Human-readable run label.
    target_selection : str
        Selection of atoms whose SASA is reported.
    context_selection : str
        Selection of atoms considered during SASA computation.
    probe_radius_nm : float
        Probe radius in nm.
    n_sphere_points : int
        Number of sphere points for Shrake-Rupley.
    start_frame : int
        First frame index (inclusive).
    stop_frame : int
        Last frame index (exclusive).
    timestep_ps : float
        Timestep in ps.
    chunk_size : int, optional
        Number of analyzed frames processed per Shrake-Rupley chunk,
        by default 100.
    stride : int, optional
        Frame stride applied before chunking (1 = analyze every frame),
        by default 1.

    Returns
    -------
    SASAComputationResult
        Raw atom-level, residue-level, and total SASA traces in Å².

    Raises
    ------
    ValueError
        If frame bounds or stride/chunk parameters are invalid.
    """
    import mdtraj as md

    if chunk_size <= 0:
        raise ValueError("chunk_size must be >= 1")
    if stride <= 0:
        raise ValueError("stride must be >= 1")

    trajectory_length = len(universe.trajectory)
    if start_frame < 0:
        raise ValueError(f"start_frame must be >= 0, got {start_frame}")
    if stop_frame > trajectory_length:
        raise ValueError(
            f"stop_frame must be <= trajectory length ({trajectory_length}), got {stop_frame}"
        )
    if start_frame >= stop_frame:
        raise ValueError(
            "start_frame must be < stop_frame, got "
            f"start_frame={start_frame}, stop_frame={stop_frame}"
        )

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
    frames = np.arange(start_frame, stop_frame, dtype=np.int64)
    analyzed_frames = frames[::stride]
    n_frames = int(analyzed_frames.size)
    if context_indices.size == 0:
        LOGGER.warning(
            "Run '%s' context selection matched zero atoms (%r); returning NaN SASA metrics",
            run_label,
            context_selection,
        )
    if target_indices.size == 0:
        LOGGER.warning(
            "Run '%s' target selection matched zero atoms (%r); returning NaN SASA metrics",
            run_label,
            target_selection,
        )

    if target_indices.size == 0 or context_indices.size == 0 or n_frames == 0:
        time_ns = (analyzed_frames.astype(np.float64) * timestep_ps) / 1000.0
        return SASAComputationResult(
            atom_sasa_a2=np.empty((n_frames, target_indices.size), dtype=np.float64),
            residue_sasa_a2=np.empty((n_frames, 0), dtype=np.float64),
            total_sasa_a2=np.full(n_frames, np.nan, dtype=np.float64),
            frames=analyzed_frames,
            time_ns=time_ns,
            target_atom_indices=target_indices,
            context_atom_indices=context_indices,
            residue_keys=[],
            residue_chainids=[],
            residue_resids=[],
            residue_resnames=[],
        )

    validate_target_subset(
        target_indices,
        context_indices,
        run_label=run_label,
        target_selection=target_selection,
        context_selection=context_selection,
    )

    context_index_to_local = {int(idx): i for i, idx in enumerate(context_indices.tolist())}
    target_local_indices = np.asarray(
        [context_index_to_local[int(idx)] for idx in target_indices.tolist()],
        dtype=np.int64,
    )

    residue_to_indices: dict[tuple[str, int, str], list[int]] = {}
    for atom_local, atom in zip(target_local_indices.tolist(), target_atoms):
        key = (str(atom.chainID), int(atom.resid), str(atom.resname))
        residue_to_indices.setdefault(key, []).append(int(atom_local))

    residue_items = list(residue_to_indices.items())
    residue_keys = [f"{chain}:{resid}:{resname}" for (chain, resid, resname), _ in residue_items]
    residue_chainids = [chain for (chain, _, _), _ in residue_items]
    residue_resids = [resid for (_, resid, _), _ in residue_items]
    residue_resnames = [resname for (_, _, resname), _ in residue_items]

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=True) as tmp_pdb:
        context_atoms.write(tmp_pdb.name)
        template = md.load(tmp_pdb.name)

    atom_sasa_target_a2 = np.empty((n_frames, target_local_indices.size), dtype=np.float64)
    residue_sasa_a2 = np.empty((n_frames, len(residue_items)), dtype=np.float64)
    total_sasa_a2 = np.empty(n_frames, dtype=np.float64)

    n_chunks = (n_frames + chunk_size - 1) // chunk_size
    for chunk_idx, chunk_start in enumerate(range(0, n_frames, chunk_size)):
        chunk_end = min(chunk_start + chunk_size, n_frames)
        LOGGER.info(
            "Computing SASA chunk %d/%d (frames %d-%d)...",
            chunk_idx + 1,
            n_chunks,
            int(analyzed_frames[chunk_start]),
            int(analyzed_frames[chunk_end - 1]),
        )

        chunk_frames = analyzed_frames[chunk_start:chunk_end]
        xyz_nm = np.empty((chunk_frames.size, len(context_atoms), 3), dtype=np.float32)
        for out_idx, frame_idx in enumerate(chunk_frames.tolist()):
            universe.trajectory[frame_idx]
            xyz_nm[out_idx] = context_atoms.positions.astype(np.float32) / 10.0

        mdtraj_traj = md.Trajectory(xyz=xyz_nm, topology=template.topology)
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

    time_ns = (analyzed_frames.astype(np.float64) * timestep_ps) / 1000.0

    return SASAComputationResult(
        atom_sasa_a2=atom_sasa_target_a2,
        residue_sasa_a2=residue_sasa_a2,
        total_sasa_a2=total_sasa_a2,
        frames=analyzed_frames,
        time_ns=time_ns,
        target_atom_indices=target_indices,
        context_atom_indices=context_indices,
        residue_keys=residue_keys,
        residue_chainids=residue_chainids,
        residue_resids=residue_resids,
        residue_resnames=residue_resnames,
    )


def save_sasa_artifacts(
    npz_path: Path,
    metadata_path: Path,
    result: SASAComputationResult,
    *,
    run_label: str,
    target_selection: str,
    context_selection: str,
    probe_radius_nm: float,
    n_sphere_points: int,
    equilibration: str,
) -> None:
    """Save raw SASA arrays plus JSON sidecar metadata.

    Parameters
    ----------
    npz_path : Path
        Output path for compressed NumPy archive.
    metadata_path : Path
        Output path for JSON metadata sidecar.
    result : SASAComputationResult
        Raw SASA arrays to persist.
    run_label : str
        Run label.
    target_selection : str
        Target selection string.
    context_selection : str
        Context selection string.
    probe_radius_nm : float
        Probe radius in nm.
    n_sphere_points : int
        Number of sphere points.
    equilibration : str
        Equilibration string.
    """
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        npz_path,
        atom_sasa_a2=result.atom_sasa_a2,
        residue_sasa_a2=result.residue_sasa_a2,
        total_sasa_a2=result.total_sasa_a2,
        frames=result.frames,
        time_ns=result.time_ns,
        target_atom_indices=result.target_atom_indices,
        context_atom_indices=result.context_atom_indices,
        residue_keys=np.asarray(result.residue_keys, dtype=str),
        residue_chainids=np.asarray(result.residue_chainids, dtype=str),
        residue_resids=np.asarray(result.residue_resids, dtype=np.int64),
        residue_resnames=np.asarray(result.residue_resnames, dtype=str),
    )

    metadata = build_sasa_artifact_metadata(
        result,
        run_label=run_label,
        target_selection=target_selection,
        context_selection=context_selection,
        probe_radius_nm=probe_radius_nm,
        n_sphere_points=n_sphere_points,
        equilibration=equilibration,
    )
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def load_sasa_artifacts(
    npz_path: Path,
    metadata_path: Path,
) -> tuple[SASAComputationResult, dict[str, Any]]:
    """Load SASA raw arrays and metadata from disk.

    Parameters
    ----------
    npz_path : Path
        Path to compressed NumPy archive.
    metadata_path : Path
        Path to JSON metadata sidecar.

    Returns
    -------
    tuple[SASAComputationResult, dict[str, Any]]
        Reconstructed raw result and parsed metadata dictionary.
    """
    with np.load(npz_path) as payload:
        result = SASAComputationResult(
            atom_sasa_a2=np.asarray(payload["atom_sasa_a2"], dtype=np.float64),
            residue_sasa_a2=np.asarray(payload["residue_sasa_a2"], dtype=np.float64),
            total_sasa_a2=np.asarray(payload["total_sasa_a2"], dtype=np.float64),
            frames=np.asarray(payload["frames"], dtype=np.int64),
            time_ns=np.asarray(payload["time_ns"], dtype=np.float64),
            target_atom_indices=(
                np.asarray(payload["target_atom_indices"], dtype=np.int64)
                if "target_atom_indices" in payload
                else np.empty((0,), dtype=np.int64)
            ),
            context_atom_indices=(
                np.asarray(payload["context_atom_indices"], dtype=np.int64)
                if "context_atom_indices" in payload
                else np.empty((0,), dtype=np.int64)
            ),
            residue_keys=[str(v) for v in payload["residue_keys"].tolist()],
            residue_chainids=[str(v) for v in payload["residue_chainids"].tolist()],
            residue_resids=[int(v) for v in payload["residue_resids"].tolist()],
            residue_resnames=[str(v) for v in payload["residue_resnames"].tolist()],
        )

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    return result, metadata


SASA_ARTIFACT_SCHEMA_NAME = "polyzymd.sasa_artifact"
SASA_ARTIFACT_SCHEMA_VERSION = 1
SASA_ARTIFACT_COMPATIBILITY_VERSION = 1
A2_TO_NM2 = 0.01
SASA_COMPAT_PROBE_RADIUS_ABS_TOL = 1e-6


@dataclass(frozen=True)
class SASAArtifactContract:
    """Canonical contract metadata for a persisted SASA artifact.

    Parameters
    ----------
    schema_name : str
        Artifact schema identifier.
    schema_version : int
        Artifact schema version.
    compatibility_version : int
        Compatibility hash payload version.
    producer : str
        Producer module name.
    engine : str
        SASA engine identifier.
    mode : str
        SASA mode. Canonical value is ``"atom"``.
    units : str
        SASA units. Canonical value is ``"A^2"``.
    run_label : str
        Human-readable run label.
    target_selection : str
        Target atom selection used for SASA reporting.
    context_selection : str
        Context atom selection used during SASA computation.
    probe_radius_nm : float
        Probe radius in nm.
    n_sphere_points : int
        Number of Shrake-Rupley sphere points.
    equilibration : str
        Equilibration descriptor.
    compatibility_hash : str
        Deterministic compatibility hash.
    """

    schema_name: str
    schema_version: int
    compatibility_version: int
    producer: str
    engine: str
    mode: str
    units: str
    run_label: str
    target_selection: str
    context_selection: str
    probe_radius_nm: float
    n_sphere_points: int
    equilibration: str
    compatibility_hash: str


@dataclass(frozen=True)
class SASAArtifactCompatibilityQuery:
    """Compatibility query metadata for sibling artifact lookup.

    Parameters
    ----------
    probe_radius_nm : float
        Probe radius expected by the consumer.
    n_sphere_points : int
        Sphere point count expected by the consumer.
    equilibration : str
        Equilibration label expected by the consumer.
    selection : str | None, optional
        Target selection string for advisory hash comparison.
    context_selection : str | None, optional
        Context selection string for advisory hash comparison.
    """

    probe_radius_nm: float
    n_sphere_points: int
    equilibration: str
    selection: str | None = None
    context_selection: str | None = None


@dataclass(frozen=True)
class SASAArtifactCompatibility:
    """Compatibility decision for a single SASA artifact metadata payload.

    Parameters
    ----------
    is_compatible : bool
        True when metadata-level required fields are compatible.
    is_legacy : bool
        True when schema version fields are absent.
    schema_version : int | None
        Parsed schema version if present.
    selection_hash_matches : bool | None
        Advisory selection hash result when query selections are supplied.
    mismatched_fields : tuple[str, ...]
        Names of fields that made the artifact incompatible.
    """

    is_compatible: bool
    is_legacy: bool
    schema_version: int | None
    selection_hash_matches: bool | None
    mismatched_fields: tuple[str, ...]


@dataclass(frozen=True)
class SASASiblingArtifactMatch:
    """A sibling SASA artifact candidate and its compatibility outcome.

    Parameters
    ----------
    sibling_analysis_dir : Path
        Directory containing sibling SASA artifacts for the replicate.
    npz_path : Path
        Path to the SASA NPZ payload.
    metadata_path : Path
        Path to the JSON metadata sidecar.
    metadata : dict[str, Any]
        Parsed metadata dictionary.
    compatibility : SASAArtifactCompatibility
        Compatibility decision for this artifact.
    """

    sibling_analysis_dir: Path
    npz_path: Path
    metadata_path: Path
    metadata: dict[str, Any]
    compatibility: SASAArtifactCompatibility


def compute_sasa_artifact_compatibility_hash(
    *,
    probe_radius_nm: float,
    n_sphere_points: int,
    selection: str,
    equilibration: str,
    context_selection: str | None = None,
) -> str:
    """Compute deterministic compatibility hash for SASA artifacts.

    Parameters
    ----------
    probe_radius_nm : float
        Probe radius in nm.
    n_sphere_points : int
        Number of Shrake-Rupley sphere points.
    selection : str
        Target selection string.
    equilibration : str
        Equilibration label.
    context_selection : str | None, optional
        Context selection string. Defaults to ``selection`` when omitted.

    Returns
    -------
    str
        First 16 characters of the SHA-256 compatibility digest.
    """
    import hashlib

    normalized_selection = selection.strip()
    normalized_context = (context_selection or selection).strip()
    payload = {
        "compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
        "mode": "atom",
        "units": "A^2",
        "selection": normalized_selection,
        "context_selection": normalized_context,
        "probe_radius_nm": round(float(probe_radius_nm), 6),
        "n_sphere_points": int(n_sphere_points),
        "equilibration": equilibration.strip(),
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:16]


def build_sasa_artifact_contract(
    *,
    run_label: str,
    target_selection: str,
    context_selection: str,
    probe_radius_nm: float,
    n_sphere_points: int,
    equilibration: str,
) -> SASAArtifactContract:
    """Build canonical SASA artifact contract metadata.

    Parameters
    ----------
    run_label : str
        Human-readable run label.
    target_selection : str
        Target atom selection used for SASA reporting.
    context_selection : str
        Context atom selection used during SASA computation.
    probe_radius_nm : float
        Probe radius in nm.
    n_sphere_points : int
        Number of Shrake-Rupley sphere points.
    equilibration : str
        Equilibration descriptor.

    Returns
    -------
    SASAArtifactContract
        Fully populated canonical metadata contract.
    """
    compatibility_hash = compute_sasa_artifact_compatibility_hash(
        probe_radius_nm=probe_radius_nm,
        n_sphere_points=n_sphere_points,
        selection=target_selection,
        context_selection=context_selection,
        equilibration=equilibration,
    )
    return SASAArtifactContract(
        schema_name=SASA_ARTIFACT_SCHEMA_NAME,
        schema_version=SASA_ARTIFACT_SCHEMA_VERSION,
        compatibility_version=SASA_ARTIFACT_COMPATIBILITY_VERSION,
        producer="polyzymd.analyses.shared.sasa",
        engine="mdtraj.shrake_rupley",
        mode="atom",
        units="A^2",
        run_label=run_label,
        target_selection=target_selection,
        context_selection=context_selection,
        probe_radius_nm=probe_radius_nm,
        n_sphere_points=n_sphere_points,
        equilibration=equilibration,
        compatibility_hash=compatibility_hash,
    )


def build_sasa_artifact_metadata(
    result: SASAComputationResult,
    *,
    run_label: str,
    target_selection: str,
    context_selection: str,
    probe_radius_nm: float,
    n_sphere_points: int,
    equilibration: str,
) -> dict[str, Any]:
    """Build flat SASA artifact metadata with schema versioning fields.

    Parameters
    ----------
    result : SASAComputationResult
        Raw SASA arrays and residue metadata.
    run_label : str
        Human-readable run label.
    target_selection : str
        Target atom selection.
    context_selection : str
        Context atom selection.
    probe_radius_nm : float
        Probe radius in nm.
    n_sphere_points : int
        Number of Shrake-Rupley sphere points.
    equilibration : str
        Equilibration descriptor.

    Returns
    -------
    dict[str, Any]
        Flat metadata dictionary for JSON sidecar persistence.
    """
    contract = build_sasa_artifact_contract(
        run_label=run_label,
        target_selection=target_selection,
        context_selection=context_selection,
        probe_radius_nm=probe_radius_nm,
        n_sphere_points=n_sphere_points,
        equilibration=equilibration,
    )
    return {
        "run_label": run_label,
        "target_selection": target_selection,
        "context_selection": context_selection,
        "units": "A^2",
        "probe_radius_nm": probe_radius_nm,
        "n_sphere_points": n_sphere_points,
        "equilibration": equilibration,
        "n_frames": int(result.frames.size),
        "n_target_atoms": int(result.target_atom_indices.size),
        "n_context_atoms": int(result.context_atom_indices.size),
        "n_target_residues": len(result.residue_keys),
        "residue_keys": result.residue_keys,
        "residue_chainids": result.residue_chainids,
        "residue_resids": result.residue_resids,
        "residue_resnames": result.residue_resnames,
        "artifact_schema": contract.schema_name,
        "artifact_schema_version": contract.schema_version,
        "artifact_compatibility_version": contract.compatibility_version,
        "artifact_producer": contract.producer,
        "sasa_engine": contract.engine,
        "sasa_mode": contract.mode,
        "compatibility_hash": contract.compatibility_hash,
    }


def check_sasa_artifact_compatibility(
    metadata: dict[str, Any],
    query: SASAArtifactCompatibilityQuery,
) -> SASAArtifactCompatibility:
    """Evaluate metadata-level compatibility for reusable SASA artifacts.

    Parameters
    ----------
    metadata : dict[str, Any]
        Parsed SASA metadata sidecar.
    query : SASAArtifactCompatibilityQuery
        Consumer compatibility query.

    Returns
    -------
    SASAArtifactCompatibility
        Compatibility decision containing definitive and advisory signals.
    """
    mismatched_fields: list[str] = []

    schema_name_raw = metadata.get("artifact_schema")
    if schema_name_raw is not None and str(schema_name_raw) != SASA_ARTIFACT_SCHEMA_NAME:
        mismatched_fields.append("artifact_schema")

    schema_raw = metadata.get("artifact_schema_version")
    schema_version: int | None
    if schema_raw is None:
        schema_version = None
    else:
        try:
            schema_version = int(schema_raw)
        except (TypeError, ValueError):
            schema_version = None
            mismatched_fields.append("artifact_schema_version")

    is_legacy = schema_raw is None
    if schema_version is not None and schema_version > SASA_ARTIFACT_SCHEMA_VERSION:
        mismatched_fields.append("artifact_schema_version")

    compat_version_raw = metadata.get("artifact_compatibility_version")
    if compat_version_raw is not None:
        try:
            compatibility_version = int(compat_version_raw)
        except (TypeError, ValueError):
            compatibility_version = None
            mismatched_fields.append("artifact_compatibility_version")
        if (
            compatibility_version is not None
            and compatibility_version > SASA_ARTIFACT_COMPATIBILITY_VERSION
        ):
            mismatched_fields.append("artifact_compatibility_version")

    probe_value = metadata.get("probe_radius_nm")
    try:
        probe_float = float(probe_value)
    except (TypeError, ValueError):
        probe_float = float("nan")
    if not np.isfinite(probe_float) or (
        abs(probe_float - float(query.probe_radius_nm)) > SASA_COMPAT_PROBE_RADIUS_ABS_TOL
    ):
        mismatched_fields.append("probe_radius_nm")

    sphere_value = metadata.get("n_sphere_points")
    try:
        sphere_points = int(sphere_value)
    except (TypeError, ValueError):
        sphere_points = -1
    if sphere_points != int(query.n_sphere_points):
        mismatched_fields.append("n_sphere_points")

    equilibration_value = str(metadata.get("equilibration", "")).strip()
    if equilibration_value != query.equilibration.strip():
        mismatched_fields.append("equilibration")

    units_value = metadata.get("units", "A^2")
    if str(units_value) != "A^2":
        mismatched_fields.append("units")

    mode_value = metadata.get("sasa_mode", "atom")
    if str(mode_value) != "atom":
        mismatched_fields.append("sasa_mode")

    selection_hash_matches: bool | None = None
    if query.selection is not None or query.context_selection is not None:
        selection_value = (
            query.selection if query.selection is not None else (query.context_selection or "")
        )
        context_value = (
            query.context_selection if query.context_selection is not None else selection_value
        )
        query_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=query.probe_radius_nm,
            n_sphere_points=query.n_sphere_points,
            selection=selection_value,
            context_selection=context_value,
            equilibration=query.equilibration,
        )
        metadata_hash = metadata.get("compatibility_hash")
        selection_hash_matches = str(metadata_hash) == query_hash

    return SASAArtifactCompatibility(
        is_compatible=len(mismatched_fields) == 0,
        is_legacy=is_legacy,
        schema_version=schema_version,
        selection_hash_matches=selection_hash_matches,
        mismatched_fields=tuple(mismatched_fields),
    )


def adapt_canonical_sasa_to_exposure(
    result: SASAComputationResult,
    *,
    exposure_threshold: float,
) -> Any:
    """Adapt canonical shared SASA result into exposure SASA trajectory format.

    Parameters
    ----------
    result : SASAComputationResult
        Canonical shared SASA computation result in Å².
    exposure_threshold : float
        Relative SASA threshold used by exposure consumers.

    Returns
    -------
    Any
        ``SASATrajectoryResult`` instance for exposure analysis.
    """
    from polyzymd.analyses.exposure._sasa_trajectory import SASATrajectoryResult
    from polyzymd.analyses.shared.aa_classification import MAX_ASA_TABLE, get_aa_class

    sasa_per_frame = result.residue_sasa_a2.astype(np.float32) * A2_TO_NM2
    resids = np.asarray(result.residue_resids, dtype=np.int32)
    resnames = [str(x).upper() for x in result.residue_resnames]
    aa_classes = [get_aa_class(name) for name in resnames]

    max_sasa_nm2 = np.asarray(
        [MAX_ASA_TABLE.get(name, 200.0) * A2_TO_NM2 for name in resnames],
        dtype=np.float32,
    )
    safe_max = np.where(max_sasa_nm2 > 0.0, max_sasa_nm2, 1.0).astype(np.float32)
    relative_sasa_per_frame = sasa_per_frame / safe_max[np.newaxis, :]

    n_frames = int(sasa_per_frame.shape[0])
    n_residues = int(sasa_per_frame.shape[1]) if sasa_per_frame.ndim == 2 else 0
    return SASATrajectoryResult(
        sasa_per_frame=sasa_per_frame,
        relative_sasa_per_frame=relative_sasa_per_frame.astype(np.float32),
        resids=resids,
        resnames=resnames,
        aa_classes=aa_classes,
        max_sasa_nm2=max_sasa_nm2,
        n_frames=n_frames,
        n_residues=n_residues,
        exposure_threshold=exposure_threshold,
    )


def find_sibling_sasa_artifacts(
    replicate_analysis_dir: Path,
    query: SASAArtifactCompatibilityQuery,
) -> list[SASASiblingArtifactMatch]:
    """Find compatible sibling SASA artifacts for a replicate analysis directory.

    Parameters
    ----------
    replicate_analysis_dir : Path
        Replicate analysis directory, for example ``analysis/<cond>/exposure/run_1``.
    query : SASAArtifactCompatibilityQuery
        Compatibility query used to pre-filter candidates.

    Returns
    -------
    list[SASASiblingArtifactMatch]
        Compatible sibling artifacts sorted by preference.
    """
    sibling_dir = replicate_analysis_dir.parent.parent / "sasa" / replicate_analysis_dir.name
    if not sibling_dir.exists():
        return []

    matches: list[SASASiblingArtifactMatch] = []
    for npz_path in sibling_dir.glob("sasa_*.npz"):
        metadata_path = npz_path.with_suffix(".json")
        if not metadata_path.exists():
            continue
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue

        compatibility = check_sasa_artifact_compatibility(metadata, query)
        if not compatibility.is_compatible:
            continue

        matches.append(
            SASASiblingArtifactMatch(
                sibling_analysis_dir=sibling_dir,
                npz_path=npz_path,
                metadata_path=metadata_path,
                metadata=metadata,
                compatibility=compatibility,
            )
        )

    def _selection_hash_rank(value: bool | None) -> int:
        if value is True:
            return 0
        if value is None:
            return 1
        return 2

    matches.sort(
        key=lambda item: (
            int(item.compatibility.is_legacy),
            _selection_hash_rank(item.compatibility.selection_hash_matches),
            item.npz_path.name,
        )
    )
    return matches
