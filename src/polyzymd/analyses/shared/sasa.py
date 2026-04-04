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
    warn_if_multi_chain_selection(
        atoms,
        selection,
        context=f"for SASA {role} selection in run '{run_label}'",
    )
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

    Returns
    -------
    SASAComputationResult
        Raw atom-level, residue-level, and total SASA traces in Å².
    """
    import mdtraj as md

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
    n_frames = max(0, stop_frame - start_frame)
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
        frames = np.arange(start_frame, stop_frame, dtype=np.int64)
        time_ns = (frames.astype(np.float64) * timestep_ps) / 1000.0
        return SASAComputationResult(
            atom_sasa_a2=np.empty((n_frames, target_indices.size), dtype=np.float64),
            residue_sasa_a2=np.empty((n_frames, 0), dtype=np.float64),
            total_sasa_a2=np.full(n_frames, np.nan, dtype=np.float64),
            frames=frames,
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

    xyz_nm = np.empty((n_frames, len(context_atoms), 3), dtype=np.float32)
    frames = np.arange(start_frame, stop_frame, dtype=np.int64)
    for out_idx, frame_idx in enumerate(frames.tolist()):
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

    atom_sasa_target_a2 = atom_sasa_nm2[:, target_local_indices] * NM2_TO_A2
    residue_sasa_a2 = np.empty((n_frames, len(residue_items)), dtype=np.float64)
    for idx, (_, atom_locals) in enumerate(residue_items):
        residue_sasa_a2[:, idx] = np.sum(atom_sasa_nm2[:, atom_locals], axis=1) * NM2_TO_A2

    total_sasa_a2 = np.sum(atom_sasa_target_a2, axis=1)
    time_ns = (frames.astype(np.float64) * timestep_ps) / 1000.0

    return SASAComputationResult(
        atom_sasa_a2=atom_sasa_target_a2,
        residue_sasa_a2=residue_sasa_a2,
        total_sasa_a2=total_sasa_a2,
        frames=frames,
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

    metadata = {
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
    }
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
            target_atom_indices=np.asarray(payload["target_atom_indices"], dtype=np.int64),
            context_atom_indices=np.asarray(payload["context_atom_indices"], dtype=np.int64),
            residue_keys=[str(v) for v in payload["residue_keys"].tolist()],
            residue_chainids=[str(v) for v in payload["residue_chainids"].tolist()],
            residue_resids=[int(v) for v in payload["residue_resids"].tolist()],
            residue_resnames=[str(v) for v in payload["residue_resnames"].tolist()],
        )

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    return result, metadata
