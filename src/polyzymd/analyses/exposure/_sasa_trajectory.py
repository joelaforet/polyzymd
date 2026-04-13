"""SASATrajectoryResult and trajectory SASA computation.

Canonical location — moved from ``analysis.sasa.trajectory`` during the
analysis → analyses migration (Phase 5, Tier 1b).

Computes per-frame, per-residue SASA using MDTraj's shrake_rupley algorithm.
SASA is always computed on a **protein-only** sub-trajectory so that the
measured exposure is intrinsic to the protein and not affected by nearby
polymer atoms.

Design decisions (from Issue #33)
----------------------------------
- MDTraj shrake_rupley is significantly faster than MDAnalysis SASA.
- Protein-only SASA: select chain A atoms before calling shrake_rupley.
- Units: MDTraj returns nm^2; we keep nm^2 internally and convert to A^2 only
  when comparing against the Tien et al. MAX_ASA_TABLE (stored in A^2).
- Caching: optional NPZ + JSON sidecar written to analysis_dir/sasa/.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses.shared.aa_classification import MAX_ASA_TABLE, get_aa_class
from polyzymd.analyses.shared.config_hash import settings_fingerprint

if TYPE_CHECKING:
    from polyzymd.analyses.exposure._sasa_config import SASAConfig

logger = logging.getLogger(__name__)

# A^2 -> nm^2 conversion factor (MAX_ASA_TABLE is in A^2)
_ANG2_TO_NM2 = 1e-2


@dataclass
class SASATrajectoryResult:
    """Per-frame, per-residue SASA for a single MD trajectory.

    SASA is computed on the **protein-only** trajectory (chain A by default)
    so that exposure is intrinsic to the protein, not affected by polymer
    proximity.

    Attributes
    ----------
    sasa_per_frame : NDArray[np.float32]
        Raw SASA values, shape (n_frames, n_residues), in nm^2.
    relative_sasa_per_frame : NDArray[np.float32]
        SASA / max_SASA_for_restype, shape (n_frames, n_residues).
        Values > 1 are possible for highly exposed loop residues.
    resids : NDArray[np.int32]
        1-indexed protein residue IDs, shape (n_residues,).
    resnames : list[str]
        3-letter residue names, length n_residues.
    aa_classes : list[str]
        Amino-acid class labels, length n_residues.
    max_sasa_nm2 : NDArray[np.float32]
        Maximum SASA from Tien et al. 2013, shape (n_residues,), in nm^2.
    n_frames : int
        Number of frames analyzed.
    n_residues : int
        Number of protein residues.
    exposure_threshold : float
        Relative SASA threshold used for classifying exposed residues (0-1).
    trajectory_path : str
        Source trajectory path (for provenance).
    topology_path : str
        Source topology path (for provenance).
    """

    sasa_per_frame: NDArray[np.float32]
    relative_sasa_per_frame: NDArray[np.float32]
    resids: NDArray[np.int32]
    resnames: list[str]
    aa_classes: list[str]
    max_sasa_nm2: NDArray[np.float32]
    n_frames: int
    n_residues: int
    exposure_threshold: float
    trajectory_path: str = ""
    topology_path: str = ""

    # ------------------------------------------------------------------ #
    # Query helpers                                                        #
    # ------------------------------------------------------------------ #

    def _resid_to_idx(self, resid: int) -> int:
        """Convert 1-indexed resid to array index.

        Raises
        ------
        KeyError
            If resid is not present.
        """
        idx = np.searchsorted(self.resids, resid)
        if idx >= len(self.resids) or self.resids[idx] != resid:
            raise KeyError(f"Residue {resid} not found in SASATrajectoryResult")
        return int(idx)

    def is_exposed(self, frame: int, resid: int) -> bool:
        """Return True if *resid* is exposed at *frame*.

        Parameters
        ----------
        frame : int
            0-indexed frame number.
        resid : int
            1-indexed protein residue ID.
        """
        idx = self._resid_to_idx(resid)
        return bool(self.relative_sasa_per_frame[frame, idx] > self.exposure_threshold)

    def exposure_fraction(self, resid: int) -> float:
        """Fraction of frames where *resid* is exposed.

        Parameters
        ----------
        resid : int
            1-indexed protein residue ID.

        Returns
        -------
        float
            Value in [0, 1].
        """
        idx = self._resid_to_idx(resid)
        col = self.relative_sasa_per_frame[:, idx]
        return float(np.mean(col > self.exposure_threshold))

    def exposure_fraction_all(self) -> NDArray[np.float64]:
        """Fraction of frames each residue is exposed.

        Returns
        -------
        NDArray[np.float64]
            Shape (n_residues,). Entry i corresponds to self.resids[i].
        """
        return np.mean(self.relative_sasa_per_frame > self.exposure_threshold, axis=0)

    def exposed_mask_per_frame(self) -> NDArray[np.bool_]:
        """Boolean mask of exposure, shape (n_frames, n_residues)."""
        return self.relative_sasa_per_frame > self.exposure_threshold

    # ------------------------------------------------------------------ #
    # Serialisation                                                        #
    # ------------------------------------------------------------------ #

    def save(self, directory: Path | str) -> Path:
        """Save SASA result to *directory*/sasa_trajectory.npz + metadata.json.

        Parameters
        ----------
        directory : Path or str
            Target directory.  Created if it does not exist.

        Returns
        -------
        Path
            The directory where files were written.
        """
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)

        # Compressed numpy arrays
        np.savez_compressed(
            directory / "sasa_trajectory.npz",
            sasa_per_frame=self.sasa_per_frame,
            relative_sasa_per_frame=self.relative_sasa_per_frame,
            resids=self.resids,
            max_sasa_nm2=self.max_sasa_nm2,
        )

        # JSON metadata (human-readable sidecar)
        meta = {
            "resnames": self.resnames,
            "aa_classes": self.aa_classes,
            "n_frames": self.n_frames,
            "n_residues": self.n_residues,
            "exposure_threshold": self.exposure_threshold,
            "trajectory_path": self.trajectory_path,
            "topology_path": self.topology_path,
        }
        (directory / "sasa_metadata.json").write_text(json.dumps(meta, indent=2))

        logger.info(f"SASA trajectory saved to {directory}")
        return directory

    @classmethod
    def load(cls, directory: Path | str) -> "SASATrajectoryResult":
        """Load from a directory previously written by :meth:`save`.

        Parameters
        ----------
        directory : Path or str
            Directory containing sasa_trajectory.npz and sasa_metadata.json.

        Returns
        -------
        SASATrajectoryResult
        """
        directory = Path(directory)
        npz = np.load(directory / "sasa_trajectory.npz")
        meta = json.loads((directory / "sasa_metadata.json").read_text())

        return cls(
            sasa_per_frame=npz["sasa_per_frame"],
            relative_sasa_per_frame=npz["relative_sasa_per_frame"],
            resids=npz["resids"],
            max_sasa_nm2=npz["max_sasa_nm2"],
            resnames=meta["resnames"],
            aa_classes=meta["aa_classes"],
            n_frames=meta["n_frames"],
            n_residues=meta["n_residues"],
            exposure_threshold=meta["exposure_threshold"],
            trajectory_path=meta.get("trajectory_path", ""),
            topology_path=meta.get("topology_path", ""),
        )

    @classmethod
    def cache_path(cls, analysis_dir: Path | str, settings_fp: str | None = None) -> Path:
        """Return the cache directory for this result."""
        base_dir = Path(analysis_dir) / "sasa"
        if settings_fp:
            return base_dir / f"fp_{settings_fp}"
        return base_dir


# --------------------------------------------------------------------------- #
# Public computation function                                                   #
# --------------------------------------------------------------------------- #


def compute_trajectory_sasa(
    topology_path: Path | str,
    trajectory_path: Path | str | list[Path] | list[str],
    config: "SASAConfig | None" = None,
    analysis_dir: Path | str | None = None,
    recompute: bool = False,
    *,
    equilibration: str = "0ns",
) -> SASATrajectoryResult:
    """Compute (or load cached) per-frame, per-residue SASA for a trajectory.

    Uses MDTraj's ``shrake_rupley`` on a protein-only sub-trajectory.

    Parameters
    ----------
    topology_path : Path or str
        Path to the topology PDB file.
    trajectory_path : Path, str, or list thereof
        Path(s) to the trajectory file(s) (DCD, XTC, etc.).
        Accepts a single path or a list of paths for multi-segment
        (daisy-chain) trajectories. MDTraj concatenates segments
        automatically when given a list.
    config : SASAConfig, optional
        SASA configuration.  Uses defaults if None.
    analysis_dir : Path or str, optional
        Analysis output directory.  SASA cache is stored under
        ``analysis_dir/sasa/`` when ``config.cache_sasa`` is True.
    recompute : bool
        Force recomputation even if cache exists.
    equilibration : str, optional
        Equilibration label used for sibling SASA artifact compatibility checks,
        by default ``"0ns"``.

    Returns
    -------
    SASATrajectoryResult
        Per-frame SASA data for all protein residues.
    """
    # Lazy imports (heavy deps)
    import mdtraj as md

    from polyzymd.analyses.exposure._sasa_config import SASAConfig as _SASAConfig

    if config is None:
        config = _SASAConfig()

    topology_path = Path(topology_path)

    # Normalize trajectory_path to a list of strings for mdtraj
    if isinstance(trajectory_path, (str, Path)):
        traj_paths = [Path(trajectory_path)]
    else:
        traj_paths = [Path(p) for p in trajectory_path]
    traj_files_str = [str(p) for p in traj_paths]

    # --- Sibling SASA artifact reuse (Phase 2) ---
    if analysis_dir is not None and not recompute:
        sibling_result = _try_load_sibling_sasa(
            topology_path=topology_path,
            config=config,
            analysis_dir=Path(analysis_dir),
            equilibration=equilibration,
        )
        if sibling_result is not None:
            # Cache the adapted result for future runs
            sasa_settings_fp = settings_fingerprint(config)
            cache_dir = SASATrajectoryResult.cache_path(analysis_dir, settings_fp=sasa_settings_fp)
            if config.cache_sasa:
                sibling_result.save(cache_dir)
            return sibling_result

    # Check cache
    sasa_settings_fp = settings_fingerprint(config)
    cache_dir = (
        SASATrajectoryResult.cache_path(analysis_dir, settings_fp=sasa_settings_fp)
        if analysis_dir is not None
        else None
    )
    if cache_dir is not None and not recompute and (cache_dir / "sasa_trajectory.npz").exists():
        logger.info(f"Loading cached SASA from {cache_dir}")
        result = SASATrajectoryResult.load(cache_dir)
        # Validate threshold matches config
        if abs(result.exposure_threshold - config.exposure_threshold) > 1e-6:
            logger.warning(
                f"Cached SASA has threshold={result.exposure_threshold}, "
                f"but config specifies {config.exposure_threshold}. "
                "Using cached SASA (relative_sasa values are threshold-independent); "
                "re-run with recompute=True if you want a fresh computation."
            )
        return result

    n_segments = len(traj_paths)
    traj_label = (
        traj_paths[0].name
        if n_segments == 1
        else f"{traj_paths[0].name} + {n_segments - 1} more segment(s)"
    )
    logger.info(
        f"Computing SASA for {traj_label} "
        f"(threshold={config.exposure_threshold}, chain={config.chain_id})"
    )

    # Load full trajectory (mdtraj natively concatenates a list of files)
    traj = md.load(traj_files_str, top=str(topology_path))
    logger.info(f"Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")

    # Select protein chain only
    protein_indices = traj.topology.select(f"chainid {_chain_letter_to_index(config.chain_id)}")
    if len(protein_indices) == 0:
        raise ValueError(
            f"Chain '{config.chain_id}' not found in topology. "
            f"Available chains: {[c.index for c in traj.topology.chains]}. "
            "Check your Settings.chain_id configuration."
        )

    protein_traj = traj.atom_slice(protein_indices)
    logger.info(
        f"Protein sub-trajectory: {protein_traj.n_atoms} atoms, {protein_traj.n_residues} residues"
    )

    # Compute SASA -- returns (n_frames, n_residues) in nm^2
    sasa_nm2 = md.shrake_rupley(
        protein_traj,
        mode="residue",
        probe_radius=config.probe_radius_nm,
        n_sphere_points=config.n_sphere_points,
    ).astype(np.float32)

    # Build residue metadata from topology
    protein_top = protein_traj.topology
    resids = []
    resnames = []
    aa_classes = []
    max_sasa_nm2_list = []

    for res in protein_top.residues:
        resids.append(res.resSeq)
        rname = res.name.upper()
        resnames.append(rname)
        aa_classes.append(get_aa_class(rname))
        # MAX_ASA_TABLE is in A^2 -> convert to nm^2
        max_ang2 = MAX_ASA_TABLE.get(rname, 200.0)
        max_sasa_nm2_list.append(max_ang2 * _ANG2_TO_NM2)

    resids_arr = np.array(resids, dtype=np.int32)
    max_sasa_arr = np.array(max_sasa_nm2_list, dtype=np.float32)

    # Relative SASA: divide each residue column by its max SASA (nm^2)
    # Guard against max_sasa == 0 (shouldn't happen with standard AA)
    safe_max = np.where(max_sasa_arr > 0, max_sasa_arr, 1.0)
    relative_sasa = sasa_nm2 / safe_max[np.newaxis, :]  # (n_frames, n_residues)

    n_frames, n_residues = sasa_nm2.shape

    result = SASATrajectoryResult(
        sasa_per_frame=sasa_nm2,
        relative_sasa_per_frame=relative_sasa.astype(np.float32),
        resids=resids_arr,
        resnames=resnames,
        aa_classes=aa_classes,
        max_sasa_nm2=max_sasa_arr,
        n_frames=n_frames,
        n_residues=n_residues,
        exposure_threshold=config.exposure_threshold,
        trajectory_path="; ".join(traj_files_str),
        topology_path=str(topology_path),
    )

    # Cache to disk
    if cache_dir is not None and config.cache_sasa:
        result.save(cache_dir)

    return result


def _resolve_protein_indices_from_topology(
    topology_path: Path,
    chain_id: str,
) -> NDArray[np.int64]:
    """Resolve sorted protein atom indices from topology without loading full trajectory.

    Uses mdtraj.load_topology() which reads only the PDB structure, not
    trajectory frames — fast even for large systems.

    Parameters
    ----------
    topology_path : Path
        Path to the topology PDB file.
    chain_id : str
        Chain letter (e.g. "A") to select.

    Returns
    -------
    NDArray[np.int64]
        Sorted atom indices for the protein chain.

    Raises
    ------
    ValueError
        If the chain is not found in the topology.
    """
    import mdtraj as md

    top = md.load_topology(str(topology_path))
    chain_idx = _chain_letter_to_index(chain_id)
    indices = top.select(f"chainid {chain_idx}")
    if len(indices) == 0:
        raise ValueError(
            f"Chain '{chain_id}' (index {chain_idx}) not found in topology {topology_path}"
        )
    return np.sort(np.asarray(indices, dtype=np.int64))


def _try_load_sibling_sasa(
    topology_path: Path,
    config: "SASAConfig",
    analysis_dir: Path,
    equilibration: str,
) -> SASATrajectoryResult | None:
    """Attempt to load a compatible SASA artifact from the sibling sasa plugin.

    Performs the full two-tier compatibility check:

    1. Metadata-level filtering via ``find_sibling_sasa_artifacts()``
    2. Atom-index verification: loads NPZ ``target_atom_indices`` and compares
       against protein indices resolved from topology.

    Parameters
    ----------
    topology_path : Path
        Path to the topology PDB file.
    config : SASAConfig
        Exposure SASA configuration (probe radius, sphere points, chain, threshold).
    analysis_dir : Path
        Per-replicate analysis directory (e.g. ``analysis/<cond>/exposure/run_1``).
    equilibration : str
        Equilibration label.

    Returns
    -------
    SASATrajectoryResult or None
        Adapted SASA result if a compatible sibling was found, else ``None``.
    """
    from polyzymd.analyses.shared.sasa import (
        SASAArtifactCompatibilityQuery,
        adapt_canonical_sasa_to_exposure,
        find_sibling_sasa_artifacts,
        load_sasa_artifacts,
    )

    # Build query — exposure wants protein-only SASA, so target and context
    # are the same selection (protein chain only, no polymer)
    # The selection strings here are advisory; definitive check is atom indices
    protein_selection = f"protein and chainid {_chain_letter_to_index(config.chain_id)}"
    query = SASAArtifactCompatibilityQuery(
        probe_radius_nm=config.probe_radius_nm,
        n_sphere_points=config.n_sphere_points,
        equilibration=equilibration,
        selection=protein_selection,
        context_selection=protein_selection,
    )

    candidates = find_sibling_sasa_artifacts(analysis_dir, query)
    if not candidates:
        logger.debug("No sibling SASA artifacts found for %s", analysis_dir)
        return None

    # Resolve expected protein atom indices from topology (fast, no trajectory load)
    try:
        expected_indices = _resolve_protein_indices_from_topology(topology_path, config.chain_id)
    except (ValueError, OSError) as exc:
        logger.debug("Cannot resolve protein indices from topology: %s", exc)
        return None

    for candidate in candidates:
        try:
            sasa_result, _meta = load_sasa_artifacts(candidate.npz_path, candidate.metadata_path)
        except (OSError, KeyError, ValueError) as exc:
            logger.debug("Corrupted sibling artifact %s: %s", candidate.npz_path, exc)
            continue

        # Tier 2: atom-index comparison
        stored_target = np.sort(sasa_result.target_atom_indices)
        stored_context = np.sort(sasa_result.context_atom_indices)

        # Exposure wants protein-only SASA: target == context == protein atoms
        if not (
            np.array_equal(stored_target, expected_indices)
            and np.array_equal(stored_context, expected_indices)
        ):
            logger.debug(
                "Sibling %s: atom-index mismatch "
                "(target: %d vs %d expected, context: %d vs %d expected)",
                candidate.npz_path.name,
                len(stored_target),
                len(expected_indices),
                len(stored_context),
                len(expected_indices),
            )
            continue

        # Match found — adapt to exposure format
        try:
            adapted = adapt_canonical_sasa_to_exposure(
                sasa_result,
                exposure_threshold=config.exposure_threshold,
            )
        except (ValueError, IndexError, TypeError, KeyError) as exc:
            logger.debug("Sibling %s: adaptation failed: %s", candidate.npz_path, exc)
            continue

        logger.info(
            "Reusing sibling SASA artifact: %s (adapted from sasa plugin)",
            candidate.npz_path,
        )
        return adapted

    logger.debug(
        "No sibling SASA artifacts matched atom-index check for %s",
        analysis_dir,
    )
    return None


def _chain_letter_to_index(chain_id: str) -> int:
    """Convert chain letter (A, B, C...) to 0-indexed integer for MDTraj.

    MDTraj's ``topology.select`` uses 0-based ``chainid`` integers.
    """
    return ord(chain_id.upper()) - ord("A")
