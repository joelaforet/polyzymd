"""SASATrajectoryResult and trajectory SASA computation.

Computes per-frame, per-residue SASA using MDTraj's shrake_rupley algorithm.
SASA is always computed on a **protein-only** sub-trajectory so that the
measured exposure is intrinsic to the protein and not affected by nearby
polymer atoms.

Design decisions (from Issue #33)
----------------------------------
- MDTraj shrake_rupley is significantly faster than MDAnalysis SASA.
- Protein-only SASA: select chain A atoms before calling shrake_rupley.
- Units: MDTraj returns nm²; we keep nm² internally and convert to Å² only
  when comparing against the Tien et al. MAX_ASA_TABLE (stored in Å²).
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

from polyzymd.analysis.common.aa_classification import MAX_ASA_TABLE, get_aa_class

if TYPE_CHECKING:
    from polyzymd.analysis.sasa.config import SASAConfig

logger = logging.getLogger(__name__)

# Ų → nm² conversion factor (MAX_ASA_TABLE is in Å²)
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
        Raw SASA values, shape (n_frames, n_residues), in nm².
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
        Maximum SASA from Tien et al. 2013, shape (n_residues,), in nm².
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
    def cache_path(cls, analysis_dir: Path | str) -> Path:
        """Return the cache directory for this result."""
        return Path(analysis_dir) / "sasa"


# --------------------------------------------------------------------------- #
# Public computation function                                                   #
# --------------------------------------------------------------------------- #


def compute_trajectory_sasa(
    topology_path: Path | str,
    trajectory_path: Path | str | list[Path] | list[str],
    config: "SASAConfig | None" = None,
    analysis_dir: Path | str | None = None,
    recompute: bool = False,
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

    Returns
    -------
    SASATrajectoryResult
        Per-frame SASA data for all protein residues.
    """
    # Lazy imports (heavy deps)
    import mdtraj as md

    from polyzymd.analysis.sasa.config import SASAConfig as _SASAConfig

    if config is None:
        config = _SASAConfig()

    topology_path = Path(topology_path)

    # Normalize trajectory_path to a list of strings for mdtraj
    if isinstance(trajectory_path, (str, Path)):
        traj_paths = [Path(trajectory_path)]
    else:
        traj_paths = [Path(p) for p in trajectory_path]
    traj_files_str = [str(p) for p in traj_paths]

    # Check cache
    cache_dir = SASATrajectoryResult.cache_path(analysis_dir) if analysis_dir is not None else None
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
        # Fallback: select all protein atoms
        protein_indices = traj.topology.select("protein")
        logger.warning(
            f"Chain '{config.chain_id}' not found; falling back to 'protein' selection "
            f"({len(protein_indices)} atoms)"
        )

    protein_traj = traj.atom_slice(protein_indices)
    logger.info(
        f"Protein sub-trajectory: {protein_traj.n_atoms} atoms, {protein_traj.n_residues} residues"
    )

    # Compute SASA — returns (n_frames, n_residues) in nm²
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
        # MAX_ASA_TABLE is in Å² → convert to nm²
        max_ang2 = MAX_ASA_TABLE.get(rname, 200.0)
        max_sasa_nm2_list.append(max_ang2 * _ANG2_TO_NM2)

    resids_arr = np.array(resids, dtype=np.int32)
    max_sasa_arr = np.array(max_sasa_nm2_list, dtype=np.float32)

    # Relative SASA: divide each residue column by its max SASA (nm²)
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


def _chain_letter_to_index(chain_id: str) -> int:
    """Convert chain letter (A, B, C...) to 0-indexed integer for MDTraj.

    MDTraj's ``topology.select`` uses 0-based ``chainid`` integers.
    """
    return ord(chain_id.upper()) - ord("A")
