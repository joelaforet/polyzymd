"""Secondary structure analysis result models.

Result classes for single-replicate and aggregated DSSP analysis.
Per-frame assignment matrices are stored as NPZ sidecars; persistence
fractions are stored in the JSON result for lightweight cross-replicate
aggregation.

Integer encoding
----------------
- 0 → coil   (DSSP ``C``)
- 1 → helix  (DSSP ``H``)
- 2 → strand (DSSP ``E``)
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
from pydantic import Field

from polyzymd.analyses._results_base import (
    AggregatedResultMixin,
    BaseAnalysisResult,
)

logger = logging.getLogger(__name__)

# Canonical mapping between DSSP letters and integer codes
SS_CHAR_TO_INT: dict[str, int] = {"C": 0, "H": 1, "E": 2}
SS_INT_TO_CHAR: dict[int, str] = {v: k for k, v in SS_CHAR_TO_INT.items()}
SS_LABELS: list[str] = ["coil", "helix", "strand"]


class SecondaryStructureResult(BaseAnalysisResult):
    """Single-replicate secondary structure result.

    The per-frame assignment matrix (shape ``(n_frames, n_residues)``,
    dtype ``int8``) is stored separately as an NPZ file to keep the JSON
    small.  The persistence fractions (per-residue fraction of frames
    spent in each SS state) are stored inline for fast aggregation.

    Attributes
    ----------
    analysis_type : str
        Always ``"secondary_structure"``.
    n_frames : int
        Number of frames analysed (after equilibration removal).
    n_residues : int
        Number of protein residues.
    residue_ids : list[int]
        1-indexed residue IDs matching the column order.
    residue_names : list[str]
        3-letter residue names, same order as ``residue_ids``.
    persistence_coil : list[float]
        Per-residue fraction of frames classified as coil.
    persistence_helix : list[float]
        Per-residue fraction of frames classified as helix.
    persistence_strand : list[float]
        Per-residue fraction of frames classified as strand.
    overall_helix_fraction : float
        Fraction of all (frame, residue) entries classified as helix.
    overall_strand_fraction : float
        Fraction classified as strand.
    overall_coil_fraction : float
        Fraction classified as coil.
    npz_path : str or None
        Path to the companion NPZ file containing the per-frame matrix.
    """

    analysis_type: ClassVar[str] = "secondary_structure"

    n_frames: int
    n_residues: int
    residue_ids: list[int]
    residue_names: list[str]

    # Per-residue persistence (fraction of frames in each state)
    persistence_coil: list[float]
    persistence_helix: list[float]
    persistence_strand: list[float]

    # Global content fractions
    overall_helix_fraction: float
    overall_strand_fraction: float
    overall_coil_fraction: float

    # Path to NPZ sidecar (set after save)
    npz_path: str | None = None

    def summary(self) -> str:
        """One-line summary for console output."""
        return (
            f"SS ({self.n_frames} frames, {self.n_residues} residues): "
            f"helix={self.overall_helix_fraction:.1%}, "
            f"strand={self.overall_strand_fraction:.1%}, "
            f"coil={self.overall_coil_fraction:.1%}"
        )

    # ------------------------------------------------------------------
    # NPZ helpers
    # ------------------------------------------------------------------

    def save_with_matrix(
        self,
        directory: Path | str,
        matrix: np.ndarray,
        filename_prefix: str = "secondary_structure",
    ) -> Path:
        """Save JSON result + NPZ per-frame matrix.

        Parameters
        ----------
        directory : Path or str
            Output directory (created if absent).
        matrix : np.ndarray
            Integer-encoded SS matrix, shape ``(n_frames, n_residues)``.
        filename_prefix : str
            Base name for the JSON and NPZ files.

        Returns
        -------
        Path
            Path to the saved JSON file.
        """
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)

        # Save NPZ
        npz_path = directory / f"{filename_prefix}_matrix.npz"
        np.savez_compressed(npz_path, ss_matrix=matrix.astype(np.int8))

        # Update internal path and save JSON
        self.npz_path = str(npz_path)
        json_path = directory / f"{filename_prefix}.json"
        self.save(json_path)
        return json_path

    @staticmethod
    def load_matrix(npz_path: Path | str) -> np.ndarray:
        """Load per-frame SS matrix from NPZ file.

        Parameters
        ----------
        npz_path : Path or str
            Path to the NPZ sidecar.

        Returns
        -------
        np.ndarray
            Integer-encoded matrix, shape ``(n_frames, n_residues)``.
        """
        data = np.load(str(npz_path))
        return data["ss_matrix"]


class SecondaryStructureAggregatedResult(AggregatedResultMixin, BaseAnalysisResult):
    """Aggregated secondary structure result across replicates.

    Stores mean and SEM of persistence fractions and overall content.

    Attributes
    ----------
    analysis_type : str
        Always ``"secondary_structure"``.
    n_replicates : int
        Number of replicates aggregated.
    replicates : list[int]
        Replicate numbers included.
    residue_ids : list[int]
        1-indexed residue IDs.
    residue_names : list[str]
        3-letter residue names.
    mean_persistence_coil : list[float]
        Mean coil persistence per residue.
    mean_persistence_helix : list[float]
        Mean helix persistence per residue.
    mean_persistence_strand : list[float]
        Mean strand persistence per residue.
    sem_persistence_coil : list[float]
        SEM of coil persistence per residue.
    sem_persistence_helix : list[float]
        SEM of helix persistence per residue.
    sem_persistence_strand : list[float]
        SEM of strand persistence per residue.
    mean_overall_helix : float
        Mean global helix fraction across replicates.
    mean_overall_strand : float
        Mean global strand fraction across replicates.
    mean_overall_coil : float
        Mean global coil fraction across replicates.
    sem_overall_helix : float
        SEM of global helix fraction.
    sem_overall_strand : float
        SEM of global strand fraction.
    sem_overall_coil : float
        SEM of global coil fraction.
    per_replicate_helix : list[float]
        Per-replicate overall helix fraction (for comparator).
    per_replicate_strand : list[float]
        Per-replicate overall strand fraction.
    per_replicate_coil : list[float]
        Per-replicate overall coil fraction.
    """

    analysis_type: ClassVar[str] = "secondary_structure"

    n_replicates: int
    replicates: list[int]
    residue_ids: list[int]
    residue_names: list[str] = Field(default_factory=list)

    # Mean persistence per residue
    mean_persistence_coil: list[float]
    mean_persistence_helix: list[float]
    mean_persistence_strand: list[float]

    # SEM persistence per residue
    sem_persistence_coil: list[float]
    sem_persistence_helix: list[float]
    sem_persistence_strand: list[float]

    # Overall content (mean ± SEM across replicates)
    mean_overall_helix: float
    mean_overall_strand: float
    mean_overall_coil: float
    sem_overall_helix: float
    sem_overall_strand: float
    sem_overall_coil: float

    # Per-replicate overall fractions (for statistical tests in comparator)
    per_replicate_helix: list[float]
    per_replicate_strand: list[float]
    per_replicate_coil: list[float]

    def summary(self) -> str:
        """One-line summary for console output."""
        return (
            f"SS aggregated ({self.n_replicates} reps): "
            f"helix={self.mean_overall_helix:.1%} ± {self.sem_overall_helix:.1%}, "
            f"strand={self.mean_overall_strand:.1%} ± {self.sem_overall_strand:.1%}, "
            f"coil={self.mean_overall_coil:.1%} ± {self.sem_overall_coil:.1%}"
        )
