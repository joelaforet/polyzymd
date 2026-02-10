"""RMSF analysis result models.

This module defines Pydantic models for storing RMSF analysis results:
- RMSFResult: Single replicate results
- RMSFAggregatedResult: Multi-replicate aggregated results

Both models support JSON serialization and include proper metadata
for reproducibility and cache validation.
"""

from __future__ import annotations

from typing import ClassVar

from pydantic import Field

from polyzymd.analysis.results.base import (
    AggregatedResultMixin,
    BaseAnalysisResult,
)


class RMSFResult(BaseAnalysisResult):
    """RMSF analysis result for a single replicate.

    Contains per-residue RMSF values and whole-protein statistics.

    Attributes
    ----------
    residue_ids : list[int]
        Residue identifiers (1-indexed, PyMOL convention)
    residue_names : list[str]
        Three-letter residue names (e.g., "ALA", "GLY")
    rmsf_values : list[float]
        RMSF value for each residue (Angstroms)
    mean_rmsf : float
        Mean RMSF over all residues
    std_rmsf : float
        Standard deviation of per-residue RMSF values
    min_rmsf : float
        Minimum per-residue RMSF
    max_rmsf : float
        Maximum per-residue RMSF

    Examples
    --------
    >>> result = RMSFResult.load("analysis/rmsf/run_1/rmsf_results.json")
    >>> print(result.summary())
    RMSF Analysis (replicate 1)
    ===========================
    Selection: protein and name CA
    Equilibration: 100ns
    Residues: 250
    Mean RMSF: 1.45 ± 0.32 Å
    Range: 0.52 - 3.21 Å
    """

    analysis_type: ClassVar[str] = "rmsf"

    # Per-residue data
    residue_ids: list[int] = Field(..., description="Residue IDs (1-indexed, PyMOL convention)")
    residue_names: list[str] = Field(..., description="Three-letter residue names")
    rmsf_values: list[float] = Field(..., description="RMSF values per residue (Angstroms)")

    # Summary statistics
    mean_rmsf: float = Field(..., description="Mean RMSF over all residues (Angstroms)")
    std_rmsf: float = Field(..., description="Std dev of per-residue RMSF values")
    min_rmsf: float = Field(..., description="Minimum per-residue RMSF")
    max_rmsf: float = Field(..., description="Maximum per-residue RMSF")

    # Reference information
    reference_frame: int | None = Field(
        default=None, description="Frame used as reference (1-indexed), None if external PDB"
    )
    reference_file: str | None = Field(
        default=None, description="Path to external reference PDB, if used"
    )

    # Trajectory info
    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(..., description="Frames used after equilibration cutoff")
    trajectory_files: list[str] = Field(
        default_factory=list, description="Trajectory files analyzed"
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            f"RMSF Analysis (replicate {self.replicate})",
            "=" * 40,
            f"Selection: {self.selection_string}",
            f"Equilibration: {self._format_equilibration()}",
            f"Correlation time: {self._format_correlation_time()}",
            f"Residues: {len(self.residue_ids)}",
            f"Frames used: {self.n_frames_used}/{self.n_frames_total}",
            f"Independent frames: {self.n_independent_frames or 'N/A'}",
            "",
            f"Mean RMSF: {self.mean_rmsf:.2f} ± {self.std_rmsf:.2f} Å",
            f"Range: {self.min_rmsf:.2f} - {self.max_rmsf:.2f} Å",
        ]
        return "\n".join(lines)

    @property
    def n_residues(self) -> int:
        """Number of residues analyzed."""
        return len(self.residue_ids)


class RMSFAggregatedResult(BaseAnalysisResult, AggregatedResultMixin):
    """Aggregated RMSF results across multiple replicates.

    Contains per-residue mean ± SEM across replicates, as well as
    whole-protein statistics with proper error propagation.

    Attributes
    ----------
    residue_ids : list[int]
        Residue identifiers (1-indexed)
    residue_names : list[str]
        Three-letter residue names
    mean_rmsf_per_residue : list[float]
        Mean RMSF per residue across replicates
    sem_rmsf_per_residue : list[float]
        SEM of RMSF per residue across replicates
    overall_mean_rmsf : float
        Mean of per-replicate whole-protein averages
    overall_sem_rmsf : float
        SEM across replicate means

    Examples
    --------
    >>> result = RMSFAggregatedResult.load("analysis/rmsf/aggregated/rmsf_reps1-5.json")
    >>> print(f"Overall RMSF: {result.overall_mean_rmsf:.2f} ± {result.overall_sem_rmsf:.2f} Å")
    Overall RMSF: 1.45 ± 0.08 Å
    """

    analysis_type: ClassVar[str] = "rmsf_aggregated"

    # Replicate info (from mixin - need to redefine for Pydantic)
    replicates: list[int] = Field(..., description="Replicate numbers included (1-indexed)")
    n_replicates: int = Field(..., description="Number of replicates")

    # Per-residue aggregated data
    residue_ids: list[int] = Field(..., description="Residue IDs (1-indexed)")
    residue_names: list[str] = Field(..., description="Three-letter residue names")
    mean_rmsf_per_residue: list[float] = Field(
        ..., description="Mean RMSF per residue across replicates"
    )
    sem_rmsf_per_residue: list[float] = Field(
        ..., description="SEM of RMSF per residue across replicates"
    )

    # Per-replicate whole-protein values (for transparency)
    per_replicate_mean_rmsf: list[float] = Field(
        ..., description="Whole-protein mean RMSF from each replicate"
    )

    # Overall statistics
    overall_mean_rmsf: float = Field(..., description="Mean of replicate whole-protein averages")
    overall_sem_rmsf: float = Field(..., description="SEM across replicate means")
    overall_min_rmsf: float = Field(..., description="Min of per-residue means")
    overall_max_rmsf: float = Field(..., description="Max of per-residue means")

    # Source files
    source_result_files: list[str] = Field(
        default_factory=list, description="Paths to individual replicate result files"
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        # Access replicate_range via the mixin property
        rep_range = self.replicate_range
        lines = [
            f"RMSF Aggregated Analysis (replicates {rep_range})",
            "=" * 50,
            f"Selection: {self.selection_string}",
            f"Equilibration: {self._format_equilibration()}",
            f"Replicates: {self.n_replicates}",
            f"Residues: {len(self.residue_ids)}",
            "",
            f"Overall RMSF: {self.overall_mean_rmsf:.2f} ± {self.overall_sem_rmsf:.2f} Å",
            f"Range of means: {self.overall_min_rmsf:.2f} - {self.overall_max_rmsf:.2f} Å",
            "",
            "Per-replicate whole-protein means:",
        ]
        for i, (rep, mean) in enumerate(zip(self.replicates, self.per_replicate_mean_rmsf)):
            lines.append(f"  Rep {rep}: {mean:.2f} Å")

        return "\n".join(lines)

    @property
    def n_residues(self) -> int:
        """Number of residues analyzed."""
        return len(self.residue_ids)

    @property
    def replicate_range(self) -> str:
        """Format replicate list as range string."""
        reps = sorted(self.replicates)
        if len(reps) == 0:
            return "none"
        if len(reps) == 1:
            return str(reps[0])
        if reps == list(range(reps[0], reps[-1] + 1)):
            return f"{reps[0]}-{reps[-1]}"
        return ",".join(map(str, reps))
