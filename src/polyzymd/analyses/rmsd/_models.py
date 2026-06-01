"""RMSD analysis result models.

This module defines Pydantic models for storing RMSD analysis results:
- RMSDRunResult: Single run (selection) within a single replicate
- RMSDReplicatePayload: All runs for one replicate
- RMSDRunAggregatePayload: One run aggregated across replicates
- RMSDConditionPayload: All runs aggregated across replicates

Supports multiple named RMSD selections (runs) per analysis, following
the multi-entry pattern established by the distances plugin.
"""

from __future__ import annotations

from typing import ClassVar

from pydantic import BaseModel, Field

from polyzymd.analyses._framework.results_base import BaseAnalysisResult


class RMSDArtifactPayloadBase(BaseModel):
    """Common artifact payload metadata for RMSD payload fragments."""

    config_hash: str
    polyzymd_version: str
    replicate: int | None = None
    equilibration_time: float
    equilibration_unit: str
    selection_string: str


class RMSDRunResult(BaseAnalysisResult):
    """RMSD result for a single named run in one replicate.

    Stores per-frame RMSD timeseries (as NPZ sidecar path) and summary
    statistics for one selection.

    Attributes
    ----------
    run_label : str
        Human-readable label for this run (e.g., "protein_backbone").
    selection : str
        MDAnalysis selection string used for RMSD calculation.
    alignment_selection : str
        MDAnalysis selection string used for trajectory alignment.
    reference_mode : str
        Reference structure mode: centroid, average, frame, or external.
    reference_frame : int | None
        1-indexed reference frame used (None for average mode).
    reference_file : str | None
        Path to external PDB used when ``reference_mode='external'``.
    mean_rmsd : float
        Mean RMSD over analyzed frames (Angstroms).
    std_rmsd : float
        Standard deviation of per-frame RMSD.
    median_rmsd : float
        Median RMSD.
    min_rmsd : float
        Minimum RMSD observed.
    max_rmsd : float
        Maximum RMSD observed.
    final_rmsd : float
        RMSD of the last frame (useful for convergence assessment).
    sem_rmsd : float | None
        Autocorrelation-corrected standard error of the mean.
    n_frames_total : int
        Total frames in trajectory.
    n_frames_used : int
        Frames used after equilibration.
    npz_path : str | None
        Path to NPZ sidecar containing per-frame RMSD timeseries.
    """

    analysis_type: ClassVar[str] = "rmsd_run"

    # Run identification
    run_label: str = Field(..., description="Human-readable run label")
    selection: str = Field(..., description="MDAnalysis selection for RMSD calculation")
    alignment_selection: str = Field(
        ..., description="MDAnalysis selection for trajectory alignment"
    )
    reference_mode: str = Field(
        ..., description="Reference structure mode: centroid, average, frame, or external"
    )
    reference_frame: int | None = Field(
        default=None, description="Reference frame (1-indexed), None for average mode"
    )
    reference_file: str | None = Field(
        default=None,
        description="Path to external PDB file used for external reference mode",
    )

    # Summary statistics
    mean_rmsd: float = Field(..., description="Mean RMSD (Angstroms)")
    std_rmsd: float = Field(..., description="Standard deviation of per-frame RMSD")
    median_rmsd: float = Field(..., description="Median RMSD")
    min_rmsd: float = Field(..., description="Minimum RMSD observed")
    max_rmsd: float = Field(..., description="Maximum RMSD observed")
    final_rmsd: float = Field(..., description="RMSD of last frame")

    # Autocorrelation-corrected uncertainty
    sem_rmsd: float | None = Field(
        default=None, description="Standard error of the mean (autocorrelation-corrected)"
    )
    correlation_time_unit: str | None = Field(default=None, description="Unit of correlation time")
    statistical_inefficiency: float | None = Field(
        default=None, description="Factor by which variance is inflated due to correlation"
    )
    autocorrelation_warning: str | None = Field(
        default=None, description="Warning if statistics may be unreliable"
    )

    # Convergence diagnostics
    converged: bool = Field(default=False, description="Whether RMSD trace converged")
    convergence_assessable: bool = Field(
        default=False,
        description="Whether convergence could be assessed from available data",
    )
    convergence_time_ns: float | None = Field(
        default=None,
        description="Detected convergence start time in ns",
    )
    convergence_message: str | None = Field(
        default=None,
        description="Convergence diagnostic message",
    )

    # Frame info
    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(..., description="Frames used after equilibration")

    # Timeseries sidecar
    npz_path: str | None = Field(
        default=None, description="Path to NPZ sidecar with per-frame RMSD timeseries"
    )

    # Time array metadata (for plotting)
    time_unit: str = Field(default="ns", description="Unit of time axis")
    timestep_ps: float | None = Field(
        default=None,
        description="Effective spacing between analyzed RMSD samples in ps",
    )
    raw_timestep_ps: float | None = Field(
        default=None,
        description="Raw trajectory frame spacing in ps before frame stride is applied",
    )
    frame_stride: int | None = Field(
        default=None,
        description="Frame stride applied when sampling this RMSD run",
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            f"RMSD Run: {self.run_label}",
            "=" * 50,
            f"Replicate: {self.replicate}",
            f"Selection: {self.selection}",
            f"Alignment: {self.alignment_selection}",
            f"Reference: {self.reference_mode}"
            + (f" (frame {self.reference_frame})" if self.reference_frame else ""),
            f"Equilibration: {self._format_equilibration()}",
            f"Frames used: {self.n_frames_used}/{self.n_frames_total}",
            "",
        ]

        if self.sem_rmsd is not None:
            lines.append(f"Mean: {self.mean_rmsd:.3f} +/- {self.sem_rmsd:.3f} A (SEM)")
            if self.n_independent_frames is not None:
                lines.append(f"  (n_independent = {self.n_independent_frames})")
        else:
            lines.append(f"Mean: {self.mean_rmsd:.3f} +/- {self.std_rmsd:.3f} A (std)")

        lines.append(f"Median: {self.median_rmsd:.3f} A")
        lines.append(f"Range: {self.min_rmsd:.3f} - {self.max_rmsd:.3f} A")
        lines.append(f"Final: {self.final_rmsd:.3f} A")

        if self.autocorrelation_warning:
            lines.append("")
            lines.append(f"WARNING: {self.autocorrelation_warning}")

        return "\n".join(lines)


class RMSDReplicatePayload(RMSDArtifactPayloadBase):
    """RMSD results for all runs in one replicate.

    Container for analyzing multiple RMSD selections simultaneously.
    """

    # Collection of run results
    run_results: list[RMSDRunResult] = Field(..., description="Results for each RMSD run")

    # Cache identity
    settings_fingerprint: str | None = Field(
        default=None,
        description="Settings fingerprint used to validate cached RMSD results",
    )

    # Trajectory info (shared across runs)
    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(..., description="Frames used after equilibration")
    trajectory_files: list[str] = Field(
        default_factory=list, description="Trajectory files analyzed"
    )

    @property
    def n_runs(self) -> int:
        """Number of RMSD runs analyzed."""
        return len(self.run_results)


class RMSDRunAggregatePayload(RMSDArtifactPayloadBase):
    """Aggregated RMSD results for one run across replicates.

    Attributes
    ----------
    run_label : str
        Human-readable label for this run.
    selection : str
        MDAnalysis selection used.
    overall_mean : float
        Mean of replicate means.
    overall_sem : float
        SEM across replicate means.
    per_replicate_means : list[float]
        Mean RMSD from each replicate.
    per_replicate_stds : list[float]
        Std dev from each replicate.
    """

    # Replicate info
    replicates: list[int] = Field(..., description="Replicate numbers included")
    n_replicates: int = Field(..., description="Number of replicates")

    # Run identification
    run_label: str = Field(..., description="Human-readable run label")
    selection: str = Field(..., description="MDAnalysis selection for RMSD")
    alignment_selection: str = Field(..., description="MDAnalysis selection for alignment")

    # Aggregated statistics
    overall_mean: float = Field(..., description="Mean of replicate means")
    overall_sem: float = Field(..., description="SEM across replicates")
    overall_median: float = Field(..., description="Median of replicate medians")

    # Per-replicate values
    per_replicate_means: list[float] = Field(..., description="Mean RMSD from each replicate")
    per_replicate_stds: list[float] = Field(..., description="Std dev from each replicate")
    per_replicate_medians: list[float] = Field(..., description="Median from each replicate")

    # Convergence summary
    per_replicate_convergence_times_ns: list[float | None] = Field(default_factory=list)
    per_replicate_convergence_assessable: list[bool] = Field(default_factory=list)
    n_converged_replicates: int = Field(default=0)
    n_assessable_replicates: int = Field(default=0)
    convergence_fraction: float = Field(default=0.0)
    all_converged: bool = Field(default=False)
    mean_convergence_time_ns: float | None = Field(default=None)
    median_convergence_time_ns: float | None = Field(default=None)


class RMSDConditionPayload(RMSDArtifactPayloadBase):
    """Aggregated RMSD results for all runs across replicates."""

    # Replicate info
    replicates: list[int] = Field(..., description="Replicate numbers included")
    n_replicates: int = Field(..., description="Number of replicates")

    # Collection of aggregated run results
    run_results: list[RMSDRunAggregatePayload] = Field(
        ..., description="Aggregated results for each run"
    )

    # Cache identity
    settings_fingerprint: str | None = Field(
        default=None,
        description="Settings fingerprint used to validate aggregated RMSD caches",
    )

    # Source files
    source_result_files: list[str] = Field(
        default_factory=list, description="Paths to individual replicate result files"
    )

    @property
    def n_runs(self) -> int:
        """Number of RMSD runs analyzed."""
        return len(self.run_results)
