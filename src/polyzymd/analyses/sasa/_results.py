"""SASA analysis result models."""

from __future__ import annotations

from typing import ClassVar

from pydantic import Field

from polyzymd.analyses._results_base import AggregatedResultMixin, BaseAnalysisResult


class SASARunResult(BaseAnalysisResult):
    """SASA result for a single named run in one replicate."""

    analysis_type: ClassVar[str] = "sasa_run"

    run_label: str = Field(..., description="Human-readable run label")
    target_selection: str = Field(..., description="MDAnalysis selection for target atoms")
    context_selection: str = Field(..., description="MDAnalysis selection for context atoms")

    mean_sasa: float = Field(..., description="Mean target SASA (A^2)")
    std_sasa: float = Field(..., description="Standard deviation of per-frame target SASA")
    median_sasa: float = Field(..., description="Median target SASA")
    min_sasa: float = Field(..., description="Minimum target SASA")
    max_sasa: float = Field(..., description="Maximum target SASA")
    final_sasa: float = Field(..., description="Target SASA of last analyzed frame")
    sem_sasa: float | None = Field(
        default=None,
        description="Standard error of the mean (autocorrelation-corrected)",
    )
    correlation_time_unit: str | None = Field(default=None, description="Unit of correlation time")
    statistical_inefficiency: float | None = Field(
        default=None,
        description="Factor by which variance is inflated due to correlation",
    )
    autocorrelation_warning: str | None = Field(
        default=None,
        description="Warning if autocorrelation-based statistics may be unreliable",
    )

    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(
        ...,
        description="Frames sampled for this run after equilibration and stride",
    )
    n_target_atoms: int = Field(..., description="Number of selected target atoms")
    n_context_atoms: int = Field(..., description="Number of selected context atoms")
    n_target_residues: int = Field(..., description="Number of selected target residues")
    zero_atom_selection: bool = Field(
        default=False,
        description="True when target or context selection matched zero atoms",
    )

    raw_npz_path: str | None = Field(
        default=None,
        description="Path to NPZ sidecar with raw SASA traces",
    )
    raw_metadata_path: str | None = Field(
        default=None,
        description="Path to JSON metadata sidecar for NPZ payload",
    )
    npz_path: str | None = Field(
        default=None,
        description="Deprecated alias for raw_npz_path",
    )
    metadata_path: str | None = Field(
        default=None,
        description="Deprecated alias for raw_metadata_path",
    )
    time_unit: str = Field(default="ns", description="Unit of time axis")
    timestep_ps: float | None = Field(
        default=None,
        description="Effective spacing between analyzed SASA samples in ps",
    )
    raw_timestep_ps: float | None = Field(
        default=None,
        description="Raw trajectory frame spacing in ps before SASA stride is applied",
    )
    frame_stride: int | None = Field(
        default=None,
        description="Frame stride applied when sampling this SASA run",
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        if self.zero_atom_selection:
            return (
                f"SASA Run: {self.run_label}\n"
                f"Selection matched zero atoms for target/context\n"
                f"Target: {self.target_selection}\n"
                f"Context: {self.context_selection}"
            )

        return (
            f"SASA Run: {self.run_label}\n"
            f"Mean SASA: {self.mean_sasa:.3f} A^2\n"
            f"Frames sampled: {self.n_frames_used}/{self.n_frames_total}"
        )


class SASAResult(BaseAnalysisResult):
    """SASA results for all runs in one replicate."""

    analysis_type: ClassVar[str] = "sasa"

    run_results: list[SASARunResult] = Field(..., description="Results for each SASA run")
    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(
        ...,
        description="Unique frames sampled across all configured runs",
    )
    trajectory_files: list[str] = Field(
        default_factory=list, description="Trajectory files analyzed"
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        return (
            f"SASA Analysis (replicate {self.replicate})\n"
            f"Runs analyzed: {len(self.run_results)}\n"
            f"Frames sampled across runs: {self.n_frames_used}/{self.n_frames_total}"
        )


class SASARunAggregatedResult(BaseAnalysisResult, AggregatedResultMixin):
    """Aggregated SASA results for one run across replicates."""

    analysis_type: ClassVar[str] = "sasa_run_aggregated"

    replicates: list[int] = Field(..., description="Replicate numbers included")
    n_replicates: int = Field(..., description="Number of replicates")

    run_label: str = Field(..., description="Human-readable run label")
    target_selection: str = Field(..., description="MDAnalysis target selection")
    context_selection: str = Field(..., description="MDAnalysis context selection")

    overall_mean: float = Field(..., description="Mean of replicate means")
    overall_sem: float = Field(..., description="SEM across replicate means")
    overall_median: float = Field(..., description="Median of replicate medians")
    overall_min: float = Field(..., description="Minimum of replicate minima")
    overall_max: float = Field(..., description="Maximum of replicate maxima")
    overall_final: float = Field(..., description="Mean of per-replicate final SASA")

    per_replicate_means: list[float] = Field(..., description="Mean SASA from each replicate")
    per_replicate_stds: list[float] = Field(..., description="Std SASA from each replicate")
    per_replicate_medians: list[float] = Field(..., description="Median SASA from each replicate")
    per_replicate_mins: list[float] = Field(..., description="Minimum SASA from each replicate")
    per_replicate_maxs: list[float] = Field(..., description="Maximum SASA from each replicate")
    per_replicate_finals: list[float] = Field(..., description="Final SASA from each replicate")

    n_target_atoms: int = Field(..., description="Target atom count")
    n_context_atoms: int | None = Field(
        ...,
        description="Context atom count, or None when counts vary across replicates",
    )
    per_replicate_context_atom_counts: list[int] = Field(
        default_factory=list,
        description="Context atom count from each contributing replicate",
    )
    n_context_atoms_variable: bool = Field(
        default=False,
        description="True when context atom counts vary across replicates",
    )
    n_target_residues: int = Field(..., description="Target residue count")
    zero_atom_selection: bool = Field(
        default=False, description="True when selection had zero atoms"
    )

    residue_keys: list[str] = Field(
        default_factory=list, description="Residue keys chain:resid:resname"
    )
    residue_chainids: list[str] = Field(default_factory=list, description="Chain IDs per residue")
    residue_resids: list[int] = Field(default_factory=list, description="Residue IDs per residue")
    residue_resnames: list[str] = Field(
        default_factory=list, description="Residue names per residue"
    )
    per_residue_mean_sasa: list[float] = Field(
        default_factory=list,
        description="Condition-level mean per-residue SASA in A^2",
    )
    per_residue_sem_sasa: list[float] = Field(
        default_factory=list,
        description="Condition-level SEM per-residue SASA in A^2",
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        return (
            f"SASA Aggregated: {self.run_label}\n"
            f"Replicates: {self.replicate_range}\n"
            f"Mean SASA: {self.overall_mean:.3f} +/- {self.overall_sem:.3f} A^2"
        )


class SASAAggregatedResult(BaseAnalysisResult, AggregatedResultMixin):
    """Aggregated SASA results for all runs across replicates."""

    analysis_type: ClassVar[str] = "sasa_aggregated"

    replicates: list[int] = Field(..., description="Replicate numbers included")
    n_replicates: int = Field(..., description="Number of replicates")
    run_results: list[SASARunAggregatedResult] = Field(
        ...,
        description="Aggregated results for each run",
    )
    settings_fingerprint: str | None = Field(
        default=None,
        description="Settings fingerprint used to validate aggregated SASA caches",
    )
    source_result_files: list[str] = Field(
        default_factory=list,
        description="Paths to individual replicate result files",
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        return (
            "SASA Aggregated Analysis\n"
            f"Replicates: {self.replicate_range}\n"
            f"Runs analyzed: {len(self.run_results)}"
        )
