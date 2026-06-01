"""Rg analysis result models.

This module defines Pydantic models for storing Rg analysis results:
- RgRunResult: Single run (selection) within a single replicate
- RgReplicatePayload: All runs for one replicate
- RgRunAggregatePayload: One run aggregated across replicates
- RgConditionPayload: All runs aggregated across replicates

Supports multiple named Rg selections (runs) per analysis, following
the multi-entry pattern established by the distances plugin.
"""

from __future__ import annotations

from typing import ClassVar

from pydantic import BaseModel, Field

from polyzymd.analyses._framework.results_base import (
    AggregatedResultMixin,
    BaseAnalysisResult,
)


class RgSkippedRunResult(BaseModel):
    """Provenance for a configured Rg run skipped in one replicate."""

    run_label: str = Field(..., description="Human-readable run label")
    selection: str = Field(..., description="MDAnalysis selection that matched no atoms")
    replicate: int = Field(..., description="Replicate where the run was skipped")
    reason: str = Field(..., description="Human-readable skip reason")
    reason_code: str = Field(default="empty_selection", description="Machine-readable skip reason")


class RgRunResult(BaseAnalysisResult):
    """Rg result for a single named run in one replicate.

    Stores per-frame Rg timeseries (as NPZ sidecar path) and summary
    statistics for one selection.

    Attributes
    ----------
    run_label : str
        Human-readable label for this run (e.g., "protein_backbone").
    selection : str
        MDAnalysis selection string used for Rg calculation.
    mean_rg : float
        Mean Rg over analyzed frames (Angstroms).
    std_rg : float
        Standard deviation of per-frame Rg.
    median_rg : float
        Median Rg.
    min_rg : float
        Minimum Rg observed.
    max_rg : float
        Maximum Rg observed.
    final_rg : float
        Rg of the last frame (useful for convergence assessment).
    sem_rg : float | None
        Autocorrelation-corrected standard error of the mean.
    n_frames_total : int
        Total frames in trajectory.
    n_frames_used : int
        Frames used after equilibration.
    npz_path : str | None
        Path to NPZ sidecar containing per-frame Rg timeseries.
    calculation_mode : str
        Calculation mode used for this run ("selection" or "fragments").
    fragment_weighting : str | None
        Fragment weighting used in fragments mode, else None.
    mean_fragments_per_frame : float | None
        Mean number of fragments observed per frame.
    min_fragments_per_frame : int | None
        Minimum fragment count observed in any frame.
    max_fragments_per_frame : int | None
        Maximum fragment count observed in any frame.
    fragment_mean_rg : float | None
        Mean of all per-fragment Rg observations.
    fragment_std_rg : float | None
        Standard deviation of all per-fragment Rg observations.
    fragment_median_rg : float | None
        Median of all per-fragment Rg observations.
    fragment_min_rg : float | None
        Minimum per-fragment Rg observed.
    fragment_max_rg : float | None
        Maximum per-fragment Rg observed.
    fragment_rg_p10 : float | None
        10th percentile of per-fragment Rg values.
    fragment_rg_p25 : float | None
        25th percentile of per-fragment Rg values.
    fragment_rg_p50 : float | None
        50th percentile of per-fragment Rg values.
    fragment_rg_p75 : float | None
        75th percentile of per-fragment Rg values.
    fragment_rg_p90 : float | None
        90th percentile of per-fragment Rg values.
    """

    analysis_type: ClassVar[str] = "rg_run"

    # Run identification
    run_label: str = Field(..., description="Human-readable run label")
    selection: str = Field(..., description="MDAnalysis selection for Rg calculation")

    # Summary statistics
    mean_rg: float = Field(..., description="Mean Rg (Angstroms)")
    std_rg: float = Field(..., description="Standard deviation of per-frame Rg")
    median_rg: float = Field(..., description="Median Rg")
    min_rg: float = Field(..., description="Minimum Rg observed")
    max_rg: float = Field(..., description="Maximum Rg observed")
    final_rg: float = Field(..., description="Rg of last frame")

    # Autocorrelation-corrected uncertainty
    sem_rg: float | None = Field(
        default=None, description="Standard error of the mean (autocorrelation-corrected)"
    )
    correlation_time_unit: str | None = Field(default=None, description="Unit of correlation time")
    statistical_inefficiency: float | None = Field(
        default=None, description="Factor by which variance is inflated due to correlation"
    )
    autocorrelation_warning: str | None = Field(
        default=None, description="Warning if statistics may be unreliable"
    )

    # Frame info
    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(..., description="Frames used after equilibration")

    # Timeseries sidecar
    npz_path: str | None = Field(
        default=None, description="Path to NPZ sidecar with per-frame Rg timeseries"
    )

    # Fragment-aware metadata (only populated in fragments mode)
    calculation_mode: str = Field(
        default="selection", description='Calculation mode: "selection" or "fragments"'
    )
    fragment_weighting: str | None = Field(
        default=None,
        description="Fragment weighting scheme used, or None for selection mode",
    )
    mean_fragments_per_frame: float | None = Field(
        default=None, description="Mean number of fragments per frame"
    )
    min_fragments_per_frame: int | None = Field(
        default=None, description="Minimum fragments in any frame"
    )
    max_fragments_per_frame: int | None = Field(
        default=None, description="Maximum fragments in any frame"
    )

    # Fragment distribution summaries (only populated in fragments mode)
    fragment_mean_rg: float | None = Field(
        default=None, description="Mean of all per-fragment Rg observations"
    )
    fragment_std_rg: float | None = Field(
        default=None, description="Std dev of all per-fragment Rg observations"
    )
    fragment_median_rg: float | None = Field(
        default=None, description="Median of all per-fragment Rg observations"
    )
    fragment_min_rg: float | None = Field(default=None, description="Min per-fragment Rg observed")
    fragment_max_rg: float | None = Field(default=None, description="Max per-fragment Rg observed")

    # Quantiles of fragment Rg distribution
    fragment_rg_p10: float | None = Field(
        default=None, description="10th percentile of fragment Rg"
    )
    fragment_rg_p25: float | None = Field(
        default=None, description="25th percentile of fragment Rg"
    )
    fragment_rg_p50: float | None = Field(
        default=None, description="50th percentile (median) of fragment Rg"
    )
    fragment_rg_p75: float | None = Field(
        default=None, description="75th percentile of fragment Rg"
    )
    fragment_rg_p90: float | None = Field(
        default=None, description="90th percentile of fragment Rg"
    )

    # Time array metadata (for plotting)
    time_unit: str = Field(default="ns", description="Unit of time axis")
    timestep_ps: float | None = Field(
        default=None,
        description="Effective spacing between analyzed Rg samples in ps",
    )
    raw_timestep_ps: float | None = Field(
        default=None,
        description="Raw trajectory frame spacing in ps before frame stride is applied",
    )
    frame_stride: int | None = Field(
        default=None,
        description="Frame stride applied when sampling this Rg run",
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            f"Rg Run: {self.run_label}",
            "=" * 50,
            f"Replicate: {self.replicate}",
            f"Selection: {self.selection}",
            f"Equilibration: {self._format_equilibration()}",
            f"Frames used: {self.n_frames_used}/{self.n_frames_total}",
            "",
        ]

        if self.sem_rg is not None:
            lines.append(f"Mean: {self.mean_rg:.3f} +/- {self.sem_rg:.3f} A (SEM)")
            if self.n_independent_frames is not None:
                lines.append(f"  (n_independent = {self.n_independent_frames})")
        else:
            lines.append(f"Mean: {self.mean_rg:.3f} +/- {self.std_rg:.3f} A (std)")

        lines.append(f"Median: {self.median_rg:.3f} A")
        lines.append(f"Range: {self.min_rg:.3f} - {self.max_rg:.3f} A")
        lines.append(f"Final: {self.final_rg:.3f} A")

        if self.autocorrelation_warning:
            lines.append("")
            lines.append(f"WARNING: {self.autocorrelation_warning}")

        return "\n".join(lines)


class RgReplicatePayload(BaseAnalysisResult):
    """Rg results for all runs in one replicate.

    Container for analyzing multiple Rg selections simultaneously.
    """

    analysis_type: ClassVar[str] = "rg"

    # Collection of run results
    run_results: list[RgRunResult] = Field(..., description="Results for each Rg run")
    skipped_runs: list[RgSkippedRunResult] = Field(
        default_factory=list,
        description="Configured Rg runs skipped with provenance",
    )

    # Cache identity
    settings_fingerprint: str | None = Field(
        default=None,
        description="Settings fingerprint used to validate cached Rg results",
    )

    # Trajectory info (shared across runs)
    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(..., description="Frames used after equilibration")
    trajectory_files: list[str] = Field(
        default_factory=list, description="Trajectory files analyzed"
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            f"Rg Analysis (replicate {self.replicate})",
            "=" * 50,
            f"Runs analyzed: {len(self.run_results)}",
            f"Runs skipped: {len(self.skipped_runs)}",
            f"Equilibration: {self._format_equilibration()}",
            f"Frames used: {self.n_frames_used}/{self.n_frames_total}",
            "",
        ]

        for rr in self.run_results:
            lines.append(f"{rr.run_label}: {rr.mean_rg:.3f} +/- {rr.std_rg:.3f} A")

        return "\n".join(lines)

    @property
    def n_runs(self) -> int:
        """Number of Rg runs analyzed."""
        return len(self.run_results)


class RgRunAggregatePayload(BaseAnalysisResult, AggregatedResultMixin):
    """Aggregated Rg results for one run across replicates.

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
        Mean Rg from each replicate.
    per_replicate_stds : list[float]
        Std dev from each replicate.
    calculation_mode : str
        Calculation mode used for this run.
    fragment_weighting : str | None
        Fragment weighting scheme used, if fragment mode was used.
    overall_mean_fragments_per_frame : float | None
        Mean fragment count per frame across replicates.
    per_replicate_mean_fragments_per_frame : list[float] | None
        Per-replicate mean fragment count per frame.
    fragment_histogram_edges : list[float] | None
        Shared histogram bin edges for fragment Rg distribution.
    fragment_histogram_density_mean : list[float] | None
        Mean normalized fragment histogram density across replicates.
    fragment_histogram_density_sem : list[float] | None
        SEM of normalized fragment histogram density across replicates.
    reduced_histogram_edges : list[float] | None
        Shared histogram bin edges for reduced frame-level Rg values.
    reduced_histogram_density_mean : list[float] | None
        Mean normalized reduced-series histogram density across replicates.
    reduced_histogram_density_sem : list[float] | None
        SEM of normalized reduced-series histogram density across replicates.
    """

    analysis_type: ClassVar[str] = "rg_run_aggregated"

    # Replicate info
    replicates: list[int] = Field(..., description="Replicate numbers included")
    n_replicates: int = Field(..., description="Number of replicates")

    # Run identification
    run_label: str = Field(..., description="Human-readable run label")
    selection: str = Field(..., description="MDAnalysis selection for Rg")

    # Aggregated statistics
    overall_mean: float = Field(..., description="Mean of replicate means")
    overall_sem: float = Field(..., description="SEM across replicates")
    overall_median: float = Field(..., description="Median of replicate medians")

    # Per-replicate values
    per_replicate_means: list[float] = Field(..., description="Mean Rg from each replicate")
    per_replicate_stds: list[float] = Field(..., description="Std dev from each replicate")
    per_replicate_medians: list[float] = Field(..., description="Median from each replicate")

    # Fragment-aware metadata
    calculation_mode: str = Field(default="selection", description="Calculation mode used")
    fragment_weighting: str | None = Field(
        default=None, description="Fragment weighting scheme used"
    )
    overall_mean_fragments_per_frame: float | None = Field(
        default=None, description="Mean fragments/frame across replicates"
    )
    per_replicate_mean_fragments_per_frame: list[float] | None = Field(
        default=None, description="Mean fragments/frame for each replicate"
    )

    # Aggregated fragment distribution histograms
    fragment_histogram_edges: list[float] | None = Field(
        default=None, description="Common histogram bin edges for fragment Rg distribution"
    )
    fragment_histogram_density_mean: list[float] | None = Field(
        default=None, description="Mean normalized density across replicates"
    )
    fragment_histogram_density_sem: list[float] | None = Field(
        default=None, description="SEM of normalized density across replicates"
    )

    # Aggregated reduced-series histograms
    reduced_histogram_edges: list[float] | None = Field(
        default=None, description="Common histogram bin edges for reduced frame-level Rg"
    )
    reduced_histogram_density_mean: list[float] | None = Field(
        default=None,
        description="Mean normalized density of reduced series across replicates",
    )
    reduced_histogram_density_sem: list[float] | None = Field(
        default=None,
        description="SEM of normalized density of reduced series across replicates",
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        rep_range = self.replicate_range
        lines = [
            f"Rg Aggregated: {self.run_label}",
            "=" * 50,
            f"Replicates: {rep_range}",
            f"Selection: {self.selection}",
            f"Equilibration: {self._format_equilibration()}",
            "",
            f"Mean: {self.overall_mean:.3f} +/- {self.overall_sem:.3f} A",
            f"Median: {self.overall_median:.3f} A",
            "",
            "Per-replicate means:",
        ]

        for rep, mean in zip(self.replicates, self.per_replicate_means):
            lines.append(f"  Rep {rep}: {mean:.3f} A")

        return "\n".join(lines)


class RgConditionPayload(BaseAnalysisResult, AggregatedResultMixin):
    """Aggregated Rg results for all runs across replicates."""

    analysis_type: ClassVar[str] = "rg_aggregated"

    # Replicate info
    replicates: list[int] = Field(..., description="Replicate numbers included")
    n_replicates: int = Field(..., description="Number of replicates")

    # Collection of aggregated run results
    run_results: list[RgRunAggregatePayload] = Field(
        ..., description="Aggregated results for each run"
    )
    skipped_runs: list[RgSkippedRunResult] = Field(
        default_factory=list,
        description="Configured Rg runs skipped in contributing replicates",
    )

    # Cache identity
    settings_fingerprint: str | None = Field(
        default=None,
        description="Settings fingerprint used to validate aggregated Rg caches",
    )

    # Source files
    source_result_files: list[str] = Field(
        default_factory=list, description="Paths to individual replicate result files"
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        rep_range = self.replicate_range
        lines = [
            "Rg Aggregated Analysis",
            "=" * 50,
            f"Replicates: {rep_range}",
            f"Runs analyzed: {len(self.run_results)}",
            f"Equilibration: {self._format_equilibration()}",
            "",
        ]

        for rr in self.run_results:
            lines.append(f"{rr.run_label}: {rr.overall_mean:.3f} +/- {rr.overall_sem:.3f} A")

        return "\n".join(lines)

    @property
    def n_runs(self) -> int:
        """Number of Rg runs analyzed."""
        return len(self.run_results)
