"""Distance analysis result models.

This module defines Pydantic models for storing distance analysis results:
- DistanceResult: Single replicate distance distributions
- DistanceAggregatedResult: Multi-replicate aggregated results

Supports both single pair and multiple pair analyses.
"""

from __future__ import annotations

from typing import ClassVar

from pydantic import Field

from polyzymd.analysis.results.base import (
    AggregatedResultMixin,
    BaseAnalysisResult,
)


class DistancePairResult(BaseAnalysisResult):
    """Distance analysis result for a single atom pair in one replicate.

    Stores the full distribution of distances as well as summary statistics.
    Includes KDE-based distribution analysis and autocorrelation-corrected
    uncertainty quantification following LiveCoMS best practices.

    Attributes
    ----------
    pair_label : str
        Human-readable label for the pair (e.g., "Ser77_OG-His133_NE2")
    selection1 : str
        MDAnalysis selection for first atom/group
    selection2 : str
        MDAnalysis selection for second atom/group
    distances : list[float]
        All distance values (for distribution analysis)
    mean_distance : float
        Mean distance
    std_distance : float
        Standard deviation
    median_distance : float
        Median distance
    min_distance : float
        Minimum distance observed
    max_distance : float
        Maximum distance observed
    fraction_below_threshold : float | None
        Fraction of frames below threshold (if threshold specified)
    threshold : float | None
        Distance threshold used (if any)
    kde_peak : float | None
        Mode (most probable distance) from KDE
    correlation_time : float | None
        Estimated correlation time from ACF
    n_independent_frames : int | None
        Effective number of independent samples
    sem_distance : float | None
        Standard error of the mean (autocorrelation-corrected)
    """

    analysis_type: ClassVar[str] = "distance_pair"

    # Pair identification
    pair_label: str = Field(..., description="Human-readable pair label")
    selection1: str = Field(..., description="MDAnalysis selection for first atom/group")
    selection2: str = Field(..., description="MDAnalysis selection for second atom/group")

    # Full distribution (optional, can be large)
    distances: list[float] | None = Field(
        default=None, description="All distance values (Angstroms)"
    )

    # Summary statistics
    mean_distance: float = Field(..., description="Mean distance (Angstroms)")
    std_distance: float = Field(..., description="Standard deviation")
    median_distance: float = Field(..., description="Median distance")
    min_distance: float = Field(..., description="Minimum distance")
    max_distance: float = Field(..., description="Maximum distance")

    # Autocorrelation-corrected uncertainty (LiveCoMS best practices)
    sem_distance: float | None = Field(
        default=None,
        description="Standard error of the mean (autocorrelation-corrected)",
    )
    correlation_time: float | None = Field(
        default=None, description="Correlation time from ACF (same unit as timestep)"
    )
    correlation_time_unit: str | None = Field(
        default=None, description="Unit of correlation time (e.g., 'ps', 'ns')"
    )
    n_independent_frames: int | None = Field(
        default=None, description="Effective number of independent samples"
    )
    statistical_inefficiency: float | None = Field(
        default=None, description="Factor by which variance is inflated due to correlation"
    )
    autocorrelation_warning: str | None = Field(
        default=None, description="Warning if statistics may be unreliable"
    )

    # Threshold analysis
    threshold: float | None = Field(default=None, description="Distance threshold (if specified)")
    fraction_below_threshold: float | None = Field(
        default=None, description="Fraction of frames below threshold"
    )

    # Histogram data (pre-computed for plotting)
    histogram_edges: list[float] | None = Field(default=None, description="Histogram bin edges")
    histogram_counts: list[int] | None = Field(default=None, description="Histogram bin counts")

    # KDE data (for smooth distribution and mode estimation)
    kde_x: list[float] | None = Field(
        default=None, description="KDE evaluation points (distance values)"
    )
    kde_y: list[float] | None = Field(default=None, description="KDE density values")
    kde_peak: float | None = Field(
        default=None, description="Mode (most probable distance) from KDE"
    )
    kde_bandwidth: float | None = Field(
        default=None, description="KDE bandwidth used (Scott's rule)"
    )

    # Trajectory info
    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(..., description="Frames used after equilibration")

    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            f"Distance Analysis: {self.pair_label}",
            "=" * 50,
            f"Replicate: {self.replicate}",
            f"Selection 1: {self.selection1}",
            f"Selection 2: {self.selection2}",
            f"Equilibration: {self._format_equilibration()}",
            f"Frames used: {self.n_frames_used}/{self.n_frames_total}",
            "",
        ]

        # Mean with proper uncertainty
        if self.sem_distance is not None:
            lines.append(f"Mean: {self.mean_distance:.2f} ± {self.sem_distance:.2f} Å (SEM)")
            if self.n_independent_frames is not None:
                lines.append(f"  (n_independent = {self.n_independent_frames})")
        else:
            lines.append(f"Mean: {self.mean_distance:.2f} ± {self.std_distance:.2f} Å (std)")

        lines.append(f"Median: {self.median_distance:.2f} Å")

        # KDE peak (mode)
        if self.kde_peak is not None:
            lines.append(f"Mode (KDE peak): {self.kde_peak:.2f} Å")

        lines.append(f"Range: {self.min_distance:.2f} - {self.max_distance:.2f} Å")

        # Threshold analysis
        if self.threshold is not None and self.fraction_below_threshold is not None:
            pct = self.fraction_below_threshold * 100
            lines.append(f"Below {self.threshold:.1f} Å: {pct:.1f}%")

        # Correlation time info
        if self.correlation_time is not None:
            unit = self.correlation_time_unit or "frames"
            lines.append(f"Correlation time: {self.correlation_time:.1f} {unit}")

        # Warning if unreliable
        if self.autocorrelation_warning:
            lines.append("")
            lines.append(f"WARNING: {self.autocorrelation_warning}")

        return "\n".join(lines)


class DistanceResult(BaseAnalysisResult):
    """Distance analysis results for multiple pairs in one replicate.

    Container for analyzing multiple distance pairs simultaneously.
    """

    analysis_type: ClassVar[str] = "distances"

    # Collection of pair results
    pair_results: list[DistancePairResult] = Field(
        ..., description="Results for each distance pair"
    )

    # Trajectory info (shared across pairs)
    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(..., description="Frames used after equilibration")
    trajectory_files: list[str] = Field(
        default_factory=list, description="Trajectory files analyzed"
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            f"Distance Analysis (replicate {self.replicate})",
            "=" * 50,
            f"Pairs analyzed: {len(self.pair_results)}",
            f"Equilibration: {self._format_equilibration()}",
            f"Frames used: {self.n_frames_used}/{self.n_frames_total}",
            "",
        ]

        for pr in self.pair_results:
            lines.append(f"{pr.pair_label}: {pr.mean_distance:.2f} ± {pr.std_distance:.2f} Å")

        return "\n".join(lines)

    @property
    def n_pairs(self) -> int:
        """Number of distance pairs analyzed."""
        return len(self.pair_results)


class DistancePairAggregatedResult(BaseAnalysisResult, AggregatedResultMixin):
    """Aggregated distance results for one pair across replicates.

    Attributes
    ----------
    pair_label : str
        Human-readable label for the pair
    overall_mean : float
        Mean of replicate means
    overall_sem : float
        SEM across replicate means
    per_replicate_means : list[float]
        Mean distance from each replicate
    per_replicate_stds : list[float]
        Std dev from each replicate
    """

    analysis_type: ClassVar[str] = "distance_pair_aggregated"

    # Replicate info
    replicates: list[int] = Field(..., description="Replicate numbers included")
    n_replicates: int = Field(..., description="Number of replicates")

    # Pair identification
    pair_label: str = Field(..., description="Human-readable pair label")
    selection1: str = Field(..., description="MDAnalysis selection for first atom/group")
    selection2: str = Field(..., description="MDAnalysis selection for second atom/group")

    # Aggregated statistics
    overall_mean: float = Field(..., description="Mean of replicate means")
    overall_sem: float = Field(..., description="SEM across replicates")
    overall_median: float = Field(..., description="Median of replicate medians")

    # Per-replicate values
    per_replicate_means: list[float] = Field(..., description="Mean distance from each replicate")
    per_replicate_stds: list[float] = Field(..., description="Std dev from each replicate")
    per_replicate_medians: list[float] = Field(..., description="Median from each replicate")

    # Threshold analysis
    threshold: float | None = Field(default=None, description="Distance threshold (if specified)")
    overall_fraction_below: float | None = Field(
        default=None, description="Mean fraction below threshold across replicates"
    )
    sem_fraction_below: float | None = Field(
        default=None, description="SEM of fraction below threshold"
    )
    per_replicate_fractions_below: list[float] | None = Field(
        default=None, description="Fraction below threshold from each replicate"
    )

    # KDE statistics (aggregated from per-replicate KDE peaks)
    overall_kde_peak: float | None = Field(
        default=None, description="Mean of per-replicate KDE peaks (mode)"
    )
    sem_kde_peak: float | None = Field(
        default=None, description="SEM of KDE peaks across replicates"
    )
    per_replicate_kde_peaks: list[float] | None = Field(
        default=None, description="KDE peak (mode) from each replicate"
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        rep_range = self.replicate_range
        lines = [
            f"Distance Aggregated: {self.pair_label}",
            "=" * 50,
            f"Replicates: {rep_range}",
            f"Equilibration: {self._format_equilibration()}",
            "",
            f"Mean: {self.overall_mean:.2f} ± {self.overall_sem:.2f} Å",
            f"Median: {self.overall_median:.2f} Å",
        ]

        # KDE peak (mode)
        if self.overall_kde_peak is not None:
            sem_str = f" ± {self.sem_kde_peak:.2f}" if self.sem_kde_peak else ""
            lines.append(f"Mode (KDE peak): {self.overall_kde_peak:.2f}{sem_str} Å")

        if self.threshold is not None and self.overall_fraction_below is not None:
            pct = self.overall_fraction_below * 100
            sem_pct = (self.sem_fraction_below or 0) * 100
            lines.append(f"Below {self.threshold:.1f} Å: {pct:.1f} ± {sem_pct:.1f}%")

        lines.append("")
        lines.append("Per-replicate means:")
        for rep, mean in zip(self.replicates, self.per_replicate_means):
            lines.append(f"  Rep {rep}: {mean:.2f} Å")

        return "\n".join(lines)


class DistanceAggregatedResult(BaseAnalysisResult, AggregatedResultMixin):
    """Aggregated distance results for multiple pairs across replicates."""

    analysis_type: ClassVar[str] = "distances_aggregated"

    # Replicate info
    replicates: list[int] = Field(..., description="Replicate numbers included")
    n_replicates: int = Field(..., description="Number of replicates")

    # Collection of aggregated pair results
    pair_results: list[DistancePairAggregatedResult] = Field(
        ..., description="Aggregated results for each pair"
    )

    # Source files
    source_result_files: list[str] = Field(
        default_factory=list, description="Paths to individual replicate result files"
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        rep_range = self.replicate_range
        lines = [
            "Distance Aggregated Analysis",
            "=" * 50,
            f"Replicates: {rep_range}",
            f"Pairs analyzed: {len(self.pair_results)}",
            f"Equilibration: {self._format_equilibration()}",
            "",
        ]

        for pr in self.pair_results:
            lines.append(f"{pr.pair_label}: {pr.overall_mean:.2f} ± {pr.overall_sem:.2f} Å")

        return "\n".join(lines)

    @property
    def n_pairs(self) -> int:
        """Number of distance pairs analyzed."""
        return len(self.pair_results)
