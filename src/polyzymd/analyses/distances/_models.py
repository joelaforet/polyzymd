"""Distance analysis result models.

This module defines Pydantic models for storing distance analysis results:
- DistanceResult: Single replicate distance distributions
- DistanceConditionModel: Multi-replicate aggregated results

Supports both single pair and multiple pair analyses.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, ClassVar, Protocol, Sequence

from pydantic import Field

from polyzymd.analyses._framework.results_base import (
    AggregatedResultMixin,
    BaseAnalysisResult,
)


@dataclass(frozen=True)
class DistanceResultMetadata:
    """Non-serialized metadata shared by distance result factories.

    Parameters
    ----------
    config_hash : str
        Hash of the simulation configuration used for cache validation.
    polyzymd_version : str
        PolyzyMD version used to compute the result.
    replicate : int or None
        Replicate number for the result being constructed.
    equilibration_time : float
        Equilibration offset value.
    equilibration_unit : str
        Equilibration offset unit.
    """

    config_hash: str
    polyzymd_version: str
    replicate: int | None
    equilibration_time: float
    equilibration_unit: str

    def with_replicate(self, replicate: int | None) -> DistanceResultMetadata:
        """Return a copy with a different replicate value.

        Parameters
        ----------
        replicate : int or None
            Replicate value for the copied metadata.

        Returns
        -------
        DistanceResultMetadata
            Metadata copy with the requested replicate.
        """

        return replace(self, replicate=replicate)

    def to_base_kwargs(self, *, selection_string: str) -> dict[str, Any]:
        """Return keyword arguments for common ``BaseAnalysisResult`` fields.

        Parameters
        ----------
        selection_string : str
            Selection string to store on the serialized result.

        Returns
        -------
        dict[str, Any]
            Flat keyword arguments accepted by distance result constructors.
        """

        return {
            "config_hash": self.config_hash,
            "polyzymd_version": self.polyzymd_version,
            "replicate": self.replicate,
            "equilibration_time": self.equilibration_time,
            "equilibration_unit": self.equilibration_unit,
            "selection_string": selection_string,
        }


class _DistancePairPayloadProtocol(Protocol):
    """Structural protocol for runner payloads used by result factories."""

    pair_label: str
    selection1: str
    selection2: str
    distances: Any
    mean_distance: float
    std_distance: float
    median_distance: float
    min_distance: float
    max_distance: float
    sem_distance: float | None
    correlation_time: float | None
    correlation_time_unit: str | None
    n_independent_frames: int | None
    statistical_inefficiency: float | None
    autocorrelation_warning: str | None
    threshold: float | None
    fraction_below_threshold: float | None
    histogram_edges: Any
    histogram_counts: Any
    kde_x: Any
    kde_y: Any
    kde_peak: float | None
    kde_bandwidth: float | None
    n_frames_total: int
    n_frames_used: int


class _DistanceAggregatedStatsProtocol(Protocol):
    """Structural protocol for per-pair aggregation statistics."""

    mean_stats: Any
    median_stats: Any
    fraction_stats: Any
    kde_peak_stats: Any
    per_rep_means: list[float]
    per_rep_stds: list[float]
    per_rep_medians: list[float]
    per_rep_fractions: list[float]
    per_rep_kde_peaks: list[float]


class _DistancePairSourceProtocol(Protocol):
    """Structural protocol for pair identity fields in aggregation."""

    pair_label: str
    selection1: str
    selection2: str


def _as_serializable_list(value: Any) -> list[Any] | None:
    """Convert array-like values to plain lists for JSON serialization.

    Parameters
    ----------
    value : Any
        Value that may be a NumPy array, tuple, list, or ``None``.

    Returns
    -------
    list[Any] or None
        Plain list representation, or ``None`` when the input is ``None``.
    """

    if value is None:
        return None
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    return [value]


def _stringify_paths(paths: Sequence[str | Path]) -> list[str]:
    """Convert path-like values to strings for serialized result fields.

    Parameters
    ----------
    paths : sequence of str or Path
        Paths to stringify.

    Returns
    -------
    list[str]
        Stringified path list.
    """

    return [str(path) for path in paths]


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

    @classmethod
    def from_runner_payload(
        cls,
        metadata: DistanceResultMetadata,
        payload: _DistancePairPayloadProtocol,
        *,
        pair_label: str | None = None,
        store_distributions: bool = True,
        selection_string: str | None = None,
    ) -> DistancePairResult:
        """Construct a pair result from a trajectory runner payload.

        Parameters
        ----------
        metadata : DistanceResultMetadata
            Common result metadata to flatten into the serialized result.
        payload : _DistancePairPayloadProtocol
            Structural runner payload for one distance pair.
        pair_label : str or None, optional
            Label override for user-facing analyses, by default ``None``.
        store_distributions : bool, optional
            If ``False``, omit the full per-frame distance distribution.
        selection_string : str or None, optional
            Selection string override, by default ``None``.

        Returns
        -------
        DistancePairResult
            Flat Pydantic result preserving the current serialized schema.
        """

        selection = selection_string or f"{payload.selection1} : {payload.selection2}"
        distances = _as_serializable_list(payload.distances) if store_distributions else None
        return cls(
            **metadata.to_base_kwargs(selection_string=selection),
            pair_label=pair_label or payload.pair_label,
            selection1=payload.selection1,
            selection2=payload.selection2,
            distances=distances,
            mean_distance=payload.mean_distance,
            std_distance=payload.std_distance,
            median_distance=payload.median_distance,
            min_distance=payload.min_distance,
            max_distance=payload.max_distance,
            sem_distance=payload.sem_distance,
            correlation_time=payload.correlation_time,
            correlation_time_unit=payload.correlation_time_unit,
            n_independent_frames=payload.n_independent_frames,
            statistical_inefficiency=payload.statistical_inefficiency,
            autocorrelation_warning=payload.autocorrelation_warning,
            threshold=payload.threshold,
            fraction_below_threshold=payload.fraction_below_threshold,
            histogram_edges=_as_serializable_list(payload.histogram_edges),
            histogram_counts=_as_serializable_list(payload.histogram_counts),
            kde_x=_as_serializable_list(payload.kde_x),
            kde_y=_as_serializable_list(payload.kde_y),
            kde_peak=payload.kde_peak,
            kde_bandwidth=payload.kde_bandwidth,
            n_frames_total=payload.n_frames_total,
            n_frames_used=payload.n_frames_used,
        )

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
    settings_fingerprint: str | None = Field(
        default=None,
        description="Settings fingerprint used to validate distance caches",
    )
    trajectory_files: list[str] = Field(
        default_factory=list, description="Trajectory files analyzed"
    )

    @classmethod
    def from_pair_results(
        cls,
        metadata: DistanceResultMetadata,
        pair_results: Sequence[DistancePairResult],
        *,
        n_frames_total: int,
        n_frames_used: int,
        trajectory_files: Sequence[str | Path] = (),
        settings_fingerprint: str | None = None,
        selection_string: str | None = None,
    ) -> DistanceResult:
        """Construct a replicate distance result from pair results.

        Parameters
        ----------
        metadata : DistanceResultMetadata
            Common result metadata to flatten into the serialized result.
        pair_results : sequence of DistancePairResult
            Per-pair distance results.
        n_frames_total : int
            Total frames available in the trajectory.
        n_frames_used : int
            Frames analyzed after equilibration and striding.
        trajectory_files : sequence of str or Path, optional
            Source trajectory files, by default ``()``.
        settings_fingerprint : str or None, optional
            Settings identity token used for cache validation, by default
            ``None``.
        selection_string : str or None, optional
            Combined selection override, by default ``None``.

        Returns
        -------
        DistanceResult
            Flat replicate result preserving the current serialized schema.
        """

        selection = selection_string or "; ".join(
            f"({pair.selection1} : {pair.selection2})" for pair in pair_results
        )
        return cls(
            **metadata.to_base_kwargs(selection_string=selection),
            pair_results=list(pair_results),
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
            settings_fingerprint=settings_fingerprint,
            trajectory_files=_stringify_paths(trajectory_files),
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

    @classmethod
    def from_aggregated_stats(
        cls,
        metadata: DistanceResultMetadata,
        source_pair: _DistancePairSourceProtocol,
        stats: _DistanceAggregatedStatsProtocol,
        *,
        replicates: Sequence[int],
        threshold: float | None,
        pair_label: str | None = None,
        selection_string: str | None = None,
    ) -> DistancePairAggregatedResult:
        """Construct an aggregated pair result from shared statistics.

        Parameters
        ----------
        metadata : DistanceResultMetadata
            Common result metadata to flatten into the serialized result.
        source_pair : _DistancePairSourceProtocol
            Pair result that provides label and selection identity.
        stats : _DistanceAggregatedStatsProtocol
            Aggregated per-pair statistics.
        replicates : sequence of int
            Replicate numbers included in the aggregation.
        threshold : float or None
            Effective threshold for this pair.
        pair_label : str or None, optional
            Label override, by default ``None``.
        selection_string : str or None, optional
            Selection string override, by default ``None``.

        Returns
        -------
        DistancePairAggregatedResult
            Flat aggregated pair result preserving the current schema.
        """

        fraction_stats = getattr(stats, "fraction_stats", None)
        kde_peak_stats = getattr(stats, "kde_peak_stats", None)
        per_rep_fractions = list(getattr(stats, "per_rep_fractions", None) or [])
        per_rep_kde_peaks = list(getattr(stats, "per_rep_kde_peaks", None) or [])
        selection = selection_string or f"{source_pair.selection1} : {source_pair.selection2}"
        replicate_list = list(replicates)

        return cls(
            **metadata.to_base_kwargs(selection_string=selection),
            replicates=replicate_list,
            n_replicates=len(replicate_list),
            pair_label=pair_label or source_pair.pair_label,
            selection1=source_pair.selection1,
            selection2=source_pair.selection2,
            overall_mean=stats.mean_stats.mean,
            overall_sem=stats.mean_stats.sem,
            overall_median=stats.median_stats.mean,
            per_replicate_means=list(stats.per_rep_means),
            per_replicate_stds=list(stats.per_rep_stds),
            per_replicate_medians=list(stats.per_rep_medians),
            threshold=threshold,
            overall_fraction_below=(fraction_stats.mean if fraction_stats else None),
            sem_fraction_below=(fraction_stats.sem if fraction_stats else None),
            per_replicate_fractions_below=per_rep_fractions or None,
            overall_kde_peak=(kde_peak_stats.mean if kde_peak_stats else None),
            sem_kde_peak=(kde_peak_stats.sem if kde_peak_stats else None),
            per_replicate_kde_peaks=per_rep_kde_peaks or None,
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


class DistanceConditionModel(BaseAnalysisResult, AggregatedResultMixin):
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
    settings_fingerprint: str | None = Field(
        default=None,
        description="Settings fingerprint used to validate aggregated distance caches",
    )
    source_result_files: list[str] = Field(
        default_factory=list, description="Paths to individual replicate result files"
    )

    @classmethod
    def from_pair_results(
        cls,
        metadata: DistanceResultMetadata,
        pair_results: Sequence[DistancePairAggregatedResult],
        *,
        replicates: Sequence[int],
        source_result_files: Sequence[str | Path] = (),
        settings_fingerprint: str | None = None,
        selection_string: str | None = None,
    ) -> DistanceConditionModel:
        """Construct an aggregated distance result from pair results.

        Parameters
        ----------
        metadata : DistanceResultMetadata
            Common result metadata to flatten into the serialized result.
        pair_results : sequence of DistancePairAggregatedResult
            Aggregated per-pair results.
        replicates : sequence of int
            Replicate numbers included in the aggregation.
        source_result_files : sequence of str or Path, optional
            Source replicate result files, by default ``()``.
        settings_fingerprint : str or None, optional
            Settings identity token used for aggregate validation, by default
            ``None``.
        selection_string : str or None, optional
            Combined selection override, by default ``None``.

        Returns
        -------
        DistanceConditionModel
            Flat aggregated result preserving the current serialized schema.
        """

        selection = selection_string or "; ".join(
            f"({pair.selection1} : {pair.selection2})" for pair in pair_results
        )
        replicate_list = list(replicates)
        return cls(
            **metadata.to_base_kwargs(selection_string=selection),
            replicates=replicate_list,
            n_replicates=len(replicate_list),
            pair_results=list(pair_results),
            settings_fingerprint=settings_fingerprint,
            source_result_files=_stringify_paths(source_result_files),
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
