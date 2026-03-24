"""Catalytic triad analysis result models.

This module defines Pydantic models for storing catalytic triad analysis results:
- TriadResult: Single replicate triad analysis
- TriadAggregatedResult: Multi-replicate aggregated results

The key metric is "simultaneous contact fraction" - the percentage of frames
where ALL pairs in the triad are below the contact threshold simultaneously.
"""

from __future__ import annotations

from typing import ClassVar

import numpy as np
from pydantic import Field

from polyzymd.analysis.results.base import (
    AggregatedResultMixin,
    BaseAnalysisResult,
)
from polyzymd.analysis.results.distances import (
    DistancePairAggregatedResult,
    DistancePairResult,
)


class TriadResult(BaseAnalysisResult):
    """Catalytic triad analysis result for a single replicate.

    Stores per-pair distance results plus the simultaneous contact fraction,
    which is the key metric for triad integrity.

    Attributes
    ----------
    triad_name : str
        Name of the catalytic triad configuration
    triad_description : str | None
        Description of the triad (e.g., enzyme name, triad residues)
    pair_results : list[DistancePairResult]
        Distance analysis results for each pair in the triad
    threshold : float
        Contact threshold in Angstroms
    simultaneous_contact_fraction : float
        Fraction of frames where ALL pairs are below threshold
    n_frames_simultaneous : int
        Number of frames with simultaneous contact

    # Autocorrelation statistics for simultaneous contact
    sim_contact_sem : float | None
        SEM of simultaneous contact fraction (autocorrelation-corrected)
    sim_contact_correlation_time : float | None
        Correlation time of simultaneous contact timeseries
    sim_contact_n_independent : int | None
        Number of independent frames for contact timeseries
    sim_contact_warning : str | None
        Warning if statistics unreliable
    """

    analysis_type: ClassVar[str] = "catalytic_triad"

    # Triad identification
    triad_name: str = Field(..., description="Name of the catalytic triad")
    triad_description: str | None = Field(default=None, description="Description of the triad")

    # Per-pair results (reuse DistancePairResult)
    pair_results: list[DistancePairResult] = Field(
        ..., description="Distance results for each pair in the triad"
    )

    # Threshold for contact
    threshold: float = Field(default=3.5, description="Contact threshold in Angstroms")

    # Simultaneous contact analysis
    simultaneous_contact_fraction: float = Field(
        ..., description="Fraction of frames where ALL pairs are below threshold"
    )
    n_frames_simultaneous: int = Field(
        ..., description="Number of frames with simultaneous contact"
    )

    # Full timeseries (optional - can be large)
    simultaneous_contact_timeseries: list[bool] | None = Field(
        default=None,
        description="Per-frame simultaneous contact (True if all pairs below threshold)",
    )

    # Autocorrelation statistics for simultaneous contact
    sim_contact_sem: float | None = Field(
        default=None,
        description="SEM of simultaneous contact fraction (autocorrelation-corrected)",
    )
    sim_contact_correlation_time: float | None = Field(
        default=None, description="Correlation time of contact timeseries"
    )
    sim_contact_correlation_time_unit: str | None = Field(
        default=None, description="Unit of correlation time"
    )
    sim_contact_n_independent: int | None = Field(
        default=None, description="Number of independent frames for contact analysis"
    )
    sim_contact_warning: str | None = Field(
        default=None, description="Warning if statistics unreliable"
    )

    # Trajectory info
    n_frames_total: int = Field(..., description="Total frames in trajectory")
    n_frames_used: int = Field(..., description="Frames used after equilibration")

    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            f"Catalytic Triad Analysis: {self.triad_name}",
            "=" * 60,
        ]

        if self.triad_description:
            lines.append(f"Description: {self.triad_description}")

        lines.extend(
            [
                f"Replicate: {self.replicate}",
                f"Equilibration: {self._format_equilibration()}",
                f"Frames used: {self.n_frames_used}/{self.n_frames_total}",
                f"Contact threshold: {self.threshold:.1f} Å",
                "",
                "Per-pair distances:",
                "-" * 40,
            ]
        )

        # Per-pair summaries
        for pr in self.pair_results:
            mean_str = f"{pr.mean_distance:.2f}"
            if pr.sem_distance is not None:
                mean_str += f" ± {pr.sem_distance:.2f}"
            else:
                mean_str += f" ± {pr.std_distance:.2f}"

            frac_str = ""
            if pr.fraction_below_threshold is not None:
                frac_str = f" ({pr.fraction_below_threshold * 100:.1f}% below threshold)"

            lines.append(f"  {pr.pair_label}: {mean_str} Å{frac_str}")

        # Simultaneous contact
        lines.extend(
            [
                "",
                "Simultaneous triad contact:",
                "-" * 40,
            ]
        )

        pct = self.simultaneous_contact_fraction * 100
        if self.sim_contact_sem is not None:
            sem_pct = self.sim_contact_sem * 100
            lines.append(f"  Contact fraction: {pct:.1f} ± {sem_pct:.1f}%")
            if self.sim_contact_n_independent is not None:
                lines.append(f"  (n_independent = {self.sim_contact_n_independent})")
        else:
            lines.append(f"  Contact fraction: {pct:.1f}%")

        lines.append(f"  Frames in contact: {self.n_frames_simultaneous}/{self.n_frames_used}")

        # Correlation time
        if self.sim_contact_correlation_time is not None:
            unit = self.sim_contact_correlation_time_unit or "frames"
            lines.append(f"  Correlation time: {self.sim_contact_correlation_time:.1f} {unit}")

        # Warning
        if self.sim_contact_warning:
            lines.extend(["", f"WARNING: {self.sim_contact_warning}"])

        return "\n".join(lines)

    @property
    def n_pairs(self) -> int:
        """Number of pairs in the triad."""
        return len(self.pair_results)

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [pr.pair_label for pr in self.pair_results]


class TriadAggregatedResult(BaseAnalysisResult, AggregatedResultMixin):
    """Aggregated catalytic triad results across multiple replicates.

    Combines per-replicate triad analyses and computes overall statistics.

    Attributes
    ----------
    triad_name : str
        Name of the catalytic triad configuration
    pair_results : list[DistancePairAggregatedResult]
        Aggregated distance results for each pair
    threshold : float
        Contact threshold in Angstroms
    overall_simultaneous_contact : float
        Mean simultaneous contact fraction across replicates
    sem_simultaneous_contact : float
        SEM of simultaneous contact across replicates
    per_replicate_simultaneous : list[float]
        Simultaneous contact fraction from each replicate
    """

    analysis_type: ClassVar[str] = "catalytic_triad_aggregated"

    # Replicate info
    replicates: list[int] = Field(..., description="Replicate numbers included")
    n_replicates: int = Field(..., description="Number of replicates")

    # Triad identification
    triad_name: str = Field(..., description="Name of the catalytic triad")
    triad_description: str | None = Field(default=None, description="Description of the triad")

    # Per-pair aggregated results
    pair_results: list[DistancePairAggregatedResult] = Field(
        ..., description="Aggregated distance results for each pair"
    )

    # Threshold
    threshold: float = Field(default=3.5, description="Contact threshold in Angstroms")

    # Aggregated simultaneous contact
    overall_simultaneous_contact: float = Field(
        ..., description="Mean simultaneous contact fraction across replicates"
    )
    sem_simultaneous_contact: float = Field(
        ..., description="SEM of simultaneous contact across replicates"
    )
    per_replicate_simultaneous: list[float] = Field(
        ..., description="Simultaneous contact fraction from each replicate"
    )

    # Source files
    source_result_files: list[str] = Field(
        default_factory=list, description="Paths to individual replicate result files"
    )

    def summary(self) -> str:
        """Return human-readable summary."""
        rep_range = self.replicate_range
        lines = [
            f"Catalytic Triad Aggregated: {self.triad_name}",
            "=" * 60,
        ]

        if self.triad_description:
            lines.append(f"Description: {self.triad_description}")

        lines.extend(
            [
                f"Replicates: {rep_range} (n={self.n_replicates})",
                f"Equilibration: {self._format_equilibration()}",
                f"Contact threshold: {self.threshold:.1f} Å",
                "",
                "Per-pair distances (mean ± SEM across replicates):",
                "-" * 50,
            ]
        )

        # Per-pair summaries
        for pr in self.pair_results:
            mean_str = f"{pr.overall_mean:.2f} ± {pr.overall_sem:.2f}"
            frac_str = ""
            if pr.overall_fraction_below is not None:
                sem_f = pr.sem_fraction_below or 0
                frac_str = f" ({pr.overall_fraction_below * 100:.1f} ± {sem_f * 100:.1f}% below)"
            lines.append(f"  {pr.pair_label}: {mean_str} Å{frac_str}")

        # Simultaneous contact
        lines.extend(
            [
                "",
                "Simultaneous triad contact:",
                "-" * 50,
                f"  Overall: {self.overall_simultaneous_contact * 100:.1f} ± {self.sem_simultaneous_contact * 100:.1f}%",
                "",
                "  Per-replicate:",
            ]
        )

        for rep, frac in zip(self.replicates, self.per_replicate_simultaneous):
            lines.append(f"    Rep {rep}: {frac * 100:.1f}%")

        return "\n".join(lines)

    @property
    def n_pairs(self) -> int:
        """Number of pairs in the triad."""
        return len(self.pair_results)

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [pr.pair_label for pr in self.pair_results]
