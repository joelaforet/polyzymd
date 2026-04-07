"""Hydrogen-bond analysis result model skeletons.

This module defines replicate and aggregated result schemas used by the
hydrogen-bonds plugin lifecycle.
"""

from __future__ import annotations

from typing import ClassVar

from pydantic import BaseModel, Field

from polyzymd.analyses._results_base import AggregatedResultMixin, BaseAnalysisResult


class ResidueRef(BaseModel):
    """Reference to a specific residue.

    Parameters
    ----------
    chain_id : str
        Chain identifier.
    resid : int
        Residue identifier.
    resname : str
        Residue name.
    """

    chain_id: str
    resid: int
    resname: str

    @property
    def label(self) -> str:
        """Return a compact residue label.

        Returns
        -------
        str
            Formatted label ``<resname><resid>:<chain_id>``.
        """

        return f"{self.resname}{self.resid}:{self.chain_id}"


class DirectedResiduePairResult(BaseModel):
    """H-bond statistics for a directed donor→acceptor residue pair."""

    donor: ResidueRef
    acceptor: ResidueRef
    frames_present: int = Field(..., description="Number of frames with at least one H-bond")
    occupancy: float = Field(..., ge=0, le=1, description="Fraction of frames with H-bond")
    event_count: int = Field(..., ge=0, description="Total H-bond events across all frames")
    mean_events_per_frame: float = Field(..., ge=0, description="Mean H-bond events per frame")


class UndirectedResiduePairResult(BaseModel):
    """H-bond statistics for an undirected (merged) residue pair."""

    residue_a: ResidueRef
    residue_b: ResidueRef
    frames_present: int
    occupancy: float = Field(..., ge=0, le=1)
    event_count: int = Field(..., ge=0)
    mean_events_per_frame: float = Field(..., ge=0)


class HydrogenBondReplicateSummary(BaseModel):
    """Per-summary results for a single replicate."""

    name: str
    mode: str = Field(..., description="'between' or 'within'")
    group_names: list[str]
    n_frames_used: int
    mean_hbonds_per_frame: float
    fraction_frames_with_any_hbond: float = Field(..., ge=0, le=1)
    counts_per_frame: list[int]
    directed_residue_pairs: list[DirectedResiduePairResult] = Field(default_factory=list)
    undirected_residue_pairs: list[UndirectedResiduePairResult] = Field(default_factory=list)


class CompositionEntry(BaseModel):
    """Composition of H-bonds between two partitions for one replicate."""

    donor_partition: str
    acceptor_partition: str
    mean_hbonds_per_frame: float
    fraction_of_total: float = Field(..., ge=0, le=1)


class HydrogenBondResult(BaseAnalysisResult):
    """Per-replicate hydrogen-bond analysis result."""

    analysis_type: ClassVar[str] = "hydrogen_bonds"

    timestep_ps: float | None = None
    summaries: list[HydrogenBondReplicateSummary]
    composition_entries: list[CompositionEntry] = Field(default_factory=list)

    def summary(self) -> str:
        """Return a compact summary string.

        Returns
        -------
        str
            Human-readable summary of replicate-level hydrogen-bond results.
        """

        lines = [f"Hydrogen bond analysis — {len(self.summaries)} summaries"]
        for item in self.summaries:
            lines.append(
                "  "
                f"{item.name}: {item.mean_hbonds_per_frame:.2f} mean H-bonds/frame "
                f"({item.n_frames_used} frames)"
            )
        return "\n".join(lines)


class DirectedPairAggregate(BaseModel):
    """Aggregated directed pair across replicates."""

    donor: ResidueRef
    acceptor: ResidueRef
    mean_occupancy: float
    sem_occupancy: float
    mean_events_per_frame: float
    sem_events_per_frame: float
    per_replicate_occupancy: list[float]


class UndirectedPairAggregate(BaseModel):
    """Aggregated undirected pair across replicates."""

    residue_a: ResidueRef
    residue_b: ResidueRef
    mean_occupancy: float
    sem_occupancy: float
    mean_events_per_frame: float
    sem_events_per_frame: float
    per_replicate_occupancy: list[float]


class HydrogenBondAggregatedSummary(BaseModel):
    """Aggregated per-summary results across replicates."""

    name: str
    mode: str
    group_names: list[str]
    n_replicates: int
    mean_hbonds_per_frame: float
    sem_hbonds_per_frame: float
    per_replicate_mean_hbonds: list[float]
    mean_fraction_with_any: float
    sem_fraction_with_any: float
    per_replicate_fraction_with_any: list[float]
    directed_pairs: list[DirectedPairAggregate] = Field(default_factory=list)
    undirected_pairs: list[UndirectedPairAggregate] = Field(default_factory=list)


class AggregatedCompositionEntry(BaseModel):
    """Aggregated composition entry across replicates."""

    donor_partition: str
    acceptor_partition: str
    mean_hbonds_per_frame: float
    sem_hbonds_per_frame: float
    per_replicate_hbonds: list[float]
    mean_fraction_of_total: float
    sem_fraction_of_total: float
    per_replicate_fraction: list[float]


class HydrogenBondAggregatedResult(BaseAnalysisResult, AggregatedResultMixin):
    """Aggregated hydrogen bond results across replicates for one condition."""

    analysis_type: ClassVar[str] = "hydrogen_bonds_aggregated"

    timestep_ps: float | None = None
    replicates: list[int]
    n_replicates: int
    summaries: list[HydrogenBondAggregatedSummary]
    composition_entries: list[AggregatedCompositionEntry] = Field(default_factory=list)

    def summary(self) -> str:
        """Return a compact summary string.

        Returns
        -------
        str
            Human-readable summary of condition-level aggregated results.
        """

        lines = [
            "Hydrogen bonds aggregated — "
            f"{self.n_replicates} replicates, {len(self.summaries)} summaries"
        ]
        for item in self.summaries:
            lines.append(
                "  "
                f"{item.name}: {item.mean_hbonds_per_frame:.2f} ± "
                f"{item.sem_hbonds_per_frame:.2f} mean H-bonds/frame"
            )
        return "\n".join(lines)
