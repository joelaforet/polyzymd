"""Hydrogen-bond analysis result model skeletons.

This module defines replicate and aggregated result schemas used by the
hydrogen-bonds plugin lifecycle.
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field


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
    mean_unique_pairs_per_frame: float = 0.0
    std_unique_pairs_per_frame: float = 0.0
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


class HydrogenBondReplicatePayload(BaseModel):
    """Canonical replicate artifact payload fragment for hydrogen bonds."""

    config_hash: str | None = None
    polyzymd_version: str | None = None
    replicate: int | None = None
    equilibration_time: float | None = None
    equilibration_unit: str | None = None
    selection_string: str | None = None
    settings_fingerprint: str | None = Field(
        default=None,
        description="Short fingerprint of hydrogen-bond settings for cache validation",
    )
    timestep_ps: float | None = Field(
        default=None,
        description="Effective spacing between analyzed hydrogen-bond samples in ps",
    )
    raw_timestep_ps: float | None = Field(
        default=None,
        description="Raw trajectory frame spacing in ps before frame stride is applied",
    )
    frame_stride: int | None = Field(
        default=None,
        description="Frame stride applied when sampling hydrogen-bond frames",
    )
    summaries: list[HydrogenBondReplicateSummary]
    composition_entries: list[CompositionEntry] = Field(default_factory=list)


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

    model_config = ConfigDict(populate_by_name=True, extra="forbid")

    name: str
    mode: str
    group_names: list[str]
    n_replicates: int
    mean_hbonds_per_frame: float
    sem_hbonds_per_frame: float
    per_replicate_mean_hbonds: list[float]
    mean_unique_pairs_per_frame: float = 0.0
    sem_unique_pairs_per_frame: float = 0.0
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


class HydrogenBondConditionPayload(BaseModel):
    """Canonical condition artifact payload fragment for hydrogen bonds."""

    config_hash: str | None = None
    polyzymd_version: str | None = None
    replicate: int | None = None
    equilibration_time: float | None = None
    equilibration_unit: str | None = None
    selection_string: str | None = None
    settings_fingerprint: str | None = Field(
        default=None,
        description="Short fingerprint of hydrogen-bond settings for cache validation",
    )
    timestep_ps: float | None = Field(
        default=None,
        description="Effective spacing between analyzed hydrogen-bond samples in ps",
    )
    raw_timestep_ps: float | None = Field(
        default=None,
        description="Raw trajectory frame spacing in ps before frame stride is applied",
    )
    frame_stride: int | None = Field(
        default=None,
        description="Frame stride applied when sampling hydrogen-bond frames",
    )
    replicates: list[int]
    n_replicates: int
    summaries: list[HydrogenBondAggregatedSummary]
    composition_entries: list[AggregatedCompositionEntry] = Field(default_factory=list)
