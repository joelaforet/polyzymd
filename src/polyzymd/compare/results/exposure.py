"""Exposure dynamics condition summary and comparison result models.

These classes inherit from the base classes in compare/core/ and add
exposure-dynamics-specific fields for chaperone event analysis.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, ClassVar

from pydantic import Field

from polyzymd import __version__
from polyzymd.compare.core.base import (
    ANOVASummary,
    BaseComparisonResult,
    BaseConditionSummary,
    PairwiseComparison,
)

# ---------------------------------------------------------------------------
# Per-condition summary
# ---------------------------------------------------------------------------


class ExposureConditionSummary(BaseConditionSummary):
    """Summary statistics for one condition in an exposure dynamics comparison.

    Attributes
    ----------
    label : str
        Display name for this condition.
    config_path : str
        Path to the simulation config file.
    n_replicates : int
        Number of replicates included.
    replicate_values : list[float]
        Per-replicate mean chaperone fraction across transient residues.
    mean_transient_fraction : float
        Mean fraction of protein residues that are transiently exposed,
        averaged across replicates.
    sem_transient_fraction : float
        Standard error of mean_transient_fraction.
    mean_chaperone_fraction : float
        Mean chaperone fraction (chaperone events / total exposed windows)
        across transient residues and replicates.
    sem_chaperone_fraction : float
        Standard error of mean_chaperone_fraction.
    mean_n_transient : float
        Mean number of transient residues across replicates.
    mean_total_chaperone_events : float
        Mean total chaperone event count across replicates.
    mean_total_unassisted_events : float
        Mean total unassisted event count across replicates.
    enrichment_by_polymer_type : dict[str, dict[str, float]]
        Nested dict: polymer_type → aa_group → mean enrichment_residue.
    polymer_types : list[str]
        Polymer types present in this condition.
    aa_groups : list[str]
        Amino-acid groups present in this condition.
    """

    mean_transient_fraction: float = Field(
        ..., description="Mean fraction of residues that are transiently exposed"
    )
    sem_transient_fraction: float = Field(..., description="SEM of transient_fraction")
    mean_chaperone_fraction: float = Field(
        ..., description="Mean chaperone fraction across transient residues"
    )
    sem_chaperone_fraction: float = Field(..., description="SEM of chaperone_fraction")
    mean_n_transient: float = Field(..., description="Mean number of transient residues")
    mean_total_chaperone_events: float = Field(
        default=0.0, description="Mean total chaperone event count"
    )
    mean_total_unassisted_events: float = Field(
        default=0.0, description="Mean total unassisted event count"
    )
    enrichment_by_polymer_type: dict[str, dict[str, float]] = Field(
        default_factory=dict,
        description="polymer_type → aa_group → mean residue-based enrichment",
    )
    polymer_types: list[str] = Field(default_factory=list, description="Polymer types present")
    aa_groups: list[str] = Field(default_factory=list, description="Amino-acid groups present")

    @property
    def primary_metric_value(self) -> float:
        """Return mean chaperone fraction as the primary metric."""
        return self.mean_chaperone_fraction

    @property
    def primary_metric_sem(self) -> float:
        """Return SEM of chaperone fraction."""
        return self.sem_chaperone_fraction


# ---------------------------------------------------------------------------
# Comparison result
# ---------------------------------------------------------------------------


class ExposureComparisonResult(BaseComparisonResult[ExposureConditionSummary, PairwiseComparison]):
    """Complete exposure dynamics comparison result.

    This is the main output from ExposureDynamicsComparator.compare().
    Contains per-condition summaries of transient exposure and chaperone
    event statistics, plus pairwise statistical comparisons.

    Attributes
    ----------
    metric : str
        Always "chaperone_fraction".
    name : str
        Comparison project name.
    control_label : str, optional
        Label of the control condition.
    conditions : list[ExposureConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[PairwiseComparison]
        Pairwise t-tests on chaperone_fraction.
    anova : ANOVASummary, optional
        One-way ANOVA across all conditions.
    ranking : list[str]
        Condition labels sorted by chaperone_fraction (highest first).
    ranking_by_transient_fraction : list[str]
        Condition labels sorted by transient_fraction (highest first).
    excluded_conditions : list[str]
        Conditions excluded (e.g., no-polymer controls).
    equilibration_time : str
        Equilibration time used.
    created_at : datetime
        When the analysis was run.
    polyzymd_version : str
        Version of polyzymd used.
    """

    comparison_type: ClassVar[str] = "exposure"

    metric: str = "chaperone_fraction"
    ranking_by_transient_fraction: list[str] = Field(
        default_factory=list,
        description="Conditions ranked by transient fraction (highest first)",
    )
    excluded_conditions: list[str] = Field(
        default_factory=list,
        description="Conditions excluded from analysis",
    )
    equilibration_time: str = "0ns"
    created_at: datetime = Field(default_factory=datetime.now)
    polyzymd_version: str = __version__
