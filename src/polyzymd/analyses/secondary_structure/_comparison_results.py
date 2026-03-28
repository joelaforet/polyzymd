"""Secondary structure condition summary and comparison result models.

These classes follow the RMSF comparison result pattern, inheriting from
the base classes in ``compare/core/``.

The primary ranking metric is **helix fraction** (most relevant for
thermal stability studies -- helix loss is a hallmark of unfolding).
Strand and coil fractions are included as additional fields.

.. versionchanged:: 1.3.0
    Migrated from ``polyzymd.compare.results.secondary_structure`` into the
    secondary structure analysis package so that contributors find all
    secondary-structure-related code in one place.
"""

from __future__ import annotations

from datetime import datetime
from typing import ClassVar

from pydantic import Field

from polyzymd import __version__
from polyzymd.compare.core.base import (
    ANOVASummary,
    BaseComparisonResult,
    BaseConditionSummary,
    PairwiseComparison,
)


class SSConditionSummary(BaseConditionSummary):
    """Summary statistics for one condition in an SS comparison.

    Attributes
    ----------
    label : str
        Display name for this condition.
    config_path : str
        Path to the simulation config file.
    n_replicates : int
        Number of replicates included.
    replicate_values : list[float]
        Per-replicate overall helix fractions (primary metric for stats).
    mean_helix : float
        Mean helix fraction across replicates.
    sem_helix : float
        SEM of helix fraction.
    mean_strand : float
        Mean strand fraction across replicates.
    sem_strand : float
        SEM of strand fraction.
    mean_coil : float
        Mean coil fraction across replicates.
    sem_coil : float
        SEM of coil fraction.
    per_replicate_helix : list[float]
        Per-replicate helix fractions.
    per_replicate_strand : list[float]
        Per-replicate strand fractions.
    per_replicate_coil : list[float]
        Per-replicate coil fractions.
    """

    mean_helix: float = Field(..., description="Mean helix fraction")
    sem_helix: float = Field(..., description="SEM of helix fraction")
    mean_strand: float = Field(..., description="Mean strand fraction")
    sem_strand: float = Field(..., description="SEM of strand fraction")
    mean_coil: float = Field(..., description="Mean coil fraction")
    sem_coil: float = Field(..., description="SEM of coil fraction")

    # Per-replicate values for all three states (for multi-metric stats)
    per_replicate_helix: list[float] = Field(default_factory=list)
    per_replicate_strand: list[float] = Field(default_factory=list)
    per_replicate_coil: list[float] = Field(default_factory=list)

    @property
    def primary_metric_value(self) -> float:
        """Return mean helix fraction as the primary metric for ranking."""
        return self.mean_helix

    @property
    def primary_metric_sem(self) -> float:
        """Return SEM of helix fraction."""
        return self.sem_helix


class SSComparisonResult(BaseComparisonResult[SSConditionSummary, PairwiseComparison]):
    """Complete secondary structure comparison analysis result.

    This is the main output from ``SecondaryStructureComparator.compare()``.
    Contains condition summaries with helix/strand/coil fractions,
    pairwise statistical comparisons (on helix fraction), and rankings.

    The primary metric is helix fraction because helix loss is the
    most common signature of thermal unfolding in globular enzymes.

    Attributes
    ----------
    comparison_type : str
        Always ``"secondary_structure"``.
    metric : str
        Always ``"helix_fraction"``.
    name : str
        Name of the comparison project.
    control_label : str, optional
        Label of the control condition.
    conditions : list[SSConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[PairwiseComparison]
        Statistical comparisons (on helix fraction).
    anova : ANOVASummary, optional
        ANOVA result if 3+ conditions.
    ranking : list[str]
        Labels sorted by helix fraction (descending = most helix first).
    equilibration_time : str
        Equilibration time used.
    created_at : datetime
        When the analysis was run.
    polyzymd_version : str
        Version of polyzymd used.
    """

    comparison_type: ClassVar[str] = "secondary_structure"

    # Override with specific types
    conditions: list[SSConditionSummary]
    pairwise_comparisons: list[PairwiseComparison]
    anova: ANOVASummary | None = None
