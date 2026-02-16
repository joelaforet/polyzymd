"""RMSF condition summary and comparison result models.

These classes inherit from the base classes in compare/core/ and add
RMSF-specific fields.
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


class RMSFConditionSummary(BaseConditionSummary):
    """Summary statistics for one condition in an RMSF comparison.

    Attributes
    ----------
    label : str
        Display name for this condition.
    config_path : str
        Path to the simulation config file.
    n_replicates : int
        Number of replicates included.
    mean_rmsf : float
        Mean RMSF across replicates (Angstroms).
    sem_rmsf : float
        Standard error of the mean.
    replicate_values : list[float]
        Per-replicate mean RMSF values.
    """

    mean_rmsf: float = Field(..., description="Mean RMSF in Angstroms")
    sem_rmsf: float = Field(..., description="Standard error of mean RMSF")

    @property
    def primary_metric_value(self) -> float:
        """Return mean RMSF as the primary metric for ranking."""
        return self.mean_rmsf

    @property
    def primary_metric_sem(self) -> float:
        """Return SEM of RMSF."""
        return self.sem_rmsf


class RMSFComparisonResult(BaseComparisonResult[RMSFConditionSummary, PairwiseComparison]):
    """Complete RMSF comparison analysis result.

    This is the main output from RMSFComparator.compare().
    Contains all condition summaries, statistical comparisons,
    and rankings.

    Attributes
    ----------
    metric : str
        Always "rmsf".
    name : str
        Name of the comparison project.
    control_label : str, optional
        Label of the control condition.
    conditions : list[RMSFConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[PairwiseComparison]
        Statistical comparisons (all vs control, or all pairs).
    anova : ANOVASummary, optional
        ANOVA result if 3+ conditions.
    ranking : list[str]
        Labels sorted by RMSF (ascending = lowest first = most stable).
    equilibration_time : str
        Equilibration time used.
    selection : str
        Atom selection used for RMSF calculation.
    created_at : datetime
        When the analysis was run.
    polyzymd_version : str
        Version of polyzymd used.
    """

    comparison_type: ClassVar[str] = "rmsf"

    # Override with specific types
    conditions: list[RMSFConditionSummary]
    pairwise_comparisons: list[PairwiseComparison]
    anova: ANOVASummary | None = None

    # RMSF-specific field
    selection: str = Field(..., description="Atom selection used for RMSF")
