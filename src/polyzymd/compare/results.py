"""Result models for comparison analysis.

This module defines Pydantic models for structured comparison results
that can be serialized to JSON and used for downstream plotting.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

from pydantic import BaseModel

from polyzymd import __version__


class ConditionSummary(BaseModel):
    """Summary statistics for one condition in a comparison.

    Attributes
    ----------
    label : str
        Display name for this condition
    config_path : str
        Path to the simulation config file
    n_replicates : int
        Number of replicates included
    mean_rmsf : float
        Mean RMSF across replicates (Angstroms)
    sem_rmsf : float
        Standard error of the mean
    replicate_values : list[float]
        Per-replicate mean RMSF values
    """

    label: str
    config_path: str
    n_replicates: int
    mean_rmsf: float
    sem_rmsf: float
    replicate_values: list[float]


class PairwiseComparison(BaseModel):
    """Statistical comparison between two conditions.

    Attributes
    ----------
    condition_a : str
        Label of first condition (typically control)
    condition_b : str
        Label of second condition (typically treatment)
    t_statistic : float
        T-test statistic
    p_value : float
        Two-tailed p-value
    cohens_d : float
        Effect size (positive = condition_a > condition_b)
    effect_size_interpretation : str
        "negligible", "small", "medium", or "large"
    direction : str
        For RMSF: "stabilizing" or "destabilizing"
    significant : bool
        Whether p < 0.05
    percent_change : float
        Percent change from condition_a to condition_b
    """

    condition_a: str
    condition_b: str
    t_statistic: float
    p_value: float
    cohens_d: float
    effect_size_interpretation: str
    direction: str
    significant: bool
    percent_change: float


class ANOVASummary(BaseModel):
    """One-way ANOVA result summary.

    Attributes
    ----------
    f_statistic : float
        F-statistic from ANOVA
    p_value : float
        P-value for the test
    significant : bool
        Whether p < 0.05
    """

    f_statistic: float
    p_value: float
    significant: bool


class ComparisonResult(BaseModel):
    """Complete comparison analysis result.

    This is the main output from RMSFComparator.compare().
    Contains all condition summaries, statistical comparisons,
    and rankings.

    Attributes
    ----------
    metric : str
        The metric being compared (e.g., "rmsf")
    name : str
        Name of the comparison project
    control_label : str, optional
        Label of the control condition
    conditions : list[ConditionSummary]
        Summary for each condition
    pairwise_comparisons : list[PairwiseComparison]
        Statistical comparisons (all vs control, or all pairs)
    anova : ANOVASummary, optional
        ANOVA result if 3+ conditions
    ranking : list[str]
        Labels sorted by metric (ascending for RMSF = lowest first)
    equilibration_time : str
        Equilibration time used
    selection : str
        Atom selection used
    created_at : datetime
        When the analysis was run
    polyzymd_version : str
        Version of polyzymd used
    """

    metric: str
    name: str
    control_label: Optional[str] = None
    conditions: list[ConditionSummary]
    pairwise_comparisons: list[PairwiseComparison]
    anova: Optional[ANOVASummary] = None
    ranking: list[str]
    equilibration_time: str
    selection: str
    created_at: datetime
    polyzymd_version: str

    def save(self, path: Path | str) -> Path:
        """Save result to JSON file.

        Parameters
        ----------
        path : Path or str
            Output path

        Returns
        -------
        Path
            Path to saved file
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path | str) -> "ComparisonResult":
        """Load result from JSON file.

        Parameters
        ----------
        path : Path or str
            Path to JSON file

        Returns
        -------
        ComparisonResult
            Loaded result
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())

    def get_condition(self, label: str) -> ConditionSummary:
        """Get a condition by label.

        Parameters
        ----------
        label : str
            Condition label

        Returns
        -------
        ConditionSummary
            The matching condition
        """
        for cond in self.conditions:
            if cond.label == label:
                return cond
        raise KeyError(f"Condition '{label}' not found")

    def get_comparison(self, label: str) -> Optional[PairwiseComparison]:
        """Get pairwise comparison for a condition vs control.

        Parameters
        ----------
        label : str
            Treatment condition label

        Returns
        -------
        PairwiseComparison or None
            The comparison, or None if not found
        """
        for comp in self.pairwise_comparisons:
            if comp.condition_b == label:
                return comp
        return None
