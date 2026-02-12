"""Result models for catalytic triad comparison analysis.

This module defines Pydantic models for structured triad comparison results
that can be serialized to JSON and used for downstream plotting.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

from pydantic import BaseModel

from polyzymd import __version__


class TriadPairSummary(BaseModel):
    """Summary statistics for one distance pair in a triad.

    Attributes
    ----------
    label : str
        Label for this pair (e.g., "Asp133-His156")
    mean_distance : float
        Mean distance in Angstroms
    sem_distance : float
        SEM of distance
    mean_fraction_below : float, optional
        Mean fraction of frames below threshold
    sem_fraction_below : float, optional
        SEM of fraction below threshold
    """

    label: str
    mean_distance: float
    sem_distance: float
    mean_fraction_below: Optional[float] = None
    sem_fraction_below: Optional[float] = None


class TriadConditionSummary(BaseModel):
    """Summary statistics for one condition in a triad comparison.

    Attributes
    ----------
    label : str
        Display name for this condition
    config_path : str
        Path to the simulation config file
    n_replicates : int
        Number of replicates included
    mean_simultaneous_contact : float
        Mean simultaneous contact fraction across replicates (0-1)
    sem_simultaneous_contact : float
        Standard error of the mean
    replicate_values : list[float]
        Per-replicate simultaneous contact fractions
    pair_summaries : list[TriadPairSummary]
        Summary for each distance pair in the triad
    """

    label: str
    config_path: str
    n_replicates: int
    mean_simultaneous_contact: float
    sem_simultaneous_contact: float
    replicate_values: list[float]
    pair_summaries: list[TriadPairSummary]


class TriadPairwiseComparison(BaseModel):
    """Statistical comparison between two conditions for triad analysis.

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
        Effect size (positive = condition_b > condition_a)
    effect_size_interpretation : str
        "negligible", "small", "medium", or "large"
    direction : str
        "improving" (more contact), "worsening" (less contact), or "unchanged"
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


class TriadANOVASummary(BaseModel):
    """One-way ANOVA result summary for triad comparison.

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


class TriadComparisonResult(BaseModel):
    """Complete triad comparison analysis result.

    This is the main output from TriadComparator.compare().
    Contains all condition summaries, statistical comparisons,
    and rankings.

    Attributes
    ----------
    metric : str
        The metric being compared ("simultaneous_contact_fraction")
    name : str
        Name of the comparison project
    triad_name : str
        Name of the catalytic triad
    triad_description : str, optional
        Description of the triad
    threshold : float
        Contact threshold in Angstroms
    n_pairs : int
        Number of distance pairs in the triad
    pair_labels : list[str]
        Labels for each pair
    control_label : str, optional
        Label of the control condition
    conditions : list[TriadConditionSummary]
        Summary for each condition
    pairwise_comparisons : list[TriadPairwiseComparison]
        Statistical comparisons (all vs control, or all pairs)
    anova : TriadANOVASummary, optional
        ANOVA result if 3+ conditions
    ranking : list[str]
        Labels sorted by metric (descending = highest contact first)
    equilibration_time : str
        Equilibration time used
    created_at : datetime
        When the analysis was run
    polyzymd_version : str
        Version of polyzymd used
    """

    metric: str
    name: str
    triad_name: str
    triad_description: Optional[str] = None
    threshold: float
    n_pairs: int
    pair_labels: list[str]
    control_label: Optional[str] = None
    conditions: list[TriadConditionSummary]
    pairwise_comparisons: list[TriadPairwiseComparison]
    anova: Optional[TriadANOVASummary] = None
    ranking: list[str]
    equilibration_time: str
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
    def load(cls, path: Path | str) -> "TriadComparisonResult":
        """Load result from JSON file.

        Parameters
        ----------
        path : Path or str
            Path to JSON file

        Returns
        -------
        TriadComparisonResult
            Loaded result
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())

    def get_condition(self, label: str) -> TriadConditionSummary:
        """Get a condition by label.

        Parameters
        ----------
        label : str
            Condition label

        Returns
        -------
        TriadConditionSummary
            The matching condition
        """
        for cond in self.conditions:
            if cond.label == label:
                return cond
        raise KeyError(f"Condition '{label}' not found")

    def get_comparison(self, label: str) -> Optional[TriadPairwiseComparison]:
        """Get pairwise comparison for a condition vs control.

        Parameters
        ----------
        label : str
            Treatment condition label

        Returns
        -------
        TriadPairwiseComparison or None
            The comparison, or None if not found
        """
        for comp in self.pairwise_comparisons:
            if comp.condition_b == label:
                return comp
        return None
