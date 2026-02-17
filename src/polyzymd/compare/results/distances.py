"""Result models for distance comparison analysis.

This module defines Pydantic models for structured distance comparison results
that can be serialized to JSON and used for downstream plotting.

The primary ranking metric is mean distance (lower = closer interactions).
Secondary metric is fraction below threshold (if threshold is specified).
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field

from polyzymd import __version__


class DistancePairSummary(BaseModel):
    """Summary statistics for one distance pair across replicates.

    Attributes
    ----------
    label : str
        Human-readable label for this pair (e.g., "Ser77-Substrate").
    selection_a : str
        MDAnalysis selection for first atom/group.
    selection_b : str
        MDAnalysis selection for second atom/group.
    mean_distance : float
        Mean distance in Angstroms (averaged across replicates).
    sem_distance : float
        SEM of distance across replicates.
    fraction_below_threshold : float, optional
        Mean fraction of frames below threshold (if threshold specified).
    sem_fraction_below : float, optional
        SEM of fraction below threshold.
    per_replicate_means : list[float]
        Mean distance from each replicate.
    """

    label: str
    selection_a: str
    selection_b: str
    mean_distance: float
    sem_distance: float
    fraction_below_threshold: Optional[float] = None
    sem_fraction_below: Optional[float] = None
    per_replicate_means: list[float] = Field(default_factory=list)


class DistanceConditionSummary(BaseModel):
    """Summary statistics for one condition in a distance comparison.

    Attributes
    ----------
    label : str
        Display name for this condition.
    config_path : str
        Path to the simulation config file.
    n_replicates : int
        Number of replicates included.
    pair_summaries : list[DistancePairSummary]
        Summary for each distance pair.
    overall_mean_distance : float
        Average mean distance across all pairs (primary metric).
    overall_sem_distance : float
        SEM of overall mean distance.
    overall_fraction_below : float, optional
        Average fraction below threshold across all pairs.
    overall_sem_fraction_below : float, optional
        SEM of overall fraction below threshold.
    replicate_values : list[float]
        Per-replicate overall mean distances (for statistical tests).
    replicate_fractions : list[float], optional
        Per-replicate overall fractions (for statistical tests).
    """

    label: str
    config_path: str
    n_replicates: int
    pair_summaries: list[DistancePairSummary]
    overall_mean_distance: float
    overall_sem_distance: float
    overall_fraction_below: Optional[float] = None
    overall_sem_fraction_below: Optional[float] = None
    replicate_values: list[float] = Field(default_factory=list)
    replicate_fractions: Optional[list[float]] = None


class DistancePairwiseComparison(BaseModel):
    """Statistical comparison between two conditions for distance analysis.

    Includes comparisons for both mean distance and fraction metrics.

    Attributes
    ----------
    condition_a : str
        Label of first condition (typically control).
    condition_b : str
        Label of second condition (typically treatment).

    # Distance metric comparisons
    distance_t_statistic : float
        T-test statistic for mean distance.
    distance_p_value : float
        Two-tailed p-value for distance comparison.
    distance_cohens_d : float
        Effect size for distance (negative = condition_b has lower distance).
    distance_effect_interpretation : str
        "negligible", "small", "medium", or "large".
    distance_direction : str
        "closer" (lower distance), "farther" (higher distance), or "unchanged".
    distance_significant : bool
        Whether p < 0.05 for distance.
    distance_percent_change : float
        Percent change in mean distance.

    # Fraction metric comparisons (optional, if threshold specified)
    fraction_t_statistic : float, optional
        T-test statistic for fraction below threshold.
    fraction_p_value : float, optional
        P-value for fraction comparison.
    fraction_cohens_d : float, optional
        Effect size for fraction.
    fraction_effect_interpretation : str, optional
        Effect size interpretation.
    fraction_direction : str, optional
        "more_contact", "less_contact", or "unchanged".
    fraction_significant : bool, optional
        Whether p < 0.05 for fraction.
    fraction_percent_change : float, optional
        Percent change in fraction below threshold.
    """

    condition_a: str
    condition_b: str

    # Distance metric
    distance_t_statistic: float
    distance_p_value: float
    distance_cohens_d: float
    distance_effect_interpretation: str
    distance_direction: str
    distance_significant: bool
    distance_percent_change: float

    # Fraction metric (optional)
    fraction_t_statistic: Optional[float] = None
    fraction_p_value: Optional[float] = None
    fraction_cohens_d: Optional[float] = None
    fraction_effect_interpretation: Optional[str] = None
    fraction_direction: Optional[str] = None
    fraction_significant: Optional[bool] = None
    fraction_percent_change: Optional[float] = None


class DistanceANOVASummary(BaseModel):
    """ANOVA result summary for distance comparison.

    Includes results for both metrics.

    Attributes
    ----------
    distance_f_statistic : float
        F-statistic for mean distance.
    distance_p_value : float
        P-value for distance ANOVA.
    distance_significant : bool
        Whether p < 0.05.
    fraction_f_statistic : float, optional
        F-statistic for fraction below threshold.
    fraction_p_value : float, optional
        P-value for fraction ANOVA.
    fraction_significant : bool, optional
        Whether fraction p < 0.05.
    """

    distance_f_statistic: float
    distance_p_value: float
    distance_significant: bool
    fraction_f_statistic: Optional[float] = None
    fraction_p_value: Optional[float] = None
    fraction_significant: Optional[bool] = None


class DistanceComparisonResult(BaseModel):
    """Complete distance comparison analysis result.

    This is the main output from DistancesComparator.compare().
    Contains all condition summaries, statistical comparisons,
    and rankings for both metrics.

    Attributes
    ----------
    metric : str
        Primary metric: "mean_distance".
    name : str
        Name of the comparison project.
    n_pairs : int
        Number of distance pairs analyzed.
    pair_labels : list[str]
        Labels for each pair.
    threshold : float, optional
        Distance threshold used (if any).
    control_label : str, optional
        Label of the control condition.
    conditions : list[DistanceConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[DistancePairwiseComparison]
        Statistical comparisons.
    anova : DistanceANOVASummary, optional
        ANOVA result if 3+ conditions.
    ranking : list[str]
        Labels sorted by mean distance (ascending = lowest first).
    ranking_by_fraction : list[str], optional
        Labels sorted by fraction below threshold (descending = highest first).
    equilibration_time : str
        Equilibration time used.
    created_at : datetime
        When the analysis was run.
    polyzymd_version : str
        Version of polyzymd used.
    """

    metric: str = "mean_distance"
    name: str
    n_pairs: int
    pair_labels: list[str]
    threshold: Optional[float] = None
    control_label: Optional[str] = None
    conditions: list[DistanceConditionSummary]
    pairwise_comparisons: list[DistancePairwiseComparison]
    anova: Optional[DistanceANOVASummary] = None
    ranking: list[str]
    ranking_by_fraction: Optional[list[str]] = None
    equilibration_time: str
    created_at: datetime
    polyzymd_version: str = __version__

    def save(self, path: Path | str) -> Path:
        """Save result to JSON file.

        Parameters
        ----------
        path : Path or str
            Output path.

        Returns
        -------
        Path
            Path to saved file.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path | str) -> "DistanceComparisonResult":
        """Load result from JSON file.

        Parameters
        ----------
        path : Path or str
            Path to JSON file.

        Returns
        -------
        DistanceComparisonResult
            Loaded result.
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())

    def get_condition(self, label: str) -> DistanceConditionSummary:
        """Get a condition by label.

        Parameters
        ----------
        label : str
            Condition label.

        Returns
        -------
        DistanceConditionSummary
            The matching condition.
        """
        for cond in self.conditions:
            if cond.label == label:
                return cond
        raise KeyError(f"Condition '{label}' not found")

    def get_comparison(self, label: str) -> Optional[DistancePairwiseComparison]:
        """Get pairwise comparison for a condition vs control.

        Parameters
        ----------
        label : str
            Treatment condition label.

        Returns
        -------
        DistancePairwiseComparison or None
            The comparison, or None if not found.
        """
        for comp in self.pairwise_comparisons:
            if comp.condition_b == label:
                return comp
        return None
