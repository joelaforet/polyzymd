"""Result models for distance comparison analysis.

This module defines Pydantic models for structured distance comparison results
that can be serialized to JSON and used for downstream plotting.

Each distance pair is compared independently across conditions - there is no
cross-pair averaging since different pairs measure fundamentally different
physical quantities (e.g., H-bond distances vs lid-opening distances).

For each pair:
- Primary ranking: mean distance (lower = closer = better)
- Secondary ranking: fraction below threshold (higher = more contact = better)
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
    threshold : float, optional
        Distance threshold used for this pair (Angstroms).
    mean_distance : float
        Mean distance in Angstroms (averaged across replicates).
    sem_distance : float
        SEM of distance across replicates.
    fraction_below_threshold : float, optional
        Mean fraction of frames below threshold (if threshold specified).
    sem_fraction_below : float, optional
        SEM of fraction below threshold.
    per_replicate_means : list[float]
        Mean distance from each replicate (for statistical tests).
    per_replicate_fractions : list[float], optional
        Fraction below threshold from each replicate.
    """

    label: str
    selection_a: str
    selection_b: str
    threshold: Optional[float] = None
    mean_distance: float
    sem_distance: float
    fraction_below_threshold: Optional[float] = None
    sem_fraction_below: Optional[float] = None
    per_replicate_means: list[float] = Field(default_factory=list)
    per_replicate_fractions: Optional[list[float]] = None


class DistanceConditionSummary(BaseModel):
    """Summary statistics for one condition in a distance comparison.

    Note: There is no "overall" distance metric across pairs, since averaging
    unrelated distances (e.g., H-bond + lid-opening) is not meaningful.
    Each pair is compared independently.

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
    """

    label: str
    config_path: str
    n_replicates: int
    pair_summaries: list[DistancePairSummary]

    def get_pair(self, label: str) -> DistancePairSummary:
        """Get a pair summary by label.

        Parameters
        ----------
        label : str
            Pair label.

        Returns
        -------
        DistancePairSummary
            The matching pair summary.
        """
        for ps in self.pair_summaries:
            if ps.label == label:
                return ps
        raise KeyError(f"Pair '{label}' not found in condition '{self.label}'")


class DistancePairwiseComparison(BaseModel):
    """Statistical comparison between two conditions for a single distance pair.

    Each comparison is specific to one pair - cross-pair comparisons are not
    meaningful since different pairs measure different physical quantities.

    Attributes
    ----------
    pair_label : str
        Label of the distance pair being compared.
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

    pair_label: str
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


class DistancePairANOVA(BaseModel):
    """ANOVA result for a single distance pair.

    Attributes
    ----------
    pair_label : str
        Label of the distance pair.
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

    pair_label: str
    distance_f_statistic: float
    distance_p_value: float
    distance_significant: bool
    fraction_f_statistic: Optional[float] = None
    fraction_p_value: Optional[float] = None
    fraction_significant: Optional[bool] = None


class DistanceComparisonResult(BaseModel):
    """Complete distance comparison analysis result.

    This is the main output from DistancesComparator.compare().

    Each distance pair is compared independently - rankings and statistics
    are computed per-pair since averaging unrelated distances is not meaningful.

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
    control_label : str, optional
        Label of the control condition.
    conditions : list[DistanceConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[DistancePairwiseComparison]
        Statistical comparisons (grouped by pair).
    anova_by_pair : list[DistancePairANOVA], optional
        ANOVA results for each pair (if 3+ conditions).
    ranking_by_pair : dict[str, list[str]]
        Condition labels sorted by mean distance, keyed by pair label.
        Lower distance = better = earlier in list.
    fraction_ranking_by_pair : dict[str, list[str]], optional
        Condition labels sorted by fraction below threshold, keyed by pair.
        Higher fraction = better = earlier in list.
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
    control_label: Optional[str] = None
    conditions: list[DistanceConditionSummary]
    pairwise_comparisons: list[DistancePairwiseComparison]
    anova_by_pair: Optional[list[DistancePairANOVA]] = None
    ranking_by_pair: dict[str, list[str]]
    fraction_ranking_by_pair: Optional[dict[str, list[str]]] = None
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

    def get_comparisons_for_pair(self, pair_label: str) -> list[DistancePairwiseComparison]:
        """Get all pairwise comparisons for a specific pair.

        Parameters
        ----------
        pair_label : str
            The pair label.

        Returns
        -------
        list[DistancePairwiseComparison]
            Comparisons for this pair.
        """
        return [c for c in self.pairwise_comparisons if c.pair_label == pair_label]

    def get_comparison(
        self, pair_label: str, condition_b: str
    ) -> Optional[DistancePairwiseComparison]:
        """Get pairwise comparison for a specific pair and condition vs control.

        Parameters
        ----------
        pair_label : str
            The pair label.
        condition_b : str
            Treatment condition label.

        Returns
        -------
        DistancePairwiseComparison or None
            The comparison, or None if not found.
        """
        for comp in self.pairwise_comparisons:
            if comp.pair_label == pair_label and comp.condition_b == condition_b:
                return comp
        return None

    def get_ranking(self, pair_label: str) -> list[str]:
        """Get ranking for a specific pair.

        Parameters
        ----------
        pair_label : str
            The pair label.

        Returns
        -------
        list[str]
            Condition labels sorted by mean distance (lowest first).
        """
        return self.ranking_by_pair.get(pair_label, [])

    def get_fraction_ranking(self, pair_label: str) -> list[str]:
        """Get fraction ranking for a specific pair.

        Parameters
        ----------
        pair_label : str
            The pair label.

        Returns
        -------
        list[str]
            Condition labels sorted by fraction below threshold (highest first).
        """
        if self.fraction_ranking_by_pair is None:
            return []
        return self.fraction_ranking_by_pair.get(pair_label, [])
