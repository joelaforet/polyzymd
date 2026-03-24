"""Result models for polymer-protein contacts comparison analysis.

This module defines Pydantic models for structured contacts comparison results
that can be serialized to JSON and used for downstream plotting.

Key features:
- Aggregate-level comparisons (coverage, mean contact fraction)
- Effect size (Cohen's d) for practical significance
- ANOVA for 3+ conditions
- Condition rankings by coverage and contact fraction

Note:
    Per-residue pairwise comparisons have been removed from this module.
    Contact data is mechanistic (explains WHY stability changes), not an
    observable. Per-residue contact-RMSF correlations are computed in
    `polyzymd compare report` which integrates contacts with RMSF data.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field

from polyzymd import __version__


class BindingPreferenceComparisonEntry(BaseModel):
    """Cross-condition binding preference comparison for one (polymer_type, protein_group) pair.

    Provides statistical comparison of enrichment ratios across conditions.

    Attributes
    ----------
    polymer_type : str
        Polymer residue type (e.g., "SBM", "EGM")
    protein_group : str
        Protein group label (e.g., "aromatic", "charged_positive")
    condition_values : dict[str, tuple[float, float]]
        Mapping of condition label to (mean_enrichment, sem_enrichment)
    pairwise_p_values : dict[str, float], optional
        P-values for pairwise comparisons (e.g., "A_vs_B": 0.05)
    highest_enrichment_condition : str
        Condition with highest mean enrichment
    lowest_enrichment_condition : str
        Condition with lowest mean enrichment
    """

    polymer_type: str
    protein_group: str
    condition_values: dict[str, tuple[float, float]] = Field(
        default_factory=dict, description="condition_label -> (mean, sem)"
    )
    pairwise_p_values: dict[str, float] = Field(
        default_factory=dict, description="Pairwise comparison p-values"
    )
    highest_enrichment_condition: Optional[str] = None
    lowest_enrichment_condition: Optional[str] = None


class BindingPreferenceComparisonSummary(BaseModel):
    """Summary of binding preference comparison across conditions.

    Contains per-condition enrichment matrices and cross-condition comparisons.

    Attributes
    ----------
    entries : list[BindingPreferenceComparisonEntry]
        Comparison data for each (polymer_type, protein_group) pair
    polymer_types : list[str]
        All polymer types found across conditions
    protein_groups : list[str]
        All protein groups analyzed
    n_conditions : int
        Number of conditions compared
    condition_labels : list[str]
        Labels of all conditions
    surface_exposure_threshold : float
        SASA threshold used for surface filtering
    """

    entries: list[BindingPreferenceComparisonEntry] = Field(default_factory=list)
    polymer_types: list[str] = Field(default_factory=list)
    protein_groups: list[str] = Field(default_factory=list)
    n_conditions: int = 0
    condition_labels: list[str] = Field(default_factory=list)
    surface_exposure_threshold: Optional[float] = None

    def get_entry(
        self, polymer_type: str, protein_group: str
    ) -> Optional[BindingPreferenceComparisonEntry]:
        """Get entry for a (polymer_type, protein_group) pair."""
        for entry in self.entries:
            if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                return entry
        return None

    def get_enrichment_matrix_for_condition(
        self, condition_label: str
    ) -> dict[str, dict[str, float]]:
        """Get enrichment matrix for a specific condition.

        Returns
        -------
        dict[str, dict[str, float]]
            {polymer_type: {protein_group: enrichment}}
        """
        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}
            values = entry.condition_values.get(condition_label)
            if values:
                result[entry.polymer_type][entry.protein_group] = values[0]  # mean
            else:
                result[entry.polymer_type][entry.protein_group] = 0.0
        return result


class AggregateComparisonResult(BaseModel):
    """Statistical comparison for aggregate metrics (coverage, mean contact fraction).

    Attributes
    ----------
    metric : str
        Name of the metric ("coverage" or "mean_contact_fraction")
    condition_a : str
        Label of condition A
    condition_b : str
        Label of condition B
    condition_a_mean : float
        Mean value for condition A
    condition_a_sem : float
        SEM for condition A
    condition_b_mean : float
        Mean value for condition B
    condition_b_sem : float
        SEM for condition B
    t_statistic : float
        T-test statistic
    p_value : float
        Two-tailed p-value
    cohens_d : float
        Effect size
    effect_size_interpretation : str
        "negligible", "small", "medium", or "large"
    significant : bool
        Whether p < alpha
    percent_change : float
        Percent change from A to B
    direction : str
        "increased", "decreased", or "unchanged"
    """

    metric: str
    condition_a: str
    condition_b: str
    condition_a_mean: float
    condition_a_sem: float
    condition_b_mean: float
    condition_b_sem: float
    t_statistic: float
    p_value: float
    cohens_d: float
    effect_size_interpretation: str
    significant: bool
    percent_change: float
    direction: str


class ContactsConditionSummary(BaseModel):
    """Summary statistics for one condition in a contacts comparison.

    Attributes
    ----------
    label : str
        Display name for this condition
    config_path : str
        Path to the simulation config file
    n_replicates : int
        Number of replicates included
    n_residues : int
        Number of protein residues
    coverage_mean : float
        Mean coverage fraction (residues contacted / total)
    coverage_sem : float
        SEM of coverage
    mean_contact_fraction : float
        Mean of mean contact fractions across replicates
    mean_contact_fraction_sem : float
        SEM of mean contact fraction
    residence_time_by_polymer_type : dict[str, tuple[float, float]]
        Mean ± SEM residence time in frames for each polymer type.
        Keys are polymer residue names (e.g., "SBM", "EGM").
        Values are (mean, sem) tuples.
    """

    label: str
    config_path: str
    n_replicates: int
    n_residues: int
    coverage_mean: float
    coverage_sem: float
    mean_contact_fraction: float
    mean_contact_fraction_sem: float
    residence_time_by_polymer_type: dict[str, tuple[float, float]] = Field(default_factory=dict)


class ContactsPairwiseComparison(BaseModel):
    """Pairwise comparison between two conditions for contacts analysis.

    Contains aggregate-level comparisons only (coverage, mean_contact_fraction).
    Per-residue contact-RMSF correlations are computed separately in the
    unified report (`polyzymd compare report`).

    Attributes
    ----------
    condition_a : str
        Label of first condition (typically control/reference)
    condition_b : str
        Label of second condition (typically treatment)
    aggregate_comparisons : list[AggregateComparisonResult]
        Aggregate-level comparisons (coverage, mean_contact_fraction)
    """

    condition_a: str
    condition_b: str
    aggregate_comparisons: list[AggregateComparisonResult]


class ContactsANOVASummary(BaseModel):
    """One-way ANOVA result summary for contacts comparison.

    Attributes
    ----------
    metric : str
        Metric being tested ("coverage" or "mean_contact_fraction")
    f_statistic : float
        F-statistic from ANOVA
    p_value : float
        P-value for the test
    significant : bool
        Whether p < 0.05
    """

    metric: str
    f_statistic: float
    p_value: float
    significant: bool


class ContactsComparisonResult(BaseModel):
    """Complete contacts comparison analysis result.

    This is the main output from ContactsComparator.compare().
    Contains all condition summaries, aggregate statistical comparisons,
    and rankings.

    Attributes
    ----------
    name : str
        Name of the comparison project
    contacts_name : str
        Name of the contacts analysis
    contacts_description : str, optional
        Description of the contacts analysis
    polymer_selection : str
        Polymer selection string used
    protein_selection : str
        Protein selection string used
    cutoff : float
        Contact cutoff distance in Angstroms
    contact_criteria : str
        Contact criteria used
    fdr_alpha : float
        FDR alpha used for Benjamini-Hochberg correction
    control_label : str, optional
        Label of the control condition
    conditions : list[ContactsConditionSummary]
        Summary for each condition
    pairwise_comparisons : list[ContactsPairwiseComparison]
        Statistical comparisons (all vs control, or all pairs)
    anova : list[ContactsANOVASummary]
        ANOVA results if 3+ conditions (one per aggregate metric)
    ranking_by_coverage : list[str]
        Labels sorted by coverage (descending)
    ranking_by_contact_fraction : list[str]
        Labels sorted by mean contact fraction (descending)
    equilibration_time : str
        Equilibration time used
    created_at : datetime
        When the analysis was run
    polyzymd_version : str
        Version of polyzymd used
    """

    name: str
    contacts_name: str
    contacts_description: Optional[str] = None
    polymer_selection: str
    protein_selection: str
    cutoff: float
    contact_criteria: str
    fdr_alpha: float
    control_label: Optional[str] = None
    conditions: list[ContactsConditionSummary]
    pairwise_comparisons: list[ContactsPairwiseComparison]
    anova: list[ContactsANOVASummary] = Field(default_factory=list)
    ranking_by_coverage: list[str] = Field(default_factory=list)
    ranking_by_contact_fraction: list[str] = Field(default_factory=list)
    excluded_conditions: list[str] = Field(
        default_factory=list,
        description="Conditions excluded from analysis (no polymer atoms found)",
    )
    binding_preference: Optional[BindingPreferenceComparisonSummary] = Field(
        default=None,
        description="Binding preference comparison across conditions (if computed)",
    )
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
    def load(cls, path: Path | str) -> "ContactsComparisonResult":
        """Load result from JSON file.

        Parameters
        ----------
        path : Path or str
            Path to JSON file

        Returns
        -------
        ContactsComparisonResult
            Loaded result
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())

    def get_condition(self, label: str) -> ContactsConditionSummary:
        """Get a condition by label.

        Parameters
        ----------
        label : str
            Condition label

        Returns
        -------
        ContactsConditionSummary
            The matching condition

        Raises
        ------
        KeyError
            If condition not found
        """
        for cond in self.conditions:
            if cond.label == label:
                return cond
        raise KeyError(f"Condition '{label}' not found")

    def get_comparison(self, label: str) -> Optional[ContactsPairwiseComparison]:
        """Get pairwise comparison for a condition vs control.

        Parameters
        ----------
        label : str
            Treatment condition label

        Returns
        -------
        ContactsPairwiseComparison or None
            The comparison, or None if not found
        """
        for comp in self.pairwise_comparisons:
            if comp.condition_b == label:
                return comp
        return None
