"""Result models for Rg comparison analysis.

This module defines Pydantic models for structured Rg comparison results
that can be serialized to JSON and used for downstream plotting.

Each Rg run is compared independently across conditions — there is no
cross-run averaging since different runs measure fundamentally different
selections (e.g., protein backbone vs polymer core).

For each run:
- Primary ranking: mean Rg (lower = more compact = better)
"""

from __future__ import annotations

from typing import ClassVar

from pydantic import BaseModel, ConfigDict, Field

from polyzymd.analyses.base import BaseComparisonResult, BaseConditionSummary


class RgRunSummary(BaseModel):
    """Summary statistics for one Rg run across replicates.

    Attributes
    ----------
    label : str
        Human-readable label for this run (e.g., "protein_backbone").
    selection : str
        MDAnalysis selection used for Rg calculation.
    mean_rg : float
        Mean Rg in Angstroms (averaged across replicates).
    sem_rg : float
        SEM of Rg across replicates.
    per_replicate_means : list[float]
        Mean Rg from each replicate (for statistical tests).
    calculation_mode : str
        Calculation mode used.
    fragment_weighting : str | None
        Fragment weighting scheme used, if applicable.
    mean_fragments_per_frame : float | None
        Mean number of fragments per frame, if applicable.
    """

    label: str
    selection: str
    mean_rg: float
    sem_rg: float
    per_replicate_means: list[float] = Field(default_factory=list)
    calculation_mode: str = Field(default="selection", description="Calculation mode used")
    fragment_weighting: str | None = Field(default=None, description="Fragment weighting scheme")
    mean_fragments_per_frame: float | None = Field(default=None, description="Mean fragments/frame")


class RgConditionSummary(BaseConditionSummary):
    """Summary statistics for one condition in an Rg comparison.

    Note: There is no "overall" Rg metric across runs, since averaging
    Rg from different selections is not meaningful. Each run is compared
    independently.

    Inherits from BaseConditionSummary for interface consistency. The
    ``replicate_values`` field defaults to an empty list because Rg
    has no single primary metric — each run is compared independently.

    Attributes
    ----------
    label : str
        Display name for this condition.
    config_path : str
        Path to the simulation config file.
    n_replicates : int
        Number of replicates included.
    replicate_values : list[float]
        Empty by default (no single primary metric for multi-run Rg).
    run_summaries : list[RgRunSummary]
        Summary for each Rg run.
    """

    replicate_values: list[float] = Field(default_factory=list)
    run_summaries: list[RgRunSummary]

    @property
    def primary_metric_value(self) -> float:
        """Return 0.0 — multi-run Rg has no single primary metric."""
        return 0.0

    @property
    def primary_metric_sem(self) -> float:
        """Return 0.0 — multi-run Rg has no single primary metric."""
        return 0.0

    def get_run(self, label: str) -> RgRunSummary:
        """Get a run summary by label.

        Parameters
        ----------
        label : str
            Run label.

        Returns
        -------
        RgRunSummary
            The matching run summary.

        Raises
        ------
        KeyError
            If run not found.
        """
        for rs in self.run_summaries:
            if rs.label == label:
                return rs
        raise KeyError(f"Run '{label}' not found in condition '{self.label}'")


class RgRunPairwiseComparison(BaseModel):
    """Statistical comparison between two conditions for a single Rg run.

    Attributes
    ----------
    run_label : str
        Label of the Rg run being compared.
    condition_a : str
        Label of first condition (typically control).
    condition_b : str
        Label of second condition (typically treatment).
    t_statistic : float
        T-test statistic for mean Rg.
    p_value : float
        Two-tailed p-value.
    p_value_adjusted : float | None
        FDR-adjusted p-value (Benjamini-Hochberg), if correction was applied.
    cohens_d : float
        Effect size (Cohen's d).
    effect_interpretation : str
        "negligible", "small", "medium", or "large".
    direction : str
        "compaction" (lower Rg), "expansion" (higher Rg),
        or "unchanged".
    significant : bool
        Whether the comparison is statistically significant (uses adjusted p-value after FDR correction).
    percent_change : float
        Percent change in mean Rg.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    run_label: str
    condition_a: str
    condition_b: str
    t_statistic: float
    p_value: float
    p_value_adjusted: float | None = None
    cohens_d: float
    effect_interpretation: str
    direction: str
    significant: bool
    percent_change: float


class RgRunANOVA(BaseModel):
    """ANOVA result for a single Rg run.

    Attributes
    ----------
    run_label : str
        Label of the Rg run.
    f_statistic : float
        F-statistic for mean Rg.
    p_value : float
        P-value for ANOVA.
    significant : bool
        Whether p < 0.05.
    """

    run_label: str
    f_statistic: float
    p_value: float
    significant: bool


class RgComparisonResult(BaseComparisonResult[RgConditionSummary, RgRunPairwiseComparison]):
    """Complete Rg comparison analysis result.

    Inherits save/load/get_condition from BaseComparisonResult.

    Each Rg run is compared independently — rankings and statistics
    are computed per-run since averaging Rg from different selections
    is not meaningful.

    Attributes
    ----------
    metric : str
        Primary metric: "mean_rg".
    name : str
        Name of the comparison project.
    n_runs : int
        Number of Rg runs analyzed.
    run_labels : list[str]
        Labels for each run.
    control_label : str | None
        Label of the control condition.
    fdr_alpha : float | None
        False discovery rate threshold used for Benjamini-Hochberg correction.
    conditions : list[RgConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[RgRunPairwiseComparison]
        Statistical comparisons (grouped by run).
    anova_by_run : list[RgRunANOVA] | None
        ANOVA results for each run (if 3+ conditions).
    ranking_by_run : dict[str, list[str]]
        Condition labels sorted by mean Rg, keyed by run label.
        Lower Rg = more compact = earlier in list.
    ranking : list[str]
        Empty by default (Rg ranks per-run via ranking_by_run).
    equilibration_time : str
        Equilibration time used.
    """

    comparison_type: ClassVar[str] = "rg"
    model_config = ConfigDict(ser_json_inf_nan="strings")

    metric: str = "mean_rg"
    name: str
    n_runs: int
    run_labels: list[str]
    control_label: str | None = None
    fdr_alpha: float | None = None
    conditions: list[RgConditionSummary]
    pairwise_comparisons: list[RgRunPairwiseComparison]
    anova: None = None  # type: ignore[assignment]
    anova_by_run: list[RgRunANOVA] | None = None
    ranking_by_run: dict[str, list[str]]
    ranking: list[str] = Field(default_factory=list)
    equilibration_time: str

    def get_comparisons_for_run(self, run_label: str) -> list[RgRunPairwiseComparison]:
        """Get all pairwise comparisons for a specific run.

        Parameters
        ----------
        run_label : str
            The run label.

        Returns
        -------
        list[RgRunPairwiseComparison]
            Comparisons for this run.
        """
        return [c for c in self.pairwise_comparisons if c.run_label == run_label]

    def get_run_comparison(
        self, run_label: str, condition_b: str
    ) -> RgRunPairwiseComparison | None:
        """Get pairwise comparison for a specific run and condition vs control.

        Parameters
        ----------
        run_label : str
            The run label.
        condition_b : str
            Treatment condition label.

        Returns
        -------
        RgRunPairwiseComparison or None
            The comparison, or None if not found.
        """
        for comp in self.pairwise_comparisons:
            if comp.run_label == run_label and comp.condition_b == condition_b:
                return comp
        return None

    def get_ranking(self, run_label: str) -> list[str]:
        """Get ranking for a specific run.

        Parameters
        ----------
        run_label : str
            The run label.

        Returns
        -------
        list[str]
            Condition labels sorted by mean Rg (lowest first).
        """
        return self.ranking_by_run.get(run_label, [])
