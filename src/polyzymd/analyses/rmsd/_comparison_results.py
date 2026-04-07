"""Result models for RMSD comparison analysis.

This module defines Pydantic models for structured RMSD comparison results
that can be serialized to JSON and used for downstream plotting.

Each RMSD run is compared independently across conditions — there is no
cross-run averaging since different runs measure fundamentally different
selections (e.g., protein backbone vs polymer core).

For each run:
- Primary ranking: mean RMSD (lower = more stable = better)

.. versionadded:: 1.3.0
"""

from __future__ import annotations

from typing import ClassVar

from pydantic import BaseModel, ConfigDict, Field

from polyzymd.compare.core.base import BaseComparisonResult, BaseConditionSummary


class RMSDRunSummary(BaseModel):
    """Summary statistics for one RMSD run across replicates.

    Attributes
    ----------
    label : str
        Human-readable label for this run (e.g., "protein_backbone").
    selection : str
        MDAnalysis selection used for RMSD calculation.
    mean_rmsd : float
        Mean RMSD in Angstroms (averaged across replicates).
    sem_rmsd : float
        SEM of RMSD across replicates.
    per_replicate_means : list[float]
        Mean RMSD from each replicate (for statistical tests).
    """

    label: str
    selection: str
    mean_rmsd: float
    sem_rmsd: float
    per_replicate_means: list[float] = Field(default_factory=list)


class RMSDConditionSummary(BaseConditionSummary):
    """Summary statistics for one condition in an RMSD comparison.

    Note: There is no "overall" RMSD metric across runs, since averaging
    RMSD from different selections is not meaningful. Each run is compared
    independently.

    Inherits from BaseConditionSummary for interface consistency. The
    ``replicate_values`` field defaults to an empty list because RMSD
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
        Empty by default (no single primary metric for multi-run RMSD).
    run_summaries : list[RMSDRunSummary]
        Summary for each RMSD run.
    """

    replicate_values: list[float] = Field(default_factory=list)
    run_summaries: list[RMSDRunSummary]

    @property
    def primary_metric_value(self) -> float:
        """Return 0.0 — multi-run RMSD has no single primary metric."""
        return 0.0

    @property
    def primary_metric_sem(self) -> float:
        """Return 0.0 — multi-run RMSD has no single primary metric."""
        return 0.0

    def get_run(self, label: str) -> RMSDRunSummary:
        """Get a run summary by label.

        Parameters
        ----------
        label : str
            Run label.

        Returns
        -------
        RMSDRunSummary
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


class RMSDRunPairwiseComparison(BaseModel):
    """Statistical comparison between two conditions for a single RMSD run.

    Attributes
    ----------
    run_label : str
        Label of the RMSD run being compared.
    condition_a : str
        Label of first condition (typically control).
    condition_b : str
        Label of second condition (typically treatment).
    t_statistic : float
        T-test statistic for mean RMSD.
    p_value : float
        Two-tailed p-value.
    cohens_d : float
        Effect size (Cohen's d).
    effect_interpretation : str
        "negligible", "small", "medium", or "large".
    direction : str
        "stabilizing" (lower RMSD), "destabilizing" (higher RMSD),
        or "unchanged".
    significant : bool
        Whether p < 0.05.
    percent_change : float
        Percent change in mean RMSD.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    run_label: str
    condition_a: str
    condition_b: str
    t_statistic: float
    p_value: float
    cohens_d: float
    effect_interpretation: str
    direction: str
    significant: bool
    percent_change: float


class RMSDRunANOVA(BaseModel):
    """ANOVA result for a single RMSD run.

    Attributes
    ----------
    run_label : str
        Label of the RMSD run.
    f_statistic : float
        F-statistic for mean RMSD.
    p_value : float
        P-value for ANOVA.
    significant : bool
        Whether p < 0.05.
    """

    run_label: str
    f_statistic: float
    p_value: float
    significant: bool


class RMSDComparisonResult(BaseComparisonResult[RMSDConditionSummary, RMSDRunPairwiseComparison]):
    """Complete RMSD comparison analysis result.

    Inherits save/load/get_condition from BaseComparisonResult.

    Each RMSD run is compared independently — rankings and statistics
    are computed per-run since averaging RMSD from different selections
    is not meaningful.

    Attributes
    ----------
    metric : str
        Primary metric: "mean_rmsd".
    name : str
        Name of the comparison project.
    n_runs : int
        Number of RMSD runs analyzed.
    run_labels : list[str]
        Labels for each run.
    control_label : str | None
        Label of the control condition.
    conditions : list[RMSDConditionSummary]
        Summary for each condition.
    pairwise_comparisons : list[RMSDRunPairwiseComparison]
        Statistical comparisons (grouped by run).
    anova_by_run : list[RMSDRunANOVA] | None
        ANOVA results for each run (if 3+ conditions).
    ranking_by_run : dict[str, list[str]]
        Condition labels sorted by mean RMSD, keyed by run label.
        Lower RMSD = better = earlier in list.
    ranking : list[str]
        Empty by default (RMSD ranks per-run via ranking_by_run).
    equilibration_time : str
        Equilibration time used.
    """

    comparison_type: ClassVar[str] = "rmsd"
    model_config = ConfigDict(ser_json_inf_nan="strings")

    metric: str = "mean_rmsd"
    name: str
    n_runs: int
    run_labels: list[str]
    control_label: str | None = None
    conditions: list[RMSDConditionSummary]
    pairwise_comparisons: list[RMSDRunPairwiseComparison]
    # Override base anova field — RMSD uses per-run ANOVA instead
    anova: None = None  # type: ignore[assignment]
    anova_by_run: list[RMSDRunANOVA] | None = None
    ranking_by_run: dict[str, list[str]]
    ranking: list[str] = Field(default_factory=list)
    equilibration_time: str

    def get_comparisons_for_run(self, run_label: str) -> list[RMSDRunPairwiseComparison]:
        """Get all pairwise comparisons for a specific run.

        Parameters
        ----------
        run_label : str
            The run label.

        Returns
        -------
        list[RMSDRunPairwiseComparison]
            Comparisons for this run.
        """
        return [c for c in self.pairwise_comparisons if c.run_label == run_label]

    def get_run_comparison(
        self, run_label: str, condition_b: str
    ) -> RMSDRunPairwiseComparison | None:
        """Get pairwise comparison for a specific run and condition vs control.

        Parameters
        ----------
        run_label : str
            The run label.
        condition_b : str
            Treatment condition label.

        Returns
        -------
        RMSDRunPairwiseComparison or None
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
            Condition labels sorted by mean RMSD (lowest first).
        """
        return self.ranking_by_run.get(run_label, [])
