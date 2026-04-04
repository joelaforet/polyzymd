"""Result models for SASA comparison analysis."""

from __future__ import annotations

from typing import ClassVar

from pydantic import BaseModel, Field

from polyzymd.compare.core.base import BaseComparisonResult, BaseConditionSummary


class SASARunSummary(BaseModel):
    """Summary statistics for one SASA run across replicates."""

    label: str
    target_selection: str
    context_selection: str
    mean_sasa: float
    sem_sasa: float
    per_replicate_means: list[float] = Field(default_factory=list)
    zero_atom_selection: bool = False


class SASAConditionSummary(BaseConditionSummary):
    """Summary statistics for one condition in a SASA comparison."""

    replicate_values: list[float] = Field(default_factory=list)
    run_summaries: list[SASARunSummary]

    @property
    def primary_metric_value(self) -> float:
        """Return 0.0 as multi-run SASA has no single primary metric."""
        return 0.0

    @property
    def primary_metric_sem(self) -> float:
        """Return 0.0 as multi-run SASA has no single primary metric."""
        return 0.0

    def get_run(self, label: str) -> SASARunSummary:
        """Get a run summary by label."""
        for run_summary in self.run_summaries:
            if run_summary.label == label:
                return run_summary
        raise KeyError(f"Run '{label}' not found in condition '{self.label}'")


class SASARunPairwiseComparison(BaseModel):
    """Statistical comparison between two conditions for one SASA run."""

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


class SASARunANOVA(BaseModel):
    """ANOVA result for a single SASA run."""

    run_label: str
    f_statistic: float
    p_value: float
    significant: bool


class SASAComparisonResult(BaseComparisonResult[SASAConditionSummary, SASARunPairwiseComparison]):
    """Complete SASA comparison analysis result."""

    comparison_type: ClassVar[str] = "sasa"

    metric: str = "mean_sasa"
    name: str
    n_runs: int
    run_labels: list[str]
    control_label: str | None = None
    conditions: list[SASAConditionSummary]
    pairwise_comparisons: list[SASARunPairwiseComparison]
    anova: None = None  # type: ignore[assignment]
    anova_by_run: list[SASARunANOVA] | None = None
    ranking_by_run: dict[str, list[str]]
    ranking: list[str] = Field(default_factory=list)
    equilibration_time: str

    def get_comparisons_for_run(self, run_label: str) -> list[SASARunPairwiseComparison]:
        """Get all pairwise comparisons for a specific run."""
        return [c for c in self.pairwise_comparisons if c.run_label == run_label]

    def get_ranking(self, run_label: str) -> list[str]:
        """Get ranking for a specific run."""
        return self.ranking_by_run.get(run_label, [])
