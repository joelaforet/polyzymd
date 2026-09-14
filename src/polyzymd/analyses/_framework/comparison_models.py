"""Comparison result models for analysis plugins."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import ClassVar, Generic, Literal, Self, TypeVar

from pydantic import BaseModel, ConfigDict, Field


class BasePlotSettings(BaseModel):
    """Base class for per-analysis plot settings.

    ``error_bar`` chooses the interval drawn on comparison bars and bands.
    ``"ci95"`` draws the 95 percent Student t interval across replicates, which
    is what Grossfield et al. (2018) ask authors to graph. ``"sem"`` draws one
    standard error, which at n = 3 is 4.3 times narrower.
    """

    error_bar: Literal["ci95", "sem"] = "ci95"


class SlurmResourceHint(BaseModel):
    """Per-plugin SLURM resource hints for HPC submission."""

    mem: str | None = None
    time: str | None = None
    cpus_per_task: int | None = None


@dataclass(frozen=True)
class MetricValue:
    """A scalar metric extracted from one aggregated condition result.

    The replicate is the sampling unit. ``unit`` names the physical unit of
    ``mean`` (``"A"``, ``"A^2"``, ``"%"``, ``"ns"``) and is ``None`` only for a
    dimensionless metric. The standard error, both limits and ``ci_method`` are
    ``None`` when a single replicate makes them inestimable.
    """

    name: str
    mean: float
    sem: float | None
    replicate_values: list[float]
    unit: str | None = None
    ci95_low: float | None = None
    ci95_high: float | None = None
    ci_method: str | None = None
    higher_is_better: bool | None = True
    direction_labels: tuple[str, str, str] = ("decreased", "unchanged", "increased")

    @classmethod
    def from_replicate_values(
        cls,
        name: str,
        replicate_values: Sequence[float],
        *,
        unit: str | None = None,
        higher_is_better: bool | None = True,
        direction_labels: tuple[str, str, str] = ("decreased", "unchanged", "increased"),
        scale: float = 1.0,
    ) -> "MetricValue":
        """Build a metric from replicate values, deriving SEM and interval.

        ``scale`` multiplies every value before the statistics are computed, for
        a unit change such as fraction to percent.
        """
        from polyzymd.analyses.shared.statistics import mean_sem_ci

        scaled = [float(value) * float(scale) for value in replicate_values]
        stats = mean_sem_ci(scaled)
        return cls(
            name=name,
            mean=stats.mean,
            sem=stats.sem,
            replicate_values=scaled,
            unit=unit,
            ci95_low=stats.ci_low,
            ci95_high=stats.ci_high,
            ci_method=stats.ci_method,
            higher_is_better=higher_is_better,
            direction_labels=direction_labels,
        )


class ConditionSummary(BaseModel):
    """Summary statistics for one condition in a scalar comparison."""

    model_config = {"extra": "allow"}

    label: str
    n_replicates: int = 0


class PairwiseResult(BaseModel):
    """Statistical comparison between two conditions for one metric.

    ``p_value`` is the raw two-tailed p-value and ``p_value_adjusted`` its
    Benjamini-Hochberg value within the run's pairwise family.
    ``hedges_g`` is Cohen's d after the Hedges (1981) small-sample
    correction. ``effect_size_interpretation`` is ``None`` when the
    combined sample is too small for a Cohen adjective to mean anything.
    ``direction`` reads "no significant change" unless ``significant`` is
    true.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    condition_a: str
    condition_b: str
    metric: str = "default"
    t_statistic: float
    p_value: float
    p_value_adjusted: float | None
    posthoc_method: str = "ttest_bh"
    cohens_d: float
    hedges_g: float | None = None
    effect_size_interpretation: str | None = None
    direction: str
    significant: bool
    percent_change: float
    testable: bool = True
    note: str | None = None


class ANOVAResult(BaseModel):
    """One-way ANOVA result for one metric.

    The ANOVA is an omnibus test outside the pairwise Benjamini-Hochberg
    family, so ``p_value_adjusted`` is always ``None`` and ``significant``
    compares the raw ``p_value`` with alpha.
    """

    metric: str = "default"
    f_statistic: float
    p_value: float
    p_value_adjusted: float | None = None
    significant: bool
    testable: bool = True
    note: str | None = None


class ComparisonResult(BaseModel):
    """Serializable result of a default scalar cross-condition comparison."""

    model_config = ConfigDict(ser_json_inf_nan="strings")

    analysis_type: str
    name: str
    control_label: str | None = None
    fdr_alpha: float | None = None
    ttest_method: str = "student"
    posthoc_method: str = "ttest_bh"
    uncertainty: dict[str, object] | None = None
    conditions: list[ConditionSummary] = Field(default_factory=list)
    pairwise_comparisons: list[PairwiseResult] = Field(default_factory=list)
    anova: list[ANOVAResult] | None = None
    ranking: list[str] = Field(default_factory=list)
    rankings_by_metric: dict[str, list[str]] | None = None
    equilibration_time: str = "0ns"
    created_at: str = ""
    polyzymd_version: str = ""

    def save(self, path: Path | str) -> Path:
        """Save the result to a JSON file.

        Parameters
        ----------
        path : Path or str
            Output path.

        Returns
        -------
        Path
            Path to the saved file.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path | str) -> Self:
        """Load the result from a JSON file.

        Parameters
        ----------
        path : Path or str
            Path to a JSON file.

        Returns
        -------
        Self
            Loaded result.
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())


class BaseConditionSummary(BaseModel, ABC):
    """Abstract base class for condition-level custom comparison summaries."""

    label: str
    config_path: str
    n_replicates: int
    replicate_values: list[float]

    @property
    @abstractmethod
    def primary_metric_value(self) -> float:
        """Return the primary metric value for ranking and comparison."""

    @property
    @abstractmethod
    def primary_metric_sem(self) -> float | None:
        """Return the SEM of the primary metric, or ``None`` for one replicate."""


TConditionSummary = TypeVar("TConditionSummary", bound=BaseConditionSummary)
TPairwiseResult = TypeVar("TPairwiseResult", bound=PairwiseResult)


class BaseComparisonResult(BaseModel, ABC, Generic[TConditionSummary, TPairwiseResult]):
    """Abstract base class for custom plugin comparison results."""

    model_config = ConfigDict(ser_json_inf_nan="strings")
    comparison_type: ClassVar[str] = "base"

    metric: str
    name: str
    control_label: str | None = None
    conditions: list[TConditionSummary]
    pairwise_comparisons: list[TPairwiseResult]
    anova: ANOVAResult | list[ANOVAResult] | None = None
    ranking: list[str]
    equilibration_time: str
    created_at: datetime
    polyzymd_version: str

    def save(self, path: Path | str) -> Path:
        """Save the result to a JSON file.

        Parameters
        ----------
        path : Path or str
            Output path.

        Returns
        -------
        Path
            Path to the saved file.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path | str) -> Self:
        """Load the result from a JSON file.

        Parameters
        ----------
        path : Path or str
            Path to a JSON file.

        Returns
        -------
        Self
            Loaded result.
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())

    def get_condition(self, label: str) -> TConditionSummary:
        """Get a condition by label.

        Parameters
        ----------
        label : str
            Condition label.

        Returns
        -------
        TConditionSummary
            The matching condition summary.

        Raises
        ------
        KeyError
            If the condition is not found.
        """
        for condition in self.conditions:
            if condition.label == label:
                return condition
        raise KeyError(f"Condition '{label}' not found")

    def get_comparison(self, label: tuple[str, str]) -> TPairwiseResult | None:
        """Get a pairwise comparison by condition pair.

        Parameters
        ----------
        label : tuple[str, str]
            Explicit ``(condition_a, condition_b)`` pair.

        Returns
        -------
        TPairwiseResult or None
            The matching comparison, or ``None`` if not found.
        """
        if not isinstance(label, tuple):
            raise TypeError("Comparison lookup requires a (condition_a, condition_b) tuple")
        condition_a, condition_b = label
        for comparison in self.pairwise_comparisons:
            if comparison.condition_a == condition_a and comparison.condition_b == condition_b:
                return comparison
        return None
