"""Comparison result models for analysis plugins."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import ClassVar, Generic, Self, TypeVar

from pydantic import BaseModel, ConfigDict, Field


class BasePlotSettings(BaseModel):
    """Base class for per-analysis plot settings."""


class SlurmResourceHint(BaseModel):
    """Per-plugin SLURM resource hints for HPC submission."""

    mem: str | None = None
    time: str | None = None
    cpus_per_task: int | None = None


@dataclass(frozen=True)
class MetricValue:
    """A scalar metric extracted from one aggregated condition result."""

    name: str
    mean: float
    sem: float
    replicate_values: list[float]
    higher_is_better: bool | None = True
    direction_labels: tuple[str, str, str] = ("decreased", "unchanged", "increased")


class ConditionSummary(BaseModel):
    """Summary statistics for one condition in a scalar comparison."""

    model_config = {"extra": "allow"}

    label: str
    n_replicates: int = 0


class PairwiseResult(BaseModel):
    """Statistical comparison between two conditions for one metric."""

    model_config = ConfigDict(ser_json_inf_nan="strings")

    condition_a: str
    condition_b: str
    metric: str = "default"
    t_statistic: float
    p_value: float
    p_value_adjusted: float | None = None
    posthoc_method: str = "ttest_bh"
    cohens_d: float
    effect_size_interpretation: str
    direction: str
    significant: bool
    percent_change: float
    testable: bool = True
    note: str | None = None


class ANOVAResult(BaseModel):
    """One-way ANOVA result for one metric."""

    metric: str = "default"
    f_statistic: float
    p_value: float
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
    def primary_metric_sem(self) -> float:
        """Return the SEM of the primary metric."""


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

    def get_comparison(self, label: str | tuple[str, str]) -> TPairwiseResult | None:
        """Get a pairwise comparison by condition pair.

        Parameters
        ----------
        label : str or tuple[str, str]
            Either a treatment label for legacy lookup or an explicit
            ``(condition_a, condition_b)`` pair.

        Returns
        -------
        TPairwiseResult or None
            The matching comparison, or ``None`` if not found.
        """
        if isinstance(label, tuple):
            condition_a, condition_b = label
            for comparison in self.pairwise_comparisons:
                if comparison.condition_a == condition_a and comparison.condition_b == condition_b:
                    return comparison
            legacy_label = condition_b
        else:
            legacy_label = label

        matches = [
            comparison
            for comparison in self.pairwise_comparisons
            if comparison.condition_b == legacy_label
        ]
        if not matches:
            return None
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous legacy comparison lookup for label '{legacy_label}': "
                f"found {len(matches)} matches; use tuple lookup (condition_a, condition_b)."
            )
        return matches[0]
