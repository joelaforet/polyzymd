"""Models the framework hands to a plugin that are not observables.

The comparison result ladder this module used to hold (``MetricValue``,
``ConditionSummary``, ``PairwiseResult``, ``ANOVAResult``, ``ComparisonResult``
and the generic ``BaseComparisonResult``) served plugins that reported a
dictionary of scalars. Every plugin reports observables now, so the contract's
``ObservableAggregate`` and ``ObservableComparison`` carry those numbers and the
ladder has gone. What is left is the two settings models a plugin class still
declares.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel

__all__ = ["BasePlotSettings", "SlurmResourceHint"]


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
