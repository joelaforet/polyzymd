"""Measurement metadata for catalytic triad scalar comparison."""

from __future__ import annotations

from typing import Any, ClassVar

from pydantic import BaseModel

from polyzymd.analyses.base import MetricSpec, ScalarMeasurement


class TriadSimultaneousContactMeasurement(ScalarMeasurement):
    """Scalar measurement metadata for simultaneous catalytic-triad contact."""

    name: ClassVar[str] = "simultaneous_contact_fraction"
    metric: ClassVar[MetricSpec] = MetricSpec(
        name="simultaneous_contact_fraction",
        higher_is_better=True,
        label="Simultaneous Contact",
        unit="%",
        direction_labels=("worsening", "unchanged", "improving"),
    )

    def measure(
        self,
        universe: Any,
        settings: BaseModel,
        *,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> float:
        """Reject direct measurement because triad uses a custom pair-distance runner.

        Parameters
        ----------
        universe : Any
            Ignored trajectory universe.
        settings : BaseModel
            Ignored triad settings.
        start : int or None, optional
            Ignored start frame.
        stop : int or None, optional
            Ignored stop frame.
        step : int or None, optional
            Ignored frame stride.

        Raises
        ------
        NotImplementedError
            Always raised because the catalytic-triad plugin keeps its legacy
            pair-distance runner and aggregation lifecycle.
        """

        del universe, settings, start, stop, step
        raise NotImplementedError(
            "Catalytic triad simultaneous contact is computed by the plugin runner."
        )

    def metric_value_from_summary(self, summary: Any) -> Any:
        """Extract the legacy percent-scaled comparison metric from a triad summary.

        Parameters
        ----------
        summary : Any
            Aggregated catalytic-triad result with simultaneous contact fields.

        Returns
        -------
        Any
            ``MetricValue`` preserving the legacy metric key, percent scaling,
            ranking direction, and direction labels.
        """

        return self.metric.to_metric_value(
            mean=float(summary.overall_simultaneous_contact) * 100,
            sem=float(summary.sem_simultaneous_contact) * 100,
            replicate_values=[float(value) * 100 for value in summary.per_replicate_simultaneous],
        )
