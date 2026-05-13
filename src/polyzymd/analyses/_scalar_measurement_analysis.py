"""Analysis adapter for scalar measurement strategies."""

from __future__ import annotations

import math
import statistics
from abc import ABC, abstractmethod
from typing import Any, ClassVar, Sequence

from pydantic import BaseModel

from polyzymd.analyses._measurement import ScalarMeasurement
from polyzymd.analyses.base import AggregateContext, Analysis, MetricValue, ReplicateContext


class _EmptyScalarMeasurementSettings(BaseModel):
    """Default empty settings for scalar measurement analyses."""


class _ScalarMeasurementRunner:
    """Runner bridge that adapts a scalar measurement to the runner seam."""

    def __init__(
        self,
        *,
        measurement: ScalarMeasurement,
        universe: Any,
        settings: BaseModel,
    ) -> None:
        """Initialize the runner bridge.

        Parameters
        ----------
        measurement : ScalarMeasurement
            Scalar measurement strategy to execute.
        universe : Any
            Trajectory universe supplied by the analysis framework.
        settings : BaseModel
            Resolved analysis settings.
        """

        self.measurement = measurement
        self.universe = universe
        self.settings = settings
        self.results: dict[str, float] = {}

    def run(
        self,
        *,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> _ScalarMeasurementRunner:
        """Execute the scalar measurement over the requested frame window.

        Parameters
        ----------
        start : int or None, optional
            First frame index included in the trajectory window.
        stop : int or None, optional
            Exclusive final frame index for the trajectory window.
        step : int or None, optional
            Frame stride for the trajectory window.

        Returns
        -------
        _ScalarMeasurementRunner
            This runner with populated ``results``.
        """

        value = self.measurement.measure(
            self.universe,
            self.settings,
            start=start,
            stop=stop,
            step=step,
        )
        self.results = {"value": float(value)}
        return self


class ScalarMeasurementAnalysis(Analysis, ABC):
    """Adapter base class for analyses backed by one scalar measurement."""

    Settings: ClassVar[type[BaseModel]] = _EmptyScalarMeasurementSettings

    @property
    @abstractmethod
    def measurement(self) -> ScalarMeasurement:
        """Return the scalar measurement strategy for this analysis."""

    def _measurement_instance(self) -> ScalarMeasurement:
        """Return the configured measurement strategy instance.

        Returns
        -------
        ScalarMeasurement
            Scalar measurement strategy configured on the analysis class.
        """

        measurement = type(self).measurement
        if isinstance(measurement, type):
            measurement = measurement()
        if not isinstance(measurement, ScalarMeasurement):
            raise TypeError(
                f"{type(self).__name__}.measurement must be a ScalarMeasurement instance or class, "
                f"got {measurement!r}."
            )
        return measurement

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> _ScalarMeasurementRunner:
        """Build a runner bridge for one scalar measurement replicate.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            One-indexed replicate number.
        universe : Any
            Loaded trajectory universe.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        _ScalarMeasurementRunner
            Runner bridge compatible with the analysis framework.
        """

        del replicate, window
        return _ScalarMeasurementRunner(
            measurement=self._measurement_instance(),
            universe=universe,
            settings=ctx.settings,
        )

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: _ScalarMeasurementRunner,
        window: Any,
    ) -> dict[str, Any]:
        """Convert runner output into a plain replicate result.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            One-indexed replicate number.
        runner : _ScalarMeasurementRunner
            Executed scalar measurement runner.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        dict[str, Any]
            Serializable per-replicate result payload.
        """

        del window
        measurement = self._measurement_instance()
        return {
            "analysis": self.name,
            "measurement": measurement.name,
            "metric": measurement.metric.name,
            "replicate": replicate,
            "value": float(runner.results["value"]),
            "cache_identity": measurement.cache_identity(
                ctx.settings,
                analysis_name=self.name,
            ).model_dump(mode="json"),
        }

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[dict[str, Any]],
    ) -> dict[str, Any]:
        """Aggregate scalar replicate results into summary statistics.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[dict[str, Any]]
            Per-replicate scalar result payloads.

        Returns
        -------
        dict[str, Any]
            Serializable aggregate result payload.
        """

        values = [float(result["value"]) for result in results]
        mean_value = statistics.mean(values) if values else math.nan
        sem_value = statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else 0.0
        measurement = self._measurement_instance()
        cache_identity = measurement.cache_identity(ctx.settings, analysis_name=self.name)
        return {
            "analysis": self.name,
            "measurement": measurement.name,
            "metric": measurement.metric.name,
            "mean_value": mean_value,
            "sem_value": sem_value,
            "replicate_values": values,
            "n_replicates": len(values),
            "replicates": list(ctx.replicates),
            "settings_fingerprint": self.aggregate_settings_fingerprint(ctx.settings),
            "cache_identity": cache_identity.model_dump(mode="json"),
        }

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract the scalar metric for default cross-condition comparison.

        Parameters
        ----------
        summary : Any
            Aggregated scalar measurement result.

        Returns
        -------
        dict[str, MetricValue]
            Mapping from metric name to comparison metric value.
        """

        metric = self._measurement_instance().metric
        if not isinstance(summary, dict):
            raise TypeError(
                f"{type(self).__name__}.extract_metrics() expects a dict aggregate result, "
                f"got {type(summary).__name__}."
            )
        metric_value = metric.to_metric_value(
            mean=float(summary["mean_value"]),
            sem=float(summary["sem_value"]),
            replicate_values=[float(value) for value in summary["replicate_values"]],
        )
        return {metric.name: metric_value}
