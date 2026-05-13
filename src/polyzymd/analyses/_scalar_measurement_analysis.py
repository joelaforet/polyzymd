"""Analysis adapter for scalar measurement strategies."""

from __future__ import annotations

import hashlib
import json
import math
import statistics
from abc import ABC
from typing import Any, ClassVar, Sequence

from pydantic import BaseModel

from polyzymd.analyses._measurement import ScalarMeasurement
from polyzymd.analyses.base import AggregateContext, Analysis, MetricValue, ReplicateContext


class _EmptyScalarMeasurementSettings(BaseModel):
    """Default empty settings for scalar measurement analyses."""


class _MeasurementClassVarContract:
    """Abstract class-variable marker for scalar measurement adapters."""

    __isabstractmethod__ = True

    def __get__(self, instance: Any, owner: type | None = None) -> Any:
        """Raise a clear error when a subclass omits ``measurement``.

        Parameters
        ----------
        instance : Any
            Ignored descriptor instance target.
        owner : type or None, optional
            Class that attempted to access the missing measurement.

        Raises
        ------
        AttributeError
            Always raised because concrete subclasses must replace this marker.
        """

        if instance is None:
            return self
        owner_name = owner.__name__ if owner is not None else "ScalarMeasurementAnalysis"
        raise AttributeError(
            f"{owner_name}.measurement must be defined as a ScalarMeasurement instance or class."
        )


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
    """Adapter base class for analyses backed by one scalar measurement.

    Concrete subclasses declare ``measurement`` as a class variable containing
    either a :class:`ScalarMeasurement` instance or a concrete
    :class:`ScalarMeasurement` subclass. This keeps simple contributor plugins
    declarative while allowing the adapter to enforce cache identity centrally.
    """

    Settings: ClassVar[type[BaseModel]] = _EmptyScalarMeasurementSettings
    measurement: ClassVar[type[ScalarMeasurement] | ScalarMeasurement] = (  # type: ignore[assignment]
        _MeasurementClassVarContract()
    )

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """Validate the scalar measurement class-variable contract."""

        super().__init_subclass__(**kwargs)
        if cls is ScalarMeasurementAnalysis:
            return
        cls._validate_measurement_contract()

    @classmethod
    def _measurement_config(cls) -> type[ScalarMeasurement] | ScalarMeasurement | None:
        """Return the nearest configured measurement class variable.

        Returns
        -------
        type[ScalarMeasurement] or ScalarMeasurement or None
            Configured measurement strategy, or ``None`` when no concrete
            subclass replaced the abstract class-variable marker.
        """

        for owner in cls.__mro__:
            value = owner.__dict__.get("measurement")
            if value is None:
                continue
            if isinstance(value, _MeasurementClassVarContract):
                return None
            return value
        return None

    @classmethod
    def _validate_measurement_contract(cls) -> None:
        """Raise a clear error for missing or invalid measurement definitions.

        Raises
        ------
        TypeError
            If ``measurement`` is missing or is not a scalar measurement
            instance or class.
        """

        measurement = cls._measurement_config()
        if measurement is None:
            raise TypeError(
                f"{cls.__name__}.measurement must be defined as a ScalarMeasurement "
                "instance or class."
            )
        if isinstance(measurement, type):
            if not issubclass(measurement, ScalarMeasurement):
                raise TypeError(
                    f"{cls.__name__}.measurement must be a ScalarMeasurement instance or "
                    f"class, got {measurement!r}."
                )
            return
        if not isinstance(measurement, ScalarMeasurement):
            raise TypeError(
                f"{cls.__name__}.measurement must be a ScalarMeasurement instance or class, "
                f"got {measurement!r}."
            )

    def _measurement_instance(self) -> ScalarMeasurement:
        """Return the configured measurement strategy instance.

        Returns
        -------
        ScalarMeasurement
            Scalar measurement strategy configured on the analysis class.
        """

        measurement = type(self)._measurement_config()
        if measurement is None:
            raise TypeError(
                f"{type(self).__name__}.measurement must be defined as a ScalarMeasurement "
                "instance or class."
            )
        if isinstance(measurement, type):
            try:
                measurement = measurement()
            except TypeError as exc:
                raise TypeError(
                    f"{type(self).__name__}.measurement must be a concrete "
                    "ScalarMeasurement class."
                ) from exc
        if not isinstance(measurement, ScalarMeasurement):
            raise TypeError(
                f"{type(self).__name__}.measurement must be a ScalarMeasurement instance or class, "
                f"got {measurement!r}."
            )
        return measurement

    def aggregate_settings_fingerprint(self, settings: BaseModel | None) -> str:
        """Return a scalar cache fingerprint that includes measurement identity.

        Parameters
        ----------
        settings : BaseModel or None
            Active analysis settings.

        Returns
        -------
        str
            Short deterministic fingerprint including settings, measurement
            name, measurement version, and metric identity.
        """

        measurement = self._measurement_instance()
        cache_identity = measurement.cache_identity(settings, analysis_name=self.name)
        payload = {
            "cache_identity": cache_identity.model_dump(mode="json"),
            "metric": {"name": measurement.metric.name},
        }
        canonical = json.dumps(payload, sort_keys=True, default=str)
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:8]

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
