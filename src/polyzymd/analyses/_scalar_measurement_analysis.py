"""Analysis adapter for scalar measurement strategies."""

from __future__ import annotations

import hashlib
import inspect
import json
import math
import statistics
from abc import ABC
from typing import Any, ClassVar, Sequence

from pydantic import BaseModel, ValidationError

from polyzymd.analyses._measurement import CacheIdentity, MetricSpec, ScalarMeasurement
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
            if inspect.isabstract(measurement):
                raise TypeError(
                    f"{cls.__name__}.measurement must be a concrete ScalarMeasurement class."
                )
            cls._validate_metric_contract(measurement)
            return
        if not isinstance(measurement, ScalarMeasurement):
            raise TypeError(
                f"{cls.__name__}.measurement must be a ScalarMeasurement instance or class, "
                f"got {measurement!r}."
            )
        cls._validate_metric_contract(measurement)

    @classmethod
    def _validate_metric_contract(
        cls, measurement: type[ScalarMeasurement] | ScalarMeasurement
    ) -> None:
        """Raise a clear error for invalid scalar measurement metric metadata.

        Parameters
        ----------
        measurement : type[ScalarMeasurement] or ScalarMeasurement
            Measurement class or instance configured on the analysis class.

        Raises
        ------
        TypeError
            If ``measurement.metric`` is missing, is not a ``MetricSpec``, or
            has a blank metric name.
        """

        metric = getattr(measurement, "metric", None)
        if not isinstance(metric, MetricSpec) or not metric.name.strip():
            raise TypeError(
                f"{cls.__name__}.measurement.metric must be a MetricSpec with a non-empty name."
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

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Run one scalar replicate with canonical cache-hit validation.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            One-indexed replicate number.

        Returns
        -------
        Any
            Cached or freshly computed replicate result.
        """

        if not ctx.recompute and ctx.result_path is not None and ctx.result_path.exists():
            cached = self._load_scalar_replicate_cache(ctx, replicate)
            if cached is not None:
                return cached
        return super().run_replicate(ctx, replicate)

    def _load_scalar_replicate_cache(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> dict[str, Any] | None:
        """Load a compatible canonical scalar replicate cache.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            One-indexed replicate number expected in the cache payload.

        Returns
        -------
        dict[str, Any] or None
            Compatible cached replicate result, or ``None`` when the existing
            file is stale or incompatible.
        """

        if ctx.result_path is None:
            return None
        try:
            cached = json.loads(ctx.result_path.read_text())
        except (OSError, UnicodeDecodeError, json.JSONDecodeError):
            return None
        if not isinstance(cached, dict):
            return None

        measurement = self._measurement_instance()
        cache_identity = measurement.cache_identity(ctx.settings, analysis_name=self.name)
        try:
            self._validate_replicate_result(
                cached,
                expected_cache_identity=cache_identity,
                expected_measurement=measurement,
                expected_replicate=replicate,
                index=0,
            )
        except (TypeError, ValueError):
            return None
        return cached

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

        if not results:
            raise ValueError(f"{type(self).__name__}.aggregate() requires at least one result.")

        measurement = self._measurement_instance()
        cache_identity = measurement.cache_identity(ctx.settings, analysis_name=self.name)
        for index, result in enumerate(results):
            try:
                expected_replicate = ctx.replicates[index]
            except IndexError as exc:
                raise ValueError(
                    f"{type(self).__name__}.aggregate() received more results than replicate IDs."
                ) from exc
            self._validate_replicate_result(
                result,
                expected_cache_identity=cache_identity,
                expected_measurement=measurement,
                expected_replicate=expected_replicate,
                index=index,
            )

        values = [float(result["value"]) for result in results]
        mean_value = statistics.mean(values)
        sem_value = statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else 0.0
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

    def _validate_replicate_result(
        self,
        result: dict[str, Any],
        *,
        expected_cache_identity: CacheIdentity,
        expected_measurement: ScalarMeasurement,
        expected_replicate: int,
        index: int,
    ) -> None:
        """Validate that a scalar replicate result matches the active identity.

        Parameters
        ----------
        result : dict[str, Any]
            Per-replicate scalar result payload to validate.
        expected_cache_identity : CacheIdentity
            Cache identity for the active analysis settings and measurement.
        expected_measurement : ScalarMeasurement
            Active measurement configured on this analysis instance.
        expected_replicate : int
            One-indexed replicate ID expected in the result payload.
        index : int
            Zero-based result index used in error messages.

        Raises
        ------
        TypeError
            If the replicate result is not a dictionary.
        ValueError
            If the result identity does not match the active analysis,
            measurement, metric, or settings payload.
        """

        if not isinstance(result, dict):
            raise TypeError(
                f"{type(self).__name__}.aggregate() expected replicate result {index} to be a "
                f"dict, got {type(result).__name__}."
            )

        expected_fields = {
            "analysis": self.name,
            "measurement": expected_measurement.name,
            "metric": expected_measurement.metric.name,
        }
        for field, expected_value in expected_fields.items():
            actual_value = result.get(field)
            if actual_value != expected_value:
                raise ValueError(
                    f"Stale scalar replicate result {index}: {field} mismatch "
                    f"({actual_value!r} != {expected_value!r})."
                )

        actual_replicate = result.get("replicate")
        if actual_replicate != expected_replicate:
            raise ValueError(
                f"Stale scalar replicate result {index}: replicate mismatch "
                f"({actual_replicate!r} != {expected_replicate!r})."
            )

        raw_identity = result.get("cache_identity")
        if not isinstance(raw_identity, dict):
            raise ValueError(f"Stale scalar replicate result {index}: cache_identity is missing.")

        for field in ("analysis_name", "measurement_name", "version", "payload"):
            if field not in raw_identity:
                raise ValueError(
                    f"Stale scalar replicate result {index}: cache_identity.{field} is missing."
                )

        try:
            actual_identity = CacheIdentity.model_validate(raw_identity)
        except ValidationError as exc:
            raise ValueError(
                f"Stale scalar replicate result {index}: cache_identity is invalid."
            ) from exc

        expected_identity_fields = {
            "analysis_name": expected_cache_identity.analysis_name,
            "measurement_name": expected_cache_identity.measurement_name,
            "version": expected_cache_identity.version,
        }
        for field, expected_value in expected_identity_fields.items():
            actual_value = getattr(actual_identity, field)
            if actual_value != expected_value:
                raise ValueError(
                    f"Stale scalar replicate result {index}: cache_identity.{field} mismatch "
                    f"({actual_value!r} != {expected_value!r})."
                )

        if actual_identity.payload != expected_cache_identity.payload:
            raise ValueError(
                f"Stale scalar replicate result {index}: cache_identity settings payload mismatch."
            )

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

        measurement = self._measurement_instance()
        metric_value = measurement.metric_value_from_summary(summary)
        return {measurement.metric.name: metric_value}
