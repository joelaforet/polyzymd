"""Tests for the lightweight measurement analysis API."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, ClassVar

import pytest
from pydantic import BaseModel

from polyzymd.analyses import ScalarMeasurement as TopLevelScalarMeasurement
from polyzymd.analyses.base import (
    AggregateContext,
    CacheIdentity,
    Condition,
    Measurement,
    MetricSpec,
    MetricValue,
    ReplicateContext,
    ScalarMeasurement,
    ScalarMeasurementAnalysis,
)


class ToyScalarSettings(BaseModel):
    """Settings for toy scalar measurement tests."""

    offset: float = 0.0


class ToyScalarMeasurement(ScalarMeasurement):
    """Toy scalar measurement that records the requested frame window."""

    name: ClassVar[str] = "toy_measurement"
    version: ClassVar[str] = "2"
    metric: ClassVar[MetricSpec] = MetricSpec(
        name="toy_value",
        higher_is_better=False,
        label="Toy value",
        unit="arb",
        direction_labels=("lower", "unchanged", "higher"),
    )

    def measure(
        self,
        universe,
        settings: BaseModel,
        *,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> float:
        """Return a deterministic value from settings and the frame window."""

        del universe
        assert isinstance(settings, ToyScalarSettings)
        return float(settings.offset + (start or 0) + (stop or 0) + (step or 0))


class ToyScalarMeasurementV3(ToyScalarMeasurement):
    """Toy scalar measurement with a changed cache version."""

    version: ClassVar[str] = "3"


class ToyScalarMeasurementRenamedMetric(ToyScalarMeasurement):
    """Toy scalar measurement with a changed scalar metric key."""

    metric: ClassVar[MetricSpec] = MetricSpec(name="renamed_toy_value")


class ToyScalarAnalysis(ScalarMeasurementAnalysis):
    """Toy scalar measurement analysis with explicit settings."""

    name: ClassVar[str] = "toy_scalar"
    Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
    measurement: ClassVar[ScalarMeasurement] = ToyScalarMeasurement()


class ClassVariableScalarAnalysis(ScalarMeasurementAnalysis):
    """Toy scalar analysis configured with a measurement class."""

    name: ClassVar[str] = "class_variable_scalar"
    Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
    measurement: ClassVar[type[ScalarMeasurement]] = ToyScalarMeasurement


class DefaultSettingsScalarAnalysis(ScalarMeasurementAnalysis):
    """Toy scalar analysis that relies on adapter default settings."""

    name: ClassVar[str] = "default_settings_scalar"
    measurement: ClassVar[ScalarMeasurement] = ToyScalarMeasurement()


class VersionedScalarAnalysis(ScalarMeasurementAnalysis):
    """Toy scalar analysis with a changed measurement version."""

    name: ClassVar[str] = "toy_scalar"
    Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
    measurement: ClassVar[type[ScalarMeasurement]] = ToyScalarMeasurementV3


class RenamedMetricScalarAnalysis(ScalarMeasurementAnalysis):
    """Toy scalar analysis with a changed metric name."""

    name: ClassVar[str] = "toy_scalar"
    Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
    measurement: ClassVar[type[ScalarMeasurement]] = ToyScalarMeasurementRenamedMetric


class FakeTrajectory:
    """Minimal trajectory with a frame count."""

    def __len__(self) -> int:
        """Return the fake frame count."""

        return 10


class FakeUniverse:
    """Minimal universe exposing a trajectory."""

    trajectory = FakeTrajectory()


class FakeWindow:
    """Trajectory window with runner keyword arguments."""

    warning_message = None

    def run_kwargs(self) -> dict[str, int]:
        """Return fake runner keyword arguments."""

        return {"start": 1, "stop": 4, "step": 2}


class FakeLoader:
    """Loader seam used by the scalar measurement adapter."""

    def __init__(self, sim_config: object) -> None:
        """Store the provided simulation configuration."""

        self.sim_config = sim_config

    def load_universe(self, replicate: int) -> FakeUniverse:
        """Return a fake universe for the requested replicate."""

        assert replicate == 1
        return FakeUniverse()


def build_aggregate_context(tmp_path: Path, settings: BaseModel) -> AggregateContext:
    """Build an aggregate context for scalar measurement tests."""

    condition = Condition(
        label="Toy",
        config_path=tmp_path / "config.yaml",
        replicates=(1,),
        sim_config=object(),
    )
    return AggregateContext(
        condition=condition,
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        equilibration="0ns",
        settings=settings,
    )


def build_scalar_replicate_result(
    analysis: ScalarMeasurementAnalysis,
    settings: BaseModel,
    *,
    replicate: int = 1,
    value: float = 17.0,
) -> dict[str, object]:
    """Build a serialized scalar replicate result for aggregation tests."""

    measurement = analysis._measurement_instance()
    return {
        "analysis": analysis.name,
        "measurement": measurement.name,
        "metric": measurement.metric.name,
        "replicate": replicate,
        "value": value,
        "cache_identity": measurement.cache_identity(
            settings,
            analysis_name=analysis.name,
        ).model_dump(mode="json"),
    }


def test_metric_spec_converts_to_metric_value() -> None:
    """Metric metadata should convert to the existing comparison model."""

    spec = ToyScalarMeasurement.metric
    metric_value = spec.to_metric_value(mean=2.0, sem=0.5, replicate_values=[1.0, 3.0])

    assert isinstance(metric_value, MetricValue)
    assert metric_value.name == "toy_value"
    assert metric_value.higher_is_better is False
    assert metric_value.direction_labels == ("lower", "unchanged", "higher")


def test_scalar_measurement_metric_value_from_dict_summary() -> None:
    """Scalar measurements should convert default dict aggregates to metrics."""

    metric_value = ToyScalarMeasurement().metric_value_from_summary(
        {
            "mean_value": 2.0,
            "sem_value": 0.5,
            "replicate_values": [1, "2.5", 3.0],
        }
    )

    assert isinstance(metric_value, MetricValue)
    assert metric_value.name == "toy_value"
    assert metric_value.mean == 2.0
    assert metric_value.sem == 0.5
    assert metric_value.replicate_values == [1.0, 2.5, 3.0]
    assert metric_value.higher_is_better is False


def test_scalar_measurement_metric_value_rejects_non_dict_summary() -> None:
    """Default scalar metric extraction should require dict aggregates."""

    with pytest.raises(TypeError, match="expects a dict aggregate result"):
        ToyScalarMeasurement().metric_value_from_summary(object())


def test_cache_identity_has_independent_payloads() -> None:
    """Cache identity payloads should not share mutable defaults."""

    left = CacheIdentity(measurement_name="left")
    right = CacheIdentity(measurement_name="right")

    left.payload["changed"] = True

    assert right.payload == {}


def test_measurement_cache_identity_includes_settings_payload() -> None:
    """Measurement cache identity should include settings metadata."""

    identity = ToyScalarMeasurement().cache_identity(
        ToyScalarSettings(offset=3.5),
        analysis_name="toy_scalar",
    )

    assert identity == CacheIdentity(
        analysis_name="toy_scalar",
        measurement_name="toy_measurement",
        version="2",
        payload={"settings": {"offset": 3.5}},
    )


def test_measurement_api_reexported_from_facades() -> None:
    """Measurement classes should be available from public facades."""

    assert issubclass(ScalarMeasurement, Measurement)
    assert TopLevelScalarMeasurement is ScalarMeasurement
    assert ScalarMeasurementAnalysis.__module__ == "polyzymd.analyses.base"
    assert MetricSpec.__module__ == "polyzymd.analyses.base"


def test_scalar_measurement_analysis_default_settings_model() -> None:
    """Scalar measurement analyses can rely on the adapter default settings."""

    settings = DefaultSettingsScalarAnalysis.Settings()

    assert isinstance(settings, BaseModel)
    assert settings.model_dump() == {}


def test_scalar_measurement_analysis_accepts_measurement_class() -> None:
    """Scalar measurement analyses should accept a class-variable measurement class."""

    analysis = ClassVariableScalarAnalysis()
    measurement = analysis._measurement_instance()

    assert isinstance(measurement, ToyScalarMeasurement)
    assert measurement.name == "toy_measurement"


def test_scalar_measurement_analysis_requires_measurement() -> None:
    """Scalar measurement analyses should fail clearly when measurement is missing."""

    with pytest.raises(TypeError, match="MissingMeasurementAnalysis.measurement must be defined"):

        class MissingMeasurementAnalysis(ScalarMeasurementAnalysis):
            name: ClassVar[str] = "missing_measurement"
            Settings: ClassVar[type[BaseModel]] = ToyScalarSettings


def test_scalar_measurement_analysis_rejects_invalid_measurement() -> None:
    """Scalar measurement analyses should fail clearly for invalid measurement values."""

    with pytest.raises(TypeError, match="InvalidMeasurementAnalysis.measurement must be"):

        class InvalidMeasurementAnalysis(ScalarMeasurementAnalysis):
            name: ClassVar[str] = "invalid_measurement"
            Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
            measurement: ClassVar[object] = object()


def test_scalar_measurement_analysis_rejects_abstract_measurement_class() -> None:
    """Scalar measurement analyses should reject abstract measurement classes."""

    class AbstractToyMeasurement(ScalarMeasurement):
        """Abstract measurement with metric metadata but no implementation."""

        name: ClassVar[str] = "abstract_toy_measurement"
        metric: ClassVar[MetricSpec] = MetricSpec(name="abstract_toy_value")

    with pytest.raises(TypeError, match="concrete ScalarMeasurement class"):

        class AbstractMeasurementAnalysis(ScalarMeasurementAnalysis):
            name: ClassVar[str] = "abstract_measurement"
            Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
            measurement: ClassVar[type[ScalarMeasurement]] = AbstractToyMeasurement


def test_scalar_measurement_analysis_requires_metric() -> None:
    """Scalar measurement analyses should fail clearly when metric is missing."""

    class MissingMetricMeasurement(ScalarMeasurement):
        """Scalar measurement without metric metadata."""

        name: ClassVar[str] = "missing_metric_measurement"

        def measure(
            self,
            universe,
            settings: BaseModel,
            *,
            start: int | None = None,
            stop: int | None = None,
            step: int | None = None,
        ) -> float:
            """Return a scalar value for the invalid metric contract test."""

            del universe, settings, start, stop, step
            return 1.0

    with pytest.raises(TypeError, match="measurement.metric must be a MetricSpec"):

        class MissingMetricAnalysis(ScalarMeasurementAnalysis):
            name: ClassVar[str] = "missing_metric"
            Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
            measurement: ClassVar[type[ScalarMeasurement]] = MissingMetricMeasurement


def test_scalar_measurement_analysis_rejects_wrong_metric_type() -> None:
    """Scalar measurement analyses should reject non-MetricSpec metric metadata."""

    class WrongMetricTypeMeasurement(ToyScalarMeasurement):
        """Scalar measurement with invalid metric metadata type."""

        metric: ClassVar[object] = object()

    with pytest.raises(TypeError, match="measurement.metric must be a MetricSpec"):

        class WrongMetricTypeAnalysis(ScalarMeasurementAnalysis):
            name: ClassVar[str] = "wrong_metric_type"
            Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
            measurement: ClassVar[type[ScalarMeasurement]] = WrongMetricTypeMeasurement


def test_scalar_measurement_analysis_rejects_blank_metric_name() -> None:
    """Scalar measurement analyses should reject blank scalar metric names."""

    class BlankMetricNameMeasurement(ToyScalarMeasurement):
        """Scalar measurement with a blank metric name."""

        metric: ClassVar[MetricSpec] = MetricSpec(name="   ")

    with pytest.raises(TypeError, match="measurement.metric must be a MetricSpec"):

        class BlankMetricNameAnalysis(ScalarMeasurementAnalysis):
            name: ClassVar[str] = "blank_metric_name"
            Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
            measurement: ClassVar[type[ScalarMeasurement]] = BlankMetricNameMeasurement


def test_scalar_measurement_fingerprint_includes_measurement_version() -> None:
    """Changing measurement version should change aggregate cache identity."""

    settings = ToyScalarSettings(offset=10.0)

    old_fingerprint = ToyScalarAnalysis().aggregate_settings_fingerprint(settings)
    new_fingerprint = VersionedScalarAnalysis().aggregate_settings_fingerprint(settings)

    assert old_fingerprint != new_fingerprint


def test_scalar_measurement_analysis_run_replicate_and_aggregate(tmp_path) -> None:
    """Scalar measurement adapter should execute and aggregate via the runner seam."""

    condition = Condition(
        label="Toy",
        config_path=tmp_path / "config.yaml",
        replicates=(1, 2),
        sim_config=object(),
    )
    settings = ToyScalarSettings(offset=10.0)
    replicate_ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = ToyScalarAnalysis()
    analysis._trajectory_loader_factory = lambda: FakeLoader  # type: ignore[method-assign]
    analysis.get_trajectory_window = lambda *args: FakeWindow()  # type: ignore[method-assign]

    replicate_result = analysis.run_replicate(replicate_ctx, replicate=1)

    assert replicate_result["value"] == 17.0
    assert replicate_result["cache_identity"]["measurement_name"] == "toy_measurement"

    aggregate_ctx = AggregateContext(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        equilibration="0ns",
        settings=settings,
    )
    aggregate_result = analysis.aggregate(
        aggregate_ctx,
        [replicate_result, {**replicate_result, "replicate": 2, "value": 19.0}],
    )

    assert aggregate_result["mean_value"] == 18.0
    assert aggregate_result["replicate_values"] == [17.0, 19.0]
    assert aggregate_result["replicates"] == [1, 2]
    assert aggregate_result["settings_fingerprint"] is not None

    metrics = analysis.extract_metrics(aggregate_result)
    assert list(metrics) == ["toy_value"]
    assert metrics["toy_value"].mean == 18.0
    assert metrics["toy_value"].replicate_values == [17.0, 19.0]


def test_scalar_measurement_analysis_cache_hit_skips_loader(tmp_path) -> None:
    """Compatible scalar canonical caches should avoid trajectory loading."""

    condition = Condition(
        label="Toy",
        config_path=tmp_path / "config.yaml",
        replicates=(1,),
        sim_config=object(),
    )
    settings = ToyScalarSettings(offset=10.0)
    result_path = tmp_path / "run_1" / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    cached_result = build_scalar_replicate_result(ToyScalarAnalysis(), settings)
    result_path.write_text(json.dumps(cached_result))
    replicate_ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=False,
        settings=settings,
        result_path=result_path,
    )

    analysis = ToyScalarAnalysis()
    analysis._trajectory_loader_factory = pytest.fail  # type: ignore[method-assign]

    result = analysis.run_replicate(replicate_ctx, replicate=1)

    assert result == cached_result


@pytest.mark.parametrize(
    ("mutation", "analysis_cls", "settings"),
    [
        ("identity", VersionedScalarAnalysis, ToyScalarSettings(offset=10.0)),
        ("settings", ToyScalarAnalysis, ToyScalarSettings(offset=12.0)),
        ("replicate", ToyScalarAnalysis, ToyScalarSettings(offset=10.0)),
    ],
)
def test_scalar_measurement_analysis_stale_cache_recomputes(
    tmp_path,
    mutation: str,
    analysis_cls: type[ScalarMeasurementAnalysis],
    settings: ToyScalarSettings,
) -> None:
    """Stale scalar canonical caches should fall through to fresh computation."""

    condition = Condition(
        label="Toy",
        config_path=tmp_path / "config.yaml",
        replicates=(1,),
        sim_config=object(),
    )
    result_path = tmp_path / "run_1" / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    cached_result = build_scalar_replicate_result(
        ToyScalarAnalysis(),
        ToyScalarSettings(offset=10.0),
        replicate=2 if mutation == "replicate" else 1,
        value=-1.0,
    )
    result_path.write_text(json.dumps(cached_result))
    replicate_ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=False,
        settings=settings,
        result_path=result_path,
    )

    analysis = analysis_cls()
    analysis._trajectory_loader_factory = lambda: FakeLoader  # type: ignore[method-assign]
    analysis.get_trajectory_window = lambda *args: FakeWindow()  # type: ignore[method-assign]

    result = analysis.run_replicate(replicate_ctx, replicate=1)

    assert result["replicate"] == 1
    assert result["value"] == settings.offset + 7.0
    assert result["value"] != cached_result["value"]


def test_scalar_measurement_analysis_extract_metrics_delegates_to_measurement() -> None:
    """Analysis metric extraction should use the measurement summary hook."""

    class CustomSummaryMeasurement(ToyScalarMeasurement):
        """Toy measurement with a custom aggregate summary schema."""

        metric: ClassVar[MetricSpec] = MetricSpec(name="custom_summary_value")

        def metric_value_from_summary(self, summary: Any) -> MetricValue:
            """Extract a metric from a custom summary object."""

            return self.metric.to_metric_value(
                mean=float(summary.custom_mean),
                sem=0.0,
                replicate_values=[float(summary.custom_mean)],
            )

    class CustomSummaryAnalysis(ScalarMeasurementAnalysis):
        """Toy scalar analysis using a custom summary hook."""

        name: ClassVar[str] = "custom_summary"
        Settings: ClassVar[type[BaseModel]] = ToyScalarSettings
        measurement: ClassVar[type[ScalarMeasurement]] = CustomSummaryMeasurement

    class CustomSummary:
        """Custom aggregate summary used by the measurement hook."""

        custom_mean = 42.0

    metrics = CustomSummaryAnalysis().extract_metrics(CustomSummary())

    assert list(metrics) == ["custom_summary_value"]
    assert metrics["custom_summary_value"].mean == 42.0


def test_scalar_measurement_aggregate_rejects_missing_cache_identity(tmp_path) -> None:
    """Aggregation should reject scalar replicates without cache identity metadata."""

    settings = ToyScalarSettings(offset=10.0)
    result = build_scalar_replicate_result(ToyScalarAnalysis(), settings)
    result.pop("cache_identity")

    with pytest.raises(ValueError, match="cache_identity is missing"):
        ToyScalarAnalysis().aggregate(build_aggregate_context(tmp_path, settings), [result])


def test_scalar_measurement_aggregate_rejects_empty_results(tmp_path) -> None:
    """Aggregation should reject empty scalar replicate collections."""

    with pytest.raises(ValueError, match="requires at least one result"):
        ToyScalarAnalysis().aggregate(
            build_aggregate_context(tmp_path, ToyScalarSettings(offset=10.0)),
            [],
        )


def test_scalar_measurement_aggregate_rejects_wrong_replicate_id(tmp_path) -> None:
    """Aggregation should reject scalar results with unexpected replicate IDs."""

    settings = ToyScalarSettings(offset=10.0)
    result = build_scalar_replicate_result(ToyScalarAnalysis(), settings, replicate=2)

    with pytest.raises(ValueError, match="replicate mismatch"):
        ToyScalarAnalysis().aggregate(build_aggregate_context(tmp_path, settings), [result])


def test_scalar_measurement_aggregate_rejects_changed_measurement_version(tmp_path) -> None:
    """Aggregation should reject replicates from an older measurement version."""

    settings = ToyScalarSettings(offset=10.0)
    stale_result = build_scalar_replicate_result(ToyScalarAnalysis(), settings)

    with pytest.raises(ValueError, match="cache_identity.version mismatch"):
        VersionedScalarAnalysis().aggregate(
            build_aggregate_context(tmp_path, settings), [stale_result]
        )


def test_scalar_measurement_aggregate_rejects_changed_metric_name(tmp_path) -> None:
    """Aggregation should reject replicates from a different metric identity."""

    settings = ToyScalarSettings(offset=10.0)
    stale_result = build_scalar_replicate_result(ToyScalarAnalysis(), settings)

    with pytest.raises(ValueError, match="metric mismatch"):
        RenamedMetricScalarAnalysis().aggregate(
            build_aggregate_context(tmp_path, settings),
            [stale_result],
        )


def test_scalar_measurement_aggregate_rejects_changed_settings_payload(tmp_path) -> None:
    """Aggregation should reject replicates computed with different settings."""

    old_settings = ToyScalarSettings(offset=10.0)
    new_settings = ToyScalarSettings(offset=12.0)
    stale_result = build_scalar_replicate_result(ToyScalarAnalysis(), old_settings)

    with pytest.raises(ValueError, match="settings payload mismatch"):
        ToyScalarAnalysis().aggregate(
            build_aggregate_context(tmp_path, new_settings), [stale_result]
        )
