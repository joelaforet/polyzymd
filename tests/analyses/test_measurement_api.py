"""Tests for the lightweight measurement analysis API."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

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


def test_metric_spec_converts_to_metric_value() -> None:
    """Metric metadata should convert to the existing comparison model."""

    spec = ToyScalarMeasurement.metric
    metric_value = spec.to_metric_value(mean=2.0, sem=0.5, replicate_values=[1.0, 3.0])

    assert isinstance(metric_value, MetricValue)
    assert metric_value.name == "toy_value"
    assert metric_value.higher_is_better is False
    assert metric_value.direction_labels == ("lower", "unchanged", "higher")


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
