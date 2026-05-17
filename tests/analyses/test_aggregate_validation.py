"""Tests for framework aggregate validation hooks."""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar

import pytest
from pydantic import BaseModel

from polyzymd.analyses._framework.cache_identity import settings_fingerprint
from polyzymd.analyses.base import AggregateValidationError, Analysis, Condition


class DummySettings(BaseModel):
    """Settings model used by aggregate validation tests."""

    cutoff: float = 4.5


class DummyAggregate(BaseModel):
    """Aggregate model with the standard validation metadata."""

    settings_fingerprint: str | None = None
    equilibration_time: float = 10.0
    equilibration_unit: str = "ns"
    replicates: list[int]
    n_replicates: int


class CountOnlyAggregate(BaseModel):
    """Aggregate model with only replicate count metadata."""

    settings_fingerprint: str | None = None
    n_replicates: int


class DummyAnalysis(Analysis):
    """Minimal analysis implementation for hook tests."""

    name: ClassVar[str] = "dummy"
    Settings: ClassVar[type] = DummySettings
    AggregatedResultClass: ClassVar[type] = DummyAggregate
    has_compute_stage: ClassVar[bool] = False
    has_aggregate_stage: ClassVar[bool] = False


class CountOnlyAnalysis(DummyAnalysis):
    """Analysis using an aggregate model without replicate IDs."""

    name: ClassVar[str] = "count_only"
    AggregatedResultClass: ClassVar[type] = CountOnlyAggregate


def _condition() -> Condition:
    """Return a lightweight condition for validation tests."""

    return Condition(
        label="condition",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1, 2, 3),
        sim_config=object(),
    )


def _payload(settings: DummySettings, **overrides: Any) -> dict[str, Any]:
    """Return a valid aggregate payload with optional overrides."""

    payload: dict[str, Any] = {
        "settings_fingerprint": settings_fingerprint(settings),
        "equilibration_time": 10.0,
        "equilibration_unit": "ns",
        "replicates": [1, 2],
        "n_replicates": 2,
    }
    payload.update(overrides)
    return payload


def test_validation_coerces_dict_through_result_class() -> None:
    """Dict aggregates should be parsed through AggregatedResultClass."""

    settings = DummySettings()
    analysis = DummyAnalysis()

    result = analysis.validate_aggregated_result(
        _payload(settings),
        condition=_condition(),
        settings=settings,
        equilibration="10ns",
        expected_replicates=(1, 2),
    )

    assert isinstance(result, DummyAggregate)
    assert result.replicates == [1, 2]


def test_validation_rejects_settings_mismatch() -> None:
    """Aggregates from different settings should be rejected."""

    settings = DummySettings()
    analysis = DummyAnalysis()

    with pytest.raises(AggregateValidationError, match="settings fingerprint mismatch"):
        analysis.validate_aggregated_result(
            _payload(settings, settings_fingerprint="deadbeef"),
            condition=_condition(),
            settings=settings,
            equilibration="10ns",
            expected_replicates=(1, 2),
        )


def test_validation_rejects_missing_fingerprint() -> None:
    """Legacy aggregates without settings identity should be rejected."""

    settings = DummySettings()
    analysis = DummyAnalysis()

    with pytest.raises(AggregateValidationError, match="missing settings fingerprint"):
        analysis.validate_aggregated_result(
            _payload(settings, settings_fingerprint=None),
            condition=_condition(),
            settings=settings,
            equilibration="10ns",
            expected_replicates=(1, 2),
        )


def test_validation_rejects_replicate_mismatch() -> None:
    """Exact replicate validation should reject stale aggregate coverage."""

    settings = DummySettings()
    analysis = DummyAnalysis()

    with pytest.raises(AggregateValidationError, match="replicate mismatch"):
        analysis.validate_aggregated_result(
            _payload(settings, replicates=[1, 3]),
            condition=_condition(),
            settings=settings,
            equilibration="10ns",
            expected_replicates=(1, 2),
        )


def test_validation_allows_replicate_subset_when_enabled() -> None:
    """Subset mode should accept successful finalized replicate subsets."""

    settings = DummySettings()
    analysis = DummyAnalysis()

    result = analysis.validate_aggregated_result(
        _payload(settings, replicates=[1, 2]),
        condition=_condition(),
        settings=settings,
        equilibration="10ns",
        expected_replicates=(1, 2, 3),
        allow_replicate_subset=True,
    )

    assert result.n_replicates == 2


def test_validation_uses_count_when_ids_are_absent() -> None:
    """Aggregates with only n_replicates should validate exact counts."""

    settings = DummySettings()
    analysis = CountOnlyAnalysis()

    result = analysis.validate_aggregated_result(
        {
            "settings_fingerprint": settings_fingerprint(settings),
            "n_replicates": 2,
        },
        condition=_condition(),
        settings=settings,
        equilibration="10ns",
        expected_replicates=(1, 2),
    )

    assert isinstance(result, CountOnlyAggregate)
