"""Tests for shared aggregation utilities."""

from __future__ import annotations

import json

import pytest

from polyzymd.analyses.shared.aggregation import collect_replicate_results


def test_collect_replicate_results_warns_and_skips_recoverable_errors(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Recoverable errors should be warned and skipped during aggregation."""

    def _compute(*, replicate: int):
        if replicate == 2:
            raise FileNotFoundError("missing trajectory")
        if replicate == 3:
            raise json.JSONDecodeError("bad json", "{", 0)
        if replicate == 4:
            raise KeyError("missing key")
        if replicate == 5:
            raise ValueError("malformed payload")
        return f"ok-{replicate}"

    with caplog.at_level("WARNING"):
        collection = collect_replicate_results(
            _compute,
            replicates=[1, 2, 3, 4, 5, 6],
            min_replicates=2,
        )

    assert collection.results == ["ok-1", "ok-6"]
    assert collection.successful_replicates == [1, 6]
    assert collection.failed_replicates == [2, 3, 4, 5]
    assert "recoverable error" in caplog.text
    assert "Aggregating 2 of 6 requested replicates" in caplog.text


def test_collect_replicate_results_propagates_unexpected_errors() -> None:
    """Unexpected runtime errors should not be silently downgraded to warnings."""

    def _compute(*, replicate: int):
        if replicate == 2:
            raise TypeError("unexpected type error")
        return f"ok-{replicate}"

    with pytest.raises(TypeError, match="unexpected type error"):
        collect_replicate_results(
            _compute,
            replicates=[1, 2, 3],
            min_replicates=2,
        )
