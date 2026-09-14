"""Tests for shared autocorrelation utilities."""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses.shared.autocorrelation import (
    ACFResult,
    compute_acf,
    estimate_correlation_time,
)


def test_compute_acf_constant_series_returns_degenerate_acf() -> None:
    """Constant timeseries should return a stable degenerate ACF."""
    series = np.ones(100, dtype=np.float64)

    result = compute_acf(series, max_lag=10, timestep=2.0, timestep_unit="ps")

    assert result.n_samples == 100
    assert result.acf[0] == 1.0
    assert np.allclose(result.acf[1:], 0.0)
    assert np.allclose(result.lags, np.arange(11, dtype=np.float64) * 2.0)


def test_acf_result_to_dict_includes_n_samples() -> None:
    """Serialized ACF payload should include n_samples metadata."""
    result = ACFResult(
        lags=np.array([0.0, 1.0], dtype=np.float64),
        acf=np.array([1.0, 0.5], dtype=np.float64),
        timestep=1.0,
        timestep_unit="frames",
        n_samples=123,
    )

    payload = result.to_dict()

    assert payload["n_samples"] == 123


def test_estimate_correlation_time_rejects_withdrawn_methods() -> None:
    """The first_zero and exponential_fit estimators were removed."""
    series = np.linspace(0.0, 1.0, 50) + np.tile([0.0, 0.1], 25)

    for method in ("first_zero", "exponential_fit"):
        with pytest.raises(ValueError, match="no longer supported"):
            estimate_correlation_time(series, method=method)


def test_estimate_correlation_time_rejects_unknown_method() -> None:
    """An unrecognised method name is still an error."""
    series = np.linspace(0.0, 1.0, 50) + np.tile([0.0, 0.1], 25)

    with pytest.raises(ValueError, match="Unknown method"):
        estimate_correlation_time(series, method="block_average")
