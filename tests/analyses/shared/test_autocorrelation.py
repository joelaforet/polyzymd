"""Tests for shared autocorrelation utilities."""

from __future__ import annotations

import numpy as np

from polyzymd.analyses.shared.autocorrelation import (
    ACFResult,
    _find_first_zero_crossing,
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


def test_estimate_correlation_time_constant_series_returns_tau_zero() -> None:
    """Degenerate ACF should map to tau=0 and g=1."""
    series = np.ones(80, dtype=np.float64)
    acf_result = compute_acf(series, max_lag=20)

    tau = estimate_correlation_time(acf_result)

    assert tau.tau == 0.0
    assert tau.statistical_inefficiency == 1.0
    assert tau.n_independent == 80


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


def test_estimate_correlation_time_uses_acf_result_n_samples() -> None:
    """n_independent should use ACFResult.n_samples, not len(acf)*4."""
    acf_result = ACFResult(
        lags=np.arange(11, dtype=np.float64),
        acf=np.array([1.0, 0.5, 0.0] + [0.0] * 8, dtype=np.float64),
        timestep=1.0,
        timestep_unit="frames",
        n_samples=40,
    )

    tau = estimate_correlation_time(acf_result, method="first_zero")

    # first_zero gives tau=2, so g=5 and n_independent=int(40/5)=8
    assert tau.n_independent == 8


def test_find_first_zero_crossing_returns_first_nonpositive_index() -> None:
    """Zero crossing index should be the first nonpositive lag."""
    acf = np.array([1.0, 0.2, -0.1, -0.2], dtype=np.float64)
    assert _find_first_zero_crossing(acf) == 2

    acf_with_exact_zero = np.array([1.0, 0.0, -0.1], dtype=np.float64)
    assert _find_first_zero_crossing(acf_with_exact_zero) == 1
