"""Tests for shared autocorrelation utilities."""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses.shared.autocorrelation import (
    ACFResult,
    _statistical_inefficiency_from_acf,
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
    assert tau.n_independent == pytest.approx(80.0)


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

    tau = estimate_correlation_time(acf_result)

    # Lag zero is excluded from the sum, so g = 1 + 2*0.5*(1 - 1/40) = 1.975
    # and tau = (g - 1)/2 = 0.4875 frames.
    assert tau.statistical_inefficiency == pytest.approx(1.975)
    assert tau.tau == pytest.approx(0.4875)
    assert tau.n_independent == pytest.approx(40.0 / 1.975)


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


def test_statistical_inefficiency_from_acf_excludes_lag_zero() -> None:
    """The sum must start at lag one; lag zero is the leading 1."""
    acf = np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float64)

    assert _statistical_inefficiency_from_acf(acf, n_frames=100) == pytest.approx(1.0)


def test_statistical_inefficiency_from_acf_truncates_at_first_nonpositive() -> None:
    """Lags past the first non-positive value beyond mintime are dropped."""
    acf = np.array([1.0, 0.4, 0.2, 0.1, -0.05, 0.9, 0.9], dtype=np.float64)

    g = _statistical_inefficiency_from_acf(acf, n_frames=1000, mintime=3)

    expected = 1.0 + 2.0 * (0.4 * 0.999 + 0.2 * 0.998 + 0.1 * 0.997)
    assert g == pytest.approx(expected)
