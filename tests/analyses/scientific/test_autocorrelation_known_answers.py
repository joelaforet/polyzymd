"""Known-answer tests for the statistical inefficiency estimator.

Two processes have an analytic statistical inefficiency. White noise has
g = 1 because successive samples are independent. A first-order
autoregressive process x[t] = phi * x[t-1] + e[t] has a geometric
autocorrelation function C(t) = phi**t, so

    g = 1 + 2 * sum_{t>=1} phi**t = (1 + phi) / (1 - phi).

Every estimator in ``polyzymd.analyses.shared.autocorrelation`` is checked
against those two answers on series drawn from a fixed seed.
"""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses.shared.autocorrelation import (
    estimate_correlation_time,
    n_effective,
    statistical_inefficiency,
)

SEED = 20240501
N_SAMPLES = 20_000


def _white_noise(n: int = N_SAMPLES) -> np.ndarray:
    """Return iid standard normal samples."""

    rng = np.random.default_rng(SEED)
    return rng.standard_normal(n)


def _ar1(phi: float, n: int = N_SAMPLES) -> np.ndarray:
    """Return a stationary AR(1) series with lag-one correlation ``phi``."""

    rng = np.random.default_rng(SEED + int(round(phi * 100)))
    noise = rng.standard_normal(n)
    series = np.empty(n, dtype=np.float64)
    series[0] = noise[0] / np.sqrt(1.0 - phi**2)
    for index in range(1, n):
        series[index] = phi * series[index - 1] + noise[index]
    return series


def _true_g(phi: float) -> float:
    """Return the analytic statistical inefficiency of an AR(1) process."""

    return (1.0 + phi) / (1.0 - phi)


def test_white_noise_integration_method_gives_g_near_one() -> None:
    """Uncorrelated data must give g close to 1, not the old floor of 3."""

    series = _white_noise()

    result = estimate_correlation_time(series, method="integration")

    assert 0.9 <= result.statistical_inefficiency <= 1.1


def test_white_noise_statistical_inefficiency_gives_g_near_one() -> None:
    """The pymbar-style estimator must also give g close to 1."""

    assert 0.9 <= statistical_inefficiency(_white_noise()) <= 1.1


def test_white_noise_tau_is_far_below_one_timestep() -> None:
    """tau must be free to fall below one timestep for uncorrelated data."""

    result = estimate_correlation_time(_white_noise(), timestep=2.0, timestep_unit="ps")

    assert result.tau < 0.5 * 2.0
    assert result.tau_unit == "ps"


@pytest.mark.parametrize("phi", [0.5, 0.9])
def test_ar1_integration_method_matches_analytic_g(phi: float) -> None:
    """The integration method must match (1+phi)/(1-phi) within 15 percent."""

    series = _ar1(phi)
    expected = _true_g(phi)

    result = estimate_correlation_time(series, method="integration")

    assert result.statistical_inefficiency == pytest.approx(expected, rel=0.15)


@pytest.mark.parametrize("phi", [0.5, 0.9])
def test_ar1_statistical_inefficiency_matches_analytic_g(phi: float) -> None:
    """The pymbar-style estimator must match the analytic answer too."""

    assert statistical_inefficiency(_ar1(phi)) == pytest.approx(_true_g(phi), rel=0.15)


@pytest.mark.parametrize("phi", [0.5, 0.9])
def test_integration_method_agrees_with_statistical_inefficiency(phi: float) -> None:
    """Both entry points must return the same number for the same series."""

    series = _ar1(phi)

    result = estimate_correlation_time(series, method="integration")

    assert result.statistical_inefficiency == pytest.approx(
        statistical_inefficiency(series), rel=1e-12
    )


@pytest.mark.parametrize("phi", [0.5, 0.9])
def test_tau_and_g_obey_the_documented_relation(phi: float) -> None:
    """tau = (g - 1) / 2 * dt must invert to g = 1 + 2 tau / dt."""

    result = estimate_correlation_time(_ar1(phi), timestep=5.0, timestep_unit="ps")

    assert 1.0 + 2.0 * result.tau / 5.0 == pytest.approx(result.statistical_inefficiency)


@pytest.mark.parametrize("phi", [0.0, 0.5, 0.9])
def test_n_independent_is_n_over_g_without_truncation(phi: float) -> None:
    """n_independent must be the float N/g, not an integer-truncated count."""

    series = _white_noise() if phi == 0.0 else _ar1(phi)

    result = estimate_correlation_time(series, method="integration")

    expected = n_effective(series.size, result.statistical_inefficiency)
    assert isinstance(result.n_independent, float)
    assert result.n_independent == pytest.approx(expected, rel=1e-12)
    assert abs(result.n_independent - expected) < 1.0
