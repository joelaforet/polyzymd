"""Known-answer tests for the statistical inefficiency estimator.

Two processes have an analytic statistical inefficiency. White noise has
g = 1 because successive samples are independent. A first-order
autoregressive process x[t] = phi * x[t-1] + e[t] has a geometric
autocorrelation function C(t) = phi**t, so

    g = 1 + 2 * sum_{t>=1} phi**t = (1 + phi) / (1 - phi).

``statistical_inefficiency`` is checked against those two answers on series
drawn from a fixed seed.
"""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses.shared.autocorrelation import statistical_inefficiency

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


def test_white_noise_statistical_inefficiency_gives_g_near_one() -> None:
    """The pymbar-style estimator must also give g close to 1."""

    assert 0.9 <= statistical_inefficiency(_white_noise()) <= 1.1


@pytest.mark.parametrize("phi", [0.5, 0.9])
def test_ar1_statistical_inefficiency_matches_analytic_g(phi: float) -> None:
    """The pymbar-style estimator must match the analytic answer too."""

    assert statistical_inefficiency(_ar1(phi)) == pytest.approx(_true_g(phi), rel=0.15)
