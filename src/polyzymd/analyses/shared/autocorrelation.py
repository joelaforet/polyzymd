"""Autocorrelation analysis for independent sampling.

MD trajectories are highly correlated in time - consecutive frames are not
independent samples. This module wraps pymbar to:

1. Compute the statistical inefficiency (g) of a timeseries
2. Convert g into an effective sample size, N_eff = N/g
3. Detect the first equilibrated sample of a timeseries

Key Concepts
------------
- **Autocorrelation function (ACF)**: Measures how correlated a signal is with
  itself at different time lags. ACF(0) = 1, and ACF decays toward 0.

- **Correlation time (τ)**: Characteristic time for decorrelation. Frames
  separated by > 2τ are approximately independent.

- **Statistical inefficiency (g)**: Factor by which variance is inflated due
  to correlation. g = 1 + 2*Σ C(t)*(1-t/N). N_eff = N/g.

- **Independent samples**: For proper SEM calculation, we need N_eff independent
  samples, not N_frames correlated observations.

Estimator
---------
The statistical inefficiency is summed directly from the normalised ACF over
positive lags,

    g = 1 + 2*Σ_{t>=1} C(t)*(1 - t/N),

with the sum truncated at the first non-positive C(t) beyond a short minimum
lag. The integrated correlation time follows from it as τ = (g - 1)/2 * Δt.

Statistical Validity
--------------------
The number of effective independent samples (N_eff) is computed as:
    N_eff = N / g = N / (1 + 2*Σ C(t)*(1-t/N))

N_eff is a real number and is reported as one. It is not rounded down to a
whole frame count.

Not implemented
---------------
Block averaging (Flyvbjerg and Petersen 1989, J. Chem. Phys. 91:461,
doi:10.1063/1.457480) is an alternative route to g. It is not implemented
anywhere in this package.

References
----------
Chodera, J. D., Swope, W. C., Pitera, J. W., Seok, C., and Dill, K. A. (2007).
    Use of the weighted histogram analysis method for the analysis of simulated
    and parallel tempering simulations. Journal of Chemical Theory and
    Computation, 3(1), 26-41. doi:10.1021/ct0502864
Shirts, M. R., and Chodera, J. D. (2008). Statistically optimal analysis of
    samples from multiple equilibrium states. The Journal of Chemical Physics,
    129(12), 124105. doi:10.1063/1.2978177. The estimator here follows the
    algorithm of pymbar's `timeseries` module, which PolyzyMD depends on.
Janke, W. (2002). Statistical analysis of simulations: data correlations and
    error estimation. In J. Grotendorst, D. Marx, and A. Muramatsu (Eds.),
    Quantum Simulations of Complex Many-Body Systems, NIC Series vol. 10,
    423-445. John von Neumann Institute for Computing.
Grossfield, A., Patrone, P. N., Roe, D. R., Schultz, A. J., Siderius, D. W.,
    and Zuckerman, D. M. (2018). Best practices for quantifying the uncertainty
    in molecular simulation. Living Journal of Computational Molecular Science,
    1(1), 5067. doi:10.33011/livecoms.1.1.5067
"""

from __future__ import annotations

import logging

import numpy as np
from numpy.typing import ArrayLike

logger = logging.getLogger(__name__)

# Lags shorter than this are always summed, so that a single noisy negative
# value near lag zero cannot truncate the sum. pymbar's default.
DEFAULT_MINTIME = 3


def _pymbar_timeseries():
    """Import ``pymbar.timeseries`` on first use.

    The import costs about a second and pymbar logs a banner about JAX being
    absent, so it is deferred to the call and the banner is silenced.
    """
    pymbar_logger = logging.getLogger("pymbar")
    if pymbar_logger.level == logging.NOTSET:
        pymbar_logger.setLevel(logging.ERROR)
    from pymbar import timeseries

    return timeseries


def statistical_inefficiency(
    timeseries: ArrayLike,
    mintime: int = DEFAULT_MINTIME,
    fft: bool = False,
) -> float:
    """Statistical inefficiency g of a timeseries, from pymbar.

    g = 1 + 2 * sum_{t>=1} C(t) * (1 - t/N), summed until the normalised ACF
    first turns non-positive beyond ``mintime`` (Chodera et al. 2007). The
    variance of the mean is Var(x) * g / N, so N_eff = N / g. A series shorter
    than three points or with no variance has g = 1.

    Parameters
    ----------
    timeseries : array_like
        Per-frame values, continuous or 0/1.
    mintime : int
        Lags summed unconditionally before the truncation rule applies.
    fft : bool
        Forwarded to pymbar; the FFT path helps for very long series.
    """
    x = np.asarray(timeseries, dtype=np.float64)
    if x.size < 3:
        logger.warning(f"Timeseries too short ({x.size} points). Returning g=1.0")
        return 1.0
    if np.var(x) < 1e-10:
        return 1.0
    return float(_pymbar_timeseries().statistical_inefficiency(x, mintime=mintime, fft=fft))


def detect_equilibration(x: ArrayLike) -> int:
    """Index of the first equilibrated sample, from pymbar ``detect_equilibration``.

    pymbar picks the start ``t0`` that maximises the number of effective samples
    in ``x[t0:]`` (Chodera 2016, doi:10.1021/acs.jctc.5b00784). Only about 100
    evenly spaced starts are tried, so a long series costs no more than a short
    one and ``t0`` is known to 1 percent of the series length.
    """
    x = np.asarray(x, dtype=np.float64)
    step = max(1, len(x) // 100)
    return int(_pymbar_timeseries().detect_equilibration(x, nskip=step)[0])


def n_effective(n_samples: int, g: float) -> float:
    """Compute number of effective independent samples.

    Parameters
    ----------
    n_samples : int
        Total number of samples
    g : float
        Statistical inefficiency

    Returns
    -------
    float
        Effective number of independent samples (N_eff = N / g)
    """
    if g <= 0:
        return float(n_samples)
    return n_samples / g
