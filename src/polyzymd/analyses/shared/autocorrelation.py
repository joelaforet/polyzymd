"""Autocorrelation analysis for independent sampling.

MD trajectories are highly correlated in time - consecutive frames are not
independent samples. This module provides tools to:

1. Compute the autocorrelation function (ACF) of an observable
2. Estimate the correlation time (τ) from the ACF
3. Compute statistical inefficiency (g) for proper uncertainty quantification

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

Method for τ estimation
-----------------------
One estimator is offered. The statistical inefficiency is summed directly
from the normalised ACF over positive lags,

    g = 1 + 2*Σ_{t>=1} C(t)*(1 - t/N),

with the sum truncated at the first non-positive C(t) beyond a short minimum
lag, and the integrated correlation time follows from it as
τ = (g - 1)/2 * Δt. Earlier releases also offered a first-zero-crossing
estimator and an exponential fit. Both were withdrawn because neither
estimates the integrated correlation time. On an AR(1) series with φ = 0.9
the first-zero method overestimated g by a factor of five and the
exponential fit by up to two.

Statistical Validity
--------------------
The number of effective independent samples (N_eff) is computed as:
    N_eff = N / g = N / (1 + 2*Σ C(t)*(1-t/N))

N_eff is a real number and is reported as one. It is not rounded down to a
whole frame count.

For multiple timeseries of different lengths (e.g., replicates), use
`statistical_inefficiency_multiple()` which correctly handles the averaging.

`MIN_RECOMMENDED_N_INDEPENDENT` is a convention of this package, not a
threshold taken from any reference below. Grossfield et al. (2018) give no
such cutoff; they discuss small sample counts in the context of plotting
every replicate rather than a summary statistic.

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
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

logger = logging.getLogger(__name__)

# Minimum recommended independent samples for reliable statistics.
# House convention of this package; see the module docstring.
MIN_RECOMMENDED_N_INDEPENDENT = 10

# Lags shorter than this are always summed, so that a single noisy negative
# value near lag zero cannot truncate the sum. pymbar's default.
DEFAULT_MINTIME = 3

_WITHDRAWN_METHODS = ("first_zero", "exponential_fit")


@dataclass
class ACFResult:
    """Result of autocorrelation function computation.

    Attributes
    ----------
    lags : NDArray[np.float64]
        Time lags in the same units as timestep
    acf : NDArray[np.float64]
        Autocorrelation values (normalized, ACF[0] = 1)
    timestep : float
        Time between frames
    timestep_unit : str
        Unit of timestep (e.g., "ps", "ns")
    n_samples : int
        Number of samples in the original timeseries
    """

    lags: NDArray[np.float64]
    acf: NDArray[np.float64]
    timestep: float
    timestep_unit: str
    n_samples: int

    def __len__(self) -> int:
        return len(self.lags)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "lags": self.lags.tolist(),
            "acf": self.acf.tolist(),
            "timestep": self.timestep,
            "timestep_unit": self.timestep_unit,
            "n_samples": self.n_samples,
        }


@dataclass
class CorrelationTimeResult:
    """Result of correlation time estimation.

    Attributes
    ----------
    tau : float
        Estimated correlation time
    tau_unit : str
        Unit of tau (same as timestep unit)
    method : str
        Method used for estimation
    n_independent : float
        Estimated number of independent samples in the trajectory, N/g. This
        is a real number and is deliberately not rounded down to whole frames.
    statistical_inefficiency : float
        g = 1 + 2*tau/dt, factor by which variance is inflated
    warning : str | None
        Warning message if statistics may be unreliable (e.g., N_ind < 10)
    """

    tau: float
    tau_unit: str
    method: str
    n_independent: float
    statistical_inefficiency: float
    warning: str | None = None

    @property
    def is_reliable(self) -> bool:
        """Return True if statistics are likely reliable (N_ind >= 10)."""
        return self.n_independent >= MIN_RECOMMENDED_N_INDEPENDENT

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "tau": self.tau,
            "tau_unit": self.tau_unit,
            "method": self.method,
            "n_independent": self.n_independent,
            "statistical_inefficiency": self.statistical_inefficiency,
            "warning": self.warning,
            "is_reliable": self.is_reliable,
        }


def compute_acf(
    timeseries: ArrayLike,
    max_lag: int | None = None,
    timestep: float = 1.0,
    timestep_unit: str = "frames",
) -> ACFResult:
    """Compute autocorrelation function of a 1D timeseries.

    Uses FFT-based computation for efficiency.

    Parameters
    ----------
    timeseries : array_like
        1D array of values (e.g., RMSD over time, distance over time)
    max_lag : int, optional
        Maximum lag to compute (in frames). Default is N//4 where N is
        the length of the timeseries.
    timestep : float, optional
        Time between frames. Default is 1.0.
    timestep_unit : str, optional
        Unit of timestep. Default is "frames".

    Returns
    -------
    ACFResult
        Container with lags, acf values, and metadata

    Examples
    --------
    >>> # Compute ACF of RMSD timeseries
    >>> rmsd = np.array([1.2, 1.3, 1.25, 1.4, ...])  # from MDAnalysis
    >>> acf_result = compute_acf(rmsd, timestep=10.0, timestep_unit="ps")
    >>> print(f"ACF at lag 100ps: {acf_result.acf[10]:.3f}")

    Notes
    -----
    The ACF is normalized so that ACF[0] = 1.

    For a stationary process: ACF(τ) = <(x(t) - μ)(x(t+τ) - μ)> / σ²

    For constant or near-constant timeseries (variance below a small epsilon),
    this function returns a defined degenerate ACF with ACF[0] = 1 and all
    positive lags set to 0.
    """
    x = np.asarray(timeseries, dtype=np.float64)
    n = len(x)

    if n < 10:
        raise ValueError(f"Timeseries too short ({n} points). Need at least 10.")

    if max_lag is None:
        max_lag = n // 4  # Reasonable default
    max_lag = min(max_lag, n - 1)

    # Remove mean
    x_centered = x - np.mean(x)
    variance = float(np.var(x_centered))

    # Define a stable degenerate ACF for near-constant timeseries
    if variance < 1e-12:
        acf = np.zeros(max_lag + 1, dtype=np.float64)
        acf[0] = 1.0
        lags = np.arange(max_lag + 1, dtype=np.float64) * timestep
        return ACFResult(
            lags=lags,
            acf=acf,
            timestep=timestep,
            timestep_unit=timestep_unit,
            n_samples=n,
        )

    # FFT-based autocorrelation (much faster than direct computation)
    # Pad to next power of 2 for FFT efficiency
    n_fft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
    fft_x = np.fft.fft(x_centered, n_fft)
    acf_full = np.fft.ifft(fft_x * np.conj(fft_x)).real[:n]

    # Normalize by decreasing sample size and variance
    acf_full = acf_full / (np.arange(n, 0, -1) * variance)

    # Take only up to max_lag
    acf = acf_full[: max_lag + 1]
    lags = np.arange(max_lag + 1) * timestep

    return ACFResult(
        lags=lags,
        acf=acf,
        timestep=timestep,
        timestep_unit=timestep_unit,
        n_samples=n,
    )


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


def estimate_correlation_time(
    timeseries: ArrayLike,
    timestep: float = 1.0,
    timestep_unit: str = "frames",
    method: Literal["integration"] = "integration",
    n_frames: int | None = None,
) -> CorrelationTimeResult:
    """Estimate the integrated correlation time of a raw timeseries.

    Parameters
    ----------
    timeseries : array_like
        Per-frame values of one observable.
    timestep : float, optional
        Time between frames.
    timestep_unit : str, optional
        Unit of ``timestep``.
    method : {"integration"}
        Only "integration" is supported: g is summed from the normalised ACF
        over positive lags by pymbar and tau = (g - 1)/2 * dt. The former
        "first_zero" and "exponential_fit" values raise ValueError.
    n_frames : int, optional
        Frame count used for ``n_independent``; defaults to the series length.

    Returns
    -------
    CorrelationTimeResult
        tau, g, N/g and a warning when N/g is below the recommended minimum.

    Raises
    ------
    ValueError
        If ``method`` is a withdrawn estimator or unknown.

    Notes
    -----
    For an uncorrelated series g approaches 1 and tau approaches 0. Neither
    is floored at one timestep, so N_eff approaches N as it should.
    """
    if method in _WITHDRAWN_METHODS:
        raise ValueError(
            f"Correlation-time method {method!r} is no longer supported because it does "
            "not estimate the integrated correlation time. Use method='integration', "
            "which sums the normalised ACF per Chodera et al. (2007)."
        )
    if method != "integration":
        raise ValueError(f"Unknown method: {method}")

    series = np.asarray(timeseries, dtype=np.float64)
    if n_frames is None:
        n_frames = int(series.size)
    g = statistical_inefficiency(series)
    tau = 0.5 * (g - 1.0) * timestep
    n_independent = n_effective(n_frames, g)

    warning = None
    if n_independent < MIN_RECOMMENDED_N_INDEPENDENT:
        warning = (
            f"Low statistical reliability: only {n_independent:.1f} independent samples "
            f"(recommended >= {MIN_RECOMMENDED_N_INDEPENDENT}). "
            f"Correlation time tau = {tau:.1f} {timestep_unit} is comparable to or longer "
            "than the trajectory sampling window. Consider: (1) extending simulation time, "
            "(2) using multiple independent trajectories, or (3) interpreting results "
            "with caution. See Grossfield et al. (2018) LiveCoMS 1:5067."
        )
        logger.warning(warning)

    return CorrelationTimeResult(
        tau=tau,
        tau_unit=timestep_unit,
        method=method,
        n_independent=n_independent,
        statistical_inefficiency=g,
        warning=warning,
    )


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


def statistical_inefficiency_multiple(timeseries_list: list[ArrayLike]) -> float:
    """Statistical inefficiency pooled over several series, from pymbar.

    The series may differ in length. pymbar computes one global mean and
    averages the ACF over every series before summing, which is what is
    needed when replicates contribute different frame counts.
    """
    arrays = [np.asarray(ts, dtype=np.float64) for ts in timeseries_list]
    if sum(a.size for a in arrays) < 3:
        return 1.0
    if np.var(np.concatenate(arrays)) < 1e-10:
        return 1.0
    return float(_pymbar_timeseries().statistical_inefficiency_multiple(arrays))


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


def check_statistical_reliability(
    n_eff: float,
    threshold: int = MIN_RECOMMENDED_N_INDEPENDENT,
) -> tuple[bool, str | None]:
    """Check if statistics are reliable based on effective sample count.

    Parameters
    ----------
    n_eff : float
        Number of effective independent samples
    threshold : int
        Minimum recommended independent samples. Default is 10.

    Returns
    -------
    is_reliable : bool
        True if n_eff >= threshold
    warning : str | None
        Warning message if not reliable, None otherwise

    Examples
    --------
    >>> g = statistical_inefficiency(contacts)
    >>> n_eff = n_effective(len(contacts), g)
    >>> is_ok, warning = check_statistical_reliability(n_eff)
    >>> if not is_ok:
    ...     print(warning)
    """
    if n_eff >= threshold:
        return True, None

    warning = (
        f"Low statistical reliability: only {n_eff:.1f} effective independent samples "
        f"(recommended >= {threshold}). Consider: (1) extending simulation time, "
        f"(2) using more independent replicates, or (3) interpreting results "
        f"with caution. See Grossfield et al. (2018) LiveCoMS 1:5067."
    )
    logger.warning(warning)

    return False, warning
