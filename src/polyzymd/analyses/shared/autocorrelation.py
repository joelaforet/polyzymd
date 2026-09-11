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
    algorithm of pymbar's MIT-licensed `timeseries` module, which accompanies
    that paper, without taking pymbar as a dependency.
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
# value near lag zero cannot truncate the sum. Same default as pymbar.
DEFAULT_MINTIME = 3

# Identifies which correlation estimator produced a stored number. Version "1"
# was the trapezoid integration that counted lag zero twice and floored tau at
# one timestep; version "2" is the pymbar-style sum below. Plugins stamp this on
# every replicate artifact whose sem_* field divides by an effective sample
# count, and refuse to aggregate artifacts that disagree, because the two
# versions give sem values that differ by a factor of about the square root of
# three for a fast observable. Bump it whenever the estimator changes a number.
AUTOCORRELATION_ESTIMATOR_VERSION = "2"

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


def estimate_correlation_time(
    acf_or_timeseries: ACFResult | ArrayLike,
    timestep: float = 1.0,
    timestep_unit: str = "frames",
    method: Literal["integration"] = "integration",
    n_frames: int | None = None,
) -> CorrelationTimeResult:
    """Estimate correlation time from ACF or raw timeseries.

    Parameters
    ----------
    acf_or_timeseries : ACFResult or array_like
        Either an ACFResult from compute_acf(), or a raw timeseries
    timestep : float, optional
        Time between frames (only used if passing raw timeseries)
    timestep_unit : str, optional
        Unit of timestep (only used if passing raw timeseries)
    method : {"integration"}
        Estimator for τ. Only "integration" is supported: g is summed from
        the normalised ACF over positive lags and τ = (g - 1)/2 * Δt. The
        former "first_zero" and "exponential_fit" values raise ValueError.
    n_frames : int, optional
        Total number of frames (for computing n_independent).
        Only needed if passing ACFResult.

    Returns
    -------
    CorrelationTimeResult
        Contains tau, method used, n_independent, statistical_inefficiency

    Raises
    ------
    ValueError
        If ``method`` is one of the withdrawn estimators, or is unknown.

    Examples
    --------
    >>> acf_result = compute_acf(rmsd, timestep=10.0, timestep_unit="ps")
    >>> tau_result = estimate_correlation_time(acf_result, method="integration")
    >>> print(f"Correlation time: {tau_result.tau:.1f} {tau_result.tau_unit}")
    >>> print(f"Independent samples: {tau_result.n_independent:.1f}")

    Notes
    -----
    Passing a raw timeseries is the accurate path, because the sum then runs
    over every available lag. Passing an ACFResult truncates the sum at the
    ACF's own ``max_lag`` (``N // 4`` by default), which underestimates g
    when the correlation time approaches a quarter of the trajectory.

    For an uncorrelated series g approaches 1 and τ approaches 0. Neither is
    floored at one timestep, so N_eff approaches N as it should.
    """
    if method in _WITHDRAWN_METHODS:
        raise ValueError(
            f"Correlation-time method {method!r} is no longer supported because it does "
            "not estimate the integrated correlation time. Use method='integration', "
            "which sums the normalised ACF per Chodera et al. (2007)."
        )
    if method != "integration":
        raise ValueError(f"Unknown method: {method}")

    # Handle input type
    if isinstance(acf_or_timeseries, ACFResult):
        acf_result = acf_or_timeseries
        dt = acf_result.timestep
        unit = acf_result.timestep_unit
        if n_frames is None:
            n_frames = acf_result.n_samples
        g = _statistical_inefficiency_from_acf(acf_result.acf, n_frames)
    else:
        series = np.asarray(acf_or_timeseries, dtype=np.float64)
        dt = timestep
        unit = timestep_unit
        if n_frames is None:
            n_frames = int(series.size)
        g = statistical_inefficiency(series)

    # Integrated correlation time implied by g = 1 + 2*tau/dt
    tau = 0.5 * (g - 1.0) * dt

    # Number of independent samples, kept as a real number
    n_independent = n_effective(n_frames, g)

    # Generate warning if statistics may be unreliable
    warning = None
    if n_independent < MIN_RECOMMENDED_N_INDEPENDENT:
        warning = (
            f"Low statistical reliability: only {n_independent:.1f} independent samples "
            f"(recommended >= {MIN_RECOMMENDED_N_INDEPENDENT}). "
            f"Correlation time τ = {tau:.1f} {unit} is comparable to or longer than "
            f"the trajectory sampling window. Consider: (1) extending simulation time, "
            f"(2) using multiple independent trajectories, or (3) interpreting results "
            f"with caution. See Grossfield et al. (2018) LiveCoMS 1:5067."
        )
        logger.warning(warning)

    return CorrelationTimeResult(
        tau=tau,
        tau_unit=unit,
        method=method,
        n_independent=n_independent,
        statistical_inefficiency=g,
        warning=warning,
    )


def _statistical_inefficiency_from_acf(
    acf: NDArray[np.float64],
    n_frames: int,
    mintime: int = DEFAULT_MINTIME,
) -> float:
    """Sum a normalised ACF into a statistical inefficiency.

    Implements g = 1 + 2*Σ_{t>=1} C(t)*(1 - t/N) with the sum truncated at
    the first non-positive C(t) at lag t >= ``mintime``. Lag zero is excluded
    from the sum; it is the leading 1.

    Parameters
    ----------
    acf : NDArray[np.float64]
        Normalised autocorrelation function, ``acf[0] == 1``, indexed by lag
        in frames.
    n_frames : int
        Number of samples the ACF was computed from, used for the (1 - t/N)
        finite-size weight.
    mintime : int
        Lags below this are summed unconditionally, so that one noisy
        negative value near lag zero cannot truncate the sum.

    Returns
    -------
    float
        Statistical inefficiency g, at least 1.0.
    """
    values = np.asarray(acf, dtype=np.float64)
    n = int(n_frames)
    if values.size < 2 or n < 3:
        return 1.0

    lags = np.arange(1, values.size, dtype=np.float64)
    tail = values[1:]

    # pymbar truncation rule: stop at the first non-positive C(t), t >= mintime
    nonpositive = np.nonzero((tail <= 0.0) & (lags >= float(mintime)))[0]
    end = int(nonpositive[0]) if nonpositive.size > 0 else tail.size

    weights = np.maximum(1.0 - lags[:end] / float(n), 0.0)
    g = 1.0 + 2.0 * float(np.sum(tail[:end] * weights))

    return max(1.0, g)


# =============================================================================
# Statistical Inefficiency Functions
# =============================================================================


def statistical_inefficiency(
    timeseries: ArrayLike,
    mintime: int = DEFAULT_MINTIME,
    fft: bool = True,
) -> float:
    """Compute statistical inefficiency g directly from a timeseries.

    The statistical inefficiency g is the factor by which the variance of
    the sample mean is increased due to correlation:

        Var(mean) = Var(x) * g / N

    This is computed as: g = 1 + 2 * Σ C(t) * (1 - t/N)

    where C(t) is the normalized autocorrelation function and the sum
    includes the finite-size correction factor (1 - t/N) per Chodera et al.
    (2007).

    Parameters
    ----------
    timeseries : array_like
        1D array of values (e.g., contact binary array, RMSD over time)
    mintime : int
        Minimum number of lags to compute before checking for zero crossing.
        Prevents early termination from noise. Default is 3.
    fft : bool
        If True, use FFT-based ACF computation (faster). Default is True.

    Returns
    -------
    float
        Statistical inefficiency g (>= 1.0). The number of effective
        independent samples is N_eff = N / g.

    Examples
    --------
    >>> # Binary contact timeseries
    >>> contacts = np.array([0, 1, 1, 1, 0, 0, 1, 1, ...])
    >>> g = statistical_inefficiency(contacts)
    >>> n_eff = len(contacts) / g
    >>> print(f"Effective samples: {n_eff:.1f}")

    >>> # Continuous observable
    >>> rmsd = np.array([1.2, 1.3, 1.25, 1.4, ...])
    >>> g = statistical_inefficiency(rmsd)

    Notes
    -----
    This implementation follows the algorithm from Chodera et al. (2007)
    J. Chem. Theory Comput. 3:26, as coded in pymbar's MIT-licensed
    `timeseries` module, with the finite-size correction.

    For binary (0/1) data, the algorithm works correctly as the variance
    of a Bernoulli random variable is p(1-p).

    References
    ----------
    Chodera, J. D., Swope, W. C., Pitera, J. W., Seok, C., and Dill, K. A.
        (2007). Journal of Chemical Theory and Computation, 3(1), 26-41.
        doi:10.1021/ct0502864
    Shirts, M. R., and Chodera, J. D. (2008). The Journal of Chemical
        Physics, 129(12), 124105. doi:10.1063/1.2978177
    """
    x = np.asarray(timeseries, dtype=np.float64)
    n = len(x)

    if n < 3:
        logger.warning(f"Timeseries too short ({n} points). Returning g=1.0")
        return 1.0

    # Compute variance
    mu = np.mean(x)
    var = np.var(x)

    if var < 1e-10:
        # Constant timeseries - no correlation
        return 1.0

    # Compute normalized fluctuations
    delta_x = x - mu

    # Compute ACF using FFT for efficiency
    if fft:
        n_fft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
        fft_x = np.fft.fft(delta_x, n_fft)
        acf_unnorm = np.fft.ifft(fft_x * np.conj(fft_x)).real[:n]
        # Normalize by decreasing sample size
        acf = acf_unnorm / (np.arange(n, 0, -1) * var)
    else:
        # Direct computation (slower but clearer)
        acf = np.zeros(n)
        for t in range(n):
            acf[t] = np.mean(delta_x[: n - t] * delta_x[t:]) / var

    return _statistical_inefficiency_from_acf(acf, n, mintime=mintime)


def statistical_inefficiency_multiple(
    timeseries_list: list[ArrayLike],
    mintime: int = DEFAULT_MINTIME,
) -> float:
    """Compute statistical inefficiency from multiple timeseries of different lengths.

    This is critical for aggregating replicates with different frame counts.
    The algorithm computes a global mean μ across all timeseries, then
    averages the ACF numerator and denominator separately before computing g.

    Parameters
    ----------
    timeseries_list : list[ArrayLike]
        List of 1D timeseries arrays (can have different lengths)
    mintime : int
        Minimum number of lags before checking for zero crossing. Default is 3.

    Returns
    -------
    float
        Statistical inefficiency g (>= 1.0)

    Examples
    --------
    >>> # Three replicates with different lengths
    >>> ts1 = np.array([0, 1, 1, 0, 0, 1])  # 6 frames
    >>> ts2 = np.array([1, 1, 0, 0, 0])      # 5 frames
    >>> ts3 = np.array([0, 0, 1, 1, 1, 0, 1])  # 7 frames
    >>> g = statistical_inefficiency_multiple([ts1, ts2, ts3])

    Notes
    -----
    This implementation follows the algorithm from PyMBAR's
    `statistical_inefficiency_multiple()`, adapted without the PyMBAR dependency.

    The algorithm:

    1. Compute global mean μ across all timeseries
    2. For each lag t:
       - Compute sum of (x - μ) products across all timeseries where t < N_k
       - Compute sum of sample counts across all timeseries where t < N_k
       - Average to get C(t)
    3. Sum with finite-size correction

    References
    ----------
    Chodera et al. (2007) J. Chem. Theory Comput. 3:26
    """
    if not timeseries_list:
        return 1.0

    # Convert to numpy arrays
    arrays = [np.asarray(ts, dtype=np.float64) for ts in timeseries_list]
    lengths = np.array([len(a) for a in arrays])
    n_total = int(np.sum(lengths))
    max_length = int(np.max(lengths))

    if n_total < 3:
        logger.warning(f"Total samples too few ({n_total}). Returning g=1.0")
        return 1.0

    # Compute global mean
    total_sum = sum(np.sum(a) for a in arrays)
    mu = total_sum / n_total

    # Compute global variance
    total_var_sum = sum(np.sum((a - mu) ** 2) for a in arrays)
    var = total_var_sum / n_total

    if var < 1e-10:
        return 1.0

    # Compute fluctuations
    deltas = [a - mu for a in arrays]

    # Compute g using averaged ACF
    g = 1.0

    for t in range(1, max_length):
        # Sum ACF contributions from all timeseries where t < N_k
        acf_numerator = 0.0
        acf_denominator = 0.0

        for k, (delta, n_k) in enumerate(zip(deltas, lengths)):
            if t < n_k:
                # This timeseries contributes at lag t
                # Number of pairs at lag t
                n_pairs = n_k - t
                # Sum of products
                product_sum = np.sum(delta[: n_k - t] * delta[t:])
                acf_numerator += product_sum
                acf_denominator += n_pairs

        if acf_denominator < 1:
            # No timeseries has this lag
            break

        # Normalized ACF at lag t
        c_t = acf_numerator / (acf_denominator * var)

        # Check for zero crossing (after mintime)
        if t >= mintime and c_t <= 0:
            break

        # Finite-size correction: use average N across contributing timeseries
        # For simplicity, use the mean length of timeseries that contribute
        contributing = lengths[lengths > t]
        if len(contributing) == 0:
            break
        mean_n = np.mean(contributing)
        weight = 1.0 - float(t) / mean_n

        g += 2.0 * c_t * weight

    # Ensure g >= 1
    g = max(1.0, g)

    return float(g)


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
