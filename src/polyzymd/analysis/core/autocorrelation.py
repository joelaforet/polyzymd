"""Autocorrelation analysis for independent sampling.

MD trajectories are highly correlated in time - consecutive frames are not
independent samples. This module provides tools to:

1. Compute the autocorrelation function (ACF) of an observable
2. Estimate the correlation time (τ) from the ACF
3. Compute statistical inefficiency (g) for proper uncertainty quantification
4. Select independent frames based on τ for proper statistics

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

Methods for τ estimation
------------------------
- **First zero crossing**: τ is lag where ACF first crosses zero
- **Exponential fit**: Fit ACF = exp(-t/τ) and extract τ
- **Integration**: τ = ∫ACF(t)dt from 0 to first zero (or cutoff)

Statistical Validity
--------------------
The number of effective independent samples (N_eff) is computed as:
    N_eff = N / g = N / (1 + 2*Σ C(t)*(1-t/N))

This matches the algorithm from Chodera et al. (2007) with the finite-size
correction factor (1-t/N). When N_eff < 10, statistical estimates (mean, SEM)
may be unreliable, and users should be warned per LiveCoMS best practices
(Grossfield et al., 2018).

For multiple timeseries of different lengths (e.g., replicates), use
`statistical_inefficiency_multiple()` which correctly handles the averaging.

References
----------
- Flyvbjerg & Petersen (1989) J. Chem. Phys. 91:461 (block averaging)
- Chodera et al. (2007) J. Chem. Theory Comput. 3:26 (statistical inefficiency)
- Grossfield et al. (2018) LiveCoMS 1:5067 (uncertainty quantification)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

logger = logging.getLogger(__name__)

# Minimum recommended independent samples for reliable statistics
MIN_RECOMMENDED_N_INDEPENDENT = 10


class CorrelationTimeMethod(str, Enum):
    """Method for estimating correlation time from ACF."""

    FIRST_ZERO = "first_zero"
    EXPONENTIAL_FIT = "exponential_fit"
    INTEGRATION = "integration"


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
    """

    lags: NDArray[np.float64]
    acf: NDArray[np.float64]
    timestep: float
    timestep_unit: str

    def __len__(self) -> int:
        return len(self.lags)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "lags": self.lags.tolist(),
            "acf": self.acf.tolist(),
            "timestep": self.timestep,
            "timestep_unit": self.timestep_unit,
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
    n_independent : int
        Estimated number of independent samples in trajectory
    statistical_inefficiency : float
        g = 1 + 2*tau/dt, factor by which variance is inflated
    warning : str | None
        Warning message if statistics may be unreliable (e.g., N_ind < 10)
    """

    tau: float
    tau_unit: str
    method: str
    n_independent: int
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

    # FFT-based autocorrelation (much faster than direct computation)
    # Pad to next power of 2 for FFT efficiency
    n_fft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
    fft_x = np.fft.fft(x_centered, n_fft)
    acf_full = np.fft.ifft(fft_x * np.conj(fft_x)).real[:n]

    # Normalize by decreasing sample size and variance
    acf_full = acf_full / (np.arange(n, 0, -1) * np.var(x_centered))

    # Take only up to max_lag
    acf = acf_full[: max_lag + 1]
    lags = np.arange(max_lag + 1) * timestep

    return ACFResult(
        lags=lags,
        acf=acf,
        timestep=timestep,
        timestep_unit=timestep_unit,
    )


def estimate_correlation_time(
    acf_or_timeseries: ACFResult | ArrayLike,
    timestep: float = 1.0,
    timestep_unit: str = "frames",
    method: Literal["first_zero", "exponential_fit", "integration"] = "integration",
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
    method : {"first_zero", "exponential_fit", "integration"}
        Method for estimating τ:
        - "first_zero": Lag where ACF first crosses zero
        - "exponential_fit": Fit ACF = exp(-t/τ)
        - "integration": τ = ∫ACF(t)dt (recommended, most robust)
    n_frames : int, optional
        Total number of frames (for computing n_independent).
        Only needed if passing ACFResult.

    Returns
    -------
    CorrelationTimeResult
        Contains tau, method used, n_independent, statistical_inefficiency

    Examples
    --------
    >>> acf_result = compute_acf(rmsd, timestep=10.0, timestep_unit="ps")
    >>> tau_result = estimate_correlation_time(acf_result, method="integration")
    >>> print(f"Correlation time: {tau_result.tau:.1f} {tau_result.tau_unit}")
    >>> print(f"Independent samples: {tau_result.n_independent}")

    Notes
    -----
    The "integration" method is most robust for noisy ACFs. It computes:
        τ = ∫₀^∞ ACF(t) dt ≈ Σ ACF[i] * dt

    Integration stops at first zero crossing to avoid noise contribution.
    """
    # Handle input type
    if isinstance(acf_or_timeseries, ACFResult):
        acf_result = acf_or_timeseries
        dt = acf_result.timestep
        unit = acf_result.timestep_unit
        acf = acf_result.acf
        lags = acf_result.lags
        if n_frames is None:
            # Estimate from ACF length (we computed up to n//4)
            n_frames = len(acf) * 4
    else:
        # Compute ACF from raw timeseries
        acf_result = compute_acf(
            acf_or_timeseries,
            timestep=timestep,
            timestep_unit=timestep_unit,
        )
        acf = acf_result.acf
        lags = acf_result.lags
        dt = timestep
        unit = timestep_unit
        n_frames = len(np.asarray(acf_or_timeseries))

    # Find first zero crossing
    zero_crossing_idx = _find_first_zero_crossing(acf)

    if method == "first_zero":
        if zero_crossing_idx is None:
            # ACF never crosses zero, use full length
            tau = lags[-1]
        else:
            tau = float(lags[zero_crossing_idx])

    elif method == "exponential_fit":
        tau = _fit_exponential_acf(lags, acf, zero_crossing_idx)

    elif method == "integration":
        tau = _integrate_acf(lags, acf, zero_crossing_idx, dt, n_frames=n_frames)

    else:
        raise ValueError(f"Unknown method: {method}")

    # Compute statistical inefficiency: g = 1 + 2*τ/dt
    # This is the factor by which variance is inflated due to correlation
    g = 1.0 + 2.0 * tau / dt

    # Number of independent samples
    n_independent = max(1, int(n_frames / g))

    # Generate warning if statistics may be unreliable
    warning = None
    if n_independent < MIN_RECOMMENDED_N_INDEPENDENT:
        warning = (
            f"Low statistical reliability: only {n_independent} independent samples "
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


def get_independent_indices(
    n_frames: int,
    correlation_time: float,
    timestep: float = 1.0,
    start_frame: int = 0,
) -> NDArray[np.int64]:
    """Get frame indices for independent samples.

    Selects frames separated by at least 2*τ (correlation time) to
    ensure approximate independence for statistical analysis.

    Parameters
    ----------
    n_frames : int
        Total number of frames in trajectory
    correlation_time : float
        Correlation time τ (in same units as timestep)
    timestep : float, optional
        Time between frames. Default is 1.0.
    start_frame : int, optional
        First frame to consider (after equilibration). Default is 0.
        Note: Frame indices are 0-indexed internally, but user-facing
        documentation uses 1-indexed (PyMOL convention).

    Returns
    -------
    NDArray[np.int64]
        Array of frame indices (0-indexed) that are approximately independent

    Examples
    --------
    >>> # Get independent frames for RMSF calculation
    >>> tau_result = estimate_correlation_time(rmsd, timestep=10.0)
    >>> indices = get_independent_indices(
    ...     n_frames=10000,
    ...     correlation_time=tau_result.tau,
    ...     timestep=10.0,
    ...     start_frame=1000,  # Skip first 1000 frames for equilibration
    ... )
    >>> print(f"Using {len(indices)} independent frames")

    Notes
    -----
    Frame indices returned are 0-indexed (for direct use with MDAnalysis).
    When displaying to users, add 1 for PyMOL convention.

    The spacing is set to 2*τ/timestep, which gives frames with
    negligible correlation (ACF < 0.05 for exponential decay).
    """
    if n_frames <= start_frame:
        raise ValueError(f"start_frame ({start_frame}) >= n_frames ({n_frames})")

    # Convert correlation time to frame spacing
    # Use 2*τ for good independence (ACF ≈ exp(-2) ≈ 0.14)
    frame_spacing = max(1, int(np.ceil(2.0 * correlation_time / timestep)))

    # Generate indices
    indices = np.arange(start_frame, n_frames, frame_spacing, dtype=np.int64)

    return indices


def _find_first_zero_crossing(acf: NDArray[np.float64]) -> int | None:
    """Find index of first zero crossing in ACF."""
    sign_changes = np.where(np.diff(np.sign(acf)))[0]
    if len(sign_changes) > 0:
        return int(sign_changes[0])
    return None


def _fit_exponential_acf(
    lags: NDArray[np.float64],
    acf: NDArray[np.float64],
    zero_crossing_idx: int | None,
) -> float:
    """Fit exponential decay to ACF and extract τ."""
    # Only fit up to first zero crossing (or halfway)
    if zero_crossing_idx is not None:
        fit_end = zero_crossing_idx
    else:
        fit_end = len(acf) // 2

    fit_end = max(3, fit_end)  # Need at least 3 points

    # ACF = exp(-t/τ) => log(ACF) = -t/τ
    # Linear fit: y = -x/τ where y = log(ACF), x = t
    acf_positive = np.maximum(acf[:fit_end], 1e-10)  # Avoid log(0)
    log_acf = np.log(acf_positive)

    # Linear regression
    try:
        slope, _ = np.polyfit(lags[:fit_end], log_acf, 1)
        tau = -1.0 / slope if slope < 0 else lags[fit_end]
    except (np.linalg.LinAlgError, ValueError):
        # Fallback to integration if fit fails
        tau = lags[fit_end]

    return float(max(tau, lags[1]))  # At least one timestep


def _integrate_acf(
    lags: NDArray[np.float64],
    acf: NDArray[np.float64],
    zero_crossing_idx: int | None,
    dt: float,
    n_frames: int | None = None,
    use_finite_size_correction: bool = True,
) -> float:
    """Estimate τ by integrating ACF.

    Parameters
    ----------
    lags : NDArray[np.float64]
        Time lags
    acf : NDArray[np.float64]
        Autocorrelation values
    zero_crossing_idx : int | None
        Index of first zero crossing
    dt : float
        Timestep
    n_frames : int | None
        Total number of frames (for finite-size correction)
    use_finite_size_correction : bool
        If True, apply (1-t/N) weighting per Chodera et al. 2007

    Returns
    -------
    float
        Estimated correlation time τ
    """
    # Integrate up to first zero crossing
    if zero_crossing_idx is not None:
        int_end = zero_crossing_idx + 1
    else:
        # Find where ACF drops below threshold
        below_threshold = np.where(acf < 0.05)[0]
        if len(below_threshold) > 0:
            int_end = below_threshold[0] + 1
        else:
            int_end = len(acf)

    # Apply finite-size correction if requested
    if use_finite_size_correction and n_frames is not None and n_frames > 0:
        # Weight ACF by (1 - t/N) per Chodera et al. 2007
        # This accounts for reduced sample size at longer lags
        lag_indices = np.arange(int_end)
        weights = 1.0 - lag_indices / n_frames
        weighted_acf = acf[:int_end] * weights
    else:
        weighted_acf = acf[:int_end]

    # Trapezoidal integration (trapezoid in numpy 2.0+, trapz in older versions)
    trapz_func = getattr(np, "trapezoid", np.trapz)
    tau = float(trapz_func(weighted_acf, lags[:int_end]))

    return max(tau, dt)  # At least one timestep


# =============================================================================
# Statistical Inefficiency Functions
# =============================================================================


def statistical_inefficiency(
    timeseries: ArrayLike,
    mintime: int = 3,
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
    J. Chem. Theory Comput. 3:26, with the finite-size correction.

    For binary (0/1) data, the algorithm works correctly as the variance
    of a Bernoulli random variable is p(1-p).

    References
    ----------
    Chodera et al. (2007) J. Chem. Theory Comput. 3:26
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

    # Compute g = 1 + 2 * sum(C(t) * (1 - t/N))
    # Start with g = 1 (for lag 0, C(0) = 1, but we don't count it in the sum)
    g = 1.0

    # Sum over positive lags with finite-size correction
    for t in range(1, n):
        # Finite-size correction factor
        weight = 1.0 - float(t) / n

        # Check for zero crossing (after mintime)
        if t >= mintime and acf[t] <= 0:
            break

        g += 2.0 * acf[t] * weight

    # Ensure g >= 1
    g = max(1.0, g)

    return float(g)


def statistical_inefficiency_multiple(
    timeseries_list: list[ArrayLike],
    mintime: int = 3,
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
