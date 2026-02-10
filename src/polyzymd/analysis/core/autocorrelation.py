"""Autocorrelation analysis for independent sampling.

MD trajectories are highly correlated in time - consecutive frames are not
independent samples. This module provides tools to:

1. Compute the autocorrelation function (ACF) of an observable
2. Estimate the correlation time (τ) from the ACF
3. Select independent frames based on τ for proper statistics

Key Concepts
------------
- **Autocorrelation function (ACF)**: Measures how correlated a signal is with
  itself at different time lags. ACF(0) = 1, and ACF decays toward 0.

- **Correlation time (τ)**: Characteristic time for decorrelation. Frames
  separated by > 2τ are approximately independent.

- **Independent samples**: For proper SEM calculation, we need N_eff independent
  samples, not N_frames correlated observations.

Methods for τ estimation
------------------------
- **First zero crossing**: τ is lag where ACF first crosses zero
- **Exponential fit**: Fit ACF = exp(-t/τ) and extract τ
- **Integration**: τ = ∫ACF(t)dt from 0 to first zero (or cutoff)

References
----------
- Flyvbjerg & Petersen (1989) J. Chem. Phys. 91:461 (block averaging)
- Chodera et al. (2007) J. Chem. Theory Comput. 3:26 (statistical inefficiency)
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray


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
    """

    tau: float
    tau_unit: str
    method: str
    n_independent: int
    statistical_inefficiency: float

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "tau": self.tau,
            "tau_unit": self.tau_unit,
            "method": self.method,
            "n_independent": self.n_independent,
            "statistical_inefficiency": self.statistical_inefficiency,
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
        tau = _integrate_acf(lags, acf, zero_crossing_idx, dt)

    else:
        raise ValueError(f"Unknown method: {method}")

    # Compute statistical inefficiency: g = 1 + 2*τ/dt
    # This is the factor by which variance is inflated due to correlation
    g = 1.0 + 2.0 * tau / dt

    # Number of independent samples
    n_independent = max(1, int(n_frames / g))

    return CorrelationTimeResult(
        tau=tau,
        tau_unit=unit,
        method=method,
        n_independent=n_independent,
        statistical_inefficiency=g,
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
) -> float:
    """Estimate τ by integrating ACF."""
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

    # Trapezoidal integration
    tau = float(np.trapezoid(acf[:int_end], lags[:int_end]))

    return max(tau, dt)  # At least one timestep
