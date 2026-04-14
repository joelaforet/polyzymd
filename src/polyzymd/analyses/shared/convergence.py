"""Convergence diagnostics for sliding-window timeseries analysis.

This module implements a sliding-window convergence heuristic adapted from a
collaborator notebook used for RMSD equilibration checks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike


@dataclass
class ConvergenceResult:
    """Container for convergence diagnostics.

    Attributes
    ----------
    converged : bool
        Whether sustained convergence was detected.
    assessable : bool
        Whether convergence could be assessed from available data.
    convergence_time_ns : float | None
        Start time of the first sustained converged period.
    window_start_times_ns : list[float]
        Start times for each sliding window.
    window_mean_values : list[float]
        Mean signal value in each sliding window.
    slope_times_ns : list[float]
        Time points associated with slope estimates.
    slopes : list[float]
        Slopes between successive window means.
    window_size_ns : float
        Sliding window width in ns.
    step_size_ns : float
        Sliding window stride in ns.
    slope_threshold : float
        Absolute slope cutoff used for convergence.
    sustained_for_ns : float
        Required sustained duration below slope threshold.
    message : str
        Human-readable status message.
    """

    converged: bool
    assessable: bool
    convergence_time_ns: float | None
    window_start_times_ns: list[float]
    window_mean_values: list[float]
    slope_times_ns: list[float]
    slopes: list[float]
    window_size_ns: float
    step_size_ns: float
    slope_threshold: float
    sustained_for_ns: float
    message: str


def find_convergence_time(
    time_ns: ArrayLike,
    values: ArrayLike,
    window_size_ns: float = 15.0,
    step_size_ns: float = 5.0,
    slope_threshold: float = 0.0005,
    sustained_for_ns: float = 15.0,
) -> ConvergenceResult:
    """Find sustained convergence time using a sliding-window slope heuristic.

    Parameters
    ----------
    time_ns : array_like
        Monotonically increasing time values in ns.
    values : array_like
        Signal values sampled at ``time_ns``.
    window_size_ns : float, optional
        Width of each averaging window in ns.
    step_size_ns : float, optional
        Sliding step between successive windows in ns.
    slope_threshold : float, optional
        Absolute slope threshold for classifying a window-to-window change as
        converged.
    sustained_for_ns : float, optional
        Required cumulative duration below threshold before declaring
        convergence.

    Returns
    -------
    ConvergenceResult
        Full convergence diagnostics, including intermediate window means and
        slope traces.

    Raises
    ------
    ValueError
        Raised when inputs are invalid.
    """
    _validate_parameters(
        window_size_ns=window_size_ns,
        step_size_ns=step_size_ns,
        slope_threshold=slope_threshold,
        sustained_for_ns=sustained_for_ns,
    )

    time_arr = np.asarray(time_ns, dtype=np.float64)
    value_arr = np.asarray(values, dtype=np.float64)

    if not np.all(np.isfinite(time_arr)) or not np.all(np.isfinite(value_arr)):
        raise ValueError("time_ns and values must contain only finite values (no NaN or inf)")

    if time_arr.ndim != 1 or value_arr.ndim != 1:
        raise ValueError("time_ns and values must be 1D arrays")
    if len(time_arr) != len(value_arr):
        raise ValueError("time_ns and values must have the same length")
    if len(time_arr) == 0:
        return ConvergenceResult(
            converged=False,
            assessable=False,
            convergence_time_ns=None,
            window_start_times_ns=[],
            window_mean_values=[],
            slope_times_ns=[],
            slopes=[],
            window_size_ns=window_size_ns,
            step_size_ns=step_size_ns,
            slope_threshold=slope_threshold,
            sustained_for_ns=sustained_for_ns,
            message="Convergence not assessable: empty timeseries",
        )

    if np.any(np.diff(time_arr) <= 0):
        raise ValueError("time_ns must be strictly increasing")

    total_duration = float(time_arr[-1] - time_arr[0])
    if total_duration < window_size_ns:
        return ConvergenceResult(
            converged=False,
            assessable=False,
            convergence_time_ns=None,
            window_start_times_ns=[],
            window_mean_values=[],
            slope_times_ns=[],
            slopes=[],
            window_size_ns=window_size_ns,
            step_size_ns=step_size_ns,
            slope_threshold=slope_threshold,
            sustained_for_ns=sustained_for_ns,
            message=(
                "Convergence not assessable: trajectory shorter than window size "
                f"({total_duration:.3f} ns < {window_size_ns:.3f} ns)"
            ),
        )

    window_means: list[float] = []
    window_start_times: list[float] = []
    start_time = float(time_arr[0])

    while start_time + window_size_ns <= float(time_arr[-1]):
        end_time = start_time + window_size_ns
        indices = np.where((time_arr >= start_time) & (time_arr < end_time))[0]
        if len(indices) > 0:
            window_means.append(float(np.mean(value_arr[indices])))
            window_start_times.append(start_time)
        start_time += step_size_ns

    if len(window_means) < 2:
        return ConvergenceResult(
            converged=False,
            assessable=False,
            convergence_time_ns=None,
            window_start_times_ns=window_start_times,
            window_mean_values=window_means,
            slope_times_ns=[],
            slopes=[],
            window_size_ns=window_size_ns,
            step_size_ns=step_size_ns,
            slope_threshold=slope_threshold,
            sustained_for_ns=sustained_for_ns,
            message="Convergence not assessable: insufficient sliding windows",
        )

    sustained_duration = 0.0
    convergence_start_time: float | None = None
    slope_times: list[float] = []
    slopes: list[float] = []

    for i in range(1, len(window_means)):
        delta_rmsd = window_means[i] - window_means[i - 1]
        delta_time = window_start_times[i] - window_start_times[i - 1]
        slope = delta_rmsd / delta_time
        slopes.append(float(slope))
        slope_times.append(float(window_start_times[i]))

        if abs(slope) < slope_threshold:
            if sustained_duration == 0:
                convergence_start_time = window_start_times[i - 1]
            sustained_duration += delta_time
        else:
            sustained_duration = 0.0
            convergence_start_time = None

        if sustained_duration >= sustained_for_ns:
            return ConvergenceResult(
                converged=True,
                assessable=True,
                convergence_time_ns=convergence_start_time,
                window_start_times_ns=window_start_times,
                window_mean_values=window_means,
                slope_times_ns=slope_times,
                slopes=slopes,
                window_size_ns=window_size_ns,
                step_size_ns=step_size_ns,
                slope_threshold=slope_threshold,
                sustained_for_ns=sustained_for_ns,
                message=f"Converged at {convergence_start_time:.3f} ns",
            )

    return ConvergenceResult(
        converged=False,
        assessable=True,
        convergence_time_ns=None,
        window_start_times_ns=window_start_times,
        window_mean_values=window_means,
        slope_times_ns=slope_times,
        slopes=slopes,
        window_size_ns=window_size_ns,
        step_size_ns=step_size_ns,
        slope_threshold=slope_threshold,
        sustained_for_ns=sustained_for_ns,
        message="Did not converge within available trajectory window",
    )


def _validate_parameters(
    *,
    window_size_ns: float,
    step_size_ns: float,
    slope_threshold: float,
    sustained_for_ns: float,
) -> None:
    """Validate sliding-window convergence parameters.

    Parameters
    ----------
    window_size_ns : float
        Window width in ns.
    step_size_ns : float
        Window stride in ns.
    slope_threshold : float
        Slope magnitude threshold.
    sustained_for_ns : float
        Required sustained converged duration.
    """
    scalar_params = {
        "window_size_ns": window_size_ns,
        "step_size_ns": step_size_ns,
        "slope_threshold": slope_threshold,
        "sustained_for_ns": sustained_for_ns,
    }
    for param_name, param in scalar_params.items():
        if not math.isfinite(param):
            raise ValueError(f"{param_name} must be finite, got {param}")

    if window_size_ns <= 0:
        raise ValueError("window_size_ns must be > 0")
    if step_size_ns <= 0:
        raise ValueError("step_size_ns must be > 0")
    if slope_threshold <= 0:
        raise ValueError("slope_threshold must be > 0")
    if sustained_for_ns <= 0:
        raise ValueError("sustained_for_ns must be > 0")
    if step_size_ns > window_size_ns:
        raise ValueError("step_size_ns must be <= window_size_ns")
