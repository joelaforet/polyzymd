"""Tests for shared sliding-window convergence diagnostics."""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses.shared.convergence import find_convergence_time


def test_find_convergence_time_detects_sustained_plateau() -> None:
    """Convergence should be detected for a trace that plateaus."""
    time_ns = np.arange(0.0, 100.0, 1.0)
    rmsd = np.concatenate(
        [
            np.linspace(0.5, 2.0, 30),
            np.full(70, 2.0),
        ]
    )

    result = find_convergence_time(
        time_ns,
        rmsd,
        window_size_ns=15.0,
        step_size_ns=5.0,
        slope_threshold=0.0005,
        sustained_for_ns=15.0,
    )

    assert result.assessable is True
    assert result.converged is True
    assert result.convergence_time_ns is not None
    assert result.convergence_time_ns >= 30.0


def test_find_convergence_time_non_converged_trace() -> None:
    """Convergence should not be detected for a drifting trace."""
    time_ns = np.arange(0.0, 100.0, 1.0)
    rmsd = np.linspace(0.5, 2.5, 100)

    result = find_convergence_time(time_ns, rmsd)

    assert result.assessable is True
    assert result.converged is False
    assert result.convergence_time_ns is None


def test_find_convergence_time_too_short_is_not_assessable() -> None:
    """Short trajectories should be marked as not assessable."""
    time_ns = np.array([0.0, 2.0, 4.0, 6.0], dtype=np.float64)
    rmsd = np.array([1.0, 1.1, 1.0, 1.05], dtype=np.float64)

    result = find_convergence_time(time_ns, rmsd, window_size_ns=15.0)

    assert result.assessable is False
    assert result.converged is False
    assert result.convergence_time_ns is None


def test_find_convergence_time_invalid_parameters_raise() -> None:
    """Invalid convergence parameters should raise ValueError."""
    time_ns = np.arange(0.0, 20.0, 1.0)
    rmsd = np.ones_like(time_ns)

    with pytest.raises(ValueError, match="step_size_ns"):
        find_convergence_time(time_ns, rmsd, window_size_ns=5.0, step_size_ns=10.0)

    with pytest.raises(ValueError, match="slope_threshold"):
        find_convergence_time(time_ns, rmsd, slope_threshold=0.0)


def test_find_convergence_time_edge_cases() -> None:
    """Constant traces can converge, single-point traces are not assessable."""
    constant_time = np.arange(0.0, 60.0, 1.0)
    constant_rmsd = np.full_like(constant_time, 1.25)
    constant_result = find_convergence_time(constant_time, constant_rmsd)
    assert constant_result.assessable is True
    assert constant_result.converged is True

    single_time = np.array([0.0], dtype=np.float64)
    single_rmsd = np.array([1.0], dtype=np.float64)
    single_result = find_convergence_time(single_time, single_rmsd)
    assert single_result.assessable is False
    assert single_result.converged is False


def test_nan_in_values_raises() -> None:
    """NaN in values should raise ValueError."""
    time_ns = np.arange(0.0, 20.0, 1.0)
    values = np.ones_like(time_ns)
    values[3] = np.nan

    with pytest.raises(ValueError, match="must contain only finite values"):
        find_convergence_time(time_ns, values)


def test_inf_in_values_raises() -> None:
    """Inf in values should raise ValueError."""
    time_ns = np.arange(0.0, 20.0, 1.0)
    values = np.ones_like(time_ns)
    values[5] = np.inf

    with pytest.raises(ValueError, match="must contain only finite values"):
        find_convergence_time(time_ns, values)


def test_nan_in_time_raises() -> None:
    """NaN in time_ns should raise ValueError."""
    time_ns = np.arange(0.0, 20.0, 1.0)
    time_ns[7] = np.nan
    values = np.ones_like(time_ns)

    with pytest.raises(ValueError, match="must contain only finite values"):
        find_convergence_time(time_ns, values)


def test_nan_parameter_raises() -> None:
    """NaN scalar parameters should raise ValueError."""
    time_ns = np.arange(0.0, 20.0, 1.0)
    values = np.ones_like(time_ns)

    with pytest.raises(ValueError, match="window_size_ns must be finite"):
        find_convergence_time(time_ns, values, window_size_ns=np.nan)
