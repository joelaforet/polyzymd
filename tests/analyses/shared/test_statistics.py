"""Tests for shared statistics helpers."""

from __future__ import annotations

import pytest

from polyzymd.analyses.shared.statistics import weighted_mean_with_sem


def test_weighted_mean_with_sem_rejects_zero_sem_in_inverse_variance_mode() -> None:
    """Inverse-variance weighting should reject zero SEM values."""
    with pytest.raises(ValueError, match="greater than zero"):
        weighted_mean_with_sem(means=[1.0, 2.0], sems=[0.0, 0.1])


def test_weighted_mean_with_sem_rejects_nonfinite_sem_in_inverse_variance_mode() -> None:
    """Inverse-variance weighting should reject nonfinite SEM values."""
    with pytest.raises(ValueError, match="finite SEM"):
        weighted_mean_with_sem(means=[1.0, 2.0], sems=[float("nan"), 0.1])


def test_weighted_mean_with_sem_allows_zero_sem_with_explicit_weights() -> None:
    """Explicit weighting mode should not enforce inverse-variance SEM constraints."""
    result = weighted_mean_with_sem(means=[1.0, 2.0], sems=[0.0, 0.1], weights=[1.0, 1.0])
    assert result.mean == pytest.approx(1.5)
