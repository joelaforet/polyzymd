"""A control value of zero gives an infinite percent change and a finite difference.

``percent_change`` returns plus or minus infinity when the control is 0.0 and
the treatment is not, and ``ReplicateValues.compare`` reports the plain
difference and its direction, never "similar".
"""

from __future__ import annotations

import math

import pytest

from polyzymd.analyses.shared.inferential_statistics import percent_change


def test_percent_change_zero_to_zero_is_zero() -> None:
    """Percent change should be zero when both values are zero."""
    assert percent_change(0.0, 0.0) == 0.0


def test_percent_change_zero_to_positive_is_inf() -> None:
    """Percent change should be +inf when control is zero and treatment is positive."""
    assert math.isinf(percent_change(0.0, 162.6))
    assert percent_change(0.0, 162.6) > 0


def test_percent_change_zero_to_negative_is_minus_inf() -> None:
    """Percent change should be -inf when control is zero and treatment is negative."""
    assert percent_change(0.0, -5.0) == -math.inf


def test_percent_change_near_zero_stays_finite() -> None:
    """Non-zero controls should use the standard finite formula."""
    pct = percent_change(1e-12, 2e-12)
    assert math.isfinite(pct)
    assert pct == 100.0


def test_percent_change_nan_input_returns_nan() -> None:
    """NaN inputs should propagate as NaN."""
    assert math.isnan(percent_change(math.nan, 5.0))


def test_function_path_zero_control_reports_a_finite_difference_and_direction() -> None:
    """A control of zero gives the plain difference and a larger verdict, never "similar".

    The function path reports ``mean(b) - mean(a)`` instead of a percent
    change, so a zero control cannot produce an infinite value.
    """
    from tests._support.analysis_testkit import replicate_values

    values = replicate_values({"Control": [0.0, 0.0, 0.0], "Treatment": [3.8, 4.0, 4.2]})
    report = values.compare()
    (row,) = report.pairwise

    assert math.isfinite(row.delta)
    assert row.delta == pytest.approx(4.0)
    assert row.testable and row.significant
    assert row.direction == "increased"
    assert report.verdict[0].startswith("Treatment larger mean_rg than Control")
