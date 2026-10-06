"""A control value of zero gives a finite difference and a direction.

``ReplicateValues.compare`` reports the plain difference and its direction,
never "similar", when the control is 0.0.
"""

from __future__ import annotations

import math

import pytest


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
