"""Regression tests for zero-control percent change and direction wording.

When the control value is 0.0 and the treatment is not, the reported percent
change is plus or minus infinity and the direction word follows its sign rather
than falling back to "similar".
"""

from __future__ import annotations

import math

from polyzymd.analyses.shared.inferential_statistics import percent_change
from polyzymd.analyses.stats import format_pct, interpret_direction


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


def test_interpret_direction_positive_infinity() -> None:
    """Positive infinity should map to the positive direction label."""
    labels = ("lower", "similar", "higher")
    assert interpret_direction(math.inf, labels) == "higher"


def test_interpret_direction_negative_infinity() -> None:
    """Negative infinity should map to the negative direction label."""
    labels = ("lower", "similar", "higher")
    assert interpret_direction(-math.inf, labels) == "lower"


def test_interpret_direction_nan_maps_to_unchanged() -> None:
    """NaN percent changes should map to the unchanged label."""
    labels = ("lower", "similar", "higher")
    assert interpret_direction(math.nan, labels) == "similar"


def test_format_pct_inf_nan_and_finite() -> None:
    """Percent formatter should render infinity and NaN safely."""
    assert format_pct(math.inf) == "new (baseline=0)"
    assert format_pct(-math.inf) == "gone (current=0)"
    assert format_pct(math.nan) == "undefined"
    assert format_pct(12.3) == "+12.3%"


def test_format_pct_normalizes_negative_zero() -> None:
    """Formatter should normalize negative zero to +0.0%."""
    assert format_pct(-0.0) == "+0.0%"
