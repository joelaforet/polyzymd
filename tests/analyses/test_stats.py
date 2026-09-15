"""Tests for the direction wording helpers in ``polyzymd.analyses.stats``."""

from __future__ import annotations

import math

import pytest

from polyzymd.analyses.stats import format_pct, interpret_direction


@pytest.mark.parametrize(
    ("change", "expected"),
    [(5.0, "increased"), (-5.0, "decreased"), (0.4, "unchanged"), (-0.4, "unchanged")],
)
def test_interpret_direction_uses_the_one_percent_threshold(change: float, expected: str) -> None:
    """A change under one percent reads as unchanged, either way."""
    assert interpret_direction(change) == expected


def test_interpret_direction_takes_custom_labels_in_order() -> None:
    """The tuple is ordered negative, unchanged, positive."""
    labels = ("lower", "similar", "higher")
    assert interpret_direction(12.0, labels) == "higher"
    assert interpret_direction(-12.0, labels) == "lower"
    assert interpret_direction(0.0, labels) == "similar"


def test_interpret_direction_honours_a_wider_threshold() -> None:
    """A caller that wants a wider dead band says so."""
    assert interpret_direction(3.0, threshold=5.0) == "unchanged"
    assert interpret_direction(7.0, threshold=5.0) == "increased"


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (math.inf, "new (baseline=0)"),
        (-math.inf, "gone (current=0)"),
        (math.nan, "undefined"),
        (12.34, "+12.3%"),
        (-12.34, "-12.3%"),
        (-0.0, "+0.0%"),
    ],
)
def test_format_pct_names_the_degenerate_cases(value: float, expected: str) -> None:
    """A zero baseline gets words rather than an infinity symbol."""
    assert format_pct(value) == expected
