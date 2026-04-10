"""Tests for Benjamini-Hochberg FDR correction in inferential statistics."""

from __future__ import annotations

import math

import pytest

from polyzymd.analyses.shared.inferential_statistics import benjamini_hochberg


def test_bh_all_significant() -> None:
    """All low p-values should remain significant after BH correction."""
    p_values = [0.001, 0.003, 0.01]
    results = benjamini_hochberg(p_values, alpha=0.05)

    assert len(results) == 3
    assert all(r.adjusted_p_value is not None for r in results)
    assert all(r.significant for r in results)


def test_bh_mixed_significance() -> None:
    """Mixed p-values should produce mixed significance after BH correction."""
    p_values = [0.001, 0.04, 0.03, 0.50]
    results = benjamini_hochberg(p_values, alpha=0.05)

    adjusted = [r.adjusted_p_value for r in results]
    assert adjusted == pytest.approx([0.004, 0.0533333333, 0.0533333333, 0.5], rel=1e-7)
    assert [r.significant for r in results] == [True, False, False, False]


def test_bh_all_non_significant() -> None:
    """High p-values should remain non-significant after BH correction."""
    p_values = [0.51, 0.72, 0.93]
    results = benjamini_hochberg(p_values, alpha=0.05)

    assert all(not r.significant for r in results)


def test_bh_single_value() -> None:
    """Single-value BH family should return adjusted equal to raw p-value."""
    p_values = [0.023]
    results = benjamini_hochberg(p_values, alpha=0.05)

    assert len(results) == 1
    assert results[0].adjusted_p_value == pytest.approx(0.023)
    assert results[0].rank == 1


def test_bh_none_handling() -> None:
    """None p-values should pass through with None adjusted values and no significance."""
    p_values = [None, 0.01, None, 0.20]
    results = benjamini_hochberg(p_values, alpha=0.05)

    assert results[0].raw_p_value is None
    assert results[0].adjusted_p_value is None
    assert not results[0].significant
    assert results[0].rank is None

    assert results[2].raw_p_value is None
    assert results[2].adjusted_p_value is None
    assert not results[2].significant
    assert results[2].rank is None

    assert results[1].rank == 1
    assert results[3].rank == 2


def test_bh_all_none() -> None:
    """All-None input should return all-None BH results."""
    p_values = [None, None, None]
    results = benjamini_hochberg(p_values, alpha=0.05)

    assert len(results) == 3
    assert all(r.raw_p_value is None for r in results)
    assert all(r.adjusted_p_value is None for r in results)
    assert all(not r.significant for r in results)
    assert all(r.rank is None for r in results)


def test_nan_treated_as_missing() -> None:
    """NaN p-values should pass through as missing and not affect others."""
    p_values = [0.01, math.nan, 0.04]
    results = benjamini_hochberg(p_values, alpha=0.05)

    assert results[1].raw_p_value is None
    assert results[1].adjusted_p_value is None
    assert results[1].rank is None
    assert not results[1].significant

    assert results[0].adjusted_p_value == pytest.approx(0.02)
    assert results[2].adjusted_p_value == pytest.approx(0.04)
    assert results[0].rank == 1
    assert results[2].rank == 2


def test_mixed_nan_and_none() -> None:
    """NaN and None should both be treated as missing values."""
    p_values = [None, math.nan, 0.03, 0.20]
    results = benjamini_hochberg(p_values, alpha=0.05)

    assert results[0].raw_p_value is None
    assert results[0].adjusted_p_value is None
    assert results[0].rank is None
    assert not results[0].significant

    assert results[1].raw_p_value is None
    assert results[1].adjusted_p_value is None
    assert results[1].rank is None
    assert not results[1].significant

    assert results[2].adjusted_p_value == pytest.approx(0.06)
    assert results[3].adjusted_p_value == pytest.approx(0.20)


def test_bh_empty_input() -> None:
    """Empty p-value family should return an empty list."""
    assert benjamini_hochberg([], alpha=0.05) == []


def test_bh_monotonicity_in_rank_order() -> None:
    """Adjusted p-values should be non-decreasing in BH rank order."""
    p_values = [0.07, 0.001, 0.8, 0.02, 0.04]
    results = benjamini_hochberg(p_values, alpha=0.05)

    ranked = sorted(
        [r for r in results if r.rank is not None],
        key=lambda r: r.rank if r.rank is not None else 0,
    )
    adjusted = [r.adjusted_p_value for r in ranked]

    assert all(a is not None for a in adjusted)
    assert adjusted == sorted(adjusted)


def test_bh_known_textbook_example() -> None:
    """Classic BH example should mark first four hypotheses significant."""
    p_values = [0.001, 0.008, 0.039, 0.041, 0.23, 0.30, 0.76]
    results = benjamini_hochberg(p_values, alpha=0.1)

    significant_flags = [r.significant for r in results]
    assert significant_flags == [True, True, True, True, False, False, False]


def test_bh_alpha_one_marks_all_significant() -> None:
    """Alpha of 1.0 should mark every non-None hypothesis significant."""
    p_values = [0.2, None, 0.9, 0.01]
    results = benjamini_hochberg(p_values, alpha=1.0)

    assert results[0].significant
    assert not results[1].significant
    assert results[2].significant
    assert results[3].significant


@pytest.mark.parametrize("alpha", [0.0, -0.1, 1.1])
def test_bh_alpha_must_be_between_zero_and_one(alpha: float) -> None:
    """Alpha outside valid range should raise a ValueError."""
    with pytest.raises(ValueError, match="alpha"):
        benjamini_hochberg([0.01, 0.02], alpha=alpha)
