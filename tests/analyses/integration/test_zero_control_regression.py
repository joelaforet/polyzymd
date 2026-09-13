"""Regression tests for zero-control percent_change and direction bug.

When the control value is 0.0 and treatment is non-zero, the system should
report ±inf percent_change and the correct direction (not "similar").
"""

from __future__ import annotations

import math

from polyzymd.analyses.base import (
    ComparisonResult,
    ConditionSummary,
    MetricValue,
    PairwiseResult,
)
from polyzymd.analyses.shared.inferential_statistics import percent_change
from polyzymd.analyses.stats import (
    format_pct,
    format_scalar_comparison,
    interpret_direction,
    pairwise_comparisons,
)


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


def test_pairwise_result_inf_round_trip_json() -> None:
    """PairwiseResult should preserve infinite percent_change across JSON round-trip."""
    result = PairwiseResult(
        condition_a="Control",
        condition_b="Treatment",
        metric="m",
        t_statistic=1.0,
        p_value=0.01,
        p_value_adjusted=0.01,
        cohens_d=2.0,
        effect_size_interpretation="large",
        direction="higher",
        significant=True,
        percent_change=math.inf,
    )

    payload = result.model_dump_json()
    loaded = PairwiseResult.model_validate_json(payload)

    assert math.isinf(loaded.percent_change)
    assert loaded.percent_change > 0


def test_comparison_result_inf_round_trip_json() -> None:
    """ComparisonResult should preserve infinite percent_change across JSON round-trip."""
    comparison = ComparisonResult(
        analysis_type="test",
        name="zero-control",
        control_label="Control",
        conditions=[
            ConditionSummary(label="Control", n_replicates=3),
            ConditionSummary(label="Treatment", n_replicates=3),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                metric="metric",
                t_statistic=1.0,
                p_value=0.01,
                p_value_adjusted=0.01,
                cohens_d=2.0,
                effect_size_interpretation="large",
                direction="higher",
                significant=True,
                percent_change=math.inf,
            )
        ],
        ranking=["Treatment", "Control"],
        equilibration_time="0ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )

    payload = comparison.model_dump_json()
    loaded = ComparisonResult.model_validate_json(payload)

    assert math.isinf(loaded.pairwise_comparisons[0].percent_change)
    assert loaded.pairwise_comparisons[0].percent_change > 0


def test_format_scalar_comparison_text_uses_semantic_infinity_label() -> None:
    """Text formatter should show semantic labels for infinite changes."""
    comparison = ComparisonResult(
        analysis_type="test",
        name="format-test",
        control_label="Control",
        conditions=[
            ConditionSummary(label="Control", n_replicates=3, metric_mean=0.0, metric_sem=0.0),
            ConditionSummary(label="Treatment", n_replicates=3, metric_mean=1.0, metric_sem=0.1),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                metric="metric",
                t_statistic=1.0,
                p_value=0.01,
                p_value_adjusted=0.01,
                cohens_d=2.0,
                effect_size_interpretation="large",
                direction="higher",
                significant=True,
                percent_change=math.inf,
            )
        ],
        ranking=["Treatment", "Control"],
        equilibration_time="0ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )

    output = format_scalar_comparison(comparison, output_format="text", metric_key="metric")
    assert "new (baseline=0)" in output
    assert "+∞%" not in output


def test_format_scalar_comparison_markdown_uses_semantic_infinity_label() -> None:
    """Markdown formatter should show semantic labels for infinite changes."""
    comparison = ComparisonResult(
        analysis_type="test",
        name="format-test",
        control_label="Control",
        conditions=[
            ConditionSummary(label="Control", n_replicates=3, metric_mean=0.0, metric_sem=0.0),
            ConditionSummary(label="Treatment", n_replicates=3, metric_mean=1.0, metric_sem=0.1),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                metric="metric",
                t_statistic=1.0,
                p_value=0.01,
                p_value_adjusted=0.01,
                cohens_d=2.0,
                effect_size_interpretation="large",
                direction="higher",
                significant=True,
                percent_change=math.inf,
            )
        ],
        ranking=["Treatment", "Control"],
        equilibration_time="0ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )

    output = format_scalar_comparison(comparison, output_format="markdown", metric_key="metric")
    assert "new (baseline=0)" in output
    assert "+∞%" not in output


def test_pairwise_comparison_zero_control_not_similar() -> None:
    """Pairwise comparison should report +inf and positive direction for zero control."""
    metrics = {
        "Control": MetricValue(
            name="metric",
            mean=0.0,
            sem=0.0,
            replicate_values=[0.0, 0.0, 0.0],
            higher_is_better=True,
            direction_labels=("lower", "similar", "higher"),
        ),
        "Treatment": MetricValue(
            name="metric",
            mean=162.6,
            sem=1.0,
            replicate_values=[160.0, 162.6, 165.0],
            higher_is_better=True,
            direction_labels=("lower", "similar", "higher"),
        ),
    }

    results = pairwise_comparisons(metrics, control_label="Control")

    assert len(results) == 1
    comp = results[0]
    assert math.isinf(comp.percent_change)
    assert comp.percent_change > 0
    assert comp.direction != "similar"
    assert comp.direction == "higher"
