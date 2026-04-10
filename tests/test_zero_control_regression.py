"""Regression tests for zero-control percent_change and direction bug.

When the control value is 0.0 and treatment is non-zero, the system should
report ±inf percent_change and the correct direction (not "similar").
"""

from __future__ import annotations

import math
from datetime import datetime
from types import SimpleNamespace

from polyzymd.analyses.base import (
    ComparisonResult,
    ConditionSummary,
    MetricValue,
    PairwiseResult,
)
from polyzymd.analyses.contacts._comparison_results import AggregateComparisonResult
from polyzymd.analyses.distances._comparison_results import (
    DistanceComparisonResult,
    DistanceConditionSummary,
    DistancePairSummary,
    DistancePairwiseComparison,
)
from polyzymd.analyses.distances._formatters import (
    format_distances_console_table,
    format_distances_markdown,
)
from polyzymd.analyses.rg import RgAnalysis
from polyzymd.analyses.rg._comparison_results import RgRunPairwiseComparison
from polyzymd.analyses.rmsd._comparison_results import RMSDRunPairwiseComparison
from polyzymd.analyses.sasa._comparison_results import SASARunPairwiseComparison
from polyzymd.analyses.stats import (
    format_pct,
    format_scalar_comparison,
    interpret_direction,
    pairwise_comparisons,
)
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


def test_rg_pairwise_inf_round_trip_json() -> None:
    """Rg pairwise model should preserve infinite percent_change across JSON round-trip."""
    comparison = RgRunPairwiseComparison(
        run_label="protein_backbone",
        condition_a="Control",
        condition_b="Treatment",
        t_statistic=1.0,
        p_value=0.01,
        cohens_d=1.5,
        effect_interpretation="large",
        direction="expansion",
        significant=True,
        percent_change=math.inf,
    )

    payload = comparison.model_dump_json()
    loaded = RgRunPairwiseComparison.model_validate_json(payload)

    assert math.isinf(loaded.percent_change)
    assert loaded.percent_change > 0


def test_rmsd_pairwise_inf_round_trip_json() -> None:
    """RMSD pairwise model should preserve infinite percent_change across JSON round-trip."""
    comparison = RMSDRunPairwiseComparison(
        run_label="protein_backbone",
        condition_a="Control",
        condition_b="Treatment",
        t_statistic=1.0,
        p_value=0.01,
        cohens_d=1.5,
        effect_interpretation="large",
        direction="destabilizing",
        significant=True,
        percent_change=math.inf,
    )

    payload = comparison.model_dump_json()
    loaded = RMSDRunPairwiseComparison.model_validate_json(payload)

    assert math.isinf(loaded.percent_change)
    assert loaded.percent_change > 0


def test_sasa_pairwise_inf_round_trip_json() -> None:
    """SASA pairwise model should preserve infinite percent_change across JSON round-trip."""
    comparison = SASARunPairwiseComparison(
        run_label="surface",
        condition_a="Control",
        condition_b="Treatment",
        t_statistic=1.0,
        p_value=0.01,
        cohens_d=1.5,
        effect_interpretation="large",
        direction="increased",
        significant=True,
        percent_change=math.inf,
    )

    payload = comparison.model_dump_json()
    loaded = SASARunPairwiseComparison.model_validate_json(payload)

    assert math.isinf(loaded.percent_change)
    assert loaded.percent_change > 0


def test_distances_pairwise_inf_round_trip_json() -> None:
    """Distances pairwise model should preserve infinite percent_change across JSON round-trip."""
    comparison = DistancePairwiseComparison(
        pair_label="Catalytic distance",
        condition_a="Control",
        condition_b="Treatment",
        distance_t_statistic=1.0,
        distance_p_value=0.01,
        distance_cohens_d=1.5,
        distance_effect_interpretation="large",
        distance_direction="farther",
        distance_significant=True,
        distance_percent_change=math.inf,
    )

    payload = comparison.model_dump_json()
    loaded = DistancePairwiseComparison.model_validate_json(payload)

    assert math.isinf(loaded.distance_percent_change)
    assert loaded.distance_percent_change > 0


def test_contacts_aggregate_inf_round_trip_json() -> None:
    """Contacts aggregate model should preserve infinite percent_change in JSON."""
    comparison = AggregateComparisonResult(
        metric="coverage",
        condition_a="Control",
        condition_b="Treatment",
        condition_a_mean=0.0,
        condition_a_sem=0.0,
        condition_b_mean=1.0,
        condition_b_sem=0.1,
        t_statistic=1.0,
        p_value=0.01,
        cohens_d=1.5,
        effect_size_interpretation="large",
        significant=True,
        percent_change=math.inf,
        direction="increased",
    )

    payload = comparison.model_dump_json()
    loaded = AggregateComparisonResult.model_validate_json(payload)

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


def test_polymer_bridging_safe_pairwise_zero_control_direction() -> None:
    """Polymer bridging safe pairwise helper should not report unchanged for zero-control inf."""
    from polyzymd.analyses.polymer_bridging import _safe_pairwise_comparisons

    metrics = {
        "Control": MetricValue(
            name="multisite_fraction",
            mean=0.0,
            sem=0.0,
            replicate_values=[0.0, 0.0, 0.0],
            higher_is_better=True,
            direction_labels=("lower", "similar", "higher"),
        ),
        "Treatment": MetricValue(
            name="multisite_fraction",
            mean=0.5,
            sem=0.1,
            replicate_values=[0.4, 0.5, 0.6],
            higher_is_better=True,
            direction_labels=("lower", "similar", "higher"),
        ),
    }

    results = _safe_pairwise_comparisons(metrics, control_label="Control")
    assert len(results) == 1
    comp = results[0]
    assert math.isinf(comp.percent_change)
    assert comp.direction == "higher"
    assert comp.direction != "similar"


def test_rg_compare_run_handles_zero_control_infinite_direction() -> None:
    """Custom Rg compare path should classify zero-control +inf as expansion."""
    run_a = SimpleNamespace(mean_rg=0.0, per_replicate_means=[0.0, 0.0, 0.0])
    run_b = SimpleNamespace(mean_rg=4.0, per_replicate_means=[3.8, 4.0, 4.2])

    comp = RgAnalysis._compare_run(
        run_label="protein",
        condition_a="Control",
        condition_b="Treatment",
        run_a=run_a,
        run_b=run_b,
    )

    assert math.isinf(comp.percent_change)
    assert comp.percent_change > 0
    assert comp.direction == "expansion"


def test_exposure_pairwise_comparison_zero_control_direction() -> None:
    """Exposure custom comparison path should classify zero-control +inf as increased."""
    from polyzymd.analyses.exposure import ExposureAnalysis

    cond_a = SimpleNamespace(
        label="Control",
        replicate_values=[0.0, 0.0, 0.0],
        mean_chaperone_fraction=0.0,
    )
    cond_b = SimpleNamespace(
        label="Treatment",
        replicate_values=[0.4, 0.5, 0.6],
        mean_chaperone_fraction=0.5,
    )

    comp = ExposureAnalysis._compare_pair(cond_a, cond_b)
    assert isinstance(comp, PairwiseResult)
    assert math.isinf(comp.percent_change)
    assert comp.percent_change > 0
    assert comp.direction == "increased"


def test_distances_formatter_zero_control_emits_infinity_with_direction() -> None:
    """Distances formatter should report direction when control distance is zero."""
    pair_label = "Catalytic distance"
    control = DistanceConditionSummary(
        label="Control",
        config_path="/tmp/control.yaml",
        n_replicates=3,
        replicate_values=[],
        pair_summaries=[
            DistancePairSummary(
                label=pair_label,
                selection_a="chainid A and name CA",
                selection_b="chainid B and name C1",
                threshold=3.5,
                mean_distance=0.0,
                sem_distance=0.0,
                fraction_below_threshold=0.0,
                sem_fraction_below=0.0,
                per_replicate_means=[0.0, 0.0, 0.0],
                per_replicate_fractions=[0.0, 0.0, 0.0],
            )
        ],
    )
    treatment = DistanceConditionSummary(
        label="Treatment",
        config_path="/tmp/treatment.yaml",
        n_replicates=3,
        replicate_values=[],
        pair_summaries=[
            DistancePairSummary(
                label=pair_label,
                selection_a="chainid A and name CA",
                selection_b="chainid B and name C1",
                threshold=3.5,
                mean_distance=-1.0,
                sem_distance=0.1,
                fraction_below_threshold=1.0,
                sem_fraction_below=0.0,
                per_replicate_means=[-1.0, -1.0, -1.0],
                per_replicate_fractions=[1.0, 1.0, 1.0],
            )
        ],
    )
    result = DistanceComparisonResult(
        name="distance-zero-control",
        n_pairs=1,
        pair_labels=[pair_label],
        control_label="Control",
        conditions=[control, treatment],
        pairwise_comparisons=[],
        anova_by_pair=None,
        ranking_by_pair={pair_label: ["Treatment", "Control"]},
        fraction_ranking_by_pair={pair_label: ["Treatment", "Control"]},
        ranking=[],
        equilibration_time="0ns",
        created_at=datetime(2026, 1, 1, 0, 0, 0),
        polyzymd_version="test",
    )

    console_output = format_distances_console_table(result, show_pairwise=False, show_anova=False)
    markdown_output = format_distances_markdown(result, show_pairwise=False, show_anova=False)

    assert "closer than control" in console_output
    assert "closer than control" in markdown_output
