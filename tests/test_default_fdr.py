"""Tests for BH/FDR behavior in default scalar comparisons."""

from __future__ import annotations

import json

import pytest

from polyzymd.analyses.base import ComparisonResult, ConditionSummary, MetricValue, PairwiseResult
from polyzymd.analyses.stats import format_scalar_comparison, pairwise_comparisons
from polyzymd.analyses.shared.inferential_statistics import EffectSize, TTestResult


def _metric(mean: float, values: list[float]) -> MetricValue:
    """Build a simple metric value for pairwise tests."""
    return MetricValue(
        name="metric",
        mean=mean,
        sem=0.1,
        replicate_values=values,
        higher_is_better=True,
        direction_labels=("lower", "unchanged", "higher"),
    )


def test_pairwise_comparisons_applies_bh_and_sets_adjusted_significance(monkeypatch) -> None:
    """Pairwise results should include BH-adjusted p-values and adjusted significance."""
    p_values = [0.01, 0.03, 0.04]

    def _fake_ttest(_group1, _group2):
        return TTestResult(t_statistic=1.0, p_value=p_values.pop(0))

    def _fake_effect(_group1, _group2):
        return EffectSize(cohens_d=0.5, interpretation="medium", direction="higher")

    monkeypatch.setattr(
        "polyzymd.analyses.shared.inferential_statistics.independent_ttest",
        _fake_ttest,
    )
    monkeypatch.setattr("polyzymd.analyses.shared.inferential_statistics.cohens_d", _fake_effect)

    metrics = {
        "A": _metric(1.0, [1.0, 1.1, 0.9]),
        "B": _metric(1.2, [1.1, 1.2, 1.3]),
        "C": _metric(1.4, [1.3, 1.4, 1.5]),
    }

    results = pairwise_comparisons(metrics, fdr_alpha=0.035)

    assert len(results) == 3
    assert all(r.p_value_adjusted is not None for r in results)
    assert any(r.p_value_adjusted != r.p_value for r in results)
    # BH-adjusted p-values are [0.03, 0.04, 0.04] for raw [0.01, 0.03, 0.04]
    assert [r.significant for r in results] == [True, False, False]


def test_comparison_result_round_trip_preserves_adjusted_pvalue() -> None:
    """ComparisonResult serialization should preserve adjusted p-values."""
    result = ComparisonResult(
        analysis_type="test",
        name="fdr-roundtrip",
        fdr_alpha=0.05,
        conditions=[ConditionSummary(label="A", n_replicates=3)],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="A",
                condition_b="B",
                metric="metric",
                t_statistic=1.0,
                p_value=0.02,
                p_value_adjusted=0.03,
                cohens_d=0.5,
                effect_size_interpretation="medium",
                direction="higher",
                significant=True,
                percent_change=10.0,
            )
        ],
    )

    loaded = ComparisonResult.model_validate_json(result.model_dump_json())

    assert loaded.fdr_alpha == pytest.approx(0.05)
    assert loaded.pairwise_comparisons[0].p_value_adjusted == pytest.approx(0.03)


def test_format_scalar_comparison_includes_adjusted_pvalues_in_text_and_markdown() -> None:
    """Formatter should render raw and adjusted p-values when available."""
    result = ComparisonResult(
        analysis_type="test",
        name="fdr-format",
        control_label="A",
        fdr_alpha=0.05,
        conditions=[
            ConditionSummary(label="A", n_replicates=3, metric_mean=1.0, metric_sem=0.1),
            ConditionSummary(label="B", n_replicates=3, metric_mean=1.2, metric_sem=0.1),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="A",
                condition_b="B",
                metric="metric",
                t_statistic=1.0,
                p_value=0.03,
                p_value_adjusted=0.04,
                cohens_d=0.5,
                effect_size_interpretation="medium",
                direction="higher",
                significant=True,
                percent_change=20.0,
            )
        ],
        ranking=["B", "A"],
        equilibration_time="0ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )

    text_output = format_scalar_comparison(result, output_format="text", metric_key="metric")
    markdown_output = format_scalar_comparison(
        result, output_format="markdown", metric_key="metric"
    )

    assert "p (adj)" in text_output
    assert "0.0300" in text_output
    assert "0.0400" in text_output
    assert "p (adj)" in markdown_output
    assert "0.0300" in markdown_output
    assert "0.0400" in markdown_output


def test_format_scalar_comparison_backward_compatible_without_adjusted_pvalues() -> None:
    """Formatter should keep legacy output when adjusted p-values are absent."""
    legacy_payload = {
        "analysis_type": "test",
        "name": "legacy",
        "control_label": "A",
        "conditions": [
            {"label": "A", "n_replicates": 3, "metric_mean": 1.0, "metric_sem": 0.1},
            {"label": "B", "n_replicates": 3, "metric_mean": 1.2, "metric_sem": 0.1},
        ],
        "pairwise_comparisons": [
            {
                "condition_a": "A",
                "condition_b": "B",
                "metric": "metric",
                "t_statistic": 1.0,
                "p_value": 0.03,
                "cohens_d": 0.5,
                "effect_size_interpretation": "medium",
                "direction": "higher",
                "significant": True,
                "percent_change": 20.0,
            }
        ],
        "anova": None,
        "ranking": ["B", "A"],
        "rankings_by_metric": None,
        "equilibration_time": "0ns",
        "created_at": "2026-01-01T00:00:00",
        "polyzymd_version": "test",
    }
    legacy_result = ComparisonResult.model_validate_json(json.dumps(legacy_payload))

    text_output = format_scalar_comparison(legacy_result, output_format="text", metric_key="metric")
    markdown_output = format_scalar_comparison(
        legacy_result,
        output_format="markdown",
        metric_key="metric",
    )

    assert "p (adj)" not in text_output
    assert "p (adj)" not in markdown_output
    assert "p-value" in text_output
    assert "| Comparison | % Change | p-value |" in markdown_output
