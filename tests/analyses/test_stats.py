"""Tests for BH/FDR behavior in default scalar comparisons."""

from __future__ import annotations

import json
import math

import pytest
from pydantic import ValidationError

from polyzymd.analyses.base import (
    ANOVAResult,
    ComparisonResult,
    ConditionSummary,
    MetricValue,
    PairwiseResult,
)
from polyzymd.analyses.shared.inferential_statistics import EffectSize, TTestResult
from polyzymd.analyses.stats import (
    anova_test,
    default_scalar_comparison,
    format_scalar_comparison,
    pairwise_comparisons,
    rank_conditions,
)


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


def test_pairwise_comparisons_returns_raw_pvalues_only(monkeypatch) -> None:
    """Pairwise results should keep raw p-values without BH adjustment."""
    p_values = [0.01, 0.03, 0.04]

    def _fake_ttest(_group1, _group2, method="student"):
        del method
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

    results = pairwise_comparisons(metrics)

    assert len(results) == 3
    assert all(r.p_value_adjusted is None for r in results)
    assert [r.significant for r in results] == [True, True, True]


def test_default_scalar_comparison_applies_bh_across_full_family(monkeypatch) -> None:
    """BH should run once across all metric pairwise tests."""
    p_values = [0.01, 0.02, 0.03, 0.04]

    def _fake_ttest(_group1, _group2, method="student"):
        del method
        return TTestResult(t_statistic=1.0, p_value=p_values.pop(0))

    def _fake_effect(_group1, _group2):
        return EffectSize(cohens_d=0.5, interpretation="medium", direction="higher")

    monkeypatch.setattr(
        "polyzymd.analyses.shared.inferential_statistics.independent_ttest",
        _fake_ttest,
    )
    monkeypatch.setattr("polyzymd.analyses.shared.inferential_statistics.cohens_d", _fake_effect)

    metrics_by_condition = {
        "Control": {
            "metric_a": _metric(1.0, [1.0, 1.1, 0.9]),
            "metric_b": _metric(2.0, [1.9, 2.0, 2.1]),
        },
        "Treatment 1": {
            "metric_a": _metric(1.2, [1.1, 1.2, 1.3]),
            "metric_b": _metric(2.2, [2.1, 2.2, 2.3]),
        },
        "Treatment 2": {
            "metric_a": _metric(1.3, [1.2, 1.3, 1.4]),
            "metric_b": _metric(2.3, [2.2, 2.3, 2.4]),
        },
    }

    result = default_scalar_comparison(
        analysis_name="test",
        project_name="Full Family BH",
        metrics_by_condition=metrics_by_condition,
        control_label="Control",
        fdr_alpha=0.035,
    )

    adjusted = [r.p_value_adjusted for r in result.pairwise_comparisons]
    assert adjusted == pytest.approx([0.04, 0.04, 0.04, 0.04])
    assert [r.significant for r in result.pairwise_comparisons] == [False, False, False, False]


def test_default_scalar_singleton_pairwise_not_testable() -> None:
    """Singleton replicate comparisons should not be treated as significant tests."""
    metrics_by_condition = {
        "Control": {"metric": _metric(1.0, [1.0])},
        "Treatment": {"metric": _metric(2.0, [2.0])},
    }

    result = default_scalar_comparison(
        analysis_name="test",
        project_name="Singleton",
        metrics_by_condition=metrics_by_condition,
        control_label="Control",
    )

    comparison = result.pairwise_comparisons[0]
    assert comparison.testable is False
    assert comparison.significant is False
    assert comparison.p_value_adjusted is None
    assert comparison.note is not None
    assert math.isnan(comparison.p_value)


def test_format_scalar_singleton_sem_rendered_as_not_available() -> None:
    """Scalar formatter should not display singleton SEM as a numeric value."""
    metrics_by_condition = {
        "Control": {"metric": _metric(1.0, [1.0])},
        "Treatment": {"metric": _metric(2.0, [2.0])},
    }
    result = default_scalar_comparison(
        analysis_name="test",
        project_name="Singleton format",
        metrics_by_condition=metrics_by_condition,
        control_label="Control",
    )

    text_output = format_scalar_comparison(result, output_format="text", metric_key="metric")
    markdown_output = format_scalar_comparison(
        result,
        output_format="markdown",
        metric_key="metric",
    )

    assert "n/a" in text_output
    assert "SEM: n/a (single replicate; not estimable)" in text_output
    assert "0.1000" not in text_output
    assert "| 1.0000 | n/a | 1 |" in markdown_output
    assert "*SEM: n/a (single replicate; not estimable).*" in markdown_output
    assert "0.1000" not in markdown_output


def test_default_scalar_singleton_anova_not_testable() -> None:
    """ANOVA should report singleton conditions as not testable."""
    metrics_by_condition = {
        "A": _metric(1.0, [1.0]),
        "B": _metric(2.0, [2.0]),
        "C": _metric(3.0, [3.0]),
    }

    anova = anova_test(metrics_by_condition, "metric")

    assert anova is not None
    assert anova.testable is False
    assert anova.significant is False
    assert anova.note is not None
    assert math.isnan(anova.p_value)


def test_default_scalar_fdr_ignores_non_testable_results(monkeypatch) -> None:
    """BH correction should skip non-testable singleton comparisons."""
    p_values = [0.01]

    def _fake_ttest(group1, group2, method="student"):
        del method
        if len(group1) < 2 or len(group2) < 2:
            return TTestResult(t_statistic=float("nan"), p_value=float("nan"))
        return TTestResult(t_statistic=3.0, p_value=p_values.pop(0))

    def _fake_effect(_group1, _group2):
        return EffectSize(cohens_d=0.5, interpretation="medium", direction="higher")

    monkeypatch.setattr(
        "polyzymd.analyses.shared.inferential_statistics.independent_ttest",
        _fake_ttest,
    )
    monkeypatch.setattr("polyzymd.analyses.shared.inferential_statistics.cohens_d", _fake_effect)

    metrics_by_condition = {
        "Control": {"metric": _metric(1.0, [1.0, 1.1])},
        "Singleton": {"metric": _metric(1.2, [1.2])},
        "Treatment": {"metric": _metric(1.4, [1.3, 1.5])},
    }

    result = default_scalar_comparison(
        analysis_name="test",
        project_name="Mixed replicate counts",
        metrics_by_condition=metrics_by_condition,
        control_label="Control",
        fdr_alpha=0.05,
    )

    singleton_comp, treatment_comp = result.pairwise_comparisons
    assert singleton_comp.testable is False
    assert singleton_comp.p_value_adjusted is None
    assert singleton_comp.significant is False
    assert treatment_comp.testable is True
    assert treatment_comp.p_value_adjusted == pytest.approx(0.01)


def test_default_scalar_comparison_threads_ttest_method(monkeypatch) -> None:
    """default_scalar_comparison should pass ttest_method to pairwise tests."""
    seen_methods: list[str] = []

    def _fake_ttest(_group1, _group2, method="student"):
        seen_methods.append(method)
        return TTestResult(t_statistic=1.0, p_value=0.5)

    def _fake_effect(_group1, _group2):
        return EffectSize(cohens_d=0.1, interpretation="small", direction="higher")

    monkeypatch.setattr(
        "polyzymd.analyses.shared.inferential_statistics.independent_ttest",
        _fake_ttest,
    )
    monkeypatch.setattr("polyzymd.analyses.shared.inferential_statistics.cohens_d", _fake_effect)

    metrics_by_condition = {
        "Control": {"metric_a": _metric(1.0, [1.0, 1.1, 0.9])},
        "Treatment": {"metric_a": _metric(1.2, [1.1, 1.2, 1.3])},
    }
    result = default_scalar_comparison(
        analysis_name="test",
        project_name="Method threading",
        metrics_by_condition=metrics_by_condition,
        control_label="Control",
        ttest_method="welch",
    )

    assert seen_methods == ["welch"]
    assert result.ttest_method == "welch"


def test_default_scalar_comparison_threads_posthoc_method(monkeypatch) -> None:
    """default_scalar_comparison should pass posthoc_method to pairwise tests."""
    seen_posthoc_methods: list[str] = []

    def _fake_pairwise(
        metrics,
        control_label=None,
        ttest_method="student",
        posthoc_method="ttest_bh",
        fdr_alpha=0.05,
    ):
        del metrics, control_label, ttest_method, fdr_alpha
        seen_posthoc_methods.append(posthoc_method)
        return [
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment 1",
                metric="metric_a",
                t_statistic=float("nan"),
                p_value=0.03,
                p_value_adjusted=None,
                posthoc_method=posthoc_method,
                cohens_d=0.5,
                effect_size_interpretation="medium",
                direction="higher",
                significant=True,
                percent_change=20.0,
            )
        ]

    monkeypatch.setattr("polyzymd.analyses.stats.pairwise_comparisons", _fake_pairwise)

    metrics_by_condition = {
        "Control": {"metric_a": _metric(1.0, [1.0, 1.1, 0.9])},
        "Treatment 1": {"metric_a": _metric(1.2, [1.1, 1.2, 1.3])},
        "Treatment 2": {"metric_a": _metric(1.3, [1.2, 1.3, 1.4])},
    }
    result = default_scalar_comparison(
        analysis_name="test",
        project_name="Posthoc threading",
        metrics_by_condition=metrics_by_condition,
        posthoc_method="tukey_hsd",
    )

    assert seen_posthoc_methods == ["tukey_hsd"]
    assert result.posthoc_method == "tukey_hsd"


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


def test_comparison_result_rejects_missing_adjusted_pvalues() -> None:
    """ComparisonResult should reject pairwise entries without adjusted p-values."""
    noncanonical_payload = {
        "analysis_type": "test",
        "name": "unstamped",
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
    with pytest.raises(ValidationError, match="p_value_adjusted"):
        ComparisonResult.model_validate_json(json.dumps(noncanonical_payload))


def test_bh_boundary_inclusive(monkeypatch) -> None:
    """Adjusted p-value exactly equal to fdr_alpha should be significant."""

    def _fake_ttest(_group1, _group2, method="student"):
        del method
        return TTestResult(t_statistic=1.0, p_value=0.05)

    def _fake_effect(_group1, _group2):
        return EffectSize(cohens_d=0.5, interpretation="medium", direction="higher")

    monkeypatch.setattr(
        "polyzymd.analyses.shared.inferential_statistics.independent_ttest",
        _fake_ttest,
    )
    monkeypatch.setattr(
        "polyzymd.analyses.shared.inferential_statistics.cohens_d",
        _fake_effect,
    )

    metrics_ctrl = {"m": _metric(1.0, [1.0, 1.1, 1.2])}
    metrics_trt = {"m": _metric(2.0, [2.0, 2.1, 2.2])}

    result = default_scalar_comparison(
        analysis_name="test",
        project_name="BH inclusive boundary",
        metrics_by_condition={"ctrl": metrics_ctrl, "trt": metrics_trt},
        control_label="ctrl",
        fdr_alpha=0.05,
        posthoc_method="ttest_bh",
    )

    assert len(result.pairwise_comparisons) == 1
    pairwise = result.pairwise_comparisons[0]
    assert pairwise.p_value_adjusted == pytest.approx(0.05)
    assert pairwise.significant is True


def test_default_scalar_comparison_tukey_mode() -> None:
    """Tukey HSD mode should produce results with NaN t-statistics."""
    metrics = {
        "ctrl": {"m": _metric(1.0, [1.0, 1.1, 1.2])},
        "trt1": {"m": _metric(2.0, [2.0, 2.1, 2.2])},
        "trt2": {"m": _metric(3.0, [3.0, 3.1, 3.2])},
    }
    result = default_scalar_comparison(
        analysis_name="test",
        project_name="Tukey mode",
        metrics_by_condition=metrics,
        control_label="ctrl",
        posthoc_method="tukey_hsd",
    )

    assert result.posthoc_method == "tukey_hsd"
    assert len(result.pairwise_comparisons) == 2
    pair_identities = {(p.condition_a, p.condition_b) for p in result.pairwise_comparisons}
    assert pair_identities == {("ctrl", "trt1"), ("ctrl", "trt2")}
    assert ("trt1", "trt2") not in pair_identities
    for pairwise in result.pairwise_comparisons:
        assert pairwise.posthoc_method == "tukey_hsd"
        assert math.isnan(pairwise.t_statistic)
        assert pairwise.p_value_adjusted == pytest.approx(pairwise.p_value)
        assert pairwise.p_value >= 0


def test_pairwise_comparisons_tukey_uses_fdr_alpha_boundary(monkeypatch) -> None:
    """Tukey significance should be inclusive at fdr_alpha boundary."""

    class _FakeTukeyResult:
        def __init__(self, group_i: int, group_j: int, p_value: float):
            self.group_i = group_i
            self.group_j = group_j
            self.p_value = p_value

    def _fake_tukey_hsd(*_groups):
        return [_FakeTukeyResult(group_i=0, group_j=1, p_value=0.05)]

    def _fake_effect(_group1, _group2):
        return EffectSize(cohens_d=0.5, interpretation="medium", direction="higher")

    monkeypatch.setattr(
        "polyzymd.analyses.shared.inferential_statistics.tukey_hsd", _fake_tukey_hsd
    )
    monkeypatch.setattr("polyzymd.analyses.shared.inferential_statistics.cohens_d", _fake_effect)

    metrics = {
        "ctrl": _metric(1.0, [1.0, 1.1, 1.2]),
        "trt": _metric(2.0, [2.0, 2.1, 2.2]),
    }
    results = pairwise_comparisons(metrics, posthoc_method="tukey_hsd", fdr_alpha=0.05)

    assert len(results) == 1
    assert results[0].p_value_adjusted == pytest.approx(0.05)
    assert results[0].significant is True


def test_pairwise_comparisons_tukey_normalizes_control_orientation(monkeypatch) -> None:
    """Tukey comparisons should keep control as condition_a baseline."""

    class _FakeTukeyResult:
        def __init__(self, group_i: int, group_j: int, p_value: float):
            self.group_i = group_i
            self.group_j = group_j
            self.p_value = p_value

    def _fake_tukey_hsd(*_groups):
        return [_FakeTukeyResult(group_i=0, group_j=1, p_value=0.01)]

    def _fake_effect(_group1, _group2):
        return EffectSize(cohens_d=0.5, interpretation="medium", direction="higher")

    monkeypatch.setattr(
        "polyzymd.analyses.shared.inferential_statistics.tukey_hsd", _fake_tukey_hsd
    )
    monkeypatch.setattr("polyzymd.analyses.shared.inferential_statistics.cohens_d", _fake_effect)

    metrics = {
        "Treatment": _metric(2.0, [2.0, 2.1, 2.2]),
        "Control": _metric(1.0, [1.0, 1.1, 1.2]),
    }
    results = pairwise_comparisons(metrics, control_label="Control", posthoc_method="tukey_hsd")

    assert len(results) == 1
    assert results[0].condition_a == "Control"
    assert results[0].condition_b == "Treatment"


def test_default_scalar_comparison_raises_on_inconsistent_metric_keys() -> None:
    """Default scalar comparison should reject mismatched metric key sets."""
    metrics_by_condition = {
        "Control": {
            "metric_a": _metric(1.0, [1.0, 1.1, 0.9]),
            "metric_b": _metric(2.0, [1.9, 2.0, 2.1]),
        },
        "Treatment": {
            "metric_a": _metric(1.2, [1.1, 1.2, 1.3]),
            "metric_c": _metric(3.2, [3.1, 3.2, 3.3]),
        },
    }

    with pytest.raises(ValueError, match="Inconsistent metric keys across conditions"):
        default_scalar_comparison(
            analysis_name="test",
            project_name="Inconsistent keys",
            metrics_by_condition=metrics_by_condition,
            control_label="Control",
        )


def test_rank_conditions_raises_on_inconsistent_higher_is_better() -> None:
    """Ranker should reject conflicting higher_is_better declarations."""
    metrics = {
        "A": MetricValue("m", 1.0, 0.1, [1.0], higher_is_better=True),
        "B": MetricValue("m", 2.0, 0.1, [2.0], higher_is_better=False),
    }

    with pytest.raises(ValueError, match="Inconsistent MetricValue.higher_is_better"):
        rank_conditions(metrics)


def test_rank_conditions_raises_on_mixed_none_and_bool_higher_is_better() -> None:
    """Ranker should reject mixed None and bool higher_is_better values."""
    metrics = {
        "A": MetricValue("m", 1.0, 0.1, [1.0], higher_is_better=None),
        "B": MetricValue("m", 2.0, 0.1, [2.0], higher_is_better=True),
    }

    with pytest.raises(ValueError, match="Inconsistent MetricValue.higher_is_better"):
        rank_conditions(metrics)


def test_formatter_uses_non_default_fdr_alpha_thresholds() -> None:
    """Formatter should display custom fdr_alpha thresholds in text and markdown."""
    result = ComparisonResult(
        analysis_type="test",
        name="alpha-format",
        control_label="A",
        fdr_alpha=0.01,
        posthoc_method="ttest_bh",
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
                p_value_adjusted=None,
                posthoc_method="ttest_bh",
                cohens_d=0.5,
                effect_size_interpretation="medium",
                direction="higher",
                significant=False,
                percent_change=20.0,
            )
        ],
        anova=[ANOVAResult(metric="metric", f_statistic=3.2, p_value=0.04, significant=False)],
        ranking=["B", "A"],
        rankings_by_metric={"metric": ["B", "A"]},
        equilibration_time="0ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )

    text_output = format_scalar_comparison(result, output_format="text", metric_key="metric")
    markdown_output = format_scalar_comparison(
        result, output_format="markdown", metric_key="metric"
    )

    assert "* p <= 0.01" in text_output
    assert "Significant: No (alpha=0.01)" in text_output
    assert "*Significance uses raw p <= 0.01.*" in markdown_output


def test_formatter_uses_non_default_fdr_alpha_for_tukey_note() -> None:
    """Markdown formatter should display custom threshold for Tukey mode."""
    result = ComparisonResult(
        analysis_type="test",
        name="alpha-format-tukey",
        control_label="A",
        fdr_alpha=0.02,
        posthoc_method="tukey_hsd",
        conditions=[
            ConditionSummary(label="A", n_replicates=3, metric_mean=1.0, metric_sem=0.1),
            ConditionSummary(label="B", n_replicates=3, metric_mean=1.2, metric_sem=0.1),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="A",
                condition_b="B",
                metric="metric",
                t_statistic=float("nan"),
                p_value=0.02,
                p_value_adjusted=0.02,
                posthoc_method="tukey_hsd",
                cohens_d=0.5,
                effect_size_interpretation="medium",
                direction="higher",
                significant=True,
                percent_change=20.0,
            )
        ],
        ranking=["B", "A"],
        rankings_by_metric={"metric": ["B", "A"]},
        equilibration_time="0ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )

    text_output = format_scalar_comparison(result, output_format="text", metric_key="metric")
    markdown_output = format_scalar_comparison(
        result, output_format="markdown", metric_key="metric"
    )

    assert "* Tukey HSD p <= 0.02" in text_output
    assert "*Significance uses Tukey HSD p <= 0.02.*" in markdown_output


class TestStats:
    """Test the shared statistical utility functions."""

    def test_interpret_direction(self):
        from polyzymd.analyses.stats import interpret_direction

        assert interpret_direction(5.0) == "increased"
        assert interpret_direction(-5.0) == "decreased"
        assert interpret_direction(0.5) == "unchanged"
        assert interpret_direction(0.5, threshold=0.3) == "increased"

    def test_interpret_direction_custom_labels(self):
        from polyzymd.analyses.stats import interpret_direction

        labels = ("stabilizing", "unchanged", "destabilizing")
        assert interpret_direction(-10.0, labels) == "stabilizing"
        assert interpret_direction(10.0, labels) == "destabilizing"

    def test_pairwise_comparisons_control(self):
        from polyzymd.analyses.stats import pairwise_comparisons

        metrics = {
            "Control": MetricValue("m", 1.0, 0.1, [0.9, 1.0, 1.1]),
            "Treatment A": MetricValue("m", 0.5, 0.05, [0.45, 0.5, 0.55]),
            "Treatment B": MetricValue("m", 0.8, 0.08, [0.75, 0.8, 0.85]),
        }
        results = pairwise_comparisons(metrics, control_label="Control")
        assert len(results) == 2
        assert all(r.condition_a == "Control" for r in results)
        assert {r.condition_b for r in results} == {"Treatment A", "Treatment B"}
        for result in results:
            assert hasattr(result, "t_statistic")
            assert hasattr(result, "p_value")
            assert hasattr(result, "cohens_d")
            assert hasattr(result, "percent_change")
            assert hasattr(result, "significant")

    def test_pairwise_comparisons_all_pairs(self):
        from polyzymd.analyses.stats import pairwise_comparisons

        metrics = {
            "A": MetricValue("m", 1.0, 0.1, [0.9, 1.0, 1.1]),
            "B": MetricValue("m", 2.0, 0.2, [1.8, 2.0, 2.2]),
            "C": MetricValue("m", 3.0, 0.3, [2.7, 3.0, 3.3]),
        }
        results = pairwise_comparisons(metrics)
        assert len(results) == 3

    def test_anova_test_with_3_conditions(self):
        from polyzymd.analyses.stats import anova_test

        metrics = {
            "A": MetricValue("m", 1.0, 0.1, [0.9, 1.0, 1.1]),
            "B": MetricValue("m", 5.0, 0.2, [4.8, 5.0, 5.2]),
            "C": MetricValue("m", 3.0, 0.3, [2.7, 3.0, 3.3]),
        }
        result = anova_test(metrics, "test_metric")
        assert result is not None
        assert hasattr(result, "f_statistic")
        assert hasattr(result, "p_value")
        assert result.significant is True

    def test_anova_test_too_few_conditions(self):
        from polyzymd.analyses.stats import anova_test

        metrics = {
            "A": MetricValue("m", 1.0, 0.1, [0.9, 1.1]),
            "B": MetricValue("m", 2.0, 0.2, [1.8, 2.2]),
        }
        result = anova_test(metrics)
        assert result is None

    def test_rank_conditions_higher_is_better(self):
        from polyzymd.analyses.stats import rank_conditions

        metrics = {
            "Low": MetricValue("m", 1.0, 0.1, [1.0], higher_is_better=True),
            "High": MetricValue("m", 3.0, 0.1, [3.0], higher_is_better=True),
            "Mid": MetricValue("m", 2.0, 0.1, [2.0], higher_is_better=True),
        }
        ranking = rank_conditions(metrics)
        assert ranking == ["High", "Mid", "Low"]

    def test_rank_conditions_lower_is_better(self):
        from polyzymd.analyses.stats import rank_conditions

        metrics = {
            "Low": MetricValue("m", 1.0, 0.1, [1.0], higher_is_better=False),
            "High": MetricValue("m", 3.0, 0.1, [3.0], higher_is_better=False),
            "Mid": MetricValue("m", 2.0, 0.1, [2.0], higher_is_better=False),
        }
        ranking = rank_conditions(metrics)
        assert ranking == ["Low", "Mid", "High"]

    def test_default_scalar_comparison(self):
        from polyzymd.analyses.stats import default_scalar_comparison

        metrics_by_condition = {
            "Control": {
                "metric_a": MetricValue("metric_a", 1.0, 0.1, [0.9, 1.0, 1.1]),
            },
            "Treatment": {
                "metric_a": MetricValue("metric_a", 0.5, 0.05, [0.45, 0.5, 0.55]),
            },
        }
        result = default_scalar_comparison(
            analysis_name="test",
            project_name="Test Project",
            metrics_by_condition=metrics_by_condition,
            control_label="Control",
            equilibration="10ns",
        )
        assert result.analysis_type == "test"
        assert result.name == "Test Project"
        assert len(result.pairwise_comparisons) == 1
        assert result.anova is None
        assert result.ranking == ["Control", "Treatment"] or result.ranking == [
            "Treatment",
            "Control",
        ]


class TestAnovaAlphaThreading:
    """ANOVA significance should respect fdr_alpha."""

    def test_anova_uses_configured_alpha(self):
        """ANOVA significant flag should reflect the pipeline's fdr_alpha."""
        # 3 conditions with clearly different means → p << 0.05
        metrics_by_condition = {
            "A": {"m": _metric(1.0, [0.9, 1.0, 1.1])},
            "B": {"m": _metric(5.0, [4.9, 5.0, 5.1])},
            "C": {"m": _metric(3.0, [2.9, 3.0, 3.1])},
        }
        result = default_scalar_comparison(
            analysis_name="test",
            project_name="ANOVA alpha",
            metrics_by_condition=metrics_by_condition,
            fdr_alpha=0.05,
        )
        assert result.anova is not None
        assert result.anova[0].significant is True

    def test_anova_not_significant_with_strict_alpha(self):
        """ANOVA p-value > very strict alpha should be not significant."""
        # 3 groups with overlapping values: p is typically moderate
        metrics_by_condition = {
            "A": {"m": _metric(1.0, [0.5, 1.0, 1.5])},
            "B": {"m": _metric(1.5, [1.0, 1.5, 2.0])},
            "C": {"m": _metric(1.2, [0.7, 1.2, 1.7])},
        }
        result = default_scalar_comparison(
            analysis_name="test",
            project_name="ANOVA strict alpha",
            metrics_by_condition=metrics_by_condition,
            fdr_alpha=0.001,  # Very strict
        )
        assert result.anova is not None
        # With these overlapping values, p should be > 0.001
        assert result.anova[0].significant is False


def test_anova_boundary_inclusive(monkeypatch) -> None:
    """ANOVA should be significant when p-value equals alpha exactly."""

    class _FakeAnovaResult:
        def __init__(self, f_statistic: float, p_value: float):
            self.f_statistic = f_statistic
            self.p_value = p_value

    def _fake_one_way_anova(*_groups):
        return _FakeAnovaResult(f_statistic=4.2, p_value=0.05)

    monkeypatch.setattr(
        "polyzymd.analyses.shared.inferential_statistics.one_way_anova",
        _fake_one_way_anova,
    )

    metrics = {
        "A": _metric(1.0, [0.9, 1.0, 1.1]),
        "B": _metric(2.0, [1.9, 2.0, 2.1]),
        "C": _metric(3.0, [2.9, 3.0, 3.1]),
    }
    result = anova_test(metrics, metric_name="m", alpha=0.05)

    assert result is not None
    assert result.p_value == pytest.approx(0.05)
    assert result.significant is True


class TestFdrAlphaValidation:
    """fdr_alpha validation edge cases."""

    def test_fdr_alpha_zero_raises(self):
        """fdr_alpha=0 should raise ValueError."""
        metrics = {
            "A": _metric(1.0, [1.0, 1.1, 1.2]),
            "B": _metric(2.0, [2.0, 2.1, 2.2]),
        }
        with pytest.raises(ValueError, match="fdr_alpha must be in"):
            pairwise_comparisons(metrics, fdr_alpha=0.0)

    def test_fdr_alpha_negative_raises(self):
        """fdr_alpha < 0 should raise ValueError."""
        metrics = {
            "A": _metric(1.0, [1.0, 1.1, 1.2]),
            "B": _metric(2.0, [2.0, 2.1, 2.2]),
        }
        with pytest.raises(ValueError, match="fdr_alpha must be in"):
            pairwise_comparisons(metrics, fdr_alpha=-0.05)

    def test_fdr_alpha_nan_raises(self):
        """fdr_alpha=NaN should raise ValueError."""
        metrics = {
            "A": _metric(1.0, [1.0, 1.1, 1.2]),
            "B": _metric(2.0, [2.0, 2.1, 2.2]),
        }
        with pytest.raises(ValueError, match="fdr_alpha must be in"):
            pairwise_comparisons(metrics, fdr_alpha=float("nan"))

    def test_fdr_alpha_one_accepted(self):
        """fdr_alpha=1.0 is valid (everything significant)."""
        metrics = {
            "A": _metric(1.0, [1.0, 1.1, 1.2]),
            "B": _metric(2.0, [2.0, 2.1, 2.2]),
        }
        results = pairwise_comparisons(metrics, fdr_alpha=1.0)
        assert len(results) == 1
        assert results[0].significant is True

    def test_fdr_alpha_greater_than_one_raises(self):
        """fdr_alpha > 1.0 should raise ValueError."""
        metrics = {
            "A": _metric(1.0, [1.0, 1.1, 1.2]),
            "B": _metric(2.0, [2.0, 2.1, 2.2]),
        }
        with pytest.raises(ValueError, match="fdr_alpha must be in"):
            pairwise_comparisons(metrics, fdr_alpha=1.5)
