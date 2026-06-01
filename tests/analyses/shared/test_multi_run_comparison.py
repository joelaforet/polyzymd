"""Tests for shared multi-run comparison helpers."""

from __future__ import annotations

import math
from datetime import datetime
from types import SimpleNamespace

import pytest

from polyzymd.analyses.base import (
    BaseComparisonResult,
    BaseConditionSummary,
    PairwiseResult,
)
from polyzymd.analyses.shared.multi_run_comparison import (
    apply_fdr_correction,
    build_condition_pairs,
    filter_summaries_with_run,
)


class _ConditionSummary(BaseConditionSummary):
    @property
    def primary_metric_value(self) -> float:
        return 0.0

    @property
    def primary_metric_sem(self) -> float:
        return 0.0


class _ComparisonResult(BaseComparisonResult[_ConditionSummary, PairwiseResult]):
    pass


def test_base_comparison_get_comparison_prefers_pair_lookup() -> None:
    """Tuple key lookup should disambiguate all-vs-all pairwise comparisons."""
    result = _ComparisonResult(
        metric="mean_value",
        name="test_project",
        control_label=None,
        conditions=[
            _ConditionSummary(
                label="A",
                config_path="/fake/a.yaml",
                n_replicates=2,
                replicate_values=[1.0, 1.1],
            ),
            _ConditionSummary(
                label="B",
                config_path="/fake/b.yaml",
                n_replicates=2,
                replicate_values=[1.2, 1.3],
            ),
            _ConditionSummary(
                label="C",
                config_path="/fake/c.yaml",
                n_replicates=2,
                replicate_values=[1.4, 1.5],
            ),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="A",
                condition_b="C",
                t_statistic=1.0,
                p_value=0.20,
                cohens_d=0.2,
                effect_size_interpretation="small",
                direction="increased",
                significant=False,
                percent_change=5.0,
            ),
            PairwiseResult(
                condition_a="B",
                condition_b="C",
                t_statistic=2.0,
                p_value=0.01,
                cohens_d=0.8,
                effect_size_interpretation="large",
                direction="decreased",
                significant=True,
                percent_change=-5.0,
            ),
        ],
        ranking=["B", "A", "C"],
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )

    exact = result.get_comparison(("B", "C"))
    assert exact is not None
    assert exact.condition_a == "B"
    assert exact.condition_b == "C"


def test_base_comparison_get_comparison_supports_noncanonical_condition_b_lookup() -> None:
    """String lookup should remain available for backward compatibility."""
    result = _ComparisonResult(
        metric="mean_value",
        name="test_project",
        conditions=[
            _ConditionSummary(
                label="Control",
                config_path="/fake/control.yaml",
                n_replicates=2,
                replicate_values=[1.0, 1.1],
            ),
            _ConditionSummary(
                label="Treatment",
                config_path="/fake/treatment.yaml",
                n_replicates=2,
                replicate_values=[0.9, 0.8],
            ),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                t_statistic=3.0,
                p_value=0.02,
                cohens_d=1.0,
                effect_size_interpretation="large",
                direction="decreased",
                significant=True,
                percent_change=-10.0,
            )
        ],
        ranking=["Treatment", "Control"],
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )

    noncanonical = result.get_comparison("Treatment")
    assert noncanonical is not None
    assert noncanonical.condition_a == "Control"
    assert noncanonical.condition_b == "Treatment"


def test_base_comparison_get_comparison_noncanonical_lookup_raises_on_ambiguity() -> None:
    """Non-canonical single-label lookup should raise when multiple matches exist."""
    result = _ComparisonResult(
        metric="mean_value",
        name="test_project",
        control_label=None,
        conditions=[
            _ConditionSummary(
                label="A",
                config_path="/fake/a.yaml",
                n_replicates=2,
                replicate_values=[1.0, 1.1],
            ),
            _ConditionSummary(
                label="B",
                config_path="/fake/b.yaml",
                n_replicates=2,
                replicate_values=[1.2, 1.3],
            ),
            _ConditionSummary(
                label="C",
                config_path="/fake/c.yaml",
                n_replicates=2,
                replicate_values=[1.4, 1.5],
            ),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="A",
                condition_b="C",
                t_statistic=1.0,
                p_value=0.20,
                cohens_d=0.2,
                effect_size_interpretation="small",
                direction="increased",
                significant=False,
                percent_change=5.0,
            ),
            PairwiseResult(
                condition_a="B",
                condition_b="C",
                t_statistic=2.0,
                p_value=0.01,
                cohens_d=0.8,
                effect_size_interpretation="large",
                direction="decreased",
                significant=True,
                percent_change=-5.0,
            ),
        ],
        ranking=["B", "A", "C"],
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )

    with pytest.raises(ValueError, match="Ambiguous non-canonical comparison lookup"):
        result.get_comparison("C")


def test_base_comparison_get_comparison_tuple_falls_back_to_condition_b_lookup() -> None:
    """Tuple lookup should fall back to non-canonical condition_b behavior when needed."""
    result = _ComparisonResult(
        metric="mean_value",
        name="test_project",
        conditions=[
            _ConditionSummary(
                label="Control",
                config_path="/fake/control.yaml",
                n_replicates=2,
                replicate_values=[1.0, 1.1],
            ),
            _ConditionSummary(
                label="Treatment",
                config_path="/fake/treatment.yaml",
                n_replicates=2,
                replicate_values=[0.9, 0.8],
            ),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                t_statistic=3.0,
                p_value=0.02,
                cohens_d=1.0,
                effect_size_interpretation="large",
                direction="decreased",
                significant=True,
                percent_change=-10.0,
            )
        ],
        ranking=["Treatment", "Control"],
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )

    fallback = result.get_comparison(("NotControl", "Treatment"))
    assert fallback is not None
    assert fallback.condition_a == "Control"
    assert fallback.condition_b == "Treatment"


def test_filter_summaries_with_run_keeps_only_matching_entries() -> None:
    """filter_summaries_with_run should keep only summaries with a given run."""

    class _Summary:
        def __init__(self, runs: set[str]) -> None:
            self._runs = runs

        def get_run(self, run_label: str) -> str:
            if run_label not in self._runs:
                raise KeyError(run_label)
            return run_label

    summaries = {
        "A": _Summary({"run_1", "run_2"}),
        "B": _Summary({"run_2"}),
        "C": _Summary({"run_1"}),
    }

    filtered = filter_summaries_with_run(
        summaries, "run_1", lambda summary, label: summary.get_run(label)
    )

    assert list(filtered.keys()) == ["A", "C"]


def test_build_condition_pairs_with_control() -> None:
    """build_condition_pairs should prefer control-vs-treatment pairs."""
    pairs = build_condition_pairs(["Control", "A", "B"], control_label="Control")
    assert pairs == [("Control", "A"), ("Control", "B")]


def test_build_condition_pairs_control_missing_falls_back_all_pairs() -> None:
    """build_condition_pairs should support all-pairs fallback on missing control."""
    pairs = build_condition_pairs(
        ["A", "B", "C"], control_label="Control", on_control_missing="all_pairs"
    )
    assert pairs == [("A", "B"), ("A", "C"), ("B", "C")]


def test_build_condition_pairs_control_missing_skip() -> None:
    """build_condition_pairs should support skip mode when control is absent."""
    pairs = build_condition_pairs(
        ["A", "B", "C"], control_label="Control", on_control_missing="skip"
    )
    assert pairs == []


def test_build_condition_pairs_rejects_invalid_control_missing_policy() -> None:
    """build_condition_pairs should reject unsupported control-missing policies."""
    with pytest.raises(ValueError, match="on_control_missing"):
        build_condition_pairs(
            ["A", "B", "C"],
            control_label="Control",
            on_control_missing="fallback",
        )


def test_apply_fdr_correction_updates_pairwise_and_anova() -> None:
    """apply_fdr_correction should write adjusted p-values and significance."""
    pairwise = [
        SimpleNamespace(p_value=0.01, p_value_adjusted=None, significant=False),
        SimpleNamespace(p_value=0.04, p_value_adjusted=None, significant=False),
    ]
    anova = [
        SimpleNamespace(p_value=0.02, p_value_adjusted=None, significant=False),
        SimpleNamespace(p_value=0.20, p_value_adjusted=None, significant=False),
    ]

    apply_fdr_correction(pairwise, anova_by_run=anova, fdr_alpha=0.05)

    assert pairwise[0].p_value_adjusted == 0.02
    assert pairwise[1].p_value_adjusted == 0.04
    assert pairwise[0].significant is True
    assert pairwise[1].significant is True
    assert anova[0].p_value_adjusted == 0.04
    assert anova[1].p_value_adjusted == 0.20
    assert anova[0].significant is True
    assert anova[1].significant is False


def test_apply_fdr_correction_supports_custom_callbacks() -> None:
    """apply_fdr_correction should support custom p-value extraction and setters."""
    pairwise = [
        SimpleNamespace(raw_p=0.03, adjusted=None, keep=False),
        SimpleNamespace(raw_p=0.50, adjusted=None, keep=False),
    ]

    def _get_p_value(result: SimpleNamespace) -> float:
        return result.raw_p

    def _set_corrected(result: SimpleNamespace, bh_result: SimpleNamespace) -> None:
        result.adjusted = bh_result.adjusted_p_value
        result.keep = bh_result.significant

    apply_fdr_correction(
        pairwise,
        anova_by_run=None,
        fdr_alpha=0.05,
        get_p_value=_get_p_value,
        set_corrected=_set_corrected,
    )

    assert pairwise[0].adjusted == 0.06
    assert pairwise[0].keep is False
    assert pairwise[1].adjusted == 0.50
    assert pairwise[1].keep is False


def test_apply_fdr_correction_handles_none_and_nan_p_values() -> None:
    """apply_fdr_correction should preserve None/NaN entries as non-significant."""
    pairwise = [
        SimpleNamespace(p_value=0.01, p_value_adjusted=None, significant=False),
        SimpleNamespace(p_value=None, p_value_adjusted=1.0, significant=True),
        SimpleNamespace(p_value=math.nan, p_value_adjusted=1.0, significant=True),
    ]

    apply_fdr_correction(pairwise, anova_by_run=None, fdr_alpha=0.05)

    assert pairwise[0].p_value_adjusted == 0.01
    assert pairwise[0].significant is True
    assert pairwise[1].p_value_adjusted is None
    assert pairwise[1].significant is False
    assert pairwise[2].p_value_adjusted is None
    assert pairwise[2].significant is False


def test_apply_fdr_correction_mixed_testable_and_non_testable_entries() -> None:
    """apply_fdr_correction should update only entries with valid p-values."""
    pairwise = [
        SimpleNamespace(
            p_value=0.02,
            p_value_adjusted=None,
            significant=False,
            testable=True,
        ),
        SimpleNamespace(
            p_value=None,
            p_value_adjusted=0.99,
            significant=True,
            testable=False,
        ),
        SimpleNamespace(
            p_value=0.20,
            p_value_adjusted=None,
            significant=False,
            testable=True,
        ),
    ]

    apply_fdr_correction(pairwise, anova_by_run=None, fdr_alpha=0.05)

    assert pairwise[0].p_value_adjusted == 0.04
    assert pairwise[0].significant is True
    assert pairwise[1].p_value_adjusted is None
    assert pairwise[1].significant is False
    assert pairwise[2].p_value_adjusted == 0.20
    assert pairwise[2].significant is False
