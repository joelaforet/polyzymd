"""Tests for shared multi-run comparison helpers."""

from __future__ import annotations

from types import SimpleNamespace

from polyzymd.analyses.shared.multi_run_comparison import (
    apply_fdr_correction,
    build_condition_pairs,
    filter_summaries_with_run,
)


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
