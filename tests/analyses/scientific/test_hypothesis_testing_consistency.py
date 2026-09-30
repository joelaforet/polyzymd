"""Hypothesis testing must behave the same way on every comparison path.

``ReplicateValues.compare`` in ``polyzymd.analyses.timeseries`` and
``default_scalar_comparison`` in ``polyzymd.analyses.stats`` run the pairwise
tests. These tests pin the properties that must not depend on which one was
asked:

1. ``test`` selects the variance assumption. With equal group sizes and
   unequal variances Welch's test gives the same t statistic as Student's but
   fewer degrees of freedom, so the Welch p-value is strictly larger.
2. Every pairwise result carries an adjusted p-value. With a single test in
   the Benjamini-Hochberg family the adjusted p-value equals the raw one, and
   all pairs and metrics of one comparison form one family.
3. Effect sizes carry Hedges' g and drop the Cohen adjectives when the
   combined sample is smaller than ten replicates.

Group A has a standard deviation of 0.1 and group B a standard deviation of
1.0, so the variances differ by a factor of 100.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

LOW_VARIANCE = (10.0, 10.1, 9.9)
HIGH_VARIANCE = (12.0, 13.0, 11.0)


def test_hedges_g_corrects_cohens_d_downward() -> None:
    """Hedges' g must apply the 1981 bias correction J to Cohen's d."""
    from polyzymd.analyses.shared.inferential_statistics import cohens_d

    effect = cohens_d(LOW_VARIANCE, HIGH_VARIANCE)

    n1 = n2 = 3
    expected_j = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    assert effect.hedges_g == pytest.approx(effect.cohens_d * expected_j)
    assert abs(effect.hedges_g) < abs(effect.cohens_d)


def test_effect_size_adjectives_dropped_for_small_samples() -> None:
    """Cohen's adjectives are noise below ten replicates, so they are withheld."""
    from polyzymd.analyses.shared.inferential_statistics import cohens_d

    small = cohens_d([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    assert small.interpretation is None

    large_group_a = [float(value) for value in range(6)]
    large_group_b = [float(value) + 10.0 for value in range(6)]
    large = cohens_d(large_group_a, large_group_b)
    assert large.interpretation == "large"


def test_public_functions_have_resolvable_annotations() -> None:
    """Every public function in the statistics module must be introspectable.

    Ruff does not flag an undefined name in an annotation here because F821
    is suppressed for this project, so a missing import only shows up when
    something resolves the annotations.
    """
    import inspect
    import typing

    from polyzymd.analyses.shared import inferential_statistics

    functions = [
        obj
        for name, obj in vars(inferential_statistics).items()
        if not name.startswith("_")
        and inspect.isfunction(obj)
        and obj.__module__ == inferential_statistics.__name__
    ]

    assert functions
    for function in functions:
        typing.get_type_hints(function)


def test_correction_family_spans_every_metric(monkeypatch: pytest.MonkeyPatch) -> None:
    """Three conditions and two metrics make one family of six tests."""
    from polyzymd.analyses.base import MetricValue
    from polyzymd.analyses.shared import inferential_statistics
    from polyzymd.analyses.stats import default_scalar_comparison

    family_sizes: list[int] = []
    real_benjamini_hochberg = inferential_statistics.benjamini_hochberg

    def _recording_benjamini_hochberg(p_values, alpha=0.05):
        family_sizes.append(len(p_values))
        return real_benjamini_hochberg(p_values, alpha=alpha)

    monkeypatch.setattr(inferential_statistics, "benjamini_hochberg", _recording_benjamini_hochberg)

    def _metrics(offset: float) -> dict[str, MetricValue]:
        first = [1.0 + offset, 1.1 + offset, 0.9 + offset]
        second = [5.0 + offset, 5.2 + offset, 4.8 + offset]
        return {
            "first": MetricValue(
                name="first",
                mean=float(np.mean(first)),
                sem=0.05,
                replicate_values=first,
            ),
            "second": MetricValue(
                name="second",
                mean=float(np.mean(second)),
                sem=0.05,
                replicate_values=second,
            ),
        }

    result = default_scalar_comparison(
        analysis_name="family",
        project_name="family",
        metrics_by_condition={
            "A": _metrics(0.0),
            "B": _metrics(1.0),
            "C": _metrics(2.0),
        },
        control_label=None,
    )

    # Three conditions give three pairs, and both metrics join the same family.
    assert len(result.pairwise_comparisons) == 6
    assert family_sizes == [6]
    assert all(comp.p_value_adjusted is not None for comp in result.pairwise_comparisons)

    # The omnibus ANOVA is outside the family and stays uncorrected.
    assert result.anova is not None
    assert len(result.anova) == 2
    assert all(anova.p_value_adjusted is None for anova in result.anova)


def _function_path(test: str) -> Any:
    """Compare the two groups through ReplicateValues.compare."""
    from tests._support.analysis_testkit import replicate_values

    values = replicate_values({"Control": list(LOW_VARIANCE), "Treated": list(HIGH_VARIANCE)})
    return values.compare(test=test)


def test_function_path_welch_gives_larger_p_than_student() -> None:
    """Welch's test widens the p-value on the function path too."""
    (welch,) = _function_path("welch").pairwise
    (student,) = _function_path("student").pairwise

    assert welch.p > student.p
    assert welch.p != pytest.approx(student.p)


def test_function_path_single_comparison_leaves_p_value_unchanged() -> None:
    """One comparison is a family of one, so its adjusted p-value is the raw one."""
    (row,) = _function_path("welch").pairwise

    assert row.p_adjusted is not None
    assert row.p_adjusted == pytest.approx(row.p)


def test_function_path_direction_requires_significance() -> None:
    """A comparison that is not significant claims no direction."""
    from polyzymd.analyses.shared.inferential_statistics import NO_SIGNIFICANT_CHANGE

    (row,) = _function_path("welch").pairwise

    assert not row.significant
    assert row.direction == NO_SIGNIFICANT_CHANGE


def test_function_path_reports_hedges_g() -> None:
    """Hedges' g is Cohen's d times J and points the same way as the difference."""
    (row,) = _function_path("welch").pairwise

    expected_j = 1.0 - 3.0 / (4.0 * 6 - 9.0)
    assert row.hedges_g == pytest.approx(row.cohens_d * expected_j)
    assert row.delta > 0 and row.cohens_d > 0
