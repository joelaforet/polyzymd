"""Properties of the pairwise tests of ``ReplicateValues.compare``.

``ReplicateValues.compare`` in ``polyzymd.analyses.timeseries`` runs the
pairwise tests. These tests pin:

1. ``test`` selects the variance assumption. With equal group sizes and
   unequal variances Welch's test gives the same t statistic as Student's but
   fewer degrees of freedom, so the Welch p-value is strictly larger.
2. Every pairwise result carries an adjusted p-value. With a single test in
   the Benjamini-Hochberg family the adjusted p-value equals the raw one.
3. Effect sizes carry Hedges' g and drop the Cohen adjectives when the
   combined sample is smaller than ten replicates.

Group A has a standard deviation of 0.1 and group B a standard deviation of
1.0, so the variances differ by a factor of 100.
"""

from __future__ import annotations

from typing import Any

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
