"""Known-answer tests for how the analyses package reports uncertainty.

Grossfield et al. (2018, LiveCoMS 1:5067) require that a simulation study
report 95 percent confidence intervals rather than bare standard errors, that
the coverage factor come from the Student t distribution with n - 1 degrees of
freedom, and that every figure describe the meaning and basis of its
uncertainties. These tests pin those three rules and the rule that a single
replicate has no estimable uncertainty at all.
"""

from __future__ import annotations

import doctest
import math

import pytest

from polyzymd.analyses.exceptions import AnalysisError, StatisticsError
from polyzymd.analyses.shared.statistics import mean_sem_ci, student_t_coverage_factor

# Two-sided Student t coverage factors at 95 percent, quoted in the article.
T_FACTOR_N3 = 4.302652729749462
T_FACTOR_N5 = 2.7764451051977934


class TestCoverageFactor:
    """The interval must use the Student t factor, not 1.96."""

    @pytest.mark.parametrize(
        ("n", "expected"),
        [(3, T_FACTOR_N3), (5, T_FACTOR_N5)],
    )
    def test_coverage_factor_matches_published_values(self, n: int, expected: float) -> None:
        """Coverage factors should match the article's table."""

        assert student_t_coverage_factor(n) == pytest.approx(expected, rel=1e-9)

    def test_coverage_factor_is_none_for_one_replicate(self) -> None:
        """One replicate has no degrees of freedom and so no coverage factor."""

        assert student_t_coverage_factor(1) is None

    def test_normal_factor_is_rejected_as_a_shortcut(self) -> None:
        """At n = 3 the factor must not collapse to the normal 1.96."""

        assert student_t_coverage_factor(3) > 4.0


class TestMeanSemCI:
    """The one interval estimator in the package."""

    def test_three_values_give_mean_plus_or_minus_t_times_sem(self) -> None:
        """A three-replicate interval is mean +/- 4.303 SEM."""

        values = [2.0, 2.2, 2.4]
        result = mean_sem_ci(values)

        expected_sem = 0.2 / math.sqrt(3.0)
        assert result.n == 3
        assert result.mean == pytest.approx(2.2)
        assert result.sem == pytest.approx(expected_sem)
        assert result.ci_low == pytest.approx(2.2 - T_FACTOR_N3 * expected_sem)
        assert result.ci_high == pytest.approx(2.2 + T_FACTOR_N3 * expected_sem)
        assert result.ci_method == "student_t"
        assert result.coverage == pytest.approx(0.95)

    def test_interval_is_wider_than_the_sem_band(self) -> None:
        """The interval must be materially wider than plus or minus one SEM."""

        result = mean_sem_ci([2.0, 2.2, 2.4])

        half_width = result.ci_high - result.mean
        assert half_width > 4.0 * result.sem

    def test_single_value_has_no_uncertainty(self) -> None:
        """One replicate yields None, never 0.0."""

        result = mean_sem_ci([1.5])

        assert result.n == 1
        assert result.mean == pytest.approx(1.5)
        assert result.sem is None
        assert result.ci_low is None
        assert result.ci_high is None
        assert result.ci_method is None

    def test_empty_input_raises(self) -> None:
        """An empty sample is an error, not a zero."""

        with pytest.raises(ValueError, match="empty array"):
            mean_sem_ci([])


class TestStatisticsErrorsAreTyped:
    """Invalid statistical input raises a typed analysis error."""

    def test_empty_sample(self) -> None:
        """An empty sample is a typed failure."""

        with pytest.raises(StatisticsError):
            mean_sem_ci([])

    def test_bad_coverage(self) -> None:
        """A coverage outside (0, 1) is a typed failure."""

        with pytest.raises(StatisticsError):
            student_t_coverage_factor(3, coverage=1.5)

    def test_errors_are_analysis_errors(self) -> None:
        """The typed error belongs to the analyses exception hierarchy."""

        assert issubclass(StatisticsError, AnalysisError)


class TestDoctests:
    """The worked example in the statistics module must stay correct."""

    def test_statistics_module_doctests(self) -> None:
        """Run the doctests in polyzymd.analyses.shared.statistics."""

        import polyzymd.analyses.shared.statistics as statistics_module

        results = doctest.testmod(statistics_module, verbose=False)

        assert results.failed == 0, f"{results.failed} doctest failures"


class TestFunctionPathUncertainty:
    """ReplicateValues.summary reports the Student t interval of the replicate values."""

    def test_three_replicates_use_the_student_t_factor(self) -> None:
        """The half-width is 4.30 standard errors at n = 3, with the unit stated."""

        from tests._support.analysis_testkit import replicate_values

        report = replicate_values({"A": [2.0, 2.2, 2.4]}).summary()
        (row,) = report.conditions

        assert report.unit == "A"
        assert row.ci_method == "student_t"
        assert row.ci95[1] - row.mean == pytest.approx(T_FACTOR_N3 * row.sem, rel=1e-12)
        assert row.sem == pytest.approx(0.2 / math.sqrt(3.0))

    def test_one_replicate_has_no_uncertainty(self) -> None:
        """A single replicate reports its value with no SEM and no interval."""

        from tests._support.analysis_testkit import replicate_values

        report = replicate_values({"Solo": [12.0]}).summary()
        (row,) = report.conditions

        assert (row.n_replicates, row.mean) == (1, 12.0)
        assert row.sem is None and row.ci95 is None and row.ci_method is None
        assert any("one replicate" in text for text in report.warnings)
