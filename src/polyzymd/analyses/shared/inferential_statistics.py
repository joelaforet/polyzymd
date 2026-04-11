"""Inferential statistical tests shared across analysis comparisons.

This module provides statistical functions for comparing analysis results
across multiple conditions, including t-tests, ANOVA, and effect sizes.

It is the canonical home for inferential statistics used by analysis plugins
and comparison utilities.

All functions use SciPy for statistical calculations.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike

logger = logging.getLogger("polyzymd.analyses")


@dataclass
class TTestResult:
    """Result of a two-sample t-test.

    Attributes
    ----------
    t_statistic : float
        The t-statistic
    p_value : float
        Two-tailed p-value
    """

    t_statistic: float
    p_value: float

    @property
    def significant(self) -> bool:
        """Whether the result is significant at p < 0.05.

        .. note::
           This uses a hardcoded alpha=0.05 threshold.  The comparison
           pipeline overrides significance with configurable thresholds
           (BH-adjusted or Tukey).  Use this property only as a
           convenience default.
        """
        return self.p_value < 0.05

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "t_statistic": self.t_statistic,
            "p_value": self.p_value,
            "significant": self.significant,
        }


@dataclass
class EffectSize:
    """Cohen's d effect size with interpretation.

    Attributes
    ----------
    cohens_d : float
        The effect size (positive = group1 > group2)
    interpretation : str
        Categorical interpretation: "negligible", "small", "medium", "large"
    direction : str
        "higher" (d > 0), "lower" (d < 0), or "unchanged" (d == 0).
    """

    cohens_d: float
    interpretation: str
    direction: str

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "cohens_d": self.cohens_d,
            "interpretation": self.interpretation,
            "direction": self.direction,
        }


@dataclass
class ANOVAResult:
    """Result of one-way ANOVA.

    Attributes
    ----------
    f_statistic : float
        The F-statistic
    p_value : float
        P-value for the test
    """

    f_statistic: float
    p_value: float

    @property
    def significant(self) -> bool:
        """Whether the result is significant at p < 0.05.

        .. note::
           This uses a hardcoded alpha=0.05 threshold.  The comparison
           pipeline overrides significance with configurable thresholds.
           Use this property only as a convenience default.
        """
        return self.p_value < 0.05

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "f_statistic": self.f_statistic,
            "p_value": self.p_value,
            "significant": self.significant,
        }


@dataclass
class BHResult:
    """Result of Benjamini-Hochberg correction for one hypothesis.

    Attributes
    ----------
    raw_p_value : float | None
        Original uncorrected p-value.
    adjusted_p_value : float | None
        BH-adjusted p-value (q-value). None if raw was None.
    significant : bool
        Whether adjusted_p_value <= alpha.
    rank : int | None
        1-based rank among non-None p-values (smallest=1). None if raw was None.
    """

    raw_p_value: float | None
    adjusted_p_value: float | None
    significant: bool
    rank: int | None


def benjamini_hochberg(
    p_values: Sequence[float | None],
    alpha: float = 0.05,
) -> list[BHResult]:
    """Apply Benjamini-Hochberg FDR correction to a family of p-values.

    Implements the Benjamini-Hochberg (1995) step-up procedure to control
    the false discovery rate. The correction adjusts p-values such that
    declaring significance at ``adjusted_p <= alpha`` controls the expected
    proportion of false discoveries at level *alpha*.

    ``None`` and ``NaN`` entries in *p_values* (e.g. cross-temperature pairs
    where statistics are suppressed, or degenerate tests with undefined
    p-values) are passed through — the corresponding ``BHResult`` has
    ``adjusted_p_value=None`` and ``significant=False``.

    Parameters
    ----------
    p_values : Sequence[float | None]
        Raw two-tailed p-values. ``None`` entries are preserved.
    alpha : float, optional
        FDR significance threshold, by default 0.05.

    Returns
    -------
    list[BHResult]
        One entry per input p-value, in the same order.

    References
    ----------
    Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery
    rate: a practical and powerful approach to multiple testing. *JRSS B*,
    57(1), 289-300.
    """
    if not 0.0 < alpha <= 1.0:
        raise ValueError(f"alpha must satisfy 0 < alpha <= 1, got {alpha}")

    if not p_values:
        return []

    results: list[BHResult] = [
        BHResult(raw_p_value=None, adjusted_p_value=None, significant=False, rank=None)
        for _ in p_values
    ]

    indexed_non_null: list[tuple[int, float]] = []
    for idx, p in enumerate(p_values):
        if p is None:
            continue
        raw_p = float(p)
        if math.isnan(raw_p):
            continue
        indexed_non_null.append((idx, raw_p))

    if not indexed_non_null:
        return results

    indexed_non_null.sort(key=lambda item: item[1])
    m = len(indexed_non_null)

    sorted_p = np.asarray([item[1] for item in indexed_non_null], dtype=np.float64)
    ranks = np.arange(1, m + 1, dtype=np.float64)

    adjusted_sorted = sorted_p * m / ranks
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)

    for rank_idx, ((original_idx, raw_p), adjusted_p) in enumerate(
        zip(indexed_non_null, adjusted_sorted, strict=False),
        start=1,
    ):
        adjusted = float(adjusted_p)
        results[original_idx] = BHResult(
            raw_p_value=raw_p,
            adjusted_p_value=adjusted,
            significant=adjusted <= alpha,
            rank=rank_idx,
        )

    return results


def independent_ttest(
    group1: ArrayLike,
    group2: ArrayLike,
    method: str = "student",
) -> TTestResult:
    """Perform a two-sample independent t-test.

    Tests the null hypothesis that two independent samples have
    identical expected values.

    The ``method`` parameter controls the variance assumption:

    - ``"student"`` uses Student's t-test (``equal_var=True``), which assumes
      equal population variances
    - ``"welch"`` uses Welch's t-test (``equal_var=False``), which does not
      assume equal variances

    Use ``"student"`` when homoscedasticity is a reasonable assumption.
    Use ``"welch"`` when variances may differ across conditions.

    Parameters
    ----------
    group1 : array_like
        First group of values (e.g., control replicate means)
    group2 : array_like
        Second group of values (e.g., treatment replicate means)
    method : str, optional
        T-test method to use: ``"student"`` or ``"welch"``, by default
        ``"student"``.

    Returns
    -------
    TTestResult
        Result containing t-statistic and p-value

    Examples
    --------
    >>> control = [0.715, 0.693, 0.696]  # No polymer RMSF
    >>> treatment = [0.517, 0.586]        # 100% SBMA RMSF
    >>> result = independent_ttest(control, treatment)
    >>> print(f"t = {result.t_statistic:.3f}, p = {result.p_value:.4f}")

    Raises
    ------
    ValueError
        If *method* is not ``"student"`` or ``"welch"``
    """
    from scipy import stats

    g1 = np.asarray(group1, dtype=np.float64)
    g2 = np.asarray(group2, dtype=np.float64)

    # Guard: need at least 2 observations per group for a t-test
    if len(g1) < 2 or len(g2) < 2:
        return TTestResult(t_statistic=float("nan"), p_value=float("nan"))

    if method == "student":
        equal_var = True
    elif method == "welch":
        equal_var = False
    else:
        raise ValueError(f"Unknown t-test method {method!r}; expected 'student' or 'welch'")

    t, p = stats.ttest_ind(g1, g2, equal_var=equal_var)

    return TTestResult(
        t_statistic=float(t),
        p_value=float(p),
    )


def cohens_d(
    group1: ArrayLike,
    group2: ArrayLike,
) -> EffectSize:
    """Compute Cohen's d effect size.

    Cohen's d is the difference between means divided by the pooled
    standard deviation. A positive d means group1 has higher values.

    Parameters
    ----------
    group1 : array_like
        First group (typically control)
    group2 : array_like
        Second group (typically treatment)

    Returns
    -------
    EffectSize
        Effect size with interpretation

    Notes
    -----
    Effect size interpretation (Cohen, 1988):
    - |d| < 0.2: negligible
    - 0.2 <= |d| < 0.5: small
    - 0.5 <= |d| < 0.8: medium
    - |d| >= 0.8: large
    """
    g1 = np.asarray(group1, dtype=np.float64)
    g2 = np.asarray(group2, dtype=np.float64)

    n1, n2 = len(g1), len(g2)

    if n1 < 2 or n2 < 2:
        # Undefined: can't compute pooled std with < 2 samples
        return EffectSize(
            cohens_d=float("nan"),
            interpretation="undefined",
            direction="undetermined",
        )
    else:
        var1 = np.var(g1, ddof=1)
        var2 = np.var(g2, ddof=1)

        # Pooled standard deviation
        pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

        if pooled_std > 0:
            d = float((np.mean(g1) - np.mean(g2)) / pooled_std)
        else:
            # Zero pooled SD with different means is undefined
            if np.mean(g1) != np.mean(g2):
                return EffectSize(
                    cohens_d=float("nan"),
                    interpretation="undefined",
                    direction="undetermined",
                )
            d = 0.0

    # Interpret magnitude
    d_abs = abs(d)
    if d_abs < 0.2:
        interpretation = "negligible"
    elif d_abs < 0.5:
        interpretation = "small"
    elif d_abs < 0.8:
        interpretation = "medium"
    else:
        interpretation = "large"

    # Interpret direction
    if d > 0:
        direction = "higher"
    elif d < 0:
        direction = "lower"
    else:
        direction = "unchanged"

    return EffectSize(
        cohens_d=d,
        interpretation=interpretation,
        direction=direction,
    )


def one_way_anova(*groups: ArrayLike) -> ANOVAResult:
    """Perform classical one-way ANOVA across multiple groups.

    Tests the null hypothesis that all groups have the same mean
    using ``scipy.stats.f_oneway`` (equal variance assumption).

    Parameters
    ----------
    *groups : array_like
        Variable number of groups to compare.  Each group must have
        at least 2 observations; groups with fewer observations cause
        the function to return NaN statistics.

    Returns
    -------
    ANOVAResult
        Result containing F-statistic and p-value.  Both are NaN if
        any group has fewer than 2 observations.

    Examples
    --------
    >>> no_poly = [0.715, 0.693, 0.696]
    >>> sbma = [0.517, 0.586]
    >>> egma = [0.558, 0.738, 0.496]
    >>> result = one_way_anova(no_poly, sbma, egma)
    >>> print(f"F = {result.f_statistic:.3f}, p = {result.p_value:.4f}")
    """
    from scipy import stats

    arrays = [np.asarray(g, dtype=np.float64) for g in groups]

    # Guard: need at least 2 groups for ANOVA
    if len(arrays) < 2:
        return ANOVAResult(f_statistic=float("nan"), p_value=float("nan"))

    # Guard: need at least 2 observations per group for ANOVA
    if any(len(a) < 2 for a in arrays):
        return ANOVAResult(f_statistic=float("nan"), p_value=float("nan"))

    stat, p = stats.f_oneway(*arrays)

    return ANOVAResult(f_statistic=float(stat), p_value=float(p))


@dataclass
class TukeyHSDResult:
    """Result of Tukey's HSD test for one pair of groups.

    Attributes
    ----------
    group_i : int
        Index of the first group.
    group_j : int
        Index of the second group.
    statistic : float
        Mean difference (group_j - group_i).
    p_value : float
        Tukey-adjusted p-value for this pair.
    """

    group_i: int
    group_j: int
    statistic: float
    p_value: float


def tukey_hsd(*groups: ArrayLike) -> list[TukeyHSDResult]:
    """Run Tukey's Honestly Significant Difference test.

    Computes family-wise-adjusted p-values for all pairwise group
    comparisons using ``scipy.stats.tukey_hsd``.

    Parameters
    ----------
    *groups : array_like
        Variable number of groups to compare.  Each group must have
        at least 2 observations.

    Returns
    -------
    list[TukeyHSDResult]
        One result per unique pair (i < j), ordered by (i, j).
        Returns an empty list if fewer than 2 groups are provided or any
        group has fewer than 2 observations.

    Examples
    --------
    >>> results = tukey_hsd([1, 2, 3], [4, 5, 6], [7, 8, 9])
    >>> for r in results:
    ...     print(f"({r.group_i}, {r.group_j}): p={r.p_value:.4f}")
    """
    from scipy import stats

    arrays = [np.asarray(g, dtype=np.float64) for g in groups]

    if len(arrays) < 2:
        return []
    if any(len(a) < 2 for a in arrays):
        logger.warning(
            "Tukey HSD requires at least 2 observations per group; "
            "got sizes %s — returning empty results",
            [len(a) for a in arrays],
        )
        return []

    result = stats.tukey_hsd(*arrays)
    pairs: list[TukeyHSDResult] = []
    n = len(arrays)
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append(
                TukeyHSDResult(
                    group_i=i,
                    group_j=j,
                    statistic=float(result.statistic[j, i]),
                    p_value=float(result.pvalue[i, j]),
                )
            )
    return pairs


def percent_change(control_mean: float, treatment_mean: float) -> float:
    """Calculate percent change from control.

    Parameters
    ----------
    control_mean : float
        Mean value of control condition
    treatment_mean : float
        Mean value of treatment condition

    Returns
    -------
    float
        Percent change: (treatment - control) / control * 100
        Negative = reduction, Positive = increase.

        Special handling for zero control values:

        - 0 -> 0 returns ``0.0``
        - 0 -> positive returns ``math.inf``
        - 0 -> negative returns ``-math.inf``

        If either input is non-finite (NaN or +/-inf), returns ``math.nan``.
    """
    if not (math.isfinite(control_mean) and math.isfinite(treatment_mean)):
        return math.nan

    if control_mean == 0:
        if treatment_mean == 0:
            return 0.0
        return math.inf if treatment_mean > 0 else -math.inf

    return (treatment_mean - control_mean) / control_mean * 100
