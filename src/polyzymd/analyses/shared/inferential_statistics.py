"""Inferential statistical tests shared across analysis comparisons.

This module provides statistical functions for comparing analysis results
across multiple conditions, including t-tests, ANOVA, and effect sizes.

It is the canonical home for inferential statistics used by analysis plugins
and comparison utilities.

All functions use SciPy for statistical calculations.

Correction family
-----------------
:func:`apply_family_correction` defines the multiple-comparison family once
for the whole package. One analysis run is one family: every pairwise test
the run produced, across all of its metrics and all of its condition pairs,
is corrected together with Benjamini-Hochberg.

One-way ANOVA sits outside that family. Its p-value is reported raw and is
labelled an omnibus test; it is never adjusted and it never gates the
pairwise tests. Significance for an ANOVA is the raw p-value compared with
the same alpha the pairwise family uses.

Direction labels ("stabilizing", "increased", and the rest) describe a
finding, so they are only assigned when the corrected test is significant.
Everything else is labelled ``"no significant change"``.

References
----------
Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate:
    a practical and powerful approach to multiple testing. *Journal of the
    Royal Statistical Society B*, 57(1), 289-300.
    doi:10.1111/j.2517-6161.1995.tb02031.x
Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences*,
    2nd edition. Lawrence Erlbaum Associates.
Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
    size and related estimators. *Journal of Educational Statistics*, 6(2),
    107-128. doi:10.3102/10769986006002107
Tukey, J. W. (1949). Comparing individual means in the analysis of variance.
    *Biometrics*, 5(2), 99-114. doi:10.2307/3001913
Welch, B. L. (1947). The generalization of "Student's" problem when several
    different population variances are involved. *Biometrika*, 34(1-2),
    28-35. doi:10.1093/biomet/34.1-2.28
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike

logger = logging.getLogger("polyzymd.analyses")

NO_SIGNIFICANT_CHANGE = "no significant change"
"""Direction label used when a comparison is not significant."""

MIN_N_FOR_EFFECT_SIZE_LABEL = 10
"""Smallest combined sample size that earns a Cohen (1988) effect adjective.

Below this the standard error of the effect size is of order one, so the
adjective would report noise. :func:`cohens_d` returns ``None`` for the
interpretation instead.
"""


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
    """Standardized mean difference with an optional interpretation.

    Attributes
    ----------
    cohens_d : float
        The effect size (positive = group1 > group2).
    interpretation : str or None
        Categorical interpretation: "negligible", "small", "medium" or
        "large". ``None`` when the combined sample is smaller than
        :data:`MIN_N_FOR_EFFECT_SIZE_LABEL`, or when the effect size is
        undefined.
    direction : str
        "higher" (d > 0), "lower" (d < 0), or "unchanged" (d == 0).
    hedges_g : float
        Cohen's d multiplied by the Hedges (1981) small-sample correction
        ``J = 1 - 3 / (4 * (n1 + n2) - 9)``. Always report this value for
        replicate counts of the size molecular dynamics produces.
    """

    cohens_d: float
    interpretation: str | None
    direction: str
    hedges_g: float = float("nan")

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "cohens_d": self.cohens_d,
            "hedges_g": self.hedges_g,
            "interpretation": self.interpretation,
            "direction": self.direction,
        }


def hedges_correction(n1: int, n2: int) -> float:
    """Return the Hedges (1981) bias correction factor J.

    Parameters
    ----------
    n1, n2 : int
        Sample sizes of the two groups.

    Returns
    -------
    float
        ``J = 1 - 3 / (4 * (n1 + n2) - 9)``. Returns NaN when the
        denominator is not positive.

    References
    ----------
    Hedges, L. V. (1981). Distribution theory for Glass's estimator of
    effect size and related estimators. *Journal of Educational
    Statistics*, 6(2), 107-128. doi:10.3102/10769986006002107
    """
    denominator = 4.0 * (n1 + n2) - 9.0
    if denominator <= 0:
        return float("nan")
    return 1.0 - 3.0 / denominator


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
        Effect size carrying both Cohen's d and Hedges' g.

    Notes
    -----
    Effect size interpretation (Cohen, 1988):

    - |d| < 0.2: negligible
    - 0.2 <= |d| < 0.5: small
    - 0.5 <= |d| < 0.8: medium
    - |d| >= 0.8: large

    The adjective is withheld (``interpretation`` is ``None``) when
    ``n1 + n2`` is below :data:`MIN_N_FOR_EFFECT_SIZE_LABEL`. At three
    replicates per condition the standard error of d is of order one, so
    the boundary between "medium" and "large" carries no information.
    ``hedges_g`` applies the Hedges (1981) correction J and is the value
    to quote at these sample sizes.

    References
    ----------
    Cohen, J. (1988). *Statistical Power Analysis for the Behavioral
    Sciences*, 2nd edition. Lawrence Erlbaum Associates.

    Hedges, L. V. (1981). Distribution theory for Glass's estimator of
    effect size and related estimators. *Journal of Educational
    Statistics*, 6(2), 107-128. doi:10.3102/10769986006002107
    """
    g1 = np.asarray(group1, dtype=np.float64)
    g2 = np.asarray(group2, dtype=np.float64)

    n1, n2 = len(g1), len(g2)
    undefined = EffectSize(
        cohens_d=float("nan"),
        interpretation=None,
        direction="undetermined",
        hedges_g=float("nan"),
    )

    if n1 < 2 or n2 < 2:
        # Undefined: can't compute pooled std with < 2 samples
        return undefined

    var1 = np.var(g1, ddof=1)
    var2 = np.var(g2, ddof=1)

    # Pooled standard deviation
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

    if pooled_std > 0:
        d = float((np.mean(g1) - np.mean(g2)) / pooled_std)
    elif np.mean(g1) != np.mean(g2):
        # Zero pooled SD with different means is undefined
        return undefined
    else:
        d = 0.0

    g = d * hedges_correction(n1, n2)

    # Interpret magnitude, but only when the sample can support an adjective
    interpretation: str | None
    if n1 + n2 < MIN_N_FOR_EFFECT_SIZE_LABEL:
        interpretation = None
    else:
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
        hedges_g=g,
    )


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


def enforce_direction_significance(
    results: Sequence[object],
    fields: Sequence[tuple[str, str]] = (("direction", "significant"),),
) -> None:
    """Replace direction labels on results that are not significant.

    A label such as "stabilizing" or "increased" is a claim about the
    system. Without a significant test there is nothing to claim, so the
    label becomes :data:`NO_SIGNIFICANT_CHANGE`.

    Parameters
    ----------
    results : Sequence[object]
        Result objects to relabel in place.
    fields : Sequence[tuple[str, str]], optional
        ``(direction_attribute, significance_attribute)`` pairs to check,
        by default the single pair ``("direction", "significant")``.
    """
    for result in results:
        for direction_attr, significant_attr in fields:
            if getattr(result, direction_attr, None) is None:
                continue
            if not getattr(result, significant_attr, False):
                setattr(result, direction_attr, NO_SIGNIFICANT_CHANGE)
