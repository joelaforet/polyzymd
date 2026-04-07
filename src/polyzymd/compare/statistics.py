"""Statistical tests for comparing simulation conditions.

This module provides statistical functions for comparing analysis results
across multiple conditions, including t-tests, ANOVA, and effect sizes.

All functions use SciPy for statistical calculations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike


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
        """Whether the result is significant at p < 0.05."""
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
        "higher" (d > 0) or "lower" (d < 0).
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
        """Whether the result is significant at p < 0.05."""
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
) -> TTestResult:
    """Perform Welch's two-sample independent t-test.

    Tests the null hypothesis that two independent samples have
    identical expected values.  Uses Welch's t-test (``equal_var=False``)
    which does **not** assume equal population variances.  This is the
    appropriate default for MD replicate data where the number of
    replicates and variance may differ between conditions.

    Parameters
    ----------
    group1 : array_like
        First group of values (e.g., control replicate means)
    group2 : array_like
        Second group of values (e.g., treatment replicate means)

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
    """
    from scipy import stats

    g1 = np.asarray(group1, dtype=np.float64)
    g2 = np.asarray(group2, dtype=np.float64)

    t, p = stats.ttest_ind(g1, g2, equal_var=False)

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
        # Can't compute pooled std with < 2 samples
        d = 0.0
    else:
        var1 = np.var(g1, ddof=1)
        var2 = np.var(g2, ddof=1)

        # Pooled standard deviation
        pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

        if pooled_std > 0:
            d = float((np.mean(g1) - np.mean(g2)) / pooled_std)
        else:
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
    direction = "higher" if d > 0 else "lower"

    return EffectSize(
        cohens_d=d,
        interpretation=interpretation,
        direction=direction,
    )


def one_way_anova(*groups: ArrayLike) -> ANOVAResult:
    """Perform one-way ANOVA across multiple groups.

    Tests the null hypothesis that all groups have the same mean.

    Parameters
    ----------
    *groups : array_like
        Variable number of groups to compare

    Returns
    -------
    ANOVAResult
        Result containing F-statistic and p-value

    Examples
    --------
    >>> no_poly = [0.715, 0.693, 0.696]
    >>> sbma = [0.517, 0.586]
    >>> egma = [0.558, 0.738, 0.496]
    >>> result = one_way_anova(no_poly, sbma, egma)
    >>> print(f"F = {result.f_statistic:.3f}, p = {result.p_value:.4f}")
    """
    from scipy import stats

    # Convert to numpy arrays
    arrays = [np.asarray(g, dtype=np.float64) for g in groups]

    f, p = stats.f_oneway(*arrays)

    return ANOVAResult(
        f_statistic=float(f),
        p_value=float(p),
    )


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
