"""Statistical tests for comparing simulation conditions.

This module provides statistical functions for comparing analysis results
across multiple conditions, including t-tests, ANOVA, and effect sizes.

All functions use SciPy for statistical calculations.
"""

from __future__ import annotations

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
        For RMSF: "stabilizing" (d > 0, lower RMSF) or "destabilizing" (d < 0)
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


def independent_ttest(
    group1: ArrayLike,
    group2: ArrayLike,
) -> TTestResult:
    """Perform two-sample independent t-test.

    Tests the null hypothesis that two independent samples have
    identical expected values.

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

    t, p = stats.ttest_ind(g1, g2)

    return TTestResult(
        t_statistic=float(t),
        p_value=float(p),
    )


def cohens_d(
    group1: ArrayLike,
    group2: ArrayLike,
    rmsf_mode: bool = True,
) -> EffectSize:
    """Compute Cohen's d effect size.

    Cohen's d is the difference between means divided by the pooled
    standard deviation. A positive d means group1 has higher values.

    For RMSF comparisons (rmsf_mode=True), direction is interpreted as:
    - d > 0 (control > treatment) = "stabilizing" (treatment reduces RMSF)
    - d < 0 (control < treatment) = "destabilizing" (treatment increases RMSF)

    Parameters
    ----------
    group1 : array_like
        First group (typically control)
    group2 : array_like
        Second group (typically treatment)
    rmsf_mode : bool, optional
        If True, interpret direction for RMSF (lower = better).
        Default is True.

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

    # Interpret direction for RMSF
    if rmsf_mode:
        # For RMSF: positive d means control > treatment, so treatment stabilizes
        direction = "stabilizing" if d > 0 else "destabilizing"
    else:
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
        Negative = reduction, Positive = increase
    """
    if control_mean == 0:
        return 0.0
    return (treatment_mean - control_mean) / control_mean * 100
