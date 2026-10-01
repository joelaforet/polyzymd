"""The mean, standard error and Student t interval of replicate values.

The replicate is the sampling unit. Every uncertainty here is computed across
replicates, never across frames. The standard uncertainty is the standard error
of the mean, ``s / sqrt(n)`` with the sample standard deviation, and the
reported interval is the two-sided Student t interval
``mean +/- t(1 - alpha/2, n - 1) * SEM``. At n = 3 that coverage factor is
4.303, so a plus-or-minus-one-SEM band covers far less than 95 percent and must
not be read as one. A single replicate supports neither, so those fields are
``None`` rather than ``0.0``.

References
----------
.. [1] Grossfield, A.; Patrone, P. N.; Roe, D. R.; Schultz, A. J.;
       Siderius, D. W.; Zuckerman, D. M. Best Practices for Quantification of
       Uncertainty and Sampling Quality in Molecular Simulations.
       Living J. Comput. Mol. Sci. 2018, 1 (1), 5067.
       https://doi.org/10.33011/livecoms.1.1.5067
.. [2] Joint Committee for Guides in Metrology. Evaluation of Measurement
       Data: Guide to the Expression of Uncertainty in Measurement,
       JCGM 100:2008; BIPM: Sevres, 2008.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from polyzymd.analyses.exceptions import StatisticsError

CI_METHOD_STUDENT_T = "student_t"
"""Name of the interval method reported alongside every confidence interval."""

DEFAULT_COVERAGE = 0.95
"""Coverage probability of the reported confidence interval."""


def student_t_coverage_factor(n: int, coverage: float = DEFAULT_COVERAGE) -> float | None:
    """Return the Student t coverage factor for *n* replicates.

    The factor is ``k = t(1 - alpha / 2, n - 1)``, which is 4.303 at n = 3 and
    2.776 at n = 5 for 95 percent coverage. Returns ``None`` when ``n < 2``, and
    raises ``StatisticsError`` when ``coverage`` is outside ``(0, 1)``.
    """
    if not 0.0 < coverage < 1.0:
        raise StatisticsError(f"coverage must be in (0, 1), got {coverage!r}")
    if n < 2:
        return None

    # Imported here, not at module level: scipy.stats costs half a second and
    # most imports of this module never take a quantile.
    from scipy.stats import t as student_t

    return float(student_t.ppf(0.5 + coverage / 2.0, n - 1))


@dataclass(frozen=True)
class MeanSemCI:
    """Mean, standard uncertainty and Student t confidence interval.

    Every field except ``mean`` and ``n`` is ``None`` for a single replicate,
    where no spread exists.
    """

    mean: float
    sem: float | None
    n: int
    ci_low: float | None
    ci_high: float | None
    ci_method: str | None
    coverage: float | None

    def to_dict(self) -> dict:
        """Convert to a dictionary for JSON serialization."""
        return {
            "mean": float(self.mean),
            "sem": None if self.sem is None else float(self.sem),
            "n": int(self.n),
            "ci95_low": None if self.ci_low is None else float(self.ci_low),
            "ci95_high": None if self.ci_high is None else float(self.ci_high),
            "ci_method": self.ci_method,
            "coverage": None if self.coverage is None else float(self.coverage),
        }


def mean_sem_ci(values: ArrayLike, coverage: float = DEFAULT_COVERAGE) -> MeanSemCI:
    """Compute the mean, SEM and Student t confidence interval of *values*.

    This is the one interval estimator in the analyses package. Everything that
    reports a condition-level uncertainty goes through it so that the coverage
    factor, the degrees of freedom and the single-replicate rule are the same
    everywhere. Raises ``StatisticsError`` on an empty sample.

    Examples
    --------
    >>> result = mean_sem_ci([2.0, 2.2, 2.4])
    >>> round(result.ci_high - result.mean, 4)
    0.4968
    """
    array = np.asarray(values, dtype=np.float64).ravel()
    n = int(array.size)
    if n == 0:
        raise StatisticsError("Cannot compute statistics on empty array")

    mean = float(np.mean(array))
    factor = student_t_coverage_factor(n, coverage)
    if factor is None:
        return MeanSemCI(
            mean=mean,
            sem=None,
            n=n,
            ci_low=None,
            ci_high=None,
            ci_method=None,
            coverage=None,
        )

    sem = float(np.std(array, ddof=1) / np.sqrt(float(n)))
    half_width = factor * sem
    return MeanSemCI(
        mean=mean,
        sem=sem,
        n=n,
        ci_low=mean - half_width,
        ci_high=mean + half_width,
        ci_method=CI_METHOD_STUDENT_T,
        coverage=float(coverage),
    )
