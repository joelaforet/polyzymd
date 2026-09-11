"""Statistical functions for replicate aggregation.

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
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import t as student_t

from polyzymd.analyses.exceptions import StatisticsError

UNCERTAINTY_KIND_REPLICATE_SEM = "sem_across_replicates"
"""Name of the only uncertainty this module estimates."""

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

    return float(student_t.ppf(0.5 + coverage / 2.0, n - 1))


def uncertainty_block(
    n: int | None,
    coverage: float = DEFAULT_COVERAGE,
    kind: str = UNCERTAINTY_KIND_REPLICATE_SEM,
) -> dict:
    """Return the serializable description of a reported uncertainty.

    Every aggregated artifact carries this block, with ``kind``, ``n``,
    ``coverage`` and ``method`` keys, so a reader never has to guess what an
    error bar or a ``sem`` field means. ``n`` is ``None`` when no replicate
    count is known.
    """
    return {
        "kind": kind,
        "n": None if n is None else int(n),
        "coverage": float(coverage),
        "method": CI_METHOD_STUDENT_T,
    }


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


def metric_summary_payload(
    name: str,
    values: Sequence[float],
    *,
    unit: str | None = None,
    coverage: float = DEFAULT_COVERAGE,
) -> dict:
    """Build the serialized metric summary for one condition-level metric.

    Every plugin writes its aggregated metrics in this shape, which matches
    ``polyzymd.analyses.mda.aggregation.AggregatedMetric``. Deriving the mean,
    the standard error, the spread and the interval from ``values`` in one place
    keeps the stored statistics consistent with the replicate values that the
    comparison layer recomputes them from. The spread, the standard error and
    both limits are ``None`` for a single replicate.
    """
    numeric = [float(value) for value in values]
    stats = mean_sem_ci(numeric, coverage=coverage)
    std = None if stats.sem is None else stats.sem * np.sqrt(float(stats.n))
    return {
        "name": name,
        "values": numeric,
        "mean": stats.mean,
        "sem": stats.sem,
        "std": None if std is None else float(std),
        "n": stats.n,
        "unit": unit,
        "ci95_low": stats.ci_low,
        "ci95_high": stats.ci_high,
        "ci_method": stats.ci_method,
    }


@dataclass
class StatResult:
    """Container for a mean with its standard uncertainty and interval.

    The standard error, both confidence limits and the method are ``None`` when
    a single replicate makes them inestimable.
    """

    mean: float
    sem: float | None
    n_samples: int
    ci95_low: float | None = None
    ci95_high: float | None = None
    ci_method: str | None = None

    def __repr__(self) -> str:
        if self.sem is None:
            return f"{self.mean:.4f} (n={self.n_samples}, SEM not estimable)"
        return f"{self.mean:.4f} ± {self.sem:.4f} (n={self.n_samples})"

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "mean": float(self.mean),
            "sem": None if self.sem is None else float(self.sem),
            "n_samples": self.n_samples,
            "ci95_low": None if self.ci95_low is None else float(self.ci95_low),
            "ci95_high": None if self.ci95_high is None else float(self.ci95_high),
            "ci_method": self.ci_method,
        }


@dataclass
class PerResidueStats:
    """Container for per-residue statistics across replicates.

    Attributes
    ----------
    residue_ids : NDArray[np.int64]
        Residue identifiers (1-indexed, following PyMOL convention)
    means : NDArray[np.float64]
        Mean value for each residue across replicates
    sems : NDArray[np.float64]
        SEM for each residue across replicates
    n_replicates : int
        Number of replicates aggregated
    """

    residue_ids: NDArray[np.int64]
    means: NDArray[np.float64]
    sems: NDArray[np.float64]
    n_replicates: int

    def __len__(self) -> int:
        return len(self.residue_ids)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "residue_ids": self.residue_ids.tolist(),
            "means": self.means.tolist(),
            "sems": self.sems.tolist(),
            "n_replicates": self.n_replicates,
        }


def compute_sem(values: ArrayLike, ddof: int = 1) -> StatResult:
    """Compute the mean, SEM and 95 percent confidence interval.

    A thin wrapper over :func:`mean_sem_ci`. ``ddof`` other than 1 is rejected
    because a confidence interval needs the sample standard deviation.

    Examples
    --------
    >>> values = [2.5, 2.7, 2.6, 2.4, 2.8]  # RMSF from 5 replicates
    >>> result = compute_sem(values)
    >>> print(f"RMSF = {result.mean:.2f} +/- {result.sem:.2f} A")
    RMSF = 2.60 +/- 0.07 A

    """
    if ddof != 1:
        raise StatisticsError(
            f"compute_sem requires the sample standard deviation (ddof=1); got ddof={ddof!r}"
        )

    result = mean_sem_ci(values)
    return StatResult(
        mean=result.mean,
        sem=result.sem,
        n_samples=result.n,
        ci95_low=result.ci_low,
        ci95_high=result.ci_high,
        ci_method=result.ci_method,
    )


def aggregate_per_residue_stats(
    per_replicate_values: Sequence[NDArray[np.float64]],
    residue_ids: NDArray[np.int64] | None = None,
) -> PerResidueStats:
    """Aggregate per-residue values across replicates.

    For each residue, computes mean +/- SEM across all replicates.
    This is the correct way to aggregate per-residue RMSF values.

    Parameters
    ----------
    per_replicate_values : sequence of arrays
        List/tuple of 1D arrays, each containing per-residue values
        from one replicate. All arrays must have the same length.
    residue_ids : array, optional
        1-indexed residue identifiers. If None, uses 1, 2, 3, ...
        Following PyMOL convention (1-indexed).

    Returns
    -------
    PerResidueStats
        Container with residue_ids, means, sems, n_replicates

    Raises
    ------
    ValueError
        If arrays have inconsistent lengths or no replicates provided
    """
    if len(per_replicate_values) == 0:
        raise ValueError("No replicate values provided")

    # Stack into 2D array: (n_replicates, n_residues)
    stacked = np.stack([np.asarray(v, dtype=np.float64) for v in per_replicate_values])
    n_replicates, n_residues = stacked.shape

    # Compute mean and SEM along replicate axis
    means = np.mean(stacked, axis=0)

    if n_replicates == 1:
        sems = np.zeros(n_residues, dtype=np.float64)
    else:
        stds = np.std(stacked, axis=0, ddof=1)
        sems = stds / np.sqrt(n_replicates)

    # Generate residue IDs if not provided (1-indexed!)
    if residue_ids is None:
        residue_ids = np.arange(1, n_residues + 1, dtype=np.int64)
    else:
        residue_ids = np.asarray(residue_ids, dtype=np.int64)
        if len(residue_ids) != n_residues:
            raise ValueError(
                f"residue_ids length ({len(residue_ids)}) doesn't match data length ({n_residues})"
            )

    return PerResidueStats(
        residue_ids=residue_ids,
        means=means,
        sems=sems,
        n_replicates=n_replicates,
    )


def aggregate_region_stats(
    per_replicate_values: Sequence[NDArray[np.float64]],
    residue_mask: NDArray[np.bool_] | None = None,
) -> StatResult:
    """Aggregate region-averaged values across replicates.

    For whole-protein or region-specific metrics, this computes the mean
    of per-replicate averages, with SEM across replicates.

    This implements the correct hierarchical aggregation:
    1. First average within each replicate (over selected residues)
    2. Then compute mean +/- SEM across replicate means

    Parameters
    ----------
    per_replicate_values : sequence of arrays
        List/tuple of 1D arrays, each containing per-residue values
        from one replicate.
    residue_mask : bool array, optional
        Boolean mask for residue selection. If None, uses all residues.

    Returns
    -------
    StatResult
        Mean +/- SEM of region-averaged values across replicates
    """
    if len(per_replicate_values) == 0:
        raise ValueError("No replicate values provided")

    # Compute per-replicate region averages
    replicate_means = []
    for values in per_replicate_values:
        arr = np.asarray(values, dtype=np.float64)
        if residue_mask is not None:
            arr = arr[residue_mask]
        replicate_means.append(float(np.mean(arr)))

    # Compute mean +/- SEM across replicates
    return compute_sem(replicate_means)


def weighted_mean_with_sem(
    means: ArrayLike,
    sems: ArrayLike,
    weights: ArrayLike | None = None,
) -> StatResult:
    """Compute weighted mean with proper error propagation.

    Useful for combining results from different conditions or
    analyses with different uncertainties.

    Parameters
    ----------
    means : array_like
        Mean values from each source
    sems : array_like
        SEM values from each source
    weights : array_like, optional
        Weights for each source. If None, uses inverse-variance weighting
        (1/sem^2), which is optimal for independent measurements.

    Returns
    -------
    StatResult
        Weighted mean with propagated uncertainty

    Notes
    -----
    For inverse-variance weighting, the combined SEM is:
        SEM_combined = 1 / sqrt(sum(1/SEM_i^2))

    For arbitrary weights, uses standard error propagation:
        SEM_combined = sqrt(sum((w_i * SEM_i)^2)) / sum(w_i)
    """
    means_arr = np.asarray(means, dtype=np.float64)
    sems_arr = np.asarray(sems, dtype=np.float64)
    n = len(means_arr)

    if n == 0:
        raise ValueError("Cannot compute weighted mean of empty arrays")

    if len(sems_arr) != n:
        raise ValueError("means and sems must have same length")

    if weights is None:
        if np.any(~np.isfinite(sems_arr)):
            raise ValueError(
                f"Inverse-variance weighting requires finite SEM values; got {sems_arr.tolist()}"
            )
        if np.any(sems_arr <= 0.0):
            raise ValueError(
                "Inverse-variance weighting requires SEM values greater than zero; "
                f"got {sems_arr.tolist()}"
            )
        # Inverse-variance weighting
        inv_var = 1.0 / (sems_arr**2)
        weights_arr = inv_var
    else:
        weights_arr = np.asarray(weights, dtype=np.float64)
        if len(weights_arr) != n:
            raise ValueError("weights must have same length as means")

    # Normalize weights
    weight_sum = np.sum(weights_arr)
    norm_weights = weights_arr / weight_sum

    # Weighted mean
    weighted_mean = float(np.sum(norm_weights * means_arr))

    # Error propagation: sqrt(sum((w_i * SEM_i)^2))
    # But for inverse-variance weighting, use the optimal formula
    if weights is None:
        combined_sem = float(1.0 / np.sqrt(np.sum(inv_var)))
    else:
        combined_sem = float(np.sqrt(np.sum((norm_weights * sems_arr) ** 2)))

    return StatResult(
        mean=weighted_mean,
        sem=combined_sem,
        n_samples=n,
    )
