"""Statistical functions for replicate aggregation.

This module provides statistical utilities for combining results across
multiple simulation replicates with proper error propagation.

Key design decisions:
- All uncertainties are reported as Standard Error of the Mean (SEM)
- SEM = std / sqrt(N) where N is the number of independent samples
- Hierarchical aggregation preserves proper statistics at each level

Functions
---------
compute_sem
    Standard error of the mean for a 1D array
aggregate_per_residue_stats
    Combine per-residue values across replicates
aggregate_region_stats
    Combine region-averaged values across replicates
weighted_mean_with_sem
    Weighted average with proper error propagation
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass
class StatResult:
    """Container for mean ± SEM results.

    Attributes
    ----------
    mean : float
        The mean value
    sem : float
        Standard error of the mean
    n_samples : int
        Number of samples used in computation
    """

    mean: float
    sem: float
    n_samples: int

    def __repr__(self) -> str:
        return f"{self.mean:.4f} ± {self.sem:.4f} (n={self.n_samples})"

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "mean": float(self.mean),
            "sem": float(self.sem),
            "n_samples": self.n_samples,
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
    """Compute mean and standard error of the mean.

    SEM = std / sqrt(N) where N is the number of samples.

    Parameters
    ----------
    values : array_like
        1D array of values (e.g., one value per replicate)
    ddof : int, optional
        Delta degrees of freedom for std calculation. Default is 1
        (Bessel's correction for sample std).

    Returns
    -------
    StatResult
        Container with mean, sem, and n_samples

    Examples
    --------
    >>> values = [2.5, 2.7, 2.6, 2.4, 2.8]  # RMSF from 5 replicates
    >>> result = compute_sem(values)
    >>> print(f"RMSF = {result.mean:.2f} ± {result.sem:.2f} Å")
    RMSF = 2.60 ± 0.07 Å

    Notes
    -----
    For a single value, SEM is undefined (returns 0.0).
    """
    arr = np.asarray(values, dtype=np.float64)
    n = len(arr)

    if n == 0:
        raise ValueError("Cannot compute statistics on empty array")

    mean = float(np.mean(arr))

    if n == 1:
        # Single sample: SEM is undefined, return 0
        return StatResult(mean=mean, sem=0.0, n_samples=1)

    std = float(np.std(arr, ddof=ddof))
    sem = std / np.sqrt(n)

    return StatResult(mean=mean, sem=sem, n_samples=n)


def aggregate_per_residue_stats(
    per_replicate_values: Sequence[NDArray[np.float64]],
    residue_ids: NDArray[np.int64] | None = None,
) -> PerResidueStats:
    """Aggregate per-residue values across replicates.

    For each residue, computes mean ± SEM across all replicates.
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

    Examples
    --------
    >>> # RMSF values for 3 residues from 4 replicates
    >>> rep1 = np.array([1.2, 2.5, 0.8])
    >>> rep2 = np.array([1.3, 2.4, 0.9])
    >>> rep3 = np.array([1.1, 2.6, 0.7])
    >>> rep4 = np.array([1.4, 2.3, 0.85])
    >>> stats = aggregate_per_residue_stats([rep1, rep2, rep3, rep4])
    >>> print(f"Residue 1 RMSF: {stats.means[0]:.2f} ± {stats.sems[0]:.2f}")
    Residue 1 RMSF: 1.25 ± 0.06

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
    2. Then compute mean ± SEM across replicate means

    Parameters
    ----------
    per_replicate_values : sequence of arrays
        List/tuple of 1D arrays, each containing per-residue values
        from one replicate.
    residue_mask : bool array, optional
        Boolean mask for residue selection. If None, uses all residues.
        E.g., to select active site residues only.

    Returns
    -------
    StatResult
        Mean ± SEM of region-averaged values across replicates

    Examples
    --------
    >>> # Whole-protein RMSF from 5 replicates
    >>> rmsf_rep1 = np.array([1.2, 2.5, 0.8, 1.5])  # 4 residues
    >>> rmsf_rep2 = np.array([1.3, 2.4, 0.9, 1.4])
    >>> # ... etc
    >>> result = aggregate_region_stats([rmsf_rep1, rmsf_rep2, ...])
    >>> print(f"Mean RMSF: {result.mean:.2f} ± {result.sem:.2f} Å")

    >>> # Active site only (residues 2 and 3)
    >>> mask = np.array([False, True, True, False])
    >>> result = aggregate_region_stats([rmsf_rep1, rmsf_rep2], mask)
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

    # Compute mean ± SEM across replicates
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

    Examples
    --------
    >>> # Combine RMSF from two different selection methods
    >>> means = [2.5, 2.7]
    >>> sems = [0.1, 0.15]
    >>> result = weighted_mean_with_sem(means, sems)
    """
    means_arr = np.asarray(means, dtype=np.float64)
    sems_arr = np.asarray(sems, dtype=np.float64)
    n = len(means_arr)

    if n == 0:
        raise ValueError("Cannot compute weighted mean of empty arrays")

    if len(sems_arr) != n:
        raise ValueError("means and sems must have same length")

    # Handle zero SEMs (replace with small value to avoid division by zero)
    sems_safe = np.where(sems_arr > 0, sems_arr, 1e-10)

    if weights is None:
        # Inverse-variance weighting
        inv_var = 1.0 / (sems_safe**2)
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
