"""Shared aggregation utilities for multi-replicate analysis.

This module extracts common patterns used by DistanceCalculator and
CatalyticTriadAnalyzer for resilient replicate collection and per-pair
statistical aggregation.

Functions
---------
collect_replicate_results
    Run analysis across replicates with error handling.
aggregate_distance_pair_stats
    Aggregate per-pair statistics across replicate results.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence, TypeVar

from polyzymd.analysis.core.statistics import StatResult, compute_sem

LOGGER = logging.getLogger(__name__)

T = TypeVar("T")


@dataclass
class ReplicateCollection:
    """Result of collecting replicate analysis results.

    Attributes
    ----------
    results : list
        Successfully computed results.
    successful_replicates : list[int]
        Replicate numbers that succeeded.
    failed_replicates : list[int]
        Replicate numbers that failed.
    """

    results: list[Any]
    successful_replicates: list[int]
    failed_replicates: list[int]


@dataclass
class PairAggregatedStats:
    """Aggregated statistics for a single distance pair.

    Attributes
    ----------
    mean_stats : StatResult
        Mean distance across replicates (mean ± SEM).
    median_stats : StatResult
        Median distance across replicates (mean of medians ± SEM).
    fraction_stats : StatResult or None
        Fraction below threshold (mean ± SEM), None if unavailable.
    kde_peak_stats : StatResult or None
        KDE peak distance (mean ± SEM), None if unavailable.
    per_rep_means : list[float]
        Per-replicate mean distances.
    per_rep_stds : list[float]
        Per-replicate standard deviations.
    per_rep_medians : list[float]
        Per-replicate median distances.
    per_rep_fractions : list[float]
        Per-replicate fractions below threshold (empty if unavailable).
    per_rep_kde_peaks : list[float]
        Per-replicate KDE peak distances (empty if unavailable).
    """

    mean_stats: StatResult
    median_stats: StatResult
    fraction_stats: StatResult | None
    kde_peak_stats: StatResult | None
    per_rep_means: list[float]
    per_rep_stds: list[float]
    per_rep_medians: list[float]
    per_rep_fractions: list[float]
    per_rep_kde_peaks: list[float]


def collect_replicate_results(
    compute_fn: Callable[..., T],
    replicates: Sequence[int],
    min_replicates: int = 2,
    output_dir_base: Path | None = None,
    **compute_kwargs: Any,
) -> ReplicateCollection:
    """Run analysis across replicates with resilient error handling.

    Iterates over replicate numbers, calling ``compute_fn(replicate=rep, **compute_kwargs)``
    for each. Failed replicates are skipped with a warning. Raises if fewer than
    ``min_replicates`` succeed.

    Parameters
    ----------
    compute_fn : callable
        Function to call for each replicate. Must accept ``replicate`` as
        a keyword argument.
    replicates : sequence of int
        Replicate numbers to process.
    min_replicates : int, optional
        Minimum number of successful replicates required. Default 2.
    output_dir_base : Path, optional
        Base directory for per-replicate output. When provided, each
        replicate call receives ``output_dir=output_dir_base / f"run_{rep}"``.
        When ``None`` (default), no ``output_dir`` kwarg is added and the
        compute function uses its own default.
    **compute_kwargs
        Additional keyword arguments passed to ``compute_fn``.

    Returns
    -------
    ReplicateCollection
        Container with results, successful replicate numbers, and failed
        replicate numbers.

    Raises
    ------
    ValueError
        If fewer than ``min_replicates`` succeed.
    """
    requested = list(replicates)
    results: list[Any] = []
    successful: list[int] = []
    failed: list[int] = []

    for rep in requested:
        try:
            call_kwargs = dict(compute_kwargs)
            if output_dir_base is not None:
                call_kwargs["output_dir"] = output_dir_base / f"run_{rep}"
            result = compute_fn(replicate=rep, **call_kwargs)
            results.append(result)
            successful.append(rep)
        except FileNotFoundError as e:
            LOGGER.warning(f"Skipping replicate {rep}: trajectory data not found. {e}")
            failed.append(rep)
        except Exception as e:
            LOGGER.warning(f"Skipping replicate {rep}: analysis failed with error: {e}")
            failed.append(rep)

    if len(results) < min_replicates:
        raise ValueError(
            f"Aggregation requires at least {min_replicates} successful replicates, "
            f"but only {len(results)} succeeded. Failed replicates: {failed}"
        )

    if failed:
        LOGGER.warning(
            f"Aggregating {len(successful)} of {len(requested)} "
            f"requested replicates. Skipped: {failed}"
        )

    return ReplicateCollection(
        results=results,
        successful_replicates=successful,
        failed_replicates=failed,
    )


def aggregate_distance_pair_stats(
    individual_results: Sequence[Any],
    pair_idx: int,
) -> PairAggregatedStats:
    """Aggregate per-pair distance statistics across replicate results.

    Collects mean, std, median, fraction_below_threshold, and kde_peak
    from each replicate's ``pair_results[pair_idx]`` and computes
    cross-replicate mean ± SEM for each.

    Parameters
    ----------
    individual_results : sequence
        List of per-replicate result objects. Each must have a
        ``pair_results`` attribute that is indexable by ``pair_idx``,
        and each pair result must have ``mean_distance``, ``std_distance``,
        ``median_distance``, and optionally ``fraction_below_threshold``
        and ``kde_peak`` attributes.
    pair_idx : int
        Index of the pair to aggregate.

    Returns
    -------
    PairAggregatedStats
        Aggregated statistics for this pair.
    """
    per_rep_means: list[float] = []
    per_rep_stds: list[float] = []
    per_rep_medians: list[float] = []
    per_rep_fractions: list[float] = []
    per_rep_kde_peaks: list[float] = []

    for result in individual_results:
        pr = result.pair_results[pair_idx]
        per_rep_means.append(pr.mean_distance)
        per_rep_stds.append(pr.std_distance)
        per_rep_medians.append(pr.median_distance)
        if pr.fraction_below_threshold is not None:
            per_rep_fractions.append(pr.fraction_below_threshold)
        if pr.kde_peak is not None:
            per_rep_kde_peaks.append(pr.kde_peak)

    mean_stats = compute_sem(per_rep_means)
    median_stats = compute_sem(per_rep_medians)

    fraction_stats = compute_sem(per_rep_fractions) if per_rep_fractions else None
    kde_peak_stats = compute_sem(per_rep_kde_peaks) if per_rep_kde_peaks else None

    return PairAggregatedStats(
        mean_stats=mean_stats,
        median_stats=median_stats,
        fraction_stats=fraction_stats,
        kde_peak_stats=kde_peak_stats,
        per_rep_means=per_rep_means,
        per_rep_stds=per_rep_stds,
        per_rep_medians=per_rep_medians,
        per_rep_fractions=per_rep_fractions,
        per_rep_kde_peaks=per_rep_kde_peaks,
    )
