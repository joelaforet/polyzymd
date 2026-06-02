"""Pair-distance ``AnalysisBase`` primitives for MDAnalysis integrations."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses.shared.statistics import StatResult, compute_sem

LOGGER = logging.getLogger(__name__)

__all__ = [
    "PairDistanceSpec",
    "build_pair_distance_analysis",
    "pair_distance_version",
]


@dataclass(frozen=True)
class PairDistanceSpec:
    """Resolved atom-group inputs for one pair-distance measurement.

    Parameters
    ----------
    label : str
        Human-readable pair label.
    selection_a : str
        Original selection string for the first atom group or point.
    selection_b : str
        Original selection string for the second atom group or point.
    atoms_a : Any
        First MDAnalysis atom group.
    atoms_b : Any
        Second MDAnalysis atom group.
    mode_a : Any
        Position-reduction mode understood by shared selection helpers.
    mode_b : Any
        Position-reduction mode understood by shared selection helpers.
    threshold : float or None, optional
        Optional distance threshold in Å for downstream state summaries.
    """

    label: str
    selection_a: str
    selection_b: str
    atoms_a: Any
    atoms_b: Any
    mode_a: Any
    mode_b: Any
    threshold: float | None = None


@dataclass
class PairAggregatedStats:
    """Aggregated statistics for a single distance pair.

    Parameters
    ----------
    mean_stats : StatResult
        Mean distance across replicates.
    median_stats : StatResult
        Median distance across replicates.
    fraction_stats : StatResult or None
        Fraction below threshold, if available.
    kde_peak_stats : StatResult or None
        KDE peak distance, if available.
    per_rep_means : list[float]
        Per-replicate mean distances.
    per_rep_stds : list[float]
        Per-replicate standard deviations.
    per_rep_medians : list[float]
        Per-replicate median distances.
    per_rep_fractions : list[float]
        Per-replicate fractions below threshold.
    per_rep_kde_peaks : list[float]
        Per-replicate KDE peak distances.
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


def aggregate_distance_pair_stats(
    individual_results: Sequence[Any],
    pair_idx: int,
) -> PairAggregatedStats:
    """Aggregate per-pair distance statistics across replicate results.

    Parameters
    ----------
    individual_results : sequence
        Per-replicate result objects with indexable ``pair_results`` entries.
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
        pair_result = result.pair_results[pair_idx]
        per_rep_means.append(pair_result.mean_distance)
        per_rep_stds.append(pair_result.std_distance)
        per_rep_medians.append(pair_result.median_distance)
        if pair_result.fraction_below_threshold is not None:
            per_rep_fractions.append(pair_result.fraction_below_threshold)
        if pair_result.kde_peak is not None:
            per_rep_kde_peaks.append(pair_result.kde_peak)

    return PairAggregatedStats(
        mean_stats=compute_sem(per_rep_means),
        median_stats=compute_sem(per_rep_medians),
        fraction_stats=compute_sem(per_rep_fractions) if per_rep_fractions else None,
        kde_peak_stats=compute_sem(per_rep_kde_peaks) if per_rep_kde_peaks else None,
        per_rep_means=per_rep_means,
        per_rep_stds=per_rep_stds,
        per_rep_medians=per_rep_medians,
        per_rep_fractions=per_rep_fractions,
        per_rep_kde_peaks=per_rep_kde_peaks,
    )


def build_pair_distance_analysis(
    *,
    universe: Any,
    pairs: Sequence[PairDistanceSpec],
    use_pbc: bool,
) -> Any:
    """Build a lazy custom ``AnalysisBase`` for pair-distance matrices.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe for one trajectory.
    pairs : sequence of PairDistanceSpec
        Resolved pair specifications.
    use_pbc : bool
        Whether to request minimum-image distances from MDAnalysis.

    Returns
    -------
    Any
        ``AnalysisBase`` instance whose ``results.distance_matrix`` has shape
        ``(n_pairs, n_frames)`` and whose results also include frame, time, and
        warning metadata.
    """

    from MDAnalysis.analysis.base import AnalysisBase
    from MDAnalysis.lib.distances import calc_bonds

    from polyzymd.analyses.shared.selections import get_position

    class PairDistanceAnalysis(AnalysisBase):  # type: ignore[misc]
        """Collect pair distances while MDAnalysis owns frame iteration."""

        def __init__(self) -> None:
            self._pairs = list(pairs)
            self._use_pbc = bool(use_pbc)
            self._warnings: list[str] = []
            super().__init__(universe.trajectory)

        def _prepare(self) -> None:
            """Initialize matrix rows before trajectory iteration."""

            self.results.distance_matrix = [[] for _ in self._pairs]
            self.results.warnings = []

        def _single_frame(self) -> None:
            """Measure all pair distances for the current frame."""

            if not self._pairs:
                return
            positions_a = np.asarray(
                [get_position(pair.atoms_a, pair.mode_a) for pair in self._pairs],
                dtype=np.float64,
            )
            positions_b = np.asarray(
                [get_position(pair.atoms_b, pair.mode_b) for pair in self._pairs],
                dtype=np.float64,
            )
            box = self._pbc_box()
            distances = calc_bonds(positions_a, positions_b, box=box).astype(np.float64)
            for pair_index, distance in enumerate(distances):
                self.results.distance_matrix[pair_index].append(float(distance))

        def _conclude(self) -> None:
            """Store arrays and metadata after frame iteration."""

            self.results.distance_matrix = np.asarray(
                self.results.distance_matrix,
                dtype=np.float64,
            )
            self.results.frames = np.asarray(getattr(self, "frames", []), dtype=np.int64)
            self.results.times_ps = np.asarray(getattr(self, "times", []), dtype=np.float64)
            self.results.warnings = list(self._warnings)

        def _pbc_box(self) -> NDArray[np.float32] | None:
            """Return the current timestep box or disable PBC with one warning."""

            if not self._use_pbc:
                return None
            dimensions = getattr(self._ts, "dimensions", None)
            if dimensions is None:
                self._warn_once(
                    "PBC requested for pair distances, but the timestep has no box dimensions; "
                    "using non-PBC distances for affected frames."
                )
                return None
            box = np.asarray(dimensions, dtype=np.float32)
            if box.shape[0] < 6 or not np.all(np.isfinite(box[:6])) or np.any(box[:3] <= 0):
                self._warn_once(
                    "PBC requested for pair distances, but the timestep box is invalid; "
                    "using non-PBC distances for affected frames."
                )
                return None
            return box[:6]

        def _warn_once(self, message: str) -> None:
            """Record and log a warning message once per analysis run."""

            if message in self._warnings:
                return
            self._warnings.append(message)
            LOGGER.warning(message)

    return PairDistanceAnalysis()


def pair_distance_version() -> str:
    """Return the pair-distance primitive schema version.

    Returns
    -------
    str
        Version string for provenance records.
    """

    return "1"
