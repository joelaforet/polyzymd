"""Runner-backed distance execution helpers.

This module keeps the trajectory-native distance work behind the runner seam so
plugins can preserve their result schemas while PolyzyMD owns replicate,
aggregation, and comparison orchestration.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Callable, Sequence

import numpy as np
from numpy.typing import NDArray

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class DistancePairPayload:
    """Trajectory-native payload for one distance pair.

    Parameters
    ----------
    pair_label : str
        Human-readable pair label derived from the selections.
    selection1 : str
        First selection string.
    selection2 : str
        Second selection string.
    distances : NDArray[np.float64]
        Per-frame distance values in Å.
    mean_distance : float
        Mean distance over analyzed frames.
    std_distance : float
        Standard deviation of the distance series.
    median_distance : float
        Median distance.
    min_distance : float
        Minimum observed distance.
    max_distance : float
        Maximum observed distance.
    sem_distance : float | None
        Autocorrelation-corrected SEM when available.
    correlation_time : float | None
        Estimated correlation time.
    correlation_time_unit : str | None
        Unit for the correlation time.
    n_independent_frames : int | None
        Effective independent frame count.
    statistical_inefficiency : float | None
        Statistical inefficiency estimate.
    autocorrelation_warning : str | None
        Reliability warning from the autocorrelation analysis.
    threshold : float | None
        Distance threshold used for contact analysis.
    fraction_below_threshold : float | None
        Fraction of frames below ``threshold``.
    histogram_edges : NDArray[np.float64]
        Histogram bin edges.
    histogram_counts : NDArray[np.int64]
        Histogram bin counts.
    kde_x : NDArray[np.float64] | None
        KDE evaluation points.
    kde_y : NDArray[np.float64] | None
        KDE density values.
    kde_peak : float | None
        KDE-derived mode.
    kde_bandwidth : float | None
        KDE bandwidth.
    n_frames_total : int
        Total number of frames in the trajectory.
    n_frames_used : int
        Number of frames analyzed for this pair.
    """

    pair_label: str
    selection1: str
    selection2: str
    distances: NDArray[np.float64]
    mean_distance: float
    std_distance: float
    median_distance: float
    min_distance: float
    max_distance: float
    sem_distance: float | None
    correlation_time: float | None
    correlation_time_unit: str | None
    n_independent_frames: int | None
    statistical_inefficiency: float | None
    autocorrelation_warning: str | None
    threshold: float | None
    fraction_below_threshold: float | None
    histogram_edges: NDArray[np.float64]
    histogram_counts: NDArray[np.int64]
    kde_x: NDArray[np.float64] | None
    kde_y: NDArray[np.float64] | None
    kde_peak: float | None
    kde_bandwidth: float | None
    n_frames_total: int
    n_frames_used: int


@dataclass(frozen=True)
class DistancesRunnerPayload:
    """Collection of distance-pair payloads for one replicate."""

    n_frames_total: int
    n_frames_used: int
    pair_payloads: list[DistancePairPayload]


@dataclass(frozen=True)
class _ResolvedDistancePair:
    """Validated atom-group inputs for one distance pair."""

    pair_label: str
    selection1: str
    selection2: str
    threshold: float | None
    atoms1: Any
    atoms2: Any
    mode1: Any
    mode2: Any


class DistancesReplicateRunner:
    """Execute distances analysis for one replicate through the runner seam."""

    def __init__(
        self,
        *,
        universe: Any,
        pairs: Sequence[tuple[str, str]],
        thresholds: Sequence[float | None],
        use_pbc: bool,
        alignment: Any,
        timestep_ps: float,
        pair_label_func: Callable[[str, str], str],
    ) -> None:
        """Store replicate-specific distance runner state.

        Parameters
        ----------
        universe : Any
            MDAnalysis universe for the replicate.
        pairs : Sequence[tuple[str, str]]
            Selection pairs to analyze.
        thresholds : Sequence[float | None]
            Per-pair thresholds for contact statistics.
        use_pbc : bool
            Whether to use the minimum-image convention.
        alignment : Any
            Alignment configuration object.
        timestep_ps : float
            Trajectory timestep in picoseconds.
        pair_label_func : Callable[[str, str], str]
            Helper used to derive legacy pair labels from selections.
        """

        self._universe = universe
        self._pairs = list(pairs)
        self._thresholds = list(thresholds)
        self._use_pbc = use_pbc
        self._alignment = alignment
        self._timestep_ps = float(timestep_ps)
        self._pair_label_func = pair_label_func
        self.results = DistancesRunnerPayload(
            n_frames_total=len(universe.trajectory),
            n_frames_used=0,
            pair_payloads=[],
        )

    def run(self, start: int, stop: int, step: int = 1) -> DistancesReplicateRunner:
        """Execute all configured distance pairs.

        Parameters
        ----------
        start : int
            Inclusive start frame.
        stop : int
            Exclusive stop frame.
        step : int, optional
            Frame stride, by default 1.

        Returns
        -------
        DistancesReplicateRunner
            Runner instance with populated ``results``.
        """

        self.results = compute_distance_payloads(
            universe=self._universe,
            pairs=self._pairs,
            thresholds=self._thresholds,
            start=start,
            stop=stop,
            step=step,
            timestep_ps=self._timestep_ps,
            use_pbc=self._use_pbc,
            alignment=self._alignment,
            pair_label_func=self._pair_label_func,
        )
        return self


def compute_distance_payloads(
    *,
    universe: Any,
    pairs: Sequence[tuple[str, str]],
    thresholds: Sequence[float | None],
    start: int,
    stop: int,
    step: int,
    timestep_ps: float,
    use_pbc: bool,
    alignment: Any,
    pair_label_func: Callable[[str, str], str],
) -> DistancesRunnerPayload:
    """Compute distance payloads for a batch of pairs.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe for the replicate.
    pairs : Sequence[tuple[str, str]]
        Selection pairs to analyze.
    thresholds : Sequence[float | None]
        Thresholds used for contact statistics.
    start : int
        Inclusive start frame.
    stop : int
        Exclusive stop frame.
    step : int
        Frame stride.
    timestep_ps : float
        Timestep between saved frames in ps.
    use_pbc : bool
        Whether to use the minimum-image convention.
    alignment : Any
        Alignment configuration object.
    pair_label_func : Callable[[str, str], str]
        Helper used to derive legacy pair labels from selections.

    Returns
    -------
    DistancesRunnerPayload
        Summarized distance payloads for the replicate.
    """

    from polyzymd.analyses.shared.alignment import align_trajectory

    if getattr(alignment, "enabled", False):
        align_trajectory(universe, alignment, start_frame=start, stop_frame=stop)

    if use_pbc:
        LOGGER.info("Using PBC-aware distance calculation (minimum image convention)")
    else:
        LOGGER.debug("PBC disabled; using simple Euclidean distances")

    resolved_pairs = _resolve_distance_pairs(
        universe=universe,
        pairs=pairs,
        thresholds=thresholds,
        pair_label_func=pair_label_func,
    )
    analysis_cls = _get_distance_analysis_base_cls()
    analysis = analysis_cls(
        universe=universe,
        resolved_pairs=resolved_pairs,
        use_pbc=use_pbc,
    ).run(start=start, stop=stop, step=step)

    pair_payloads = [
        _summarize_distance_series(
            distances=np.asarray(distances, dtype=np.float64),
            pair_label=resolved_pair.pair_label,
            selection1=resolved_pair.selection1,
            selection2=resolved_pair.selection2,
            threshold=resolved_pair.threshold,
            timestep_ps=timestep_ps,
            n_frames_total=len(universe.trajectory),
        )
        for resolved_pair, distances in zip(
            resolved_pairs,
            analysis.results.distance_arrays,
            strict=True,
        )
    ]
    n_frames_used = len(range(start, stop, step))
    return DistancesRunnerPayload(
        n_frames_total=len(universe.trajectory),
        n_frames_used=n_frames_used,
        pair_payloads=pair_payloads,
    )


def _resolve_distance_pairs(
    *,
    universe: Any,
    pairs: Sequence[tuple[str, str]],
    thresholds: Sequence[float | None],
    pair_label_func: Callable[[str, str], str],
) -> list[_ResolvedDistancePair]:
    """Resolve selections and validate atom groups before trajectory iteration."""

    from polyzymd.analyses.shared.diagnostics import (
        get_selection_diagnostics,
        warn_if_multi_chain_selection,
    )
    from polyzymd.analyses.shared.selections import parse_selection_string

    resolved_pairs: list[_ResolvedDistancePair] = []
    for (selection1, selection2), threshold in zip(pairs, thresholds, strict=True):
        parsed1 = parse_selection_string(selection1)
        parsed2 = parse_selection_string(selection2)
        atoms1 = universe.select_atoms(parsed1.selection)
        atoms2 = universe.select_atoms(parsed2.selection)

        if len(atoms1) == 0:
            diagnostics = get_selection_diagnostics(universe, selection1)
            raise ValueError(f"Selection '{selection1}' matched no atoms.\n\n{diagnostics}")
        if len(atoms2) == 0:
            diagnostics = get_selection_diagnostics(universe, selection2)
            raise ValueError(f"Selection '{selection2}' matched no atoms.\n\n{diagnostics}")

        pair_label = pair_label_func(selection1, selection2)
        warn_if_multi_chain_selection(atoms1, selection1, f"for distance pair '{pair_label}'")
        warn_if_multi_chain_selection(atoms2, selection2, f"for distance pair '{pair_label}'")

        resolved_pairs.append(
            _ResolvedDistancePair(
                pair_label=pair_label,
                selection1=selection1,
                selection2=selection2,
                threshold=threshold,
                atoms1=atoms1,
                atoms2=atoms2,
                mode1=parsed1.mode,
                mode2=parsed2.mode,
            )
        )
    return resolved_pairs


def _get_distance_analysis_base_cls() -> type:
    """Return a cached ``AnalysisBase`` subclass for pair-distance collection."""

    cached = getattr(_get_distance_analysis_base_cls, "_cls", None)
    if cached is not None:
        return cached

    from MDAnalysis.analysis.base import AnalysisBase

    class _DistanceAnalysisBase(AnalysisBase):  # type: ignore[misc]
        """Collect per-pair distances while MDAnalysis owns trajectory iteration."""

        def __init__(
            self,
            *,
            universe: Any,
            resolved_pairs: Sequence[_ResolvedDistancePair],
            use_pbc: bool,
        ) -> None:
            super().__init__(universe.trajectory)
            self._resolved_pairs = list(resolved_pairs)
            self._use_pbc = use_pbc

        def _prepare(self) -> None:
            self.results.distance_arrays = [[] for _ in self._resolved_pairs]

        def _single_frame(self) -> None:
            from polyzymd.analyses.shared.pbc import minimum_image_distance
            from polyzymd.analyses.shared.selections import get_position

            box = self._ts.dimensions
            for pair_index, pair in enumerate(self._resolved_pairs):
                position1 = get_position(pair.atoms1, pair.mode1)
                position2 = get_position(pair.atoms2, pair.mode2)
                if self._use_pbc:
                    distance = minimum_image_distance(position1, position2, box)
                else:
                    distance = float(np.linalg.norm(position2 - position1))
                self.results.distance_arrays[pair_index].append(distance)

        def _conclude(self) -> None:
            self.results.distance_arrays = [
                np.asarray(values, dtype=np.float64) for values in self.results.distance_arrays
            ]

    _get_distance_analysis_base_cls._cls = _DistanceAnalysisBase  # type: ignore[attr-defined]
    return _DistanceAnalysisBase


def _summarize_distance_series(
    *,
    distances: NDArray[np.float64],
    pair_label: str,
    selection1: str,
    selection2: str,
    threshold: float | None,
    timestep_ps: float,
    n_frames_total: int,
) -> DistancePairPayload:
    """Summarize one distance series into the legacy payload schema."""

    from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

    n_frames_used = len(distances)
    mean_distance = float(np.mean(distances))
    std_distance = float(np.std(distances))
    median_distance = float(np.median(distances))
    min_distance = float(np.min(distances))
    max_distance = float(np.max(distances))

    fraction_below_threshold: float | None = None
    if threshold is not None:
        fraction_below_threshold = float(np.mean(distances < threshold))

    histogram_counts, histogram_edges = np.histogram(distances, bins=50)

    kde_x: NDArray[np.float64] | None = None
    kde_y: NDArray[np.float64] | None = None
    kde_peak: float | None = None
    kde_bandwidth: float | None = None
    try:
        from scipy.stats import gaussian_kde

        if len(distances) > 10:
            kde = gaussian_kde(distances)
            std_value = float(np.std(distances))
            kde_bandwidth = float(kde.factor) * std_value
            x_min = max(0.0, min_distance - 0.5)
            x_max = max_distance + 0.5
            kde_x = np.linspace(x_min, x_max, 200, dtype=np.float64)
            kde_y = np.asarray(kde(kde_x), dtype=np.float64)
            kde_peak = float(kde_x[int(np.argmax(kde_y))])
            LOGGER.debug("KDE peak (mode) for %s: %.2f Å", pair_label, kde_peak)
    except ImportError:
        pass
    except (ValueError, np.linalg.LinAlgError, RuntimeError) as exc:
        LOGGER.warning("KDE computation failed for %s: %s", pair_label, exc)

    sem_distance: float | None = None
    correlation_time: float | None = None
    correlation_time_unit: str | None = None
    n_independent_frames: int | None = None
    statistical_inefficiency: float | None = None
    autocorrelation_warning: str | None = None
    if len(distances) >= 20:
        try:
            tau_result = estimate_correlation_time(
                distances,
                timestep=timestep_ps,
                timestep_unit="ps",
                method="integration",
                n_frames=n_frames_used,
            )
            correlation_time = tau_result.tau
            correlation_time_unit = tau_result.tau_unit
            n_independent_frames = tau_result.n_independent
            statistical_inefficiency = tau_result.statistical_inefficiency
            autocorrelation_warning = tau_result.warning
            if n_independent_frames > 0:
                sem_distance = float(std_distance / np.sqrt(float(n_independent_frames)))
        except (ValueError, np.linalg.LinAlgError, RuntimeError) as exc:
            LOGGER.warning("Autocorrelation analysis failed for %s: %s", pair_label, exc)
            sem_distance = float(std_distance / np.sqrt(float(n_frames_used)))

    return DistancePairPayload(
        pair_label=pair_label,
        selection1=selection1,
        selection2=selection2,
        distances=np.asarray(distances, dtype=np.float64),
        mean_distance=mean_distance,
        std_distance=std_distance,
        median_distance=median_distance,
        min_distance=min_distance,
        max_distance=max_distance,
        sem_distance=sem_distance,
        correlation_time=correlation_time,
        correlation_time_unit=correlation_time_unit,
        n_independent_frames=n_independent_frames,
        statistical_inefficiency=statistical_inefficiency,
        autocorrelation_warning=autocorrelation_warning,
        threshold=threshold,
        fraction_below_threshold=fraction_below_threshold,
        histogram_edges=np.asarray(histogram_edges, dtype=np.float64),
        histogram_counts=np.asarray(histogram_counts, dtype=np.int64),
        kde_x=kde_x,
        kde_y=kde_y,
        kde_peak=kde_peak,
        kde_bandwidth=kde_bandwidth,
        n_frames_total=n_frames_total,
        n_frames_used=n_frames_used,
    )
