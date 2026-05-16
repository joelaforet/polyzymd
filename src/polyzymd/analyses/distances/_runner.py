"""Compatibility wrappers for distance runner consumers.

The distances plugin itself uses MDAnalysis-native jobs. This module remains as
a thin adapter for the catalytic-triad plugin until its dedicated migration
replaces that runner path.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from polyzymd.analyses.distances._mda import (
    DistancePairPayload,
    DistancesRunnerPayload,
    compute_distance_payloads,
)


class DistancesReplicateRunner:
    """Execute pair distances through the shared MDAnalysis primitive.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe for the replicate.
    pairs : sequence of tuple[str, str]
        Selection pairs to analyze.
    thresholds : sequence of float or None
        Per-pair thresholds for contact statistics.
    use_pbc : bool
        Whether to use minimum-image distances.
    alignment : Any
        Alignment configuration object.
    timestep_ps : float
        Trajectory timestep in picoseconds.
    pair_label_func : Callable[[str, str], str]
        Helper used to derive pair labels from selections.
    """

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
        self._universe = universe
        self._pairs = list(pairs)
        self._thresholds = list(thresholds)
        self._use_pbc = bool(use_pbc)
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


__all__ = [
    "DistancePairPayload",
    "DistancesReplicateRunner",
    "DistancesRunnerPayload",
    "compute_distance_payloads",
]
