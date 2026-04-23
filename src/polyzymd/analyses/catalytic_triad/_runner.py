"""Runner-backed catalytic-triad execution helpers."""

from __future__ import annotations

from typing import Any, Sequence

from polyzymd.analyses.distances._runner import DistancesReplicateRunner, DistancesRunnerPayload


class CatalyticTriadReplicateRunner:
    """Delegate catalytic-triad pair distances to the distances runner."""

    def __init__(
        self,
        *,
        universe: Any,
        pair_selections: Sequence[tuple[str, str]],
        threshold: float,
        timestep_ps: float,
        pair_label_func: Any,
    ) -> None:
        """Store replicate-specific triad runner state.

        Parameters
        ----------
        universe : Any
            MDAnalysis universe for the replicate.
        pair_selections : Sequence[tuple[str, str]]
            Triad pair selections.
        threshold : float
            Shared simultaneous-contact threshold in Å.
        timestep_ps : float
            Trajectory timestep in picoseconds.
        pair_label_func : Any
            Helper used to derive legacy pair labels.
        """

        from polyzymd.analyses.shared.alignment import AlignmentConfig

        self._delegate = DistancesReplicateRunner(
            universe=universe,
            pairs=list(pair_selections),
            thresholds=[float(threshold)] * len(pair_selections),
            use_pbc=True,
            alignment=AlignmentConfig(),
            timestep_ps=timestep_ps,
            pair_label_func=pair_label_func,
        )
        self.results = DistancesRunnerPayload(
            n_frames_total=len(universe.trajectory),
            n_frames_used=0,
            pair_payloads=[],
        )

    def run(self, start: int, stop: int, step: int = 1) -> CatalyticTriadReplicateRunner:
        """Execute all configured triad pairs.

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
        CatalyticTriadReplicateRunner
            Runner instance with populated ``results``.
        """

        self._delegate.run(start=start, stop=stop, step=step)
        self.results = self._delegate.results
        return self
