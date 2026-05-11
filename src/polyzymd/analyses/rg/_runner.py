"""Runner-backed Rg execution helpers.

This module keeps trajectory-native radius-of-gyration work behind the runner
seam while preserving the plugin's legacy result schemas and NPZ sidecars.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Callable

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from polyzymd.analyses.rg import RgRunSettings

LOGGER = logging.getLogger(__name__)


@dataclass
class RgRunPayload:
    """Trajectory-native Rg data for one configured run."""

    run_label: str
    selection: str
    calculation_mode: str
    fragment_weighting: str | None
    rg_values: NDArray[np.float64]
    frames: NDArray[np.int64]
    time_ns: NDArray[np.float64]
    raw_timestep_ps: float
    frame_stride: int
    effective_timestep_ps: float
    mean_rg: float
    std_rg: float
    median_rg: float
    min_rg: float
    max_rg: float
    final_rg: float
    sem_rg: float | None
    correlation_time: float | None
    correlation_time_unit: str | None
    n_independent_frames: int | None
    statistical_inefficiency: float | None
    autocorrelation_warning: str | None
    fragment_rg_values: NDArray[np.float64] | None = None
    fragment_counts_per_frame: NDArray[np.int64] | None = None
    fragment_masses: NDArray[np.float64] | None = None
    frag_metadata: dict[str, float | int] = field(default_factory=dict)


@dataclass
class RgRunnerPayload:
    """Collection of Rg run payloads for one replicate."""

    n_frames_total: int
    run_payloads: list[RgRunPayload]
    skipped_run_payloads: list[RgSkippedRunPayload] = field(default_factory=list)


@dataclass
class RgSkippedRunPayload:
    """Provenance for one skipped Rg run in a replicate."""

    run_label: str
    selection: str
    replicate: int
    reason: str
    reason_code: str = "empty_selection"


class RgEmptySelectionError(ValueError):
    """Raised when an Rg run selection matches no atoms."""

    def __init__(self, *, run_label: str, selection: str, replicate: int) -> None:
        """Build a diagnostic empty-selection error.

        Parameters
        ----------
        run_label : str
            Configured run label.
        selection : str
            MDAnalysis selection string that matched no atoms.
        replicate : int
            Replicate number where the selection was empty.
        """

        self.run_label = run_label
        self.selection = selection
        self.replicate = replicate
        super().__init__(
            f"Run '{run_label}' selection matched no atoms for replicate "
            f"{replicate}: {selection!r}"
        )


class RgReplicateRunner:
    """Execute configured Rg runs for one replicate."""

    def __init__(
        self,
        *,
        sim_config: Any,
        replicate: int,
        runs: list["RgRunSettings"],
        loader_factory: Callable[[Any], Any],
        n_frames_total: int,
        timestep_ps: float,
    ) -> None:
        """Store replicate-specific Rg runner state.

        Parameters
        ----------
        sim_config : Any
            Simulation configuration used to build fresh universes.
        replicate : int
            Replicate number.
        runs : list[RgRunSettings]
            Configured Rg runs.
        loader_factory : Callable[[Any], Any]
            Factory used to create trajectory loaders.
        n_frames_total : int
            Total number of frames in the replicate trajectory.
        timestep_ps : float
            Trajectory timestep in picoseconds.
        """

        self._sim_config = sim_config
        self._replicate = replicate
        self._runs = list(runs)
        self._loader_factory = loader_factory
        self._n_frames_total = int(n_frames_total)
        self._timestep_ps = float(timestep_ps)
        self.results = RgRunnerPayload(
            n_frames_total=0,
            run_payloads=[],
            skipped_run_payloads=[],
        )

    def run(self, start: int, stop: int, step: int = 1) -> RgReplicateRunner:
        """Execute all configured Rg runs.

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
        RgReplicateRunner
            Runner instance with populated ``results``.
        """

        loader = self._loader_factory(self._sim_config)
        run_payloads: list[RgRunPayload] = []
        skipped_run_payloads: list[RgSkippedRunPayload] = []
        for run in self._runs:
            try:
                payload = compute_rg_run(
                    universe=loader.load_universe(self._replicate, cache=False),
                    run=run,
                    replicate=self._replicate,
                    start=start,
                    stop=stop,
                    step=step,
                    timestep_ps=self._timestep_ps,
                )
            except RgEmptySelectionError as exc:
                LOGGER.warning(
                    "Skipping Rg run '%s' for replicate %s because the selection matched no atoms",
                    exc.run_label,
                    exc.replicate,
                )
                skipped_run_payloads.append(
                    RgSkippedRunPayload(
                        run_label=exc.run_label,
                        selection=exc.selection,
                        replicate=exc.replicate,
                        reason=str(exc),
                    )
                )
                continue
            run_payloads.append(payload)
        self.results = RgRunnerPayload(
            n_frames_total=self._n_frames_total,
            run_payloads=run_payloads,
            skipped_run_payloads=skipped_run_payloads,
        )
        return self


def compute_rg_run(
    *,
    universe: Any,
    run: "RgRunSettings",
    replicate: int,
    start: int,
    stop: int,
    step: int,
    timestep_ps: float,
) -> RgRunPayload:
    """Compute one Rg run on a fresh universe.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe for the replicate.
    run : RgRunSettings
        Run configuration.
    replicate : int
        Replicate number used for warning messages.
    start : int
        Inclusive start frame.
    stop : int
        Exclusive stop frame.
    step : int
        Frame stride.
    timestep_ps : float
        Trajectory timestep in picoseconds.

    Returns
    -------
    RgRunPayload
        Trajectory-native payload for later result-model serialization.
    """

    from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

    atom_group = universe.select_atoms(run.selection)
    if len(atom_group) == 0:
        raise RgEmptySelectionError(
            run_label=run.label,
            selection=run.selection,
            replicate=replicate,
        )

    frame_indices = np.arange(start, stop, step, dtype=np.int64)
    frag_masses: NDArray[np.float64] | None = None
    frag_metadata: dict[str, float | int] = {}

    if run.calculation_mode == "fragments":
        fragments = list(atom_group.fragments) if atom_group.fragments else [atom_group]
        if len(fragments) == 1:
            LOGGER.warning(
                "Run '%s' selection produced only 1 fragment; fragment mode will behave "
                "identically to selection mode",
                run.label,
            )

        if run.fragment_weighting == "mass":
            frag_masses = np.asarray([frag.total_mass() for frag in fragments], dtype=np.float64)
            if np.any(frag_masses <= 0) or np.any(np.isnan(frag_masses)):
                raise ValueError(
                    f"Run '{run.label}': fragment masses contain zero or NaN values. "
                    "This suggests a problem with the MDAnalysis universe topology. "
                    f"Fragment masses: {frag_masses.tolist()}"
                )

        rg_values, fragment_counts_per_frame, fragment_rg_values = _run_fragment_rg_analysis(
            atom_group=atom_group,
            fragments=fragments,
            save_fragment_distribution=run.save_fragment_distribution,
            frag_masses=frag_masses,
            start=start,
            stop=stop,
            step=step,
        )

        if fragment_rg_values is not None and fragment_rg_values.size > 0:
            frag_metadata = {
                "fragment_mean_rg": float(np.mean(fragment_rg_values)),
                "fragment_std_rg": float(np.std(fragment_rg_values, ddof=0)),
                "fragment_median_rg": float(np.median(fragment_rg_values)),
                "fragment_min_rg": float(np.min(fragment_rg_values)),
                "fragment_max_rg": float(np.max(fragment_rg_values)),
                "fragment_rg_p10": float(np.percentile(fragment_rg_values, 10)),
                "fragment_rg_p25": float(np.percentile(fragment_rg_values, 25)),
                "fragment_rg_p50": float(np.percentile(fragment_rg_values, 50)),
                "fragment_rg_p75": float(np.percentile(fragment_rg_values, 75)),
                "fragment_rg_p90": float(np.percentile(fragment_rg_values, 90)),
            }

        frag_metadata["mean_fragments_per_frame"] = float(np.mean(fragment_counts_per_frame))
        frag_metadata["min_fragments_per_frame"] = int(np.min(fragment_counts_per_frame))
        frag_metadata["max_fragments_per_frame"] = int(np.max(fragment_counts_per_frame))
    else:
        rg_values = _run_selection_rg_analysis(
            atom_group=atom_group,
            start=start,
            stop=stop,
            step=step,
        )
        fragment_counts_per_frame = None
        fragment_rg_values = None

    time_ns = (frame_indices.astype(np.float64) * timestep_ps) / 1000.0
    effective_timestep_ps = float(timestep_ps) * float(step)

    mean_rg = float(np.mean(rg_values))
    std_rg = float(np.std(rg_values, ddof=0))
    median_rg = float(np.median(rg_values))
    min_rg = float(np.min(rg_values))
    max_rg = float(np.max(rg_values))
    final_rg = float(rg_values[-1])

    sem_rg: float | None = None
    correlation_time: float | None = None
    correlation_time_unit: str | None = None
    n_independent_frames: int | None = None
    statistical_inefficiency: float | None = None
    autocorrelation_warning: str | None = None

    if len(rg_values) >= 20:
        tau_result = estimate_correlation_time(
            rg_values,
            timestep=effective_timestep_ps,
            timestep_unit="ps",
            method="integration",
            n_frames=len(rg_values),
        )
        correlation_time = tau_result.tau
        correlation_time_unit = tau_result.tau_unit
        n_independent_frames = tau_result.n_independent
        statistical_inefficiency = tau_result.statistical_inefficiency
        autocorrelation_warning = tau_result.warning
        if n_independent_frames > 0:
            sem_rg = float(std_rg / np.sqrt(float(n_independent_frames)))

    return RgRunPayload(
        run_label=run.label,
        selection=run.selection,
        calculation_mode=run.calculation_mode,
        fragment_weighting=(
            run.fragment_weighting if run.calculation_mode == "fragments" else None
        ),
        rg_values=rg_values,
        frames=frame_indices,
        time_ns=time_ns,
        raw_timestep_ps=float(timestep_ps),
        frame_stride=int(step),
        effective_timestep_ps=effective_timestep_ps,
        mean_rg=mean_rg,
        std_rg=std_rg,
        median_rg=median_rg,
        min_rg=min_rg,
        max_rg=max_rg,
        final_rg=final_rg,
        sem_rg=sem_rg,
        correlation_time=correlation_time,
        correlation_time_unit=correlation_time_unit,
        n_independent_frames=n_independent_frames,
        statistical_inefficiency=statistical_inefficiency,
        autocorrelation_warning=autocorrelation_warning,
        fragment_rg_values=fragment_rg_values,
        fragment_counts_per_frame=fragment_counts_per_frame,
        fragment_masses=frag_masses,
        frag_metadata=frag_metadata,
    )


def _run_selection_rg_analysis(
    *,
    atom_group: Any,
    start: int,
    stop: int,
    step: int,
) -> NDArray[np.float64]:
    """Run selection-mode Rg through an MDAnalysis analysis object."""

    from MDAnalysis.analysis.base import AnalysisBase

    class _SelectionRgAnalysis(AnalysisBase):
        """Collect selection-mode Rg values for each trajectory frame."""

        def __init__(self, group: Any) -> None:
            super().__init__(group.universe.trajectory)
            self._group = group

        def _prepare(self) -> None:
            self.results.rg_values = []

        def _single_frame(self) -> None:
            self.results.rg_values.append(float(self._group.radius_of_gyration()))

        def _conclude(self) -> None:
            self.results.rg_values = np.asarray(self.results.rg_values, dtype=np.float64)

    analysis = _SelectionRgAnalysis(atom_group).run(start=start, stop=stop, step=step)
    return np.asarray(analysis.results.rg_values, dtype=np.float64)


def _run_fragment_rg_analysis(
    *,
    atom_group: Any,
    fragments: list[Any],
    save_fragment_distribution: bool,
    frag_masses: NDArray[np.float64] | None,
    start: int,
    stop: int,
    step: int,
) -> tuple[NDArray[np.float64], NDArray[np.int64], NDArray[np.float64] | None]:
    """Run fragment-mode Rg through an MDAnalysis analysis object."""

    from MDAnalysis.analysis.base import AnalysisBase

    class _FragmentRgAnalysis(AnalysisBase):
        """Collect fragment-mode Rg values for each trajectory frame."""

        def __init__(
            self,
            *,
            group: Any,
            fragment_groups: list[Any],
            record_distribution: bool,
            weights: NDArray[np.float64] | None,
        ) -> None:
            super().__init__(group.universe.trajectory)
            self._fragment_groups = fragment_groups
            self._record_distribution = record_distribution
            self._weights = weights

        def _prepare(self) -> None:
            self.results.rg_values = []
            self.results.fragment_counts = []
            self.results.fragment_rg_values = []

        def _single_frame(self) -> None:
            frag_rg = np.asarray(
                [fragment.radius_of_gyration() for fragment in self._fragment_groups],
                dtype=np.float64,
            )
            self.results.fragment_counts.append(len(self._fragment_groups))
            if self._record_distribution:
                self.results.fragment_rg_values.append(frag_rg.copy())
            reduced_rg = (
                float(np.average(frag_rg, weights=self._weights))
                if self._weights is not None
                else float(np.mean(frag_rg))
            )
            self.results.rg_values.append(reduced_rg)

        def _conclude(self) -> None:
            self.results.rg_values = np.asarray(self.results.rg_values, dtype=np.float64)
            self.results.fragment_counts = np.asarray(self.results.fragment_counts, dtype=np.int64)
            if self.results.fragment_rg_values:
                self.results.fragment_rg_values = np.concatenate(
                    self.results.fragment_rg_values,
                    axis=0,
                ).astype(np.float64)
            else:
                self.results.fragment_rg_values = None

    analysis = _FragmentRgAnalysis(
        group=atom_group,
        fragment_groups=fragments,
        record_distribution=save_fragment_distribution,
        weights=frag_masses,
    ).run(start=start, stop=stop, step=step)
    return (
        np.asarray(analysis.results.rg_values, dtype=np.float64),
        np.asarray(analysis.results.fragment_counts, dtype=np.int64),
        (
            None
            if analysis.results.fragment_rg_values is None
            else np.asarray(analysis.results.fragment_rg_values, dtype=np.float64)
        ),
    )
