"""Runner-backed RMSD execution helpers.

This module keeps the trajectory-native RMSD work behind the runner seam so
the plugin can preserve its legacy cache layout while PolyzyMD owns the outer
replicate workflow.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory

if TYPE_CHECKING:
    from polyzymd.analyses.rmsd import RMSDRunSettings

LOGGER = logging.getLogger(__name__)


@dataclass
class RMSDRunPayload:
    """Trajectory-native RMSD data for one configured run.

    Parameters
    ----------
    run_label : str
        Human-readable run label.
    selection : str
        MDAnalysis selection used for RMSD calculation.
    alignment_selection : str
        Selection used for trajectory alignment.
    reference_mode : str
        Reference-structure mode used for alignment and RMSD.
    reference_frame : int | None
        1-indexed reference frame used for the run, when applicable.
    reference_file : str | None
        External reference path for ``reference_mode="external"``.
    rmsd_values : NDArray[np.float64]
        Per-frame RMSD values.
    frames : NDArray[np.int64]
        Absolute frame indices analyzed.
    time_ns : NDArray[np.float64]
        Time axis in nanoseconds.
    mean_rmsd : float
        Mean RMSD over analyzed frames.
    std_rmsd : float
        Standard deviation of RMSD values.
    median_rmsd : float
        Median RMSD.
    min_rmsd : float
        Minimum RMSD.
    max_rmsd : float
        Maximum RMSD.
    final_rmsd : float
        Final RMSD value.
    sem_rmsd : float | None
        Autocorrelation-corrected SEM when available.
    correlation_time : float | None
        Estimated correlation time.
    correlation_time_unit : str | None
        Unit for correlation time.
    n_independent_frames : int | None
        Effective independent frame count.
    statistical_inefficiency : float | None
        Statistical inefficiency estimate.
    autocorrelation_warning : str | None
        Reliability warning from autocorrelation analysis.
    convergence_result : Any
        Sliding-window convergence result object.
    """

    run_label: str
    selection: str
    alignment_selection: str
    reference_mode: str
    reference_frame: int | None
    reference_file: str | None
    rmsd_values: NDArray[np.float64]
    frames: NDArray[np.int64]
    time_ns: NDArray[np.float64]
    mean_rmsd: float
    std_rmsd: float
    median_rmsd: float
    min_rmsd: float
    max_rmsd: float
    final_rmsd: float
    sem_rmsd: float | None
    correlation_time: float | None
    correlation_time_unit: str | None
    n_independent_frames: int | None
    statistical_inefficiency: float | None
    autocorrelation_warning: str | None
    convergence_result: Any


@dataclass
class RMSDRunnerPayload:
    """Collection of RMSD run payloads for one replicate."""

    n_frames_total: int
    run_payloads: list[RMSDRunPayload]


class RMSDReplicateRunner:
    """Execute RMSD runs for one replicate through the runner seam."""

    def __init__(
        self,
        *,
        sim_config: Any,
        replicate: int,
        runs: list["RMSDRunSettings"],
        loader_factory: Callable[[Any], Any],
        n_frames_total: int,
        timestep_ps: float,
    ) -> None:
        """Store replicate-specific RMSD runner state.

        Parameters
        ----------
        sim_config : Any
            Simulation configuration used to materialize fresh universes.
        replicate : int
            Replicate number.
        runs : list[RMSDRunSettings]
            Configured RMSD runs.
        loader_factory : Callable[[Any], Any]
            Factory used to build trajectory loaders.
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
        self.results = RMSDRunnerPayload(n_frames_total=0, run_payloads=[])

    def run(self, start: int, stop: int, step: int = 1) -> RMSDReplicateRunner:
        """Execute all configured RMSD runs.

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
        RMSDReplicateRunner
            Runner instance with populated ``results``.
        """

        loader = self._loader_factory(self._sim_config)
        run_payloads: list[RMSDRunPayload] = []
        for run in self._runs:
            payload = compute_rmsd_run(
                universe=loader.load_universe(self._replicate, cache=False),
                run=run,
                start=start,
                stop=stop,
                step=step,
                timestep_ps=self._timestep_ps,
            )
            run_payloads.append(payload)
        self.results = RMSDRunnerPayload(
            n_frames_total=self._n_frames_total,
            run_payloads=run_payloads,
        )
        return self


def compute_rmsd_run(
    *,
    universe: Any,
    run: "RMSDRunSettings",
    start: int,
    stop: int,
    step: int,
    timestep_ps: float,
) -> RMSDRunPayload:
    """Compute one RMSD run on a fresh universe.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe for the replicate.
    run : RMSDRunSettings
        Run configuration.
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
    RMSDRunPayload
        Trajectory-native payload for later result-model serialization.
    """

    from MDAnalysis.analysis.rms import RMSD

    from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time
    from polyzymd.analyses.shared.convergence import find_convergence_time

    centroid_selection = run.centroid_selection
    if run.reference_mode == "centroid" and centroid_selection is None:
        centroid_selection = run.alignment_selection
        LOGGER.info(
            "Run '%s': centroid_selection not set, using alignment_selection='%s'",
            run.label,
            centroid_selection,
        )

    reference_frame_1indexed: int | None
    if run.reference_mode == "frame":
        reference_frame_1indexed = run.reference_frame + 1
    else:
        reference_frame_1indexed = None

    alignment_config = AlignmentConfig(
        enabled=True,
        reference_mode=run.reference_mode,
        reference_frame=reference_frame_1indexed,
        selection=run.alignment_selection,
        centroid_selection=centroid_selection or run.alignment_selection,
        reference_file=(Path(run.reference_file) if run.reference_file is not None else None),
    )
    ref_frame_idx = align_trajectory(
        universe,
        alignment_config,
        start_frame=start,
        stop_frame=stop,
    )

    atom_group = universe.select_atoms(run.selection)
    if len(atom_group) == 0:
        raise ValueError(f"Run '{run.label}' selection matched no atoms: {run.selection!r}")

    _reference_universe, reference_atom_group = _build_reference_structure(
        universe=universe,
        atom_group=atom_group,
        run=run,
        start=start,
        stop=stop,
        step=step,
        ref_frame_idx=ref_frame_idx,
    )

    rmsd_analysis = RMSD(
        atom_group,
        reference=reference_atom_group,
        select="all",
        ref_frame=0,
    ).run(start=start, stop=stop, step=step)

    rmsd_values = np.asarray(rmsd_analysis.results.rmsd[:, 2], dtype=np.float64)
    frames = np.arange(start, stop, step, dtype=np.int64)
    time_ns = (frames.astype(np.float64) * timestep_ps) / 1000.0

    mean_rmsd = float(np.mean(rmsd_values))
    std_rmsd = float(np.std(rmsd_values, ddof=0))
    median_rmsd = float(np.median(rmsd_values))
    min_rmsd = float(np.min(rmsd_values))
    max_rmsd = float(np.max(rmsd_values))
    final_rmsd = float(rmsd_values[-1])

    sem_rmsd: float | None = None
    correlation_time: float | None = None
    correlation_time_unit: str | None = None
    n_independent_frames: int | None = None
    statistical_inefficiency: float | None = None
    autocorrelation_warning: str | None = None

    if len(rmsd_values) >= 20:
        tau_result = estimate_correlation_time(
            rmsd_values,
            timestep=timestep_ps,
            timestep_unit="ps",
            method="integration",
            n_frames=len(rmsd_values),
        )
        correlation_time = tau_result.tau
        correlation_time_unit = tau_result.tau_unit
        n_independent_frames = tau_result.n_independent
        statistical_inefficiency = tau_result.statistical_inefficiency
        autocorrelation_warning = tau_result.warning
        if n_independent_frames > 0:
            sem_rmsd = float(std_rmsd / np.sqrt(float(n_independent_frames)))

    convergence_result = find_convergence_time(
        time_ns,
        rmsd_values,
        window_size_ns=run.convergence_window_size_ns,
        step_size_ns=run.convergence_step_size_ns,
        slope_threshold=run.convergence_slope_threshold,
        sustained_for_ns=run.convergence_sustained_for_ns,
    )

    return RMSDRunPayload(
        run_label=run.label,
        selection=run.selection,
        alignment_selection=run.alignment_selection,
        reference_mode=run.reference_mode,
        reference_frame=(ref_frame_idx + 1 if ref_frame_idx is not None else None),
        reference_file=run.reference_file,
        rmsd_values=rmsd_values,
        frames=frames,
        time_ns=time_ns,
        mean_rmsd=mean_rmsd,
        std_rmsd=std_rmsd,
        median_rmsd=median_rmsd,
        min_rmsd=min_rmsd,
        max_rmsd=max_rmsd,
        final_rmsd=final_rmsd,
        sem_rmsd=sem_rmsd,
        correlation_time=correlation_time,
        correlation_time_unit=correlation_time_unit,
        n_independent_frames=n_independent_frames,
        statistical_inefficiency=statistical_inefficiency,
        autocorrelation_warning=autocorrelation_warning,
        convergence_result=convergence_result,
    )


def _build_reference_structure(
    *,
    universe: Any,
    atom_group: Any,
    run: "RMSDRunSettings",
    start: int,
    stop: int,
    step: int,
    ref_frame_idx: int | None,
) -> tuple[Any, Any]:
    """Materialize the reference structure for RMSD calculations.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe.
    atom_group : Any
        Selected atoms to compare.
    run : RMSDRunSettings
        Run configuration.
    start : int
        Inclusive start frame.
    stop : int
        Exclusive stop frame.
    step : int
        Frame stride.
    ref_frame_idx : int | None
        Reference frame returned by the alignment helper.

    Returns
    -------
    tuple[Any, Any]
        Reference universe and matching atom group.
    """

    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    if run.reference_mode in {"centroid", "frame"}:
        if ref_frame_idx is None:
            raise ValueError(
                f"Run '{run.label}' expected a reference frame for mode '{run.reference_mode}'"
            )
        universe.trajectory[ref_frame_idx]
        ref_positions = atom_group.positions.copy().astype(np.float64)
    elif run.reference_mode == "average":
        mean_analysis = _build_mean_position_analysis(atom_group).run(
            start=start,
            stop=stop,
            step=step,
        )
        ref_positions = np.asarray(mean_analysis.results.mean_positions, dtype=np.float64)
    elif run.reference_mode == "external":
        if run.reference_file is None:
            raise ValueError(
                f"Run '{run.label}' requires reference_file when reference_mode='external'"
            )

        ref_path = Path(run.reference_file)
        LOGGER.info("Using external reference from: %s", ref_path)

        ref_universe = mda.Universe(str(ref_path))
        ref_atoms = ref_universe.select_atoms(run.selection)

        if len(ref_atoms) == 0:
            raise ValueError(
                f"Run '{run.label}' external PDB '{ref_path.name}' has no atoms matching "
                f"selection {run.selection!r}."
            )

        if len(ref_atoms) != len(atom_group):
            LOGGER.warning(
                "Run '%s' external reference atom count mismatch for selection %r "
                "(trajectory=%d, external=%d)",
                run.label,
                run.selection,
                len(atom_group),
                len(ref_atoms),
            )
            raise ValueError(
                f"Run '{run.label}' atom count mismatch between trajectory "
                f"({len(atom_group)}) and external PDB ({len(ref_atoms)}) for "
                f"selection {run.selection!r}."
            )

        ref_positions = ref_atoms.positions.copy().astype(np.float64)
    else:
        raise ValueError(f"Unsupported RMSD reference_mode: {run.reference_mode!r}")

    reference_universe = mda.Merge(atom_group)
    reference_universe.load_new(ref_positions[np.newaxis, :, :], format=MemoryReader)
    reference_atom_group = reference_universe.atoms
    return reference_universe, reference_atom_group


def _build_mean_position_analysis(atoms: Any) -> Any:
    """Build an MDAnalysis analysis object for mean-position accumulation."""

    from MDAnalysis.analysis.base import AnalysisBase

    class _MeanPositionAnalysis(AnalysisBase):
        """Accumulate average coordinates over the analyzed frames."""

        def __init__(self, atom_group: Any) -> None:
            super().__init__(atom_group.universe.trajectory)
            self._atom_group = atom_group

        def _prepare(self) -> None:
            self._position_sum = np.zeros_like(self._atom_group.positions, dtype=np.float64)
            self._n_frames = 0

        def _single_frame(self) -> None:
            self._position_sum += self._atom_group.positions.astype(np.float64)
            self._n_frames += 1

        def _conclude(self) -> None:
            if self._n_frames == 0:
                raise ValueError("Mean-position analysis requires at least one frame")
            self.results.mean_positions = self._position_sum / float(self._n_frames)

    return _MeanPositionAnalysis(atoms)
