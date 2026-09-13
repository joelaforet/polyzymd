"""Root mean square deviation, written against the observable contract.

One ``mean_of_timeseries`` observable per configured run, in angstrom. Each run
names a selection, an alignment selection and a reference structure, and the
observable name carries the reference mode, so a comparison table never shows
two runs measured against different references under one name.

The per-frame value is the minimised RMSD between the selection and the
reference structure, computed by MDAnalysis with the quaternion characteristic
polynomial solution of the Kabsch problem. MDAnalysis superimposes every frame
itself, so no separate alignment pass runs over the trajectory.

``alignment_selection`` chooses the atoms that superposition minimises over.
When it equals ``selection``, the reported value is the minimised RMSD of those
atoms. When it differs, each frame is superimposed on ``alignment_selection``
and the deviation is reported for ``selection``, which is how a flexible loop
or a bound ligand is measured against a rigid core. ``average`` mode is the one
mode that superimposes the trajectory in place, because a mean structure means
nothing until every frame shares a frame of reference.

Aggregation over replicates, uncertainty, cross-condition tests, persistence and
formatting belong to the framework.

References
----------
Kabsch, W. (1976). A solution for the best rotation to relate two sets of
vectors. *Acta Crystallographica Section A*, 32(5), 922-923.
doi:10.1107/S0567739476001873

Theobald, D. L. (2005). Rapid calculation of RMSDs using a quaternion-based
characteristic polynomial. *Acta Crystallographica Section A*, 61(4), 478-480.
doi:10.1107/S0108767305015266

Michaud-Agrawal, N., Denning, E. J., Woolf, T. B. & Beckstein, O. (2011).
MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
*Journal of Computational Chemistry*, 32(10), 2319-2327. doi:10.1002/jcc.21787
"""

from __future__ import annotations

import logging
import re
import warnings
from pathlib import Path
from typing import Any, ClassVar, Literal, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.contract import Observable, contract_analysis, iter_frames
from polyzymd.analyses.exceptions import ReplicateError, SelectionError
from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory

#: Settings of the dropped convergence heuristic, accepted for one more release.
DEPRECATED_SETTINGS = (
    "convergence_window_size_ns",
    "convergence_step_size_ns",
    "convergence_slope_threshold",
    "convergence_sustained_for_ns",
)


class RMSDRunSettings(BaseModel):
    """One named RMSD measurement.

    Parameters
    ----------
    label : str
        Run label. The observable is named ``rmsd_<label>_ref_<reference_mode>``.
    selection : str
        Atoms whose deviation is measured.
    alignment_selection : str
        Atoms used to superimpose the trajectory before the reference is built.
    reference_mode : {"centroid", "average", "frame", "external"}
        Structure the deviation is measured from. ``centroid`` is the frame
        closest to the aligned mean, ``average`` the aligned mean itself,
        ``frame`` a chosen frame, ``external`` a PDB file. Only ``external``
        compares across conditions, because every other mode gives each
        replicate its own reference.
    reference_frame : int
        Frame index, 0-indexed, used when ``reference_mode="frame"``.
    reference_file : Path or None
        PDB file used when ``reference_mode="external"``.
    centroid_selection : str or None
        Atoms used to find the representative frame in ``centroid`` mode.
        Defaults to ``alignment_selection``.
    """

    label: str = Field(min_length=1, description="Run label")
    selection: str = Field(
        default="protein and name CA", description="Atoms whose deviation is measured"
    )
    alignment_selection: str = Field(
        default="protein and name CA", description="Atoms used for superposition"
    )
    reference_mode: Literal["centroid", "average", "frame", "external"] = Field(
        default="centroid", description="Reference structure mode"
    )
    reference_frame: int = Field(default=0, ge=0, description="0-indexed reference frame")
    reference_file: Path | None = Field(default=None, description="External reference PDB")
    centroid_selection: str | None = Field(
        default=None, description="Atoms used to find the representative frame"
    )

    @model_validator(mode="before")
    @classmethod
    def _warn_on_dropped_settings(cls, data: Any) -> Any:
        """Accept the convergence settings of the old plugin once more, with a warning."""
        if not isinstance(data, dict):
            return data
        dropped = [key for key in DEPRECATED_SETTINGS if key in data]
        if not dropped:
            return data
        message = (
            f"rmsd: {', '.join(dropped)} no longer has an effect and will be rejected in the "
            "next release. The sliding-window convergence flag was removed because its slope "
            "threshold sat below the noise of successive window means."
        )
        warnings.warn(message, UserWarning, stacklevel=2)
        logging.getLogger("polyzymd.analyses").warning(message)
        return {key: value for key, value in data.items() if key not in dropped}

    @model_validator(mode="after")
    def _check_external_reference(self) -> RMSDRunSettings:
        """Require a reference file in external mode, and warn if it is not here.

        The file is read on the machine that runs the trajectory, so a missing
        path is a warning at parse time and an error at compute time. That way a
        comparison file naming a cluster path still validates on a laptop.
        """
        if self.reference_mode != "external":
            return self
        if self.reference_file is None:
            raise ValueError(f"run {self.label!r}: reference_mode='external' needs reference_file")
        if not Path(self.reference_file).exists():
            warnings.warn(
                f"rmsd run {self.label!r}: reference_file {self.reference_file} is not on this "
                "machine; it must exist where the analysis runs",
                UserWarning,
                stacklevel=2,
            )
        return self


class RMSDSettings(BaseModel):
    """Runs measured by the RMSD analysis."""

    runs: list[RMSDRunSettings] = Field(min_length=1, description="RMSD runs to compute")

    @field_validator("runs", mode="after")
    @classmethod
    def _unique_observable_names(cls, runs: list[RMSDRunSettings]) -> list[RMSDRunSettings]:
        """Reject two runs that would write one observable name.

        Two labels that differ only in case or punctuation, such as ``Core
        Frame`` and ``core_frame``, reach the same slug, and the second would
        then overwrite the first in the sidecar.
        """
        names = [observable_name(run) for run in runs]
        if len(set(names)) != len(names):
            raise ValueError(
                "rmsd run labels must give unique observable names, got "
                f"{[run.label for run in runs]} for {names}"
            )
        return runs


class RMSD:
    """Deviation of a selection from a reference structure, per frame."""

    name: ClassVar[str] = "rmsd"
    Settings: ClassVar[type[BaseModel]] = RMSDSettings
    references: ClassVar[tuple[str, ...]] = (
        "Kabsch 1976, Acta Cryst A32:922, doi:10.1107/S0567739476001873",
        "Theobald 2005, Acta Cryst A61:478, doi:10.1107/S0108767305015266",
        "Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787",
    )

    def compute(self, universe: Any, frames: Any, settings: RMSDSettings) -> Sequence[Observable]:
        """Measure every configured run on one replicate.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework. Its coordinates are superimposed
            in place, which leaves every RMSD value unchanged because the
            measurement itself minimises over rigid motions.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : RMSDSettings
            Runs to measure.

        Returns
        -------
        Sequence[Observable]
            One ``mean_of_timeseries`` observable per run, in angstrom.
        """
        return [_measure(universe, frames, run) for run in settings.runs]


def observable_name(run: RMSDRunSettings) -> str:
    """Name of the observable one run reports, reference mode included.

    Parameters
    ----------
    run : RMSDRunSettings
        Run whose observable name is wanted.

    Returns
    -------
    str
        ``rmsd_<label slug>_ref_<reference mode>``.
    """
    slug = re.sub(r"[^a-z0-9]+", "_", run.label.lower()).strip("_")
    return f"rmsd_{slug}_ref_{run.reference_mode}"


def _measure(universe: Any, frames: Any, run: RMSDRunSettings) -> Observable:
    """Per-frame RMSD of one run against its reference structure."""
    from MDAnalysis.analysis.rms import RMSD as MDAnalysisRMSD

    start, stop, step = _window(universe, frames)
    scope = _scope(universe, run)
    reference = _reference_universe(universe, frames, run, scope, start, stop)
    on_core = run.alignment_selection != run.selection
    analysis = MDAnalysisRMSD(
        universe,
        reference,
        select=run.alignment_selection if on_core else run.selection,
        groupselections=[run.selection] if on_core else None,
    )
    analysis.run(start=start, stop=stop, step=step)
    return Observable(
        name=observable_name(run),
        kind="mean_of_timeseries",
        unit="A",
        values=np.asarray(analysis.results.rmsd[:, 3 if on_core else 2], dtype=np.float64),
        higher_is_better=False,
    )


def _scope(universe: Any, run: RMSDRunSettings) -> Any:
    """Atoms the reference structure must carry, measured and superposed alike.

    Also rejects a superposition group of fewer than three atoms, which leaves
    the rotation undetermined and makes MDAnalysis return NaN.
    """
    measured = _select(universe, run.selection, run.label)
    superposed = measured
    if run.alignment_selection != run.selection:
        superposed = _select(universe, run.alignment_selection, run.label)
        measured = measured | superposed
    if len(superposed) < 3:
        raise SelectionError(
            f"rmsd run {run.label!r}: superposition needs at least three atoms, but "
            f"{run.alignment_selection!r} matches {len(superposed)}"
        )
    return measured


def _reference_universe(
    universe: Any, frames: Any, run: RMSDRunSettings, scope: Any, start: int, stop: int
) -> Any:
    """Universe holding the reference structure of one run in a single frame."""
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    if run.reference_mode == "external":
        return _external_universe(universe, run)
    if run.reference_mode == "average":
        positions = _average_positions(universe, frames, run, scope, start, stop)
    else:
        universe.trajectory[_reference_frame(universe, run, start, stop)]
        positions = scope.positions.astype(np.float64)
    reference = mda.Merge(scope)
    reference.load_new(positions[np.newaxis, :, :], format=MemoryReader)
    return reference


def _reference_frame(universe: Any, run: RMSDRunSettings, start: int, stop: int) -> int:
    """Index of the frame that ``frame`` and ``centroid`` mode measure against."""
    if run.reference_mode == "centroid":
        from polyzymd.analyses.shared.centroid import find_centroid_frame

        return int(
            find_centroid_frame(
                universe,
                selection=run.centroid_selection or run.alignment_selection,
                start_frame=start,
                stop_frame=stop,
                verbose=False,
            )
        )
    if not 0 <= run.reference_frame < len(universe.trajectory):
        raise ReplicateError(
            f"rmsd run {run.label!r}: reference_frame {run.reference_frame} is outside the "
            f"trajectory, which holds {len(universe.trajectory)} frames"
        )
    return int(run.reference_frame)


def _external_universe(universe: Any, run: RMSDRunSettings) -> Any:
    """Universe of the external PDB, checked against the trajectory atom by atom."""
    import MDAnalysis as mda

    path = Path(str(run.reference_file))
    if not path.exists():
        raise ReplicateError(
            f"rmsd run {run.label!r}: reference_file {path} does not exist on this machine"
        )
    reference = mda.Universe(str(path))
    for selection in dict.fromkeys((run.alignment_selection, run.selection)):
        atoms = _select(reference, selection, run.label, source=str(path))
        expected = len(universe.select_atoms(selection))
        if len(atoms) != expected:
            raise SelectionError(
                f"rmsd run {run.label!r}: selection {selection!r} matches {expected} atoms in "
                f"the trajectory but {len(atoms)} in {path}"
            )
    return reference


def _average_positions(
    universe: Any, frames: Any, run: RMSDRunSettings, scope: Any, start: int, stop: int
) -> np.ndarray:
    """Mean position of every atom in ``scope``, taken in a common frame of reference.

    This is the one place the trajectory is superimposed in place. A mean taken
    over frames that still carry rigid-body drift is not a structure.
    """
    align_trajectory(
        universe,
        AlignmentConfig(enabled=True, reference_mode="average", selection=run.alignment_selection),
        start_frame=start,
        stop_frame=stop,
        step_frame=int(frames.step or 1),
    )
    total = np.zeros_like(scope.positions, dtype=np.float64)
    n_frames = 0
    for _ in iter_frames(universe, frames):
        total += scope.positions
        n_frames += 1
    if n_frames == 0:
        raise ReplicateError(
            f"rmsd run {run.label!r}: the production window holds no frames, so there is no "
            "average structure to measure against"
        )
    return total / float(n_frames)


def _select(universe: Any, selection: str, label: str, source: str = "the trajectory") -> Any:
    """Select atoms, or say which run and which structure came up empty."""
    group = universe.select_atoms(selection)
    if len(group) == 0:
        raise SelectionError(
            f"rmsd run {label!r}: selection {selection!r} matched no atoms in {source}"
        )
    return group


def _window(universe: Any, frames: Any) -> tuple[int, int, int]:
    """Start, stop and step of the production window."""
    if frames.frames is not None:
        raise ReplicateError(
            "rmsd needs a contiguous production window because its reference structure is "
            "defined on one, but the frame selection lists explicit frames"
        )
    n_total = frames.n_frames_total or len(universe.trajectory)
    return (
        int(frames.start or 0),
        int(n_total if frames.stop is None else frames.stop),
        int(frames.step or 1),
    )


RMSDAnalysis = contract_analysis(RMSD)
