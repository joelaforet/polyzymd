"""Root mean square deviation, written against the observable contract.

One ``mean_of_timeseries`` observable per configured run, in angstrom. Each run
names a selection, an alignment selection and a reference structure, and the
observable name carries the reference mode, so a comparison table never shows
two runs measured against different references under one name.

The per-frame value is the minimised RMSD between the selection and the
reference structure, computed by MDAnalysis with the quaternion characteristic
polynomial solution of the Kabsch problem. The trajectory is superimposed on
``alignment_selection`` first because the centroid and average references are
defined in aligned space.

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

import re
import warnings
from pathlib import Path
from typing import Any, ClassVar, Literal, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.contract import Observable, iter_frames
from polyzymd.analyses.contract_runner import contract_analysis
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
        warnings.warn(
            f"rmsd: {', '.join(dropped)} no longer has an effect and will be rejected in the "
            "next release. The sliding-window convergence flag was removed because its slope "
            "threshold sat below the noise of successive window means.",
            DeprecationWarning,
            stacklevel=2,
        )
        return {key: value for key, value in data.items() if key not in dropped}

    @model_validator(mode="after")
    def _check_external_reference(self) -> RMSDRunSettings:
        """Require an existing PDB file in external mode."""
        if self.reference_mode != "external":
            return self
        if self.reference_file is None or not Path(self.reference_file).exists():
            raise ValueError(
                f"run {self.label!r}: reference_mode='external' needs reference_file to name an "
                f"existing PDB file, got {self.reference_file}"
            )
        return self


class RMSDSettings(BaseModel):
    """Runs measured by the RMSD analysis."""

    runs: list[RMSDRunSettings] = Field(min_length=1, description="RMSD runs to compute")

    @field_validator("runs", mode="after")
    @classmethod
    def _unique_labels(cls, runs: list[RMSDRunSettings]) -> list[RMSDRunSettings]:
        """Reject two runs that would write one observable name."""
        labels = [run.label for run in runs]
        if len(set(labels)) != len(labels):
            raise ValueError(f"rmsd run labels must be unique, got {labels}")
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
    reference_frame = align_trajectory(
        universe,
        _alignment_config(run),
        start_frame=start,
        stop_frame=stop,
        step_frame=step,
    )
    group = _select(universe, run.selection, run.label)
    reference = _reference_atoms(universe, frames, run, group, reference_frame)
    analysis = MDAnalysisRMSD(group, reference=reference, select="all", ref_frame=0)
    analysis.run(start=start, stop=stop, step=step)
    return Observable(
        name=observable_name(run),
        kind="mean_of_timeseries",
        unit="A",
        values=np.asarray(analysis.results.rmsd[:, 2], dtype=np.float64),
        higher_is_better=False,
    )


def _alignment_config(run: RMSDRunSettings) -> AlignmentConfig:
    """Alignment configuration for one run, with the 1-indexed frame the helper expects."""
    return AlignmentConfig(
        enabled=True,
        reference_mode=run.reference_mode,
        reference_frame=run.reference_frame + 1 if run.reference_mode == "frame" else None,
        selection=run.alignment_selection,
        centroid_selection=run.centroid_selection or run.alignment_selection,
        reference_file=run.reference_file,
    )


def _reference_atoms(
    universe: Any, frames: Any, run: RMSDRunSettings, group: Any, reference_frame: int | None
) -> Any:
    """Materialise the reference structure as a one-frame atom group."""
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    if run.reference_mode == "external":
        positions = _external_positions(run, group)
    elif run.reference_mode == "average":
        positions = _average_positions(universe, frames, group)
    else:
        if reference_frame is None:
            raise ReplicateError(
                f"rmsd run {run.label!r}: alignment returned no reference frame for mode "
                f"{run.reference_mode!r}"
            )
        universe.trajectory[reference_frame]
        positions = group.positions.astype(np.float64)
    reference_universe = mda.Merge(group)
    reference_universe.load_new(positions[np.newaxis, :, :], format=MemoryReader)
    return reference_universe.atoms


def _external_positions(run: RMSDRunSettings, group: Any) -> np.ndarray:
    """Reference positions read from the external PDB file of one run."""
    import MDAnalysis as mda

    reference = mda.Universe(str(run.reference_file))
    atoms = _select(reference, run.selection, run.label, source=str(run.reference_file))
    if len(atoms) != len(group):
        raise SelectionError(
            f"rmsd run {run.label!r}: selection {run.selection!r} matches {len(group)} atoms in "
            f"the trajectory but {len(atoms)} in {run.reference_file}"
        )
    return atoms.positions.astype(np.float64)


def _average_positions(universe: Any, frames: Any, group: Any) -> np.ndarray:
    """Mean position of every selected atom over the aligned production window."""
    total = np.zeros_like(group.positions, dtype=np.float64)
    n_frames = 0
    for _ in iter_frames(universe, frames):
        total += group.positions
        n_frames += 1
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
