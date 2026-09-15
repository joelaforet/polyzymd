"""Per-residue root mean square fluctuation, written against the observable contract.

The plugin reports one ``profile`` observable, ``rmsf``, holding the fluctuation
of every selected residue about the mean structure of the aligned production
window. Every frame in the window is used. The mean over residues is reported as
the scalar ``rmsf_mean`` through the profile's ``mean_over_index`` reduction, and
that scalar is what the framework tests across conditions.

In ``external`` reference mode the trajectory is superposed on a structure that
is not drawn from the simulation, so a second profile,
``rmsd_about_reference_per_residue``, reports the per-residue deviation from that
structure. The two answer different questions: ``rmsf`` measures spread about the
trajectory's own mean, the second measures distance from the external model. The
pre-port plugin stored the second under the name of the first.

References
----------
Michaud-Agrawal, N., Denning, E. J., Woolf, T. B. & Beckstein, O. (2011).
MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
*Journal of Computational Chemistry*, 32(10), 2319-2327. doi:10.1002/jcc.21787

Kuzmanic, A. & Zagrovic, B. (2010). Determination of ensemble-average pairwise
root mean-square deviation from experimental B-factors. *Biophysical Journal*,
98(5), 861-871. doi:10.1016/j.bpj.2009.11.011

Grossfield, A., Patrone, P. N., Roe, D. R., Schultz, A. J., Siderius, D. W. &
Zuckerman, D. M. (2018). Best practices for quantifying the uncertainty in
molecular simulations. *Living Journal of Computational Molecular Science*,
1(1), 5067. doi:10.33011/livecoms.1.1.5067
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar, Literal, Sequence

import numpy as np
from pydantic import BaseModel, Field, model_validator

from polyzymd.analyses.contract import Observable, contract_analysis, iter_frames
from polyzymd.analyses.exceptions import PluginContractError, ReplicateError, SelectionError
from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory

ReferenceMode = Literal["centroid", "average", "frame", "external"]


class RMSFSettings(BaseModel):
    """Settings for the per-residue RMSF analysis."""

    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection whose residues carry the profile",
    )
    alignment_selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection superposed before the fluctuation is measured",
    )
    centroid_selection: str = Field(
        default="protein",
        description="MDAnalysis selection used to pick the representative frame in centroid mode",
    )
    reference_mode: ReferenceMode = Field(
        default="centroid",
        description="Alignment reference: centroid, average, frame or external",
    )
    reference_frame: int | None = Field(
        default=None,
        description="One-indexed frame used when reference_mode is 'frame'",
    )
    reference_file: str | None = Field(
        default=None,
        description="Structure file used when reference_mode is 'external'",
    )

    @model_validator(mode="after")
    def _check_reference(self) -> RMSFSettings:
        """Reject a reference mode whose input is missing."""
        if self.reference_mode == "frame" and self.reference_frame is None:
            raise ValueError("reference_frame is required when reference_mode is 'frame'")
        if self.reference_mode == "external" and self.reference_file is None:
            raise ValueError("reference_file is required when reference_mode is 'external'")
        return self


class RMSF:
    """Per-residue fluctuation about the mean structure of the aligned window."""

    name: ClassVar[str] = "rmsf"
    Settings: ClassVar[type[BaseModel]] = RMSFSettings
    references: ClassVar[tuple[str, ...]] = (
        "Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787",
        "Kuzmanic and Zagrovic 2010, Biophys J 98:861, doi:10.1016/j.bpj.2009.11.011",
        "Grossfield et al. 2018, LiveCoMS 1:5067, doi:10.33011/livecoms.1.1.5067",
    )

    @staticmethod
    def identity_files(settings: RMSFSettings) -> Sequence[Path]:
        """External reference structure, whose contents change the answer."""
        if settings.reference_mode != "external" or settings.reference_file is None:
            return ()
        return (Path(settings.reference_file).expanduser(),)

    def compute(self, universe: Any, frames: Any, settings: RMSFSettings) -> Sequence[Observable]:
        """Measure the per-residue fluctuation of one replicate.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework. Alignment rewrites its
            coordinates in memory.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : RMSFSettings
            Selections and alignment reference.

        Returns
        -------
        Sequence[Observable]
            The ``rmsf`` profile in angstrom, plus
            ``rmsd_about_reference_per_residue`` in external mode.

        Raises
        ------
        SelectionError
            If the RMSF selection matches no atoms, or an external reference is
            missing or does not match the selected atoms.
        PluginContractError
            If an explicit frame list is given in a reference mode that builds
            its reference from a contiguous slice.
        ReplicateError
            If the production window holds no frames.
        """
        atoms = universe.select_atoms(settings.selection)
        if len(atoms) == 0:
            raise SelectionError(f"rmsf: selection {settings.selection!r} matched no atoms")
        reference_file = _reference_file(settings)
        start, stop, step = _window(universe, frames, settings)
        align_trajectory(
            universe,
            AlignmentConfig(
                enabled=True,
                reference_mode=settings.reference_mode,
                reference_frame=settings.reference_frame,
                selection=settings.alignment_selection,
                centroid_selection=settings.centroid_selection,
                reference_file=reference_file,
            ),
            start_frame=start,
            stop_frame=stop,
            step_frame=step,
        )
        reference = _reference_positions(atoms, settings, reference_file)

        n_frames = 0
        mean = np.zeros((len(atoms), 3), dtype=np.float64)
        sum_squares = np.zeros((len(atoms), 3), dtype=np.float64)
        about_reference = np.zeros(len(atoms), dtype=np.float64)
        for _ in iter_frames(universe, frames):
            positions = atoms.positions.astype(np.float64)
            n_frames += 1
            delta = positions - mean
            mean += delta / float(n_frames)
            sum_squares += delta * (positions - mean)
            if reference is not None:
                about_reference += np.sum((positions - reference) ** 2, axis=1)
        if n_frames == 0:
            raise ReplicateError("rmsf: the production window holds no frames")

        observables = [
            _profile(
                atoms,
                "rmsf",
                np.sqrt(np.sum(sum_squares / n_frames, axis=1)),
                n_frames,
                reduced_kind="fluctuation",
            )
        ]
        if reference is not None:
            observables.append(
                _profile(
                    atoms,
                    "rmsd_about_reference_per_residue",
                    np.sqrt(about_reference / n_frames),
                    n_frames,
                )
            )
        return observables


def _profile(
    atoms: Any,
    name: str,
    per_atom: np.ndarray,
    n_frames: int,
    *,
    reduced_kind: str | None = None,
) -> Observable:
    """Average per-atom values inside each residue and label them by residue ID."""
    _, inverse = np.unique(np.asarray(atoms.resindices), return_inverse=True)
    per_residue = np.bincount(inverse, weights=per_atom) / np.bincount(inverse)
    return Observable(
        name=name,
        kind="profile",
        unit="A",
        values=per_residue,
        index=np.asarray(atoms.residues.resids, dtype=np.float64),
        higher_is_better=False,
        n_frames=n_frames,
        reduce="mean_over_index",
        reduced_kind=reduced_kind,
    )


def _window(universe: Any, frames: Any, settings: RMSFSettings) -> tuple[int, int, int]:
    """Contiguous frame bounds the alignment reference is built from.

    Raises
    ------
    PluginContractError
        If the framework passed an explicit frame list in a reference mode that
        builds its reference from a contiguous slice.
    """
    if frames.frames is not None:
        if settings.reference_mode in {"centroid", "average"}:
            raise PluginContractError(
                f"rmsf: reference_mode={settings.reference_mode!r} builds its reference from a "
                "contiguous trajectory slice, so it cannot be used with an explicit frame list. "
                "Use reference_mode 'frame' or 'external', or a start/stop/step window."
            )
        selected = np.asarray(list(frames.frames))
        if selected.dtype == bool:
            selected = np.flatnonzero(selected)
        return int(selected.min()), int(selected.max()) + 1, 1
    stop = len(universe.trajectory) if frames.stop is None else int(frames.stop)
    return int(frames.start or 0), stop, int(frames.step or 1)


def _reference_file(settings: RMSFSettings) -> Path | None:
    """Resolved external reference path, checked before the trajectory is touched."""
    if settings.reference_mode != "external" or settings.reference_file is None:
        return None
    path = Path(settings.reference_file).expanduser()
    if not path.exists():
        raise SelectionError(
            f"rmsf: reference_file {path} does not exist; reference_mode is 'external', so "
            "the analysis needs the structure it measures deviation from"
        )
    return path


def _reference_positions(
    atoms: Any, settings: RMSFSettings, reference_file: Path | None
) -> np.ndarray | None:
    """Positions of the external reference structure, or None in the other modes."""
    if reference_file is None:
        return None
    import MDAnalysis as mda

    reference = mda.Universe(str(reference_file))
    selected = reference.select_atoms(settings.selection)
    if len(selected) != len(atoms) or list(selected.residues.resids) != list(atoms.residues.resids):
        raise SelectionError(
            f"rmsf: external reference {reference_file} gives "
            f"{len(selected)} atoms over {len(selected.residues)} residues for selection "
            f"{settings.selection!r}, the trajectory gives {len(atoms)} over "
            f"{len(atoms.residues)}; use a reference with the same selected atoms in the "
            "same order"
        )
    return selected.positions.astype(np.float64)


RMSFAnalysis = contract_analysis(RMSF)
