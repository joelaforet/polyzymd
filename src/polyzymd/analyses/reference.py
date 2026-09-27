"""Reference structures for per-frame functions, built once per replicate.

:func:`reference` stands for an ``AtomGroup`` that holds the reference
coordinates of a selection. :func:`build_reference` builds it for one
replicate and copies it into a separate one-frame universe with
``MDAnalysis.Merge``, so nothing that runs later can move it and the
replicate's trajectory is never modified.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from polyzymd.analyses.exceptions import ProtocolError

if TYPE_CHECKING:
    import numpy as np

MODES = ("external", "frame", "average", "centroid")


@dataclass(frozen=True)
class Reference:
    """Reference coordinates of ``selection``, built per replicate."""

    mode: str
    selection: str
    alignment: str
    frame: int | None = None
    file: str | None = None


def reference(
    mode: str,
    selection: str,
    *,
    frame: int | None = None,
    file: str | Path | None = None,
    alignment: str | None = None,
) -> Reference:
    """Stand for reference coordinates of ``selection`` in each replicate.

    Parameters
    ----------
    mode : {"external", "frame", "average", "centroid"}
        ``"external"`` reads the atoms of ``selection`` from ``file``.
        ``"frame"`` takes them from production frame ``frame``, counted
        from 1 at the first frame after the equilibration window.
        ``"average"`` superposes the ``alignment`` atoms of every production
        frame on the first production frame, averages, superposes every
        frame again on that average and takes the mean positions.
        ``"centroid"`` takes the production frame whose ``alignment`` atoms
        have the smallest RMSD to their iterative average structure, from
        :func:`~polyzymd.analyses.shared.centroid.find_centroid_frame`.
    selection : str
        MDAnalysis selection string of the atoms the reference holds.
    frame : int, optional
        Production frame, from 1, for ``mode="frame"``.
    file : str or Path, optional
        Structure file for ``mode="external"``. Its SHA-256 hash is recorded.
    alignment : str, optional
        Selection that is superposed in ``"average"`` and ``"centroid"``
        modes. Defaults to ``selection``.

    Returns
    -------
    Reference
        Placeholder that :func:`~polyzymd.analyses.timeseries.run_timeseries`
        replaces with an ``AtomGroup`` of a one-frame universe.

    Raises
    ------
    ProtocolError
        If the mode is unknown, ``frame`` is missing or below 1 in frame
        mode, or ``file`` is missing in external mode.
    """
    if mode not in MODES:
        raise ProtocolError(f"Unknown reference mode {mode!r}.", hint=f"Use one of {MODES}.")
    if mode == "frame" and (isinstance(frame, bool) or not isinstance(frame, int) or frame < 1):
        raise ProtocolError(
            f"Reference mode 'frame' needs a production frame from 1, got {frame!r}.",
            hint="Pass frame=1 for the first frame after the equilibration window.",
        )
    path = Path(file).expanduser().resolve() if file is not None else None
    if mode == "external" and (path is None or not path.is_file()):
        raise ProtocolError(
            f"Reference mode 'external' needs an existing structure file, got {file!r}.",
            hint="Pass file='reference.pdb'.",
        )
    return Reference(
        mode,
        str(selection),
        str(alignment or selection),
        frame if mode == "frame" else None,
        str(path) if mode == "external" else None,
    )


def build_reference(
    ref: Reference, universe: Any, frames: np.ndarray
) -> tuple[Any, dict[str, int]]:
    """Build the reference atoms of one replicate.

    Parameters
    ----------
    ref : Reference
        What to build.
    universe : MDAnalysis.Universe
        The replicate's universe. Its coordinates are read, never changed.
    frames : numpy.ndarray
        The replicate's production frame indices, consecutive.

    Returns
    -------
    atoms : MDAnalysis.core.groups.AtomGroup
        All atoms of a one-frame universe holding the reference coordinates.
    chosen : dict
        For ``"frame"`` and ``"centroid"``, the production frame from 1 as
        ``frame`` and its trajectory index as ``trajectory_frame``.

    Raises
    ------
    ProtocolError
        If a selection matches no atoms, the external file has a different
        number of atoms for ``selection``, or ``frame`` is past the last
        production frame.
    """
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    from polyzymd.analyses.timeseries import Select, _build

    atoms = _build(Select(ref.selection), universe)
    chosen: dict[str, int] = {}
    if ref.mode == "external":
        positions = mda.Universe(ref.file).select_atoms(ref.selection).positions
        if len(positions) != len(atoms):
            raise ProtocolError(
                f"{ref.file} has {len(positions)} atoms for {ref.selection!r} and the "
                f"trajectory has {len(atoms)}.",
                hint="Use a reference file with the same atoms as the simulation topology.",
            )
    elif ref.mode == "average":
        fit = _build(Select(ref.alignment), universe)
        group = atoms | fit
        coordinates = np.array([group.positions for _ in universe.trajectory[frames]], float)
        where = np.searchsorted(group.indices, fit.indices)
        first = _fitted_mean(coordinates, where, coordinates[0][where])
        positions = _fitted_mean(coordinates, where, first[where])
        positions = positions[np.searchsorted(group.indices, atoms.indices)]
    else:
        if ref.mode == "centroid":
            from polyzymd.analyses.shared.centroid import find_centroid_frame

            index = find_centroid_frame(
                universe, ref.alignment, int(frames[0]), int(frames[-1]) + 1, verbose=False
            )
        elif ref.frame > len(frames):
            raise ProtocolError(
                f"Reference frame {ref.frame} is past the last of {len(frames)} production frames.",
                hint="Pass a smaller frame, counted from 1 after the equilibration window.",
            )
        else:
            index = int(frames[ref.frame - 1])
        universe.trajectory[index]
        positions = atoms.positions
        chosen = {"frame": index - int(frames[0]) + 1, "trajectory_frame": index}
    copy = mda.Merge(atoms)
    copy.load_new(np.asarray(positions, np.float32)[np.newaxis], format=MemoryReader)
    return copy.atoms, chosen


def _fitted_mean(coordinates: np.ndarray, fit: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Average the frames after superposing each frame's ``fit`` atoms on ``target``.

    Each frame is centred on the centre of geometry of its ``fit`` atoms,
    rotated by ``MDAnalysis.analysis.align.rotation_matrix`` onto ``target``
    and moved to the centre of ``target``, as ``AlignTraj`` does.
    """
    import numpy as np
    from MDAnalysis.analysis.align import rotation_matrix

    center = target.mean(axis=0)
    total = np.zeros(coordinates.shape[1:])
    for frame in coordinates:
        origin = frame[fit].mean(axis=0)
        rotation, _ = rotation_matrix(frame[fit] - origin, target - center)
        total += (frame - origin) @ rotation.T + center
    return total / len(coordinates)
