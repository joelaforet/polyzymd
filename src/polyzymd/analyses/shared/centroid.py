"""Find the representative frame of a trajectory for a centroid reference.

:func:`find_centroid_frame` returns the frame closest, after optimal
superposition, to the iterative average structure of
``MDAnalysis.analysis.align.iterative_average``.
:func:`polyzymd.analyses.reference.reference` calls it for the
``"centroid"`` reference mode.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    import MDAnalysis as mda
    from MDAnalysis.core.universe import Universe

LOGGER = logging.getLogger(__name__)


def _find_frame_closest_to_aligned_mean(coordinates: NDArray[np.float64]) -> tuple[int, float]:
    """Find the frame closest to the iterative average structure.

    The average comes from ``MDAnalysis.analysis.align.iterative_average``,
    which superposes every frame on the current average and averages again
    until the average moves by less than 1e-4 Å. The MDAnalysis default of
    1e-6 Å is below the float32 rounding of coordinates near 50 Å, where it
    does not converge. Each frame's distance to it
    is ``MDAnalysis.analysis.rms.rmsd`` after optimal superposition, the
    square root of the mean squared distance per atom.

    Parameters
    ----------
    coordinates : NDArray[np.float64]
        Raw coordinates with shape (n_frames, n_atoms, 3).

    Returns
    -------
    tuple[int, float]
        Relative index of the frame closest to the average structure and its
        RMSD to that average in Å.
    """
    if coordinates.ndim != 3:
        raise ValueError("coordinates must have shape (n_frames, n_atoms, 3)")

    if coordinates.shape[0] == 1:
        return 0, 0.0

    import contextlib
    import io
    import warnings

    import MDAnalysis as mda
    from MDAnalysis.analysis import align, rms
    from MDAnalysis.coordinates.memory import MemoryReader

    frames = mda.Universe.empty(coordinates.shape[1], trajectory=True)
    frames.add_TopologyAttr("masses", np.ones(coordinates.shape[1]))
    frames.load_new(coordinates.astype(np.float32), format=MemoryReader)
    # iterative_average always draws a progress bar on stderr and warns that
    # this universe has no atom types, so both are silenced here.
    with warnings.catch_warnings(), contextlib.redirect_stderr(io.StringIO()):
        warnings.simplefilter("ignore")
        mean_coordinates = align.iterative_average(frames, eps=1e-4).results.positions
    rmsd_to_mean = np.array(
        [
            rms.rmsd(frame, mean_coordinates, center=True, superposition=True)
            for frame in coordinates
        ]
    )
    relative_idx = int(np.argmin(rmsd_to_mean))
    return relative_idx, float(rmsd_to_mean[relative_idx])


def find_centroid_frame(
    universe: "Universe",
    selection: str = "protein",
    start_frame: int = 0,
    stop_frame: int | None = None,
    verbose: bool = True,
) -> int:
    """Find a representative aligned frame.

    This function identifies the frame closest to the aligned mean structure.
    The mean is the iterative average structure of
    ``MDAnalysis.analysis.align.iterative_average``, and the frame returned is
    the one with the smallest RMSD to it after optimal superposition.

    This approach avoids contamination from translation/rotation and provides a
    scientifically defensible representative frame for downstream alignment and
    RMSF calculations.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory to analyze.
    selection : str, optional
        MDAnalysis selection string for atoms used to find the representative
        frame.
        Default is "protein" (all protein atoms) to capture both backbone
        and side chain conformations.
    start_frame : int, optional
        First frame to include in analysis (0-indexed). Default is 0.
        Use this to skip equilibration frames.
    stop_frame : int, optional
        Last frame to include (exclusive). Default is None (all frames).
    verbose : bool, optional
        If True, log progress messages. Default is True.

    Returns
    -------
    int
        Index of the representative frame (0-indexed, relative to full
        trajectory, not to ``start_frame``).

    Notes
    -----
    Superposing every frame on the iterative average and averaging again
    gives the same average back. A single pass that superposes every frame on
    one frame gives an average that depends on which frame that was.

    Using all protein atoms (default) rather than just CA atoms captures
    the full conformational state including side chain rotamers.

    Examples
    --------
    >>> import MDAnalysis as mda
    >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
    >>> # Find centroid after 100 frames of equilibration
    >>> centroid_idx = find_centroid_frame(u, start_frame=100)
    >>> print(f"Representative aligned frame: {centroid_idx}")

    >>> # Use only backbone atoms
    >>> centroid_idx = find_centroid_frame(u, selection="protein and backbone")
    """
    # Select atoms for representative-frame calculation
    atoms = universe.select_atoms(selection)
    if len(atoms) == 0:
        from polyzymd.analyses.shared.diagnostics import get_selection_diagnostics

        diag = get_selection_diagnostics(universe, selection)
        raise ValueError(f"Selection '{selection}' matched no atoms.\n\n{diag}")

    if verbose:
        LOGGER.info(
            f"Finding representative aligned frame using {len(atoms)} atoms from '{selection}'"
        )

    # Determine frame range
    n_frames_total = len(universe.trajectory)
    if stop_frame is None:
        stop_frame = n_frames_total

    # Validate frame range
    if start_frame < 0 or start_frame >= n_frames_total:
        raise ValueError(f"start_frame={start_frame} is out of range [0, {n_frames_total})")
    if stop_frame <= start_frame:
        raise ValueError(f"stop_frame={stop_frame} must be greater than start_frame={start_frame}")

    n_frames = stop_frame - start_frame

    if verbose:
        LOGGER.info(f"Analyzing frames {start_frame} to {stop_frame - 1} ({n_frames} frames)")

    # Collect coordinates for all frames in range
    if verbose:
        LOGGER.info("Collecting coordinates...")

    coordinates = np.empty((n_frames, len(atoms), 3), dtype=np.float64)
    for i, _ in enumerate(universe.trajectory[start_frame:stop_frame]):
        coordinates[i] = atoms.positions

    if verbose:
        LOGGER.info("Aligning frames and selecting representative frame...")

    relative_idx, rmsd_to_mean = _find_frame_closest_to_aligned_mean(coordinates)

    # Convert to absolute frame index
    representative_frame_idx = relative_idx + start_frame

    if verbose:
        LOGGER.info(
            f"Representative aligned frame: {representative_frame_idx} "
            f"(RMSD to aligned mean: {rmsd_to_mean:.3f} Å)"
        )

    return representative_frame_idx
