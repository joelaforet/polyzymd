"""Representative frame finding utilities.

This module provides functions to find representative frames from MD trajectories
using different methods. The representative frame is commonly used as a reference
for trajectory alignment before RMSF calculations.

Methods
-------
centroid (aligned-mean representative frame)
    Finds the frame closest to the aligned mean structure.
    Uses all protein atoms by default to capture side chain conformations.
    Best for: Finding a representative equilibrium conformation while
    removing rigid-body translation and rotation effects.

average
    Aligns to an average structure computed from all frames.
    Note: The average structure is synthetic and may have unphysical geometry
    (e.g., distorted bond lengths/angles).
    Best for: Pure mathematical measure of thermal fluctuations around the mean.

frame
    Uses a specific frame as the reference (user-specified).
    Best for: Analyzing fluctuations relative to a known functional state,
    such as a catalytically competent conformation.

See Also
--------
The documentation at docs/analysis/reference_selection.md provides detailed
guidance on when to use each method.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    import MDAnalysis as mda
    from MDAnalysis.core.universe import Universe

from polyzymd.analyses.shared.alignment import ReferenceMode

LOGGER = logging.getLogger(__name__)


def _kabsch_align_to_reference(
    coordinates: NDArray[np.float64],
    reference: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Align coordinates to a reference with Kabsch superposition.

    Parameters
    ----------
    coordinates : NDArray[np.float64]
        Coordinates to align with shape (n_atoms, 3).
    reference : NDArray[np.float64]
        Reference coordinates with shape (n_atoms, 3).

    Returns
    -------
    NDArray[np.float64]
        Coordinates aligned to the centered reference.
    """
    centered = coordinates - np.mean(coordinates, axis=0)
    ref_centered = reference - np.mean(reference, axis=0)

    covariance = centered.T @ ref_centered
    u_matrix, _, v_t = np.linalg.svd(covariance)

    rotation = u_matrix @ v_t
    if np.linalg.det(rotation) < 0:
        u_matrix[:, -1] *= -1.0
        rotation = u_matrix @ v_t

    return centered @ rotation


def _find_frame_closest_to_aligned_mean(coordinates: NDArray[np.float64]) -> tuple[int, float]:
    """Find frame closest to aligned mean coordinates.

    Parameters
    ----------
    coordinates : NDArray[np.float64]
        Raw coordinates with shape (n_frames, n_atoms, 3).

    Returns
    -------
    tuple[int, float]
        Relative index of frame closest to the aligned mean structure and the
        corresponding RMSD to aligned mean.
    """
    if coordinates.ndim != 3:
        raise ValueError("coordinates must have shape (n_frames, n_atoms, 3)")

    n_frames = coordinates.shape[0]
    if n_frames == 1:
        return 0, 0.0

    reference = coordinates[0]
    aligned = np.empty_like(coordinates)
    for frame_idx in range(n_frames):
        aligned[frame_idx] = _kabsch_align_to_reference(coordinates[frame_idx], reference)

    mean_coordinates = np.mean(aligned, axis=0)
    rmsd_to_mean = np.sqrt(np.mean((aligned - mean_coordinates) ** 2, axis=(1, 2)))
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
    It first performs rigid-body alignment of each frame to a common reference,
    computes the mean coordinates in aligned space, then returns the trajectory
    frame with minimum RMSD to that aligned mean.

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
    The algorithm aligns all candidate frames to a common reference frame,
    computes the aligned mean structure, then selects the frame with minimum
    RMSD to that mean.

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

    See Also
    --------
    find_reference_frame : High-level function supporting multiple methods
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


def find_reference_frame(
    universe: "Universe",
    mode: ReferenceMode = "centroid",
    selection: str = "protein",
    start_frame: int = 0,
    stop_frame: int | None = None,
    specific_frame: int | None = None,
    verbose: bool = True,
) -> int | None:
    """Find a reference frame for trajectory alignment.

    This is the high-level interface for selecting a reference structure
    for RMSF calculations. It supports multiple methods for choosing the
    reference, each appropriate for different scientific questions.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory.
    mode : {"centroid", "average", "frame", "external"}, optional
        Method for selecting the reference. Default is "centroid".

        - "centroid": Representative aligned frame mode.
          Returns the frame index closest to the aligned mean structure.
        - "average": Use average structure as reference.
          Returns None (caller should use AverageStructure).
        - "frame": Use a specific frame specified by `specific_frame`.
          Returns the specified frame index (converted to 0-indexed).
        - "external": Use an external PDB file as reference.
          Returns None (alignment handled by align_trajectory).

    selection : str, optional
        MDAnalysis selection string for atoms to use. Default is "protein".
        Only used for "centroid" mode.
    start_frame : int, optional
        First frame for analysis (0-indexed). Default is 0.
    stop_frame : int, optional
        Last frame for analysis (exclusive). Default is None.
    specific_frame : int, optional
        Frame index to use when mode="frame" (1-indexed, PyMOL convention).
        Required when mode="frame".
    verbose : bool, optional
        Log progress messages. Default is True.

    Returns
    -------
    int or None
        Frame index (0-indexed) to use as reference, or None if mode="average"
        (indicating the caller should compute an average structure).

    Raises
    ------
    ValueError
        If mode="frame" but specific_frame is not provided.
        If specific_frame is out of range.

    Examples
    --------
    >>> # Find a representative aligned frame (equilibrium conformation)
    >>> ref_frame = find_reference_frame(u, mode="centroid", start_frame=100)

    >>> # Use average structure
    >>> ref_frame = find_reference_frame(u, mode="average")
    >>> # ref_frame is None, use AverageStructure instead

    >>> # Use a specific frame (e.g., catalytically competent state)
    >>> ref_frame = find_reference_frame(u, mode="frame", specific_frame=500)

    See Also
    --------
    find_centroid_frame : Low-level representative frame selection
    """
    if mode == "centroid":
        return find_centroid_frame(
            universe,
            selection=selection,
            start_frame=start_frame,
            stop_frame=stop_frame,
            verbose=verbose,
        )

    elif mode == "average":
        if verbose:
            LOGGER.info("Using average structure as reference (will be computed during alignment)")
        return None  # Caller should use AverageStructure

    elif mode == "frame":
        if specific_frame is None:
            raise ValueError("specific_frame is required when mode='frame'")

        # Convert from 1-indexed (PyMOL) to 0-indexed
        frame_idx = specific_frame - 1

        # Validate
        n_frames = len(universe.trajectory)
        if frame_idx < 0 or frame_idx >= n_frames:
            raise ValueError(
                f"specific_frame={specific_frame} (1-indexed) is out of range. "
                f"Valid range: 1 to {n_frames}"
            )

        if verbose:
            LOGGER.info(f"Using user-specified frame {specific_frame} (0-indexed: {frame_idx})")

        return frame_idx

    elif mode == "external":
        # External PDB reference — no trajectory frame to return.
        # Alignment to external PDB is handled by align_trajectory() in alignment.py.
        if verbose:
            LOGGER.info("Using external PDB as reference (alignment handled by align_trajectory)")
        return None

    else:
        raise ValueError(
            f"Unknown mode: {mode}. Must be 'centroid', 'average', 'frame', or 'external'"
        )


def get_reference_mode_description(mode: ReferenceMode) -> str:
    """Get a human-readable description of a reference mode.

    Parameters
    ----------
    mode : {"centroid", "average", "frame", "external"}
        The reference mode.

    Returns
    -------
    str
        Description of what this mode represents.
    """
    descriptions = {
        "centroid": (
            "Representative aligned frame (closest to aligned mean) - "
            "measures flexibility around a representative equilibrium conformation"
        ),
        "average": ("Average structure - pure thermal fluctuations around the mathematical mean"),
        "frame": (
            "Specific frame - "
            "fluctuations relative to a user-defined reference (e.g., functional state)"
        ),
        "external": (
            "External PDB structure - "
            "deviations from a condition-independent reference geometry "
            "(e.g., catalytically competent crystal structure)"
        ),
    }
    return descriptions.get(mode, f"Unknown mode: {mode}")
