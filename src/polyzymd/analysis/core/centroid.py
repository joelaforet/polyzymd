"""Centroid/representative frame finding utilities.

This module provides functions to find representative frames from MD trajectories
using different methods. The representative frame is commonly used as a reference
for trajectory alignment before RMSF calculations.

Methods
-------
centroid (K-Means clustering)
    Finds the frame closest to the center of the most populated cluster.
    Uses all protein atoms by default to capture side chain conformations.
    Best for: Finding the equilibrium conformation where the protein spends
    the most time during the simulation.

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

from polyzymd.analysis.core.alignment import ReferenceMode

LOGGER = logging.getLogger(__name__)


def find_centroid_frame(
    universe: "Universe",
    selection: str = "protein",
    start_frame: int = 0,
    stop_frame: int | None = None,
    verbose: bool = True,
) -> int:
    """Find the most representative frame using K-Means clustering.

    This function identifies the frame that best represents the most populated
    conformational state in the trajectory. It uses K-Means clustering with
    a single cluster to find the centroid of all conformations, then returns
    the actual frame closest to that centroid.

    The approach finds where the protein spends most of its time during the
    simulation, making it suitable as a reference for RMSF calculations when
    you want to measure flexibility around the equilibrium state.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory to analyze.
    selection : str, optional
        MDAnalysis selection string for atoms to use in clustering.
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
        Index of the most representative frame (0-indexed, relative to full
        trajectory, not to start_frame).

    Notes
    -----
    The algorithm:
    1. Extract coordinates for selected atoms across all frames in range
    2. Reshape to 2D array: (n_frames, n_atoms * 3)
    3. Perform K-Means clustering with k=1 to find the centroid
    4. Find the frame with minimum Euclidean distance to the centroid
    5. Return the index of that frame (adjusted for start_frame offset)

    Using all protein atoms (default) rather than just CA atoms captures
    the full conformational state including side chain rotamers.

    Examples
    --------
    >>> import MDAnalysis as mda
    >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
    >>> # Find centroid after 100 frames of equilibration
    >>> centroid_idx = find_centroid_frame(u, start_frame=100)
    >>> print(f"Most representative frame: {centroid_idx}")

    >>> # Use only backbone atoms
    >>> centroid_idx = find_centroid_frame(u, selection="protein and backbone")

    See Also
    --------
    find_reference_frame : High-level function supporting multiple methods
    """
    try:
        from sklearn.cluster import KMeans
    except ImportError:
        raise ImportError(
            "scikit-learn is required for centroid frame finding.\n"
            "Install with: pip install scikit-learn"
        )

    # Select atoms for clustering
    atoms = universe.select_atoms(selection)
    if len(atoms) == 0:
        from polyzymd.analysis.core.diagnostics import get_selection_diagnostics

        diag = get_selection_diagnostics(universe, selection)
        raise ValueError(f"Selection '{selection}' matched no atoms.\n\n{diag}")

    if verbose:
        LOGGER.info(f"Finding centroid frame using {len(atoms)} atoms from '{selection}'")

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
    for i, ts in enumerate(universe.trajectory[start_frame:stop_frame]):
        coordinates[i] = atoms.positions

    # Reshape to 2D for K-Means: (n_frames, n_atoms * 3)
    n_atoms = len(atoms)
    coordinates_reshaped = coordinates.reshape(n_frames, n_atoms * 3)

    if verbose:
        LOGGER.info("Performing K-Means clustering...")

    # K-Means with 1 cluster to find centroid
    kmeans = KMeans(n_clusters=1, random_state=42, n_init=10)
    kmeans.fit(coordinates_reshaped)

    # Get cluster center and reshape back
    cluster_center = kmeans.cluster_centers_[0]
    cluster_center_reshaped = cluster_center.reshape(n_atoms, 3)

    if verbose:
        LOGGER.info("Finding frame closest to centroid...")

    # Find frame closest to centroid (Euclidean distance across all atoms)
    distances = np.linalg.norm(coordinates - cluster_center_reshaped, axis=(1, 2))
    relative_idx = int(np.argmin(distances))

    # Convert to absolute frame index
    centroid_frame_idx = relative_idx + start_frame

    if verbose:
        LOGGER.info(
            f"Most representative frame: {centroid_frame_idx} "
            f"(distance to centroid: {distances[relative_idx]:.3f} Å)"
        )

    return centroid_frame_idx


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
    mode : {"centroid", "average", "frame"}, optional
        Method for selecting the reference. Default is "centroid".

        - "centroid": K-Means clustering to find most populated state.
          Returns the frame index closest to the cluster center.
        - "average": Use average structure as reference.
          Returns None (caller should use AverageStructure).
        - "frame": Use a specific frame specified by `specific_frame`.
          Returns the specified frame index (converted to 0-indexed).

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
    >>> # Find the most populated state (equilibrium conformation)
    >>> ref_frame = find_reference_frame(u, mode="centroid", start_frame=100)

    >>> # Use average structure
    >>> ref_frame = find_reference_frame(u, mode="average")
    >>> # ref_frame is None, use AverageStructure instead

    >>> # Use a specific frame (e.g., catalytically competent state)
    >>> ref_frame = find_reference_frame(u, mode="frame", specific_frame=500)

    See Also
    --------
    find_centroid_frame : Low-level centroid finding with K-Means
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

    else:
        raise ValueError(f"Unknown mode: {mode}. Must be 'centroid', 'average', or 'frame'")


def get_reference_mode_description(mode: ReferenceMode) -> str:
    """Get a human-readable description of a reference mode.

    Parameters
    ----------
    mode : {"centroid", "average", "frame"}
        The reference mode.

    Returns
    -------
    str
        Description of what this mode represents.
    """
    descriptions = {
        "centroid": (
            "Most populated state (K-Means centroid) - "
            "measures flexibility around the equilibrium conformation"
        ),
        "average": ("Average structure - pure thermal fluctuations around the mathematical mean"),
        "frame": (
            "Specific frame - "
            "fluctuations relative to a user-defined reference (e.g., functional state)"
        ),
    }
    return descriptions.get(mode, f"Unknown mode: {mode}")
