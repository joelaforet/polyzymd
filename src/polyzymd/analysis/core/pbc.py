"""Periodic boundary condition (PBC) utilities for distance calculations.

This module provides PBC-aware distance functions using the minimum image
convention. All analysis modules should import distance utilities from here
to ensure consistent handling across the codebase.

Supported Box Types
-------------------
- **Orthorhombic boxes** (cubic, rectangular): Fully supported
- **Triclinic boxes**: Not supported; a warning is logged and PBC correction
  is skipped (falls back to simple Euclidean distance)

Usage
-----
>>> from polyzymd.analysis.core.pbc import minimum_image_distance
>>>
>>> # Single distance calculation
>>> pos1 = np.array([1.0, 2.0, 3.0])
>>> pos2 = np.array([99.0, 2.0, 3.0])  # Near box boundary
>>> box = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])
>>>
>>> # Without PBC: distance would be 98.0 Å
>>> # With PBC: distance is 2.0 Å (minimum image)
>>> dist = minimum_image_distance(pos1, pos2, box)

Notes
-----
The minimum image convention assumes that the distance of interest is less
than half the box length. For very small boxes or extended molecules, this
assumption may break down.

References
----------
- Allen & Tildesley, "Computer Simulation of Liquids", Chapter 1.5
- MDAnalysis: https://docs.mdanalysis.org/stable/documentation_pages/lib/distances.html
"""

from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    pass

LOGGER = logging.getLogger(__name__)

# Track whether we've warned about triclinic boxes (warn once)
_TRICLINIC_WARNING_ISSUED = False


def is_orthorhombic(box: NDArray[np.floating]) -> bool:
    """Check if box is orthorhombic (all angles approximately 90°).

    Parameters
    ----------
    box : NDArray
        Box dimensions in MDAnalysis format: [Lx, Ly, Lz, alpha, beta, gamma]
        or just lengths [Lx, Ly, Lz].

    Returns
    -------
    bool
        True if box is orthorhombic (angles within 0.01° of 90°).
    """
    if box is None:
        return False

    # If only lengths provided, assume orthorhombic
    if len(box) == 3:
        return True

    if len(box) >= 6:
        alpha, beta, gamma = box[3:6]
        # Check if all angles are approximately 90 degrees
        return all(abs(angle - 90.0) < 0.01 for angle in [alpha, beta, gamma])

    return False


def minimum_image_distance(
    pos1: NDArray[np.floating],
    pos2: NDArray[np.floating],
    box: NDArray[np.floating] | None = None,
) -> float:
    """Calculate minimum image distance between two positions.

    Uses the minimum image convention for periodic boundary conditions.
    For orthorhombic boxes, this finds the shortest distance considering
    all periodic images.

    Parameters
    ----------
    pos1, pos2 : NDArray
        Position vectors of shape (3,) in Angstroms.
    box : NDArray, optional
        Box dimensions in MDAnalysis format: [Lx, Ly, Lz, alpha, beta, gamma]
        or just [Lx, Ly, Lz]. If None, no PBC correction is applied.

    Returns
    -------
    float
        Minimum image distance in Angstroms.

    Warnings
    --------
    If a triclinic box is detected (non-90° angles), a warning is logged
    once per session and PBC correction is skipped.

    Examples
    --------
    >>> pos1 = np.array([1.0, 50.0, 50.0])
    >>> pos2 = np.array([99.0, 50.0, 50.0])
    >>> box = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])
    >>> minimum_image_distance(pos1, pos2, box)
    2.0  # Not 98.0, because pos2 is closer via the periodic image
    """
    global _TRICLINIC_WARNING_ISSUED

    diff = pos2 - pos1

    if box is not None:
        if is_orthorhombic(box):
            # Apply minimum image convention for orthorhombic box
            box_lengths = box[:3]
            diff = diff - box_lengths * np.round(diff / box_lengths)
        else:
            # Triclinic box - warn and skip PBC
            if not _TRICLINIC_WARNING_ISSUED:
                warnings.warn(
                    "Triclinic box detected. PBC correction is only implemented for "
                    "orthorhombic boxes. Falling back to simple Euclidean distance. "
                    "This warning is shown once per session.",
                    UserWarning,
                    stacklevel=2,
                )
                LOGGER.warning(
                    "Triclinic box detected (angles: %.1f, %.1f, %.1f). PBC correction skipped.",
                    box[3] if len(box) > 3 else 90,
                    box[4] if len(box) > 4 else 90,
                    box[5] if len(box) > 5 else 90,
                )
                _TRICLINIC_WARNING_ISSUED = True

    return float(np.linalg.norm(diff))


def pairwise_distances_pbc(
    positions1: NDArray[np.floating],
    positions2: NDArray[np.floating],
    box: NDArray[np.floating] | None = None,
) -> NDArray[np.float64]:
    """Calculate pairwise distances between two sets of positions.

    Uses minimum image convention if an orthorhombic box is provided.

    Parameters
    ----------
    positions1 : NDArray
        First set of positions, shape (N, 3) in Angstroms.
    positions2 : NDArray
        Second set of positions, shape (M, 3) in Angstroms.
    box : NDArray, optional
        Box dimensions in MDAnalysis format: [Lx, Ly, Lz, alpha, beta, gamma]
        or just [Lx, Ly, Lz]. If None, no PBC correction is applied.

    Returns
    -------
    NDArray[np.float64]
        Distance matrix of shape (N, M) in Angstroms.

    Notes
    -----
    For large arrays, this can be memory-intensive as it creates an
    intermediate array of shape (N, M, 3). For very large systems,
    consider using MDAnalysis.lib.distances.distance_array() which
    is optimized for performance.

    Examples
    --------
    >>> pos1 = np.array([[1.0, 50.0, 50.0], [50.0, 50.0, 50.0]])
    >>> pos2 = np.array([[99.0, 50.0, 50.0]])
    >>> box = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])
    >>> pairwise_distances_pbc(pos1, pos2, box)
    array([[2.], [49.]])
    """
    global _TRICLINIC_WARNING_ISSUED

    # Compute difference vectors: (N, 1, 3) - (1, M, 3) -> (N, M, 3)
    diff = positions1[:, np.newaxis, :] - positions2[np.newaxis, :, :]

    if box is not None:
        if is_orthorhombic(box):
            box_lengths = box[:3]
            diff = diff - box_lengths * np.round(diff / box_lengths)
        else:
            # Triclinic box - warn and skip PBC
            if not _TRICLINIC_WARNING_ISSUED:
                warnings.warn(
                    "Triclinic box detected. PBC correction is only implemented for "
                    "orthorhombic boxes. Falling back to simple Euclidean distance. "
                    "This warning is shown once per session.",
                    UserWarning,
                    stacklevel=2,
                )
                LOGGER.warning(
                    "Triclinic box detected (angles: %.1f, %.1f, %.1f). PBC correction skipped.",
                    box[3] if len(box) > 3 else 90,
                    box[4] if len(box) > 4 else 90,
                    box[5] if len(box) > 5 else 90,
                )
                _TRICLINIC_WARNING_ISSUED = True

    # Compute distances
    distances = np.linalg.norm(diff, axis=2)
    return distances.astype(np.float64)


def reset_triclinic_warning() -> None:
    """Reset the triclinic warning flag.

    This is primarily useful for testing. In production, the warning
    should only be shown once per session.
    """
    global _TRICLINIC_WARNING_ISSUED
    _TRICLINIC_WARNING_ISSUED = False
