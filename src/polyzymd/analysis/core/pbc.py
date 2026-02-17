"""Periodic boundary condition (PBC) utilities for distance calculations.

This module provides PBC-aware distance functions using the minimum image
convention. All analysis modules should import distance utilities from here
to ensure consistent handling across the codebase.

Supported Box Types
-------------------
Both **orthorhombic** (cubic, rectangular) and **triclinic** (e.g., truncated
dodecahedron, rhombic dodecahedron) boxes are fully supported. This is achieved
by delegating to MDAnalysis distance functions which handle all box types.

Box Format
----------
The box parameter uses MDAnalysis format: ``[Lx, Ly, Lz, alpha, beta, gamma]``

- ``Lx, Ly, Lz``: Box edge lengths in Angstroms
- ``alpha``: Angle between edges b and c (degrees)
- ``beta``: Angle between edges a and c (degrees)
- ``gamma``: Angle between edges a and b (degrees)

For orthorhombic boxes, all angles are 90°. For triclinic boxes (like truncated
dodecahedron), angles differ from 90°.

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

import numpy as np
from numpy.typing import NDArray


def minimum_image_distance(
    pos1: NDArray[np.floating],
    pos2: NDArray[np.floating],
    box: NDArray[np.floating] | None = None,
) -> float:
    """Calculate minimum image distance between two positions.

    Uses the minimum image convention for periodic boundary conditions.
    Both orthorhombic and triclinic boxes are fully supported.

    Parameters
    ----------
    pos1, pos2 : NDArray
        Position vectors of shape (3,) in Angstroms.
    box : NDArray, optional
        Box dimensions in MDAnalysis format: [Lx, Ly, Lz, alpha, beta, gamma]
        or just [Lx, Ly, Lz] (assumes orthorhombic). If None, no PBC
        correction is applied.

    Returns
    -------
    float
        Minimum image distance in Angstroms.

    Examples
    --------
    >>> pos1 = np.array([1.0, 50.0, 50.0])
    >>> pos2 = np.array([99.0, 50.0, 50.0])
    >>> box = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])
    >>> minimum_image_distance(pos1, pos2, box)
    2.0  # Not 98.0, because pos2 is closer via the periodic image
    """
    from MDAnalysis.lib.distances import calc_bonds

    # Reshape to (1, 3) for calc_bonds which expects arrays
    p1 = np.asarray(pos1, dtype=np.float32).reshape(1, 3)
    p2 = np.asarray(pos2, dtype=np.float32).reshape(1, 3)

    if box is not None:
        box_arr = np.asarray(box, dtype=np.float32)
        # Expand 3-element box to 6-element (assume orthorhombic)
        if len(box_arr) == 3:
            box_arr = np.array(
                [box_arr[0], box_arr[1], box_arr[2], 90.0, 90.0, 90.0],
                dtype=np.float32,
            )
        result = calc_bonds(p1, p2, box=box_arr)
    else:
        result = calc_bonds(p1, p2)

    return float(result[0])


def pairwise_distances_pbc(
    positions1: NDArray[np.floating],
    positions2: NDArray[np.floating],
    box: NDArray[np.floating] | None = None,
) -> NDArray[np.float64]:
    """Calculate pairwise distances between two sets of positions.

    Uses the minimum image convention for periodic boundary conditions.
    Both orthorhombic and triclinic boxes are fully supported.

    Parameters
    ----------
    positions1 : NDArray
        First set of positions, shape (N, 3) in Angstroms.
    positions2 : NDArray
        Second set of positions, shape (M, 3) in Angstroms.
    box : NDArray, optional
        Box dimensions in MDAnalysis format: [Lx, Ly, Lz, alpha, beta, gamma]
        or just [Lx, Ly, Lz] (assumes orthorhombic). If None, no PBC
        correction is applied.

    Returns
    -------
    NDArray[np.float64]
        Distance matrix of shape (N, M) in Angstroms.

    Examples
    --------
    >>> pos1 = np.array([[1.0, 50.0, 50.0], [50.0, 50.0, 50.0]])
    >>> pos2 = np.array([[99.0, 50.0, 50.0]])
    >>> box = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])
    >>> pairwise_distances_pbc(pos1, pos2, box)
    array([[2.], [49.]])
    """
    from MDAnalysis.lib.distances import distance_array

    p1 = np.asarray(positions1, dtype=np.float32)
    p2 = np.asarray(positions2, dtype=np.float32)

    if box is not None:
        box_arr = np.asarray(box, dtype=np.float32)
        # Expand 3-element box to 6-element (assume orthorhombic)
        if len(box_arr) == 3:
            box_arr = np.array(
                [box_arr[0], box_arr[1], box_arr[2], 90.0, 90.0, 90.0],
                dtype=np.float32,
            )
        result = distance_array(p1, p2, box=box_arr)
    else:
        result = distance_array(p1, p2)

    return result.astype(np.float64)
