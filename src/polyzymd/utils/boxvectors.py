"""
Box vector utilities for simulation setup.

This module provides functions for sizing the periodic box of a simulation:

- get_topology_positions: Coordinates of all atoms in a topology
- get_topology_bbox_bounds: Per-axis min/max coordinate bounds
- max_pairwise_distance: Largest atom-to-atom distance (solute diameter)
- plan_box: Box vectors for a solute, a box shape and a padding
- get_box_volume: Volume enclosed by box vectors

Adapted from the Polymerist package by Timotej Bernat, used under the MIT License.
Original source: https://github.com/timbernat/polymerist
Copyright (c) 2024 Timotej Bernat
"""

from typing import TYPE_CHECKING, Union

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from openff.toolkit import Topology
    from openff.units import Quantity


def get_topology_bbox_bounds(topology: "Topology") -> tuple[NDArray, NDArray]:
    """Return the min and max coordinate bounds of all atoms in a topology.

    This is needed when the absolute position of the bounding box matters
    (e.g. computing an exclusion zone for Packmol).

    Parameters
    ----------
    topology : openff.toolkit.Topology
        Topology with at least one molecule that has conformer coordinates.

    Returns
    -------
    min_coords : NDArray
        1-D array of shape (3,) with [xmin, ymin, zmin] in Angstrom.
    max_coords : NDArray
        1-D array of shape (3,) with [xmax, ymax, zmax] in Angstrom.

    Raises
    ------
    ValueError
        If the topology has no molecules with coordinates.
    """
    all_coords = get_topology_positions(topology)
    min_coords: NDArray = np.min(all_coords, axis=0)
    max_coords: NDArray = np.max(all_coords, axis=0)

    return min_coords, max_coords


def get_topology_positions(topology: "Topology") -> NDArray:
    """Return the coordinates of every atom with a conformer, in Angstrom.

    Raises
    ------
    ValueError
        If the topology has no molecules with coordinates.
    """
    all_positions = [
        molecule.conformers[0].m_as("angstrom")
        for molecule in topology.molecules
        if molecule.n_conformers > 0
    ]
    if not all_positions:
        raise ValueError("Cannot size the box: no molecules in topology have coordinates")
    return np.vstack(all_positions)


def max_pairwise_distance(positions: NDArray) -> float:
    """Return the largest distance between two points (the solute diameter).

    Only the points on the convex hull can be the farthest pair, so the
    distances are computed between those points.
    """
    from scipy.spatial import ConvexHull, QhullError
    from scipy.spatial.distance import pdist

    points = np.asarray(positions, dtype=float)
    if len(points) < 2:
        return 0.0
    try:
        points = points[ConvexHull(points).vertices]
    except (QhullError, ValueError):
        pass  # flat or tiny point sets: use all points
    return float(pdist(points).max())


def plan_box(
    positions_nm: NDArray,
    shape_matrix: NDArray,
    padding_nm: float,
    margin_nm: float,
) -> dict:
    """Size a periodic cell for a solute.

    The cell edge is the solute diameter plus ``2 * padding_nm``. Every
    lattice vector of a cube or a rhombic dodecahedron is
    at least one edge long, so the solute is at least ``2 * padding_nm`` from
    each of its periodic copies, in any orientation.

    The solute is packed in the rectangular brick of the cell (the diagonal of
    the box vectors). The edge also grows, if needed, so that the solute's
    bounding box fits in that brick with ``margin_nm`` to every face.

    Parameters
    ----------
    positions_nm : NDArray
        Solute coordinates, shape (N, 3), in nm.
    shape_matrix : NDArray
        3x3 box-shape matrix with unit-length rows (e.g. ``openff.packmol.UNIT_CUBE``).
    padding_nm : float
        Distance from the solute to the cell edge, in nm.
    margin_nm : float
        Smallest distance from the solute bounding box to a brick face, in nm.

    Returns
    -------
    dict
        ``box_vectors`` (3x3, nm), ``padding``, ``edge``, ``diameter``, ``extent`` (bounding
        box, 3 values), ``brick`` (3 values) and ``clearance`` (solute bounding
        box to each brick face, 3 values), all in nm.
    """
    positions = np.asarray(positions_nm, dtype=float)
    shape = np.asarray(shape_matrix, dtype=float)
    extent = positions.max(axis=0) - positions.min(axis=0)
    diameter = max_pairwise_distance(positions)
    edge = max(
        diameter + 2.0 * padding_nm, float(np.max((extent + 2.0 * margin_nm) / np.diagonal(shape)))
    )
    box_vectors = edge * shape
    brick = np.diagonal(box_vectors).copy()
    return {
        "box_vectors": box_vectors,
        "padding": float(padding_nm),
        "edge": edge,
        "diameter": diameter,
        "extent": extent,
        "brick": brick,
        "clearance": (brick - extent) / 2.0,
    }


def describe_box_plan(plan: dict) -> str:
    """Return a one-line description of a :func:`plan_box` result."""
    return (
        "edge %.2f nm (solute diameter %.2f nm + 2 x %.2f nm padding); "
        "brick %.2f x %.2f x %.2f nm; solute bounding box %.2f x %.2f x %.2f nm; "
        "clearance to the brick faces %.2f / %.2f / %.2f nm"
        % (
            plan["edge"],
            plan["diameter"],
            plan["padding"],
            *plan["brick"],
            *plan["extent"],
            *plan["clearance"],
        )
    )


def get_box_volume(
    box_vectors: "Quantity",
    units_as_openmm: bool = False,
) -> "Quantity":
    """Calculate the volume enclosed by box vectors.

    This function computes the volume of the parallelopiped defined by
    the three box vectors using the scalar triple product formula:
    V = |a . (b x c)|

    Args:
        box_vectors: A 3x3 matrix of box vectors (with length units).
            Each row represents one box vector [a, b, c].
        units_as_openmm: If True, return volume with OpenMM units.
            If False (default), return with OpenFF units.
            Note: This parameter is named for API compatibility with
            Polymerist but the actual return type is always OpenFF Quantity.

    Returns:
        The box volume with appropriate cubic length units.

    Example:
        >>> volume = get_box_volume(box_vectors)
        >>> print(volume.to("nanometer**3"))
    """
    from openff.units import Quantity

    # Get the magnitude in consistent units
    box_angstrom = box_vectors.m_as("angstrom")

    # Calculate volume using the determinant (equivalent to scalar triple product)
    # For a matrix where rows are the box vectors, |det(M)| gives the volume
    volume = abs(np.linalg.det(box_angstrom))

    # Return as Quantity with cubic angstrom units
    return Quantity(volume, "angstrom**3")
