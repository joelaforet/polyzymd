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

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Union

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from openff.toolkit import Topology
    from openff.units import Quantity

LOGGER = logging.getLogger(__name__)


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


def plan_rectangular_box(
    positions_nm: NDArray,
    padding_xyz_nm: "tuple[float, float, float]",
    margin_nm: float,
) -> dict:
    """Size a rectangular (orthorhombic) cell for a solute in its current orientation.

    Each edge is the solute's bounding-box extent along that axis plus twice
    that axis's padding, so the solute is ``2 * padding`` from its periodic
    copy along each axis *as long as it keeps this orientation*. Unlike a cube
    or a rhombic dodecahedron, the cell is not sized for any orientation: a
    solute that rotates during the run can come closer to its images.

    Parameters
    ----------
    positions_nm : NDArray
        Solute coordinates, shape (N, 3), in nm.
    padding_xyz_nm : tuple of three floats
        Padding along x, y and z, in nm.
    margin_nm : float
        Smallest distance from the solute bounding box to a face, in nm.

    Returns
    -------
    dict
        The keys of :func:`plan_box`, plus ``shape`` = ``"rectangular"`` and
        ``padding_xyz`` (3 values).
    """
    positions = np.asarray(positions_nm, dtype=float)
    padding = np.asarray(padding_xyz_nm, dtype=float)
    extent = positions.max(axis=0) - positions.min(axis=0)
    brick = np.maximum(extent + 2.0 * padding, extent + 2.0 * margin_nm)
    return {
        "shape": "rectangular",
        "box_vectors": np.diag(brick),
        "padding": float(padding.min()),
        "padding_xyz": [float(v) for v in padding],
        "edge": float(brick.max()),
        "diameter": max_pairwise_distance(positions),
        "extent": extent,
        "brick": brick.copy(),
        "clearance": (brick - extent) / 2.0,
    }


#: Adjacent principal moments closer than this fraction are treated as equal.
DEGENERATE_MOMENT_FRACTION = 0.01


def _canonical_sign(vector: NDArray) -> NDArray:
    """``vector`` with its component of largest magnitude made positive."""
    return vector if vector[int(np.argmax(np.abs(vector)))] >= 0 else -vector


def principal_axes_rotation(positions_nm: NDArray, masses: NDArray | None = None) -> NDArray:
    """Rotation matrix that puts a solute's principal axes on x, y and z.

    The longest axis (smallest moment of inertia) goes on z and the shortest
    (largest moment) on x. The matrix is a proper rotation (determinant +1);
    apply it as ``(positions - center) @ R.T + center``.

    The result does not depend on atom order or on noise far below the
    solute's size: the sign of each axis is fixed so that its component of
    largest magnitude is positive, and the shortest axis (x) is flipped when
    needed to make the rotation proper. When two adjacent moments are within
    1% (:data:`DEGENERATE_MOMENT_FRACTION`), the axes in their plane are
    arbitrary, so they are not rotated: the input axis that belongs on the
    later of the two rows is projected onto the plane, which keeps the input
    orientation within it. A solute with three equal moments is not rotated.

    Parameters
    ----------
    positions_nm : NDArray
        Coordinates, shape (N, 3); PolyzyMD passes the heavy atoms of the fixed
        solutes (see :func:`fixed_solute_alignment`).
    masses : NDArray or None
        Per-atom masses; equal masses when None.

    Returns
    -------
    NDArray
        3x3 rotation matrix.
    """
    positions = np.asarray(positions_nm, dtype=float)
    weights = np.ones(len(positions)) if masses is None else np.asarray(masses, dtype=float)
    center = (positions * weights[:, None]).sum(axis=0) / weights.sum()
    r = positions - center
    inertia = np.eye(3) * (weights * (r * r).sum(axis=1)).sum() - (r * weights[:, None]).T @ r
    moments, axes = np.linalg.eigh(inertia)  # ascending moments
    # Rows of R are the new x, y, z directions: largest moment -> x, smallest -> z.
    rows = [2, 1, 0]
    rotation = np.stack([_canonical_sign(axes[:, i]) for i in rows])
    order = moments[rows]  # moments of the x, y and z rows
    close = [
        abs(order[i] - order[i + 1])
        <= DEGENERATE_MOMENT_FRACTION * max(order[i], order[i + 1], np.finfo(float).tiny)
        for i in (0, 1)
    ]
    if all(close):
        LOGGER.info("The solute's three principal moments are within 1%%: it is not rotated")
        return np.eye(3)
    if any(close):
        first = 0 if close[0] else 1  # the degenerate rows are (first, first + 1)
        second, fixed = first + 1, 2 if first == 0 else 0
        normal = rotation[fixed]
        LOGGER.info(
            "Principal moments %.4g and %.4g are within 1%%: the %s-%s plane keeps the "
            "input orientation",
            order[first],
            order[second],
            "xyz"[first],
            "xyz"[second],
        )
        identity = np.eye(3)
        projected = identity[second] - identity[second] @ normal * normal
        if np.linalg.norm(projected) > 0.1:
            rotation[second] = projected / np.linalg.norm(projected)
            rotation[first] = np.cross(rotation[second], normal)
            target = first
        else:  # that input axis is almost along the normal: keep the other one
            projected = identity[first] - identity[first] @ normal * normal
            rotation[first] = projected / np.linalg.norm(projected)
            rotation[second] = np.cross(normal, rotation[first])
            target = second
        if np.linalg.det(rotation) < 0:
            rotation[target] *= -1.0
        return rotation
    if np.linalg.det(rotation) < 0:
        rotation[0] *= -1.0
    return rotation


def fixed_solute_alignment(
    positions_nm: NDArray,
    masses: NDArray,
    atomic_numbers: NDArray,
    exclude: NDArray | None = None,
) -> tuple[NDArray, NDArray]:
    """Rotation and centre that align the fixed solutes for a rectangular box.

    The inertia is that of the heavy atoms, mass-weighted, leaving out the
    atoms flagged in ``exclude``: PolyzyMD excludes a modified protein's
    attachment segments (``G###``, ``P###``, ``M###``), whose grafted chains
    can be drawn per replicate, so every replicate of a condition is turned
    the same way and shares one box.

    Returns
    -------
    (rotation, center_nm)
        Apply as ``(positions - center) @ rotation.T + center``.
    """
    positions = np.asarray(positions_nm, dtype=float)
    keep = np.asarray(atomic_numbers) > 1
    if exclude is not None:
        keep &= ~np.asarray(exclude, dtype=bool)
    weights = np.asarray(masses, dtype=float)[keep]
    rotation = principal_axes_rotation(positions[keep], weights)
    center = (positions[keep] * weights[:, None]).sum(axis=0) / weights.sum()
    return rotation, center


def describe_box_plan(plan: dict) -> str:
    """Return a one-line description of a :func:`plan_box` result."""
    if plan.get("shape") == "rectangular":
        return (
            "rectangular %.2f x %.2f x %.2f nm (solute bounding box %.2f x %.2f x %.2f nm "
            "+ 2 x padding %.2f / %.2f / %.2f nm)"
            % (*plan["brick"], *plan["extent"], *plan["padding_xyz"])
        )
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
