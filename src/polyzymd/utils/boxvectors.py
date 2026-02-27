"""
Box vector utilities for simulation setup.

This module provides functions for calculating and manipulating periodic box
vectors used in molecular dynamics simulations:

- get_topology_bbox: Calculate the bounding box of all molecules in a topology
- pad_box_vectors_uniform: Add uniform padding to box vectors
- get_box_volume: Calculate the volume enclosed by box vectors

These functions are essential for setting up properly-sized simulation boxes
that accommodate the solute with appropriate padding for periodic boundary
conditions.

Example:
    >>> from openff.toolkit import Topology
    >>> from openff.units import Quantity
    >>> from polyzymd.utils.boxvectors import (
    ...     get_topology_bbox,
    ...     pad_box_vectors_uniform,
    ...     get_box_volume,
    ... )
    >>>
    >>> topology = Topology.from_molecules([molecule])
    >>> bbox = get_topology_bbox(topology)
    >>> padded = pad_box_vectors_uniform(bbox, Quantity(1.0, "nanometer"))
    >>> volume = get_box_volume(padded)

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

    Unlike :func:`get_topology_bbox`, which returns a 3x3 diagonal matrix of
    box *dimensions*, this function returns the raw per-axis minimum and
    maximum coordinates.  This is needed when the absolute position of the
    bounding box matters (e.g. computing an exclusion zone for Packmol).

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
    all_positions: list[NDArray] = []

    for molecule in topology.molecules:
        if molecule.n_conformers > 0:
            coords_angstrom = molecule.conformers[0].m_as("angstrom")
            all_positions.append(coords_angstrom)

    if not all_positions:
        raise ValueError(
            "Cannot calculate bounding box bounds: no molecules in topology have coordinates"
        )

    all_coords = np.vstack(all_positions)
    min_coords: NDArray = np.min(all_coords, axis=0)
    max_coords: NDArray = np.max(all_coords, axis=0)

    return min_coords, max_coords


def get_topology_bbox(topology: "Topology") -> "Quantity":
    """Calculate the bounding box of all molecules in a topology.

    This function finds the minimum axis-aligned bounding box that contains
    all atoms in the topology. The result is returned as a 3x3 diagonal
    matrix representing box vectors.

    Args:
        topology: OpenFF Topology with at least one molecule that has
            conformer coordinates.

    Returns:
        A 3x3 diagonal matrix (as an OpenFF Quantity with length units)
        where the diagonal elements are the box dimensions [Lx, Ly, Lz].

    Raises:
        ValueError: If the topology has no molecules with coordinates.

    Example:
        >>> bbox = get_topology_bbox(topology)
        >>> print(bbox.magnitude)
        [[5.2, 0.0, 0.0],
         [0.0, 4.8, 0.0],
         [0.0, 0.0, 6.1]]
    """
    from openff.units import Quantity

    # Collect all atomic positions from all molecules
    all_positions = []

    for molecule in topology.molecules:
        if molecule.n_conformers > 0:
            # Get the first conformer's coordinates
            coords = molecule.conformers[0]
            # Convert to angstrom for consistent processing
            coords_angstrom = coords.m_as("angstrom")
            all_positions.append(coords_angstrom)

    if not all_positions:
        raise ValueError("Cannot calculate bounding box: no molecules in topology have coordinates")

    # Stack all positions into a single array
    all_coords = np.vstack(all_positions)

    # Find min and max along each axis
    min_coords = np.min(all_coords, axis=0)
    max_coords = np.max(all_coords, axis=0)

    # Calculate box dimensions
    box_dims = max_coords - min_coords

    # Create a diagonal box vector matrix
    # OpenFF uses row vectors, so each row is a box vector
    box_vectors = np.diag(box_dims)

    return Quantity(box_vectors, "angstrom")


def pad_box_vectors_uniform(
    box_vectors: "Quantity",
    padding: "Quantity",
) -> "Quantity":
    """Add uniform padding to box vectors.

    This function increases the box size by adding the same padding amount
    to each dimension. This is useful for ensuring adequate space between
    the solute and the periodic boundary.

    Args:
        box_vectors: A 3x3 matrix of box vectors (with length units).
        padding: The amount of padding to add to each dimension (with
            length units). This value is added to EACH side, so the total
            increase in each dimension is 2 * padding.

    Returns:
        Padded box vectors as a 3x3 matrix with the same units as input.

    Example:
        >>> bbox = get_topology_bbox(topology)
        >>> padded = pad_box_vectors_uniform(bbox, Quantity(1.0, "nm"))
        >>> # If bbox was 5x4x6 nm, padded will be 7x6x8 nm
    """
    from openff.units import Quantity

    # Convert both to the same units (angstrom for precision)
    box_angstrom = box_vectors.m_as("angstrom")
    padding_angstrom = padding.m_as("angstrom")

    # For a diagonal box, add 2*padding to each diagonal element
    # (padding on each side of each dimension)
    padded = box_angstrom.copy()

    # Add padding to diagonal elements
    for i in range(3):
        padded[i, i] += 2 * padding_angstrom

    # Return in the original units
    return Quantity(padded, "angstrom")


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


def box_vectors_from_lengths_angles(
    a: float,
    b: float,
    c: float,
    alpha: float = 90.0,
    beta: float = 90.0,
    gamma: float = 90.0,
    unit: str = "nanometer",
) -> "Quantity":
    """Create box vectors from lengths and angles.

    This function constructs the 3x3 box vector matrix from crystallographic
    parameters (lengths a, b, c and angles alpha, beta, gamma).

    Args:
        a: Length of first box vector.
        b: Length of second box vector.
        c: Length of third box vector.
        alpha: Angle between b and c vectors in degrees (default 90).
        beta: Angle between a and c vectors in degrees (default 90).
        gamma: Angle between a and b vectors in degrees (default 90).
        unit: Length unit for a, b, c (default "nanometer").

    Returns:
        A 3x3 box vector matrix as an OpenFF Quantity.

    Example:
        >>> # Create a cubic box of 5 nm
        >>> box = box_vectors_from_lengths_angles(5.0, 5.0, 5.0)
        >>>
        >>> # Create a rhombic dodecahedron
        >>> box = box_vectors_from_lengths_angles(
        ...     5.0, 5.0, 5.0,
        ...     alpha=60.0, beta=60.0, gamma=90.0
        ... )
    """
    from openff.units import Quantity

    # Convert angles to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)

    # Calculate box vectors using standard crystallographic convention
    # Vector a is along x-axis
    # Vector b is in the xy-plane
    # Vector c has components along all three axes

    cos_alpha = np.cos(alpha_rad)
    cos_beta = np.cos(beta_rad)
    cos_gamma = np.cos(gamma_rad)
    sin_gamma = np.sin(gamma_rad)

    # First vector: along x
    ax = a
    ay = 0.0
    az = 0.0

    # Second vector: in xy plane
    bx = b * cos_gamma
    by = b * sin_gamma
    bz = 0.0

    # Third vector: general
    cx = c * cos_beta
    cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz = np.sqrt(c * c - cx * cx - cy * cy)

    box_vectors = np.array(
        [
            [ax, ay, az],
            [bx, by, bz],
            [cx, cy, cz],
        ]
    )

    return Quantity(box_vectors, unit)


def reduce_box_vectors(box_vectors: "Quantity") -> "Quantity":
    """Reduce box vectors to a standard form.

    This function applies lattice reduction to convert the box vectors
    to a form where the box is as close to orthogonal as possible while
    maintaining the same periodic images. This is useful for visualization
    and for some simulation algorithms.

    Args:
        box_vectors: A 3x3 matrix of box vectors (with length units).

    Returns:
        Reduced box vectors in standard form.

    Note:
        The reduced form ensures:
        - c_z > 0
        - b_z = 0
        - a_z = a_y = 0
        - -b_x/2 <= c_x < b_x/2
        - -a_x/2 <= b_x < a_x/2
        - -a_x/2 <= c_x < a_x/2
    """
    from openff.units import Quantity

    # Get magnitude and work in consistent units
    box = box_vectors.m_as("angstrom").copy()

    # Ensure the matrix is in lower triangular form
    # (standard convention for periodic boxes)
    a = box[0]
    b = box[1]
    c = box[2]

    # Reduce c with respect to b
    if b[1] != 0:
        c = c - b * np.round(c[1] / b[1])

    # Reduce c with respect to a
    if a[0] != 0:
        c = c - a * np.round(c[0] / a[0])

    # Reduce b with respect to a
    if a[0] != 0:
        b = b - a * np.round(b[0] / a[0])

    reduced = np.array([a, b, c])

    return Quantity(reduced, "angstrom")


def get_box_lengths(box_vectors: "Quantity") -> "Quantity":
    """Extract the lengths of the three box vectors.

    Args:
        box_vectors: A 3x3 matrix of box vectors (with length units).

    Returns:
        A 1D array of three lengths [|a|, |b|, |c|] with length units.

    Example:
        >>> lengths = get_box_lengths(box_vectors)
        >>> print(lengths.to("nanometer"))
        [5.0, 5.0, 5.0] nanometer
    """
    from openff.units import Quantity

    box = box_vectors.m_as("angstrom")

    lengths = np.array(
        [
            np.linalg.norm(box[0]),
            np.linalg.norm(box[1]),
            np.linalg.norm(box[2]),
        ]
    )

    return Quantity(lengths, "angstrom")


def get_box_angles(box_vectors: "Quantity") -> NDArray[np.float64]:
    """Extract the angles between box vectors in degrees.

    Args:
        box_vectors: A 3x3 matrix of box vectors (with length units).

    Returns:
        A 1D numpy array of three angles [alpha, beta, gamma] in degrees,
        where:
        - alpha is the angle between b and c
        - beta is the angle between a and c
        - gamma is the angle between a and b

    Example:
        >>> angles = get_box_angles(box_vectors)
        >>> print(angles)
        [90.0, 90.0, 90.0]
    """
    box = box_vectors.magnitude  # Just need the numbers, not units

    a = box[0]
    b = box[1]
    c = box[2]

    # Calculate angles using dot product formula
    # cos(theta) = (u . v) / (|u| |v|)
    def angle_between(v1: NDArray, v2: NDArray) -> float:
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        # Clip to handle numerical errors
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))

    alpha = angle_between(b, c)  # angle between b and c
    beta = angle_between(a, c)  # angle between a and c
    gamma = angle_between(a, b)  # angle between a and b

    return np.array([alpha, beta, gamma])
