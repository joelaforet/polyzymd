"""
Residue-based periodic boundary condition (PBC) unwrapping.

This module provides mathematically rigorous tools to "unwrap" molecules that have
been split across periodic boundaries in molecular dynamics trajectories. Unlike
bond-based approaches (which fail when connectivity information is unavailable or
corrupted), this module uses **residue membership** to determine which atoms belong
to the same molecule.

Theory of Periodic Boundary Conditions
======================================

In MD simulations, the simulation box is replicated infinitely in all directions
(periodic images). When a molecule moves across a box boundary, it reappears on
the opposite side. This creates visual artifacts where molecules appear "split"
across the box.

The simulation box is defined by three lattice vectors **a**, **b**, **c** and
angles α (between b and c), β (between a and c), γ (between a and b). These are
stored in MDAnalysis as::

    dimensions = [a, b, c, alpha, beta, gamma]  # lengths in Å, angles in degrees

Box Matrix Construction
-----------------------

The box matrix **H** transforms fractional coordinates to Cartesian coordinates::

    r = H · s

where:
- **r** is a position in Cartesian coordinates (Å)
- **s** is a position in fractional coordinates (dimensionless, range [0, 1))
- **H** is the 3×3 box matrix (columns are the lattice vectors)

For a general triclinic box with crystallographic parameters (a, b, c, α, β, γ),
the box matrix in the standard orientation (a along x-axis, b in xy-plane) is::

    H = | a_x  b_x  c_x |   | a    b·cos(γ)    c·cos(β)                        |
        | 0    b_y  c_y | = | 0    b·sin(γ)    c·[cos(α)-cos(β)cos(γ)]/sin(γ)  |
        | 0    0    c_z |   | 0    0           c·√(1 - cos²α - cos²β - cos²γ   |
                            |                      + 2·cos(α)cos(β)cos(γ)) / sin(γ) |

For orthogonal boxes (α = β = γ = 90°), this simplifies to::

    H = | a  0  0 |
        | 0  b  0 |
        | 0  0  c |

Minimum Image Convention
------------------------

To find the shortest vector between two atoms across periodic boundaries, we use
the Minimum Image Convention (MIC). Given a displacement vector Δr:

1. Convert to fractional coordinates: Δs = H⁻¹ · Δr
2. Apply MIC: Δs_mic = Δs - round(Δs)    # wraps to [-0.5, 0.5)
3. Convert back: Δr_mic = H · Δs_mic

The result Δr_mic is the shortest displacement vector considering all periodic
images.

Residue-Based Unwrapping Algorithm
----------------------------------

For each residue (which corresponds to one molecule in PolyzyMD's PDB output):

1. Use the first atom as the reference position
2. For each subsequent atom:
   a. Calculate displacement from previous atom: Δr = r_i - r_{i-1}
   b. Apply MIC to get shortest displacement: Δr_mic = MIC(Δr)
   c. Place atom at: r_i_unwrapped = r_{i-1}_unwrapped + Δr_mic

This "chains" atoms together, ensuring the molecule is continuous even if it
spans multiple periodic images.

Why Residue-Based Works for PolyzyMD
------------------------------------

PolyzyMD's `SystemBuilder._assign_pdb_identifiers()` ensures that each molecule
(protein, ligand, individual polymer, water, DMSO) gets its own unique residue
number. This means:

- All atoms of a water molecule share the same resid
- All atoms of a DMSO molecule share the same resid
- All atoms of a polymer chain share the same resid

Even when CONECT records are corrupted (atom indices > 99,999 overflow PDB format),
the residue assignments remain valid, allowing correct unwrapping.

References
----------
.. [1] Allen, M.P. & Tildesley, D.J. (2017). Computer Simulation of Liquids.
       Oxford University Press. Chapter 1.5: Periodic Boundary Conditions.
.. [2] Tuckerman, M.E. (2010). Statistical Mechanics: Theory and Molecular
       Simulation. Oxford University Press. Section 3.3.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    import MDAnalysis as mda
    from MDAnalysis.core.groups import AtomGroup

LOGGER = logging.getLogger(__name__)


def box_dimensions_to_matrix(dimensions: np.ndarray) -> np.ndarray:
    """
    Convert MDAnalysis box dimensions to a 3×3 box matrix.

    Transforms crystallographic parameters [a, b, c, α, β, γ] into the box matrix
    **H** whose columns are the lattice vectors. This matrix converts fractional
    coordinates to Cartesian coordinates: r = H · s.

    Parameters
    ----------
    dimensions : np.ndarray
        Box dimensions as [a, b, c, alpha, beta, gamma] where:
        - a, b, c are box edge lengths in Ångströms
        - alpha, beta, gamma are box angles in degrees:
          - α (alpha): angle between b and c
          - β (beta): angle between a and c
          - γ (gamma): angle between a and b

    Returns
    -------
    np.ndarray
        3×3 box matrix H with lattice vectors as columns::

            H = | a_x  b_x  c_x |
                | a_y  b_y  c_y |
                | a_z  b_z  c_z |

    Notes
    -----
    The box matrix is constructed in the standard crystallographic orientation:

    - Vector **a** lies along the positive x-axis
    - Vector **b** lies in the xy-plane (positive y component)
    - Vector **c** has components in all three directions

    The explicit formulas are::

        a_x = a
        b_x = b · cos(γ)
        b_y = b · sin(γ)
        c_x = c · cos(β)
        c_y = c · [cos(α) - cos(β)·cos(γ)] / sin(γ)
        c_z = √(c² - c_x² - c_y²)

    For orthogonal boxes (α = β = γ = 90°), this reduces to a diagonal matrix.

    Examples
    --------
    >>> # Cubic box with 50 Å edges
    >>> dims = np.array([50.0, 50.0, 50.0, 90.0, 90.0, 90.0])
    >>> H = box_dimensions_to_matrix(dims)
    >>> H
    array([[50.,  0.,  0.],
           [ 0., 50.,  0.],
           [ 0.,  0., 50.]])

    >>> # Triclinic box (truncated octahedron)
    >>> dims = np.array([66.433, 76.808, 78.786, 60.83, 65.06, 90.0])
    >>> H = box_dimensions_to_matrix(dims)
    """
    a, b, c = dimensions[:3]
    # Convert angles from degrees to radians
    alpha, beta, gamma = np.radians(dimensions[3:6])

    # Precompute trigonometric values
    cos_alpha = np.cos(alpha)
    cos_beta = np.cos(beta)
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)

    # Build box matrix (columns are lattice vectors a, b, c)
    # Standard orientation: a along x, b in xy-plane
    ax = a
    bx = b * cos_gamma
    by = b * sin_gamma
    cx = c * cos_beta
    cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz = np.sqrt(c**2 - cx**2 - cy**2)

    # Construct matrix with lattice vectors as columns
    H = np.array(
        [
            [ax, bx, cx],
            [0.0, by, cy],
            [0.0, 0.0, cz],
        ]
    )

    return H


def minimum_image_displacement(
    dr: np.ndarray,
    box_matrix: np.ndarray,
    box_matrix_inv: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Apply the Minimum Image Convention to a displacement vector.

    Given a displacement vector Δr between two atoms, returns the shortest
    equivalent displacement considering all periodic images. This is essential
    for correctly handling atoms that may be on opposite sides of the box.

    Parameters
    ----------
    dr : np.ndarray
        Displacement vector(s) in Cartesian coordinates (Å).
        Shape: (3,) for single vector, or (N, 3) for N vectors.
    box_matrix : np.ndarray
        3×3 box matrix from `box_dimensions_to_matrix()`.
    box_matrix_inv : np.ndarray, optional
        Pre-computed inverse of box_matrix for efficiency when processing
        many displacements with the same box. If None, computed internally.

    Returns
    -------
    np.ndarray
        Minimum image displacement vector(s), same shape as input.
        Each component is in the range [-L/2, L/2) where L is the box length
        in that direction.

    Notes
    -----
    The algorithm works in fractional coordinates where the MIC is trivial:

    1. Convert to fractional: Δs = H⁻¹ · Δr
    2. Apply MIC: Δs_mic = Δs - round(Δs)
       - This maps each component to the range [-0.5, 0.5)
    3. Convert back: Δr_mic = H · Δs_mic

    The rounding operation `round(Δs)` finds the nearest integer, which
    corresponds to the nearest periodic image. Subtracting it gives the
    shortest path.

    For a single component in an orthogonal box::

        Δx_mic = Δx - L_x · round(Δx / L_x)

    Examples
    --------
    >>> H = np.diag([50.0, 50.0, 50.0])  # 50 Å cubic box
    >>> dr = np.array([48.0, 2.0, -3.0])  # atom nearly across box in x
    >>> minimum_image_displacement(dr, H)
    array([-2.,  2., -3.])  # shortest path wraps around
    """
    if box_matrix_inv is None:
        box_matrix_inv = np.linalg.inv(box_matrix)

    # Handle both single vectors and arrays of vectors
    single_vector = dr.ndim == 1
    if single_vector:
        dr = dr.reshape(1, 3)

    # Transform to fractional coordinates: s = H^{-1} · r
    # For array of vectors, use transpose: (N, 3) @ (3, 3).T = (N, 3)
    ds = dr @ box_matrix_inv.T

    # Apply minimum image convention in fractional space
    # round() maps to nearest integer, subtracting wraps to [-0.5, 0.5)
    ds_mic = ds - np.round(ds)

    # Transform back to Cartesian: r = H · s
    dr_mic = ds_mic @ box_matrix.T

    if single_vector:
        return dr_mic[0]
    return dr_mic


def unwrap_residue(
    positions: np.ndarray,
    box_matrix: np.ndarray,
    box_matrix_inv: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Make a single residue/molecule whole by unwrapping across periodic boundaries.

    Takes the atomic positions of a single residue and adjusts them so the molecule
    is continuous (not split across the box). The first atom is used as the reference
    point; subsequent atoms are placed at the minimum image distance from the
    previous atom.

    Parameters
    ----------
    positions : np.ndarray
        Atomic positions of shape (N_atoms, 3) in Cartesian coordinates (Å).
    box_matrix : np.ndarray
        3×3 box matrix from `box_dimensions_to_matrix()`.
    box_matrix_inv : np.ndarray, optional
        Pre-computed inverse of box_matrix for efficiency.

    Returns
    -------
    np.ndarray
        Unwrapped positions of shape (N_atoms, 3). The first atom position is
        unchanged; subsequent atoms are shifted to be continuous with the molecule.

    Notes
    -----
    The algorithm "chains" atoms together:

    1. Keep first atom at original position: r₀_unwrapped = r₀
    2. For each subsequent atom i:
       a. Calculate displacement from previous: Δr = rᵢ - rᵢ₋₁_unwrapped
       b. Apply MIC: Δr_mic = MIC(Δr)
       c. Place atom: rᵢ_unwrapped = rᵢ₋₁_unwrapped + Δr_mic

    This ensures each atom is at the minimum image distance from its predecessor,
    resulting in a continuous molecule.

    **Why chain from previous atom?**

    For small molecules (water, DMSO), chaining from the first atom would work.
    But for polymers with many atoms, an atom far down the chain might be more
    than half a box length from the first atom (legitimate, since the polymer
    is large). Chaining from the previous atom handles arbitrarily long molecules.

    Examples
    --------
    >>> # Water molecule split across periodic boundary
    >>> H = np.diag([50.0, 50.0, 50.0])
    >>> # Oxygen at x=1, hydrogens wrapped to x=49 (should be at x=-1)
    >>> positions = np.array([
    ...     [1.0, 25.0, 25.0],   # O
    ...     [49.0, 25.5, 25.0],  # H1 (wrapped, should be ~-1)
    ...     [49.0, 24.5, 25.0],  # H2 (wrapped, should be ~-1)
    ... ])
    >>> unwrapped = unwrap_residue(positions, H)
    >>> unwrapped
    array([[ 1.,  25.,  25.],
           [-1.,  25.5, 25.],
           [-1.,  24.5, 25.]])
    """
    if box_matrix_inv is None:
        box_matrix_inv = np.linalg.inv(box_matrix)

    n_atoms = positions.shape[0]
    if n_atoms == 0:
        return positions.copy()

    unwrapped = np.empty_like(positions)
    unwrapped[0] = positions[0]  # First atom unchanged

    # Chain each atom to the previous one
    for i in range(1, n_atoms):
        dr = positions[i] - unwrapped[i - 1]
        dr_mic = minimum_image_displacement(dr, box_matrix, box_matrix_inv)
        unwrapped[i] = unwrapped[i - 1] + dr_mic

    return unwrapped


def unwrap_by_residue(
    atom_group: "AtomGroup",
    box_matrix: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Unwrap all molecules in an AtomGroup using residue membership.

    Iterates over all residues in the AtomGroup and unwraps each one independently.
    This is the main entry point for residue-based unwrapping.

    Parameters
    ----------
    atom_group : MDAnalysis.AtomGroup
        The atoms to unwrap. Must have valid residue assignments (each molecule
        should have a unique resid). The AtomGroup's Universe must have valid
        box dimensions.
    box_matrix : np.ndarray, optional
        Pre-computed 3×3 box matrix. If None, computed from the Universe's
        current box dimensions. Providing this is useful when processing
        multiple frames with the same box.

    Returns
    -------
    np.ndarray
        Unwrapped positions of shape (N_atoms, 3), in the same atom order as
        the input AtomGroup.

    Warns
    -----
    Logs a warning and skips residues that contain NaN positions.

    Notes
    -----
    This function assumes that PolyzyMD's residue assignment convention is used:

    - Each molecule (water, DMSO, polymer, etc.) has a unique residue number
    - All atoms of a molecule share the same resid

    This allows correct unwrapping even when bond connectivity information is
    unavailable or corrupted (e.g., PDB CONECT records with atom indices > 99,999).

    The function modifies positions in-place conceptually but returns a new array.
    To apply to a trajectory, use with MDAnalysis transformations.

    See Also
    --------
    ResidueUnwrapTransform : MDAnalysis transformation wrapper for trajectories.

    Examples
    --------
    >>> import MDAnalysis as mda
    >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
    >>> # Unwrap all atoms in current frame
    >>> unwrapped_pos = unwrap_by_residue(u.atoms)
    >>> u.atoms.positions = unwrapped_pos
    """
    # Get box matrix from universe if not provided
    if box_matrix is None:
        dimensions = atom_group.universe.dimensions
        if dimensions is None:
            raise ValueError("Universe has no box dimensions. Cannot unwrap without periodic box.")
        box_matrix = box_dimensions_to_matrix(dimensions)

    box_matrix_inv = np.linalg.inv(box_matrix)

    # Get current positions
    positions = atom_group.positions.copy()

    # Build index mapping: for each atom in atom_group, what's its index in positions?
    # This handles the case where atom_group is a subset of the universe
    atom_indices = {atom.ix: i for i, atom in enumerate(atom_group)}

    # Process each residue
    for residue in atom_group.residues:
        # Get atoms of this residue that are in our atom_group
        res_atoms = residue.atoms.intersection(atom_group)
        if len(res_atoms) == 0:
            continue

        # Get indices into our positions array
        indices = [atom_indices[atom.ix] for atom in res_atoms]
        res_positions = positions[indices]

        # Check for NaN positions
        if np.any(np.isnan(res_positions)):
            LOGGER.warning(
                f"Skipping residue {residue.resname}:{residue.resid} - contains NaN positions"
            )
            continue

        # Unwrap this residue
        unwrapped = unwrap_residue(res_positions, box_matrix, box_matrix_inv)

        # Write back
        positions[indices] = unwrapped

    return positions


class ResidueUnwrapTransform:
    """
    MDAnalysis transformation for residue-based PBC unwrapping.

    This class wraps `unwrap_by_residue()` as an MDAnalysis on-the-fly transformation,
    allowing it to be added to a trajectory's transformation pipeline.

    Parameters
    ----------
    atom_group : MDAnalysis.AtomGroup
        The atoms to unwrap. Typically `universe.atoms` for all atoms.

    Notes
    -----
    Unlike MDAnalysis's built-in `unwrap()` transformation which requires bond
    information, this transformation uses residue membership to determine molecular
    boundaries. This is essential for PolyzyMD trajectories where:

    1. PDB CONECT records may be corrupted (atom indices > 99,999)
    2. Each molecule is guaranteed to have a unique residue number

    Examples
    --------
    >>> import MDAnalysis as mda
    >>> from polyzymd.analysis.unwrap import ResidueUnwrapTransform
    >>>
    >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
    >>> transform = ResidueUnwrapTransform(u.atoms)
    >>> u.trajectory.add_transformations(transform)
    >>>
    >>> # Now iterating through frames automatically applies unwrapping
    >>> for ts in u.trajectory:
    ...     # Positions are already unwrapped
    ...     print(u.atoms.positions)
    """

    def __init__(self, atom_group: "AtomGroup") -> None:
        self.atom_group = atom_group
        self._box_matrix: Optional[np.ndarray] = None
        self._last_dimensions: Optional[np.ndarray] = None

    def __call__(self, ts):
        """
        Apply residue-based unwrapping to the current timestep.

        Parameters
        ----------
        ts : MDAnalysis.coordinates.base.Timestep
            The current timestep object.

        Returns
        -------
        ts : MDAnalysis.coordinates.base.Timestep
            The modified timestep with unwrapped positions.
        """
        # Cache box matrix if dimensions haven't changed (common for NVT/NPT)
        dimensions = ts.dimensions
        if (
            self._box_matrix is None
            or self._last_dimensions is None
            or not np.allclose(dimensions, self._last_dimensions)
        ):
            self._box_matrix = box_dimensions_to_matrix(dimensions)
            self._last_dimensions = dimensions.copy()

        # Unwrap all atoms
        unwrapped_positions = unwrap_by_residue(self.atom_group, self._box_matrix)

        # Update positions in-place
        self.atom_group.positions = unwrapped_positions

        return ts
