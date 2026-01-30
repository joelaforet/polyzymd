"""
Cluster-based trajectory visualization for periodic boundary conditions.

This module provides tools to create "droplet" visualizations of solvated MD
trajectories, where all molecules cluster around a central reference group
(typically protein + ligand). This approach is designed for **visualization**,
not for quantitative analysis of diffusion or other PBC-sensitive properties.

The Problem
===========

In MD simulations with periodic boundary conditions (PBC), molecules can:
1. Diffuse across box boundaries, appearing on the opposite side
2. Become "split" across boundaries (part of molecule on each side)
3. Scatter throughout periodic images, making visualization confusing

Standard unwrapping approaches rely on **bond connectivity** to determine which
atoms belong to the same molecule. However, PDB CONECT records have a 5-digit
atom index limit (99,999). For large solvated systems, atom indices overflow,
corrupting the CONECT records and breaking bond-based methods.

The Solution: Cluster Around Protein
====================================

PolyzyMD's approach leverages two key facts:
1. **Residue membership is preserved**: Each molecule has a unique residue number
2. **We want a "droplet" visualization**: Molecules should cluster around protein

Algorithm:
1. Calculate reference center of mass (protein + ligand)
2. Center the reference at the box center
3. For each molecule (identified by residue):
   a. Compute the molecule's center of mass
   b. Find the periodic image that places the COM closest to the reference
   c. Translate the entire molecule to that image

This creates a visualization where:
- The protein is centered and whole
- All solvent/polymer molecules are whole (not split across boundaries)
- Everything clusters around the protein in a "droplet" shape
- Box boundaries are effectively ignored for visualization purposes

Mathematical Background
=======================

Minimum Image Convention
------------------------

For a displacement vector Δr between two points in a periodic box:

1. Convert to fractional coordinates: Δs = H⁻¹ · Δr
2. Apply minimum image: Δs_mic = Δs - round(Δs)
3. Convert back to Cartesian: Δr_mic = H · Δs_mic

Where H is the 3×3 box matrix with lattice vectors as columns.

Box Matrix from Crystallographic Parameters
-------------------------------------------

For a box with edge lengths (a, b, c) and angles (α, β, γ)::

    H = | a    b·cos(γ)    c·cos(β)                        |
        | 0    b·sin(γ)    c·[cos(α)-cos(β)cos(γ)]/sin(γ)  |
        | 0    0           c·√(1 - cos²α - cos²β - cos²γ   |
                               + 2·cos(α)cos(β)cos(γ)) / sin(γ) |

Vectorized Center of Mass Calculation
-------------------------------------

For N atoms with positions r_i, masses m_i, and residue indices k_i, the
center of mass of residue K is::

    COM_K = Σ_{i: k_i = K} (m_i · r_i) / Σ_{i: k_i = K} m_i

This is computed efficiently using NumPy's bincount for the weighted sums.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional, Tuple

import numpy as np

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup

LOGGER = logging.getLogger(__name__)


def box_dimensions_to_matrix(dimensions: np.ndarray) -> np.ndarray:
    """
    Convert MDAnalysis box dimensions to a 3×3 box matrix.

    Parameters
    ----------
    dimensions : np.ndarray
        Box dimensions as [a, b, c, alpha, beta, gamma] where lengths are in
        Ångströms and angles are in degrees.

    Returns
    -------
    np.ndarray
        3×3 box matrix H with lattice vectors as columns.

    Notes
    -----
    The box matrix transforms fractional coordinates to Cartesian: r = H · s
    """
    a, b, c = dimensions[:3]
    alpha, beta, gamma = np.radians(dimensions[3:6])

    cos_alpha = np.cos(alpha)
    cos_beta = np.cos(beta)
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)

    ax = a
    bx = b * cos_gamma
    by = b * sin_gamma
    cx = c * cos_beta
    cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz = np.sqrt(c**2 - cx**2 - cy**2)

    return np.array(
        [
            [ax, bx, cx],
            [0.0, by, cy],
            [0.0, 0.0, cz],
        ]
    )


def minimum_image_shift(
    displacement: np.ndarray,
    box_matrix: np.ndarray,
    box_matrix_inv: np.ndarray,
) -> np.ndarray:
    """
    Compute the shift needed to place a point at its minimum image position.

    Given a displacement from a reference point, computes how much to translate
    so the point is at the minimum image distance from the reference.

    Parameters
    ----------
    displacement : np.ndarray
        Displacement vector(s) from reference. Shape (3,) or (N, 3).
    box_matrix : np.ndarray
        3×3 box matrix.
    box_matrix_inv : np.ndarray
        Inverse of box matrix (precomputed for efficiency).

    Returns
    -------
    np.ndarray
        Shift vector(s) to add to positions. Same shape as input.
    """
    # Convert to fractional coordinates
    if displacement.ndim == 1:
        ds = box_matrix_inv @ displacement
    else:
        ds = displacement @ box_matrix_inv.T

    # Find which periodic image to use
    image_shift = -np.round(ds)

    # Convert back to Cartesian
    if displacement.ndim == 1:
        return box_matrix @ image_shift
    else:
        return image_shift @ box_matrix.T


def compute_residue_coms(
    positions: np.ndarray,
    resindices: np.ndarray,
    masses: Optional[np.ndarray] = None,
    n_residues: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute center of mass for each residue using vectorized operations.

    Parameters
    ----------
    positions : np.ndarray
        Atomic positions, shape (N_atoms, 3).
    resindices : np.ndarray
        Residue index for each atom, shape (N_atoms,).
    masses : np.ndarray, optional
        Atomic masses, shape (N_atoms,). If None, uses geometric center.
    n_residues : int, optional
        Number of residues. If None, computed from max(resindices) + 1.

    Returns
    -------
    coms : np.ndarray
        Center of mass for each residue, shape (N_residues, 3).
    total_masses : np.ndarray
        Total mass of each residue, shape (N_residues,).
        If masses not provided, this is the atom count per residue.

    Notes
    -----
    Uses np.bincount for O(N_atoms) complexity instead of O(N_atoms × N_residues).
    """
    if n_residues is None:
        n_residues = resindices.max() + 1

    if masses is None:
        # Geometric center - equal weights
        masses = np.ones(len(positions))
        LOGGER.debug("No masses provided, using geometric center")

    # Weighted position sums for each residue: Σ(m_i * r_i) for each residue
    weighted_pos = positions * masses[:, np.newaxis]  # (N_atoms, 3)

    # Sum weighted positions by residue index
    com_x = np.bincount(resindices, weights=weighted_pos[:, 0], minlength=n_residues)
    com_y = np.bincount(resindices, weights=weighted_pos[:, 1], minlength=n_residues)
    com_z = np.bincount(resindices, weights=weighted_pos[:, 2], minlength=n_residues)

    # Total mass per residue
    total_masses = np.bincount(resindices, weights=masses, minlength=n_residues)

    # Avoid division by zero for empty residues (shouldn't happen, but be safe)
    total_masses = np.where(total_masses == 0, 1.0, total_masses)

    # COM = Σ(m_i * r_i) / Σ(m_i)
    coms = np.column_stack([com_x, com_y, com_z]) / total_masses[:, np.newaxis]

    return coms, total_masses


def cluster_around_point(
    positions: np.ndarray,
    resindices: np.ndarray,
    reference_point: np.ndarray,
    box_matrix: np.ndarray,
    masses: Optional[np.ndarray] = None,
    exclude_resindices: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Translate each residue to its minimum image position relative to a reference point.

    This is the core algorithm for creating "droplet" visualizations. Each residue
    (molecule) is translated as a unit so its center of mass is at the minimum
    image distance from the reference point.

    Parameters
    ----------
    positions : np.ndarray
        Atomic positions, shape (N_atoms, 3). Modified in-place.
    resindices : np.ndarray
        Residue index for each atom, shape (N_atoms,).
    reference_point : np.ndarray
        The point to cluster around (e.g., protein COM), shape (3,).
    box_matrix : np.ndarray
        3×3 box matrix.
    masses : np.ndarray, optional
        Atomic masses for COM calculation. If None, uses geometric center.
    exclude_resindices : np.ndarray, optional
        Residue indices to skip (e.g., protein residues already at reference).

    Returns
    -------
    np.ndarray
        Modified positions array (same object as input, modified in-place).

    Notes
    -----
    Algorithm complexity: O(N_atoms) - single pass over all atoms.

    The algorithm:
    1. Compute COM of each residue: O(N_atoms)
    2. Compute displacement from reference: O(N_residues)
    3. Compute minimum image shift: O(N_residues)
    4. Apply shifts to atoms via broadcasting: O(N_atoms)
    """
    n_residues = resindices.max() + 1
    box_matrix_inv = np.linalg.inv(box_matrix)

    # Step 1: Compute COM for each residue
    residue_coms, _ = compute_residue_coms(positions, resindices, masses, n_residues)

    # Step 2: Displacement of each residue COM from reference
    displacements = residue_coms - reference_point  # (N_residues, 3)

    # Step 3: Compute shift needed to bring each residue to minimum image
    shifts = minimum_image_shift(displacements, box_matrix, box_matrix_inv)

    # Step 4: Zero out shifts for excluded residues (protein/ligand)
    if exclude_resindices is not None and len(exclude_resindices) > 0:
        shifts[exclude_resindices] = 0.0

    # Step 5: Apply shifts to atoms (broadcast from residue to atoms)
    atom_shifts = shifts[resindices]  # (N_atoms, 3)
    positions += atom_shifts

    return positions


class ClusterAroundProtein:
    """
    MDAnalysis transformation that clusters all molecules around protein+ligand.

    Creates a "droplet" visualization where:
    - Protein and ligand are centered at the box center
    - All other molecules are whole and clustered around the protein
    - Box boundaries are effectively ignored

    This transformation is designed for **visualization only**, not for
    quantitative analysis of diffusion or other PBC-sensitive properties.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The Universe containing the trajectory.
    protein_selection : str
        MDAnalysis selection string for protein atoms (e.g., "protein").
    ligand_selection : str, optional
        MDAnalysis selection string for ligand atoms (e.g., "resname LIG").
        If provided, ligand is grouped with protein for centering.

    Examples
    --------
    >>> import MDAnalysis as mda
    >>> from polyzymd.analysis.unwrap import ClusterAroundProtein
    >>>
    >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
    >>> transform = ClusterAroundProtein(u, "protein", "resname SUB")
    >>> u.trajectory.add_transformations(transform)
    >>>
    >>> # Now iterate - positions are automatically transformed
    >>> for ts in u.trajectory:
    ...     # Positions show protein-centered droplet
    ...     pass
    """

    def __init__(
        self,
        universe,
        protein_selection: str = "protein",
        ligand_selection: Optional[str] = None,
    ) -> None:
        self.universe = universe
        self.all_atoms = universe.atoms

        # Select reference group (protein + optional ligand)
        self.protein_ag = universe.select_atoms(protein_selection)
        if len(self.protein_ag) == 0:
            raise ValueError(f"Protein selection '{protein_selection}' matched no atoms")

        if ligand_selection:
            self.ligand_ag = universe.select_atoms(ligand_selection)
            if len(self.ligand_ag) == 0:
                LOGGER.warning(f"Ligand selection '{ligand_selection}' matched no atoms")
                self.reference_ag = self.protein_ag
            else:
                self.reference_ag = self.protein_ag | self.ligand_ag
        else:
            self.ligand_ag = None
            self.reference_ag = self.protein_ag

        # Get residue indices for the reference group (to exclude from clustering)
        self.reference_resindices = np.unique(self.reference_ag.resindices)

        # Check if masses are available
        try:
            _ = self.all_atoms.masses
            self._has_masses = True
        except Exception:
            self._has_masses = False
            LOGGER.warning(
                "Masses not available in topology. Using geometric centers instead of "
                "center of mass. This may slightly affect clustering for asymmetric molecules."
            )

        # Cache for box matrix (recomputed if box changes)
        self._cached_box_matrix: Optional[np.ndarray] = None
        self._cached_dimensions: Optional[np.ndarray] = None

        LOGGER.info(
            f"ClusterAroundProtein initialized: "
            f"{len(self.protein_ag)} protein atoms, "
            f"{len(self.ligand_ag) if self.ligand_ag else 0} ligand atoms, "
            f"{len(self.reference_resindices)} reference residues"
        )

    def _get_box_matrix(self, dimensions: np.ndarray) -> np.ndarray:
        """Get box matrix, using cache if dimensions unchanged."""
        if self._cached_dimensions is None or not np.allclose(dimensions, self._cached_dimensions):
            self._cached_box_matrix = box_dimensions_to_matrix(dimensions)
            self._cached_dimensions = dimensions.copy()
        return self._cached_box_matrix

    def __call__(self, ts):
        """
        Apply clustering transformation to the current timestep.

        Parameters
        ----------
        ts : MDAnalysis.coordinates.base.Timestep
            The current timestep.

        Returns
        -------
        ts : MDAnalysis.coordinates.base.Timestep
            The modified timestep.
        """
        dimensions = ts.dimensions
        if dimensions is None:
            raise ValueError("Trajectory has no box dimensions. Cannot apply PBC clustering.")

        box_matrix = self._get_box_matrix(dimensions)
        box_center = np.sum(box_matrix, axis=1) / 2  # Center of the box

        # Get current positions and residue indices
        positions = self.all_atoms.positions
        resindices = self.all_atoms.resindices
        masses = self.all_atoms.masses if self._has_masses else None

        # Step 1: Compute reference (protein+ligand) center of mass
        ref_positions = self.reference_ag.positions
        if self._has_masses:
            ref_masses = self.reference_ag.masses
            ref_com = np.average(ref_positions, axis=0, weights=ref_masses)
        else:
            ref_com = np.mean(ref_positions, axis=0)

        # Step 2: Center the entire system so reference COM is at box center
        translation = box_center - ref_com
        positions += translation

        # Step 3: Cluster all molecules around the (now centered) reference
        # The reference is now at box_center
        cluster_around_point(
            positions=positions,
            resindices=resindices,
            reference_point=box_center,
            box_matrix=box_matrix,
            masses=masses,
            exclude_resindices=self.reference_resindices,
        )

        # Update positions in the AtomGroup
        self.all_atoms.positions = positions

        return ts
