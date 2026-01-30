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

The Solution: Unwrapped Trajectories + NoJump + Cluster
=======================================================

PolyzyMD writes trajectories with ``enforcePeriodicBox=False`` in OpenMM's
DCDReporter. This means:

1. **Molecules start whole**: PACKMOL creates continuous molecules
2. **No per-frame wrapping**: OpenMM writes raw coordinates without applying
   PBC wrapping on each frame, so molecules stay whole
3. **Atoms may drift**: Over long simulations, atoms can move far from origin

The post-processing transformation then applies:

**Phase 1: NoJump (all frames)**
    Prevent atoms from jumping more than half a box between frames due to
    any residual discontinuities. This maintains molecular continuity.

**Phase 2: Cluster Around Protein**
    Translate each whole molecule to be at minimum image distance from the
    protein center of mass, creating a "droplet" visualization.

.. note::

   This module requires trajectories written with ``enforcePeriodicBox=False``.
   Trajectories written with the default OpenMM behavior (per-frame wrapping)
   will NOT work correctly with this post-processor.

Mathematical Background
=======================

Minimum Image Convention
------------------------

For a displacement vector Δr between two points in a periodic box:

1. Convert to fractional coordinates: Δs = H⁻¹ · Δr
2. Apply minimum image: Δs_mic = Δs - round(Δs)
3. Convert back to Cartesian: Δr_mic = H · Δs_mic

Where H is the 3×3 box matrix with lattice vectors as columns.

NoJump Algorithm (Kulke & Vermaas 2022)
---------------------------------------

For each atom, compare current position to previous frame::

    f_current = positions @ box_inverse  # fractional coordinates
    f_new = f_current - round(f_current - f_prev)  # remove jumps
    positions = f_new @ box  # back to Cartesian

This prevents atoms from moving more than half a box between frames.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

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


def make_molecules_whole(
    positions: np.ndarray,
    resindices: np.ndarray,
    box_matrix: np.ndarray,
    box_matrix_inv: np.ndarray,
) -> None:
    """
    Make all molecules whole by chaining atoms within each residue.

    For each residue (molecule), atoms are repositioned so each atom is at the
    minimum image distance from the previous atom in the residue. This "chains"
    atoms together, ensuring the molecule is continuous even if it spans
    periodic boundaries.

    Parameters
    ----------
    positions : np.ndarray
        Atomic positions, shape (N_atoms, 3). **Modified in-place.**
    resindices : np.ndarray
        Residue index for each atom, shape (N_atoms,).
    box_matrix : np.ndarray
        3×3 box matrix (triclinic dimensions).
    box_matrix_inv : np.ndarray
        Inverse of box matrix (precomputed for efficiency).

    Notes
    -----
    Performance is optimized by batching residues of the same size together.
    Water molecules (3 atoms) are the most common and processed in a single
    vectorized operation.

    The algorithm for each residue:
    1. Keep first atom at original position
    2. For each subsequent atom:
       - Compute displacement from previous atom
       - Apply minimum image convention
       - Place at previous + minimum_image_displacement
    """
    n_atoms = len(positions)
    n_residues = resindices.max() + 1

    # Build mapping: for each residue, list of atom indices in order
    # We need atoms in their original order within each residue
    residue_atom_lists: Dict[int, List[int]] = {i: [] for i in range(n_residues)}
    for atom_idx in range(n_atoms):
        residue_atom_lists[resindices[atom_idx]].append(atom_idx)

    # Group residues by size for batched processing
    size_to_residues: Dict[int, List[int]] = {}
    for res_idx, atom_list in residue_atom_lists.items():
        size = len(atom_list)
        if size not in size_to_residues:
            size_to_residues[size] = []
        size_to_residues[size].append(res_idx)

    # Process each size group
    for size, residue_list in size_to_residues.items():
        if size <= 1:
            # Single-atom residues (ions) - nothing to chain
            continue

        n_res = len(residue_list)

        # Gather atom indices for all residues of this size
        # Shape: (n_residues_of_size, atoms_per_residue)
        atom_indices = np.array(
            [residue_atom_lists[res_idx] for res_idx in residue_list], dtype=np.int64
        )

        # Extract positions for these residues
        # Shape: (n_res, size, 3)
        res_positions = positions[atom_indices]

        # Chain atoms: each atom placed at MIC distance from previous
        for i in range(1, size):
            # Displacement from atom i-1 to atom i
            displacement = res_positions[:, i, :] - res_positions[:, i - 1, :]

            # Convert to fractional, apply MIC, convert back
            disp_frac = displacement @ box_matrix_inv.T
            disp_mic_frac = disp_frac - np.round(disp_frac)
            disp_mic = disp_mic_frac @ box_matrix.T

            # Place atom i at (atom i-1) + minimum_image_displacement
            res_positions[:, i, :] = res_positions[:, i - 1, :] + disp_mic

        # Write back to positions array
        positions[atom_indices.ravel()] = res_positions.reshape(-1, 3)

    LOGGER.debug(f"Made {n_residues} molecules whole")


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

    Notes
    -----
    Uses np.bincount for O(N_atoms) complexity.
    """
    if n_residues is None:
        n_residues = resindices.max() + 1

    if masses is None:
        masses = np.ones(len(positions))

    # Weighted position sums
    weighted_pos = positions * masses[:, np.newaxis]

    com_x = np.bincount(resindices, weights=weighted_pos[:, 0], minlength=n_residues)
    com_y = np.bincount(resindices, weights=weighted_pos[:, 1], minlength=n_residues)
    com_z = np.bincount(resindices, weights=weighted_pos[:, 2], minlength=n_residues)

    total_masses = np.bincount(resindices, weights=masses, minlength=n_residues)
    total_masses = np.where(total_masses == 0, 1.0, total_masses)

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

    Each residue (molecule) is translated as a unit so its center of mass is at
    the minimum image distance from the reference point.

    Parameters
    ----------
    positions : np.ndarray
        Atomic positions, shape (N_atoms, 3). **Modified in-place.**
    resindices : np.ndarray
        Residue index for each atom, shape (N_atoms,).
    reference_point : np.ndarray
        The point to cluster around (e.g., protein COM), shape (3,).
    box_matrix : np.ndarray
        3×3 box matrix.
    masses : np.ndarray, optional
        Atomic masses for COM calculation.
    exclude_resindices : np.ndarray, optional
        Residue indices to skip (e.g., protein residues).

    Returns
    -------
    np.ndarray
        Modified positions array (same object as input).
    """
    n_residues = resindices.max() + 1
    box_matrix_inv = np.linalg.inv(box_matrix)

    # Compute COM for each residue
    residue_coms, _ = compute_residue_coms(positions, resindices, masses, n_residues)

    # Displacement from reference
    displacements = residue_coms - reference_point

    # Minimum image shift
    shifts = minimum_image_shift(displacements, box_matrix, box_matrix_inv)

    # Zero out excluded residues
    if exclude_resindices is not None and len(exclude_resindices) > 0:
        shifts[exclude_resindices] = 0.0

    # Apply to atoms
    positions += shifts[resindices]

    return positions


class ClusterAroundProtein:
    """
    MDAnalysis transformation that creates a "droplet" visualization.

    This transformation applies two phases to create a visualization where
    all molecules are whole and clustered around the protein:

    1. **NoJump (all frames)**: Prevent atoms from jumping more than half a box
       between frames, maintaining molecular continuity from the initial state.

    2. **Cluster**: Translate each molecule to be at minimum image distance
       from the protein center of mass.

    The result is a "droplet" visualization where:
    - Protein and ligand are centered
    - All molecules are whole (not split across boundaries)
    - Solvent clusters around the protein
    - Some atoms may be outside the primary unit cell (that's intentional!)

    .. important::

       This transformation requires trajectories written with OpenMM's
       ``enforcePeriodicBox=False`` option. PolyzyMD's simulation runner
       writes trajectories this way by default. Trajectories written with
       per-frame PBC wrapping will NOT work correctly.

    This transformation is for **visualization only**, not for diffusion analysis.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The Universe containing the trajectory.
    protein_selection : str
        MDAnalysis selection string for protein atoms (default: "protein").
    ligand_selection : str, optional
        MDAnalysis selection string for ligand atoms (e.g., "resname LIG").

    Examples
    --------
    >>> import MDAnalysis as mda
    >>> from polyzymd.analysis.unwrap import ClusterAroundProtein
    >>>
    >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
    >>> transform = ClusterAroundProtein(u, "protein", "resname SUB")
    >>> u.trajectory.add_transformations(transform)
    >>>
    >>> for ts in u.trajectory:
    ...     # Positions show protein-centered droplet with whole molecules
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

        # Get residue indices for reference group
        self.reference_resindices = np.unique(self.reference_ag.resindices)

        # Check if masses are available
        try:
            _ = self.all_atoms.masses
            self._has_masses = True
        except Exception:
            self._has_masses = False
            LOGGER.warning("Masses not available in topology. Using geometric centers.")

        # State for NoJump algorithm
        self._prev_frac: Optional[np.ndarray] = None  # Previous frame fractional coords
        self._is_first_frame = True

        # Cache for box matrix
        self._cached_box_matrix: Optional[np.ndarray] = None
        self._cached_box_matrix_inv: Optional[np.ndarray] = None
        self._cached_dimensions: Optional[np.ndarray] = None

        LOGGER.info(
            f"ClusterAroundProtein initialized: "
            f"{len(self.protein_ag)} protein atoms, "
            f"{len(self.ligand_ag) if self.ligand_ag else 0} ligand atoms, "
            f"{len(self.reference_resindices)} reference residues, "
            f"{len(self.all_atoms)} total atoms"
        )

    def _get_box_matrices(self, dimensions: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Get box matrix and inverse, using cache if dimensions unchanged."""
        if self._cached_dimensions is None or not np.allclose(dimensions, self._cached_dimensions):
            self._cached_box_matrix = box_dimensions_to_matrix(dimensions)
            self._cached_box_matrix_inv = np.linalg.inv(self._cached_box_matrix)
            self._cached_dimensions = dimensions.copy()
        return self._cached_box_matrix, self._cached_box_matrix_inv

    def __call__(self, ts):
        """
        Apply the NoJump + Cluster transformation to the current timestep.

        Phase 1: NoJump - prevent atoms from jumping >half box between frames
        Phase 2: Cluster around protein COM

        Parameters
        ----------
        ts : MDAnalysis.coordinates.timestep.Timestep
            The current timestep.

        Returns
        -------
        ts : MDAnalysis.coordinates.timestep.Timestep
            The modified timestep.

        Notes
        -----
        This transformation assumes trajectories were written with
        ``enforcePeriodicBox=False``, so molecules are already whole.
        """
        dimensions = ts.dimensions
        if dimensions is None:
            raise ValueError("Trajectory has no box dimensions.")

        box_matrix, box_matrix_inv = self._get_box_matrices(dimensions)

        # Get positions and metadata
        positions = self.all_atoms.positions
        resindices = self.all_atoms.resindices
        masses = self.all_atoms.masses if self._has_masses else None

        # =====================================================================
        # PHASE 1: NoJump - maintain molecular continuity between frames
        # =====================================================================
        # Trajectories are written with enforcePeriodicBox=False, so molecules
        # start whole. NoJump ensures they stay whole by preventing any atom
        # from jumping more than half a box length between frames.

        if self._is_first_frame:
            # Frame 0: Just store fractional coordinates for NoJump on next frame
            # Molecules are already whole from PACKMOL/initial configuration
            LOGGER.debug("Frame 0: Storing initial fractional coordinates")
            self._prev_frac = positions @ box_matrix_inv.T
            self._is_first_frame = False
        else:
            # Frame 1+: Apply NoJump algorithm (Kulke & Vermaas 2022)
            # This prevents atoms from jumping more than half a box between frames
            pos_frac = positions @ box_matrix_inv.T
            pos_frac_unwrapped = pos_frac - np.round(pos_frac - self._prev_frac)
            positions[:] = pos_frac_unwrapped @ box_matrix.T

            # Update stored fractional coords for next frame
            self._prev_frac = pos_frac_unwrapped

        # =====================================================================
        # PHASE 2: Cluster around protein
        # =====================================================================

        # Compute protein+ligand center of mass
        ref_positions = positions[self.reference_ag.ix]
        if self._has_masses:
            ref_masses = self.reference_ag.masses
            ref_com = np.average(ref_positions, axis=0, weights=ref_masses)
        else:
            ref_com = np.mean(ref_positions, axis=0)

        # Cluster all molecules around the protein COM
        # Molecules will be translated to their minimum image position relative to protein
        cluster_around_point(
            positions=positions,
            resindices=resindices,
            reference_point=ref_com,
            box_matrix=box_matrix,
            masses=masses,
            exclude_resindices=self.reference_resindices,
        )

        # Update positions in the Universe
        self.all_atoms.positions = positions

        return ts

    def reset(self) -> None:
        """
        Reset the transformation state.

        Call this if you need to reprocess the trajectory from the beginning,
        or if you're jumping to a non-sequential frame.
        """
        self._prev_frac = None
        self._is_first_frame = True
        LOGGER.debug("ClusterAroundProtein state reset")
