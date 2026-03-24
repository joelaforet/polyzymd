"""Distance-based contact criteria.

This module provides the most common contact criteria based on
distance measurements:

- AnyAtomWithinCutoff: Contact if ANY atom pair is within cutoff
- AnyAtomToCOM: Contact if any query atom is within cutoff of target COM
- COMToCOM: Contact if COMs are within cutoff
- MinimumDistance: Returns minimum distance for continuous analysis

These cover the most common use cases. For specialized criteria
(e.g., hydrogen bonds), subclass ContactCriteria.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from polyzymd.analysis.contacts.criteria.base import ContactCheckResult, ContactCriteria
from polyzymd.analysis.core.pbc import minimum_image_distance, pairwise_distances_pbc

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup


class AnyAtomWithinCutoff(ContactCriteria):
    """Contact if any atom pair is within the cutoff distance.

    This is the default and most common criteria. A contact exists if
    ANY atom from the query group is within the cutoff of ANY atom
    in the target group.

    Parameters
    ----------
    cutoff : float
        Distance cutoff in Angstroms. Default 4.0.
    use_pbc : bool
        Whether to apply periodic boundary conditions. Default True.

    Examples
    --------
    >>> criteria = AnyAtomWithinCutoff(cutoff=4.0)
    >>> result = criteria.check_contact(polymer_atoms, protein_atoms, box)
    >>> if result.is_contact:
    ...     print(f"Contact! Minimum distance: {result.distance:.2f} A")
    """

    def __init__(self, cutoff: float = 4.0, use_pbc: bool = True):
        self._cutoff = cutoff
        self.use_pbc = use_pbc

    def check_contact(
        self,
        query_atoms: "AtomGroup",
        target_atoms: "AtomGroup",
        box: NDArray[np.float64] | None = None,
    ) -> ContactCheckResult:
        """Check if any atom pair is within cutoff."""
        if len(query_atoms) == 0 or len(target_atoms) == 0:
            return ContactCheckResult(
                is_contact=False,
                distance=float("inf"),
                metadata={"reason": "empty_selection"},
            )

        # Get positions
        pos1 = query_atoms.positions
        pos2 = target_atoms.positions

        # Apply PBC if requested
        pbc_box = box if self.use_pbc else None

        # Calculate pairwise distances
        distances = pairwise_distances_pbc(pos1, pos2, pbc_box)

        # Find minimum distance
        min_dist = float(np.min(distances))

        # Find which atoms form the closest contact
        min_idx = np.unravel_index(np.argmin(distances), distances.shape)

        return ContactCheckResult(
            is_contact=min_dist <= self._cutoff,
            distance=min_dist,
            metadata={
                "query_atom_idx": int(min_idx[0]),
                "target_atom_idx": int(min_idx[1]),
            },
        )

    @property
    def label(self) -> str:
        return f"any_atom_{self._cutoff:.1f}A"

    @property
    def cutoff(self) -> float:
        return self._cutoff


class AnyAtomToCOM(ContactCriteria):
    """Contact if any query atom is within cutoff of target center of mass.

    This is the criteria described in the requirements: "any atom of
    polymer residue within cutoff of COM of protein residue."

    Useful when you want to detect if a polymer residue is "near"
    a protein residue's center, rather than near any specific atom.

    Parameters
    ----------
    cutoff : float
        Distance cutoff in Angstroms. Default 4.0.
    use_pbc : bool
        Whether to apply periodic boundary conditions. Default True.
    use_center_of_geometry : bool
        If True, use center of geometry instead of center of mass.
        Default False (use COM).

    Examples
    --------
    >>> # Contact if any polymer atom is within 4A of protein residue COM
    >>> criteria = AnyAtomToCOM(cutoff=4.0)
    """

    def __init__(
        self,
        cutoff: float = 4.0,
        use_pbc: bool = True,
        use_center_of_geometry: bool = False,
    ):
        self._cutoff = cutoff
        self.use_pbc = use_pbc
        self.use_center_of_geometry = use_center_of_geometry

    def check_contact(
        self,
        query_atoms: "AtomGroup",
        target_atoms: "AtomGroup",
        box: NDArray[np.float64] | None = None,
    ) -> ContactCheckResult:
        """Check if any query atom is within cutoff of target COM."""
        if len(query_atoms) == 0 or len(target_atoms) == 0:
            return ContactCheckResult(
                is_contact=False,
                distance=float("inf"),
                metadata={"reason": "empty_selection"},
            )

        # Get target center
        if self.use_center_of_geometry:
            target_center = target_atoms.center_of_geometry()
        else:
            target_center = target_atoms.center_of_mass()

        # Get query positions
        query_positions = query_atoms.positions

        # Calculate distances to center
        pbc_box = box if self.use_pbc else None

        min_dist = float("inf")
        closest_idx = 0

        for i, pos in enumerate(query_positions):
            dist = minimum_image_distance(pos, target_center, pbc_box)
            if dist < min_dist:
                min_dist = dist
                closest_idx = i

        return ContactCheckResult(
            is_contact=min_dist <= self._cutoff,
            distance=min_dist,
            metadata={
                "closest_query_atom_idx": closest_idx,
                "target_center": target_center.tolist(),
            },
        )

    @property
    def label(self) -> str:
        center_type = "cog" if self.use_center_of_geometry else "com"
        return f"atom_to_{center_type}_{self._cutoff:.1f}A"

    @property
    def cutoff(self) -> float:
        return self._cutoff


class COMToCOM(ContactCriteria):
    """Contact if centers of mass are within cutoff.

    Simple criteria based on the distance between the COMs of
    the two groups. Useful for coarse-grained analysis.

    Parameters
    ----------
    cutoff : float
        Distance cutoff in Angstroms. Default 8.0 (larger since COM-to-COM).
    use_pbc : bool
        Whether to apply periodic boundary conditions. Default True.
    use_center_of_geometry : bool
        If True, use center of geometry instead of center of mass.
        Default False.
    """

    def __init__(
        self,
        cutoff: float = 8.0,
        use_pbc: bool = True,
        use_center_of_geometry: bool = False,
    ):
        self._cutoff = cutoff
        self.use_pbc = use_pbc
        self.use_center_of_geometry = use_center_of_geometry

    def check_contact(
        self,
        query_atoms: "AtomGroup",
        target_atoms: "AtomGroup",
        box: NDArray[np.float64] | None = None,
    ) -> ContactCheckResult:
        """Check if COMs are within cutoff."""
        if len(query_atoms) == 0 or len(target_atoms) == 0:
            return ContactCheckResult(
                is_contact=False,
                distance=float("inf"),
                metadata={"reason": "empty_selection"},
            )

        # Get centers
        if self.use_center_of_geometry:
            query_center = query_atoms.center_of_geometry()
            target_center = target_atoms.center_of_geometry()
        else:
            query_center = query_atoms.center_of_mass()
            target_center = target_atoms.center_of_mass()

        # Calculate distance
        pbc_box = box if self.use_pbc else None
        dist = minimum_image_distance(query_center, target_center, pbc_box)

        return ContactCheckResult(
            is_contact=dist <= self._cutoff,
            distance=dist,
            metadata={
                "query_center": query_center.tolist(),
                "target_center": target_center.tolist(),
            },
        )

    @property
    def label(self) -> str:
        center_type = "cog" if self.use_center_of_geometry else "com"
        return f"{center_type}_to_{center_type}_{self._cutoff:.1f}A"

    @property
    def cutoff(self) -> float:
        return self._cutoff


class MinimumDistance(ContactCriteria):
    """Returns minimum distance between groups (for continuous analysis).

    Unlike other criteria that return binary contact yes/no, this
    is useful when you want to track the actual distance over time
    and apply cutoffs later.

    The is_contact field is always True if groups are non-empty
    (since we're interested in the distance value itself).

    Parameters
    ----------
    use_pbc : bool
        Whether to apply periodic boundary conditions. Default True.

    Examples
    --------
    >>> criteria = MinimumDistance()
    >>> for ts in universe.trajectory:
    ...     result = criteria.check_contact(polymer, protein, ts.dimensions)
    ...     distances.append(result.distance)
    """

    def __init__(self, use_pbc: bool = True):
        self.use_pbc = use_pbc
        self._cutoff = float("inf")  # No cutoff - always report distance

    def check_contact(
        self,
        query_atoms: "AtomGroup",
        target_atoms: "AtomGroup",
        box: NDArray[np.float64] | None = None,
    ) -> ContactCheckResult:
        """Calculate minimum distance between groups."""
        if len(query_atoms) == 0 or len(target_atoms) == 0:
            return ContactCheckResult(
                is_contact=False,
                distance=float("inf"),
                metadata={"reason": "empty_selection"},
            )

        pos1 = query_atoms.positions
        pos2 = target_atoms.positions
        pbc_box = box if self.use_pbc else None

        distances = pairwise_distances_pbc(pos1, pos2, pbc_box)
        min_dist = float(np.min(distances))
        min_idx = np.unravel_index(np.argmin(distances), distances.shape)

        return ContactCheckResult(
            is_contact=True,  # Always "in contact" for distance tracking
            distance=min_dist,
            metadata={
                "query_atom_idx": int(min_idx[0]),
                "target_atom_idx": int(min_idx[1]),
                "query_atom_name": query_atoms[min_idx[0]].name,
                "target_atom_name": target_atoms[min_idx[1]].name,
            },
        )

    @property
    def label(self) -> str:
        return "min_distance"

    @property
    def cutoff(self) -> float:
        return self._cutoff
