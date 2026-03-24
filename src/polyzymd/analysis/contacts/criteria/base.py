"""Base class for contact criteria.

Defines the abstract interface for determining whether two
molecular groups are "in contact".
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup


@dataclass
class ContactCheckResult:
    """Result of checking if two groups are in contact.

    Attributes
    ----------
    is_contact : bool
        Whether the groups are in contact
    distance : float
        The distance used for the contact check (interpretation depends on criteria)
    metadata : dict
        Additional information (e.g., which atoms formed the contact)
    """

    is_contact: bool
    distance: float
    metadata: dict | None = None


class ContactCriteria(ABC):
    """Abstract base class for contact criteria.

    Defines what constitutes a "contact" between two molecular groups.
    Different criteria can be used depending on the analysis needs:

    - Distance-based (any atom, COM, minimum distance)
    - Hydrogen bonding
    - Hydrophobic contacts
    - Custom criteria

    The Strategy pattern allows swapping criteria without changing
    the analysis code.

    Examples
    --------
    >>> class HBondCriteria(ContactCriteria):
    ...     '''Hydrogen bond contact criteria.'''
    ...
    ...     def __init__(self, distance_cutoff=3.5, angle_cutoff=30.0):
    ...         self.distance_cutoff = distance_cutoff
    ...         self.angle_cutoff = angle_cutoff
    ...
    ...     def check_contact(self, query_atoms, target_atoms, box=None):
    ...         # Custom H-bond detection logic
    ...         ...
    ...
    ...     @property
    ...     def label(self):
    ...         return f"hbond_{self.distance_cutoff}A_{self.angle_cutoff}deg"
    """

    @abstractmethod
    def check_contact(
        self,
        query_atoms: "AtomGroup",
        target_atoms: "AtomGroup",
        box: NDArray[np.float64] | None = None,
    ) -> ContactCheckResult:
        """Check if two atom groups are in contact.

        Parameters
        ----------
        query_atoms : AtomGroup
            The "query" group (e.g., polymer residue atoms)
        target_atoms : AtomGroup
            The "target" group (e.g., protein residue atoms)
        box : NDArray[np.float64], optional
            Box dimensions for periodic boundary conditions.
            Shape (3,) for orthorhombic or (3,3) for triclinic.

        Returns
        -------
        ContactCheckResult
            Contains is_contact bool, distance, and optional metadata
        """
        ...

    @property
    @abstractmethod
    def label(self) -> str:
        """Short label for this criteria (for filenames/logging)."""
        ...

    @property
    @abstractmethod
    def cutoff(self) -> float:
        """The cutoff distance used by this criteria."""
        ...

    def to_dict(self) -> dict[str, Any]:
        """Serialize criteria parameters to dictionary."""
        return {
            "type": self.__class__.__name__,
            "label": self.label,
            "cutoff": self.cutoff,
        }

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(cutoff={self.cutoff})"
