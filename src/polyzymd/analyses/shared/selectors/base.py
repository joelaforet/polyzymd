"""Base class for molecular selectors.

This module defines the abstract base class for all molecular selectors,
providing a consistent interface for selecting atoms, residues, or groups
from an MDAnalysis Universe.

The Strategy pattern allows users to define custom selection logic by
subclassing MolecularSelector and implementing the select() method.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup, ResidueGroup
    from MDAnalysis.core.universe import Universe


@dataclass
class SelectionResult:
    """Container for selection results with metadata.

    Attributes
    ----------
    atoms : AtomGroup
        The selected atoms
    residues : ResidueGroup
        The residues containing the selected atoms
    label : str
        Human-readable label for this selection
    metadata : dict
        Additional metadata about the selection (e.g., selection string used,
        cutoff values, etc.)
    """

    atoms: "AtomGroup"
    residues: "ResidueGroup"
    label: str
    metadata: dict = field(default_factory=dict)

    @property
    def n_atoms(self) -> int:
        """Number of selected atoms."""
        return len(self.atoms)

    @property
    def n_residues(self) -> int:
        """Number of selected residues."""
        return len(self.residues)

    @property
    def residue_ids(self) -> NDArray[np.int64]:
        """1-indexed residue IDs (PyMOL convention)."""
        return np.array(self.residues.resids, dtype=np.int64)

    @property
    def residue_names(self) -> list[str]:
        """Residue names for each residue."""
        return list(self.residues.resnames)

    def __repr__(self) -> str:
        return (
            f"SelectionResult(label='{self.label}', "
            f"n_atoms={self.n_atoms}, n_residues={self.n_residues})"
        )


class MolecularSelector(ABC):
    """Abstract base class for molecular selections.

    Subclasses must implement the select() method to define how atoms
    or residues are selected from a Universe.

    This follows the Strategy pattern - different selectors can be
    swapped in to change selection behavior without modifying the
    analysis code.

    Examples
    --------
    >>> class ActiveSiteSelector(MolecularSelector):
    ...     def __init__(self, active_site_resids: list[int]):
    ...         self.resids = active_site_resids
    ...
    ...     def select(self, universe: Universe) -> SelectionResult:
    ...         resid_str = " ".join(str(r) for r in self.resids)
    ...         atoms = universe.select_atoms(f"resid {resid_str}")
    ...         return SelectionResult(
    ...             atoms=atoms,
    ...             residues=atoms.residues,
    ...             label="active_site",
    ...             metadata={"resids": self.resids}
    ...         )
    ...
    ...     @property
    ...     def label(self) -> str:
    ...         return "active_site"
    """

    @abstractmethod
    def select(self, universe: "Universe") -> SelectionResult:
        """Select atoms/residues from a Universe.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe to select from

        Returns
        -------
        SelectionResult
            Container with selected atoms, residues, and metadata
        """
        ...

    @property
    @abstractmethod
    def label(self) -> str:
        """Short label identifying this selector (for filenames/logging)."""
        ...

    def validate(self, universe: "Universe") -> dict[str, Any]:
        """Validate the selector against a Universe.

        Returns diagnostic information about whether the selection
        would succeed and what it would select.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe to validate against

        Returns
        -------
        dict
            Validation results with keys:
            - valid: bool
            - n_atoms: int
            - n_residues: int
            - error: str (if invalid)
            - warnings: list[str]
        """
        try:
            result = self.select(universe)
            warnings = []

            if result.n_atoms == 0:
                return {
                    "valid": False,
                    "error": "Selection matched no atoms",
                    "n_atoms": 0,
                    "n_residues": 0,
                    "warnings": [],
                }

            return {
                "valid": True,
                "n_atoms": result.n_atoms,
                "n_residues": result.n_residues,
                "warnings": warnings,
            }

        except Exception as e:
            return {
                "valid": False,
                "error": str(e),
                "n_atoms": 0,
                "n_residues": 0,
                "warnings": [],
            }

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(label='{self.label}')"


class MDAnalysisSelector(MolecularSelector):
    """Simple selector using an MDAnalysis selection string.

    This is the most flexible selector - it allows arbitrary MDAnalysis
    selection syntax. Use this when you need direct control over the
    selection or when the specialized selectors don't fit your needs.

    Parameters
    ----------
    selection : str
        MDAnalysis selection string (e.g., "protein", "resname SBM EGM",
        "resid 1-50 and name CA")
    label : str, optional
        Human-readable label. If not provided, uses a sanitized version
        of the selection string.

    Examples
    --------
    >>> # Select polymer residues by name
    >>> selector = MDAnalysisSelector("resname SBM EGM")
    >>> result = selector.select(universe)
    >>>
    >>> # Select protein backbone near ligand
    >>> selector = MDAnalysisSelector(
    ...     "protein and backbone and around 5.0 resname LIG",
    ...     label="protein_near_ligand"
    ... )
    """

    def __init__(self, selection: str, label: str | None = None):
        self.selection = selection
        self._label = label

    def select(self, universe: "Universe") -> SelectionResult:
        """Select atoms using the MDAnalysis selection string."""
        atoms = universe.select_atoms(self.selection)
        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={"selection": self.selection},
        )

    @property
    def label(self) -> str:
        if self._label:
            return self._label
        # Create a sanitized label from selection string
        sanitized = self.selection.replace(" ", "_").replace(".", "_")
        # Truncate if too long
        if len(sanitized) > 30:
            sanitized = sanitized[:27] + "..."
        return sanitized


class CompositeSelector(MolecularSelector):
    """Combines multiple selectors with AND/OR logic.

    Useful for complex selections like "protein residues that are
    both aromatic AND within 5A of the active site".

    Parameters
    ----------
    selectors : list[MolecularSelector]
        List of selectors to combine
    mode : {"union", "intersection"}
        How to combine selections:
        - "union": Include atoms selected by ANY selector (OR)
        - "intersection": Include only atoms selected by ALL selectors (AND)
    label : str, optional
        Custom label. If not provided, generates from component labels.
    """

    def __init__(
        self,
        selectors: list[MolecularSelector],
        mode: str = "union",
        label: str | None = None,
    ):
        if mode not in ("union", "intersection"):
            raise ValueError(f"mode must be 'union' or 'intersection', got {mode}")

        self.selectors = selectors
        self.mode = mode
        self._label = label

    def select(self, universe: "Universe") -> SelectionResult:
        """Select atoms using combined selectors."""
        if not self.selectors:
            raise ValueError("CompositeSelector requires at least one selector")

        results = [s.select(universe) for s in self.selectors]

        if self.mode == "union":
            # OR: combine all atoms
            combined_atoms = results[0].atoms
            for r in results[1:]:
                combined_atoms = combined_atoms | r.atoms
        else:
            # AND: intersect all atoms
            combined_atoms = results[0].atoms
            for r in results[1:]:
                combined_atoms = combined_atoms & r.atoms

        return SelectionResult(
            atoms=combined_atoms,
            residues=combined_atoms.residues,
            label=self.label,
            metadata={
                "mode": self.mode,
                "component_labels": [s.label for s in self.selectors],
            },
        )

    @property
    def label(self) -> str:
        if self._label:
            return self._label
        op = "_or_" if self.mode == "union" else "_and_"
        return op.join(s.label for s in self.selectors)
