"""Protein residue selectors.

This module provides selectors for protein residues:

- ProteinResidues: Select all protein residues
- ProteinResiduesByGroup: Select protein residues by amino acid classification
- ProteinResiduesNearReference: Select residues within cutoff of reference atoms

Examples
--------
>>> # Select all protein residues
>>> selector = ProteinResidues()
>>> result = selector.select(universe)
>>>
>>> # Select aromatic residues only
>>> selector = ProteinResiduesByGroup(
...     grouping=ProteinAAClassification(),
...     groups=["aromatic"]
... )
>>>
>>> # Select residues near catalytic triad
>>> selector = ProteinResiduesNearReference(
...     reference_selection="resid 77 133 156",
...     cutoff=5.0
... )
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from polyzymd.analysis.common.selectors.base import MolecularSelector, SelectionResult

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from polyzymd.analysis.common.groupings.base import ResidueGrouping


class ProteinResidues(MolecularSelector):
    """Select all protein residues.

    Uses MDAnalysis "protein" selection keyword which matches standard
    amino acid residues.

    Parameters
    ----------
    selection_modifier : str, optional
        Additional selection criteria to AND with "protein".
        E.g., "and not name H*" to exclude hydrogens.
    """

    def __init__(self, selection_modifier: str | None = None):
        self.selection_modifier = selection_modifier

    def select(self, universe: "Universe") -> SelectionResult:
        """Select all protein atoms/residues."""
        selection = "protein"
        if self.selection_modifier:
            selection = f"({selection}) and ({self.selection_modifier})"

        atoms = universe.select_atoms(selection)

        if len(atoms) == 0:
            raise ValueError(
                f"Selection '{selection}' matched no atoms. "
                "Check that the topology contains protein residues."
            )

        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={"selection": selection},
        )

    @property
    def label(self) -> str:
        return "protein"


class ProteinResiduesByGroup(MolecularSelector):
    """Select protein residues by amino acid group classification.

    Uses a ResidueGrouping to classify amino acids (e.g., aromatic, charged,
    polar, nonpolar) and selects only residues in the specified groups.

    Parameters
    ----------
    grouping : ResidueGrouping
        Classification scheme for amino acids
    groups : list[str]
        Names of groups to include (e.g., ["aromatic", "charged_positive"])
    exclude : bool, optional
        If True, select residues NOT in the specified groups. Default False.

    Examples
    --------
    >>> from polyzymd.analysis.common.groupings import ProteinAAClassification
    >>>
    >>> # Select aromatic residues
    >>> grouping = ProteinAAClassification()
    >>> selector = ProteinResiduesByGroup(grouping, groups=["aromatic"])
    >>>
    >>> # Select all charged residues
    >>> selector = ProteinResiduesByGroup(
    ...     grouping,
    ...     groups=["charged_positive", "charged_negative"]
    ... )
    """

    def __init__(
        self,
        grouping: "ResidueGrouping",
        groups: list[str],
        exclude: bool = False,
    ):
        self.grouping = grouping
        self.groups = groups
        self.exclude = exclude

    def select(self, universe: "Universe") -> SelectionResult:
        """Select protein residues matching the specified groups."""
        # First get all protein residues
        protein_atoms = universe.select_atoms("protein")
        if len(protein_atoms) == 0:
            raise ValueError("No protein atoms found in universe")

        protein_residues = protein_atoms.residues

        # Filter by group membership
        matching_resids = []
        for res in protein_residues:
            group = self.grouping.classify(res.resname)
            in_group = group in self.groups

            if self.exclude:
                if not in_group:
                    matching_resids.append(res.resid)
            else:
                if in_group:
                    matching_resids.append(res.resid)

        if not matching_resids:
            group_str = ", ".join(self.groups)
            mode = "excluding" if self.exclude else "in"
            raise ValueError(
                f"No protein residues found {mode} groups: {group_str}. "
                f"Available groups: {self.grouping.available_groups}"
            )

        # Select the matching residues
        resid_str = " ".join(str(r) for r in matching_resids)
        atoms = universe.select_atoms(f"protein and resid {resid_str}")

        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={
                "groups": self.groups,
                "exclude": self.exclude,
                "n_matching": len(matching_resids),
            },
        )

    @property
    def label(self) -> str:
        prefix = "not_" if self.exclude else ""
        return f"protein_{prefix}{'_'.join(self.groups)}"


class ProteinResiduesNearReference(MolecularSelector):
    """Select protein residues within a cutoff distance of reference atoms.

    Useful for selecting residues near active sites, binding pockets,
    or other regions of interest.

    Parameters
    ----------
    reference_selection : str
        MDAnalysis selection string for reference atoms (e.g., "resid 77 133 156")
    cutoff : float
        Distance cutoff in Angstroms. Residues with any atom within this
        distance of any reference atom are selected.
    include_reference : bool, optional
        Whether to include the reference residues themselves. Default True.
    frame : int, optional
        Frame to use for distance calculation. Default is current frame (0).

    Examples
    --------
    >>> # Select residues within 5A of catalytic triad
    >>> selector = ProteinResiduesNearReference(
    ...     reference_selection="resid 77 133 156",
    ...     cutoff=5.0,
    ... )
    >>>
    >>> # Select residues near substrate binding site (not including the site itself)
    >>> selector = ProteinResiduesNearReference(
    ...     reference_selection="resname LIG",
    ...     cutoff=4.0,
    ...     include_reference=False,
    ... )
    """

    def __init__(
        self,
        reference_selection: str,
        cutoff: float,
        include_reference: bool = True,
        frame: int = 0,
    ):
        self.reference_selection = reference_selection
        self.cutoff = cutoff
        self.include_reference = include_reference
        self.frame = frame

    def select(self, universe: "Universe") -> SelectionResult:
        """Select protein residues near the reference atoms."""
        # Go to the specified frame
        universe.trajectory[self.frame]

        # Select reference atoms
        ref_atoms = universe.select_atoms(self.reference_selection)
        if len(ref_atoms) == 0:
            raise ValueError(f"Reference selection '{self.reference_selection}' matched no atoms")

        # Select protein atoms within cutoff of reference
        # MDAnalysis "around" selection finds atoms within cutoff
        nearby_selection = f"protein and around {self.cutoff} ({self.reference_selection})"
        nearby_atoms = universe.select_atoms(nearby_selection)

        if self.include_reference:
            # Include reference residues if they are protein
            ref_protein = universe.select_atoms(
                f"protein and same residue as ({self.reference_selection})"
            )
            nearby_atoms = nearby_atoms | ref_protein

        if len(nearby_atoms) == 0:
            raise ValueError(
                f"No protein atoms found within {self.cutoff}A of '{self.reference_selection}'"
            )

        return SelectionResult(
            atoms=nearby_atoms,
            residues=nearby_atoms.residues,
            label=self.label,
            metadata={
                "reference_selection": self.reference_selection,
                "cutoff": self.cutoff,
                "include_reference": self.include_reference,
                "frame": self.frame,
                "n_reference_atoms": len(ref_atoms),
            },
        )

    @property
    def label(self) -> str:
        return f"near_ref_{self.cutoff:.1f}A"
