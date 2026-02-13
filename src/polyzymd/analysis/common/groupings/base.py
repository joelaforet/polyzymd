"""Base classes for residue grouping/classification.

This module provides the abstract base class for residue classification
schemes and concrete implementations for protein amino acids.

The Strategy pattern allows users to define custom classification schemes
for polymers, modified residues, or other systems.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class ResidueGrouping(ABC):
    """Abstract base class for residue classification schemes.

    Subclasses must implement classify() to map residue names to group labels.

    Examples
    --------
    >>> class MyPolymerGrouping(ResidueGrouping):
    ...     def classify(self, resname: str) -> str:
    ...         if resname in ["SBM", "SBMA"]:
    ...             return "zwitterionic"
    ...         elif resname in ["EGP", "EGMA"]:
    ...             return "hydrophilic"
    ...         return "unknown"
    ...
    ...     @property
    ...     def available_groups(self) -> list[str]:
    ...         return ["zwitterionic", "hydrophilic", "unknown"]
    """

    @abstractmethod
    def classify(self, resname: str) -> str:
        """Classify a residue name into a group.

        Parameters
        ----------
        resname : str
            Residue name (3-letter code for amino acids)

        Returns
        -------
        str
            Group label for this residue type
        """
        ...

    @property
    @abstractmethod
    def available_groups(self) -> list[str]:
        """List of all group labels in this classification scheme."""
        ...

    def get_residues_in_group(self, group: str) -> list[str]:
        """Get all residue names that belong to a group.

        Parameters
        ----------
        group : str
            Group label

        Returns
        -------
        list[str]
            Residue names in this group

        Raises
        ------
        ValueError
            If group is not in available_groups
        """
        if group not in self.available_groups:
            raise ValueError(f"Unknown group '{group}'. Available: {self.available_groups}")
        return [r for r in self._all_resnames if self.classify(r) == group]

    @property
    def _all_resnames(self) -> list[str]:
        """All residue names known to this grouping (override in subclass)."""
        return []

    def to_dict(self) -> dict[str, Any]:
        """Serialize grouping scheme to dictionary."""
        return {
            "type": self.__class__.__name__,
            "groups": {group: self.get_residues_in_group(group) for group in self.available_groups},
        }


class ProteinAAClassification(ResidueGrouping):
    """Standard amino acid classification.

    Groups amino acids into:
    - aromatic: PHE, TRP, TYR, HIS
    - charged_positive: ARG, LYS
    - charged_negative: ASP, GLU
    - polar: ASN, CYS, GLN, SER, THR
    - nonpolar: ALA, GLY, ILE, LEU, MET, PRO, VAL

    This classification matches the scaffold notebooks and common
    biochemistry conventions.

    Parameters
    ----------
    include_his_aromatic : bool, optional
        Whether to classify HIS as aromatic (default True).
        Some classifications put HIS with charged_positive.

    Examples
    --------
    >>> grouping = ProteinAAClassification()
    >>> grouping.classify("PHE")
    'aromatic'
    >>> grouping.classify("LYS")
    'charged_positive'
    >>> grouping.get_residues_in_group("aromatic")
    ['PHE', 'TRP', 'TYR', 'HIS']
    """

    # Standard amino acid classification
    # Based on scaffold notebooks: Contact_Analysis_SBMA_EGMA_EGPMA_per_Residue.ipynb
    _CLASSIFICATION = {
        # Aromatic
        "PHE": "aromatic",
        "TRP": "aromatic",
        "TYR": "aromatic",
        "HIS": "aromatic",  # Can be overridden
        # Charged positive
        "ARG": "charged_positive",
        "LYS": "charged_positive",
        # Charged negative
        "ASP": "charged_negative",
        "GLU": "charged_negative",
        # Polar
        "ASN": "polar",
        "CYS": "polar",
        "GLN": "polar",
        "SER": "polar",
        "THR": "polar",
        # Nonpolar
        "ALA": "nonpolar",
        "GLY": "nonpolar",
        "ILE": "nonpolar",
        "LEU": "nonpolar",
        "MET": "nonpolar",
        "PRO": "nonpolar",
        "VAL": "nonpolar",
    }

    _ALL_GROUPS = [
        "aromatic",
        "charged_positive",
        "charged_negative",
        "polar",
        "nonpolar",
        "unknown",
    ]

    def __init__(self, include_his_aromatic: bool = True):
        self.include_his_aromatic = include_his_aromatic
        self._classification = self._CLASSIFICATION.copy()

        if not include_his_aromatic:
            # Move HIS to charged_positive (alternative classification)
            self._classification["HIS"] = "charged_positive"

    def classify(self, resname: str) -> str:
        """Classify amino acid by residue name."""
        # Normalize: uppercase, handle common variants
        resname = resname.upper().strip()

        # Handle common protonation state variants
        variants = {
            "HIE": "HIS",
            "HID": "HIS",
            "HIP": "HIS",  # Protonated histidine
            "HSE": "HIS",
            "HSD": "HIS",
            "HSP": "HIS",
            "CYSH": "CYS",
            "CYX": "CYS",  # Disulfide bonded
            "ASH": "ASP",  # Protonated aspartate
            "GLH": "GLU",  # Protonated glutamate
            "LYN": "LYS",  # Neutral lysine
        }
        resname = variants.get(resname, resname)

        return self._classification.get(resname, "unknown")

    @property
    def available_groups(self) -> list[str]:
        return self._ALL_GROUPS

    @property
    def _all_resnames(self) -> list[str]:
        return list(self._CLASSIFICATION.keys())

    def get_charged_groups(self) -> list[str]:
        """Convenience: get both charged group names."""
        return ["charged_positive", "charged_negative"]

    def get_hydrophobic_groups(self) -> list[str]:
        """Convenience: groups typically considered hydrophobic."""
        return ["aromatic", "nonpolar"]

    def get_hydrophilic_groups(self) -> list[str]:
        """Convenience: groups typically considered hydrophilic."""
        return ["charged_positive", "charged_negative", "polar"]


class CustomGrouping(ResidueGrouping):
    """User-defined residue classification.

    Allows arbitrary mapping from residue names to group labels.

    Parameters
    ----------
    classification : dict[str, str]
        Mapping from residue name to group label.
    default_group : str, optional
        Group label for unclassified residues. Default "other".

    Examples
    --------
    >>> # Custom polymer classification
    >>> grouping = CustomGrouping({
    ...     "SBM": "zwitterionic",
    ...     "SBMA": "zwitterionic",
    ...     "EGP": "peg_like",
    ...     "EGMA": "peg_like",
    ... }, default_group="unknown")
    >>> grouping.classify("SBM")
    'zwitterionic'
    """

    def __init__(
        self,
        classification: dict[str, str],
        default_group: str = "other",
    ):
        self._classification = {k.upper(): v for k, v in classification.items()}
        self.default_group = default_group

        # Compute available groups
        self._groups = sorted(set(self._classification.values()))
        if default_group not in self._groups:
            self._groups.append(default_group)

    def classify(self, resname: str) -> str:
        """Classify residue by name using custom mapping."""
        return self._classification.get(resname.upper(), self.default_group)

    @property
    def available_groups(self) -> list[str]:
        return self._groups

    @property
    def _all_resnames(self) -> list[str]:
        return list(self._classification.keys())

    @classmethod
    def from_groups(
        cls, groups: dict[str, list[str]], default_group: str = "other"
    ) -> "CustomGrouping":
        """Create grouping from group -> residue list mapping.

        Parameters
        ----------
        groups : dict[str, list[str]]
            Mapping from group name to list of residue names
        default_group : str
            Group for unlisted residues

        Returns
        -------
        CustomGrouping

        Examples
        --------
        >>> grouping = CustomGrouping.from_groups({
        ...     "zwitterionic": ["SBM", "SBMA"],
        ...     "peg_like": ["EGP", "EGMA", "OEGMA"],
        ... })
        """
        classification = {}
        for group_name, resnames in groups.items():
            for resname in resnames:
                classification[resname] = group_name
        return cls(classification, default_group=default_group)
