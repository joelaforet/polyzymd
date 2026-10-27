"""Residue groupings: :class:`ProteinAAClassification` and its interface :class:`ResidueGrouping`.

Examples
--------
>>> from polyzymd.analyses.shared.groupings import ProteinAAClassification
>>>
>>> # Classify amino acids
>>> grouping = ProteinAAClassification()
>>> print(grouping.classify("PHE"))  # "aromatic"
>>> print(grouping.classify("LYS"))  # "charged_positive"
>>>
>>> # Get all residues in a group
>>> aromatics = grouping.get_residues_in_group("aromatic")
>>> # Returns: ["PHE", "TRP", "TYR", "HIS"]
"""

from polyzymd.analyses.shared.groupings.base import ProteinAAClassification, ResidueGrouping

__all__ = [
    "ResidueGrouping",
    "ProteinAAClassification",
]
