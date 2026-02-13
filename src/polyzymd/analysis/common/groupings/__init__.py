"""Residue grouping abstractions for analysis modules.

This module provides classification systems for residues:

- ResidueGrouping: Abstract base class for residue classification
- ProteinAAClassification: Standard amino acid classification
- CustomGrouping: User-defined classification scheme

Examples
--------
>>> from polyzymd.analysis.common.groupings import ProteinAAClassification
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

from polyzymd.analysis.common.groupings.base import (
    ResidueGrouping,
    ProteinAAClassification,
    CustomGrouping,
)

__all__ = [
    "ResidueGrouping",
    "ProteinAAClassification",
    "CustomGrouping",
]
