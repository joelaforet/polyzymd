"""Common infrastructure for analysis modules.

This package provides shared components used across multiple analysis modules:

- **selectors**: Molecular selection abstractions (protein, polymer, solvent, etc.)
- **groupings**: Residue classification systems (amino acid types, custom groups)

These components enable consistent, extensible analysis across different
molecular systems while maintaining separation of concerns.
"""

from polyzymd.analysis.common.selectors import (
    MolecularSelector,
    ProteinResidues,
    ProteinResiduesByGroup,
    ProteinResiduesNearReference,
    PolymerChains,
    PolymerResiduesByType,
    SubstrateMolecule,
    SolventMolecules,
    CosolventMolecules,
)
from polyzymd.analysis.common.groupings import (
    ResidueGrouping,
    ProteinAAClassification,
    CustomGrouping,
)

__all__ = [
    # Selectors
    "MolecularSelector",
    "ProteinResidues",
    "ProteinResiduesByGroup",
    "ProteinResiduesNearReference",
    "PolymerChains",
    "PolymerResiduesByType",
    "SubstrateMolecule",
    "SolventMolecules",
    "CosolventMolecules",
    # Groupings
    "ResidueGrouping",
    "ProteinAAClassification",
    "CustomGrouping",
]
