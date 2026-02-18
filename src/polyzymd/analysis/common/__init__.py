"""Common infrastructure for analysis modules.

This package provides shared components used across multiple analysis modules:

- **selectors**: Molecular selection abstractions (protein, polymer, solvent, etc.)
- **groupings**: Residue classification systems (amino acid types, custom groups)
- **aa_classification**: Amino acid properties (maxASA, classifications, selections)

These components enable consistent, extensible analysis across different
molecular systems while maintaining separation of concerns.
"""

from polyzymd.analysis.common.aa_classification import (
    AAClass,
    AA_CLASSIFICATION_TABLE,
    AA_CLASS_RESIDUES,
    DEFAULT_AA_CLASS_SELECTIONS,
    MAX_ASA_TABLE,
    STANDARD_AA_CODES,
    get_aa_class,
    get_max_asa,
    get_residues_for_class,
    get_selection_for_class,
)
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
    # AA Classification & SASA
    "AAClass",
    "AA_CLASSIFICATION_TABLE",
    "AA_CLASS_RESIDUES",
    "DEFAULT_AA_CLASS_SELECTIONS",
    "MAX_ASA_TABLE",
    "STANDARD_AA_CODES",
    "get_aa_class",
    "get_max_asa",
    "get_residues_for_class",
    "get_selection_for_class",
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
