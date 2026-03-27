"""Common infrastructure for analysis modules.

.. deprecated::
    This package has moved to ``polyzymd.analyses.shared``.
    Sub-packages:
    - ``aa_classification`` -> ``polyzymd.analyses.shared.aa_classification``
    - ``selectors`` -> ``polyzymd.analyses.shared.selectors``
    - ``groupings`` -> ``polyzymd.analyses.shared.groupings``

    This shim re-exports all symbols for backward compatibility with legacy
    ``analysis/`` code. New code should import from the canonical locations.
"""

from polyzymd.analyses.shared.aa_classification import (
    AA_CLASS_RESIDUES,
    AA_CLASSIFICATION_TABLE,
    DEFAULT_AA_CLASS_SELECTIONS,
    MAX_ASA_TABLE,
    STANDARD_AA_CODES,
    AAClass,
    get_aa_class,
    get_max_asa,
    get_residues_for_class,
    get_selection_for_class,
)
from polyzymd.analyses.shared.groupings import (
    CustomGrouping,
    ProteinAAClassification,
    ResidueGrouping,
)
from polyzymd.analyses.shared.selectors import (
    CosolventMolecules,
    MolecularSelector,
    PolymerChains,
    PolymerResiduesByType,
    ProteinResidues,
    ProteinResiduesByGroup,
    ProteinResiduesNearReference,
    SolventMolecules,
    SubstrateMolecule,
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
