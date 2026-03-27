"""Molecular selector abstractions for analysis modules.

.. deprecated::
    This package has moved to ``polyzymd.analyses.shared.selectors``.
    This shim re-exports all symbols for backward compatibility with legacy
    ``analysis/`` code. New code should import from the canonical location.
"""

from polyzymd.analyses.shared.selectors.base import (  # noqa: F401
    CompositeSelector,
    MDAnalysisSelector,
    MolecularSelector,
    SelectionResult,
)
from polyzymd.analyses.shared.selectors.polymer import (  # noqa: F401
    PolymerChains,
    PolymerResiduesByType,
    PolymerSegments,
)
from polyzymd.analyses.shared.selectors.protein import (  # noqa: F401
    ProteinResidues,
    ProteinResiduesByGroup,
    ProteinResiduesNearReference,
)
from polyzymd.analyses.shared.selectors.solvent import (  # noqa: F401
    CosolventMolecules,
    IonSelector,
    SolventMolecules,
    SubstrateMolecule,
)

__all__ = [
    # Base
    "MolecularSelector",
    "MDAnalysisSelector",
    "SelectionResult",
    "CompositeSelector",
    # Protein
    "ProteinResidues",
    "ProteinResiduesByGroup",
    "ProteinResiduesNearReference",
    # Polymer
    "PolymerChains",
    "PolymerResiduesByType",
    "PolymerSegments",
    # Solvent/Substrate
    "SolventMolecules",
    "CosolventMolecules",
    "SubstrateMolecule",
    "IonSelector",
]
