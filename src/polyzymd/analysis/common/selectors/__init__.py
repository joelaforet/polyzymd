"""Molecular selector abstractions for analysis modules.

This module provides a unified interface for selecting atoms, residues, or
molecular groups from an MDAnalysis Universe. Selectors enable:

1. Consistent selection logic across different analysis types
2. Configurable protein/polymer/solvent selection
3. Support for arbitrary user-defined selections
4. Proximity-based selections (e.g., "residues near active site")

The Strategy pattern is used so users can define custom selectors by
subclassing MolecularSelector.

Examples
--------
>>> from polyzymd.analysis.common.selectors import ProteinResidues, PolymerChains
>>>
>>> # Select all protein residues
>>> protein_selector = ProteinResidues()
>>> protein_residues = protein_selector.select(universe)
>>>
>>> # Select polymer chains by type
>>> polymer_selector = PolymerResiduesByType(residue_names=["SBM", "EGP"])
>>> polymer_residues = polymer_selector.select(universe)
>>>
>>> # Select protein residues near catalytic triad
>>> triad_selector = ProteinResiduesNearReference(
...     reference_selection="resid 77 133 156",
...     cutoff=5.0
... )
>>> nearby_residues = triad_selector.select(universe)
"""

from polyzymd.analysis.common.selectors.base import (
    MolecularSelector,
    MDAnalysisSelector,
    SelectionResult,
)
from polyzymd.analysis.common.selectors.protein import (
    ProteinResidues,
    ProteinResiduesByGroup,
    ProteinResiduesNearReference,
)
from polyzymd.analysis.common.selectors.polymer import (
    PolymerChains,
    PolymerResiduesByType,
)
from polyzymd.analysis.common.selectors.solvent import (
    SolventMolecules,
    CosolventMolecules,
    SubstrateMolecule,
)

__all__ = [
    # Base
    "MolecularSelector",
    "MDAnalysisSelector",
    "SelectionResult",
    # Protein
    "ProteinResidues",
    "ProteinResiduesByGroup",
    "ProteinResiduesNearReference",
    # Polymer
    "PolymerChains",
    "PolymerResiduesByType",
    # Solvent/Substrate
    "SolventMolecules",
    "CosolventMolecules",
    "SubstrateMolecule",
]
