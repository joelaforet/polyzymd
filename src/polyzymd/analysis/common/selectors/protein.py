"""Protein residue selectors.

.. deprecated::
    This module has moved to ``polyzymd.analyses.shared.selectors.protein``.
    This shim re-exports all symbols for backward compatibility.
"""

from polyzymd.analyses.shared.selectors.protein import (  # noqa: F401
    ProteinResidues,
    ProteinResiduesByGroup,
    ProteinResiduesNearReference,
)
