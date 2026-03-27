"""Base class for molecular selectors.

.. deprecated::
    This module has moved to ``polyzymd.analyses.shared.selectors.base``.
    This shim re-exports all symbols for backward compatibility.
"""

from polyzymd.analyses.shared.selectors.base import (  # noqa: F401
    CompositeSelector,
    MDAnalysisSelector,
    MolecularSelector,
    SelectionResult,
)
