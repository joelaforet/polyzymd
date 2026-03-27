"""Polymer chain and residue selectors.

.. deprecated::
    This module has moved to ``polyzymd.analyses.shared.selectors.polymer``.
    This shim re-exports all symbols for backward compatibility.
"""

from polyzymd.analyses.shared.selectors.polymer import (  # noqa: F401
    DEFAULT_POLYMER_RESNAMES,
    POLYZYMD_POLYMER_CHAIN_ID,
    PolymerChains,
    PolymerResiduesByType,
    PolymerSegments,
)
