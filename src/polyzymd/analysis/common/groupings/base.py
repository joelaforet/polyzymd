"""Base classes for residue grouping/classification.

.. deprecated::
    This module has moved to ``polyzymd.analyses.shared.groupings.base``.
    This shim re-exports all symbols for backward compatibility.
"""

from polyzymd.analyses.shared.groupings.base import (  # noqa: F401
    CustomGrouping,
    ProteinAAClassification,
    ResidueGrouping,
)
