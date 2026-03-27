"""Residue grouping abstractions for analysis modules.

.. deprecated::
    This package has moved to ``polyzymd.analyses.shared.groupings``.
    This shim re-exports all symbols for backward compatibility.
"""

from polyzymd.analyses.shared.groupings.base import (  # noqa: F401
    CustomGrouping,
    ProteinAAClassification,
    ResidueGrouping,
)

__all__ = [
    "ResidueGrouping",
    "ProteinAAClassification",
    "CustomGrouping",
]
