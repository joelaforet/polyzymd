"""Solvent, cosolvent, and substrate selectors.

.. deprecated::
    This module has moved to ``polyzymd.analyses.shared.selectors.solvent``.
    This shim re-exports all symbols for backward compatibility.
"""

from polyzymd.analyses.shared.selectors.solvent import (  # noqa: F401
    DEFAULT_COSOLVENT_RESNAMES,
    DEFAULT_WATER_RESNAMES,
    CosolventMolecules,
    IonSelector,
    SolventMolecules,
    SubstrateMolecule,
)
