"""Amino acid classification and SASA reference data.

.. deprecated::
    This module has moved to ``polyzymd.analyses.shared.aa_classification``.
    This shim re-exports all symbols for backward compatibility with legacy
    ``analysis/`` code. New code should import from the canonical location.
"""

from polyzymd.analyses.shared.aa_classification import (  # noqa: F401
    AA_CLASS_RESIDUES,
    AA_CLASSIFICATION_TABLE,
    CANONICAL_AA_CLASS_ORDER,
    DEFAULT_AA_CLASS_SELECTIONS,
    MAX_ASA_TABLE,
    STANDARD_AA_CODES,
    AAClass,
    get_aa_class,
    get_max_asa,
    get_residues_for_class,
    get_selection_for_class,
)
