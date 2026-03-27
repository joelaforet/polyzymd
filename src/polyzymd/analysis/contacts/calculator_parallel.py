"""Backward-compatibility shim — real implementation lives in analyses.contacts.

.. deprecated::
    Import from :mod:`polyzymd.analyses.contacts` instead::

        from polyzymd.analyses.contacts import ParallelContactAnalyzer
"""

from __future__ import annotations

# Re-export from canonical location
from polyzymd.analyses.contacts import (  # noqa: F401
    ParallelContactAnalyzer,
    _get_contact_analysis_base_cls,
    identify_polymer_chains,
)

__all__ = [
    "ParallelContactAnalyzer",
    "_get_contact_analysis_base_cls",
    "identify_polymer_chains",
]
