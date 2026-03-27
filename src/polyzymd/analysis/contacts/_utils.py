"""Backward-compatibility shim — real implementation lives in analyses.contacts.

.. deprecated::
    Import from :mod:`polyzymd.analyses.contacts` instead::

        from polyzymd.analyses.contacts import identify_polymer_chains
"""

from __future__ import annotations

# Re-export from canonical location
from polyzymd.analyses.contacts import identify_polymer_chains  # noqa: F401

__all__ = ["identify_polymer_chains"]
