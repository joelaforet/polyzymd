"""Output formatters for contacts comparison results.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    contacts formatters is now :mod:`polyzymd.analyses.contacts._formatters`.

All symbols are re-exported so that existing ``from polyzymd.compare.contacts_formatters import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.contacts._formatters import (  # noqa: F401
    contacts_to_json,
    format_contacts_console_table,
    format_contacts_markdown,
    format_contacts_result,
)

__all__ = [
    "contacts_to_json",
    "format_contacts_console_table",
    "format_contacts_markdown",
    "format_contacts_result",
]
