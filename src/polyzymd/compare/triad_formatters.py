"""Output formatters for triad comparison results.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    triad formatters is now :mod:`polyzymd.analyses.catalytic_triad._formatters`.

All symbols are re-exported so that existing ``from polyzymd.compare.triad_formatters import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.catalytic_triad._formatters import (  # noqa: F401
    format_triad_console_table,
    format_triad_markdown,
    format_triad_result,
    triad_to_json,
)

__all__ = [
    "format_triad_console_table",
    "format_triad_markdown",
    "format_triad_result",
    "triad_to_json",
]
