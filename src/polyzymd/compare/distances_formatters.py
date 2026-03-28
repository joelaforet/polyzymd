"""Output formatters for distance comparison results.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    distances formatters is now :mod:`polyzymd.analyses.distances._formatters`.

All symbols are re-exported so that existing ``from polyzymd.compare.distances_formatters import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.distances._formatters import (  # noqa: F401
    distances_to_json,
    format_distances_console_table,
    format_distances_markdown,
    format_distances_result,
)

__all__ = [
    "distances_to_json",
    "format_distances_console_table",
    "format_distances_markdown",
    "format_distances_result",
]
