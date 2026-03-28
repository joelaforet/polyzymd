"""Output formatters for polymer affinity score comparison results.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    polymer affinity formatters is now
    :mod:`polyzymd.analyses.polymer_affinity._formatters`.

All symbols are re-exported so that existing ``from polyzymd.compare.polymer_affinity_formatters import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.polymer_affinity._formatters import (  # noqa: F401
    format_affinity_console_table,
    format_affinity_json,
    format_affinity_markdown,
    format_affinity_result,
)

__all__ = [
    "format_affinity_console_table",
    "format_affinity_json",
    "format_affinity_markdown",
    "format_affinity_result",
]
