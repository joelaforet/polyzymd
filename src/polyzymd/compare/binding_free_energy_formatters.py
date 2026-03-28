"""Output formatters for binding free energy comparison results.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    binding free energy formatters is now
    :mod:`polyzymd.analyses.binding_free_energy._formatters`.

All symbols are re-exported so that existing ``from polyzymd.compare.binding_free_energy_formatters import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.binding_free_energy._formatters import (  # noqa: F401
    format_bfe_console_table,
    format_bfe_json,
    format_bfe_markdown,
    format_bfe_result,
)

__all__ = [
    "format_bfe_console_table",
    "format_bfe_json",
    "format_bfe_markdown",
    "format_bfe_result",
]
