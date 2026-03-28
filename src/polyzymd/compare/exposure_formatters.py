"""Output formatters for exposure dynamics comparison results.

.. deprecated:: 1.3.0
    This module is a backward-compatibility shim.  The canonical location for
    exposure formatters is now :mod:`polyzymd.analyses.exposure._formatters`.

All symbols are re-exported so that existing ``from polyzymd.compare.exposure_formatters import …``
statements continue to work.
"""

from __future__ import annotations

# Re-export everything from the new canonical location
from polyzymd.analyses.exposure._formatters import (  # noqa: F401
    exposure_to_json,
    format_exposure_console_table,
    format_exposure_markdown,
    format_exposure_result,
)

__all__ = [
    "exposure_to_json",
    "format_exposure_console_table",
    "format_exposure_markdown",
    "format_exposure_result",
]
