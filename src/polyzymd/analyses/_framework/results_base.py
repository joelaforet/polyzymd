"""Version helpers for analysis artifact metadata.

The implementation lives in :mod:`polyzymd.utils.version` so that build and
simulation provenance records can share it; this module re-exports it for
backwards compatibility.
"""

from __future__ import annotations

from polyzymd.utils.version import get_polyzymd_version

__all__ = ["get_polyzymd_version"]
