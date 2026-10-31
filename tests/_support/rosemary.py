"""The newest OpenFF Rosemary candidate installed, for tests that build on it.

OpenFF publishes Rosemary candidates as ``openff_no_water-3.0.0-alpha*``; the
tests take the newest one installed instead of naming a file, so a new
candidate needs no test changes.
"""

from __future__ import annotations

import re

from polyzymd.config.schema import is_rosemary

#: The oldest candidate, which PolyzyMD keeps accepting.
ROSEMARY_ALPHA0 = "openff_no_water-3.0.0-alpha0.offxml"


def _order(name: str) -> tuple:
    """Sort key of a candidate: alpha number, then its letter suffix (alpha2b > alpha2)."""
    match = re.search(r"alpha(\d+)([a-z]*)\.offxml$", name)
    return (int(match.group(1)), match.group(2)) if match else (-1, "")


def newest_rosemary() -> str:
    """File name of the newest installed Rosemary candidate (water-free, constrained)."""
    from openff.toolkit.typing.engines.smirnoff.forcefield import get_available_force_fields

    candidates = [
        name
        for name in get_available_force_fields()
        if is_rosemary([name])
        and name.startswith("openff_no_water-")
        and "unconstrained" not in name
    ]
    if not candidates:
        raise RuntimeError("no OpenFF Rosemary candidate is installed")
    return max(candidates, key=_order)
