"""Per-frame measurements shipped with PolyzyMD, as plain functions.

Each function takes MDAnalysis ``AtomGroup`` arguments positioned at one frame
and returns one number, so it runs through
:meth:`polyzymd.analyses.study.Study.timeseries` like any function you write.
"""

from __future__ import annotations

from typing import Any


def radius_of_gyration(atoms: Any) -> float:
    """Return the mass-weighted radius of gyration of ``atoms`` at the current frame.

    This calls ``AtomGroup.radius_of_gyration()``, which weights each atom by
    its mass, uses the coordinates as loaded without unwrapping molecules
    split across periodic boundaries, and returns Å.

    Parameters
    ----------
    atoms : MDAnalysis.core.groups.AtomGroup
        Atoms to measure.

    Returns
    -------
    float
        Radius of gyration in Å.
    """
    return float(atoms.radius_of_gyration())
