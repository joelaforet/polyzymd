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


def rmsd(atoms: Any, reference: Any) -> float:
    """Return the RMSD of ``atoms`` from ``reference`` after optimal superposition.

    This calls ``MDAnalysis.analysis.rms.rmsd`` with ``center=True`` and
    ``superposition=True``, which moves both sets of coordinates to their
    centres of geometry, rotates ``atoms`` onto ``reference`` and returns the
    square root of the mean squared distance per atom, in Å, with every atom
    weighted equally. This is what ``MDAnalysis.analysis.rms.RMSD`` computes
    for one selection, as the legacy rmsd plugin did.

    Parameters
    ----------
    atoms : MDAnalysis.core.groups.AtomGroup
        Atoms to measure.
    reference : MDAnalysis.core.groups.AtomGroup
        The same number of atoms at their reference positions, usually from
        :func:`polyzymd.analyses.reference.reference`.

    Returns
    -------
    float
        RMSD in Å.
    """
    from MDAnalysis.analysis.rms import rmsd as _rmsd

    return float(_rmsd(atoms.positions, reference.positions, center=True, superposition=True))
