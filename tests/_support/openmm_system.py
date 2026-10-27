"""Write serialized OpenMM systems, as PolyzyMD saves beside each trajectory segment.

A PolyzyMD OpenMM segment ``production_0/production_0_trajectory.dcd``
has its system in ``production_0/production_0_system.xml``, written by
``openmm.XmlSerializer``. :func:`write_openmm_system` builds such a system
with OpenMM itself, so the file has the layout OpenMM writes.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Sequence


def openmm_system_xml(
    charges: Sequence[float],
    bonds: Sequence[tuple[int, int]] = (),
    constraints: Sequence[tuple[int, int]] = (),
    *,
    decoys: bool = False,
) -> str:
    """Return the XML of an OpenMM system with one particle per charge.

    Parameters
    ----------
    charges : sequence of float
        Partial charge of each particle, in elementary charges, in its
        ``NonbondedForce``. Every particle has mass 12.
    bonds : sequence of (int, int)
        Pairs of particle indices in a ``HarmonicBondForce``.
    constraints : sequence of (int, int)
        Pairs of particle indices held at a fixed distance.
    decoys : bool, optional
        Also add a ``CustomBondForce`` with a bond between particles 0 and
        the last one and a ``GBSAOBCForce`` whose particles have charge 9,
        which a reader of bonds and charges must ignore.

    Returns
    -------
    str
        The serialized system.
    """
    mm = importlib.import_module("openmm")

    system = mm.System()
    for _ in charges:
        system.addParticle(12.0)
    for a, b in constraints:
        system.addConstraint(int(a), int(b), 0.1)
    harmonic = mm.HarmonicBondForce()
    for a, b in bonds:
        harmonic.addBond(int(a), int(b), 0.1, 1000.0)
    system.addForce(harmonic)
    if decoys:
        custom = mm.CustomBondForce("k*r^2")
        custom.addPerBondParameter("k")
        custom.addBond(0, len(charges) - 1, [1.0])
        system.addForce(custom)
        gbsa = mm.GBSAOBCForce()
        for _ in charges:
            gbsa.addParticle(9.0, 0.15, 1.0)
        system.addForce(gbsa)
    nonbonded = mm.NonbondedForce()
    for charge in charges:
        nonbonded.addParticle(float(charge), 0.3, 0.5)
    system.addForce(nonbonded)
    return mm.XmlSerializer.serialize(system)


def write_openmm_system(
    run_dir: Path,
    charges: Sequence[float],
    bonds: Sequence[tuple[int, int]] = (),
    constraints: Sequence[tuple[int, int]] = (),
    *,
    segment: str = "production_0",
    decoys: bool = False,
) -> Path:
    """Write ``<run_dir>/<segment>/<segment>_system.xml`` for :func:`openmm_system_xml`.

    Returns
    -------
    Path
        The written file.
    """
    path = Path(run_dir) / segment / f"{segment}_system.xml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(openmm_system_xml(charges, bonds, constraints, decoys=decoys))
    return path
