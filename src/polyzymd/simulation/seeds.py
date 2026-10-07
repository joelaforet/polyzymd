"""Random seeds of the dynamics, derived from the replicate number.

The replicate number seeds a replicate's starting structure (Packmol and
polymer draws) and its dynamics: the initial velocities, the thermostat
noise and the Monte Carlo barostat moves. Each phase of a run (an
equilibration stage, a production segment) takes its own seed, so no phase
repeats another's noise, and the same replicate run again draws the same
numbers.

The same seeds do not give the same trajectory. OpenMM's CPU threads, PME
and GPUs add up forces in a different order on each run, so two runs of a
replicate differ from the first minimization on and agree statistically,
not frame by frame.
"""

from __future__ import annotations

import hashlib
from typing import Any

#: Largest seed both engines accept (a positive 32-bit integer).
MAX_SEED = 2**31 - 1


def dynamics_seed(replicate: int, phase: str) -> int:
    """Return the seed of ``phase`` of replicate ``replicate``, from 1 to :data:`MAX_SEED`.

    ``phase`` names one draw, such as ``"velocities"``,
    ``"equilibration:2"`` or ``"production:5"``. The seed is the SHA-256 of
    the two, so it never depends on the machine or the engine. It is never
    0, which OpenMM reads as "choose a random seed" (GROMACS reads -1 so).
    """
    digest = hashlib.sha256(f"{int(replicate)}:{phase}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % MAX_SEED + 1


def openmm_seeds(simulation: Any, velocities: int | None = None) -> dict[str, int]:
    """Return the random seeds an OpenMM ``Simulation`` runs with, for ``progress.json``.

    ``integrator`` is the thermostat noise seed and ``barostat`` the seed of
    the Monte Carlo volume moves, read back from the integrator and the
    system; 0 means OpenMM chose a random one. ``velocities`` is the seed the
    phase drew its starting velocities with, given only when it drew them.
    """
    seeds = {"integrator": int(simulation.integrator.getRandomNumberSeed())}
    for force in simulation.system.getForces():
        if "Barostat" in type(force).__name__:
            seeds["barostat"] = int(force.getRandomNumberSeed())
    if velocities is not None:
        seeds["velocities"] = int(velocities)
    return seeds
