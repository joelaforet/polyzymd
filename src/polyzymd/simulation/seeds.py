"""Random seeds of the dynamics, derived from the replicate number.

The replicate number seeds a replicate's starting structure (Packmol and
polymer draws) and its dynamics: the initial velocities, the thermostat
noise and the Monte Carlo barostat moves. Each phase of a run (an
equilibration stage, a production segment) takes its own seed, so no phase
repeats another's noise, and the same replicate run again draws the same
numbers.

The same seeds give the same trajectory only with ``openmm.deterministic:
true`` on the same platform, precision and software versions. Otherwise CPU
threads and PME add up forces in a different order on each run, so two runs
of a replicate differ from the first minimization on and agree only
statistically.
"""

from __future__ import annotations

import hashlib

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
