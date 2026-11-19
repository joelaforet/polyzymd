"""Run OpenMM dynamics on stock forces or on NVIDIA's cuPotential plugin.

The System PolyzyMD builds, keeps and writes to disk is always stock OpenMM.
For the ``openmm_cupotential`` backend, :func:`create_simulation` converts a private
copy right before the Context is made, so ``system.xml``, segment system
files and saved States never depend on the backend. Energy minimization
always runs on stock forces.

``openmm_cupotential`` is all or nothing: ``nv.NonbondedForce``,
``nv.MonteCarloBarostat`` and ``nv.LangevinMiddleIntegrator``, with the
CMMotionRemover removed because the nv integrator removes centre-of-mass
motion itself. ``openmm`` is stock OpenMM for everything.
"""

from __future__ import annotations

import logging
from typing import Any

LOGGER = logging.getLogger(__name__)

#: Ewald error tolerance cuPotential 0.3.5 requires on the replaced NonbondedForce.
CUPOTENTIAL_EWALD_TOLERANCE = 5e-4

BACKENDS = ("openmm", "openmm_cupotential")

PLUGIN_MISSING = (
    "openmm.backend 'openmm_cupotential' needs NVIDIA's nv_openmm_plugin, which is not installed in "
    "this environment. PolyzyMD expects it from NVIDIA's conda package nv-openmm-plugin 0.3.5 "
    "(https://github.com/NVIDIA/cuPotential/releases, file "
    "nv-openmm-plugin-conda-channel-0.3.5.tar.gz, unpacked as a local conda channel). It needs "
    "glibc >= 2.38 and a CUDA 13.4 driver; on CURC, run inside the openmm-cupotential Apptainer "
    "image. See docs/source/how_to/run_with_cupotential.md."
)


def _plugin() -> Any:
    """Import ``nv_openmm_plugin``, or raise an ImportError that says how to get it."""
    try:
        import nv_openmm_plugin
    except ImportError as exc:
        raise ImportError(PLUGIN_MISSING) from exc
    return nv_openmm_plugin


def refuse_backend_change(recorded: str | None, requested: str, where: str) -> None:
    """Raise when ``where`` ran another backend; a record without one ran ``openmm``.

    A run never mixes stock OpenMM and cuPotential dynamics.
    """
    recorded = recorded or "openmm"
    if recorded != requested:
        raise RuntimeError(
            f"{where} ran on backend {recorded!r}, but this job would run on {requested!r}. "
            f"A run cannot change force backend; set openmm.backend to {recorded!r}, "
            "or start a new run."
        )


def com_frequency(stock: Any) -> int:
    """Return the stock System's CMMotionRemover frequency, or 0 (no COM removal) without one."""
    import openmm

    removers = [f for f in stock.getForces() if isinstance(f, openmm.CMMotionRemover)]
    if len(removers) > 1:
        raise ValueError(f"Found {len(removers)} CMMotionRemovers; cuPotential needs at most one")
    return removers[0].getFrequency() if removers else 0


def stock_system(simulation: Any) -> Any:
    """Return the stock System behind ``simulation``, the only System PolyzyMD writes."""
    # vars(), not getattr(): a Mock invents any attribute it is asked for.
    stock = vars(simulation).get("polyzymd_stock_system")
    return simulation.system if stock is None else stock


def simulation_backend(simulation: Any) -> str:
    """Return the backend ``simulation`` was created with (``openmm`` by default)."""
    return vars(simulation).get("polyzymd_backend", "openmm")


def runtime_system(stock: Any, backend: str = "openmm") -> Any:
    """Return the System a Context should run for ``backend``.

    ``stock`` itself for ``openmm``; for ``openmm_cupotential`` a converted private
    copy (``stock`` is not modified). Force groups and the barostat pressure,
    temperature, frequency and random seed are kept.
    """
    if backend == "openmm":
        return stock
    if backend not in BACKENDS:
        raise ValueError(f"Unknown OpenMM backend {backend!r}")
    nv = _plugin()
    import openmm

    system = openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(stock))
    forces = system.getForces()
    nb_indices = [i for i, f in enumerate(forces) if isinstance(f, openmm.NonbondedForce)]
    if len(nb_indices) != 1:
        raise ValueError(f"cuPotential needs exactly one NonbondedForce, found {len(nb_indices)}")
    nb = forces[nb_indices[0]]
    method = nb.getNonbondedMethod()
    if method not in (openmm.NonbondedForce.PME, openmm.NonbondedForce.LJPME):
        raise ValueError("cuPotential accepts only a PME or LJPME NonbondedForce")
    if not nb.getIncludeDirectSpace():
        raise ValueError("cuPotential does not support NonbondedForce includeDirectSpace=False")
    if nb.getEwaldErrorTolerance() != CUPOTENTIAL_EWALD_TOLERANCE:
        LOGGER.warning(
            f"cuPotential: Ewald tolerance {nb.getEwaldErrorTolerance()} -> "
            f"{CUPOTENTIAL_EWALD_TOLERANCE} (required by the plugin)"
        )
        nb.setEwaldErrorTolerance(CUPOTENTIAL_EWALD_TOLERANCE)
    nb.setPMEParameters(0, 0, 0, 0)
    if method == openmm.NonbondedForce.LJPME:
        nb.setLJPMEParameters(0, 0, 0, 0)
    replacement = nv.NonbondedForce(nb)
    replacement.setForceGroup(nb.getForceGroup())
    remove, add = [nb_indices[0]], [replacement]

    for i, force in enumerate(forces):
        if isinstance(force, openmm.CMMotionRemover):
            remove.append(i)
        elif isinstance(force, openmm.MonteCarloBarostat):
            barostat = nv.MonteCarloBarostat(
                force.getDefaultPressure(), force.getDefaultTemperature(), force.getFrequency()
            )
            barostat.setRandomNumberSeed(force.getRandomNumberSeed())
            barostat.setForceGroup(force.getForceGroup())
            remove.append(i)
            add.append(barostat)
        elif "Barostat" in type(force).__name__:
            raise ValueError(f"{type(force).__name__} has no cuPotential equivalent")
    for i in sorted(set(remove), reverse=True):
        system.removeForce(i)
    for force in add:
        system.addForce(force)
    return system


def make_integrator(
    thermostat: str,
    temperature: Any,
    friction: Any,
    timestep: Any,
    backend: str = "openmm",
    stock: Any = None,
) -> Any:
    """Return the dynamics integrator: nv.LangevinMiddleIntegrator for openmm_cupotential.

    ``thermostat`` is ``"LangevinMiddle"`` or ``"Langevin"``; ``temperature``,
    ``friction`` and ``timestep`` are OpenMM quantities. For openmm_cupotential,
    ``stock`` is the stock System: the nv integrator removes centre-of-mass
    motion at its CMMotionRemover frequency (:func:`com_frequency`), since
    :func:`runtime_system` removes that force.
    """
    import openmm

    if backend == "openmm_cupotential":
        if thermostat != "LangevinMiddle":
            raise ValueError(
                f"The cuPotential backend supports only the LangevinMiddle thermostat, "
                f"not {thermostat!r}"
            )
        if stock is None:
            raise ValueError("The cuPotential integrator needs the stock System")
        frequency = com_frequency(stock)
        integrator = _plugin().LangevinMiddleIntegrator(
            temperature, friction, timestep, comFrequency=frequency
        )
        if integrator.getComFrequency() != frequency:
            raise RuntimeError(
                f"nv.LangevinMiddleIntegrator removes COM every {integrator.getComFrequency()} "
                f"steps, not {frequency} like the stock System"
            )
        return integrator
    if thermostat == "Langevin":
        return openmm.LangevinIntegrator(temperature, friction, timestep)
    if thermostat == "LangevinMiddle":
        return openmm.LangevinMiddleIntegrator(temperature, friction, timestep)
    raise ValueError(f"Unknown thermostat {thermostat!r}")


def create_simulation(
    topology: Any,
    stock: Any,
    integrator: Any,
    platform: Any,
    properties: Any,
    backend: str = "openmm",
) -> Any:
    """Create a Simulation for ``backend`` that remembers its stock System and backend."""
    from openmm.app import Simulation

    simulation = Simulation(
        topology, runtime_system(stock, backend), integrator, platform, properties
    )
    simulation.polyzymd_stock_system = stock
    simulation.polyzymd_backend = backend
    return simulation
