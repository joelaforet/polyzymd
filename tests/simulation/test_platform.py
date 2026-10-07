"""Tests for fail-fast OpenMM platform selection."""

from unittest.mock import MagicMock

import pytest

from polyzymd.simulation.platform import resolve_platform


def test_cuda_unavailable_never_falls_back(monkeypatch):
    """A missing requested CUDA platform must be fatal."""
    import openmm

    lookup = MagicMock(side_effect=openmm.OpenMMException("no CUDA"))
    monkeypatch.setattr(openmm.Platform, "getPlatformByName", lookup)

    with pytest.raises(RuntimeError, match="will not fall back to CPU"):
        resolve_platform("CUDA")

    lookup.assert_called_once_with("CUDA")


def test_cuda_properties_are_explicit(monkeypatch):
    """CUDA precision and device selection should reach Context creation."""
    import openmm

    platform = object()
    monkeypatch.setattr(openmm.Platform, "getPlatformByName", lambda name: platform)
    selection = resolve_platform("CUDA", precision="double", device_index="1")
    assert selection.platform is platform
    assert selection.properties == {"Precision": "double", "DeviceIndex": "1"}


def test_explicit_cpu_caps_threads_from_slurm(monkeypatch):
    """Explicit CPU selection should respect the scheduler allocation."""
    import openmm

    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setattr(openmm.Platform, "getPlatformByName", lambda name: object())
    assert resolve_platform("CPU").properties == {"Threads": "8"}


def test_deterministic_properties_for_each_platform(monkeypatch):
    """Deterministic runs set DeterministicForces on every GPU/CPU platform, and one CPU thread."""
    import openmm

    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setattr(openmm.Platform, "getPlatformByName", lambda name: object())
    deterministic = {"DeterministicForces": "true"}
    assert resolve_platform("CPU", deterministic=True).properties == {
        "Threads": "1",
        **deterministic,
    }
    assert resolve_platform("CUDA", deterministic=True).properties == {
        "Precision": "mixed",
        **deterministic,
    }
    assert resolve_platform("OpenCL", deterministic=True).properties == deterministic
    assert resolve_platform("Reference", deterministic=True).properties == {}


def _water_box():
    import openmm
    from openmm import app, unit

    forcefield = app.ForceField("amber14/tip3p.xml")
    modeller = app.Modeller(app.Topology(), [])
    modeller.addSolvent(forcefield, boxSize=openmm.Vec3(2.0, 2.0, 2.0) * unit.nanometer)
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=0.9 * unit.nanometer
    )
    return system, modeller.positions


def _minimize_and_run(system, positions):
    import numpy as np
    import openmm
    from openmm import unit

    from polyzymd.simulation.platform import platform_record

    selection = resolve_platform("CPU", deterministic=True)
    integrator = openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtosecond
    )
    integrator.setRandomNumberSeed(12345)
    context = openmm.Context(system, integrator, selection.platform, selection.properties)
    context.setPositions(positions)
    forces = context.getState(getForces=True).getForces(asNumpy=True)._value
    openmm.LocalEnergyMinimizer.minimize(context, 10.0, 200)
    context.setVelocitiesToTemperature(300 * unit.kelvin, 777)
    integrator.step(100)
    final = context.getState(getPositions=True).getPositions(asNumpy=True)._value
    return np.array(forces), np.array(final), platform_record(context)


def test_deterministic_cpu_runs_repeat_bitwise():
    """Two deterministic CPU contexts give identical forces, minimum and dynamics."""
    import numpy as np

    system, positions = _water_box()
    forces_a, final_a, record = _minimize_and_run(system, positions)
    forces_b, final_b, _ = _minimize_and_run(system, positions)
    assert np.array_equal(forces_a, forces_b)
    assert np.array_equal(final_a, final_b)
    assert record == {"name": "CPU", "properties": {"Threads": "1", "DeterministicForces": "true"}}
