"""The cuPotential backend never changes what PolyzyMD writes, and a chain never mixes backends."""

from __future__ import annotations

import json
import sys

import pytest

openmm = pytest.importorskip("openmm")

from polyzymd.simulation import backends  # noqa: E402
from polyzymd.simulation.continuation import ContinuationManager  # noqa: E402


def _tiny_system():
    system = openmm.System()
    for _ in range(2):
        system.addParticle(1.0)
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(3, 0, 0), openmm.Vec3(0, 3, 0), openmm.Vec3(0, 0, 3)
    )
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.PME)
    nb.setCutoffDistance(0.9)
    for charge in (1.0, -1.0):
        nb.addParticle(charge, 0.3, 0.5)
    nb.setForceGroup(3)
    system.addForce(nb)
    system.addForce(openmm.CMMotionRemover())
    barostat = openmm.MonteCarloBarostat(1.0, 300.0, 25)
    barostat.setRandomNumberSeed(1234)
    system.addForce(barostat)
    return system


def test_openmm_backend_returns_the_stock_system_itself():
    system = _tiny_system()
    assert backends.runtime_system(system, "openmm") is system


def test_interrupt_and_restart_files_hold_the_stock_system(tmp_path, monkeypatch):
    """A Context running a converted System still writes and seeds from the stock one."""
    from openmm.app import Topology

    from polyzymd.simulation.seeds import openmm_seeds
    from polyzymd.simulation.signals import save_interrupted_state, save_restart_checkpoint

    stock = _tiny_system()

    def converted(system, backend):  # stands in for the nv conversion
        copy = openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system))
        copy.addForce(openmm.CustomExternalForce("0"))
        return copy

    monkeypatch.setattr(backends, "runtime_system", converted)
    topology = Topology()
    chain = topology.addChain()
    for _ in range(2):
        topology.addAtom("X", None, topology.addResidue("X", chain))
    simulation = backends.create_simulation(
        topology,
        stock,
        openmm.LangevinMiddleIntegrator(300, 1, 0.001),
        openmm.Platform.getPlatformByName("Reference"),
        {},
        backend="openmm_cupotential",
    )
    simulation.context.setPositions([openmm.Vec3(0, 0, 0), openmm.Vec3(1, 1, 1)])
    assert simulation.system.getNumForces() == stock.getNumForces() + 1
    assert backends.simulation_backend(simulation) == "openmm_cupotential"

    save_interrupted_state(simulation, tmp_path, 1, 0, 10)
    save_restart_checkpoint(simulation, tmp_path)
    expected = openmm.XmlSerializer.serialize(stock)
    assert (tmp_path / "interrupted_system.xml").read_text() == expected
    assert (tmp_path / "restart_system.xml").read_text() == expected
    assert openmm_seeds(simulation)["barostat"] == 1234


@pytest.mark.parametrize(
    ("recorded", "backend", "refused"),
    [
        (None, "openmm", False),
        (None, "openmm_cupotential", True),
        ("openmm_cupotential", "openmm_cupotential", False),
        ("openmm_cupotential", "openmm", True),
    ],
)
def test_continuation_refuses_a_backend_change(tmp_path, recorded, backend, refused):
    manager = ContinuationManager(tmp_path, 2, backend=backend)
    provenance = {"openmm_version": "8.4"} if recorded is None else {"backend": recorded}
    manager._param_dict = json.loads(json.dumps({"provenance": provenance}))
    if refused:
        with pytest.raises(RuntimeError, match="cannot change force backend"):
            manager._refuse_backend_change()
    else:
        manager._refuse_backend_change()


def test_cupotential_conversion_keeps_groups_and_barostat():
    pytest.importorskip("nv_openmm_plugin")
    stock = _tiny_system()
    before = openmm.XmlSerializer.serialize(stock)
    system = backends.runtime_system(stock, "openmm_cupotential")
    assert openmm.XmlSerializer.serialize(stock) == before
    names = [type(force).__name__ for force in system.getForces()]
    assert "CMMotionRemover" not in names
    nonbonded = next(f for f in system.getForces() if "Nonbonded" in type(f).__name__)
    assert nonbonded.getForceGroup() == 3
    barostat = next(f for f in system.getForces() if "Barostat" in type(f).__name__)
    assert barostat.getRandomNumberSeed() == 1234
    assert barostat.getFrequency() == 25


class _FakeNvIntegrator:
    def __init__(self, temperature, friction, step, comFrequency=100):  # noqa: N803
        self._com = comFrequency

    def getComFrequency(self):  # noqa: N802
        return self._com


@pytest.mark.parametrize("removers", [0, 1, 2])
def test_cupotential_integrator_removes_com_like_the_stock_system(monkeypatch, removers):
    import types

    monkeypatch.setitem(
        sys.modules,
        "nv_openmm_plugin",
        types.SimpleNamespace(LangevinMiddleIntegrator=_FakeNvIntegrator),
    )
    stock = _tiny_system()
    stock.removeForce(1)  # the CMMotionRemover (frequency 1, as Interchange writes)
    for _ in range(removers):
        stock.addForce(openmm.CMMotionRemover(1))
    make = lambda: backends.make_integrator(  # noqa: E731
        "LangevinMiddle", 300, 1, 0.002, "openmm_cupotential", stock
    )
    if removers == 2:
        with pytest.raises(ValueError, match="CMMotionRemovers"):
            make()
    else:
        assert make().getComFrequency() == removers


def test_missing_plugin_says_how_to_install_it(monkeypatch):
    monkeypatch.setitem(sys.modules, "nv_openmm_plugin", None)
    with pytest.raises(ImportError, match="not installed in this environment"):
        backends.make_integrator(
            "LangevinMiddle", 300, 1, 0.002, "openmm_cupotential", _tiny_system()
        )


@pytest.mark.parametrize("backend", ["openmm", "openmm_cupotential"])
def test_runner_refuses_a_backend_other_than_the_run_directory_recorded(tmp_path, backend):
    from polyzymd.simulation.phase_state import write_phase_record
    from polyzymd.simulation.runner import SimulationRunner

    write_phase_record(
        tmp_path / "minimization" / "phase.json", phase="minimization", status="completed"
    )
    write_phase_record(tmp_path / "heating" / "phase.json", phase="heating", status="recovery")
    runner = SimulationRunner.__new__(SimulationRunner)
    runner._working_dir, runner._backend = tmp_path, backend
    if backend == "openmm":
        runner._refuse_backend_change()  # a legacy record without backend ran openmm
    else:
        with pytest.raises(RuntimeError, match="heating ran on backend 'openmm'"):
            runner._refuse_backend_change()
