"""Tests for build bundle identity and replicate serialization."""

from __future__ import annotations

import json
from multiprocessing import Process, Queue
from pathlib import Path

import pytest

from polyzymd.simulation.artifact_integrity import (
    ArtifactIntegrityError,
    publish_build_bundle,
    replicate_lock,
    validate_build_bundle,
    validate_openmm_identity,
)


class _Config:
    def __init__(self, name: str = "campaign") -> None:
        self.name = name

    def model_dump(self, *, mode: str) -> dict[str, str]:
        assert mode == "json"
        return {"name": self.name}


def _tiny_openmm_bundle(particles: int = 2):
    from openmm import System, Vec3, unit
    from openmm.app import Element, Topology

    topology = Topology()
    chain = topology.addChain("A")
    residue = topology.addResidue("HOH", chain)
    system = System()
    for index in range(particles):
        topology.addAtom(f"H{index}", Element.getByAtomicNumber(1), residue)
        system.addParticle(1.0)
    positions = [Vec3(float(index), 0, 0) for index in range(particles)] * unit.nanometer
    return topology, system, positions


def _contend_for_lock(path: str, queue: Queue) -> None:
    try:
        with replicate_lock(Path(path)):
            queue.put("acquired")
    except ArtifactIntegrityError:
        queue.put("blocked")


def test_manifest_rejects_hash_and_config_drift(tmp_path):
    topology, system, positions = _tiny_openmm_bundle()
    publish_build_bundle(tmp_path, topology, system, positions, _Config())
    validate_build_bundle(tmp_path, _Config())

    (tmp_path / "system.xml").write_text("stale")
    with pytest.raises(ArtifactIntegrityError, match="hash mismatch"):
        validate_build_bundle(tmp_path, _Config())

    publish_build_bundle(tmp_path, topology, system, positions, _Config())
    with pytest.raises(ArtifactIntegrityError, match="Configuration does not match"):
        validate_build_bundle(tmp_path, _Config("changed"))


def test_failed_publication_restores_previous_bundle(tmp_path, monkeypatch):
    import polyzymd.simulation.artifact_integrity as integrity

    topology, system, positions = _tiny_openmm_bundle()
    original = publish_build_bundle(tmp_path, topology, system, positions, _Config())
    real_replace = integrity.os.replace
    calls = 0

    def fail_between_artifact_replacements(source, destination):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected publication failure")
        return real_replace(source, destination)

    monkeypatch.setattr(integrity.os, "replace", fail_between_artifact_replacements)
    with pytest.raises(OSError, match="injected"):
        publish_build_bundle(tmp_path, topology, system, positions, _Config("replacement"))

    assert json.loads((tmp_path / "build_manifest.json").read_text()) == original
    validate_build_bundle(tmp_path, _Config())


def test_legacy_bundle_requires_consistent_counts(tmp_path):
    topology, system, positions = _tiny_openmm_bundle()
    publish_build_bundle(tmp_path, topology, system, positions, _Config())
    (tmp_path / "build_manifest.json").unlink()
    with pytest.warns(RuntimeWarning, match="Legacy build"):
        validate_build_bundle(tmp_path, _Config())

    _, wrong_system, _ = _tiny_openmm_bundle(3)
    from openmm import XmlSerializer

    (tmp_path / "system.xml").write_text(XmlSerializer.serialize(wrong_system))
    with (
        pytest.warns(RuntimeWarning),
        pytest.raises(ArtifactIntegrityError, match="particle-count"),
    ):
        validate_build_bundle(tmp_path, _Config())


def test_state_position_velocity_mismatch_is_rejected(tmp_path):
    topology, system, _ = _tiny_openmm_bundle(2)

    class State:
        def getPositions(self):
            return [object(), object()]

        def getVelocities(self):
            return [object()]

    with pytest.raises(ArtifactIntegrityError) as error:
        validate_openmm_identity(
            topology,
            system,
            topology_path=tmp_path / "topology.pdb",
            system_path=tmp_path / "system.xml",
            state=State(),
            state_path=tmp_path / "state.xml",
        )
    assert "topology.pdb=2" in str(error.value)
    assert "system.xml=2" in str(error.value)
    assert "state.xml velocities=1" in str(error.value)


def test_replicate_lock_blocks_concurrent_process(tmp_path):
    queue: Queue = Queue()
    with replicate_lock(tmp_path):
        process = Process(target=_contend_for_lock, args=(str(tmp_path), queue))
        process.start()
        process.join(5)
    assert process.exitcode == 0
    assert queue.get(timeout=1) == "blocked"
