"""Tests for build bundle identity and replicate serialization."""

from __future__ import annotations

import json
from multiprocessing import Process, Queue
from pathlib import Path

import pytest

from polyzymd.simulation.artifact_integrity import (
    MANIFEST_NAME,
    ArtifactIntegrityError,
    _absolute_path_config_hash,
    assert_rebuild_allowed,
    config_hash,
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
        # Only the bundle's own files count, not entries of the file-hash cache.
        if Path(destination).parent == tmp_path:
            calls += 1
            if calls == 2:
                raise OSError("injected publication failure")
        return real_replace(source, destination)

    monkeypatch.setattr(integrity.os, "replace", fail_between_artifact_replacements)
    with pytest.raises(OSError, match="injected"):
        publish_build_bundle(tmp_path, topology, system, positions, _Config("replacement"))

    assert json.loads((tmp_path / "build_manifest.json").read_text()) == original
    validate_build_bundle(tmp_path, _Config())


def test_bundle_without_manifest_is_refused_with_rebuild_advice(tmp_path):
    topology, system, positions = _tiny_openmm_bundle()
    publish_build_bundle(tmp_path, topology, system, positions, _Config())
    (tmp_path / "build_manifest.json").unlink()
    with pytest.raises(ArtifactIntegrityError, match="Build manifest is missing") as excinfo:
        validate_build_bundle(tmp_path, _Config())
    message = str(excinfo.value)
    assert "polyzymd build -c <config> -r <replicate>" in message
    assert "run the same command again" in message


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


def test_manifest_records_provenance_and_versions(tmp_path):
    topology, system, positions = _tiny_openmm_bundle()
    manifest = publish_build_bundle(
        tmp_path,
        topology,
        system,
        positions,
        _Config(),
        provenance={"packmol_seed": 3, "polymer_seed": 3},
    )
    assert manifest["provenance"] == {"packmol_seed": 3, "polymer_seed": 3}
    assert isinstance(manifest["polyzymd_version"], str)
    assert manifest["openmm_version"]
    on_disk = json.loads((tmp_path / MANIFEST_NAME).read_text())
    assert on_disk["provenance"]["packmol_seed"] == 3


def test_manifest_provenance_defaults_to_empty(tmp_path):
    topology, system, positions = _tiny_openmm_bundle()
    manifest = publish_build_bundle(tmp_path, topology, system, positions, _Config())
    assert manifest["provenance"] == {}


QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart"


def _copied_quickstart_config(folder: Path):
    import shutil

    from polyzymd.config.loader import load_config

    folder.mkdir()
    for name in ("config.yaml", "trpcage.pdb"):
        shutil.copy(QUICKSTART / name, folder / name)
    return load_config(folder / "config.yaml")


def test_config_hash_is_the_same_for_copies_in_different_folders(tmp_path):
    first = _copied_quickstart_config(tmp_path / "a")
    second = _copied_quickstart_config(tmp_path / "b")
    assert config_hash(first) == config_hash(second)
    assert _absolute_path_config_hash(first) != _absolute_path_config_hash(second)

    (tmp_path / "b" / "trpcage.pdb").write_text("REMARK changed\n")
    assert config_hash(first) != config_hash(second)


def test_build_validates_with_a_copied_config_and_with_the_absolute_path_hash(tmp_path):
    first = _copied_quickstart_config(tmp_path / "a")
    second = _copied_quickstart_config(tmp_path / "b")
    run = tmp_path / "run"
    topology, system, positions = _tiny_openmm_bundle()
    publish_build_bundle(run, topology, system, positions, first)
    validate_build_bundle(run, second)

    # A build written before the portable hash recorded the absolute-path hash.
    manifest = json.loads((run / MANIFEST_NAME).read_text())
    manifest["config_hash"] = _absolute_path_config_hash(first)
    (run / MANIFEST_NAME).write_text(json.dumps(manifest))
    validate_build_bundle(run, first)
    with pytest.raises(ArtifactIntegrityError, match="Configuration does not match"):
        validate_build_bundle(run, second)


def test_a_moved_build_validates(tmp_path):
    """Artifacts are recorded relative to the run folder, so a moved run folder still validates."""
    import shutil

    topology, system, positions = _tiny_openmm_bundle()
    manifest = publish_build_bundle(tmp_path / "run", topology, system, positions, _Config())
    assert manifest["artifacts"]["system.xml"]["path"] == "system.xml"
    moved = Path(shutil.move(tmp_path / "run", tmp_path / "moved"))
    validate_build_bundle(moved, _Config())

    # A manifest written before records absolute paths, which a move leaves behind.
    recorded = json.loads((moved / MANIFEST_NAME).read_text())
    for name, artifact in recorded["artifacts"].items():
        artifact["path"] = str(tmp_path / "run" / name)
    (moved / MANIFEST_NAME).write_text(json.dumps(recorded))
    validate_build_bundle(moved, _Config())

    recorded["artifacts"]["system.xml"]["path"] = "other.xml"
    (moved / MANIFEST_NAME).write_text(json.dumps(recorded))
    with pytest.raises(ArtifactIntegrityError, match="path mismatch"):
        validate_build_bundle(moved, _Config())


def test_validation_reads_an_artifact_replaced_with_the_same_size_and_time(tmp_path, monkeypatch):
    import os

    from polyzymd.analyses.shared.file_hashes import cache_dir

    monkeypatch.setenv("POLYZYMD_CACHE_DIR", str(tmp_path / "cache"))
    topology, system, positions = _tiny_openmm_bundle()
    publish_build_bundle(tmp_path, topology, system, positions, _Config())
    validate_build_bundle(tmp_path, _Config())
    system_xml = tmp_path / "system.xml"
    stat = system_xml.stat()
    text = system_xml.read_text()
    system_xml.write_text(text.replace("<", "[", 1))
    os.utime(system_xml, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    with pytest.raises(ArtifactIntegrityError, match="hash mismatch"):
        validate_build_bundle(tmp_path, _Config())
    # Staging files are deleted after the build, so their hashes are not cached.
    cached = [json.loads(p.read_text())["path"] for p in cache_dir().glob("*.json")]
    assert not any(".build-bundle-" in path for path in cached)


def test_gromacs_build_manifest_records_the_exported_inputs(tmp_path):
    """A GROMACS build writes build_manifest.json beside the replicate, as an OpenMM build does."""
    from polyzymd.analyses.shared.file_hashes import file_sha256
    from polyzymd.simulation.artifact_integrity import write_gromacs_build_manifest

    gromacs = tmp_path / "gromacs"
    gromacs.mkdir()
    (tmp_path / "solvated_system.pdb").write_text("pdb")
    (gromacs / "system.gro").write_text("gro")
    (gromacs / "system.top").write_text('#include "system_MOL0.itp"\n')
    (gromacs / "system_MOL0.itp").write_text("itp")
    (gromacs / "prod.mdp").write_text("ld_seed = 7\n")
    (gromacs / "backup.top").write_text("not an input")

    write_gromacs_build_manifest(tmp_path, gromacs, _Config(), 3, {"packmol_seed": 1})

    manifest = json.loads((tmp_path / MANIFEST_NAME).read_text())
    assert manifest["config_hash"] == config_hash(_Config())
    assert manifest["particle_count"] == 3
    assert manifest["openmm_version"] is None
    assert manifest["polyzymd_version"]
    assert manifest["provenance"] == {"packmol_seed": 1}
    assert sorted(manifest["artifacts"]) == [
        "gromacs/prod.mdp",
        "gromacs/system.gro",
        "gromacs/system.top",
        "gromacs/system_MOL0.itp",
        "solvated_system.pdb",
    ]
    assert manifest["artifacts"]["gromacs/system.gro"] == {
        "path": "gromacs/system.gro",
        "sha256": file_sha256(gromacs / "system.gro", use_cache=False),
    }


@pytest.mark.parametrize(
    "marker", ["gromacs/em.log", "gromacs/eq_01.tpr", "gromacs/prod.cpt", "gromacs/state.cpt"]
)
def test_rebuild_is_refused_once_a_gromacs_run_has_started(tmp_path, marker):
    (tmp_path / "gromacs").mkdir()
    (tmp_path / marker).write_text("")

    with pytest.raises(ArtifactIntegrityError, match="Refusing to rebuild"):
        assert_rebuild_allowed(tmp_path)


def test_rebuild_is_allowed_over_gromacs_build_outputs(tmp_path):
    gromacs = tmp_path / "gromacs"
    gromacs.mkdir()
    for name in (
        "system.gro",
        "system.top",
        "system_MOL0.itp",
        "em.mdp",
        "eq_01_nvt.mdp",
        "prod.mdp",
        "run_system_gromacs.sh",
    ):
        (gromacs / name).write_text("")
    (tmp_path / MANIFEST_NAME).write_text("{}")
    (tmp_path / "progress.json").write_text("{}")

    assert_rebuild_allowed(tmp_path)


def test_rebuild_is_refused_once_an_openmm_run_has_started(tmp_path):
    (tmp_path / "equilibration_0").mkdir()

    with pytest.raises(ArtifactIntegrityError, match="equilibration_0"):
        assert_rebuild_allowed(tmp_path)
