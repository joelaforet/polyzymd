"""Integrity checks and transactional publication for OpenMM build artifacts."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import tempfile
import uuid
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

MANIFEST_NAME = "build_manifest.json"
_BUILD_ARTIFACTS = ("solvated_system.pdb", "system.xml")


class ArtifactIntegrityError(RuntimeError):
    """Raised when molecular artifacts cannot belong to the same build."""


def config_hash(config: Any) -> str:
    """Return a stable SHA-256 identity for a validated simulation config."""
    payload = config.model_dump(mode="json")
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def publish_build_bundle(
    working_dir: Path,
    topology: Any,
    system: Any,
    positions: Any,
    config: Any,
) -> dict[str, Any]:
    """Atomically publish a complete PDB/System bundle, committing manifest last."""
    from openmm import XmlSerializer, version
    from openmm.app import PDBFile

    working_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".build-bundle-", dir=working_dir) as staging_text:
        staging = Path(staging_text)
        pdb_path = staging / _BUILD_ARTIFACTS[0]
        system_path = staging / _BUILD_ARTIFACTS[1]
        with pdb_path.open("w") as stream:
            PDBFile.writeFile(topology, positions, stream, keepIds=True)
        system_path.write_text(XmlSerializer.serialize(system))

        manifest = {
            "schema_version": 1,
            "build_uuid": str(uuid.uuid4()),
            "config_hash": config_hash(config),
            "particle_count": int(system.getNumParticles()),
            "openmm_version": version.full_version,
            "artifacts": {
                name: {
                    "path": str((working_dir / name).resolve()),
                    "sha256": _file_hash(staging / name),
                }
                for name in _BUILD_ARTIFACTS
            },
        }
        previous = {
            name: (working_dir / name).read_bytes()
            for name in (*_BUILD_ARTIFACTS, MANIFEST_NAME)
            if (working_dir / name).exists()
        }
        try:
            # The manifest is the final commit marker for the new pair.
            for name in _BUILD_ARTIFACTS:
                os.replace(staging / name, working_dir / name)
            _atomic_write(
                working_dir / MANIFEST_NAME,
                (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode(),
            )
        except Exception:
            for name in (*_BUILD_ARTIFACTS, MANIFEST_NAME):
                path = working_dir / name
                if name in previous:
                    _atomic_write(path, previous[name])
                else:
                    path.unlink(missing_ok=True)
            raise
    return manifest


def validate_build_bundle(working_dir: Path, config: Any, *, allow_legacy: bool = True) -> None:
    """Validate a prebuilt bundle against its manifest and current configuration."""
    manifest_path = working_dir / MANIFEST_NAME
    paths = {name: working_dir / name for name in _BUILD_ARTIFACTS}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ArtifactIntegrityError(
            "Pre-built bundle is incomplete; missing: " + ", ".join(missing)
        )
    manifest = None
    if not manifest_path.is_file():
        if not allow_legacy:
            raise ArtifactIntegrityError(f"Build manifest is missing: {manifest_path}")
        warnings.warn(
            f"Legacy build has no {MANIFEST_NAME}; accepting only after particle-count validation",
            RuntimeWarning,
            stacklevel=2,
        )
    else:
        try:
            manifest = json.loads(manifest_path.read_text())
            artifacts = manifest["artifacts"]
            expected_config = manifest["config_hash"]
            int(manifest["particle_count"])
        except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
            raise ArtifactIntegrityError(f"Invalid build manifest {manifest_path}: {exc}") from exc
        if expected_config != config_hash(config):
            raise ArtifactIntegrityError(
                f"Configuration does not match {manifest_path}: expected {expected_config}, "
                f"got {config_hash(config)}"
            )
        for name, path in paths.items():
            recorded_path = artifacts.get(name, {}).get("path")
            if recorded_path != str(path.resolve()):
                raise ArtifactIntegrityError(
                    f"Artifact path mismatch for {name}: manifest={recorded_path!r}, "
                    f"current={str(path.resolve())!r}"
                )
            expected = artifacts.get(name, {}).get("sha256")
            actual = _file_hash(path)
            if expected != actual:
                raise ArtifactIntegrityError(
                    f"Artifact hash mismatch for {path}: manifest={expected!r}, actual={actual}"
                )

    from openmm import XmlSerializer
    from openmm.app import PDBFile

    pdb = PDBFile(str(paths["solvated_system.pdb"]))
    system = XmlSerializer.deserialize(paths["system.xml"].read_text())
    validate_openmm_identity(
        pdb.topology,
        system,
        topology_path=paths["solvated_system.pdb"],
        system_path=paths["system.xml"],
    )
    if manifest is not None and int(manifest["particle_count"]) != system.getNumParticles():
        raise ArtifactIntegrityError(
            f"Manifest particle count mismatch for {manifest_path}: "
            f"manifest={manifest['particle_count']}, system={system.getNumParticles()}"
        )


def validate_openmm_identity(
    topology: Any,
    system: Any,
    *,
    topology_path: Path,
    system_path: Path,
    state: Any | None = None,
    state_path: Path | None = None,
) -> None:
    """Fail with path-and-count diagnostics when OpenMM artifact counts disagree."""
    topology_count = sum(1 for _ in topology.atoms())
    system_count = int(system.getNumParticles())
    counts = {str(topology_path): topology_count, str(system_path): system_count}
    if state is not None:
        positions = state.getPositions()
        velocities = state.getVelocities()
        label = str(state_path or "state")
        counts[f"{label} positions"] = len(positions)
        counts[f"{label} velocities"] = len(velocities)
    if len(set(counts.values())) != 1:
        detail = "; ".join(f"{path}={count}" for path, count in counts.items())
        raise ArtifactIntegrityError(
            f"OpenMM artifact particle-count mismatch: {detail}. Rebuild into a new output "
            "directory or restore artifacts from the same completed build."
        )


def assert_rebuild_allowed(working_dir: Path) -> None:
    """Refuse replacement once any simulation phase or production artifact exists."""
    markers = [working_dir / "simulation_progress.json"]
    markers.extend(working_dir.glob("minimization*"))
    markers.extend(working_dir.glob("equilibration_*"))
    markers.extend(working_dir.glob("production*"))
    active = [path for path in markers if path.exists()]
    if active:
        raise ArtifactIntegrityError(
            f"Refusing to rebuild started campaign {working_dir}; found {active[0]}. "
            "Use a new output directory."
        )


@contextmanager
def replicate_lock(working_dir: Path) -> Iterator[None]:
    """Hold a non-blocking per-replicate build/run lock."""
    working_dir.mkdir(parents=True, exist_ok=True)
    lock_path = working_dir / ".polyzymd.lock"
    with lock_path.open("a+") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ArtifactIntegrityError(
                f"Another PolyzyMD build or run holds the replicate lock: {lock_path}"
            ) from exc
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
