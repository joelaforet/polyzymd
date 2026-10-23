"""Integrity checks and transactional publication for build artifacts."""

from __future__ import annotations

import errno
import fcntl
import hashlib
import json
import os
import tempfile
import uuid
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Sequence

from polyzymd.analyses.shared.file_hashes import file_sha256
from polyzymd.simulation.analysis_topology import ANALYSIS_TOPOLOGY_NAME, write_analysis_topology

#: The build record: provenance, components, physics audit and artifact hashes.
MANIFEST_NAME = "build_provenance.json"
#: The same record for people (Markdown).
PROVENANCE_MARKDOWN = "build_provenance.md"
#: Name of the record before 1.3.0rc5, still read for existing runs.
LEGACY_MANIFEST_NAME = "build_manifest.json"
_BUILD_ARTIFACTS = ("solvated_system.pdb", "system.xml")
# Written when ParmEd can convert the system; its absence is not an error.
_OPTIONAL_ARTIFACTS = (ANALYSIS_TOPOLOGY_NAME,)


class ArtifactIntegrityError(RuntimeError):
    """Raised when molecular artifacts cannot belong to the same build."""


def build_record_path(working_dir: Path) -> Path:
    """The run folder's build record: ``build_provenance.json``, or the older
    ``build_manifest.json`` of a run built before 1.3.0rc5."""
    current = Path(working_dir) / MANIFEST_NAME
    legacy = Path(working_dir) / LEGACY_MANIFEST_NAME
    return legacy if legacy.is_file() and not current.is_file() else current


def _digest(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _portable(value: Any, key: str | None = None) -> Any:
    from polyzymd.config.loader import OUTPUT_PATH_KEYS, PATH_KEYS

    if isinstance(value, dict):
        return {k: _portable(v, k) for k, v in value.items() if k not in OUTPUT_PATH_KEYS}
    if isinstance(value, list):
        return [_portable(v, key) for v in value]
    if key in PATH_KEYS and isinstance(value, str) and Path(value).is_absolute():
        path = Path(value)
        return file_sha256(path) if path.is_file() else path.name
    return value


def config_hash(config: Any) -> str:
    """Return a SHA-256 identity for a validated simulation config that does not depend on its folder.

    The hash covers the config's fields. An input file (such as ``enzyme.pdb_path``) counts by its SHA-256, a
    folder by its name, and the output folders (``projects_directory``,
    ``scratch_directory``) are left out. So copies of one config.yaml in two
    folders give the same hash.
    """
    return _digest(_portable(config.model_dump(mode="json")))


def _absolute_path_config_hash(config: Any) -> str:
    """Return the hash that earlier builds recorded, which includes absolute paths."""
    return _digest(config.model_dump(mode="json"))


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


def _manifest(
    config: Any,
    particle_count: int,
    openmm_version: str | None,
    provenance: dict[str, Any] | None,
    artifacts: dict[str, Path],
    engine: str = "openmm",
    components: Sequence[Any] | None = None,
    audit: Any = None,
) -> dict[str, Any]:
    """Return the build record; ``artifacts`` maps each recorded path to the file to hash.

    ``components`` (their records) and ``audit`` (the physics audit of the files
    the engine will run) are recorded when given.
    """
    from polyzymd.utils.version import build_versions, get_polyzymd_version

    return {
        "schema_version": 2,
        "engine": engine,
        "build_uuid": str(uuid.uuid4()),
        "config_hash": config_hash(config),
        "particle_count": int(particle_count),
        "openmm_version": openmm_version,
        **build_versions(),
        "polyzymd_version": get_polyzymd_version(),
        "provenance": dict(provenance or {}),
        "artifacts": {
            name: {"path": name, "sha256": file_sha256(path, use_cache=False)}
            for name, path in artifacts.items()
        },
        "components": [c.to_dict() for c in components or []],
        "audit": audit.to_dict() if audit is not None else None,
    }


def _audit(realize: Any, components: Sequence[Any] | None, config: Any) -> Any:
    """The physics audit of the system ``realize()`` reads, or None without components.

    Reading the system is part of the audit: if it cannot be read, or any check
    cannot run, the audit is one FAIL saying so, never a pass.
    """
    if not components:
        return None
    from polyzymd.audit.findings import Finding
    from polyzymd.audit.report import AuditReport, audit

    try:
        realized = realize() if callable(realize) else realize
        return audit(realized, components, config)
    except Exception as error:  # an audit that cannot run must not pass silently
        return AuditReport(
            [
                Finding(
                    "audit",
                    "system",
                    "FAIL",
                    f"the physics audit could not run: {type(error).__name__}: {error}",
                )
            ]
        )


def provenance_markdown(record: dict[str, Any], audit: Any = None) -> str:
    """The build record for people: versions, components and the physics audit."""
    lines = [
        "# Build provenance",
        "",
        f"- engine: {record.get('engine')}",
        f"- PolyzyMD {record.get('polyzymd_version')}; OpenMM {record.get('openmm_version')}",
        f"- config hash `{record.get('config_hash')}`; {record.get('particle_count')} atoms",
        f"- build {record.get('build_uuid')}",
        "",
        "## Components",
        "",
        "| Component | Type | Segment | Molecules | Force field | Charges |",
        "|---|---|---|---|---|---|",
    ]
    for c in record.get("components") or []:
        fields = sorted({", ".join(d["force_field"]) for d in c["domains"]})
        charges = sorted(
            {
                f"{d['charges']['method']}"
                + (f" ({d['charges']['model']})" if d["charges"].get("model") else "")
                for d in c["domains"]
            }
        )
        lines.append(
            f"| {c['name']} | {c['type']} | {c['segid']} | {c.get('n_molecules', len(c['molecules']))} | "
            f"{'; '.join(fields)} | {'; '.join(charges)} |"
        )
    if audit is not None:
        lines += ["", audit.to_markdown()]
    elif record.get("audit") is None:
        lines += ["", "No physics audit was recorded for this build."]
    seeds = record.get("provenance") or {}
    if seeds:
        lines += ["", "## Inputs and seeds", ""] + [
            f"- {k}: `{v}`" for k, v in sorted(seeds.items())
        ]
    return "\n".join(lines) + "\n"


def write_gromacs_build_manifest(
    replicate_dir: Path,
    gromacs_dir: Path,
    config: Any,
    particle_count: int,
    provenance: dict[str, Any] | None = None,
    components: Sequence[Any] | None = None,
) -> dict[str, Any]:
    """Write ``build_provenance.json`` (and ``.md``) for a GROMACS build into the replicate folder.

    With ``components``, the record carries the physics audit of the exported
    ``.top``/``.gro`` read under the production ``prod.mdp``, the files GROMACS
    will run.

    It has the keys of an OpenMM build's manifest, with ``openmm_version``
    null and the OpenFF versions that parameterized the system. The artifacts are ``solvated_system.pdb``, the exported
    ``<prefix>.gro``, and the topology, the files it includes and the MDP
    files (:func:`polyzymd.analyses.shared.gromacs.run_input_files`), each
    by its path relative to the replicate folder.
    """
    from polyzymd.analyses.shared.gromacs import run_input_files, system_prefix

    candidates = [
        replicate_dir / _BUILD_ARTIFACTS[0],
        gromacs_dir / _BUILD_ARTIFACTS[0],
        gromacs_dir / f"{system_prefix(config)}.gro",
        *run_input_files(gromacs_dir, config),
    ]
    artifacts = {
        path.relative_to(replicate_dir).as_posix(): path for path in candidates if path.is_file()
    }
    report = None
    if components:
        from polyzymd.audit.realized import RealizedSystem

        prefix = system_prefix(config)
        report = _audit(
            lambda: RealizedSystem.from_gromacs(
                gromacs_dir / f"{prefix}.top",
                gromacs_dir / f"{prefix}.gro",
                gromacs_dir / "prod.mdp",
            ),
            components,
            config,
        )
    manifest = _manifest(
        config, particle_count, None, provenance, artifacts, "gromacs", components, report
    )
    # The Markdown first, the JSON last: the JSON is the record readers trust.
    _atomic_write(
        replicate_dir / PROVENANCE_MARKDOWN, provenance_markdown(manifest, report).encode()
    )
    _atomic_write(
        replicate_dir / MANIFEST_NAME,
        (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode(),
    )
    return manifest


def publish_build_bundle(
    working_dir: Path,
    topology: Any,
    system: Any,
    positions: Any,
    config: Any,
    provenance: dict[str, Any] | None = None,
    segids: Sequence[str] | None = None,
    components: Sequence[Any] | None = None,
) -> dict[str, Any]:
    """Atomically publish a complete PDB/System bundle, committing its record last.

    With ``components``, the staged ``system.xml`` and PDB are audited
    (:mod:`polyzymd.audit`) before the bundle is committed, and the record
    (``build_provenance.json``, with ``build_provenance.md`` for people) carries
    the components and the audit.

    ``segids`` (optional) is each atom's component segment ID, written into
    columns 73-76 of the PDB.

    ``provenance`` (optional) is recorded verbatim under the manifest's
    ``"provenance"`` key alongside the PolyzyMD version; builders use it for
    the Packmol/polymer seeds that determine the starting coordinates.
    """
    from openmm import XmlSerializer, version
    from openmm.app import PDBFile

    working_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".build-bundle-", dir=working_dir) as staging_text:
        staging = Path(staging_text)
        pdb_path = staging / _BUILD_ARTIFACTS[0]
        system_path = staging / _BUILD_ARTIFACTS[1]
        with pdb_path.open("w") as stream:
            PDBFile.writeFile(topology, positions, stream, keepIds=True)
        if segids is not None:
            from polyzymd.core.components import write_pdb_segids

            write_pdb_segids(pdb_path, segids)
        system_path.write_text(XmlSerializer.serialize(system))
        written = list(_BUILD_ARTIFACTS)
        if write_analysis_topology(topology, system, positions, staging / ANALYSIS_TOPOLOGY_NAME):
            written.append(ANALYSIS_TOPOLOGY_NAME)

        report = None
        if components:
            from polyzymd.audit.realized import RealizedSystem

            step = config.simulation_phases.production.time_step
            report = _audit(
                lambda: RealizedSystem.from_openmm(system_path, pdb_path, step), components, config
            )
        manifest = _manifest(
            config,
            system.getNumParticles(),
            version.full_version,
            provenance,
            {name: staging / name for name in written},
            "openmm",
            components,
            report,
        )
        (staging / PROVENANCE_MARKDOWN).write_text(provenance_markdown(manifest, report))
        written.append(PROVENANCE_MARKDOWN)
        previous = {
            name: (working_dir / name).read_bytes()
            for name in (
                *_BUILD_ARTIFACTS,
                *_OPTIONAL_ARTIFACTS,
                PROVENANCE_MARKDOWN,
                MANIFEST_NAME,
            )
            if (working_dir / name).exists()
        }
        try:
            # The manifest is the final commit marker for the new pair.
            for name in written:
                os.replace(staging / name, working_dir / name)
            _atomic_write(
                working_dir / MANIFEST_NAME,
                (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode(),
            )
        except Exception:
            for name in (
                *_BUILD_ARTIFACTS,
                *_OPTIONAL_ARTIFACTS,
                PROVENANCE_MARKDOWN,
                MANIFEST_NAME,
            ):
                path = working_dir / name
                if name in previous:
                    _atomic_write(path, previous[name])
                else:
                    path.unlink(missing_ok=True)
            raise
    return manifest


def validate_build_bundle(working_dir: Path, config: Any) -> None:
    """Validate a prebuilt bundle against its build record and current configuration."""
    manifest_path = build_record_path(working_dir)
    paths = {name: working_dir / name for name in _BUILD_ARTIFACTS}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ArtifactIntegrityError(
            "Pre-built bundle is incomplete; missing: " + ", ".join(missing)
        )
    if not manifest_path.is_file():
        raise ArtifactIntegrityError(
            f"Build record is missing: {manifest_path}. Rebuild with "
            "'polyzymd build -c <config> -r <replicate>', then run the same command again."
        )
    try:
        manifest = json.loads(manifest_path.read_text())
        artifacts = manifest["artifacts"]
        expected_config = manifest["config_hash"]
        int(manifest["particle_count"])
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        raise ArtifactIntegrityError(f"Invalid build record {manifest_path}: {exc}") from exc
    # Builds written before the portable hash recorded the absolute-path hash.
    if expected_config not in (config_hash(config), _absolute_path_config_hash(config)):
        raise ArtifactIntegrityError(
            f"Configuration does not match {manifest_path}: expected {expected_config}, "
            f"got {config_hash(config)}"
        )
    checked = dict(paths)
    for name in _OPTIONAL_ARTIFACTS:
        if name in artifacts:
            checked[name] = working_dir / name
    for name, path in checked.items():
        if name in _OPTIONAL_ARTIFACTS and not path.is_file():
            raise ArtifactIntegrityError(
                f"Manifest lists {name} but {path} is missing; rerun the build or "
                "polyzymd analysis-topology"
            )
        # Paths are relative to the run folder, so a moved run folder validates.
        # Earlier manifests recorded absolute paths, compared by file name.
        recorded_path = artifacts.get(name, {}).get("path")
        if recorded_path is None or Path(recorded_path).name != name:
            raise ArtifactIntegrityError(
                f"Artifact path mismatch for {name}: manifest={recorded_path!r}, current={name!r}"
            )
        expected = artifacts.get(name, {}).get("sha256")
        actual = file_sha256(path, use_cache=False)
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
    if int(manifest["particle_count"]) != system.getNumParticles():
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
    """Refuse replacement once an OpenMM or GROMACS run has written output.

    GROMACS run outputs live in ``working_dir/gromacs``; the ``.mdp`` files
    there are build outputs.
    """
    markers = [*working_dir.glob("minimization*")]
    markers.extend(working_dir.glob("equilibration_*"))
    markers.extend(working_dir.glob("production*"))
    for pattern in ("em.*", "eq_[0-9][0-9].*", "prod*", "state.cpt"):
        markers.extend(
            path for path in (working_dir / "gromacs").glob(pattern) if path.suffix != ".mdp"
        )
    active = sorted(markers)
    if active:
        raise ArtifactIntegrityError(
            f"Refusing to rebuild started campaign {working_dir}; found {active[0]}. "
            "Use a new output directory."
        )


# Lock files this process holds through replicate_lock().  A POSIX lock is
# dropped when the process closes *any* descriptor of the file, so this
# process must never open a lock file it holds.
_HELD_LOCKS: set[Path] = set()


# lockf() fails with one of these when another process holds the lock.
_LOCK_HELD_ERRNOS = (errno.EAGAIN, errno.EACCES)


def _lock_path(working_dir: Path) -> Path:
    # resolve(): another spelling of the same folder must find the held lock.
    return (Path(working_dir) / ".polyzymd.lock").resolve()


@contextmanager
def replicate_lock(working_dir: Path) -> Iterator[None]:
    """Hold a non-blocking per-replicate build/run lock.

    The lock is a POSIX record lock (``fcntl.lockf``).  GPFS and NFS share
    these between nodes, while ``flock`` locks are seen only on the node that
    took them.  The kernel releases the lock when the holding process dies,
    however it is killed.
    """
    working_dir.mkdir(parents=True, exist_ok=True)
    lock_path = _lock_path(working_dir)
    if lock_path in _HELD_LOCKS:
        raise ArtifactIntegrityError(
            f"Another PolyzyMD build or run holds the replicate lock: {lock_path}"
        )
    with lock_path.open("a+") as stream:
        try:
            fcntl.lockf(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            if exc.errno not in _LOCK_HELD_ERRNOS:
                raise
            raise ArtifactIntegrityError(
                f"Another PolyzyMD build or run holds the replicate lock: {lock_path}"
            ) from exc
        _HELD_LOCKS.add(lock_path)
        try:
            yield
        finally:
            _HELD_LOCKS.discard(lock_path)
            fcntl.lockf(stream.fileno(), fcntl.LOCK_UN)


def replicate_lock_held_elsewhere(working_dir: Path) -> bool | None:
    """Return whether another process holds the replicate lock of *working_dir*.

    A run holds this lock for as long as its process lives, so the answer
    cannot be stale.  Returns ``False`` when this process holds the lock or
    no lock file exists, and ``None`` when the lock file cannot be opened.
    """
    lock_path = _lock_path(working_dir)
    if lock_path in _HELD_LOCKS or not lock_path.exists():
        return False
    try:
        with lock_path.open("a+") as stream:
            try:
                fcntl.lockf(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as exc:
                return True if exc.errno in _LOCK_HELD_ERRNOS else None
            fcntl.lockf(stream.fileno(), fcntl.LOCK_UN)
            return False
    except OSError:
        return None
