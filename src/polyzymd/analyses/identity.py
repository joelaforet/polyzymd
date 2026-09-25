"""What a replicate result was computed from, and whether that still holds.

Every replicate artifact carries an identity block: the input files and their
sizes and modification times, the equilibration window, the settings, the
source of the plugin and of the framework code it ran through, and the
simulation config. :func:`identity_mismatch` compares a stored block with the
current one and is the only test the runner applies before it reuses a cached
replicate or aggregates one from disk, so a change to any of those recomputes
the number.

The PolyzyMD version and the git commit are recorded as provenance but not
compared. The code hashes already change whenever code that can change a
number changes, and a version bump on its own should not throw away a
campaign's cached replicates. A version difference is logged instead.
"""

from __future__ import annotations

import ast
import hashlib
import importlib
import inspect
import json
import logging
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping, Sequence

from pydantic import BaseModel

from polyzymd.analyses.exceptions import PluginContractError

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig

logger = logging.getLogger("polyzymd.analyses")

#: Identity fields that must match before a cached replicate is reused.
COMPARED_KEYS = (
    "inputs",
    "equilibration",
    "settings_fingerprint",
    "plugin_code_hash",
    "framework_code_hash",
    "config_hash",
    "settings_files",
)

#: Framework modules every plugin computes through, whatever it imports.
_FRAMEWORK_MODULES = (
    "polyzymd.analyses.contract",
    "polyzymd.analyses.orchestrator",
    "polyzymd.analyses.loading",
)

#: Package prefixes the framework hash walks into. ``mda`` is included because
#: ``mda/frame_selection.py`` decides which frames a plugin sees and
#: ``mda/universe.py`` decides what it reads them from.
_HASHED_PREFIXES = (
    "polyzymd.analyses.shared.",
    "polyzymd.analyses.mda.",
)

#: Modules inside a hashed prefix that only affect figures. A change to a figure
#: must not throw away every cached replicate.
_FIGURE_MODULES = frozenset({"polyzymd.analyses.shared.plotting"})


def replicate_identity(
    plugin: Any,
    sim_config: SimulationConfig,
    settings: BaseModel,
    equilibration: str,
    inputs: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Build the identity block of one replicate.

    Parameters
    ----------
    plugin : AnalysisProtocol
        Plugin instance that computes the replicate.
    sim_config : SimulationConfig
        Configuration of the condition.
    settings : BaseModel
        Plugin settings.
    equilibration : str
        Equilibration window, for example ``"10ns"``.
    inputs : sequence of mapping
        Identities of the topology and trajectory files, as
        :func:`input_files` returns them.

    Returns
    -------
    dict
        JSON-compatible identity block.
    """
    from polyzymd import __version__

    return {
        "polyzymd_version": __version__,
        "plugin": plugin.name,
        "git_commit": _git_commit(),
        "git_dirty": _git_dirty(),
        "plugin_code_hash": code_hash(plugin),
        "framework_code_hash": framework_code_hash(type(plugin).__module__),
        "settings_fingerprint": settings_fingerprint(settings),
        "config_hash": compute_config_hash(sim_config),
        "equilibration": str(equilibration),
        "inputs": [dict(entry) for entry in inputs],
        "settings_files": settings_file_identity(plugin, settings),
    }


def identity_mismatch(stored: Mapping[str, Any] | None, current: Mapping[str, Any]) -> str | None:
    """Say why a stored identity no longer describes the current run.

    Parameters
    ----------
    stored : mapping or None
        Identity block read from a cached artifact.
    current : mapping
        Identity block of the run that wants to reuse it.

    Returns
    -------
    str or None
        ``None`` when every compared field matches, otherwise one sentence
        naming the first field that differs. An empty or missing stored block
        cannot be shown to match, so it is a mismatch.
    """
    if not stored or not stored.get("inputs"):
        return "it records no input file identity"
    for key in COMPARED_KEYS:
        if key == "inputs":
            reason = _inputs_mismatch(stored.get("inputs") or [], current.get("inputs") or [])
            if reason is not None:
                return reason
        elif stored.get(key) != current.get(key):
            return f"{key} changed: cached {stored.get(key)!r}, current {current.get(key)!r}"
    return None


def input_files(provenance: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    """Topology and trajectory identities from a universe provider's provenance."""
    if not isinstance(provenance, Mapping):
        return []
    topology = provenance.get("topology")
    files = ([topology] if isinstance(topology, Mapping) else []) + list(
        provenance.get("trajectories") or []
    )
    return [dict(entry) for entry in files if isinstance(entry, Mapping)]


def restat(inputs: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """The recorded input files as they are on disk now.

    Used where the runner has only the artifact to go on, when aggregating
    replicate results a worker wrote. A file that is gone is recorded as
    missing, so the comparison reports it.
    """
    from polyzymd.analyses.mda.universe import FileIdentity

    current: list[dict[str, Any]] = []
    for entry in inputs:
        path = Path(str(entry.get("path", "")))
        if path.is_file():
            current.append(FileIdentity.from_path(path, entry.get("format")).as_dict())
        else:
            current.append({"path": str(path), "missing": True})
    return current


def warn_on_version_mismatch(recorded: Any, source: Path | str) -> str | None:
    """Log when a cached artifact was written by another PolyzyMD version.

    Returns the message, or ``None`` when the versions agree or none is
    recorded.
    """
    from polyzymd import __version__ as running_version

    if not isinstance(recorded, str) or not recorded or recorded == running_version:
        return None
    message = (
        f"{source} was written by PolyzyMD {recorded} but this process runs "
        f"{running_version}; the cached numbers may not match what this version computes"
    )
    logger.warning("%s", message)
    return message


def compute_config_hash(config: SimulationConfig) -> str:
    """Hash the simulation config fields that change how a trajectory is read.

    The enzyme, substrate and polymer definitions, the thermodynamic state and
    the output locations are included. Simulation phases and force fields are
    left out, because they are already baked into a finished trajectory.

    Returns
    -------
    str
        First 16 hex characters of a SHA-256 over the canonical JSON.
    """
    hash_data: dict[str, Any] = {
        "name": config.name,
        "enzyme": {"name": config.enzyme.name, "pdb_path": str(config.enzyme.pdb_path)},
        "thermodynamics": {
            "temperature": config.thermodynamics.temperature,
            "pressure": config.thermodynamics.pressure,
        },
        "output": {
            "projects_directory": str(config.output.projects_directory),
            "scratch_directory": str(config.output.effective_scratch_directory),
            "naming_template": config.output.naming_template,
        },
    }
    if config.substrate is not None:
        hash_data["substrate"] = {
            "name": config.substrate.name,
            "sdf_path": str(config.substrate.sdf_path),
        }
    if config.polymers is not None and config.polymers.enabled:
        hash_data["polymers"] = {
            "type_prefix": config.polymers.type_prefix,
            "length": config.polymers.length,
            "count": config.polymers.count,
            "monomers": [
                {"label": m.label, "probability": m.probability, "name": m.name}
                for m in config.polymers.monomers
            ],
        }
    serialized = json.dumps(hash_data, sort_keys=True, default=str)
    return hashlib.sha256(serialized.encode()).hexdigest()[:16]


def settings_fingerprint(settings: BaseModel | None) -> str | None:
    """First 8 hex characters of a SHA-256 over the canonical settings JSON."""
    if settings is None:
        return None
    serialized = json.dumps(settings.model_dump(mode="json"), sort_keys=True)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:8]


def settings_file_identity(plugin: Any, settings: BaseModel) -> list[dict[str, Any]]:
    """Identity of the extra input files a plugin's settings name.

    A plugin whose answer depends on a file the framework does not load, such
    as an external reference structure, declares
    ``identity_files(settings) -> Sequence[Path]``. A file that is missing is
    recorded as missing rather than skipped, so it appearing later also
    invalidates the replicate.
    """
    from polyzymd.analyses.mda.universe import FileIdentity

    declare = getattr(plugin, "identity_files", None)
    if declare is None:
        return []
    entries: list[dict[str, Any]] = []
    for path in declare(settings):
        resolved = Path(path).expanduser()
        entries.append(
            FileIdentity.from_path(resolved).as_dict()
            if resolved.exists()
            else {"path": str(resolved), "missing": True}
        )
    return entries


def code_hash(plugin: Any) -> str:
    """Short SHA-256 of the module defining the plugin, so any edit invalidates."""
    plugin_type = type(plugin)
    for target in (inspect.getmodule(plugin_type), plugin_type):
        try:
            source = inspect.getsource(target)
        except (OSError, TypeError):
            continue
        return hashlib.sha256(source.encode("utf-8")).hexdigest()[:16]
    return "unknown"


@lru_cache(maxsize=None)
def framework_code_hash(plugin_module: str) -> str:
    """Hash the framework and shared code a plugin's answer depends on.

    Covers the contract, the runner, the loader, and every module under
    ``analyses/shared/`` and ``analyses/mda/`` the plugin reaches, so a fix in
    alignment, frame selection or the reduction rules invalidates the replicate.
    Plotting is excluded.

    Raises
    ------
    PluginContractError
        If the source of a module in the walk cannot be read. A hash that
        silently skipped it would compare equal to one taken before the module
        changed.
    """
    names = sorted(set(_FRAMEWORK_MODULES) | _hashed_imports(plugin_module))
    digest = hashlib.sha256()
    for name in names:
        source = _module_source(name)
        if source is None:
            raise PluginContractError(
                f"cannot read the source of {name}, so the cache identity of "
                f"{plugin_module} cannot be computed. Run from a source checkout or an "
                "installed package that ships its .py files, not a zipped or frozen build."
            )
        digest.update(name.encode("utf-8"))
        digest.update(source.encode("utf-8"))
    return digest.hexdigest()[:16]


def _inputs_mismatch(
    stored: Sequence[Mapping[str, Any]], current: Sequence[Mapping[str, Any]]
) -> str | None:
    """Name the first input file that changed, appeared or vanished."""
    stored_by_path = {str(entry.get("path")): entry for entry in stored}
    current_by_path = {str(entry.get("path")): entry for entry in current}
    for path, entry in stored_by_path.items():
        now = current_by_path.get(path)
        if now is None:
            return f"{path} is no longer among the inputs"
        if now.get("missing"):
            return f"{path} is missing or unreadable"
        if now.get("size_bytes") != entry.get("size_bytes"):
            return (
                f"{path} changed size: recorded {entry.get('size_bytes')} bytes, "
                f"found {now.get('size_bytes')}"
            )
        if now.get("mtime_ns") != entry.get("mtime_ns"):
            return (
                f"{path} changed modification time: recorded {entry.get('mtime_ns')} ns, "
                f"found {now.get('mtime_ns')} ns"
            )
    for path in current_by_path:
        if path not in stored_by_path:
            return f"{path} is a new input the cached result never read"
    return None


def _hashed_imports(module_name: str) -> set[str]:
    """Every hashed module reachable from one module's imports, transitively."""
    seen: set[str] = set()
    queue = [module_name, *_FRAMEWORK_MODULES]
    visited: set[str] = set()
    while queue:
        name = queue.pop()
        if name in visited:
            continue
        visited.add(name)
        source = _module_source(name)
        if source is None:
            continue
        for imported in _imported_names(source, name):
            if not imported.startswith(_HASHED_PREFIXES):
                continue
            if imported in seen or imported in _FIGURE_MODULES:
                continue
            if _module_source(imported) is None:
                continue
            seen.add(imported)
            queue.append(imported)
    return seen


def _imported_names(source: str, module_name: str) -> set[str]:
    """Dotted module names a source file imports, relative imports resolved."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return set()
    package = module_name.rsplit(".", 1)[0]
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                base = f"{package}.{base}" if base else package
            names.add(base)
            names.update(f"{base}.{alias.name}" for alias in node.names)
    return {name for name in names if name}


def _module_source(module_name: str) -> str | None:
    """Source of an importable module, or ``None`` when it has none on disk."""
    try:
        module = importlib.import_module(module_name)
    except Exception:
        return None
    try:
        return inspect.getsource(module)
    except (OSError, TypeError):
        return None


@lru_cache(maxsize=1)
def _git_commit() -> str | None:
    """Commit of the checkout PolyzyMD runs from, or ``None`` for an install."""
    return _git("rev-parse", "HEAD") or None


@lru_cache(maxsize=1)
def _git_dirty() -> bool | None:
    """Whether that checkout had uncommitted changes, or ``None`` for an install."""
    if _git_commit() is None:
        return None
    return bool(_git("status", "--porcelain"))


def _git(*args: str) -> str | None:
    """Run one git command in the source tree, or return ``None`` without git."""
    try:
        completed = subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parent), *args],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None
