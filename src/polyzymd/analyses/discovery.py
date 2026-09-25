"""Discovery of analysis plugins, built in and user written.

Built-in analyses are found by walking ``polyzymd.analyses``: every
non-infrastructure top-level module or package is imported and its concrete
:class:`~polyzymd.analyses.base.Analysis` subclasses are collected.

Analyses written outside PolyzyMD live in a study's ``analyses/`` folder (see
:mod:`polyzymd.config.study`). :func:`load_analysis_directory` imports each
module there and registers every analysis it defines, whether the module wraps
it with :func:`~polyzymd.analyses.contract.contract_analysis` or only defines
the plugin class. :func:`register_analysis` does the same for one analysis
defined in a script. A registered analysis may not take the name of a built-in
one, so a name in a comparison always means the same code.

Name collisions between two plugins raise immediately.
"""

from __future__ import annotations

import hashlib
import importlib
import importlib.util
import inspect
import logging
import pkgutil
import sys
import types
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any

from polyzymd.analyses.exceptions import PluginContractError

if TYPE_CHECKING:
    from polyzymd.analyses.base import Analysis

logger = logging.getLogger("polyzymd.analyses")

# Modules that are infrastructure, not plugins
_SKIP_MODULES = frozenset(
    {
        "shared",
        "base",
        "stats",
        "discovery",
        "orchestrator",
        "exceptions",
        "protocols",
        "mda",
        "runner",
        "config",
    }
)

# Heavy optional dependencies that may not be installed in all environments
# ImportError for these is expected and should be silently skipped
_OPTIONAL_HEAVY_DEPS = frozenset(
    {
        "openmm",
        "openff",
        "MDAnalysis",
        "mdanalysis",
        "parmed",
        "pdbfixer",
        "espaloma_charge",
        "dgl",
        "torch",
        "ambertools",
    }
)


def _is_concrete_analysis(obj: type) -> bool:
    """Return True if *obj* is a concrete (non-abstract) Analysis subclass."""
    from polyzymd.analyses.base import Analysis

    return (
        inspect.isclass(obj)
        and issubclass(obj, Analysis)
        and obj is not Analysis
        and not getattr(obj, "__abstractmethods__", None)
    )


def _should_skip_module(modname: str, package_prefix: str) -> bool:
    """Return True when module path includes skipped components.

    Parameters
    ----------
    modname : str
        Fully qualified module name discovered by ``pkgutil``.
    package_prefix : str
        Base package prefix including trailing dot, for example
        ``"polyzymd.analyses."``.

    Returns
    -------
    bool
        True if any path component is private (starts with ``"_"``)
        or listed in ``_SKIP_MODULES``.
    """
    relative_name = modname
    if modname.startswith(package_prefix):
        relative_name = modname[len(package_prefix) :]
    components = relative_name.split(".")
    return any(component.startswith("_") or component in _SKIP_MODULES for component in components)


def _is_top_level_module(modname: str, package_prefix: str) -> bool:
    """Return whether *modname* is a direct module under ``polyzymd.analyses``.

    Parameters
    ----------
    modname : str
        Fully qualified module name discovered by ``pkgutil``.
    package_prefix : str
        Base package prefix including trailing dot, for example
        ``"polyzymd.analyses."``.

    Returns
    -------
    bool
        ``True`` when the relative module name has no package separator.
    """

    relative_name = modname
    if modname.startswith(package_prefix):
        relative_name = modname[len(package_prefix) :]
    return "." not in relative_name


def _discover_plugins() -> dict[str, type["Analysis"]]:
    """Import all analysis modules and collect concrete Analysis subclasses.

    Returns
    -------
    dict[str, type[Analysis]]
        Mapping from canonical analysis name to Analysis subclass.

    Raises
    ------
    RuntimeError
        If two plugins register the same ``name``.
    """
    import polyzymd.analyses as analyses_pkg

    registry: dict[str, type[Analysis]] = {}

    # Import only direct plugin packages and simple single-file plugin modules
    package_path = analyses_pkg.__path__
    package_prefix = analyses_pkg.__name__ + "."

    for _, modname, is_pkg in pkgutil.iter_modules(package_path, prefix=package_prefix):
        del is_pkg
        # Skip infrastructure modules
        if _should_skip_module(modname, package_prefix):
            continue

        try:
            module = importlib.import_module(modname)
        except ImportError as exc:
            # Distinguish optional-dep failures (skip) from plugin bugs (re-raise)
            failing_module = getattr(exc, "name", None) or ""
            is_optional_dep = any(
                failing_module == dep or failing_module.startswith(dep + ".")
                for dep in _OPTIONAL_HEAVY_DEPS
            )
            if is_optional_dep:
                logger.info(
                    "Skipping analysis module %s: optional dependency %r not available",
                    modname,
                    failing_module,
                )
            else:
                logger.error(
                    "Failed to import analysis module %s: %s",
                    modname,
                    exc,
                    exc_info=True,
                )
                raise
            continue

        for attr_name in dir(module):
            try:
                obj = getattr(module, attr_name)
            except AttributeError:
                logger.debug(
                    "Could not access attribute %s.%s — skipping.",
                    modname,
                    attr_name,
                )
                continue  # Module __getattr__ raised; skip this attribute
            if not _is_concrete_analysis(obj):
                continue

            name = obj.name
            if not name or not name.strip():
                logger.warning(
                    "Analysis class %s.%s has empty name — skipping.",
                    obj.__module__,
                    obj.__qualname__,
                )
                continue
            name = name.strip()
            if name in registry:
                existing = registry[name]
                if existing is obj:
                    continue  # Same class found in multiple imports (sub-package re-export)
                raise RuntimeError(
                    f"Analysis name collision: both {existing.__module__}.{existing.__qualname__} "
                    f"and {obj.__module__}.{obj.__qualname__} use name={name!r}."
                )

            registry[name] = obj
            logger.debug(f"Discovered analysis plugin: {name} ({obj.__qualname__})")

        # A module may define only the plugin class; wrap it as a study module would be
        for plugin in _analyses_defined_in(module):
            if _is_concrete_analysis(plugin):
                continue
            from polyzymd.analyses.contract import contract_analysis

            wrapped = contract_analysis(plugin)
            if wrapped.name in registry:
                raise RuntimeError(
                    f"Analysis name collision: {modname}.{plugin.__qualname__} and "
                    f"{_source_module(registry[wrapped.name])} use name={wrapped.name!r}."
                )
            registry[wrapped.name] = wrapped

    return registry


@lru_cache(maxsize=1)
def _cached_registry() -> dict[str, type["Analysis"]]:
    """Return canonical analysis registry with caching.

    The cache is invalidated only by :func:`clear_cache`.
    """
    return _discover_plugins()


#: Analyses defined outside the package, by name.
_registered: dict[str, type["Analysis"]] = {}

#: Analysis folders already imported, with the names each one registered.
_loaded_directories: dict[Path, list[str]] = {}

#: Entries of an analysis folder that are never imported as analyses.
_SKIP_ENTRY_PREFIXES = ("_", ".", "test_")
_SKIP_ENTRY_NAMES = frozenset({"tests", "conftest.py"})


def clear_cache() -> None:
    """Forget discovered and registered analyses.  Useful in tests."""
    _cached_registry.cache_clear()
    _registered.clear()
    _loaded_directories.clear()


def register_analysis(analysis: Any) -> type["Analysis"]:
    """Make an analysis defined outside PolyzyMD available by name.

    Parameters
    ----------
    analysis : type[Analysis] or plugin class or plugin instance
        An ``Analysis`` subclass, or a plugin that satisfies
        :class:`~polyzymd.analyses.contract.AnalysisProtocol`, which is wrapped
        with :func:`~polyzymd.analyses.contract.contract_analysis`.

    Returns
    -------
    type[Analysis]
        The registered analysis class.

    Raises
    ------
    PluginContractError
        If the plugin breaks the contract, takes the name of a built-in
        analysis, or takes the name of a registered analysis defined in another
        module.
    """
    from polyzymd.analyses.contract import contract_analysis

    cls = analysis if _is_concrete_analysis(analysis) else contract_analysis(analysis)
    name = cls.name.strip()
    source = _source_module(cls)
    if name in _cached_registry():
        raise PluginContractError(
            f"{source} defines an analysis named {name!r}, which is a built-in PolyzyMD "
            "analysis. Give it another name so a comparison that names it cannot mean "
            "two different pieces of code."
        )
    existing = _registered.get(name)
    if existing is not None and _source_module(existing) != source:
        raise PluginContractError(
            f"Analysis name collision: {_source_module(existing)} and {source} both "
            f"define {name!r}."
        )
    _registered[name] = cls
    return cls


def load_analysis_directory(directory: Path | str) -> list[str]:
    """Import every analysis module in a folder and register what it defines.

    Each ``*.py`` file and each package in the folder is imported, except
    names starting with ``_``, ``.`` or ``test_`` and a ``tests`` folder. Those
    can hold helpers and tests; a module in the folder imports a helper with a
    relative import such as ``from ._common import f``. A folder is imported
    once per process; call :func:`clear_cache` to import it again after
    editing it.

    Parameters
    ----------
    directory : Path or str
        Folder to import, usually a study's ``analyses/``.

    Returns
    -------
    list[str]
        Names of the analyses the folder defines.

    Raises
    ------
    PluginContractError
        If a module fails to import or defines an analysis that cannot be
        registered.
    """
    root = Path(directory).resolve()
    if root in _loaded_directories:
        return list(_loaded_directories[root])
    package = _folder_package(root)
    names: list[str] = []
    for entry in sorted(root.iterdir()):
        if entry.name.startswith(_SKIP_ENTRY_PREFIXES) or entry.name in _SKIP_ENTRY_NAMES:
            continue
        if entry.is_file() and entry.suffix == ".py":
            module = _import_file(f"{package}.{entry.stem}", entry, package_dir=None)
        elif entry.is_dir() and (entry / "__init__.py").is_file():
            module = _import_file(
                f"{package}.{entry.name}", entry / "__init__.py", package_dir=entry
            )
        else:
            continue
        names += [register_analysis(item).name for item in _analyses_defined_in(module)]
    _loaded_directories[root] = names
    return list(names)


def _folder_package(root: Path) -> str:
    """Create the package an analysis folder's modules are imported under.

    The name carries a hash of the folder path, so two studies loaded in one
    process never share module names.
    """
    package = "polyzymd_study_" + hashlib.sha256(str(root).encode()).hexdigest()[:10]
    if package not in sys.modules:
        module = types.ModuleType(package)
        module.__path__ = [str(root)]
        sys.modules[package] = module
    return package


def _import_file(modname: str, path: Path, *, package_dir: Path | None) -> types.ModuleType:
    """Import one file of an analysis folder under ``modname``."""
    search = [str(package_dir)] if package_dir is not None else None
    spec = importlib.util.spec_from_file_location(modname, path, submodule_search_locations=search)
    if spec is None or spec.loader is None:
        raise PluginContractError(f"Cannot import analysis file {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[modname] = module
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        del sys.modules[modname]
        raise PluginContractError(
            f"Importing analysis file {path} failed: {type(exc).__name__}: {exc}"
        ) from exc
    return module


def _analyses_defined_in(module: types.ModuleType) -> list[Any]:
    """Analyses a module defines itself, wrapped or as bare plugin classes.

    Classes the module only imports are left out, so importing a built-in
    analysis to reuse it does not register it a second time.
    """
    found: list[Any] = []
    wrapped: set[type] = set()
    for obj in vars(module).values():
        if _is_concrete_analysis(obj) and _source_module(obj) == module.__name__:
            found.append(obj)
            plugin = getattr(obj, "plugin", None)
            if plugin is not None:
                wrapped.add(type(plugin))
    for obj in vars(module).values():
        if (
            inspect.isclass(obj)
            and obj.__module__ == module.__name__
            and obj not in wrapped
            and not _is_concrete_analysis(obj)
            and isinstance(getattr(obj, "name", None), str)
            and inspect.isclass(getattr(obj, "Settings", None))
            and callable(getattr(obj, "compute", None))
        ):
            found.append(obj)
    return found


def _source_module(cls: type) -> str:
    """Module that defines an analysis, looking through the contract wrapper."""
    plugin = getattr(cls, "plugin", None)
    return type(plugin).__module__ if plugin is not None else cls.__module__


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def get_analysis(name: str) -> type["Analysis"]:
    """Look up an Analysis class by canonical name.

    Parameters
    ----------
    name : str
        Canonical analysis name, for example ``"rmsf"``.

    Returns
    -------
    type[Analysis]
        The concrete Analysis subclass.

    Raises
    ------
    KeyError
        If no analysis matches *name*.
    """
    registry = {**_cached_registry(), **_registered}
    if name in registry:
        return registry[name]

    available = sorted(registry.keys())
    raise KeyError(f"Unknown analysis {name!r}.  Available: {', '.join(available)}")


def list_analyses() -> dict[str, type["Analysis"]]:
    """Return all discovered analyses.

    Returns
    -------
    dict[str, type[Analysis]]
        Mapping ``canonical_name -> Analysis subclass``, sorted by name.
    """
    registry = {**_cached_registry(), **_registered}
    return dict(sorted(registry.items()))


def list_all_names() -> list[str]:
    """Return all canonical analysis names, sorted.

    Returns
    -------
    list[str]
        All canonical names.
    """
    registry = {**_cached_registry(), **_registered}
    return sorted(registry.keys())
