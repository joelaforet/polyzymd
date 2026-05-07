"""Automatic discovery of Analysis plugins via ``pkgutil``.

Scans ``src/polyzymd/analyses/`` for modules and sub-packages, imports
them, and collects all concrete :class:`Analysis` subclasses.  No
bootstrap files, no ``__init__.py`` edits, no decorators needed.

How Discovery Works
-------------------
1. ``pkgutil.iter_modules()`` yields every ``.py`` file and sub-package
   inside ``polyzymd.analyses``.
2. Each module is imported via ``importlib.import_module()``.
3. All module-level names are inspected; concrete subclasses of
   :class:`~polyzymd.analyses.base.Analysis` are collected.
4. Name collisions (two plugins with the same ``name``) raise immediately.

Contributor Impact
------------------
To add a new analysis, create a file in ``src/polyzymd/analyses/`` (or a
sub-package with ``__init__.py``), define a class inheriting from
``Analysis``, and set ``name`` as a ``ClassVar[str]``.  That's it.
"""

from __future__ import annotations

import importlib
import inspect
import logging
import pkgutil
from functools import lru_cache
from typing import TYPE_CHECKING

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
        "_results_base",
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


def _discover_plugins() -> tuple[dict[str, type["Analysis"]], dict[str, str]]:
    """Import all analysis modules and collect concrete Analysis subclasses.

    Returns
    -------
    tuple[dict[str, type[Analysis]], dict[str, str]]
        ``(name_registry, alias_to_name)`` mappings.

    Raises
    ------
    RuntimeError
        If two plugins register the same ``name`` or alias.
    """
    import polyzymd.analyses as analyses_pkg

    registry: dict[str, type[Analysis]] = {}
    alias_registry: dict[str, str] = {}  # alias -> canonical name

    # Walk all modules in the analyses package (one level deep for files,
    # recurse into sub-packages)
    package_path = analyses_pkg.__path__
    package_prefix = analyses_pkg.__name__ + "."

    for _, modname, _ in pkgutil.walk_packages(package_path, prefix=package_prefix):
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
            if name in alias_registry:
                raise RuntimeError(
                    f"Analysis name {name!r} (from {obj.__module__}.{obj.__qualname__}) "
                    f"conflicts with existing alias for {alias_registry[name]!r}."
                )
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

            # Register aliases
            for alias in getattr(obj, "aliases", ()):
                if not alias or not alias.strip():
                    logger.warning(
                        "Analysis %r has empty alias — skipping.",
                        name,
                    )
                    continue
                alias = alias.strip()
                if alias in alias_registry:
                    raise RuntimeError(
                        f"Analysis alias collision: {alias!r} is claimed by both "
                        f"{alias_registry[alias]!r} and {name!r}."
                    )
                if alias in registry:
                    raise RuntimeError(
                        f"Analysis alias {alias!r} (from {name!r}) conflicts with "
                        f"existing analysis name {alias!r}."
                    )
                alias_registry[alias] = name

    return registry, alias_registry


@lru_cache(maxsize=1)
def _cached_registry() -> tuple[dict[str, type["Analysis"]], dict[str, str]]:
    """Return (name_registry, alias_to_name) with caching.

    The cache is invalidated only by :func:`clear_cache`.
    """
    return _discover_plugins()


def clear_cache() -> None:
    """Clear the discovery cache.  Useful in tests."""
    _cached_registry.cache_clear()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def get_analysis(name: str) -> type["Analysis"]:
    """Look up an Analysis class by name or alias.

    Parameters
    ----------
    name : str
        Analysis name (e.g. ``"rmsf"``) or alias (e.g. ``"triad"``).

    Returns
    -------
    type[Analysis]
        The concrete Analysis subclass.

    Raises
    ------
    KeyError
        If no analysis matches *name*.
    """
    registry, aliases = _cached_registry()
    if name in registry:
        return registry[name]
    canonical = aliases.get(name)
    if canonical is not None:
        return registry[canonical]
    from polyzymd.core.archived_features import (
        format_archived_analysis_message,
        get_archived_analysis_plugin,
    )

    if get_archived_analysis_plugin(name) is not None:
        raise KeyError(format_archived_analysis_message(name, context="analysis lookup"))

    available = sorted(set(registry.keys()) | set(aliases.keys()))
    raise KeyError(f"Unknown analysis {name!r}.  Available: {', '.join(available)}")


def list_analyses() -> dict[str, type["Analysis"]]:
    """Return all discovered analyses.

    Returns
    -------
    dict[str, type[Analysis]]
        Mapping ``canonical_name -> Analysis subclass``, sorted by name.
    """
    registry, _ = _cached_registry()
    return dict(sorted(registry.items()))


def list_all_names() -> list[str]:
    """Return all names and aliases, sorted.

    Returns
    -------
    list[str]
        All canonical names and aliases.
    """
    registry, aliases = _cached_registry()
    return sorted(set(registry.keys()) | set(aliases.keys()))
