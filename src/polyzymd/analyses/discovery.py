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
        "base",
        "stats",
        "discovery",
        "orchestrator",
        "runner",
        "config",
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


def _discover_plugins() -> dict[str, type["Analysis"]]:
    """Import all analysis modules and collect concrete Analysis subclasses.

    Returns
    -------
    dict[str, type[Analysis]]
        Mapping ``analysis_name -> Analysis subclass``.

    Raises
    ------
    RuntimeError
        If two plugins register the same ``name``.
    """
    import polyzymd.analyses as analyses_pkg

    registry: dict[str, type[Analysis]] = {}
    alias_registry: dict[str, str] = {}  # alias -> canonical name

    # Walk all modules in the analyses package (one level deep for files,
    # recurse into sub-packages)
    package_path = analyses_pkg.__path__
    package_prefix = analyses_pkg.__name__ + "."

    for importer, modname, ispkg in pkgutil.walk_packages(package_path, prefix=package_prefix):
        # Skip infrastructure modules
        short_name = modname.rsplit(".", 1)[-1]
        if short_name.startswith("_") or short_name in _SKIP_MODULES:
            continue

        try:
            module = importlib.import_module(modname)
        except Exception:
            logger.warning(f"Failed to import analysis module {modname}", exc_info=True)
            continue

        for attr_name in dir(module):
            obj = getattr(module, attr_name)
            if not _is_concrete_analysis(obj):
                continue

            name = obj.name
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

    return registry


@lru_cache(maxsize=1)
def _cached_registry() -> tuple[dict[str, type["Analysis"]], dict[str, str]]:
    """Return (name_registry, alias_to_name) with caching.

    The cache is invalidated only by :func:`clear_cache`.
    """
    registry = _discover_plugins()

    # Build alias map
    alias_to_name: dict[str, str] = {}
    for name, cls in registry.items():
        for alias in getattr(cls, "aliases", ()):
            alias_to_name[alias] = name

    return registry, alias_to_name


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
