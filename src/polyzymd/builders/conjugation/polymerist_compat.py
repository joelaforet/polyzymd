"""Compatibility helpers for importing Polymerist on Python 3.12."""

from __future__ import annotations

import importlib
import sys
from types import ModuleType
from typing import Any


def ensure_polymerist_py312_compat() -> bool:
    """Install the temporary Polymerist Python 3.12 compatibility shim.

    Polymerist currently imports the private
    :mod:`importlib.resources._common.get_package` helper. Python 3.12 removed
    that helper, so importing Polymerist can fail before PolyzyMD reaches any
    conjugation-specific code. This function restores the small Python 3.11
    behavior Polymerist expects by injecting a local ``get_package`` function
    only when the helper is missing.

    The shim is temporary and should be removed after Polymerist switches to a
    public :mod:`importlib.resources` API. It intentionally does not import
    Polymerist itself.

    Returns
    -------
    bool
        ``True`` when a compatibility function was injected, otherwise
        ``False``.
    """
    import importlib.resources._common as resources_common

    if hasattr(resources_common, "get_package"):
        return False

    resources_common.get_package = _get_package_compat  # type: ignore[attr-defined]
    return True


def import_polymerist_building() -> ModuleType:
    """Import :mod:`polymerist.polymers.building` behind the compatibility shim.

    Returns
    -------
    types.ModuleType
        Imported Polymerist building module.

    Raises
    ------
    ImportError
        If Polymerist or one of its runtime dependencies is unavailable.
    """
    ensure_polymerist_py312_compat()
    return importlib.import_module("polymerist.polymers.building")


def polymerist_py312_compat_status() -> dict[str, Any]:
    """Return diagnostic details for the temporary Polymerist shim.

    Returns
    -------
    dict[str, Any]
        Compatibility status suitable for inclusion in conjugation diagnostics.
    """
    import importlib.resources._common as resources_common

    missing_get_package = not hasattr(resources_common, "get_package")
    relevant = sys.version_info >= (3, 12)
    return {
        "relevant": relevant,
        "python_version": sys.version.split()[0],
        "get_package_present": not missing_get_package,
        "shim_required": relevant and missing_get_package,
        "rationale": (
            "Polymerist imports private importlib.resources._common.get_package; "
            "PolyzyMD injects a Python 3.11-compatible helper on Python 3.12 when needed"
        ),
    }


def _get_package_compat(package: str | ModuleType) -> ModuleType:
    """Resolve and validate a package using Python 3.11-compatible behavior.

    Parameters
    ----------
    package : str or types.ModuleType
        Package name or imported package module.

    Returns
    -------
    types.ModuleType
        Imported package module.

    Raises
    ------
    TypeError
        If the resolved object is not a package.
    """
    module = importlib.import_module(package) if isinstance(package, str) else package
    spec = getattr(module, "__spec__", None)

    if spec is None or not _spec_has_package_locations(spec):
        raise TypeError(f"{module!r} is not a package")

    return module


def _spec_has_package_locations(spec: Any) -> bool:
    """Return whether an import spec represents a package."""
    if getattr(spec, "submodule_search_locations", None) is not None:
        return True

    wrapped_spec = getattr(spec, "_unwrapped", None)
    if wrapped_spec is not None:
        return getattr(wrapped_spec, "submodule_search_locations", None) is not None

    loader_spec = getattr(getattr(spec, "loader", None), "spec", None)
    if loader_spec is not None:
        return getattr(loader_spec, "submodule_search_locations", None) is not None

    return False
