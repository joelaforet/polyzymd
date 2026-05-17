"""Subclass contract validation for the analysis facade."""

from __future__ import annotations

import warnings
from typing import Any

from pydantic import BaseModel

_LEGACY_LIFECYCLE_HOOKS = (
    "run_replicate",
    "aggregate",
    "compare",
    "plot",
    "filter_conditions",
)

_REMOVED_RUNNER_HOOKS = (
    "build_runner",
    "summarize_replicate",
    "_run_replicate_via_runner",
)

_BUILTIN_PLUGIN_MODULE_PREFIXES = (
    "polyzymd.analyses.catalytic_triad",
    "polyzymd.analyses.contacts",
    "polyzymd.analyses.distances",
    "polyzymd.analyses.hydrogen_bonds",
    "polyzymd.analyses.rg",
    "polyzymd.analyses.rmsd",
    "polyzymd.analyses.rmsf",
    "polyzymd.analyses.sasa",
    "polyzymd.analyses.secondary_structure",
)


def validate_analysis_subclass(cls: type, *, base_cls: type, kwargs: dict[str, Any]) -> None:
    """Validate an analysis subclass after ABC initialization.

    Parameters
    ----------
    cls : type
        Subclass being initialized.
    base_cls : type
        Public facade base class.
    kwargs : dict[str, Any]
        Original subclass keyword arguments, retained for compatibility with
        the previous validation seam.
    """
    del kwargs
    if cls is base_cls:
        return
    if getattr(cls, "__abstractmethods__", False) or any(
        getattr(getattr(cls, name, None), "__isabstractmethod__", False) for name in dir(cls)
    ):
        return
    if not hasattr(cls, "name") or not isinstance(cls.name, str):
        raise TypeError(f"Analysis subclass {cls.__name__} must define 'name' as a ClassVar[str].")
    if not hasattr(cls, "Settings"):
        raise TypeError(
            f"Analysis subclass {cls.__name__} must define 'Settings' as a ClassVar[type]."
        )

    settings_cls = cls.Settings
    if not (isinstance(settings_cls, type) and issubclass(settings_cls, BaseModel)):
        raise TypeError(
            f"Analysis subclass {cls.__name__}.Settings must be a "
            f"pydantic BaseModel subclass, got {settings_cls!r}."
        )
    if not cls.has_compute_stage and cls.has_aggregate_stage:
        raise TypeError(
            f"Analysis subclass {cls.__name__} cannot set has_aggregate_stage=True "
            "when has_compute_stage=False."
        )

    removed_hooks = [name for name in _REMOVED_RUNNER_HOOKS if name in cls.__dict__]
    if removed_hooks:
        raise TypeError(
            f"Analysis subclass {cls.__name__} defines removed runner hook(s): "
            f"{', '.join(removed_hooks)}. Implement build_mda_jobs() for the "
            "MDAnalysis job lifecycle or override run_replicate() explicitly."
        )

    uses_run_replicate = cls.run_replicate is not base_cls.run_replicate
    uses_mda_jobs = cls.build_mda_jobs is not base_cls.build_mda_jobs

    if cls.has_compute_stage and not uses_run_replicate:
        if not uses_mda_jobs:
            raise TypeError(
                f"Analysis subclass {cls.__name__} public plugins must implement "
                "build_mda_jobs() or override run_replicate() when has_compute_stage=True; "
                "set has_compute_stage=False for compare-only plugins."
            )

    _warn_on_direct_lifecycle_overrides(cls)


def _warn_on_direct_lifecycle_overrides(cls: type) -> None:
    """Warn when a concrete external plugin defines lifecycle hooks directly.

    Parameters
    ----------
    cls : type
        Concrete analysis subclass being validated.
    """

    if _is_shipped_plugin_module(cls.__module__):
        return
    hooks = [name for name in _LEGACY_LIFECYCLE_HOOKS if name in cls.__dict__]
    if not hooks:
        return
    warnings.warn(
        f"{cls.__name__} directly defines analysis lifecycle hook(s): "
        f"{', '.join(hooks)}. Direct lifecycle overrides remain supported for "
        "advanced/internal integrations; prefer build_mda_jobs() for trajectory-native "
        "compute-stage plugins.",
        DeprecationWarning,
        stacklevel=3,
    )


def _is_shipped_plugin_module(module_name: str) -> bool:
    """Return whether a module belongs to an existing shipped plugin.

    Parameters
    ----------
    module_name : str
        Module path defining an analysis subclass.

    Returns
    -------
    bool
        ``True`` for built-in plugin modules that should not emit import-time
        deprecation warnings during the transition.
    """

    return any(
        module_name == prefix or module_name.startswith(prefix + ".")
        for prefix in _BUILTIN_PLUGIN_MODULE_PREFIXES
    )
