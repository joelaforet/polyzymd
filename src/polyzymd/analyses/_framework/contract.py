"""Subclass contract validation for the analysis facade."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel

_REMOVED_RUNNER_HOOKS = (
    "build_runner",
    "summarize_replicate",
    "_run_replicate_via_runner",
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
        Original subclass keyword arguments supplied during class creation.
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
            "MDAnalysis job lifecycle with build_mda_collector(), or set "
            "has_compute_stage=False for compare-only plugins."
        )

    run_replicate_owner = _find_removed_run_replicate_owner(cls, base_cls=base_cls)
    if run_replicate_owner is cls:
        raise TypeError(
            f"Analysis subclass {cls.__name__} defines removed hook run_replicate(). "
            "Plugins must not override or inherit run_replicate(); implement "
            "build_mda_jobs() with build_mda_collector(), or set has_compute_stage=False "
            "for compare-only plugins."
        )
    if run_replicate_owner is not None:
        raise TypeError(
            f"Analysis subclass {cls.__name__} inherits removed hook run_replicate() "
            f"from {run_replicate_owner.__name__}. Plugins must not override or inherit "
            "run_replicate(); implement build_mda_jobs() with build_mda_collector(), or set "
            "has_compute_stage=False for compare-only plugins."
        )

    uses_mda_jobs = cls.build_mda_jobs is not base_cls.build_mda_jobs

    if cls.has_compute_stage and not uses_mda_jobs:
        raise TypeError(
            f"Analysis subclass {cls.__name__} public plugins must implement "
            "build_mda_jobs() when has_compute_stage=True; set has_compute_stage=False "
            "for compare-only plugins."
        )


def _find_removed_run_replicate_owner(cls: type, *, base_cls: type) -> type | None:
    """Return the MRO class that provides the removed ``run_replicate`` hook.

    Parameters
    ----------
    cls : type
        Subclass being validated.
    base_cls : type
        Public analysis base class whose own definitions are framework-owned.

    Returns
    -------
    type | None
        Class in the effective MRO that defines ``run_replicate``, excluding
        the public analysis base class and ``object``.
    """

    for mro_cls in cls.__mro__:
        if mro_cls in {base_cls, object}:
            continue
        if "run_replicate" in mro_cls.__dict__:
            return mro_cls
    return None
