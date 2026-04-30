"""Subclass contract validation for the analysis facade."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel


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

    uses_run_replicate = cls.run_replicate is not base_cls.run_replicate
    uses_runner_build = cls.build_runner is not base_cls.build_runner
    uses_runner_summary = cls.summarize_replicate is not base_cls.summarize_replicate

    if cls.has_compute_stage and not uses_run_replicate:
        if uses_runner_build != uses_runner_summary:
            raise TypeError(
                f"Analysis subclass {cls.__name__} must implement both build_runner() and "
                "summarize_replicate() when using the runner-based path."
            )
        if not uses_runner_build:
            raise TypeError(
                f"Analysis subclass {cls.__name__} public plugins must implement both "
                "build_runner() and summarize_replicate() when has_compute_stage=True; "
                "direct run_replicate() overrides are advanced/internal only."
            )
