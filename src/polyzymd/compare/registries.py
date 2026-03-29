"""Registry pattern for plot settings.

This module provides the plot settings registry and base class used by the
``compare`` workflow to discover per-analysis plot customization at parse time.

Registries
----------
- ``PlotSettingsRegistry`` — Per-analysis plot customisation

Base classes
------------
- ``BasePlotSettings``  — Base class for per-analysis plot settings

Note
----
Analysis plugins are discovered via ``analyses.discovery`` and their
``Settings`` inner classes provide the authoritative configuration.
"""

from __future__ import annotations

import logging
from typing import Type

from pydantic import BaseModel

logger = logging.getLogger(__name__)


# ============================================================================
# Plot Settings Registry (active — used by compare/config.py)
# ============================================================================


class BasePlotSettings(BaseModel):
    """Abstract base class for per-analysis plot settings.

    Each analysis type that has plot customization should subclass this
    and register with ``PlotSettingsRegistry``.  The class is intentionally
    minimal — it exists only so the registry can enforce a common type.

    Examples
    --------
    >>> @PlotSettingsRegistry.register("rmsf")
    ... class RMSFPlotSettings(BasePlotSettings):
    ...     show_error: bool = True
    ...     highlight_residues: list[int] = Field(default_factory=list)
    """


class PlotSettingsRegistry:
    """Registry for per-analysis plot settings types.

    Allows new analysis types to register their plot-customization models
    without modifying the central ``PlotSettings`` container. At parse time,
    ``PlotSettings`` discovers registered types and constructs them from
    the YAML dict; at access time, ``PlotSettings.__getattr__`` returns a
    default-constructed instance for any registered type that was not
    explicitly configured, so ``self.settings.rmsf.figsize_profile`` always
    works even if the user omitted the ``rmsf:`` block.

    Examples
    --------
    >>> @PlotSettingsRegistry.register("rmsf")
    ... class RMSFPlotSettings(BasePlotSettings):
    ...     show_error: bool = True
    >>>
    >>> cls = PlotSettingsRegistry.get("rmsf")
    >>> PlotSettingsRegistry.list_available()
    ['rmsf']
    """

    _registry: dict[str, Type[BasePlotSettings]] = {}

    @classmethod
    def register(cls, name: str):
        """Decorator to register a plot settings class.

        Parameters
        ----------
        name : str
            Registry key (must match the analysis type identifier,
            e.g. ``"rmsf"``, ``"contacts"``).

        Returns
        -------
        Callable
            Decorator function.
        """

        def decorator(settings_class: Type[BasePlotSettings]):
            key = name.lower()
            if key in cls._registry:
                logger.warning(f"Overwriting existing plot settings registration: {key}")
            cls._registry[key] = settings_class
            logger.debug(f"Registered plot settings: {key}")
            return settings_class

        return decorator

    @classmethod
    def get(cls, name: str) -> Type[BasePlotSettings]:
        """Get plot settings class by name.

        Parameters
        ----------
        name : str
            Analysis type identifier.

        Returns
        -------
        Type[BasePlotSettings]
            The registered plot settings class.

        Raises
        ------
        ValueError
            If the analysis type is not registered.
        """
        key = name.lower()
        if key not in cls._registry:
            available = ", ".join(sorted(cls._registry.keys()))
            raise ValueError(f"Unknown plot settings type: '{name}'. Available: {available}")
        return cls._registry[key]

    @classmethod
    def list_available(cls) -> list[str]:
        """List all registered plot settings types.

        Returns
        -------
        list[str]
            Sorted list of registered type names.
        """
        return sorted(cls._registry.keys())

    @classmethod
    def is_registered(cls, name: str) -> bool:
        """Check if a plot settings type is registered.

        Parameters
        ----------
        name : str
            Analysis type identifier.

        Returns
        -------
        bool
            True if registered, False otherwise.
        """
        return name.lower() in cls._registry

    @classmethod
    def clear(cls) -> None:
        """Clear the registry (for testing purposes)."""
        cls._registry.clear()
