"""Registry for comparator types.

This module provides extensible infrastructure for registering comparator types
following the Open-Closed Principle (OCP). New comparators can be added by
registering with the ComparatorRegistry without modifying core code.

Example
-------
Registering a new comparator:

>>> from polyzymd.compare.core.registry import ComparatorRegistry
>>> from polyzymd.compare.core.base import BaseComparator
>>>
>>> @ComparatorRegistry.register("my_metric")
... class MyComparator(BaseComparator):
...     @classmethod
...     def comparison_type_name(cls) -> str:
...         return "my_metric"
...     ...
>>>
>>> # Create comparator instance via registry
>>> comparator = ComparatorRegistry.create("my_metric", config, settings)
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Type

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig
    from polyzymd.compare.core.base import BaseComparator

logger = logging.getLogger("polyzymd.compare")


class ComparatorRegistry:
    """Registry for comparator implementations.

    Allows new comparators to be registered without modifying core code.
    Use the `register` decorator to add new comparator classes.

    Examples
    --------
    >>> @ComparatorRegistry.register("rmsf")
    ... class RMSFComparator(BaseComparator):
    ...     ...
    >>>
    >>> # List available comparators
    >>> ComparatorRegistry.list_available()
    ['contacts', 'rmsf', 'triad']
    >>>
    >>> # Create comparator instance
    >>> comparator = ComparatorRegistry.create("rmsf", config, settings)
    """

    _registry: dict[str, Type["BaseComparator"]] = {}

    @classmethod
    def register(cls, name: str | None = None):
        """Decorator to register a comparator class.

        Parameters
        ----------
        name : str, optional
            Registry key. If None, uses the class's comparison_type_name().

        Returns
        -------
        Callable
            Decorator function.

        Examples
        --------
        >>> @ComparatorRegistry.register("rmsf")
        ... class RMSFComparator(BaseComparator):
        ...     @classmethod
        ...     def comparison_type_name(cls) -> str:
        ...         return "rmsf"
        """

        def decorator(comparator_class: Type["BaseComparator"]):
            key = name if name is not None else comparator_class.comparison_type_name()
            key = key.lower()
            if key in cls._registry:
                logger.warning(f"Overwriting existing comparator registration: {key}")
            cls._registry[key] = comparator_class
            logger.debug(f"Registered comparator: {key}")
            return comparator_class

        return decorator

    @classmethod
    def get(cls, name: str) -> Type["BaseComparator"]:
        """Get comparator class by name.

        Parameters
        ----------
        name : str
            Comparator type identifier.

        Returns
        -------
        Type[BaseComparator]
            The registered comparator class.

        Raises
        ------
        ValueError
            If the comparator type is not registered.
        """
        key = name.lower()
        if key not in cls._registry:
            available = ", ".join(sorted(cls._registry.keys()))
            raise ValueError(f"Unknown comparator type: '{name}'. Available: {available}")
        return cls._registry[key]

    @classmethod
    def list_available(cls) -> list[str]:
        """List all registered comparator types.

        Returns
        -------
        list[str]
            Sorted list of registered type names.
        """
        return sorted(cls._registry.keys())

    @classmethod
    def is_registered(cls, name: str) -> bool:
        """Check if a comparator type is registered.

        Parameters
        ----------
        name : str
            Comparator type identifier.

        Returns
        -------
        bool
            True if registered, False otherwise.
        """
        return name.lower() in cls._registry

    @classmethod
    def create(
        cls,
        name: str,
        config: "ComparisonConfig",
        analysis_settings,
        equilibration: str | None = None,
        **kwargs,
    ) -> "BaseComparator":
        """Factory to create a comparator instance.

        Parameters
        ----------
        name : str
            Comparator type identifier.
        config : ComparisonConfig
            Comparison configuration.
        analysis_settings
            Analysis-specific settings.
        equilibration : str, optional
            Equilibration time override.
        **kwargs
            Additional comparator-specific arguments.

        Returns
        -------
        BaseComparator
            Configured comparator instance.
        """
        comparator_class = cls.get(name)
        return comparator_class(
            config=config,
            analysis_settings=analysis_settings,
            equilibration=equilibration,
            **kwargs,
        )

    @classmethod
    def get_for_analysis_type(cls, settings_key: str) -> Type["BaseComparator"]:
        """Look up comparator by analysis_settings YAML key.

        First checks if ``settings_key`` is a direct registry key.
        Then checks all registered comparators for a matching
        ``config_settings_key()``.

        Parameters
        ----------
        settings_key : str
            Analysis settings key (e.g., ``"catalytic_triad"``).

        Returns
        -------
        Type[BaseComparator]
            The registered comparator class.

        Raises
        ------
        ValueError
            If no comparator matches the given settings key.
        """
        key = settings_key.lower()
        if key in cls._registry:
            return cls._registry[key]

        for comp_cls in cls._registry.values():
            if comp_cls.config_settings_key() == key:
                return comp_cls

        available = ", ".join(sorted(cls._registry.keys()))
        raise ValueError(
            f"No comparator for analysis settings key '{settings_key}'. Available: {available}"
        )

    @classmethod
    def create_from_config(
        cls,
        settings_key: str,
        config: "ComparisonConfig",
        equilibration: str | None = None,
        **kwargs,
    ) -> "BaseComparator":
        """Config-driven factory: resolve comparator and build from config.

        Automatically:
        1. Resolves the comparator class from the settings key
        2. Extracts analysis_settings from config
        3. Injects comparison_settings if the comparator opts in

        Parameters
        ----------
        settings_key : str
            Analysis settings YAML key (e.g., ``"rmsf"``,
            ``"catalytic_triad"``).
        config : ComparisonConfig
            Comparison configuration.
        equilibration : str, optional
            Equilibration time override.
        **kwargs
            Additional keyword arguments forwarded to the comparator
            constructor (e.g., selection_override for RMSF).

        Returns
        -------
        BaseComparator
            Configured comparator instance.

        Raises
        ------
        ValueError
            If no comparator is registered for the settings key, or
            the analysis_settings section is missing from config.
        """
        comparator_cls = cls.get_for_analysis_type(settings_key)

        analysis_settings = config.analysis_settings.get(comparator_cls.config_settings_key())
        if analysis_settings is None:
            raise ValueError(
                f"No '{comparator_cls.config_settings_key()}' in "
                "analysis_settings section of config."
            )

        init_kwargs: dict[str, object] = {
            "config": config,
            "analysis_settings": analysis_settings,
            "equilibration": equilibration,
            **kwargs,
        }

        if comparator_cls.uses_comparison_settings:
            comparison_settings = config.comparison_settings.get(
                comparator_cls.config_settings_key()
            )
            init_kwargs.setdefault("comparison_settings", comparison_settings)

        return comparator_cls(**init_kwargs)

    @classmethod
    def clear(cls) -> None:
        """Clear the registry (for testing purposes)."""
        cls._registry.clear()
