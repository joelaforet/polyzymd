"""Registry pattern for analysis types.

This module provides extensible infrastructure for registering analysis types
following the Open-Closed Principle (OCP). New analysis types can be added
by registering with the appropriate registry without modifying core code.

Registries:
- AnalysisSettingsRegistry: For analysis parameter configs (WHAT to analyze)
- ComparisonSettingsRegistry: For comparison parameter configs (HOW to compare)
- AnalyzerRegistry: For analyzer implementations (future use)

Example
-------
Registering a new analysis settings type:

>>> from polyzymd.analysis.core.registry import (
...     BaseAnalysisSettings,
...     AnalysisSettingsRegistry,
... )
>>>
>>> @AnalysisSettingsRegistry.register("my_analysis")
... class MyAnalysisSettings(BaseAnalysisSettings):
...     my_param: str = "default"
...
...     @classmethod
...     def analysis_type(cls) -> str:
...         return "my_analysis"
...
...     def to_analysis_yaml_dict(self) -> dict:
...         return {"enabled": True, "my_param": self.my_param}
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Type

from pydantic import BaseModel

if TYPE_CHECKING:
    from collections.abc import Sequence

    from polyzymd.analysis.results.base import BaseAnalysisResult
    from polyzymd.config.schema import SimulationConfig

logger = logging.getLogger(__name__)


# ============================================================================
# Abstract Base Classes
# ============================================================================


class BaseAnalysisSettings(BaseModel, ABC):
    """Abstract base class for analysis settings configurations.

    Analysis settings define WHAT to analyze (parameters shared across
    conditions in a comparison). Subclasses must implement:
    - analysis_type(): Returns the unique identifier for this analysis
    - to_analysis_yaml_dict(): Converts to analysis.yaml format

    Attributes
    ----------
    None defined here; subclasses define their own attributes.

    Notes
    -----
    Unlike AnalysisConfig models which have `enabled: bool`, these settings
    models imply enabled by their presence (no `enabled` field).
    """

    @classmethod
    @abstractmethod
    def analysis_type(cls) -> str:
        """Return the unique identifier for this analysis type.

        Returns
        -------
        str
            Analysis type identifier (e.g., "rmsf", "contacts").
        """
        ...

    @abstractmethod
    def to_analysis_yaml_dict(self) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary.

        The returned dict should include `enabled: True` and all
        analysis-specific parameters.

        Returns
        -------
        dict[str, Any]
            Dictionary suitable for writing to analysis.yaml.
        """
        ...


class BaseComparisonSettings(BaseModel, ABC):
    """Abstract base class for comparison settings configurations.

    Comparison settings define HOW to compare across conditions
    (statistical parameters, thresholds, etc.). Some analysis types
    may have no comparison-specific parameters (empty class).

    Subclasses must implement:
    - analysis_type(): Returns the unique identifier (must match analysis settings)
    """

    @classmethod
    @abstractmethod
    def analysis_type(cls) -> str:
        """Return the unique identifier for this analysis type.

        Returns
        -------
        str
            Analysis type identifier (must match corresponding analysis settings).
        """
        ...


class BaseAnalyzer(ABC):
    """Abstract base class for all analyzers.

    Analyzers perform the actual computation on trajectories. They follow
    a consistent interface for initialization, single-replicate computation,
    and multi-replicate aggregation.

    This ABC is provided for future refactoring of existing analyzers
    (RMSFCalculator, ContactAnalyzer, etc.) to a unified interface.

    Attributes
    ----------
    None defined here; subclasses define their own attributes.

    Notes
    -----
    Current analyzers do not inherit from this class. This is provided
    for future extensibility and to define the expected interface.
    """

    @classmethod
    @abstractmethod
    def analysis_type(cls) -> str:
        """Return the unique identifier for this analyzer.

        Returns
        -------
        str
            Analysis type identifier (e.g., "rmsf", "contacts").
        """
        ...

    @classmethod
    @abstractmethod
    def from_config(
        cls,
        analysis_settings: BaseAnalysisSettings,
        sim_config: "SimulationConfig",
        equilibration: str = "0ns",
    ) -> "BaseAnalyzer":
        """Factory method to create analyzer from config.

        Parameters
        ----------
        analysis_settings : BaseAnalysisSettings
            Analysis-specific settings.
        sim_config : SimulationConfig
            Simulation configuration for trajectory loading.
        equilibration : str
            Equilibration time to skip.

        Returns
        -------
        BaseAnalyzer
            Configured analyzer instance.
        """
        ...

    @abstractmethod
    def compute(self, replicate: int, **kwargs) -> "BaseAnalysisResult":
        """Run analysis for a single replicate.

        Parameters
        ----------
        replicate : int
            Replicate number to analyze.
        **kwargs
            Additional analyzer-specific parameters.

        Returns
        -------
        BaseAnalysisResult
            Analysis result for this replicate.
        """
        ...

    @abstractmethod
    def compute_aggregated(self, replicates: "Sequence[int]", **kwargs) -> "BaseAnalysisResult":
        """Run aggregated analysis across multiple replicates.

        Parameters
        ----------
        replicates : Sequence[int]
            List of replicate numbers to analyze.
        **kwargs
            Additional analyzer-specific parameters.

        Returns
        -------
        BaseAnalysisResult
            Aggregated analysis result.
        """
        ...

    @property
    @abstractmethod
    def label(self) -> str:
        """Human-readable label for this analyzer.

        Returns
        -------
        str
            Display label for reports and logs.
        """
        ...


# ============================================================================
# Registry Classes
# ============================================================================


class AnalysisSettingsRegistry:
    """Registry for analysis settings types.

    Allows new analysis types to be registered without modifying core code.
    Use the `register` decorator to add new analysis settings classes.

    Examples
    --------
    >>> @AnalysisSettingsRegistry.register("my_analysis")
    ... class MyAnalysisSettings(BaseAnalysisSettings):
    ...     ...
    >>>
    >>> # Get registered class
    >>> cls = AnalysisSettingsRegistry.get("my_analysis")
    >>>
    >>> # List all available types
    >>> types = AnalysisSettingsRegistry.list_available()
    """

    _registry: dict[str, Type[BaseAnalysisSettings]] = {}

    @classmethod
    def register(cls, name: str | None = None):
        """Decorator to register an analysis settings class.

        Parameters
        ----------
        name : str, optional
            Registry key. If None, uses the class's analysis_type().

        Returns
        -------
        Callable
            Decorator function.

        Examples
        --------
        >>> @AnalysisSettingsRegistry.register("rmsf")
        ... class RMSFAnalysisSettings(BaseAnalysisSettings):
        ...     ...
        """

        def decorator(settings_class: Type[BaseAnalysisSettings]):
            key = name if name is not None else settings_class.analysis_type()
            key = key.lower()
            if key in cls._registry:
                logger.warning(f"Overwriting existing analysis settings registration: {key}")
            cls._registry[key] = settings_class
            logger.debug(f"Registered analysis settings: {key}")
            return settings_class

        return decorator

    @classmethod
    def get(cls, name: str) -> Type[BaseAnalysisSettings]:
        """Get analysis settings class by name.

        Parameters
        ----------
        name : str
            Analysis type identifier.

        Returns
        -------
        Type[BaseAnalysisSettings]
            The registered settings class.

        Raises
        ------
        ValueError
            If the analysis type is not registered.
        """
        key = name.lower()
        if key not in cls._registry:
            available = ", ".join(sorted(cls._registry.keys()))
            raise ValueError(f"Unknown analysis type: '{name}'. Available: {available}")
        return cls._registry[key]

    @classmethod
    def list_available(cls) -> list[str]:
        """List all registered analysis types.

        Returns
        -------
        list[str]
            Sorted list of registered analysis type names.
        """
        return sorted(cls._registry.keys())

    @classmethod
    def is_registered(cls, name: str) -> bool:
        """Check if an analysis type is registered.

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


class ComparisonSettingsRegistry:
    """Registry for comparison settings types.

    Allows new comparison settings to be registered without modifying core code.
    Use the `register` decorator to add new comparison settings classes.

    Examples
    --------
    >>> @ComparisonSettingsRegistry.register("contacts")
    ... class ContactsComparisonSettings(BaseComparisonSettings):
    ...     fdr_alpha: float = 0.05
    ...     ...
    """

    _registry: dict[str, Type[BaseComparisonSettings]] = {}

    @classmethod
    def register(cls, name: str | None = None):
        """Decorator to register a comparison settings class.

        Parameters
        ----------
        name : str, optional
            Registry key. If None, uses the class's analysis_type().

        Returns
        -------
        Callable
            Decorator function.
        """

        def decorator(settings_class: Type[BaseComparisonSettings]):
            key = name if name is not None else settings_class.analysis_type()
            key = key.lower()
            if key in cls._registry:
                logger.warning(f"Overwriting existing comparison settings registration: {key}")
            cls._registry[key] = settings_class
            logger.debug(f"Registered comparison settings: {key}")
            return settings_class

        return decorator

    @classmethod
    def get(cls, name: str) -> Type[BaseComparisonSettings]:
        """Get comparison settings class by name.

        Parameters
        ----------
        name : str
            Analysis type identifier.

        Returns
        -------
        Type[BaseComparisonSettings]
            The registered settings class.

        Raises
        ------
        ValueError
            If the analysis type is not registered.
        """
        key = name.lower()
        if key not in cls._registry:
            available = ", ".join(sorted(cls._registry.keys()))
            raise ValueError(f"Unknown comparison settings type: '{name}'. Available: {available}")
        return cls._registry[key]

    @classmethod
    def list_available(cls) -> list[str]:
        """List all registered comparison settings types.

        Returns
        -------
        list[str]
            Sorted list of registered type names.
        """
        return sorted(cls._registry.keys())

    @classmethod
    def is_registered(cls, name: str) -> bool:
        """Check if a comparison settings type is registered.

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


class AnalyzerRegistry:
    """Registry for analyzer implementations.

    Allows new analyzers to be registered without modifying core code.
    This registry is provided for future use when analyzers are refactored
    to inherit from BaseAnalyzer.

    Examples
    --------
    >>> @AnalyzerRegistry.register("rmsf")
    ... class RMSFCalculator(BaseAnalyzer):
    ...     ...
    >>>
    >>> # Create analyzer instance
    >>> analyzer = AnalyzerRegistry.create("rmsf", settings, sim_config)
    """

    _registry: dict[str, Type[BaseAnalyzer]] = {}

    @classmethod
    def register(cls, name: str | None = None):
        """Decorator to register an analyzer class.

        Parameters
        ----------
        name : str, optional
            Registry key. If None, uses the class's analysis_type().

        Returns
        -------
        Callable
            Decorator function.
        """

        def decorator(analyzer_class: Type[BaseAnalyzer]):
            key = name if name is not None else analyzer_class.analysis_type()
            key = key.lower()
            if key in cls._registry:
                logger.warning(f"Overwriting existing analyzer registration: {key}")
            cls._registry[key] = analyzer_class
            logger.debug(f"Registered analyzer: {key}")
            return analyzer_class

        return decorator

    @classmethod
    def get(cls, name: str) -> Type[BaseAnalyzer]:
        """Get analyzer class by name.

        Parameters
        ----------
        name : str
            Analysis type identifier.

        Returns
        -------
        Type[BaseAnalyzer]
            The registered analyzer class.

        Raises
        ------
        ValueError
            If the analysis type is not registered.
        """
        key = name.lower()
        if key not in cls._registry:
            available = ", ".join(sorted(cls._registry.keys()))
            raise ValueError(f"Unknown analyzer type: '{name}'. Available: {available}")
        return cls._registry[key]

    @classmethod
    def list_available(cls) -> list[str]:
        """List all registered analyzer types.

        Returns
        -------
        list[str]
            Sorted list of registered type names.
        """
        return sorted(cls._registry.keys())

    @classmethod
    def is_registered(cls, name: str) -> bool:
        """Check if an analyzer type is registered.

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
    def create(
        cls,
        name: str,
        analysis_settings: BaseAnalysisSettings,
        sim_config: "SimulationConfig",
        equilibration: str = "0ns",
    ) -> BaseAnalyzer:
        """Factory to create an analyzer instance.

        Parameters
        ----------
        name : str
            Analysis type identifier.
        analysis_settings : BaseAnalysisSettings
            Analysis-specific settings.
        sim_config : SimulationConfig
            Simulation configuration for trajectory loading.
        equilibration : str
            Equilibration time to skip.

        Returns
        -------
        BaseAnalyzer
            Configured analyzer instance.
        """
        analyzer_class = cls.get(name)
        return analyzer_class.from_config(analysis_settings, sim_config, equilibration)

    @classmethod
    def clear(cls) -> None:
        """Clear the registry (for testing purposes)."""
        cls._registry.clear()


# ============================================================================
# Plot Settings Registry
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
