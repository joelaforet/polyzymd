"""Registry pattern for analysis and comparison settings.

This module is the authoritative home for the base classes and registries
used by the ``compare`` workflow to discover analysis/comparison/plot
settings types at parse time.

Registries
----------
- ``AnalysisSettingsRegistry``  — WHAT to analyze (parameter configs)
- ``ComparisonSettingsRegistry`` — HOW to compare (statistical thresholds)
- ``AnalyzerRegistry``          — Analyzer implementations (legacy, will be removed)
- ``PlotSettingsRegistry``      — Per-analysis plot customisation

Base classes
------------
- ``BaseAnalysisSettings``     — ABC for analysis parameter configs
- ``BaseComparisonSettings``   — ABC for comparison parameter configs
- ``BaseAnalyzer``             — ABC for analyzer implementations (legacy)
- ``BasePlotSettings``         — Base class for per-analysis plot settings

Note
----
Previously these lived in ``polyzymd.analysis.core.registry``.  They were
relocated here as part of the OCP-compliance refactor (Phase 3) so that
the ``compare`` package — the sole consumer — owns the infrastructure and
the dependency direction flows ``compare → analyses`` rather than
``compare → analysis``.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Type

from pydantic import BaseModel

if TYPE_CHECKING:
    from collections.abc import Sequence

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
    Unlike AnalysisConfig models which have ``enabled: bool``, these settings
    models imply enabled by their presence (no ``enabled`` field).
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

        The returned dict should include ``enabled: True`` and all
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
    """Abstract base class defining the expected interface for analyzers.

    Analyzers perform computation on MD trajectories, producing structured
    results.

    .. deprecated::
        This class is part of the legacy ``analysis/`` infrastructure and
        will be removed once all calculators are migrated to the
        ``analyses/`` plugin system.
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
    def compute(self, replicate: int, **kwargs) -> Any:
        """Run analysis for a single replicate.

        Parameters
        ----------
        replicate : int
            Replicate number to analyze.
        **kwargs
            Additional analyzer-specific parameters.

        Returns
        -------
        Any
            Analysis result for this replicate.
        """
        ...

    @abstractmethod
    def compute_aggregated(self, replicates: "Sequence[int]", **kwargs) -> Any:
        """Run aggregated analysis across multiple replicates.

        Parameters
        ----------
        replicates : Sequence[int]
            List of replicate numbers to analyze.
        **kwargs
            Additional analyzer-specific parameters.

        Returns
        -------
        Any
            Aggregated analysis result.
        """
        ...

    @property
    def label(self) -> str:
        """Human-readable label for this analyzer.

        Returns
        -------
        str
            Display label for reports and logs.
        """
        return self.analysis_type()


# ============================================================================
# Registry Classes
# ============================================================================


class AnalysisSettingsRegistry:
    """Registry for analysis settings types.

    Allows new analysis types to be registered without modifying core code.
    Use the ``register`` decorator to add new analysis settings classes.

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
    Use the ``register`` decorator to add new comparison settings classes.

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

    .. deprecated::
        This registry is part of the legacy ``analysis/`` infrastructure.
        New analysis types should use the ``analyses/`` plugin system instead.
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
