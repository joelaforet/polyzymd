"""Config-driven plotting for comparison results.

This module provides extensible infrastructure for generating plots from
comparison.yaml configurations. It follows the Open-Closed Principle (OCP)
using a registry pattern - new plot types can be added by registering
with PlotterRegistry without modifying core code.

Architecture
------------
- BasePlotter: Abstract base class for all plotters
- PlotterRegistry: Registry for plot type implementations
- ComparisonPlotter: Main orchestrator that uses the registry

Adding New Plot Types
---------------------
Contributors can add new plot types by:

1. Create a plotter class inheriting from BasePlotter
2. Implement the required methods (plot_type, can_plot, plot)
3. Register with the PlotterRegistry decorator

Example::

    from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

    @PlotterRegistry.register("my_plot")
    class MyPlotter(BasePlotter):
        @classmethod
        def plot_type(cls) -> str:
            return "my_plot"

        def can_plot(self, comparison_config, analysis_type: str) -> bool:
            return analysis_type == "my_analysis"

        def plot(self, data, labels, settings, output_path) -> Path:
            # Generate plot and save
            ...
            return output_path
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence, Type

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from polyzymd.compare.config import ComparisonConfig, PlotSettings

logger = logging.getLogger(__name__)


# ============================================================================
# Abstract Base Class
# ============================================================================


class BasePlotter(ABC):
    """Abstract base class for all plotters.

    Plotters generate specific plot types from analysis results.
    Each plotter handles one or more related plot types and knows
    how to load the required data and generate the visualization.

    Subclasses must implement:
    - plot_type(): Returns the unique identifier for this plotter
    - can_plot(): Determines if this plotter can handle given data
    - plot(): Generates and saves the plot

    Attributes
    ----------
    settings : PlotSettings
        Global plot settings from comparison.yaml
    """

    def __init__(self, settings: "PlotSettings"):
        """Initialize plotter with settings.

        Parameters
        ----------
        settings : PlotSettings
            Global plot settings from comparison.yaml
        """
        self.settings = settings

    @classmethod
    @abstractmethod
    def plot_type(cls) -> str:
        """Return the unique identifier for this plotter.

        Returns
        -------
        str
            Plot type identifier (e.g., "triad_kde_panel", "rmsf_comparison")
        """
        ...

    @abstractmethod
    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the given analysis type.

        Parameters
        ----------
        comparison_config : ComparisonConfig
            Full comparison configuration
        analysis_type : str
            Analysis type to check (e.g., "rmsf", "triad", "distances")

        Returns
        -------
        bool
            True if this plotter can generate plots for the analysis type
        """
        ...

    @abstractmethod
    def plot(
        self,
        data: Any,
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate and save plot(s).

        This method receives pre-loaded condition metadata from
        `ComparisonPlotter._load_analysis_data()` and must load its own
        analysis results from the filesystem.

        Parameters
        ----------
        data : dict[str, dict]
            Mapping of condition_label -> condition data dict containing:
            - "condition": ConditionConfig object with condition metadata
            - "sim_config": SimulationConfig object with full config
            - "analysis_dir": Path to analysis/{analysis_type}/ directory
            - "aggregated_dir": Path to analysis/{analysis_type}/aggregated/
            - "replicates": list[int] of replicate numbers

            **IMPORTANT**: Plotters must load their own analysis results from
            `analysis_dir` or `aggregated_dir`. The orchestrator does NOT pass
            pre-loaded results via kwargs.
        labels : sequence of str
            Condition labels in desired display order
        output_dir : Path
            Directory to save plot files
        **kwargs
            Reserved for future use. Do NOT expect analysis results here.

        Returns
        -------
        list[Path]
            Paths to generated plot files (may be empty if no data available)

        Examples
        --------
        Correct pattern for loading data in a plotter:

        >>> def plot(self, data, labels, output_dir, **kwargs):
        ...     for label in labels:
        ...         analysis_dir = Path(data[label]["analysis_dir"])
        ...         result_file = analysis_dir / "my_result.json"
        ...         result = MyResult.load(result_file)
        ...         # ... generate plot from result ...
        """
        ...

    def _get_output_path(self, output_dir: Path, name: str) -> Path:
        """Generate output file path with correct format extension.

        Parameters
        ----------
        output_dir : Path
            Output directory
        name : str
            Base filename (without extension)

        Returns
        -------
        Path
            Full output path with extension
        """
        return output_dir / f"{name}.{self.settings.format}"

    def _save_figure(self, fig: "Figure", output_path: Path) -> Path:
        """Save figure with settings.

        Parameters
        ----------
        fig : Figure
            Matplotlib figure to save
        output_path : Path
            Output file path

        Returns
        -------
        Path
            Path to saved figure
        """
        import matplotlib.pyplot as plt

        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            output_path,
            dpi=self.settings.dpi,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
        plt.close(fig)
        logger.info(f"Saved plot: {output_path}")
        return output_path


# ============================================================================
# Plotter Registry
# ============================================================================


class PlotterRegistry:
    """Registry for plotter implementations.

    Allows new plot types to be registered without modifying core code.
    Use the `register` decorator to add new plotter classes.

    Examples
    --------
    >>> @PlotterRegistry.register("triad_kde_panel")
    ... class TriadKDEPanelPlotter(BasePlotter):
    ...     @classmethod
    ...     def plot_type(cls) -> str:
    ...         return "triad_kde_panel"
    ...     ...
    >>>
    >>> # Get all plotters that can handle "catalytic_triad" analysis
    >>> plotters = PlotterRegistry.get_plotters_for_analysis("catalytic_triad", config)
    """

    _registry: dict[str, Type[BasePlotter]] = {}

    @classmethod
    def register(cls, name: str | None = None):
        """Decorator to register a plotter class.

        Parameters
        ----------
        name : str, optional
            Registry key. If None, uses the class's plot_type().

        Returns
        -------
        Callable
            Decorator function.

        Examples
        --------
        >>> @PlotterRegistry.register("rmsf_comparison")
        ... class RMSFComparisonPlotter(BasePlotter):
        ...     ...
        """

        def decorator(plotter_class: Type[BasePlotter]):
            key = name if name is not None else plotter_class.plot_type()
            key = key.lower()
            if key in cls._registry:
                logger.warning(f"Overwriting existing plotter registration: {key}")
            cls._registry[key] = plotter_class
            logger.debug(f"Registered plotter: {key}")
            return plotter_class

        return decorator

    @classmethod
    def get(cls, name: str) -> Type[BasePlotter]:
        """Get plotter class by name.

        Parameters
        ----------
        name : str
            Plot type identifier.

        Returns
        -------
        Type[BasePlotter]
            The registered plotter class.

        Raises
        ------
        ValueError
            If the plot type is not registered.
        """
        key = name.lower()
        if key not in cls._registry:
            available = ", ".join(sorted(cls._registry.keys()))
            raise ValueError(f"Unknown plot type: '{name}'. Available: {available}")
        return cls._registry[key]

    @classmethod
    def list_available(cls) -> list[str]:
        """List all registered plot types.

        Returns
        -------
        list[str]
            Sorted list of registered plot type names.
        """
        return sorted(cls._registry.keys())

    @classmethod
    def is_registered(cls, name: str) -> bool:
        """Check if a plot type is registered.

        Parameters
        ----------
        name : str
            Plot type identifier.

        Returns
        -------
        bool
            True if registered, False otherwise.
        """
        return name.lower() in cls._registry

    @classmethod
    def get_plotters_for_analysis(
        cls,
        analysis_type: str,
        comparison_config: "ComparisonConfig",
        settings: "PlotSettings",
    ) -> list[BasePlotter]:
        """Get all plotters that can handle an analysis type.

        Parameters
        ----------
        analysis_type : str
            Analysis type (e.g., "rmsf", "catalytic_triad")
        comparison_config : ComparisonConfig
            Full comparison configuration
        settings : PlotSettings
            Plot settings for instantiating plotters

        Returns
        -------
        list[BasePlotter]
            Instantiated plotters that can handle this analysis type
        """
        compatible = []
        for plotter_class in cls._registry.values():
            plotter = plotter_class(settings)
            if plotter.can_plot(comparison_config, analysis_type):
                compatible.append(plotter)
        return compatible

    @classmethod
    def clear(cls) -> None:
        """Clear the registry (for testing purposes)."""
        cls._registry.clear()


# ============================================================================
# Main Orchestrator
# ============================================================================


class ComparisonPlotter:
    """Config-driven plot generator for comparison results.

    This class orchestrates plot generation based on comparison.yaml
    settings. It uses the PlotterRegistry to find appropriate plotters
    for each enabled analysis type.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration loaded from comparison.yaml

    Examples
    --------
    >>> from polyzymd.compare.config import ComparisonConfig
    >>> from polyzymd.compare.plotter import ComparisonPlotter
    >>>
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> plotter = ComparisonPlotter(config)
    >>> generated_files = plotter.plot_all()
    >>> print(f"Generated {len(generated_files)} plots")
    """

    def __init__(self, config: "ComparisonConfig"):
        """Initialize plotter with comparison config.

        Parameters
        ----------
        config : ComparisonConfig
            Comparison configuration
        """
        self.config = config
        self.settings = config.plot_settings
        self._output_dir = self._resolve_output_dir()

    def _resolve_output_dir(self) -> Path:
        """Resolve the output directory relative to comparison.yaml location.

        Returns
        -------
        Path
            Absolute path to output directory
        """
        output_dir = self.settings.output_dir
        if not output_dir.is_absolute():
            # Use comparison.yaml parent directory as base
            if self.config.source_path is not None:
                base = self.config.source_path.parent
            elif self.config.conditions:
                # Fallback: use first condition's config parent
                base = self.config.conditions[0].config.parent
            else:
                base = Path.cwd()
            output_dir = base / output_dir
        return output_dir.resolve()

    @property
    def output_dir(self) -> Path:
        """Get the resolved output directory."""
        return self._output_dir

    def list_available_plots(self) -> dict[str, list[str]]:
        """List all available plots for enabled analyses.

        Returns
        -------
        dict[str, list[str]]
            Mapping of analysis_type -> list of plot_type names
        """
        result = {}
        enabled = self.config.analysis_settings.get_enabled_analyses()

        for analysis_type in enabled:
            plotters = PlotterRegistry.get_plotters_for_analysis(
                analysis_type, self.config, self.settings
            )
            result[analysis_type] = [p.plot_type() for p in plotters]

        return result

    def plot_all(self) -> list[Path]:
        """Generate all plots for all enabled analyses.

        Returns
        -------
        list[Path]
            Paths to all generated plot files

        Notes
        -----
        This method iterates through all enabled analyses in the config
        and generates all registered plot types for each.
        """
        generated = []
        enabled = self.config.analysis_settings.get_enabled_analyses()

        logger.info(f"Generating plots for {len(enabled)} analysis types")
        logger.info(f"Output directory: {self._output_dir}")

        for analysis_type in enabled:
            try:
                paths = self.plot_analysis(analysis_type)
                generated.extend(paths)
            except Exception as e:
                logger.error(f"Failed to generate plots for {analysis_type}: {e}")

        logger.info(f"Generated {len(generated)} total plots")
        return generated

    def plot_analysis(self, analysis_type: str) -> list[Path]:
        """Generate all plots for a specific analysis type.

        Parameters
        ----------
        analysis_type : str
            Analysis type (e.g., "rmsf", "catalytic_triad")

        Returns
        -------
        list[Path]
            Paths to generated plot files
        """
        plotters = PlotterRegistry.get_plotters_for_analysis(
            analysis_type, self.config, self.settings
        )

        if not plotters:
            logger.warning(f"No plotters registered for analysis type: {analysis_type}")
            return []

        generated = []
        labels = [c.label for c in self.config.conditions]

        for plotter in plotters:
            try:
                logger.info(f"Running plotter: {plotter.plot_type()}")

                # Load data for this analysis type
                data = self._load_analysis_data(analysis_type)

                # Generate plots
                output_subdir = self._output_dir / analysis_type
                paths = plotter.plot(data, labels, output_subdir)
                generated.extend(paths)

            except Exception as e:
                logger.error(f"Plotter {plotter.plot_type()} failed: {e}")

        return generated

    def plot_single(self, plot_type: str, analysis_type: str) -> list[Path]:
        """Generate a specific plot type.

        Parameters
        ----------
        plot_type : str
            Plot type identifier (e.g., "triad_kde_panel")
        analysis_type : str
            Analysis type this plot is for

        Returns
        -------
        list[Path]
            Paths to generated plot files
        """
        plotter_class = PlotterRegistry.get(plot_type)
        plotter = plotter_class(self.settings)

        if not plotter.can_plot(self.config, analysis_type):
            raise ValueError(f"Plotter '{plot_type}' cannot handle analysis type '{analysis_type}'")

        labels = [c.label for c in self.config.conditions]
        data = self._load_analysis_data(analysis_type)
        output_subdir = self._output_dir / analysis_type

        return plotter.plot(data, labels, output_subdir)

    def _load_analysis_data(self, analysis_type: str) -> dict[str, Any]:
        """Load analysis results for all conditions.

        Parameters
        ----------
        analysis_type : str
            Analysis type to load

        Returns
        -------
        dict[str, Any]
            Mapping of condition_label -> loaded data, plus a special
            ``"__meta__"`` key with comparison-level metadata.

            Per-condition entries contain:
            - ``condition``: ConditionConfig object
            - ``sim_config``: SimulationConfig object
            - ``analysis_dir``: Path to analysis/{analysis_type}/ directory
            - ``aggregated_dir``: Path to analysis/{analysis_type}/aggregated/
            - ``replicates``: list[int] of replicate numbers

            The ``__meta__`` entry contains:
            - ``comparison_source_path``: Path to comparison.yaml (or None)
            - ``results_dir``: Path to results/ adjacent to comparison.yaml (or None)
        """
        from polyzymd.config.schema import SimulationConfig

        data: dict[str, Any] = {}

        # Provide comparison-level metadata so cross-condition plotters
        # (e.g., BFE) can locate results/ without heuristic path guessing.
        source_path = self.config.source_path
        results_dir = source_path.parent / "results" if source_path is not None else None
        data["__meta__"] = {
            "comparison_source_path": source_path,
            "results_dir": results_dir,
        }

        for condition in self.config.conditions:
            try:
                # Load the simulation config to find the projects directory
                sim_config = SimulationConfig.from_yaml(condition.config)
                projects_dir = sim_config.output.projects_directory

                # Build path to analysis results - with fallback to config parent
                analysis_dir = projects_dir / "analysis" / analysis_type
                if not analysis_dir.exists() and condition.config is not None:
                    # Fallback: check relative to condition config file
                    fallback_dir = condition.config.parent / "analysis" / analysis_type
                    if fallback_dir.exists():
                        analysis_dir = fallback_dir

                # Load aggregated results if available
                aggregated_dir = analysis_dir / "aggregated"

                data[condition.label] = {
                    "condition": condition,
                    "sim_config": sim_config,
                    "analysis_dir": analysis_dir,
                    "aggregated_dir": aggregated_dir,
                    "replicates": condition.replicates,
                }

            except Exception as e:
                logger.warning(f"Failed to load data for condition '{condition.label}': {e}")
                data[condition.label] = None

        return data


# ============================================================================
# Built-in Plotters - Import triggers registration
# ============================================================================


# Import built-in plotters to register them
# These are imported at the end to avoid circular imports
def _register_builtin_plotters():
    """Register all built-in plotters.

    This is called when the module is imported, but can also be called
    explicitly to ensure all plotters are registered.
    """
    # Import submodules that contain @PlotterRegistry.register decorators
    try:
        from polyzymd.compare import plotters  # noqa: F401
    except ImportError:
        logger.debug("Built-in plotters module not found - skipping registration")


# Register built-in plotters on import
_register_builtin_plotters()
