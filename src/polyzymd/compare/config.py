"""Configuration schema for comparison projects.

This module defines the YAML schema for comparison.yaml files that
specify which simulation conditions to compare.

The schema has two main sections:
- analysis_settings: Defines WHAT analyses to run (shared across conditions)
- comparison_settings: Defines HOW to compare (statistical parameters)

Both sections use a registry-based approach for extensibility. New analysis
types can be added by registering with AnalysisSettingsRegistry and
ComparisonSettingsRegistry (see polyzymd.compare.settings).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, ClassVar

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analysis.config import AnalysisDefaults
from polyzymd.analysis.core.registry import (
    AnalysisSettingsRegistry,
    BaseAnalysisSettings,
    BaseComparisonSettings,
    BasePlotSettings,
    ComparisonSettingsRegistry,
    PlotSettingsRegistry,
)

# Import settings to trigger registration
from polyzymd.compare.settings import (  # noqa: F401
    BindingFreeEnergyAnalysisSettings,
    BindingFreeEnergyComparisonSettings,
    CatalyticTriadAnalysisSettings,
    CatalyticTriadComparisonSettings,
    ContactsAnalysisSettings,
    ContactsComparisonSettings,
    DistancePairSettings,
    DistancesAnalysisSettings,
    DistancesComparisonSettings,
    RMSFAnalysisSettings,
    RMSFComparisonSettings,
    SecondaryStructureAnalysisSettings,
    SecondaryStructureComparisonSettings,
    TriadPairSettings,
)

# Backward-compatible aliases for analysis module
# The analysis/triad module still imports these old names
CatalyticTriadConfig = CatalyticTriadAnalysisSettings
TriadPairConfig = TriadPairSettings

logger = logging.getLogger(__name__)


# ============================================================================
# Condition Configuration
# ============================================================================


class ConditionConfig(BaseModel):
    """Configuration for one condition in a comparison.

    Attributes
    ----------
    label : str
        Display name for this condition (e.g., "No Polymer", "100% SBMA")
    config : Path
        Path to the simulation's config.yaml file
    replicates : list[int]
        List of replicate numbers to include in the analysis
    """

    label: str
    config: Path
    replicates: list[int]

    @field_validator("config", mode="before")
    @classmethod
    def resolve_path(cls, v: str | Path) -> Path:
        """Convert string paths to Path objects."""
        return Path(v)

    @field_validator("replicates", mode="before")
    @classmethod
    def ensure_list(cls, v: list[int] | int) -> list[int]:
        """Ensure replicates is a list."""
        if isinstance(v, int):
            return [v]
        return list(v)


# ============================================================================
# Dynamic Settings Containers
# ============================================================================


class AnalysisSettingsContainer(BaseModel):
    """Container for analysis settings (WHAT to analyze).

    Uses dynamic attribute access to support any registered analysis type
    without hardcoding field names.
    """

    model_config = {"extra": "allow"}

    def __init__(self, **data: Any):
        """Initialize with dynamic analysis settings.

        Parameters
        ----------
        **data : Any
            Analysis settings keyed by analysis type name.
        """
        # Parse each setting using the registry
        parsed_settings: dict[str, BaseAnalysisSettings] = {}
        for key, value in data.items():
            if value is None:
                continue
            key_lower = key.lower()
            if AnalysisSettingsRegistry.is_registered(key_lower):
                settings_class = AnalysisSettingsRegistry.get(key_lower)
                if isinstance(value, dict):
                    parsed_settings[key_lower] = settings_class(**value)
                elif isinstance(value, BaseAnalysisSettings):
                    parsed_settings[key_lower] = value
                else:
                    raise ValueError(
                        f"Invalid value for {key}: expected dict or {settings_class.__name__}"
                    )
            else:
                logger.warning(f"Unknown analysis type '{key}' - skipping")

        super().__init__(**parsed_settings)

    def get(self, analysis_type: str) -> BaseAnalysisSettings | None:
        """Get settings for a specific analysis type.

        Parameters
        ----------
        analysis_type : str
            Analysis type identifier (e.g., "rmsf", "contacts").

        Returns
        -------
        BaseAnalysisSettings or None
            Settings for the analysis type, or None if not configured.
        """
        return getattr(self, analysis_type.lower(), None)

    def get_enabled_analyses(self) -> list[str]:
        """Get list of enabled analysis types.

        Returns
        -------
        list[str]
            Names of configured analyses (presence implies enabled).

        Notes
        -----
        Uses actual model data from comparison.yaml rather than relying on
        a registry. This makes comparison.yaml the source of truth for which
        analyses are enabled.
        """
        # Use actual model data - if an analysis section exists and has a value, it's enabled
        return [key for key, value in self.model_dump().items() if value is not None]

    def to_analysis_yaml_dict(self, replicates: list[int], eq_time: str) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary.

        Parameters
        ----------
        replicates : list[int]
            Replicate numbers for the analysis.yaml.
        eq_time : str
            Equilibration time for the analysis.yaml.

        Returns
        -------
        dict[str, Any]
            Dictionary suitable for writing to analysis.yaml.
        """
        result: dict[str, Any] = {
            "replicates": replicates,
            "defaults": {"equilibration_time": eq_time},
        }
        for analysis_type in self.get_enabled_analyses():
            settings = self.get(analysis_type)
            if settings is not None:
                result[analysis_type] = settings.to_analysis_yaml_dict()
        return result


class ComparisonSettingsContainer(BaseModel):
    """Container for comparison settings (HOW to compare).

    Uses dynamic attribute access to support any registered comparison type.
    Each analysis type in analysis_settings must have a corresponding entry
    here (can be empty dict) to enable comparison.
    """

    model_config = {"extra": "allow"}

    def __init__(self, **data: Any):
        """Initialize with dynamic comparison settings.

        Parameters
        ----------
        **data : Any
            Comparison settings keyed by analysis type name.
        """
        # Parse each setting using the registry
        parsed_settings: dict[str, BaseComparisonSettings] = {}
        for key, value in data.items():
            if value is None:
                continue
            key_lower = key.lower()
            if ComparisonSettingsRegistry.is_registered(key_lower):
                settings_class = ComparisonSettingsRegistry.get(key_lower)
                if isinstance(value, dict):
                    parsed_settings[key_lower] = settings_class(**value)
                elif isinstance(value, BaseComparisonSettings):
                    parsed_settings[key_lower] = value
                else:
                    raise ValueError(
                        f"Invalid value for {key}: expected dict or {settings_class.__name__}"
                    )
            else:
                logger.warning(f"Unknown comparison type '{key}' - skipping")

        super().__init__(**parsed_settings)

    def get(self, analysis_type: str) -> BaseComparisonSettings | None:
        """Get settings for a specific comparison type.

        Parameters
        ----------
        analysis_type : str
            Analysis type identifier (e.g., "rmsf", "contacts").

        Returns
        -------
        BaseComparisonSettings or None
            Comparison settings, or None if not configured.
        """
        return getattr(self, analysis_type.lower(), None)

    def get_enabled_comparisons(self) -> list[str]:
        """Get list of enabled comparison types.

        Returns
        -------
        list[str]
            Names of configured comparisons.
        """
        enabled = []
        for analysis_type in ComparisonSettingsRegistry.list_available():
            if self.get(analysis_type) is not None:
                enabled.append(analysis_type)
        return enabled


# ============================================================================
# Plot Settings Configuration
# ============================================================================


@PlotSettingsRegistry.register("rmsf")
class RMSFPlotSettings(BasePlotSettings):
    """RMSF-specific plot customization.

    Attributes
    ----------
    show_error : bool
        Show error bands/bars on plots (default True)
    highlight_residues : list[int]
        Residue numbers to highlight with vertical lines (e.g., active site)
    figsize_profile : tuple[float, float]
        Figure size for per-residue profile plots
    figsize_comparison : tuple[float, float]
        Figure size for bar comparison plots
    """

    show_error: bool = True
    highlight_residues: list[int] = Field(default_factory=list)
    figsize_profile: tuple[float, float] = (14, 4)
    figsize_comparison: tuple[float, float] = (8, 6)


@PlotSettingsRegistry.register("triad")
class TriadPlotSettings(BasePlotSettings):
    """Triad-specific plot customization.

    Attributes
    ----------
    generate_kde_panel : bool
        Generate multi-row KDE panel plot (default True)
    generate_bars : bool
        Generate grouped threshold bar chart (default True)
    generate_2d_kde : bool
        Generate 2D joint KDE plot (default False, more specialized)
    threshold_line_color : str
        Color for threshold vertical line
    kde_fill_alpha : float
        Transparency for KDE fill (0-1)
    kde_xlim : tuple[float, float]
        X-axis limits for KDE panel in Angstroms (default ``(0, 7)``).
    figsize_kde_panel : tuple[float, float] | None
        Figure size for KDE panel (auto-calculated if None)
    figsize_bars : tuple[float, float]
        Figure size for bar chart
    """

    generate_kde_panel: bool = True
    generate_bars: bool = True
    generate_2d_kde: bool = False
    threshold_line_color: str = "red"
    kde_fill_alpha: float = 0.7
    kde_xlim: tuple[float, float] = (0.0, 7.0)
    figsize_kde_panel: tuple[float, float] | None = None
    figsize_bars: tuple[float, float] = (10, 6)


@PlotSettingsRegistry.register("distances")
class DistancesPlotSettings(BasePlotSettings):
    """Distance analysis plot customization.

    Attributes
    ----------
    show_threshold : bool
        Show threshold line on distribution plots
    use_kde : bool
        Use KDE instead of histogram for distributions
    generate_state_bars : bool
        Generate per-pair state bar charts (above/below threshold).
        Each pair gets its own figure showing the fraction of frames
        in each state per condition. Default True.
    figsize : tuple[float, float]
        Default figure size for distance plots
    """

    show_threshold: bool = True
    use_kde: bool = True
    generate_state_bars: bool = True
    figsize: tuple[float, float] = (10, 6)


@PlotSettingsRegistry.register("contacts")
class ContactsPlotSettings(BasePlotSettings):
    """Contacts analysis plot customization.

    Attributes
    ----------
    figsize : tuple[float, float]
        Default figure size for contact plots
    generate_enrichment_heatmap : bool
        Generate binding preference enrichment heatmap (default True)
    generate_enrichment_bars : bool
        Generate binding preference bar charts (default True)
    figsize_enrichment_heatmap : tuple[float, float] | None
        Figure size for enrichment heatmap (auto-calculated if None)
    figsize_enrichment_bars : tuple[float, float]
        Figure size for enrichment bar charts
    enrichment_colormap : str
        Colormap for enrichment heatmap (diverging recommended)
    show_enrichment_error : bool
        Show error bars on enrichment bar charts (default True)
    generate_system_coverage_heatmap : bool
        Generate system coverage enrichment heatmap (default True)
    generate_system_coverage_bars : bool
        Generate system coverage bar charts (default True)
    figsize_system_coverage_heatmap : tuple[float, float] | None
        Figure size for system coverage heatmap (auto-calculated if None)
    figsize_system_coverage_bars : tuple[float, float]
        Figure size for system coverage bar charts
    show_system_coverage_error : bool
        Show error bars on system coverage bar charts (default True)
    generate_user_partition_bars : bool
        Generate user-defined partition bar charts (default True)
    figsize_user_partition_bars : tuple[float, float]
        Figure size for user-defined partition bar charts
    show_user_partition_error : bool
        Show error bars on user-defined partition bar charts (default True)
    generate_contact_fraction_profile : bool
        Generate per-residue contact fraction line plot (default True)
    figsize_contact_fraction_profile : tuple[float, float]
        Figure size for contact fraction profile plot
    show_contact_fraction_profile_error : bool
        Show SEM fill_between bands on contact fraction profile (default True)
    contact_fraction_profile_threshold : float or None
        If set, draw a horizontal threshold line on the contact fraction
        profile. Residues above this value are considered "high contact".
    generate_residence_time_profile : bool
        Generate per-residue mean residence time line plot (default True)
    figsize_residence_time_profile : tuple[float, float]
        Figure size for residence time profile plot
    show_residence_time_profile_error : bool
        Show SEM fill_between bands on residence time profile (default True)
    generate_cf_by_aa_class_bars : bool
        Generate contact fraction by AA class grouped bar chart (default True)
    figsize_cf_by_aa_class_bars : tuple[float, float]
        Figure size for contact fraction by AA class bar chart
    show_cf_by_aa_class_error : bool
        Show error bars on contact fraction by AA class bar chart (default True)
    generate_cf_by_partition_bars : bool
        Generate contact fraction by user-defined partition bar charts (default True)
    figsize_cf_by_partition_bars : tuple[float, float]
        Figure size for contact fraction by partition bar charts
    show_cf_by_partition_error : bool
        Show error bars on contact fraction by partition bar charts (default True)
    generate_rt_by_aa_class_bars : bool
        Generate residence time by AA class grouped bar chart (default True)
    figsize_rt_by_aa_class_bars : tuple[float, float]
        Figure size for residence time by AA class bar chart
    show_rt_by_aa_class_error : bool
        Show error bars on residence time by AA class bar chart (default True)
    generate_rt_by_partition_bars : bool
        Generate residence time by user-defined partition bar charts (default True)
    figsize_rt_by_partition_bars : tuple[float, float]
        Figure size for residence time by partition bar charts
    show_rt_by_partition_error : bool
        Show error bars on residence time by partition bar charts (default True)
    highlight_residues : list[int]
        Residue IDs to mark with vertical dashed lines on profile plots.
        Useful for highlighting active-site residues or known anchor points.
    """

    figsize: tuple[float, float] = (10, 8)
    generate_enrichment_heatmap: bool = True
    generate_enrichment_bars: bool = True
    figsize_enrichment_heatmap: tuple[float, float] | None = None
    figsize_enrichment_bars: tuple[float, float] = (10, 6)
    enrichment_colormap: str = "RdBu_r"  # Diverging: red=high, blue=low
    show_enrichment_error: bool = True

    # System coverage plot settings
    generate_system_coverage_heatmap: bool = True
    generate_system_coverage_bars: bool = True
    figsize_system_coverage_heatmap: tuple[float, float] | None = None
    figsize_system_coverage_bars: tuple[float, float] = (10, 6)
    show_system_coverage_error: bool = True

    # User-defined partition plot settings
    generate_user_partition_bars: bool = True
    figsize_user_partition_bars: tuple[float, float] = (10, 6)
    show_user_partition_error: bool = True

    # Contact fraction profile plot settings
    generate_contact_fraction_profile: bool = True
    figsize_contact_fraction_profile: tuple[float, float] = (14, 5)
    show_contact_fraction_profile_error: bool = True
    contact_fraction_profile_threshold: float | None = None

    # Residence time profile plot settings
    generate_residence_time_profile: bool = True
    figsize_residence_time_profile: tuple[float, float] = (14, 5)
    show_residence_time_profile_error: bool = True

    # Contact fraction by AA class bar chart settings
    generate_cf_by_aa_class_bars: bool = True
    figsize_cf_by_aa_class_bars: tuple[float, float] = (10, 6)
    show_cf_by_aa_class_error: bool = True

    # Contact fraction by user partition bar chart settings
    generate_cf_by_partition_bars: bool = True
    figsize_cf_by_partition_bars: tuple[float, float] = (10, 6)
    show_cf_by_partition_error: bool = True

    # Residence time by AA class bar chart settings
    generate_rt_by_aa_class_bars: bool = True
    figsize_rt_by_aa_class_bars: tuple[float, float] = (10, 6)
    show_rt_by_aa_class_error: bool = True

    # Residence time by user partition bar chart settings
    generate_rt_by_partition_bars: bool = True
    figsize_rt_by_partition_bars: tuple[float, float] = (10, 6)
    show_rt_by_partition_error: bool = True

    # Shared profile plot settings
    highlight_residues: list[int] = Field(default_factory=list)


@PlotSettingsRegistry.register("binding_free_energy")
class BFEPlotSettings(BasePlotSettings):
    """Binding free energy plot customization.

    Attributes
    ----------
    generate_heatmap : bool
        Generate ΔG_sel heatmap (rows = AA groups, columns = conditions). Default True.
    generate_bars : bool
        Generate ΔG_sel grouped bar chart (one bar per condition per AA group). Default True.
    figsize_heatmap : tuple[float, float] | None
        Figure size for ΔG_sel heatmap (auto-calculated if None).
    figsize_bars : tuple[float, float]
        Figure size for ΔG_sel bar charts.
    colormap : str
        Diverging colormap for heatmap (default "RdBu_r": red = avoidance, blue = preference).
    show_error_bars : bool
        Show SEM error bars on bar charts. Default True.
    annotate_heatmap : bool
        Annotate each heatmap cell with its ΔG_sel value. Default True.
    """

    generate_heatmap: bool = True
    generate_bars: bool = True
    figsize_heatmap: tuple[float, float] | None = None
    figsize_bars: tuple[float, float] = (10, 6)
    colormap: str = "RdBu_r"
    show_error_bars: bool = True
    annotate_heatmap: bool = True


@PlotSettingsRegistry.register("polymer_affinity")
class AffinityPlotSettings(BasePlotSettings):
    """Polymer affinity score plot customization.

    Attributes
    ----------
    generate_stacked_bars : bool
        Generate stacked bar chart of total score by condition, broken
        down by polymer type. Default True.
    generate_group_bars : bool
        Generate grouped bar chart showing per-group contributions
        across conditions. Default True.
    figsize_stacked : tuple[float, float]
        Figure size for stacked bar chart.
    figsize_group_bars : tuple[float, float]
        Figure size for grouped bar charts.
    show_error_bars : bool
        Show SEM error bars on plots. Default True.
    """

    generate_stacked_bars: bool = True
    generate_group_bars: bool = True
    figsize_stacked: tuple[float, float] = (10, 6)
    figsize_group_bars: tuple[float, float] = (10, 6)
    show_error_bars: bool = True


@PlotSettingsRegistry.register("secondary_structure")
class SSPlotSettings(BasePlotSettings):
    """Secondary structure plot customization.

    Attributes
    ----------
    generate_timeline : bool
        Generate per-condition residue x time SS heatmap. Default True.
    generate_content_bars : bool
        Generate grouped bar chart of helix/strand/coil fractions. Default True.
    generate_individual_bars : bool
        Generate one bar chart per SS type (helix, beta-sheet, no-SS). Default True.
    generate_diff_heatmap : bool
        Generate condition x residue persistence difference heatmap. Default True.
    figsize_timeline : tuple[float, float]
        Figure size for timeline heatmap.
    figsize_content_bars : tuple[float, float]
        Figure size for content bar chart.
    figsize_diff_heatmap : tuple[float, float] | None
        Figure size for difference heatmap (auto-calculated if None).
    diff_colormap : str
        Diverging colormap for difference heatmap.
    """

    generate_timeline: bool = True
    generate_content_bars: bool = True
    generate_individual_bars: bool = True
    generate_diff_heatmap: bool = True
    figsize_timeline: tuple[float, float] = (14, 6)
    figsize_content_bars: tuple[float, float] = (10, 6)
    figsize_diff_heatmap: tuple[float, float] | None = None
    diff_colormap: str = "RdBu_r"


class PlotTheme(BaseModel):
    """Centralized visual defaults for all comparison plots.

    Replaces ~219 hardcoded style values (font sizes, alphas, line widths,
    marker sizes, spine visibility, etc.) across all plotter files with a
    single configurable Pydantic model.

    Three presets are available via class methods:

    - ``PlotTheme.publication()`` — default; print-ready sizes and weights.
    - ``PlotTheme.presentation()`` — ~1.3x larger fonts/dots/lines for slides.
    - ``PlotTheme.minimal()`` — no dots, no bar edges, thinner lines.

    Users can override individual values in ``comparison.yaml``::

        plot_settings:
          style: "publication"
          theme:
            title_fontsize: 16
            dot_size: 24

    Parameters
    ----------
    title_fontsize : int
        Font size for axes titles.
    suptitle_fontsize : int
        Font size for figure suptitles.
    label_fontsize : int
        Font size for axis labels (xlabel/ylabel).
    tick_fontsize : int
        Font size for tick labels.
    legend_fontsize : int
        Font size for legend entries.
    annotation_fontsize : int
        Font size for heatmap cell annotations and inline text.
    small_fontsize : int
        Font size for secondary annotations (e.g. SEM ± labels).
    tiny_fontsize : int
        Font size for fine-grained annotations (e.g. residue IDs).
    bar_alpha : float
        Opacity for bar chart fill.
    bar_edgecolor : str
        Edge colour for bar outlines.
    bar_linewidth : float
        Edge line width for bars.
    bar_capsize : int
        Error bar cap size in points.
    dot_size : int
        Marker size for replicate dot overlays (``s=`` in ``scatter``).
    dot_alpha : float
        Opacity for replicate dots.
    dot_color : str
        Colour for replicate dots.
    line_alpha : float
        Opacity for line plots (e.g. RMSF profiles).
    fill_alpha : float
        Opacity for fill_between bands (e.g. SEM regions).
    reference_line_color : str
        Colour for horizontal/vertical reference lines.
    reference_line_style : str
        Linestyle for reference lines (e.g. ``"--"``).
    reference_line_width : float
        Line width for reference lines.
    highlight_line_alpha : float
        Opacity for highlight / vertical reference lines.
    hide_top_spine : bool
        Whether to hide the top axis spine.
    hide_right_spine : bool
        Whether to hide the right axis spine.
    title_fontweight : str
        Font weight for titles (e.g. ``"bold"``, ``"normal"``).
    legend_loc : str
        Matplotlib legend location string (e.g. ``"center left"``).
        Used with ``legend_bbox`` to place the legend outside the axes.
    legend_bbox : tuple of float
        ``bbox_to_anchor`` for legend placement, relative to axes.
        Default ``(1.02, 0.5)`` places it just outside the right edge,
        vertically centred.
    show_watermark : bool
        Whether to render a subtle "Made by PolyzyMD" watermark in the
        bottom-right corner of every saved figure.  Default ``True``.
    """

    # Font sizes by semantic role
    title_fontsize: int = 13
    suptitle_fontsize: int = 14
    label_fontsize: int = 11
    tick_fontsize: int = 9
    legend_fontsize: int = 9
    annotation_fontsize: int = 9
    small_fontsize: int = 8
    tiny_fontsize: int = 7

    # Bar chart defaults
    bar_alpha: float = 0.85
    bar_edgecolor: str = "black"
    bar_linewidth: float = 0.5
    bar_capsize: int = 4

    # Replicate dot overlay
    dot_size: int = 18
    dot_alpha: float = 0.7
    dot_color: str = "black"

    # Line defaults
    line_alpha: float = 0.8
    fill_alpha: float = 0.25
    reference_line_color: str = "black"
    reference_line_style: str = "--"
    reference_line_width: float = 1.5
    highlight_line_alpha: float = 0.5

    # Axes chrome
    hide_top_spine: bool = True
    hide_right_spine: bool = True

    # Title style
    title_fontweight: str = "bold"

    # Legend placement
    legend_loc: str = "center left"
    legend_bbox: tuple[float, float] = (1.02, 0.5)

    # Watermark
    show_watermark: bool = True

    @classmethod
    def publication(cls) -> PlotTheme:
        """Publication preset — default values, print-ready."""
        return cls()

    @classmethod
    def presentation(cls) -> PlotTheme:
        """Presentation preset — ~1.3x larger fonts/dots/lines for slides."""
        return cls(
            title_fontsize=18,
            suptitle_fontsize=20,
            label_fontsize=15,
            tick_fontsize=12,
            legend_fontsize=12,
            annotation_fontsize=12,
            small_fontsize=10,
            tiny_fontsize=9,
            dot_size=30,
            bar_linewidth=0.8,
            bar_capsize=5,
            reference_line_width=2.0,
            fill_alpha=0.3,
        )

    @classmethod
    def minimal(cls) -> PlotTheme:
        """Minimal preset — no dots, no bar edges, thinner lines."""
        return cls(
            dot_size=0,
            dot_alpha=0.0,
            bar_edgecolor="none",
            bar_linewidth=0.0,
            bar_capsize=3,
            reference_line_width=1.0,
            fill_alpha=0.15,
        )


class PlotSettings(BaseModel):
    """Global plot settings for comparison.yaml.

    Controls plot generation for all analyses. Per-analysis plot settings
    are discovered via ``PlotSettingsRegistry`` — any key in the YAML that
    matches a registered analysis type is parsed into the corresponding
    settings class.  Unrecognised keys that are not global fields are
    logged and skipped.

    Attributes
    ----------
    output_dir : Path
        Directory for generated plots (relative to comparison.yaml)
    format : str
        Image format: "png", "pdf", or "svg"
    dpi : int
        Resolution for raster formats (PNG)
    style : str
        Plot style preset: "publication", "presentation", or "minimal"
    color_palette : str
        Seaborn/matplotlib color palette name
    theme : PlotTheme
        Resolved visual theme.  Built from the ``style`` preset and
        any user overrides in the ``theme:`` YAML block.

    Notes
    -----
    Attribute access for any registered analysis type always succeeds:
    if the user did not provide that section in YAML, a default-constructed
    settings instance is returned.  This means ``self.settings.rmsf.show_error``
    is always safe, even when the YAML has no ``rmsf:`` block.

    Examples
    --------
    In comparison.yaml:

    .. code-block:: yaml

        plot_settings:
          output_dir: "figures/"
          format: "png"
          dpi: 300
          style: "publication"

          rmsf:
            highlight_residues: [77, 133, 156]

          triad:
            generate_2d_kde: true
    """

    model_config = {"extra": "allow"}

    _GLOBAL_FIELDS: ClassVar[set[str]] = {
        "output_dir",
        "format",
        "dpi",
        "style",
        "color_palette",
        "theme",
    }

    output_dir: Path = Field(default=Path("figures/"))
    format: str = Field(default="png", pattern="^(png|pdf|svg)$")
    dpi: int = Field(default=300, ge=50, le=600)
    style: str = Field(default="publication", pattern="^(publication|presentation|minimal)$")
    color_palette: str = "tab10"
    theme: PlotTheme = Field(default_factory=PlotTheme)

    def __init__(self, **data: Any):
        """Initialize with global fields and registry-discovered per-analysis settings.

        Theme resolution: the ``style`` field selects a preset (publication,
        presentation, or minimal) and then any user-supplied ``theme:``
        overrides are merged on top.  This allows ``style: presentation``
        with ``theme: {dot_size: 40}`` to use the presentation preset but
        override just the dot size.

        Parameters
        ----------
        **data : Any
            Plot settings from YAML.  Keys matching registered analysis
            types are parsed into their settings classes; global keys are
            handled by Pydantic; unknown keys are logged and skipped.
        """
        global_data: dict[str, Any] = {}
        per_analysis: dict[str, BasePlotSettings] = {}

        for key, value in data.items():
            if key in PlotSettings._GLOBAL_FIELDS:
                global_data[key] = value
            elif PlotSettingsRegistry.is_registered(key):
                settings_class = PlotSettingsRegistry.get(key)
                if isinstance(value, dict):
                    per_analysis[key] = settings_class(**value)
                elif isinstance(value, BasePlotSettings):
                    per_analysis[key] = value
                else:
                    raise ValueError(
                        f"Invalid value for plot settings '{key}': "
                        f"expected dict or {settings_class.__name__}"
                    )
            else:
                logger.warning(f"Unknown plot settings key '{key}' — skipping")

        # ── Resolve theme from style preset + user overrides ──
        style = global_data.get("style", "publication")
        theme_overrides = global_data.pop("theme", None)

        _THEME_PRESETS = {
            "publication": PlotTheme.publication,
            "presentation": PlotTheme.presentation,
            "minimal": PlotTheme.minimal,
        }
        preset_factory = _THEME_PRESETS.get(style, PlotTheme.publication)

        if theme_overrides is None or (isinstance(theme_overrides, dict) and not theme_overrides):
            # No user overrides — use preset as-is
            global_data["theme"] = preset_factory()
        elif isinstance(theme_overrides, dict):
            # Merge user overrides on top of preset defaults
            preset = preset_factory()
            merged = {**preset.model_dump(), **theme_overrides}
            global_data["theme"] = PlotTheme(**merged)
        elif isinstance(theme_overrides, PlotTheme):
            # Already a PlotTheme instance (programmatic usage)
            global_data["theme"] = theme_overrides

        super().__init__(**global_data, **per_analysis)

    def __getattr__(self, name: str) -> Any:
        """Fall back to default-constructed settings for registered types.

        This ensures ``self.settings.rmsf.show_error`` works even when
        the user omitted the ``rmsf:`` block from their YAML.

        Parameters
        ----------
        name : str
            Attribute name.

        Returns
        -------
        BasePlotSettings
            Default-constructed settings if *name* is a registered type.

        Raises
        ------
        AttributeError
            If *name* is not a registered plot settings type.
        """
        if PlotSettingsRegistry.is_registered(name):
            settings_class = PlotSettingsRegistry.get(name)
            return settings_class()
        raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'")

    @field_validator("output_dir", mode="before")
    @classmethod
    def resolve_output_dir(cls, v: str | Path) -> Path:
        """Convert string paths to Path objects."""
        return Path(v)


# ============================================================================
# Main Comparison Configuration
# ============================================================================


class ComparisonConfig(BaseModel):
    """Schema for comparison.yaml configuration files.

    A comparison config defines multiple simulation conditions to compare,
    along with analysis settings and comparison-specific parameters.

    The schema follows a three-section pattern:
    - analysis_settings: WHAT to analyze (shared across conditions)
    - comparison_settings: HOW to compare (statistical parameters)
    - plot_settings: HOW to visualize (plot customization)

    Attributes
    ----------
    name : str
        Name of the comparison project
    description : str, optional
        Description of what is being compared
    control : str, optional
        Label of the control condition for relative comparisons
    conditions : list[ConditionConfig]
        List of conditions to compare
    defaults : AnalysisDefaults
        Default analysis parameters (equilibration_time)
    analysis_settings : AnalysisSettingsContainer
        Analysis parameters (WHAT to analyze)
    comparison_settings : ComparisonSettingsContainer
        Comparison parameters (HOW to compare)
    plot_settings : PlotSettings
        Plot customization (HOW to visualize)

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> print(config.name)
    "Polymer Stabilization Study"
    >>> for cond in config.conditions:
    ...     print(f"{cond.label}: {cond.config}")
    >>> print("Enabled analyses:", config.analysis_settings.get_enabled_analyses())
    >>> rmsf_settings = config.analysis_settings.get("rmsf")
    >>> if rmsf_settings:
    ...     print(f"RMSF selection: {rmsf_settings.selection}")
    """

    name: str
    description: str | None = None
    control: str | None = None
    conditions: list[ConditionConfig]
    defaults: AnalysisDefaults = Field(default_factory=AnalysisDefaults)
    analysis_settings: AnalysisSettingsContainer = Field(default_factory=AnalysisSettingsContainer)
    comparison_settings: ComparisonSettingsContainer = Field(
        default_factory=ComparisonSettingsContainer
    )
    plot_settings: PlotSettings = Field(default_factory=PlotSettings)
    source_path: Path | None = Field(default=None, exclude=True)

    @field_validator("analysis_settings", mode="before")
    @classmethod
    def parse_analysis_settings(cls, v: Any) -> AnalysisSettingsContainer:
        """Parse analysis_settings from dict or container."""
        if v is None:
            return AnalysisSettingsContainer()
        if isinstance(v, dict):
            return AnalysisSettingsContainer(**v)
        return v

    @field_validator("comparison_settings", mode="before")
    @classmethod
    def parse_comparison_settings(cls, v: Any) -> ComparisonSettingsContainer:
        """Parse comparison_settings from dict or container."""
        if v is None:
            return ComparisonSettingsContainer()
        if isinstance(v, dict):
            return ComparisonSettingsContainer(**v)
        return v

    @model_validator(mode="after")
    def validate_comparison_coverage(self) -> "ComparisonConfig":
        """Validate that comparison_settings covers all analysis_settings.

        Each analysis type in analysis_settings must have a corresponding
        entry in comparison_settings (can be empty {}).
        """
        enabled_analyses = self.analysis_settings.get_enabled_analyses()
        enabled_comparisons = self.comparison_settings.get_enabled_comparisons()

        missing = set(enabled_analyses) - set(enabled_comparisons)
        if missing:
            raise ValueError(
                f"Missing comparison_settings for: {sorted(missing)}. "
                f"Add 'comparison_settings.{list(missing)[0]}: {{}}' to enable comparison."
            )
        return self

    @classmethod
    def from_yaml(cls, path: Path | str) -> "ComparisonConfig":
        """Load comparison config from YAML file.

        Parameters
        ----------
        path : Path or str
            Path to comparison.yaml file

        Returns
        -------
        ComparisonConfig
            Loaded and validated configuration

        Raises
        ------
        FileNotFoundError
            If the config file doesn't exist
        ValidationError
            If the config is invalid
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Comparison config not found: {path}")

        with open(path) as f:
            data = yaml.safe_load(f)

        # Resolve relative paths relative to the config file location
        config_dir = path.parent.resolve()
        if "conditions" in data:
            for cond in data["conditions"]:
                if "config" in cond:
                    cond_path = Path(cond["config"])
                    if not cond_path.is_absolute():
                        cond["config"] = str(config_dir / cond_path)

        config = cls(**data)
        config.source_path = path.resolve()
        return config

    def to_yaml(self, path: Path | str) -> None:
        """Save comparison config to YAML file.

        Parameters
        ----------
        path : Path or str
            Output path for comparison.yaml
        """
        path = Path(path)

        # Convert to dict, handling Path objects and nested containers
        data = self.model_dump(mode="json")
        for cond in data["conditions"]:
            cond["config"] = str(cond["config"])

        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def get_condition(self, label: str) -> ConditionConfig:
        """Get a condition by its label.

        Parameters
        ----------
        label : str
            The condition label to find

        Returns
        -------
        ConditionConfig
            The matching condition

        Raises
        ------
        KeyError
            If no condition with that label exists
        """
        for cond in self.conditions:
            if cond.label == label:
                return cond
        raise KeyError(f"Condition '{label}' not found in: {[c.label for c in self.conditions]}")

    def validate_config(self) -> list[str]:
        """Validate the comparison configuration.

        Returns
        -------
        list[str]
            List of error messages (empty if valid)
        """
        errors = []

        # Check minimum conditions
        if len(self.conditions) < 2:
            errors.append("Need at least 2 conditions to compare")

        # Check for duplicate labels
        labels = [c.label for c in self.conditions]
        if len(labels) != len(set(labels)):
            errors.append("Duplicate condition labels found")

        # Check control label exists
        if self.control and self.control not in labels:
            errors.append(f"Control '{self.control}' not in conditions: {labels}")

        # Check config files exist
        for cond in self.conditions:
            if not cond.config.exists():
                errors.append(f"Config not found for '{cond.label}': {cond.config}")

        # Check analysis/comparison coverage
        enabled_analyses = self.analysis_settings.get_enabled_analyses()
        enabled_comparisons = self.comparison_settings.get_enabled_comparisons()
        missing = set(enabled_analyses) - set(enabled_comparisons)
        if missing:
            errors.append(
                f"Missing comparison_settings for: {sorted(missing)}. "
                f"Add comparison_settings entries for these analyses."
            )

        return errors

    def generate_analysis_yaml(self, condition: ConditionConfig) -> str:
        """Generate analysis.yaml content for a specific condition.

        Parameters
        ----------
        condition : ConditionConfig
            The condition to generate analysis.yaml for.

        Returns
        -------
        str
            YAML content for the analysis.yaml file.
        """
        data = self.analysis_settings.to_analysis_yaml_dict(
            replicates=condition.replicates,
            eq_time=self.defaults.equilibration_time,
        )
        return yaml.dump(data, default_flow_style=False, sort_keys=False)

    def generate_analysis_yaml_for_all(self) -> dict[str, str]:
        """Generate analysis.yaml content for all conditions.

        Returns
        -------
        dict[str, str]
            Dictionary mapping condition labels to analysis.yaml content.
        """
        return {cond.label: self.generate_analysis_yaml(cond) for cond in self.conditions}


def generate_comparison_template(name: str, eq_time: str = "10ns") -> str:
    """Generate a template comparison.yaml file.

    Parameters
    ----------
    name : str
        Project name
    eq_time : str
        Default equilibration time

    Returns
    -------
    str
        YAML template content
    """
    return f"""\
# ============================================================================
# PolyzyMD Comparison Configuration
# ============================================================================
# Compare analyses across simulation conditions (e.g., polymer, temperature).
# Docs: https://polyzymd.readthedocs.io/en/latest/
# ============================================================================

name: "{name}"
description: "Comparison of simulation conditions"

# Control condition for relative comparisons.
# Must match one of the 'label' values in 'conditions' below, or null if none.
control: null

# ============================================================================
# Conditions
# ============================================================================
# Each condition points to a simulation's config.yaml file.
conditions:
  - label: "Condition A"
    config: "../path/to/condition_a/config.yaml"
    replicates: [1, 2, 3]

  - label: "Condition B"
    config: "../path/to/condition_b/config.yaml"
    replicates: [1, 2, 3]

# ============================================================================
# Defaults
# ============================================================================
defaults:
  equilibration_time: "{eq_time}"

# ============================================================================
# Analysis Settings (WHAT to analyze - applied to all conditions)
# ============================================================================
# Define which analyses to run. Presence of a section enables that analysis.
# Running `polyzymd compare analyze` will run these for each condition.

analysis_settings:
  # RMSF Analysis
  rmsf:
    selection: "protein and name CA"
    reference_mode: "centroid"  # centroid, average, or frame
    # reference_frame: 500      # Required if reference_mode is "frame"
    # reference_file: "structures/enzyme.pdb"  # Crystal/input PDB for SS annotation bar

  # Secondary Structure (DSSP) Analysis
  # Per-residue and per-frame secondary structure via mdtraj DSSP.
  # Produces timeline heatmaps, content bar charts, and persistence profiles.
  #
  # secondary_structure:
  #   chain_id: "A"              # chain letter for the protein to analyze

  # Catalytic Triad / Active Site Distances
  #
  # IMPORTANT: Always use "protein and resid X" for protein residues!
  # Residue numbers restart per chain. Without "protein and", your selection
  # may match atoms from polymer or water chains, causing incorrect distances.
  #
  # catalytic_triad:
  #   name: "enzyme_catalytic_triad"
  #   threshold: 3.5  # Angstroms (H-bond cutoff)
  #   pairs:
  #     - label: "Asp-His"
  #       selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
  #       selection_b: "protein and resid 156 and name ND1"
  #     - label: "His-Ser"
  #       selection_a: "protein and resid 156 and name NE2"
  #       selection_b: "protein and resid 77 and name OG"

  # Distance Analysis (general inter-atomic distances)
  #
  # IMPORTANT: Always use "protein and resid X" for protein residues!
  # See warning above in catalytic_triad section.
  #
  # Each pair can have its own threshold and display labels for above/below
  # states.  If threshold is omitted, the global threshold (default 3.5 Å)
  # is used.  If below_label / above_label are omitted, defaults are
  # "Below {threshold}Å" / "Above {threshold}Å".
  #
  # distances:
  #   threshold: 3.5                         # global default threshold
  #   pairs:
  #     - label: "Ser77-Substrate"
  #       selection_a: "protein and resid 77 and name OG"
  #       selection_b: "resname RBY and name C1"
  #       threshold: 3.5                     # per-pair override (optional)
  #       below_label: "Bound"               # d <= threshold
  #       above_label: "Unbound"             # d > threshold
  #     - label: "Lid Domain"
  #       selection_a: "com(resid 141:148 and chainID A)"
  #       selection_b: "com(resid 281:289 and chainID A)"
  #       threshold: 15.0
  #       below_label: "Closed"
  #       above_label: "Open"

  # Polymer-Protein Contact Analysis
  # contacts:
  #   polymer_selection: "chainID C"
  #   protein_selection: "protein"
  #   cutoff: 4.5
  #   grouping: "aa_class"  # aa_class, secondary_structure, or none
  #   compute_residence_times: true
  #
  #   # EXPERIMENTAL: Binding Preference Analysis (enrichment by residue group)
  #   # --------------------------------------------------------
  #   # Computes which residue types (aromatic, polar, etc.) are preferentially
  #   # contacted by the polymer, normalized by surface exposure.
  #   # Definitions and interpretation may change after the presentation release.
  #   #
  #   # IMPORTANT: Place your enzyme PDB in the structures/ directory!
  #   # The path is relative to this comparison.yaml file.
  #   #
  #   # compute_binding_preference: true
  #   # surface_exposure_threshold: 0.2  # 20% relative SASA = surface exposed
  #   # enzyme_pdb_for_sasa: "structures/enzyme.pdb"
  #   # include_default_aa_groups: true  # aromatic, polar, nonpolar, charged
  #   #
  #   # Custom protein groups (residue IDs from your enzyme):
  #   # protein_groups:
  #   #   catalytic_triad: [77, 133, 156]
  #   #   lid_helix_5: [141, 142, 143, 144, 145]
  #   #   lid_helix_10: [281, 282, 283, 284, 285]
  #   #
  #   # User-defined partitions for system coverage plots
  #   # -------------------------------------------------
  #   # Groups within a partition MUST be mutually exclusive (no overlapping
  #   # residues). A "rest_of_protein" element is automatically added if the
  #   # groups don't cover all surface-exposed residues. One plot is generated
  #   # per partition.
  #   #
  #   # protein_partitions:
  #   #   lid_helices:          # partition name (becomes plot title)
  #   #     - lid_helix_5       # must be defined in protein_groups above
  #   #     - lid_helix_10

  # EXPERIMENTAL: Exposure Dynamics Analysis (chaperone-like polymer activity)
  # Requires contacts analysis to be run first for each condition.
  # Definitions and interpretation may change after the presentation release.
  # Run: polyzymd compare exposure
  #
  # exposure:
  #   exposure_threshold: 0.20     # fraction SASA defining 'exposed' residue
  #   transient_lower: 0.20        # lower bound: residue must reach this to be transient
  #   transient_upper: 0.80        # upper bound: residue must also reach this threshold
  #   min_event_length: 1          # minimum consecutive frames to count as an event
  #   protein_chain: "A"           # chain ID for the protein
  #   protein_selection: "protein" # MDAnalysis selection for protein
  #   polymer_selection: "chainID C"  # MDAnalysis selection for polymer
  #   # polymer_resnames: [SBM, EGM]  # optional: residue names for enrichment analysis
  #   probe_radius_nm: 0.14        # SASA probe radius (nm)
  #   n_sphere_points: 960         # number of sphere points for SASA computation

  # EXPERIMENTAL: Binding Free Energy Analysis (ΔG_sel via Boltzmann inversion)
  # Requires contacts analysis with compute_binding_preference: true to be run first.
  # Converts binding preference probabilities into ΔG_sel = -k_B·T·ln(contact_share / expected_share).
  # Definitions and interpretation may change after the presentation release.
  # Run: polyzymd compare binding-free-energy
  #
  # binding_free_energy:
  #   units: kT                    # energy units: kT (default), kcal/mol, or kJ/mol
  #   surface_exposure_threshold: 0.2  # minimum relative SASA to be considered surface-exposed
  #   # protein_partitions: null   # optional: restrict to user-defined AA partitions

  # EXPERIMENTAL: Polymer Affinity Score (composite selectivity metric)
  # Requires contacts analysis with compute_binding_preference: true.
  # Computes S = Σ N_pg × ΔG_sel_pg for each polymer composition.
  # Definitions and interpretation may change after the presentation release.
  # Run: polyzymd compare polymer-affinity
  #
  # polymer_affinity:
  #   surface_exposure_threshold: 0.2
  #   enzyme_pdb_for_sasa: "structures/enzyme.pdb"
  #   include_default_aa_groups: true
  #   # protein_groups:            # same format as contacts.protein_groups
  #   #   catalytic_triad: [77, 133, 156]
  #   # protein_partitions:        # same format as contacts.protein_partitions
  #   #   lid_helices:
  #   #     - lid_helix_5

# ============================================================================
# Comparison Settings (HOW to compare - statistical parameters)
# ============================================================================
# Each analysis in analysis_settings MUST have a corresponding entry here.
# Use empty {{}} for analyses with no comparison-specific parameters.

comparison_settings:
  rmsf: {{}}  # No comparison-specific parameters

  # secondary_structure: {{}}     # No comparison-specific parameters

  # catalytic_triad: {{}}

  # distances: {{}}

  # contacts:
  #   fdr_alpha: 0.05           # FDR for Benjamini-Hochberg correction
  #   min_effect_size: 0.5      # Cohen's d threshold (0.2=small, 0.5=medium, 0.8=large)
  #   top_residues: 10          # Number of top residues to show in console

  # exposure: {{}}              # No comparison-specific parameters for exposure

  # binding_free_energy:
  #   fdr_alpha: 0.05           # FDR for Benjamini-Hochberg correction

  # polymer_affinity:
  #   fdr_alpha: 0.05           # FDR for Benjamini-Hochberg correction

# ============================================================================
# Plot Settings (HOW to visualize - figure customization)
# ============================================================================
# Controls plot generation for all analyses. Per-analysis sections override
# defaults. Run `polyzymd compare plot` to generate all configured plots.

plot_settings:
  output_dir: "figures/"         # relative to this file
  format: "png"                  # png, pdf, or svg
  dpi: 300                       # resolution for raster formats
  style: "publication"           # publication, presentation, or minimal
  color_palette: "tab10"         # seaborn/matplotlib color palette

  # ── Visual theme (all fields optional — defaults come from style preset) ──
  # Uncomment individual lines to override the preset values.
  # theme:
  #   # Font sizes by semantic role
  #   title_fontsize: 13           # axes titles
  #   suptitle_fontsize: 14        # figure suptitles
  #   label_fontsize: 11           # axis labels (xlabel/ylabel)
  #   tick_fontsize: 9             # tick labels
  #   legend_fontsize: 9           # legend entries
  #   annotation_fontsize: 9       # heatmap cell annotations
  #   small_fontsize: 8            # secondary annotations (SEM labels)
  #   tiny_fontsize: 7             # fine-grained annotations (residue IDs)
  #
  #   # Bar chart defaults
  #   bar_alpha: 0.85              # bar fill opacity
  #   bar_edgecolor: "black"       # bar edge colour
  #   bar_linewidth: 0.5           # bar edge line width
  #   bar_capsize: 4               # error bar cap size (points)
  #
  #   # Replicate dot overlay
  #   dot_size: 18                 # scatter marker size (s=)
  #   dot_alpha: 0.7               # dot opacity
  #   dot_color: "black"           # dot colour
  #
  #   # Line defaults
  #   line_alpha: 0.8              # line plot opacity
  #   fill_alpha: 0.25             # fill_between band opacity
  #   reference_line_color: "black"   # reference/threshold line colour
  #   reference_line_style: "--"      # reference line style
  #   reference_line_width: 1.5       # reference line width
  #   highlight_line_alpha: 0.5       # vertical highlight line opacity
  #
  #   # Axes chrome
  #   hide_top_spine: true         # hide top axis spine
  #   hide_right_spine: true       # hide right axis spine
  #
  #   # Title style
  #   title_fontweight: "bold"     # title font weight
  #
  #   # Legend placement
  #   legend_loc: "center left"           # matplotlib legend loc string
  #   legend_bbox: [1.02, 0.5]            # bbox_to_anchor (outside right)

  # Per-analysis plot customization (uncomment sections as needed):

  # rmsf:
  #   show_error: true             # show SEM fill_between bands
  #   highlight_residues: []       # residue IDs for vertical reference lines
  #   figsize_profile: [14, 4]     # per-residue profile figure size
  #   figsize_comparison: [8, 6]   # bar comparison figure size

  # triad:
  #   generate_kde_panel: true     # multi-row KDE panel
  #   generate_bars: true          # threshold bar chart
  #   generate_2d_kde: false       # 2D joint KDE (specialized)
  #   kde_xlim: [0, 7]             # x-axis range for KDE panel (Angstroms)

  # distances:
  #   show_threshold: true         # threshold line on distributions
  #   use_kde: true                # KDE vs histogram
  #   generate_state_bars: true    # per-pair above/below threshold bars

  # contacts:
  #   generate_enrichment_heatmap: true
  #   generate_enrichment_bars: true
  #   generate_system_coverage_heatmap: true
  #   generate_system_coverage_bars: true
  #   generate_contact_fraction_profile: true
  #   generate_residence_time_profile: true

  # binding_free_energy:
  #   generate_heatmap: true       # ΔG_sel heatmap (AA groups × conditions)
  #   generate_bars: true          # ΔG_sel grouped bar chart
  #   colormap: "RdBu_r"           # diverging colormap for heatmap

  # polymer_affinity:
  #   generate_stacked_bars: true  # total score by condition
  #   generate_group_bars: true    # per-group contributions

  # secondary_structure:
  #   generate_timeline: true      # per-condition residue × time SS heatmap
  #   generate_content_bars: true  # helix/strand/coil fraction bars
  #   generate_individual_bars: true  # one bar chart per SS type
  #   generate_diff_heatmap: true  # Δ(helix persistence) vs control
  #   diff_colormap: "RdBu_r"     # diverging colormap for diff heatmap
"""
