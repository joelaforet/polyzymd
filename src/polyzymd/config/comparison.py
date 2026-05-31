"""Comparison configuration schema in ``polyzymd.config``."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, ClassVar, Literal

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import BasePlotSettings
from polyzymd.core.archived_features import (
    format_archived_analysis_message,
    get_archived_analysis_plugin,
)
from polyzymd.core.branding import prepend_file_header
from polyzymd.utils.templates import render_package_template

logger = logging.getLogger(__name__)

CANONICAL_PLOT_STYLES: tuple[str, ...] = ("compact", "large_elements", "low_ink")
_PLOT_STYLE_ALLOWED_MESSAGE = "allowed values are 'compact', 'large_elements', and 'low_ink'"


def _normalize_plot_style(value: Any) -> str:
    """Validate plot style presets before Pydantic stores them.

    Parameters
    ----------
    value : Any
        Raw ``plot_settings.style`` value from YAML or programmatic usage.

    Returns
    -------
    str
        Canonical style name.

    Raises
    ------
    TypeError
        If the style is not a string.
    ValueError
        If the style string is not a canonical style.
    """
    if not isinstance(value, str):
        raise TypeError(
            f"plot_settings.style must be a string, got {type(value).__name__}; "
            f"{_PLOT_STYLE_ALLOWED_MESSAGE}"
        )

    if value in CANONICAL_PLOT_STYLES:
        return value

    raise ValueError(f"Invalid plot_settings.style '{value}'; {_PLOT_STYLE_ALLOWED_MESSAGE}")


class AnalysisDefaults(BaseModel):
    """Default parameters applied to every configured analysis.

    The comparison config owns these project-level defaults because they are
    parsed from ``comparison.yaml`` before plugin lifecycle contexts are built.
    """

    equilibration_time: str = "10ns"
    fdr_alpha: float = Field(default=0.05, gt=0.0, le=1.0)
    ttest_method: str = Field(default="student", pattern="^(welch|student)$")
    posthoc_method: str = Field(default="ttest_bh", pattern="^(ttest_bh|tukey_hsd)$")


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

    @model_validator(mode="after")
    def _validate_replicates(self) -> "ConditionConfig":
        """Reject empty, duplicate, or non-positive replicate lists."""
        if not self.replicates:
            raise ValueError(f"Condition '{self.label}': replicates list must not be empty")
        if any(r <= 0 for r in self.replicates):
            bad = [r for r in self.replicates if r <= 0]
            raise ValueError(
                f"Condition '{self.label}': replicate numbers must be positive, got {bad}"
            )
        if len(self.replicates) != len(set(self.replicates)):
            dupes = [r for r in self.replicates if self.replicates.count(r) > 1]
            raise ValueError(
                f"Condition '{self.label}': duplicate replicate numbers: {sorted(set(dupes))}"
            )
        return self


# ============================================================================
# Dynamic Settings Containers
# ============================================================================


class PluginSettingsContainer(BaseModel):
    """Container for unified plugin settings keyed by analysis name."""

    model_config = {"extra": "allow"}

    def __init__(self, **data: Any):
        """Initialize plugin settings using discovered analysis classes."""
        from polyzymd.analyses.discovery import get_analysis

        parsed_settings: dict[str, Any] = {}
        for key, value in data.items():
            key_lower = str(key).lower()
            if get_archived_analysis_plugin(key_lower) is not None:
                raise ValueError(format_archived_analysis_message(key, context="plugins section"))
            if value is None:
                continue
            try:
                analysis_cls = get_analysis(key_lower)
            except KeyError:
                from polyzymd.analyses.discovery import list_analyses

                available = sorted(list_analyses().keys())
                raise ValueError(
                    f"Unknown analysis plugin '{key}' in plugins section. "
                    f"Available plugins: {available}"
                ) from None

            # Always store under the canonical name so that
            # _resolve_settings() (which queries by analysis.name) finds it,
            # even when the YAML key is an alias like "triad".
            canonical_name = analysis_cls.name

            settings_class = analysis_cls.Settings
            if isinstance(value, settings_class):
                parsed_settings[canonical_name] = value
            elif isinstance(value, dict):
                parsed_settings[canonical_name] = settings_class.model_validate(value)
            elif hasattr(value, "model_dump"):
                parsed_settings[canonical_name] = settings_class.model_validate(value.model_dump())
            else:
                raise ValueError(
                    f"Invalid value for {key}: expected dict or {settings_class.__name__}"
                )

        super().__init__(**parsed_settings)

    def get(self, analysis_type: str) -> Any | None:
        """Get settings for a specific analysis type."""
        return getattr(self, analysis_type.lower(), None)

    def get_enabled_plugins(self) -> list[str]:
        """Get the configured plugin names."""
        return [key for key, value in self.model_dump().items() if value is not None]

    def to_analysis_yaml_dict(self, replicates: list[int], eq_time: str) -> dict[str, Any]:
        """Convert unified plugin settings to analysis.yaml-compatible data."""
        result: dict[str, Any] = {
            "replicates": replicates,
            "defaults": {"equilibration_time": eq_time},
        }
        for analysis_type in self.get_enabled_plugins():
            settings = self.get(analysis_type)
            if settings is not None:
                if hasattr(settings, "to_analysis_yaml_dict"):
                    result[analysis_type] = settings.to_analysis_yaml_dict()
                elif hasattr(settings, "model_dump"):
                    result[analysis_type] = settings.model_dump(exclude_none=True)
                else:
                    result[analysis_type] = settings
        return result


# ============================================================================
# Plot Settings Configuration
# ============================================================================

# Per-analysis plot settings classes (RMSFPlotSettings, TriadPlotSettings, etc.)
# live in their respective plugin packages at analyses/<name>/_plot_settings.py.
# Each plugin exposes its plot settings model via
# Analysis.PlotSettingsModel. PlotSettings.__init__ discovers plugins and
# builds a mapping from analysis name to plot settings model.


class SemanticConditionColorConfig(BaseModel):
    """Semantic metadata for one plotted condition.

    The model is intentionally chemistry-agnostic. Projects can describe any
    condition family, numeric or ordinal value, display order, direct color,
    or role without core PolyzyMD knowing domain-specific condition names.
    """

    color: Any | None = None
    family: str | None = None
    value: Any | None = None
    order: int | None = None
    role: str | None = None

    @field_validator("family")
    @classmethod
    def strip_optional_family(cls, value: str | None) -> str | None:
        """Normalize optional family text while rejecting empty values."""
        if value is None:
            return None
        stripped = value.strip()
        if not stripped:
            raise ValueError("semantic condition family must not be empty")
        return stripped

    @field_validator("role")
    @classmethod
    def strip_optional_role(cls, value: str | None) -> str | None:
        """Normalize optional role text while rejecting empty values."""
        if value is None:
            return None
        stripped = value.strip().lower()
        if not stripped:
            raise ValueError("semantic condition role must not be empty")
        return stripped


class SemanticFamilyColorConfig(BaseModel):
    """Color mapping rules for a semantic family of conditions.

    A family can map condition values through a matplotlib colormap using a
    linear or ordinal scale, or through explicit ``value_colors`` for selected
    values. Invalid color and colormap names are handled by plotting helpers so
    figure generation can fall back with warnings.
    """

    colormap: str = "viridis"
    scale: Literal["linear", "ordinal"] = "linear"
    value_order: list[Any] = Field(default_factory=list)
    vmin: float | None = None
    vmax: float | None = None
    colormap_range: tuple[float, float] = (0.0, 1.0)
    reverse: bool = False
    value_colors: dict[str, Any] = Field(default_factory=dict)

    @field_validator("colormap")
    @classmethod
    def validate_colormap_name(cls, value: str) -> str:
        """Reject empty colormap names while deferring existence checks."""
        colormap = value.strip()
        if not colormap:
            raise ValueError("colormap must not be empty")
        return colormap

    @field_validator("colormap_range")
    @classmethod
    def validate_colormap_range(cls, value: tuple[float, float]) -> tuple[float, float]:
        """Ensure colormap sampling bounds are ordered fractions."""
        low, high = value
        if not 0.0 <= low <= 1.0 or not 0.0 <= high <= 1.0:
            raise ValueError("colormap_range values must be between 0.0 and 1.0")
        if low > high:
            raise ValueError("colormap_range lower bound must be <= upper bound")
        return value

    @model_validator(mode="after")
    def validate_linear_bounds(self) -> "SemanticFamilyColorConfig":
        """Validate optional linear scale bounds."""
        if self.vmin is not None and self.vmax is not None and self.vmin > self.vmax:
            raise ValueError("vmin must be <= vmax")
        return self


class SemanticColorSettings(BaseModel):
    """Semantic condition color and ordering settings for comparison plots.

    Defaults keep semantic behavior disabled, preserving legacy palette order
    and colors until users opt in with ``enabled: true``.
    """

    enabled: bool = False
    order: list[str] = Field(default_factory=list)
    conditions: dict[str, SemanticConditionColorConfig] = Field(default_factory=dict)
    families: dict[str, SemanticFamilyColorConfig] = Field(default_factory=dict)
    manual_colors: dict[str, Any] = Field(default_factory=dict)
    control_color: Any = "black"
    missing_color: Any = "lightgray"
    default_color: Any | None = None

    @field_validator("order")
    @classmethod
    def validate_order_labels(cls, value: list[str]) -> list[str]:
        """Reject duplicate or empty labels in explicit plot order."""
        normalized: list[str] = []
        for label in value:
            stripped = label.strip()
            if not stripped:
                raise ValueError("semantic color order labels must not be empty")
            normalized.append(stripped)
        if len(normalized) != len(set(normalized)):
            raise ValueError("semantic color order labels must be unique")
        return normalized


class PlotTheme(BaseModel):
    """Centralized visual defaults for all comparison plots.

    Replaces ~219 hardcoded style values (font sizes, alphas, line widths,
    marker sizes, spine visibility, etc.) across all plugin ``plot()`` methods
    with a single configurable Pydantic model.

    Canonical presets are available via class methods:

    - ``PlotTheme.compact()`` — default; print-ready sizes and weights.
    - ``PlotTheme.large_elements()`` — ~1.3x larger fonts/dots/lines for slides.
    - ``PlotTheme.low_ink()`` — no dots, no bar edges, thinner lines.

    Users can override individual values in ``comparison.yaml``::

        plot_settings:
          style: "compact"
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
    def compact(cls) -> PlotTheme:
        """Compact preset with print-ready sizes and weights."""
        return cls()

    @classmethod
    def large_elements(cls) -> PlotTheme:
        """Large-elements preset with slide-oriented fonts, dots, and lines."""
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
    def low_ink(cls) -> PlotTheme:
        """Low-ink preset with no dots, no bar edges, and thinner lines."""
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
    are discovered via analysis plugin metadata — any key in the YAML that
    matches a discovered analysis name with a ``PlotSettingsModel`` is parsed
    into the corresponding settings class. Unrecognized keys that are not
    global fields raise ``ValueError``.

    Attributes
    ----------
    output_dir : Path
        Directory for generated plots (relative to comparison.yaml)
    format : str
        Image format: "png", "pdf", or "svg"
    dpi : int
        Resolution for raster formats (PNG)
    style : str
        Canonical plot style preset: "compact", "large_elements", or
        "low_ink".
    color_palette : str
        Seaborn/matplotlib color palette name
    theme : PlotTheme
        Resolved visual theme.  Built from the ``style`` preset and
        any user overrides in the ``theme:`` YAML block.

    Notes
    -----
    Attribute access for any discovered analysis type with a
    ``PlotSettingsModel`` always succeeds:
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
          style: "compact"

          rmsf:
            highlight_residues: [77, 133, 156]

          catalytic_triad:
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
        "semantic_colors",
    }

    output_dir: Path = Field(default=Path("figures/"))
    format: str = Field(default="png", pattern="^(png|pdf|svg)$")
    dpi: int = Field(default=300, ge=50, le=600)
    style: str = "compact"
    color_palette: str = "tab10"
    theme: PlotTheme = Field(default_factory=PlotTheme)
    semantic_colors: SemanticColorSettings = Field(default_factory=SemanticColorSettings)

    def __init__(self, **data: Any):
        """Initialize with global fields and plugin-discovered per-analysis settings.

        Theme resolution: the ``style`` field selects a canonical preset
        (compact, large_elements, or low_ink) and then any user-supplied
        ``theme:`` overrides are merged on top. This allows
        ``style: large_elements`` with ``theme: {dot_size: 40}`` to use
        the large-elements preset but override just the dot size.

        Parameters
        ----------
        **data : Any
            Plot settings from YAML.  Keys matching discovered analysis
            types are parsed into their settings classes; global keys are
            handled by Pydantic; unknown keys are logged and skipped.
        """
        from polyzymd.analyses.discovery import list_analyses

        plugin_registry = list_analyses()
        _plot_models: dict[str, type[BasePlotSettings]] = {}
        for analysis_name, analysis_cls in plugin_registry.items():
            if (
                hasattr(analysis_cls, "PlotSettingsModel")
                and analysis_cls.PlotSettingsModel is not None
            ):
                _plot_models[analysis_name] = analysis_cls.PlotSettingsModel

        # Check for deprecated "triad" key and give helpful warning
        if "triad" in data and "triad" not in _plot_models:
            import warnings

            warnings.warn(
                "plot_settings key 'triad' has been renamed to 'catalytic_triad'. "
                "Please update your comparison.yaml to use 'catalytic_triad' instead.",
                FutureWarning,
                stacklevel=2,
            )
            # Do NOT silently remap — raise the normal unknown-key error below

        for key in data:
            if (
                key not in PlotSettings._GLOBAL_FIELDS
                and get_archived_analysis_plugin(key) is not None
            ):
                raise ValueError(
                    format_archived_analysis_message(key, context="plot_settings section")
                )

        global_data: dict[str, Any] = {}
        per_analysis: dict[str, BasePlotSettings] = {}

        for key, value in data.items():
            if key in PlotSettings._GLOBAL_FIELDS:
                global_data[key] = value
            elif key in _plot_models:
                settings_class = _plot_models[key]
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
                raise ValueError(
                    f"Unknown plot settings key '{key}'. "
                    f"Expected a global field ({sorted(PlotSettings._GLOBAL_FIELDS)}) "
                    f"or a discovered analysis type with PlotSettingsModel."
                )

        # ── Resolve theme from style preset + user overrides ──
        style = _normalize_plot_style(global_data.get("style", "compact"))
        global_data["style"] = style
        theme_overrides = global_data.pop("theme", None)

        _THEME_PRESETS = {
            "compact": PlotTheme.compact,
            "large_elements": PlotTheme.large_elements,
            "low_ink": PlotTheme.low_ink,
        }
        preset_factory = _THEME_PRESETS[style]

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
        else:
            raise TypeError(
                f"Invalid 'theme' value: expected None, dict, or PlotTheme, "
                f"got {type(theme_overrides).__name__}"
            )

        super().__init__(**global_data, **per_analysis)

    def __getattr__(self, name: str) -> Any:
        """Fall back to default-constructed settings for discovered types.

        This ensures ``self.settings.rmsf.show_error`` works even when
        the user omitted the ``rmsf:`` block from their YAML.

        Parameters
        ----------
        name : str
            Attribute name.

        Returns
        -------
        BasePlotSettings
            Default-constructed settings if *name* is a discovered type.

        Raises
        ------
        AttributeError
            If *name* is not a discovered plot settings type.
        """
        from polyzymd.analyses.discovery import get_analysis

        try:
            analysis_cls = get_analysis(name)
        except KeyError as exc:
            raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'") from exc

        if (
            hasattr(analysis_cls, "PlotSettingsModel")
            and analysis_cls.PlotSettingsModel is not None
        ):
            return analysis_cls.PlotSettingsModel()

        raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'")

    @field_validator("output_dir", mode="before")
    @classmethod
    def resolve_output_dir(cls, v: str | Path) -> Path:
        """Convert string paths to Path objects."""
        return Path(v)


# ============================================================================
# MDAnalysis Backend Policy Configuration
# ============================================================================


class MDABackendPolicyConfig(BaseModel):
    """Configuration for optional MDAnalysis internal backend execution.

    The default configuration forwards no backend-related keyword arguments to
    MDAnalysis. This keeps PolyzyMD replicate-level scheduling as the default
    execution model and avoids nested oversubscription unless users explicitly
    opt in.
    """

    model_config = {"extra": "forbid"}

    backend: str | None = None
    n_workers: int | None = Field(default=None, ge=1)
    n_parts: int | None = Field(default=None, ge=1)

    @field_validator("backend")
    @classmethod
    def validate_backend(cls, value: str | None) -> str | None:
        """Reject empty backend names while preserving default serial behavior."""
        if value is None:
            return None
        backend = value.strip()
        if not backend:
            raise ValueError("backend must not be empty")
        return backend

    @field_validator("n_workers", "n_parts", mode="before")
    @classmethod
    def reject_bool_counts(cls, value: Any) -> Any:
        """Reject booleans before integer coercion handles worker counts."""
        if isinstance(value, bool):
            raise ValueError("worker and part counts must be positive integers")
        return value

    @model_validator(mode="after")
    def validate_backend_required(self) -> "MDABackendPolicyConfig":
        """Require an explicit backend before forwarding worker controls."""
        if self.backend is None and (self.n_workers is not None or self.n_parts is not None):
            raise ValueError("n_workers and n_parts require an explicit backend")
        return self

    def to_policy(self) -> Any:
        """Convert this config to an MDAnalysis backend policy.

        Returns
        -------
        MDABackendPolicy
            Import-light policy object consumed by ``MDAAnalysisJob``.
        """
        from polyzymd.analyses.mda import MDABackendPolicy

        return MDABackendPolicy(
            backend=self.backend,
            n_workers=self.n_workers,
            n_parts=self.n_parts,
        )


# ============================================================================
# Main Comparison Configuration
# ============================================================================


class ComparisonConfig(BaseModel):
    """Schema for comparison.yaml configuration files.

    A comparison config defines multiple simulation conditions to compare,
    along with unified plugin settings and plot customization.

    The schema follows a two-section pattern:
    - plugins: unified analysis settings keyed by plugin name
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
    plugins : PluginSettingsContainer
        Unified plugin parameters for compute, compare, and plot-aware metadata.
    mda_backend_policy : MDABackendPolicyConfig
        Optional MDAnalysis internal backend policy. Defaults to no backend
        keyword arguments.
    plot_settings : PlotSettings
        Plot customization (HOW to visualize)

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> print(config.name)
    "Polymer Stabilization Study"
    >>> for cond in config.conditions:
    ...     print(f"{cond.label}: {cond.config}")
    >>> print("Enabled analyses:", config.plugins.get_enabled_plugins())
    >>> rmsf_settings = config.plugins.get("rmsf")
    >>> if rmsf_settings:
    ...     print(f"RMSF selection: {rmsf_settings.selection}")
    """

    name: str
    description: str | None = None
    control: str | None = None
    conditions: list[ConditionConfig]
    defaults: AnalysisDefaults = Field(default_factory=AnalysisDefaults)
    plugins: PluginSettingsContainer = Field(default_factory=PluginSettingsContainer)
    mda_backend_policy: MDABackendPolicyConfig = Field(default_factory=MDABackendPolicyConfig)
    plot_settings: PlotSettings = Field(default_factory=PlotSettings)
    source_path: Path | None = Field(default=None, exclude=True)

    @field_validator("plugins", mode="before")
    @classmethod
    def parse_plugin_settings(cls, v: Any) -> PluginSettingsContainer:
        """Parse unified plugin settings from dict or container."""
        if v is None:
            return PluginSettingsContainer()
        if isinstance(v, dict):
            return PluginSettingsContainer(**v)
        return v

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

        if data is None:
            data = {}

        if "plugins" not in data and "analysis_settings" in data:
            logger.warning(
                "comparison.yaml uses legacy 'analysis_settings:'; treating it as 'plugins:' "
                "for backward compatibility."
            )
            data["plugins"] = data.pop("analysis_settings")

        if "analysis_settings" in data:
            logger.warning(
                "comparison.yaml contains both 'plugins' and legacy 'analysis_settings'; "
                "ignoring 'analysis_settings'."
            )
            data.pop("analysis_settings")

        allowed_keys = {
            "name",
            "description",
            "control",
            "conditions",
            "defaults",
            "plugins",
            "mda_backend_policy",
            "plot_settings",
        }
        unknown_keys = sorted(set(data.keys()) - allowed_keys)
        if unknown_keys:
            raise ValueError(
                f"comparison.yaml contains unknown top-level key(s): {unknown_keys}. "
                f"Valid keys: {', '.join(sorted(allowed_keys))}"
            )

        if "plugins" not in data:
            data["plugins"] = {}

        # Resolve relative paths relative to the config file location
        config_dir = path.parent.resolve()
        if "conditions" in data:
            for cond in data["conditions"]:
                if "config" in cond:
                    cond_path = Path(cond["config"])
                    if not cond_path.is_absolute():
                        cond["config"] = str(config_dir / cond_path)

        # Resolve plugin-declared settings path fields relative to comparison.yaml
        plugins = data.get("plugins")
        if isinstance(plugins, dict):
            from polyzymd.analyses.discovery import get_analysis

            yaml_dir = path.parent
            for plugin_name, plugin_settings in plugins.items():
                if not isinstance(plugin_settings, dict):
                    continue
                try:
                    analysis_cls = get_analysis(str(plugin_name).lower())
                except KeyError:
                    continue

                for field_name in analysis_cls.settings_path_fields:
                    raw_value = plugin_settings.get(field_name)
                    if not isinstance(raw_value, str):
                        continue
                    candidate = Path(raw_value)
                    if candidate.is_absolute():
                        continue
                    plugin_settings[field_name] = str(yaml_dir / candidate)

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

    @model_validator(mode="after")
    def _validate_conditions(self) -> "ComparisonConfig":
        """Validate conditions at construction time.

        Checks:
        - No duplicate labels
        - Control label exists in conditions (if specified)

        Note: The minimum-2-conditions check is deferred to
        :meth:`validate_config` because a partially-built config
        (e.g. during template generation) may have fewer.
        """
        labels = [c.label for c in self.conditions]
        if len(labels) != len(set(labels)):
            dupes = [lbl for lbl in labels if labels.count(lbl) > 1]
            raise ValueError(f"Duplicate condition labels: {set(dupes)}")
        if self.control and self.control not in labels:
            raise ValueError(f"Control '{self.control}' not in conditions: {labels}")
        return self

    def validate_config(self) -> list[str]:
        """Validate the comparison configuration (runtime checks).

        This method performs checks that cannot be done at construction time
        (e.g. filesystem existence checks).

        Returns
        -------
        list[str]
            List of error messages (empty if valid)
        """
        errors = []

        # Check minimum conditions
        if len(self.conditions) < 1:
            errors.append("Need at least 1 condition")

        # Check config files exist
        for cond in self.conditions:
            if not cond.config.exists():
                errors.append(f"Config not found for '{cond.label}': {cond.config}")

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
        data = self.plugins.to_analysis_yaml_dict(
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
    rendered = render_package_template(
        "polyzymd.templates",
        "comparison_template.yaml.jinja",
        {"name": name, "eq_time": eq_time},
    )
    return prepend_file_header(rendered, comment_prefix="#")
