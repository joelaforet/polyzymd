"""Comparison configuration schema in ``polyzymd.config``."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, ClassVar

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import BasePlotSettings
from polyzymd.analyses.shared.defaults import AnalysisDefaults
from polyzymd.core.branding import prepend_file_header

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
            if value is None:
                continue
            key_lower = key.lower()
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


class PlotTheme(BaseModel):
    """Centralized visual defaults for all comparison plots.

    Replaces ~219 hardcoded style values (font sizes, alphas, line widths,
    marker sizes, spine visibility, etc.) across all plugin ``plot()`` methods
    with a single configurable Pydantic model.

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
        Plot style preset: "publication", "presentation", or "minimal"
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
          style: "publication"

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
    }

    output_dir: Path = Field(default=Path("figures/"))
    format: str = Field(default="png", pattern="^(png|pdf|svg)$")
    dpi: int = Field(default=300, ge=50, le=600)
    style: str = Field(default="publication", pattern="^(publication|presentation|minimal)$")
    color_palette: str = "tab10"
    theme: PlotTheme = Field(default_factory=PlotTheme)

    def __init__(self, **data: Any):
        """Initialize with global fields and plugin-discovered per-analysis settings.

        Theme resolution: the ``style`` field selects a preset (publication,
        presentation, or minimal) and then any user-supplied ``theme:``
        overrides are merged on top.  This allows ``style: presentation``
        with ``theme: {dot_size: 40}`` to use the presentation preset but
        override just the dot size.

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
    return prepend_file_header(f"""\
# PolyzyMD Comparison Configuration
# Compare analyses across simulation conditions (e.g., polymer, temperature).
# Docs: https://polyzymd.readthedocs.io/en/latest/

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
# Plugins (unified analysis settings)
# ============================================================================
# Define which analyses to run. Presence of a section enables that analysis.

plugins:
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

  # SASA (Solvent Accessible Surface Area) Analysis
  #
  # sasa:
  #   runs:
  #     - label: "protein_total"
  #       target_selection: "protein"
  #       context_selection: "protein"  # optional, defaults to target_selection
  #       stride: 1
  #   probe_radius_nm: 0.14
  #   n_sphere_points: 960
  #   chunk_size: 100

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
  # "Below {{threshold}}Å" / "Above {{threshold}}Å".
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

  # Hydrogen-bond Analysis
  #
  # hydrogen_bonds:
  #   distance_cutoff: 3.0
  #   angle_cutoff: 150
  #   groups:
  #     protein_all: "protein"
  #     substrate: "resname LIG"     # adjust to your ligand
  #     polymer_all: "chainid C"
  #   summaries:
  #     protein_polymer:
  #       between: [protein_all, polymer_all]
  #     protein_substrate:
  #       between: [protein_all, substrate]
  #     protein_internal:
  #       within: protein_all
  #     polymer_internal:
  #       within: polymer_all
  #   composition:
  #     partitions:
  #       protein: "protein"
  #       substrate: "resname LIG"
  #       polymer: "chainid C"

  # EXPERIMENTAL: Exposure Dynamics Analysis (chaperone-like polymer activity)
  # Requires contacts analysis to be run first for each condition.
  # Definitions and interpretation may change after the presentation release.
  # Run: polyzymd compare run exposure
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
  # Run: polyzymd compare run binding_free_energy
  #
  # binding_free_energy:
  #   units: kT                    # energy units: kT (default), kcal/mol, or kJ/mol
  #   surface_exposure_threshold: 0.2  # minimum relative SASA to be considered surface-exposed
  #   # protein_partitions: null   # optional: restrict to user-defined AA partitions

  # EXPERIMENTAL: Polymer Affinity Score (composite selectivity metric)
  # Requires contacts analysis with compute_binding_preference: true.
  # Computes S = Σ N_pg × ΔG_sel_pg for each polymer composition.
  # Definitions and interpretation may change after the presentation release.
  # Run: polyzymd compare run polymer_affinity
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
# Plot Settings (HOW to visualize - figure customization)
# ============================================================================
# Controls plot generation for all analyses. Per-analysis sections override
# defaults. Run `polyzymd compare plot-all` to generate all configured plots.

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

  # catalytic_triad:
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
""")
