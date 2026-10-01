"""Plot settings and the default equilibration window of ``polyzymd analyze``.

:class:`PlotSettings` holds the figure format, resolution, style preset,
palette, resolved :class:`PlotTheme` and :class:`SemanticColorSettings` that
the functions in :mod:`polyzymd.analyses.figures` and
:mod:`polyzymd.analyses.shared.plotting` read. :class:`AnalysisDefaults`
gives the equilibration window that ``polyzymd analyze`` discards when
``--eq`` is not given.
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

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
    """Package defaults of ``polyzymd analyze``.

    ``equilibration_time`` is the window discarded from the start of every
    replicate when ``--eq`` (or ``equilibration=``) is not given.
    ``fdr_alpha``, ``ttest_method`` and ``posthoc_method`` record the
    significance level, t test and post hoc method defaults.
    """

    equilibration_time: str = "10ns"
    fdr_alpha: float = Field(default=0.05, gt=0.0, le=1.0)
    ttest_method: str = Field(default="student", pattern="^(welch|student)$")
    posthoc_method: str = Field(default="ttest_bh", pattern="^(ttest_bh|tukey_hsd)$")


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
    """Condition colours and plotting order, from families, values, roles or explicit colours.

    Defaults define the canonical v1.3 plot behavior: semantic color mapping is
    disabled until users opt in with ``enabled: true``.
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
    """Font sizes, alphas, line widths, marker sizes and spines of every figure.

    Canonical presets are available via class methods:

    - ``PlotTheme.compact()`` — default; print-ready sizes and weights.
    - ``PlotTheme.large_elements()`` — ~1.3x larger fonts/dots/lines for slides.
    - ``PlotTheme.low_ink()`` — no dots, no bar edges, thinner lines.

    Override individual values through :class:`PlotSettings`::

        PlotSettings(style="compact", theme={"title_fontsize": 16, "dot_size": 24})

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
    """Figure format, resolution, style preset, palette, theme and condition colours.

    The functions in :mod:`polyzymd.analyses.figures` take one as
    ``plot_settings``; ``PlotSettings()`` is the default. Unknown keys raise a
    validation error.

    Attributes
    ----------
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
        any overrides passed as ``theme``.

    Examples
    --------
    >>> from polyzymd.config.analysis_settings import PlotSettings
    >>> settings = PlotSettings(format="pdf", style="large_elements", theme={"dot_size": 40})
    >>> ts.plot(plot_settings=settings)  # doctest: +SKIP
    """

    model_config = {"extra": "forbid"}

    format: str = Field(default="png", pattern="^(png|pdf|svg)$")
    dpi: int = Field(default=300, ge=50, le=600)
    style: str = "compact"
    color_palette: str = "tab10"
    theme: PlotTheme = Field(default_factory=PlotTheme)
    semantic_colors: SemanticColorSettings = Field(default_factory=SemanticColorSettings)

    def __init__(self, **data: Any):
        """Resolve ``theme`` from the ``style`` preset and any ``theme`` overrides.

        ``style`` selects the compact, large_elements or low_ink preset, and
        the keys of a ``theme`` mapping replace that preset's values, so
        ``style="large_elements", theme={"dot_size": 40}`` keeps the preset
        and changes only the dot size. A :class:`PlotTheme` passed as
        ``theme`` is used as it is.
        """
        style = _normalize_plot_style(data.get("style", "compact"))
        data["style"] = style
        theme_overrides = data.pop("theme", None)
        preset_factory = {
            "compact": PlotTheme.compact,
            "large_elements": PlotTheme.large_elements,
            "low_ink": PlotTheme.low_ink,
        }[style]
        if theme_overrides is None or (isinstance(theme_overrides, dict) and not theme_overrides):
            data["theme"] = preset_factory()
        elif isinstance(theme_overrides, dict):
            data["theme"] = PlotTheme(**{**preset_factory().model_dump(), **theme_overrides})
        elif isinstance(theme_overrides, PlotTheme):
            data["theme"] = theme_overrides
        else:
            raise TypeError(
                f"Invalid 'theme' value: expected None, dict, or PlotTheme, "
                f"got {type(theme_overrides).__name__}"
            )
        super().__init__(**data)
