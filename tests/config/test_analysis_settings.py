"""Tests for the plot settings of polyzymd.config.analysis_settings."""

import pytest

from polyzymd.config.analysis_settings import (
    PlotSettings,
    PlotTheme,
    SemanticColorSettings,
    SemanticFamilyColorConfig,
)


class TestPlotThemeValidation:
    """PlotSettings rejects a theme that is not a mapping or a PlotTheme."""

    def test_string_theme_raises(self):
        with pytest.raises(TypeError, match="Invalid 'theme' value"):
            PlotSettings(theme="bad")

    def test_int_theme_raises(self):
        with pytest.raises(TypeError, match="Invalid 'theme' value"):
            PlotSettings(theme=123)

    def test_list_theme_raises(self):
        with pytest.raises(TypeError, match="Invalid 'theme' value"):
            PlotSettings(theme=["bad"])

    def test_dict_theme_accepted(self):
        """Dict overrides should still work."""
        ps = PlotSettings(theme={"dot_size": 42})
        assert ps.theme.dot_size == 42

    def test_none_theme_uses_default(self):
        """None theme should use compact preset."""
        ps = PlotSettings(theme=None)
        assert ps.theme is not None


class TestPlotStylePresets:
    """Plot style presets should use canonical names only."""

    def test_default_style_is_compact(self) -> None:
        """Omitted style should use the compact preset values."""
        settings = PlotSettings()

        assert settings.style == "compact"
        assert settings.theme == PlotTheme.compact()

    @pytest.mark.parametrize(
        ("style", "expected_theme"),
        [
            ("compact", PlotTheme.compact()),
            ("large_elements", PlotTheme.large_elements()),
            ("low_ink", PlotTheme.low_ink()),
        ],
    )
    def test_canonical_styles_are_accepted(
        self,
        style: str,
        expected_theme: PlotTheme,
    ) -> None:
        """Canonical style names should select the expected theme."""
        settings = PlotSettings(style=style)

        assert settings.style == style
        assert settings.theme == expected_theme

    @pytest.mark.parametrize(
        "style",
        [
            "publication",
            "presentation",
            "minimal",
        ],
    )
    def test_old_aliases_raise_value_error(self, style: str) -> None:
        """Old style aliases should be rejected as invalid strings."""
        with pytest.raises(ValueError, match="compact.*large_elements.*low_ink"):
            PlotSettings(style=style)

    def test_theme_overrides_merge_after_canonical_style(self) -> None:
        """Theme overrides should layer on top of canonical presets."""
        settings = PlotSettings(style="large_elements", theme={"dot_size": 40})

        assert settings.style == "large_elements"
        assert settings.theme.dot_size == 40
        assert settings.theme.title_fontsize == PlotTheme.large_elements().title_fontsize

    def test_invalid_style_string_raises_value_error(self) -> None:
        """Invalid style strings should list allowed canonical values."""
        with pytest.raises(ValueError, match="compact.*large_elements.*low_ink"):
            PlotSettings(style="poster")

    @pytest.mark.parametrize("style", [123, None, ["compact"]])
    def test_non_string_style_raises_type_error(self, style: object) -> None:
        """Non-string style values should produce a clear type error."""
        with pytest.raises(TypeError, match="plot_settings.style must be a string"):
            PlotSettings(style=style)


class TestSemanticColorSettings:
    """Semantic color settings should parse canonical defaults."""

    def test_defaults_disable_semantic_mapping(self) -> None:
        """Semantic colors should be opt-in by default."""
        settings = PlotSettings()

        assert settings.semantic_colors.enabled is False
        assert settings.semantic_colors.order == []
        assert settings.semantic_colors.conditions == {}
        assert settings.semantic_colors.families == {}
        assert settings.semantic_colors.default_color is None

    def test_plot_settings_accepts_semantic_colors(self) -> None:
        """PlotSettings should parse the global semantic_colors block."""
        settings = PlotSettings(
            semantic_colors={
                "enabled": True,
                "order": ["Control", "Condition B"],
                "manual_colors": {"Control": "#222222"},
                "conditions": {
                    "Control": {"role": "Control", "order": 0},
                    "Condition B": {"family": "composition", "value": 0.5},
                },
                "families": {
                    "composition": {
                        "scale": "linear",
                        "colormap": "viridis",
                        "vmin": 0.0,
                        "vmax": 1.0,
                        "colormap_range": [0.2, 0.8],
                    }
                },
            }
        )

        semantic = settings.semantic_colors
        assert semantic.enabled is True
        assert semantic.order == ["Control", "Condition B"]
        assert semantic.conditions["Control"].role == "control"
        assert semantic.conditions["Condition B"].family == "composition"
        assert semantic.families["composition"].colormap_range == (0.2, 0.8)

    def test_semantic_colors_is_global_plot_field(self) -> None:
        """semantic_colors takes a SemanticColorSettings instance as well as a mapping."""
        settings = PlotSettings(semantic_colors=SemanticColorSettings(enabled=True))

        assert settings.semantic_colors.enabled is True

    def test_duplicate_explicit_order_raises(self) -> None:
        """Explicit plot order labels should be unique."""
        with pytest.raises(ValueError, match="order labels must be unique"):
            SemanticColorSettings(order=["A", "A"])

    def test_invalid_family_scale_raises(self) -> None:
        """Family color scale should be linear or ordinal."""
        with pytest.raises(ValueError):
            SemanticFamilyColorConfig(scale="bad")

    def test_invalid_colormap_range_raises(self) -> None:
        """Family colormap ranges should be ordered fractions."""
        with pytest.raises(ValueError, match="colormap_range"):
            SemanticFamilyColorConfig(colormap_range=(0.8, 0.2))


class TestUnknownKeys:
    """PlotSettings takes no per-analysis block and no output_dir."""

    @pytest.mark.parametrize("key", ["rmsf", "output_dir", "contacts"])
    def test_unknown_key_is_rejected(self, key: str) -> None:
        with pytest.raises(ValueError, match=key):
            PlotSettings(**{key: {}})
