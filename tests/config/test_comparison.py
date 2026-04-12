"""Tests for comparison config validation."""

from pathlib import Path

import pytest
import yaml

from polyzymd.config.comparison import ComparisonConfig, ConditionConfig, PlotSettings


class TestUnknownTopLevelKeys:
    """C8-H1: Unknown comparison.yaml keys should raise ValueError."""

    def test_typo_key_raises(self, tmp_path: Path) -> None:
        """Typo like 'contol' instead of 'control' should error."""
        yaml_path = tmp_path / "comparison.yaml"
        yaml_path.write_text(
            yaml.dump(
                {
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                        {"label": "B", "config": "/fake/b.yaml", "replicates": [1]},
                    ],
                    "contol": "A",
                }
            )
        )
        with pytest.raises(ValueError, match="unknown top-level key.*contol"):
            ComparisonConfig.from_yaml(yaml_path)

    def test_multiple_unknown_keys_raises(self, tmp_path: Path) -> None:
        """Multiple unknown keys listed in error."""
        yaml_path = tmp_path / "comparison.yaml"
        yaml_path.write_text(
            yaml.dump(
                {
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                    ],
                    "foo": 1,
                    "bar": 2,
                }
            )
        )
        with pytest.raises(ValueError, match="unknown top-level key"):
            ComparisonConfig.from_yaml(yaml_path)

    def test_valid_keys_accepted(self, tmp_path: Path) -> None:
        """All valid keys should be accepted without error."""
        yaml_path = tmp_path / "comparison.yaml"
        yaml_path.write_text(
            yaml.dump(
                {
                    "name": "test",
                    "description": "test comparison",
                    "control": "A",
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                        {"label": "B", "config": "/fake/b.yaml", "replicates": [1]},
                    ],
                    "defaults": {"equilibration_time": "10ns"},
                    "plugins": {},
                    "plot_settings": {},
                }
            )
        )
        config = ComparisonConfig.from_yaml(yaml_path)
        assert config.name == "test"

    def test_legacy_analysis_settings_still_migrates(self, tmp_path: Path) -> None:
        """Legacy analysis_settings key should warn but migrate to plugins."""
        yaml_path = tmp_path / "comparison.yaml"
        yaml_path.write_text(
            yaml.dump(
                {
                    "name": "legacy-migration",
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                    ],
                    "analysis_settings": {"rmsf": {"selection": "protein"}},
                }
            )
        )

        config = ComparisonConfig.from_yaml(yaml_path)
        assert config is not None


class TestConditionReplicateValidation:
    """C8-H2: ConditionConfig should reject bad replicate lists."""

    def test_scalar_coerced_to_list(self):
        """Single int replicate should become [int]."""
        cond = ConditionConfig(label="A", config=Path("/fake/a.yaml"), replicates=3)
        assert cond.replicates == [3]

    def test_empty_replicates_raises(self):
        """Empty replicate list should raise."""
        with pytest.raises(ValueError, match="must not be empty"):
            ConditionConfig(label="A", config=Path("/fake/a.yaml"), replicates=[])

    def test_duplicate_replicates_raises(self):
        """Duplicate replicate numbers should raise."""
        with pytest.raises(ValueError, match="duplicate"):
            ConditionConfig(label="A", config=Path("/fake/a.yaml"), replicates=[1, 2, 1])

    def test_zero_replicate_raises(self):
        """Zero replicate number should raise."""
        with pytest.raises(ValueError, match="positive"):
            ConditionConfig(label="A", config=Path("/fake/a.yaml"), replicates=[0, 1])

    def test_negative_replicate_raises(self):
        """Negative replicate number should raise."""
        with pytest.raises(ValueError, match="positive"):
            ConditionConfig(label="A", config=Path("/fake/a.yaml"), replicates=[-1])

    def test_valid_replicates_accepted(self):
        """Normal replicate list should work."""
        cond = ConditionConfig(label="A", config=Path("/fake/a.yaml"), replicates=[1, 2, 3])
        assert cond.replicates == [1, 2, 3]


class TestPlotThemeValidation:
    """C8-H3: PlotSettings should reject invalid theme types."""

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
        """None theme should use publication preset."""
        ps = PlotSettings(theme=None)
        assert ps.theme is not None
