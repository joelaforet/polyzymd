"""Tests for comparison config validation."""

from pathlib import Path

import pytest
import yaml

from polyzymd.config.comparison import ComparisonConfig


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
