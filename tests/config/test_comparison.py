"""Tests for comparison config validation."""

from pathlib import Path
from typing import ClassVar

import pytest
import yaml
from pydantic import BaseModel

from polyzymd.config.comparison import (
    ComparisonConfig,
    ConditionConfig,
    MDABackendPolicyConfig,
    PlotSettings,
)


class _FakePathSettings(BaseModel):
    """Settings model with a plugin-declared path field."""

    enzyme_pdb_for_sasa: str | None = None
    units: str = "kT"


class _FakePathAnalysis:
    """Fake plugin used to test config path resolution."""

    name: ClassVar[str] = "fake_paths"
    Settings: ClassVar[type[_FakePathSettings]] = _FakePathSettings
    settings_path_fields: ClassVar[tuple[str, ...]] = ("enzyme_pdb_for_sasa",)


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
                    "mda_backend_policy": {},
                    "plot_settings": {},
                }
            )
        )
        config = ComparisonConfig.from_yaml(yaml_path)
        assert config.name == "test"


class TestMDABackendPolicyConfig:
    """MDAnalysis backend policy config should be opt-in and validated."""

    def test_default_policy_forwards_no_backend_kwargs(self) -> None:
        """Default config should keep per-replicate execution serial."""
        policy = MDABackendPolicyConfig().to_policy()

        assert policy.run_kwargs() == {}
        assert policy.is_default()

    def test_explicit_backend_accepts_workers_and_parts(self) -> None:
        """Explicit backend should convert to an MDA job policy."""
        policy = MDABackendPolicyConfig(
            backend="multiprocessing",
            n_workers=2,
            n_parts=4,
        ).to_policy()

        assert policy.run_kwargs() == {
            "backend": "multiprocessing",
            "n_workers": 2,
            "n_parts": 4,
        }

    @pytest.mark.parametrize("field_name", ["n_workers", "n_parts"])
    def test_worker_controls_require_backend(self, field_name: str) -> None:
        """Worker and part controls should not imply backend opt-in."""
        with pytest.raises(ValueError, match="require an explicit backend"):
            MDABackendPolicyConfig(**{field_name: 2})

    @pytest.mark.parametrize("field_name", ["n_workers", "n_parts"])
    def test_worker_controls_must_be_positive_integers(self, field_name: str) -> None:
        """Worker and part counts should reject invalid values."""
        with pytest.raises(ValueError):
            MDABackendPolicyConfig(backend="multiprocessing", **{field_name: 0})

        with pytest.raises(ValueError, match="positive integers"):
            MDABackendPolicyConfig(backend="multiprocessing", **{field_name: True})

    def test_backend_name_must_not_be_empty(self) -> None:
        """Empty backend strings should not be treated as backend opt-in."""
        with pytest.raises(ValueError, match="backend must not be empty"):
            MDABackendPolicyConfig(backend="  ")

    def test_from_yaml_accepts_top_level_policy(self, tmp_path: Path) -> None:
        """comparison.yaml should expose a top-level MDA backend policy."""
        yaml_path = tmp_path / "comparison.yaml"
        yaml_path.write_text(
            yaml.dump(
                {
                    "name": "backend-test",
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                    ],
                    "mda_backend_policy": {
                        "backend": "multiprocessing",
                        "n_workers": 2,
                    },
                }
            )
        )

        config = ComparisonConfig.from_yaml(yaml_path)

        assert config.mda_backend_policy.to_policy().run_kwargs() == {
            "backend": "multiprocessing",
            "n_workers": 2,
        }

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


class TestArchivedAnalysisDiagnostics:
    """Archived analyses should fail with recovery diagnostics."""

    @pytest.mark.parametrize(
        "analysis_name",
        [
            "exposure",
            "binding_free_energy",
            "bfe",
            "polymer_affinity",
            "pa",
            "polymer_bridging",
            "bridging",
        ],
    )
    def test_plugins_archived_analysis_reports_archive_location(
        self,
        analysis_name: str,
        tmp_path: Path,
    ) -> None:
        """Archived plugin settings should explain where the archived code lives."""
        yaml_path = tmp_path / "comparison.yaml"
        yaml_path.write_text(
            yaml.dump(
                {
                    "name": "archive-test",
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                    ],
                    "plugins": {analysis_name: {}},
                }
            )
        )

        with pytest.raises(ValueError) as excinfo:
            ComparisonConfig.from_yaml(yaml_path)

        message = str(excinfo.value)
        assert "plugins section" in message
        assert "archive_experimental_analysis" in message
        assert "feature/mda-analysis-migration" in message

    @pytest.mark.parametrize(
        "analysis_name",
        [
            "exposure",
            "binding_free_energy",
            "bfe",
            "polymer_affinity",
            "pa",
            "polymer_bridging",
            "bridging",
        ],
    )
    def test_plot_settings_archived_analysis_reports_archive_location(
        self,
        analysis_name: str,
        tmp_path: Path,
    ) -> None:
        """Archived plot settings should explain where the archived code lives."""
        yaml_path = tmp_path / "comparison.yaml"
        yaml_path.write_text(
            yaml.dump(
                {
                    "name": "archive-test",
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                    ],
                    "plugins": {},
                    "plot_settings": {analysis_name: {}},
                }
            )
        )

        with pytest.raises(ValueError) as excinfo:
            ComparisonConfig.from_yaml(yaml_path)

        message = str(excinfo.value)
        assert "plot_settings section" in message
        assert "archive_experimental_analysis" in message
        assert "feature/mda-analysis-migration" in message


class TestPluginSettingsPathResolution:
    """Plugin path-like settings should resolve relative to comparison.yaml."""

    def test_relative_plugin_path_resolves_from_yaml_parent(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        """Relative plugin path fields should be rebased to yaml directory."""
        monkeypatch.setattr(
            "polyzymd.analyses.discovery.get_analysis",
            lambda name: _FakePathAnalysis,
        )
        yaml_path = tmp_path / "comparison.yaml"
        yaml_path.write_text(
            yaml.dump(
                {
                    "name": "path-test",
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                    ],
                    "plugins": {
                        "fake_paths": {
                            "enzyme_pdb_for_sasa": "structures/enzyme.pdb",
                        }
                    },
                }
            )
        )

        config = ComparisonConfig.from_yaml(yaml_path)
        settings = config.plugins.get("fake_paths")
        assert settings is not None
        assert settings.enzyme_pdb_for_sasa == str(tmp_path / "structures" / "enzyme.pdb")

    def test_absolute_plugin_path_is_unchanged(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        """Absolute plugin path fields should remain unchanged."""
        monkeypatch.setattr(
            "polyzymd.analyses.discovery.get_analysis",
            lambda name: _FakePathAnalysis,
        )
        yaml_path = tmp_path / "comparison.yaml"
        absolute_path = tmp_path / "already_absolute" / "enzyme.pdb"
        yaml_path.write_text(
            yaml.dump(
                {
                    "name": "path-test",
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                    ],
                    "plugins": {
                        "fake_paths": {
                            "enzyme_pdb_for_sasa": str(absolute_path),
                        }
                    },
                }
            )
        )

        config = ComparisonConfig.from_yaml(yaml_path)
        settings = config.plugins.get("fake_paths")
        assert settings is not None
        assert settings.enzyme_pdb_for_sasa == str(absolute_path)

    def test_missing_declared_path_field_is_ignored(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        """Declared path fields missing from YAML should be ignored."""
        monkeypatch.setattr(
            "polyzymd.analyses.discovery.get_analysis",
            lambda name: _FakePathAnalysis,
        )
        yaml_path = tmp_path / "comparison.yaml"
        yaml_path.write_text(
            yaml.dump(
                {
                    "name": "path-test",
                    "conditions": [
                        {"label": "A", "config": "/fake/a.yaml", "replicates": [1]},
                    ],
                    "plugins": {
                        "fake_paths": {
                            "units": "kT",
                        }
                    },
                }
            )
        )

        config = ComparisonConfig.from_yaml(yaml_path)
        settings = config.plugins.get("fake_paths")
        assert settings is not None
        assert settings.enzyme_pdb_for_sasa is None
