"""Tests for config loader path handling and YAML dumping."""

from pathlib import Path

import pytest


class TestSaveConfigYamlDumper:
    """save_config uses a local Dumper subclass, not the global yaml.Dumper."""

    def test_global_dumper_unmodified_after_save(self, tmp_path):
        """Calling save_config must not add representers to yaml.Dumper."""
        import yaml

        from polyzymd.config.loader import save_config
        from polyzymd.config.schema import SimulationConfig

        original_representer = yaml.Dumper.yaml_representers.get(str)

        config = SimulationConfig.model_construct(
            name="test",
            _fields_set=set(),
        )
        out = tmp_path / "test_config.yaml"
        save_config(config, out)

        after_representer = yaml.Dumper.yaml_representers.get(str)
        message = "save_config must not mutate yaml.Dumper.yaml_representers"
        assert after_representer is original_representer, message

    def test_multiline_strings_use_block_style(self, tmp_path):
        """Multiline strings in saved YAML use literal block (|) style."""
        from polyzymd.config.loader import save_config
        from polyzymd.config.schema import SimulationConfig

        config = SimulationConfig.model_construct(
            name="multi\nline\ntest",
            _fields_set=set(),
        )
        out = tmp_path / "test_config.yaml"
        save_config(config, out)

        content = out.read_text()
        assert "|-" in content or "|\n" in content


class TestReactionPathResolution:
    """_expand_paths and _convert_paths_to_relative handle reaction template keys."""

    def test_expand_paths_resolves_initiation(self):
        """Relative initiation path should be expanded to absolute."""
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"initiation": "templates/init.rxn"}}
        base = Path("/configs")
        result = _expand_paths(data, base)
        assert result["reactions"]["initiation"] == "/configs/templates/init.rxn"

    def test_expand_paths_resolves_polymerization(self):
        """Relative polymerization path should be expanded to absolute."""
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"polymerization": "templates/poly.rxn"}}
        base = Path("/configs")
        result = _expand_paths(data, base)
        assert result["reactions"]["polymerization"] == "/configs/templates/poly.rxn"

    def test_expand_paths_resolves_termination(self):
        """Relative termination path should be expanded to absolute."""
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"termination": "templates/term.rxn"}}
        base = Path("/configs")
        result = _expand_paths(data, base)
        assert result["reactions"]["termination"] == "/configs/templates/term.rxn"

    def test_expand_paths_preserves_absolute_reaction_path(self):
        """Absolute reaction paths should not be modified."""
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"initiation": "/absolute/init.rxn"}}
        base = Path("/configs")
        result = _expand_paths(data, base)
        assert result["reactions"]["initiation"] == "/absolute/init.rxn"

    def test_convert_paths_relativizes_reaction_paths(self):
        """Absolute reaction paths should be converted to relative for saving."""
        from polyzymd.config.loader import _convert_paths_to_relative

        data = {"reactions": {"initiation": "/configs/templates/init.rxn"}}
        base = Path("/configs")
        result = _convert_paths_to_relative(data, base)
        assert result["reactions"]["initiation"] == "templates/init.rxn"


class TestReactionDefaultSentinel:
    """Default reaction sentinel should pass through path transforms."""

    def test_expand_paths_preserves_default_initiation(self):
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"initiation": "default"}}
        result = _expand_paths(data, Path("/some/project"))
        assert result["reactions"]["initiation"] == "default"

    def test_expand_paths_preserves_default_polymerization(self):
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"polymerization": "default"}}
        result = _expand_paths(data, Path("/some/project"))
        assert result["reactions"]["polymerization"] == "default"

    def test_expand_paths_preserves_default_termination(self):
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"termination": "default"}}
        result = _expand_paths(data, Path("/some/project"))
        assert result["reactions"]["termination"] == "default"

    def test_expand_paths_preserves_default_case_insensitive(self):
        """'Default', 'DEFAULT', ' default ' should all be preserved."""
        from polyzymd.config.loader import _expand_paths

        for variant in ("Default", "DEFAULT", " default ", "  Default  "):
            data = {"reactions": {"initiation": variant}}
            result = _expand_paths(data, Path("/some/project"))
            message = f"Sentinel variant {variant!r} was not preserved"
            assert result["reactions"]["initiation"] == variant, message

    def test_expand_paths_still_expands_non_sentinel_reaction_paths(self):
        """Regular .rxn paths must still be expanded to absolute."""
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"initiation": "templates/init.rxn"}}
        result = _expand_paths(data, Path("/configs"))
        assert result["reactions"]["initiation"] == "/configs/templates/init.rxn"

    def test_relativize_preserves_default_initiation(self):
        from polyzymd.config.loader import _convert_paths_to_relative

        data = {"reactions": {"initiation": "default"}}
        result = _convert_paths_to_relative(data, Path("/some/project"))
        assert result["reactions"]["initiation"] == "default"

    def test_relativize_preserves_default_all_keys(self):
        from polyzymd.config.loader import _convert_paths_to_relative

        data = {
            "reactions": {
                "initiation": "default",
                "polymerization": "default",
                "termination": "default",
            }
        }
        result = _convert_paths_to_relative(data, Path("/some/project"))
        for key in ("initiation", "polymerization", "termination"):
            assert result["reactions"][key] == "default"

    def test_schema_preserves_default_markers(self):
        """ReactionConfig must preserve native default markers."""
        from polyzymd.config.schema import ReactionConfig

        cfg = ReactionConfig(
            initiation="default",
            polymerization="default",
            termination="default",
        )
        assert cfg.initiation == "default"
        assert cfg.polymerization == "default"
        assert cfg.termination == "default"

    def test_expand_then_validate_default_end_to_end(self):
        """Simulate the full load path: _expand_paths then model_validate."""
        from polyzymd.config.loader import _expand_paths
        from polyzymd.config.schema import ReactionConfig

        raw = {
            "initiation": "default",
            "polymerization": "default",
            "termination": "default",
        }
        expanded = _expand_paths({"reactions": raw}, Path("/any/base/dir"))
        cfg = ReactionConfig.model_validate(expanded["reactions"])
        for field in ("initiation", "polymerization", "termination"):
            assert getattr(cfg, field) == "default"

    def test_custom_rxn_validate_end_to_end_rejects(self):
        """Custom reaction paths expanded by the loader should still be rejected."""
        from polyzymd.config.loader import _expand_paths
        from polyzymd.config.schema import ReactionConfig

        raw = {
            "initiation": "custom_init.rxn",
            "polymerization": "custom_poly.rxn",
            "termination": "custom_term.rxn",
        }
        expanded = _expand_paths({"reactions": raw}, Path("/any/base/dir"))
        with pytest.raises(ValueError, match="provided_molecules"):
            ReactionConfig.model_validate(expanded["reactions"])
