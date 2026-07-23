"""Tests for config loader path handling and YAML dumping."""

from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError


def _write_config(path: Path, conjugation: dict) -> None:
    """Write a minimal simulation config containing conjugation data."""
    data = {
        "name": "csv-attachments",
        "engine": "openmm",
        "enzyme": {"name": "TestEnzyme", "pdb_path": "enzyme.pdb"},
        "thermodynamics": {"temperature": 300.0},
        "simulation_phases": {
            "equilibration_stages": [
                {
                    "name": "eq1",
                    "duration": 0.1,
                    "temperature": 300.0,
                    "ensemble": "NVT",
                }
            ],
            "production": {
                "ensemble": "NPT",
                "duration": 1.0,
                "samples": 10,
                "checkpoint_interval": 60.0,
            },
        },
        "conjugation": conjugation,
    }
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


class TestAttachmentCsvLoading:
    """Attachment CSV files are materialized before schema validation."""

    def test_loads_bom_rows_in_order_and_resolves_csv_relative_moiety_path(
        self, tmp_path: Path
    ) -> None:
        """CSV rows preserve order and use the CSV directory for moiety paths."""
        from polyzymd.config.loader import load_config, save_config

        csv_dir = tmp_path / "tables"
        csv_dir.mkdir()
        csv_path = csv_dir / "attachments.csv"
        csv_path.write_text(
            "\ufeffname,moiety.name,moiety.input_path,site.residue_number,site.atom_name\n"
            "first,modifier-one,structures/one.pdb,10,NZ\n"
            "second,modifier-two,,11,SG\n",
            encoding="utf-8",
        )
        yaml_path = tmp_path / "config.yaml"
        _write_config(yaml_path, {"attachments_file": "tables/attachments.csv"})

        config = load_config(yaml_path)

        assert [item.name for item in config.conjugation.attachments] == ["first", "second"]
        assert config.conjugation.attachments[0].moiety.input_path == (
            csv_dir / "structures/one.pdb"
        )
        assert config.conjugation.attachments[1].moiety.input_path is None

        saved_path = tmp_path / "saved.yaml"
        save_config(config, saved_path)
        saved = yaml.safe_load(saved_path.read_text(encoding="utf-8"))
        assert "attachments_file" not in saved["conjugation"]
        assert [item["name"] for item in saved["conjugation"]["attachments"]] == [
            "first",
            "second",
        ]

    @pytest.mark.parametrize(
        ("csv_text", "message"),
        [
            ("name,name,moiety.name\na,b,c\n", "Duplicate attachments CSV header"),
            ("name,unknown\na,b\n", "Unsupported attachments CSV header"),
            (
                "name,mechanism.leaving_atoms.site\na,OW\n",
                "does not support list-valued header",
            ),
            ("name,moiety.name\na,b,extra\n", "row 2 has extra cells"),
        ],
    )
    def test_rejects_malformed_csv_shapes(
        self, tmp_path: Path, csv_text: str, message: str
    ) -> None:
        """Ambiguous or non-scalar CSV shapes fail with targeted messages."""
        from polyzymd.config.loader import load_config

        (tmp_path / "attachments.csv").write_text(csv_text, encoding="utf-8")
        yaml_path = tmp_path / "config.yaml"
        _write_config(yaml_path, {"attachments_file": "attachments.csv"})

        with pytest.raises(ValidationError, match=message):
            load_config(yaml_path)

    def test_pydantic_reports_inflated_attachment_field_location(self, tmp_path: Path) -> None:
        """Semantic cell errors retain the ordinary attachment field location."""
        from polyzymd.config.loader import load_config

        (tmp_path / "attachments.csv").write_text(
            "name,moiety.name,site.residue_number\ninvalid,modifier,not-an-integer\n",
            encoding="utf-8",
        )
        yaml_path = tmp_path / "config.yaml"
        _write_config(yaml_path, {"attachments_file": "attachments.csv"})

        with pytest.raises(ValidationError) as exc_info:
            load_config(yaml_path)

        assert any(
            error["loc"] == ("conjugation", "attachments", 0, "site", "residue_number")
            for error in exc_info.value.errors()
        )


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

    def test_schema_resolves_default_to_bundled_paths(self):
        """ReactionConfig must resolve 'default' to real .rxn file paths."""
        from polyzymd.config.schema import ReactionConfig

        cfg = ReactionConfig(
            initiation="default",
            polymerization="default",
            termination="default",
        )
        assert cfg.initiation.suffix == ".rxn"
        assert cfg.polymerization.suffix == ".rxn"
        assert cfg.termination.suffix == ".rxn"
        assert cfg.initiation.exists()
        assert cfg.polymerization.exists()
        assert cfg.termination.exists()

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
            path = getattr(cfg, field)
            assert path.suffix == ".rxn"
            assert "atrp" in path.name
            assert path.exists()
