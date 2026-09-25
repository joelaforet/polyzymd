"""The ``polyzymd study`` commands and ``polyzymd new-analysis`` inside a study."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from polyzymd.analyses import discovery
from polyzymd.cli.main import cli
from polyzymd.cli.study import STUDY_FOLDERS


@pytest.fixture(autouse=True)
def _forget_registered_analyses():
    """Leave discovery and ``sys.modules`` as each test found them."""
    discovery.clear_cache()
    yield
    discovery.clear_cache()
    for name in [module for module in sys.modules if module.startswith("polyzymd_study_")]:
        del sys.modules[name]


def _init(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str = "thermal") -> Path:
    """Run ``polyzymd study init`` in ``tmp_path`` and return the study root."""
    monkeypatch.chdir(tmp_path)
    result = CliRunner().invoke(cli, ["study", "init", "-n", name])
    assert result.exit_code == 0, result.output
    return tmp_path / name


class TestStudyInit:
    """``polyzymd study init`` lays out a study that git can keep whole."""

    def test_creates_the_layout(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        root = _init(tmp_path, monkeypatch)

        assert yaml.safe_load((root / "study.yaml").read_text()) == {
            "name": "thermal",
            "description": None,
            "analyses": "analyses",
        }
        for folder in STUDY_FOLDERS:
            assert (root / folder / ".gitkeep").is_file()
        assert "*.dcd" in (root / ".gitignore").read_text()
        assert "polyzymd new-analysis" in (root / "README.md").read_text()

    def test_refuses_an_existing_directory(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        (tmp_path / "thermal").mkdir()
        monkeypatch.chdir(tmp_path)

        result = CliRunner().invoke(cli, ["study", "init", "-n", "thermal"])

        assert result.exit_code == 1
        assert "already exists" in result.output


class TestNewAnalysisInAStudy:
    """Inside a study the scaffold writes into the study, not the source tree."""

    def test_writes_the_plugin_and_its_test_into_analyses(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        root = _init(tmp_path, monkeypatch)
        (root / "conditions" / "a").mkdir()
        monkeypatch.chdir(root / "conditions" / "a")

        result = CliRunner().invoke(cli, ["new-analysis", "lid_opening"])

        assert result.exit_code == 0, result.output
        plugin = root / "analyses" / "lid_opening.py"
        test = root / "analyses" / "test_lid_opening.py"
        assert plugin.is_file() and test.is_file()
        assert "from lid_opening import LidOpening, LidOpeningSettings" in test.read_text()
        assert "plugins:" in result.output

    def test_the_scaffold_is_found_by_a_comparison(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from polyzymd.config.comparison import ComparisonConfig

        root = _init(tmp_path, monkeypatch)
        monkeypatch.chdir(root)
        assert CliRunner().invoke(cli, ["new-analysis", "lid_opening"]).exit_code == 0
        comparison = root / "comparisons" / "c1" / "comparison.yaml"
        comparison.parent.mkdir()
        comparison.write_text(
            "name: c1\n"
            "conditions:\n"
            "  - label: A\n"
            "    config: ../../conditions/a/config.yaml\n"
            "    replicates: [1]\n"
            "plugins:\n"
            "  lid_opening:\n"
            "    selection: protein\n"
        )

        config = ComparisonConfig.from_yaml(comparison)

        assert config.plugins.get_enabled_plugins() == ["lid_opening"]

    def test_refuses_a_builtin_name(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        root = _init(tmp_path, monkeypatch)
        monkeypatch.chdir(root)

        result = CliRunner().invoke(cli, ["new-analysis", "rg"])

        assert result.exit_code != 0
        assert not (root / "analyses" / "rg.py").exists()

    def test_refuses_an_importable_module_name(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        root = _init(tmp_path, monkeypatch)
        monkeypatch.chdir(root)

        result = CliRunner().invoke(cli, ["new-analysis", "json"])

        assert result.exit_code != 0
        assert "importable Python module" in result.output

    def test_builtin_flag_writes_into_the_source_tree(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        root = _init(tmp_path, monkeypatch)
        repo = tmp_path / "repo"
        (repo / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (repo / "tests").mkdir()
        monkeypatch.chdir(root)

        result = CliRunner().invoke(
            cli, ["new-analysis", "lid_opening", "--builtin", "--project-root", str(repo)]
        )

        assert result.exit_code == 0, result.output
        assert (repo / "src" / "polyzymd" / "analyses" / "lid_opening.py").is_file()
        assert not (root / "analyses" / "lid_opening.py").exists()
