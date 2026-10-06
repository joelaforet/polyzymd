"""polyzymd study locate: mapping downloaded runs to a study's conditions."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from polyzymd.analyses.study_file import load_study_file
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_committed_study

pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore"), pytest.mark.usefixtures("git_identity")]


def test_alike_runs_are_told_apart_by_folder_name_or_refused(tmp_path: Path) -> None:
    """Two conditions whose runs share a name are never mapped to one folder."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    download = tmp_path / "download"
    shutil.copytree(tmp_path / "scratch", download / "by_condition")
    result = CliRunner().invoke(cli, ["study", "locate", str(download), "--study", str(root)])
    assert result.exit_code == 0, result.output
    data = load_study_file(root).data
    assert data["Polymer"].name == "polymer" and data["No polymer"].name == "no_polymer"
    mixed = tmp_path / "mixed"
    shutil.copytree(tmp_path / "scratch" / "polymer", mixed / "a")
    shutil.copytree(tmp_path / "scratch" / "no_polymer", mixed / "b")
    refused = CliRunner().invoke(cli, ["study", "locate", str(mixed), "--study", str(root)])
    assert refused.exit_code == 2 and "cannot be told apart" in refused.output


def _write_manifest(root: Path, folders: dict[str, Path]) -> None:
    """Write a manifest.json listing the size of every run file of each condition."""
    import json

    conditions = {
        label: {
            "replicates": {
                "1": {
                    "files": [
                        {"path": str(path.relative_to(folder)), "size": path.stat().st_size}
                        for path in sorted(folder.rglob("*"))
                        if path.is_file()
                    ]
                }
            }
        }
        for label, folder in folders.items()
    }
    (root / "manifest.json").write_text(json.dumps({"conditions": conditions}))


def test_runs_of_equal_size_go_to_the_folder_named_for_each_condition(tmp_path: Path) -> None:
    """Folders named for the conditions decide before file sizes; one folder never serves two."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    scratch = tmp_path / "scratch"
    _write_manifest(root, {"No polymer": scratch / "no_polymer", "Polymer": scratch / "polymer"})
    download = tmp_path / "download"
    shutil.copytree(scratch, download)
    result = CliRunner().invoke(cli, ["study", "locate", str(download), "--study", str(root)])
    assert result.exit_code == 0, result.output
    data = load_study_file(root).data
    assert data["Polymer"].name == "polymer" and data["No polymer"].name == "no_polymer"
    tied = tmp_path / "tied"
    shutil.copytree(scratch / "no_polymer", tied / "a")
    shutil.copytree(scratch / "polymer", tied / "b")
    refused = CliRunner().invoke(cli, ["study", "locate", str(tied), "--study", str(root)])
    assert refused.exit_code == 2, refused.output
    data = load_study_file(root).data
    assert data["Polymer"] != data["No polymer"]


def test_locate_keeps_the_data_file_when_nothing_is_located(tmp_path: Path) -> None:
    """Entries written by hand survive a locate that finds no condition, and nothing is written."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    text = "# mine\nNo polymer: ../scratch/no_polymer\nPolymer: ../scratch/polymer\n"
    (root / "data.local.yaml").write_text(text)
    empty = tmp_path / "empty"
    empty.mkdir()
    result = CliRunner().invoke(cli, ["study", "locate", str(empty), "--study", str(root)])
    assert result.exit_code == 2
    assert "wrote" not in result.output
    assert (root / "data.local.yaml").read_text() == text


def test_add_condition_names_where_its_runs_are(tmp_path: Path) -> None:
    """add-condition prints the data folder it records and warns when it holds no runs."""
    from tests._support.analysis_testkit import write_simulation_config

    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    config = write_simulation_config(tmp_path / "runs" / "new", scratch=tmp_path / "nothing")
    (config.parent / "test.pdb").write_text("REMARK input\nEND\n")
    (tmp_path / "nothing").mkdir()
    result = CliRunner().invoke(
        cli, ["study", "add-condition", "New", "--config", str(config), "--study", str(root)]
    )
    assert result.exit_code == 0, result.output
    where = (tmp_path / "nothing").resolve()
    assert f"data New: {where} (from the config's scratch_directory)" in result.output
    assert "warning" in result.output and "no run" in result.output
    help_text = CliRunner().invoke(cli, ["study", "add-condition", "--help"]).output
    assert "data.local.yaml" in help_text
