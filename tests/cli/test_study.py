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
