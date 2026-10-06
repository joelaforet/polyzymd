"""project freeze: what each study's manifest in a project lists."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.analyses.test_project import project  # noqa: F401  (fixture)

pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore"), pytest.mark.usefixtures("git_identity")]


def test_project_study_manifests_list_what_git_tracks(project: Path) -> None:  # noqa: F811
    """A study's manifest in a project lists the files the deposit holds."""
    from polyzymd.analyses.project_freeze import freeze_project
    from polyzymd.analyses.study_git import init_repository

    (project / ".gitignore").write_text("__pycache__/\ndata.local.yaml\n")
    init_repository(project, "start")
    (project / "lipa" / "__pycache__").mkdir()
    (project / "lipa" / "__pycache__" / "x.pyc").write_bytes(b"x")
    freeze_project(project)
    manifest = json.loads((project / "lipa" / "manifest.json").read_text())
    assert not any("__pycache__" in name for name in manifest["files"])
    assert manifest["git"]["commit"]
