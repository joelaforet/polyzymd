"""project init and project add-study: the folders and files they write."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.project_scaffold import add_study, create_project


def test_add_study_git_ignores_runs_in_a_project_made_before_runs_existed(tmp_path: Path) -> None:
    """A project .gitignore without runs/ gets it when a study is added."""
    root = create_project(tmp_path / "paper", ["lipa"], git=False).root
    gitignore = root / ".gitignore"
    gitignore.write_text("data.local.yaml\n")

    add_study(root, "calb")

    assert gitignore.read_text().splitlines() == ["data.local.yaml", "runs/"]


@pytest.mark.parametrize("label", ["runs", "slurm", "stats", "analyses", "environment"])
def test_a_study_may_not_take_a_folder_name_polyzymd_uses(tmp_path: Path, label: str) -> None:
    """A study named runs/ or stats/ would collide with the project's own folders."""
    with pytest.raises(ProtocolError, match="reserved"):
        create_project(tmp_path / "new", [label], git=False)
    root = create_project(tmp_path / "paper", ["lipa"], git=False).root
    with pytest.raises(ProtocolError, match="reserved"):
        add_study(root, label)
