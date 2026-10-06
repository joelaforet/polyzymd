"""project.yaml: where a project's studies may sit and what an analysis entry must be."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.project_file import load_project_file


@pytest.mark.parametrize("where", ["nested/lipa", "../outside/lipa"])
def test_a_study_must_sit_directly_in_the_project(tmp_path: Path, where: str) -> None:
    """A study elsewhere never finds its paper, or is not published."""
    paper = tmp_path / "Paper"
    paper.mkdir()
    folder = (paper / where).resolve()
    folder.mkdir(parents=True)
    (folder / "study.yaml").write_text("conditions: {}\n")
    (paper / "project.yaml").write_text(f"studies: {{lipa: {where}}}\n")
    with pytest.raises(ProtocolError, match="directly inside the project"):
        load_project_file(paper)


def test_an_analysis_entry_must_be_a_mapping(tmp_path: Path) -> None:
    """rg: notamapping gets a message, not a traceback."""
    paper = tmp_path / "Paper"
    (paper / "lipa").mkdir(parents=True)
    (paper / "lipa" / "study.yaml").write_text("conditions: {}\n")
    (paper / "project.yaml").write_text("studies: {lipa: lipa}\nanalyses: {rg: notamapping}\n")
    with pytest.raises(ProtocolError, match="must be a mapping"):
        load_project_file(paper)
