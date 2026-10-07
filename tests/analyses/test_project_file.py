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


def test_region_is_resolved_only_in_selections() -> None:
    """A label that says 'region' is left as written; a selection's region is resolved."""
    from polyzymd.analyses.project_file import resolve_names

    entry = {
        "label": "catalytic region distance",
        "selection": "region core",
        "pairs": [{"label": "lid region gap", "selection_a": "region core and name CA"}],
    }
    resolved = resolve_names(entry, {"core": "resid 1-5"}, {}, "test")
    assert resolved == {
        "label": "catalytic region distance",
        "selection": "(resid 1-5)",
        "pairs": [{"label": "lid region gap", "selection_a": "(resid 1-5) and name CA"}],
    }


@pytest.mark.parametrize("key", ["donors", "hydrogens", "acceptors"])
def test_region_is_resolved_in_hydrogen_bond_selections(key: str) -> None:
    """donors, hydrogens and acceptors of hydrogen_bonds are atom selections."""
    from polyzymd.analyses.project_file import resolve_names

    resolved = resolve_names({key: "region core and name N"}, {"core": "resid 1-5"}, {}, "test")
    assert resolved == {key: "(resid 1-5) and name N"}
