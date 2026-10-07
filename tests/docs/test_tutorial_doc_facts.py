"""Facts that the tutorials and the pages around them must state correctly."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs" / "source"


def _read(relative: str) -> str:
    return (DOCS / relative).read_text(encoding="utf-8")


def _all_pages() -> dict[str, str]:
    return {
        str(path.relative_to(DOCS)): path.read_text(encoding="utf-8")
        for path in sorted(DOCS.rglob("*.md"))
    }


@pytest.mark.parametrize(
    "page", ["first_analysis.md", "analysis_complete_workflow.md", "sasa_analysis.md"]
)
def test_tutorial_analyze_commands_use_the_study(page: str) -> None:
    text = _read(f"tutorials/{page}")
    # Only the quick-look section of the first analysis lesson may use -c.
    text = text.split("## When a quick look with `-c` is enough")[0]
    commands = [line for line in text.splitlines() if line.startswith("polyzymd analyze")]
    assert commands, page
    for command in commands:
        assert "--study" in command or "--project" in command, command


def test_quickstart_is_in_the_tutorials_toctree() -> None:
    toctree = _read("tutorials/index.md").split("```{toctree}")[1]
    assert "../get_started/quickstart" in toctree


def test_pages_say_that_build_runs_analyze() -> None:
    assert "The `build` environment includes the analysis tools" in _read("get_started/index.md")
    assert "`polyzymd analyze` runs in the `build` environment" in _read(
        "get_started/installation.md"
    )
    skill = (ROOT / ".claude" / "skills" / "polyzymd-analyze" / "SKILL.md").read_text()
    assert "only the `analysis` and" not in skill


def test_lipase_catalytic_residues_use_built_topology_numbers() -> None:
    crystal = re.compile(
        r"(Ser|Asp|His) ?(77|133|156)\b|resid (77|133|156) and name (OG|OD|NE2|ND1)"
    )
    for name, text in _all_pages().items():
        for line in text.splitlines():
            if crystal.search(line):
                assert "1ISP" in line, f"{name}: {line}"


def test_project_reference_lists_its_methods() -> None:
    from polyzymd.analyses.project import Project

    reference = _read("reference/study_api.md")
    public = [name for name in vars(Project) if not name.startswith("_")]
    assert public
    for name in public:
        assert f"project.{name}" in reference, name


def test_stats_script_finds_the_project_from_its_own_path() -> None:
    project = _read("how_to/project.md")
    assert "root = Path(__file__).parents[1]" in project
    assert 'to_csv("stats/' not in project


@pytest.mark.parametrize("page", ["reference/configuration.md", "how_to/equilibration.md"])
def test_equilibration_cannot_be_skipped_is_stated(page: str) -> None:
    text = " ".join(_read(page).split())
    assert "Equilibration cannot be skipped" in text
    assert "with a duration above 0 ns" in text


def test_hydrogen_bond_how_to_links_to_the_rules() -> None:
    how_to = _read("how_to/hydrogen_bonds.md")
    for moved in ("| Donor hydrogen |", "| Engine |", "hydrogen_bonds(group_a, group_b=None"):
        assert moved not in how_to
    assert "copolymer" in _read("explanation/analysis_hydrogen_bonds_verification.md")


def test_docs_say_replicates_for_replicate_folders() -> None:
    for name, text in _all_pages().items():
        assert re.search(r": runs \[\d", text) is None, name
    assert "a project laid out like this" not in _read("tutorials/analysis_complete_workflow.md")


def test_docs_do_not_describe_comparison_configs() -> None:
    for name, text in _all_pages().items():
        assert "comparison.yaml" not in text, name
        assert "comparison configuration" not in text, name


def test_how_to_pages_do_not_call_themselves_tutorials() -> None:
    for path in sorted((DOCS / "how_to").glob("*.md")):
        assert "This tutorial" not in path.read_text(encoding="utf-8"), path.name


def test_unit_examples_use_the_report_spelling() -> None:
    for name, text in _all_pages().items():
        assert re.search(r'unit="Å', text) is None, name
        assert re.search(r"unit: Å", text) is None, name


def test_contributor_material_is_in_the_contributor_guide() -> None:
    assert not (DOCS / "explanation" / "architecture.md").exists()
    assert (DOCS / "contributor_guide" / "architecture.md").exists()
    assert "Regenerate `pixi.lock`" not in _read("how_to/hardware_platforms.md")
    assert "Regenerate `pixi.lock`" in _read("contributor_guide/packaging.md")
