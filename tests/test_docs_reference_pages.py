"""Facts in the CLI, SLURM and GROMACS pages of the docs agree with the code.

Each test reads the Markdown source under ``docs/source`` and compares one
stated fact with the code or with the other pages, so a page that drifts from
the code fails here.
"""

from __future__ import annotations

import re
import typing
from pathlib import Path

import pytest

from polyzymd.workflow.slurm import PresetType, SlurmConfig

DOCS = Path(__file__).resolve().parents[1] / "docs" / "source"


def _pages() -> dict[Path, str]:
    return {path: path.read_text() for path in DOCS.rglob("*.md")}


def _preset_rows(text: str) -> dict[str, list[str]]:
    rows = {}
    for line in text.splitlines():
        cells = [cell.strip().strip("`") for cell in line.strip().strip("|").split("|")]
        if line.startswith("| `") and cells[0] in typing.get_args(PresetType):
            rows[cells[0]] = cells
    return rows


def test_slurm_preset_table_matches_the_presets() -> None:
    rows = _preset_rows((DOCS / "how_to" / "hpc_slurm.md").read_text())
    assert set(rows) == set(typing.get_args(PresetType))
    for name, (_, partition, qos, time_limit, _use) in rows.items():
        preset = SlurmConfig.from_preset(name)
        assert partition == preset.partition, name
        assert qos == (preset.qos or "none"), name
        assert time_limit == preset.time_limit, name


def test_cli_reference_has_no_second_preset_table() -> None:
    assert _preset_rows((DOCS / "reference" / "cli_reference.md").read_text()) == {}


def test_study_check_table_has_no_paragraph_inside() -> None:
    text = (DOCS / "reference" / "cli_reference.md").read_text()
    section = text.split("### polyzymd study check", 1)[1]
    table = section[section.index("| Line | Fields |") :].split("\n\n", 1)[0]
    assert "| citation |" in table and "| warning |" in table


def test_no_page_says_gromacs_job_scripts_ignore_stop() -> None:
    for path, text in _pages().items():
        flat = " ".join(text.split())
        assert not re.search(r"GROMACS job scripts,? do not check", flat), path


@pytest.mark.parametrize(
    "site_fact", ["slurm/blanca", "bgpu-bortz1", "bgpu-g4-u2", "acceptance test"]
)
def test_cu_boulder_site_notes_are_on_one_page(site_fact: str) -> None:
    pages = [path.name for path, text in _pages().items() if site_fact in text]
    assert pages == ["site_cu_boulder.md"]
