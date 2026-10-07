"""The SLURM pages of the docs agree with the code.

The preset table in ``docs/source/how_to/hpc_slurm.md`` is compared with the
presets in :mod:`polyzymd.workflow.slurm`. One docs lint keeps the CU Boulder
site notes on their own page.
"""

from __future__ import annotations

import typing
from pathlib import Path

import pytest

from polyzymd.workflow.slurm import PresetType, SlurmConfig

DOCS = Path(__file__).resolve().parents[1] / "docs" / "source"


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


@pytest.mark.parametrize(
    "site_fact", ["slurm/blanca", "bgpu-bortz1", "bgpu-g4-u2", "acceptance test"]
)
def test_docs_lint_cu_boulder_site_notes_are_on_one_page(site_fact: str) -> None:
    pages = [path.name for path in DOCS.rglob("*.md") if site_fact in path.read_text()]
    assert pages == ["site_cu_boulder.md"]
