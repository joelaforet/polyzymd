"""Execution cost hints an analysis exposes to the CLI and the SLURM workflow."""

from __future__ import annotations

import pytest

from polyzymd.analyses.base import Analysis


def test_execution_cost_hint_default() -> None:
    """An analysis that states no cost is treated as medium."""
    assert Analysis.execution_cost_hint == "medium"


@pytest.mark.parametrize("analysis_name", ["sasa", "contacts"])
def test_execution_cost_hint_high(analysis_name: str) -> None:
    """The expensive built-ins say so, which contract_analysis copies from the plugin."""
    from polyzymd.analyses.discovery import get_analysis

    assert get_analysis(analysis_name).execution_cost_hint == "high"
