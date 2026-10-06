"""Regressions for the 1.3 audit round, wave D (friction and console output).

Each test names its finding in
``PAPERS/polyzymd_v1.3_refactor_handoff/audit_2026-10-06/AUDIT_LOG.md``.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner


def test_identical_warnings_for_several_conditions_are_one_line() -> None:
    """WS-5: eight conditions with one problem give one line naming them all."""
    from polyzymd.analyses.study_freeze import group_warnings

    warnings = [
        "A: no hashes; run x",
        "B: no hashes; run x",
        "metadata.doi is not set",
        "A replicate 1: odd",
    ]
    assert group_warnings(warnings, ["A", "B"]) == [
        "A, B: no hashes; run x",
        "metadata.doi is not set",
        "A replicate 1: odd",
    ]


@pytest.mark.parametrize(
    ("limits", "text"),
    [((1.9999995, 2.0000004), "1.999999 to 2"), ((0.7733, 1.227), "0.7733 to 1.227"), (None, "na")],
)
def test_a_narrow_interval_is_not_printed_as_one_number(limits, text) -> None:
    """WF-5: 'ci95 2 to 2' hid a real, narrow interval."""
    from polyzymd.analyses.protocols import _interval

    assert _interval(limits) == text


def test_short_runs_show_progress_in_ps() -> None:
    """SIM-9: a 4 ps run showed 0.0/0 ns."""
    from polyzymd.cli.status_report import _progress_text

    assert _progress_text(0.004, 0.004) == "   4.0/4ps"
    assert _progress_text(12.5, 100.0) == "  12.5/100ns"


def test_status_without_slurm_says_so_instead_of_warning() -> None:
    """SIM-9: --no-slurm printed 'squeue unavailable'."""
    from polyzymd.cli.status_report import render_agent

    text = render_agent([], slurm_available=False, slurm_queried=False)
    assert "SLURM not queried" in text and "squeue unavailable" not in text


def test_project_check_takes_production() -> None:
    """WF-1: the skill's advice to read production lengths works on a project too."""
    from polyzymd.cli.main import cli

    assert "--production" in CliRunner().invoke(cli, ["project", "check", "--help"]).output


def test_the_untestable_reason_is_the_real_one() -> None:
    """NOV-12: two replicates with no variance were said to lack replicates."""
    from polyzymd.analyses.protocols import ConditionReport, PairwiseReport, _verdict

    conditions = [
        ConditionReport(label=label, n_replicates=2, mean=0.0, replicate_values=[0.0, 0.0])
        for label in ("A", "B")
    ]
    pair = PairwiseReport(a="A", b="B", delta=0.0, testable=False)
    text = " ".join(_verdict("m", None, conditions, [pair]))
    assert "same value in every replicate" in text and "at least two" not in text


def test_a_missing_packmol_is_named(monkeypatch) -> None:
    """SIM-6: a missing Packmol gave 'Unexpected error (<class 'TypeError'>)'."""
    import shutil

    from polyzymd.utils import packmol

    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(OSError, match="Packmol is not on PATH"):
        packmol._require_packmol()
