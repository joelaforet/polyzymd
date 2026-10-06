"""Study.replicate_table: one row per replicate and quantity, filled as the report fills it."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses.study_scaffold import create_study
from polyzymd.cli.main import cli

pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore"), pytest.mark.usefixtures("git_identity")]


def _analyze(*arguments: str):
    return CliRunner().invoke(cli, ["analyze", *arguments], catch_exceptions=False)


def _json(result) -> dict:
    return json.loads(result.output[result.output.index("{") :])


def _hbond_study(tmp_path: Path, analyses: str, never_egm: bool = False) -> Path:
    from tests.analyses.test_hydrogen_bonds_analyze import _schedule, _write

    schedules = {
        (label, r): _schedule(10 * r + offset) for label, offset in (("A", 1), ("B", 5)) for r in (1, 2, 3)
    }
    if never_egm:
        schedules[("A", 3)] = _schedule(31, never_egm=True)
    configs = _write(tmp_path / "runs", schedules)
    root = tmp_path / "st"
    create_study(root, conditions=configs, equilibration="0ns")
    text = (root / "study.yaml").read_text().replace("analyses: {}", "analyses:\n" + analyses)
    (root / "study.yaml").write_text(text)
    subprocess.run(["git", "-C", str(root), "commit", "-qam", "analyses"], capture_output=True)
    return root


def test_replicate_table_keeps_each_quantity(tmp_path: Path) -> None:
    """Two hydrogen-bond summaries in one run are two quantities, not two replicates."""
    root = _hbond_study(
        tmp_path,
        "  hbonds:\n    analysis: hydrogen_bonds\n"
        "    groups: {protein: chainid A, polymer: chainid C}\n"
        "    summaries:\n      protein_polymer: {between: [protein, polymer]}\n"
        "      protein_protein: {within: protein}\n",
    )
    options = ["--study", str(root), "--no-plots", "--no-eq-check"]
    assert _analyze("hbonds", *options).exit_code == 0
    assert _analyze("hbonds", *options, "--run", "protein_protein_mean_hbonds").exit_code == 0
    table = pz.Study(root).replicate_table("hbonds")
    keys = ["condition", "replicate", "name", "part", "label"]
    assert table.groupby(keys, dropna=False).size().max() == 1
    assert set(table["name"]) == {
        "hydrogen_bonds_protein_polymer",
        "hydrogen_bonds_protein_protein",
    }


def test_replicate_table_fills_missing_labels_as_the_report(tmp_path: Path) -> None:
    """A pair one replicate never formed counts as 0 in the table, as in the report."""
    root = _hbond_study(
        tmp_path,
        "  hbonds: {analysis: hydrogen_bonds, groups: {protein: chainid A, polymer: chainid C}}\n",
        never_egm=True,
    )
    report = _json(
        _analyze(
            "hbonds", "--study", str(root), "--no-plots", "--no-eq-check",
            "--run", "protein_polymer_pairs", "--format", "json",
        )
    )
    table = pz.Study(root).replicate_table("hbonds")
    for row in report["conditions"]:
        sub = table[(table.condition == row["label"]) & (table.label.astype(str) == str(row["entry"]))]
        assert len(sub) == row["n_replicates"], row["entry"]
        assert sub.value.mean() == pytest.approx(row["mean"]), row["entry"]


def test_equal_condition_means_read_as_no_trend() -> None:
    """A flat trend says the means are the same, not that it cannot be tested."""
    from types import SimpleNamespace

    from polyzymd.analyses.protocols import ConditionReport, _trend_line
    from polyzymd.analyses.study_statistics import trend_sentence, trend_tests

    conditions = [
        ConditionReport(label=f"c{i}", n_replicates=2, mean=1.0, replicate_values=[0.9, 1.1])
        for i in range(3)
    ]
    factors = {"c0": {"x": 0.1}, "c1": {"x": 0.5}, "c2": {"x": 0.9}}
    (flat,) = trend_tests(SimpleNamespace(conditions=conditions), factors)
    assert trend_sentence("m", None, flat).startswith(
        "no trend of m with x: every condition mean is the same"
    )
    assert _trend_line(flat).startswith("trend x  no trend: every condition mean is the same")
