"""A comparison.yaml that still configures rg keeps working for every other plugin.

rg left the plugin system for ``polyzymd analyze rg``. Its ``plugins.rg`` block
is ignored with one warning, and ``compare run rg`` prints the equivalent
``polyzymd analyze rg`` command.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from polyzymd.cli.compare import compare
from polyzymd.config.comparison import ComparisonConfig
from tests._support.analysis_testkit import write_simulation_config


@pytest.fixture()
def comparison_file(tmp_path: Path) -> Path:
    """A comparison.yaml with an rg block and an rmsd block over two conditions."""
    conditions = []
    for label in ("No Polymer", "SBMA 50"):
        config = write_simulation_config(tmp_path / label.replace(" ", "_"), scratch=tmp_path)
        conditions.append({"label": label, "config": str(config), "replicates": [1, 2]})
    data = {
        "name": "retired_rg",
        "control": "No Polymer",
        "conditions": conditions,
        "defaults": {"equilibration_time": "200ns"},
        "plugins": {
            "rg": {"runs": [{"label": "Protein", "selection": "protein"}]},
            "rmsd": {"runs": [{"label": "CA", "selection": "name CA"}]},
        },
    }
    path = tmp_path / "comparison.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def test_rg_block_is_ignored_with_one_warning(comparison_file: Path) -> None:
    """The file loads, rmsd keeps its settings, and rg is gone with a warning."""
    with pytest.warns(UserWarning, match="plugins.rg block, which is ignored") as record:
        config = ComparisonConfig.from_yaml(comparison_file)

    assert len([item for item in record if "plugins.rg" in str(item.message)]) == 1
    assert "polyzymd analyze rg" in str(record[0].message)
    assert config.plugins.get_enabled_plugins() == ["rmsd"]
    assert config.plugins.get("rmsd").runs[0].selection == "name CA"
    assert config.validate_config() == []


def test_compare_validate_and_run_rmsd_succeed(comparison_file: Path, monkeypatch) -> None:
    """compare validate passes and compare run rmsd reaches the pipeline with rmsd settings."""
    seen = {}

    def _pipeline(analysis, config, **kwargs):
        seen["analysis"], seen["plugins"] = analysis.name, config.plugins.get_enabled_plugins()
        raise ValueError("stop after the pipeline was reached")

    monkeypatch.setattr("polyzymd.analyses.orchestrator.run_comparison", _pipeline)
    runner = CliRunner()

    validated = runner.invoke(compare, ["validate", "-f", str(comparison_file)])
    assert validated.exit_code == 0, validated.output

    runner.invoke(compare, ["run", "rmsd", "-f", str(comparison_file)])
    assert seen == {"analysis": "rmsd", "plugins": ["rmsd"]}


def test_compare_run_rg_prints_the_analyze_command(comparison_file: Path) -> None:
    """compare run rg exits 1 with the polyzymd analyze rg command for this file."""
    result = CliRunner().invoke(compare, ["run", "rg", "-f", str(comparison_file)])

    assert result.exit_code == 1
    assert "no longer runs through 'polyzymd compare run'" in result.stderr
    assert "Unknown comparison type" not in result.stderr
    fix = next(line for line in result.stderr.splitlines() if line.startswith("fix: "))
    assert fix.startswith("fix: polyzymd analyze rg -c ")
    assert "--label 'No Polymer'" in fix and "--label 'SBMA 50'" in fix
    assert fix.endswith("--replicates 1,2 --eq 200ns")


def test_compare_run_all_skips_rg(comparison_file: Path, monkeypatch, tmp_path: Path) -> None:
    """compare run-all warns about rg once and runs only the plugins that remain."""
    seen = {}

    def _run_all(config, **kwargs):
        seen["plugins"] = config.plugins.get_enabled_plugins()
        return {"rmsd": {"comparison": {"ok": True}, "comparison_path": tmp_path / "r.json"}}

    monkeypatch.setattr("polyzymd.analyses.orchestrator.run_all_comparisons", _run_all)
    with pytest.warns(UserWarning, match="plugins.rg block, which is ignored"):
        result = CliRunner().invoke(compare, ["run-all", "-f", str(comparison_file)])

    assert result.exit_code == 0, result.output
    assert seen == {"plugins": ["rmsd"]}
