"""A comparison.yaml that still configures rg or rmsd keeps working for every other plugin.

rg and rmsd left the plugin system for ``polyzymd analyze``. Each of their
``plugins`` blocks is ignored with one warning, and ``compare run rg`` or
``compare run rmsd`` prints the equivalent ``polyzymd analyze`` command.
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
    """A comparison.yaml with rg, rmsd and rmsf blocks over two conditions."""
    conditions = []
    for label in ("No Polymer", "SBMA 50"):
        config = write_simulation_config(tmp_path / label.replace(" ", "_"), scratch=tmp_path)
        conditions.append({"label": label, "config": str(config), "replicates": [1, 2]})
    data = {
        "name": "retired_plugins",
        "control": "No Polymer",
        "conditions": conditions,
        "defaults": {"equilibration_time": "200ns"},
        "plugins": {
            "rg": {"runs": [{"label": "Protein", "selection": "protein"}]},
            "rmsd": {"runs": [{"label": "CA", "selection": "name CA"}]},
            "rmsf": {"selection": "name CA"},
        },
    }
    path = tmp_path / "comparison.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def test_retired_blocks_are_ignored_with_one_warning_each(comparison_file: Path) -> None:
    """The file loads, rmsf keeps its settings, and rg and rmsd are gone with a warning each."""
    with pytest.warns(UserWarning, match="block, which is ignored") as record:
        config = ComparisonConfig.from_yaml(comparison_file)

    messages = [str(item.message) for item in record]
    for name, function in (("rg", "radius_of_gyration"), ("rmsd", "rmsd")):
        (message,) = [text for text in messages if f"plugins.{name} block" in text]
        assert f"polyzymd analyze {name} " in message
        assert f"polyzymd.analyses.functions.{function}." in message
    assert config.plugins.get_enabled_plugins() == ["rmsf"]
    assert config.plugins.get("rmsf").selection == "name CA"
    assert config.validate_config() == []


def test_compare_validate_and_run_rmsf_succeed(comparison_file: Path, monkeypatch) -> None:
    """compare validate passes and compare run rmsf reaches the pipeline with rmsf settings."""
    seen = {}

    def _pipeline(analysis, config, **kwargs):
        seen["analysis"], seen["plugins"] = analysis.name, config.plugins.get_enabled_plugins()
        raise ValueError("stop after the pipeline was reached")

    monkeypatch.setattr("polyzymd.analyses.orchestrator.run_comparison", _pipeline)
    runner = CliRunner()

    validated = runner.invoke(compare, ["validate", "-f", str(comparison_file)])
    assert validated.exit_code == 0, validated.output

    runner.invoke(compare, ["run", "rmsf", "-f", str(comparison_file)])
    assert seen == {"analysis": "rmsf", "plugins": ["rmsf"]}


@pytest.mark.parametrize("name", ["rg", "rmsd"])
def test_compare_run_prints_the_analyze_command(comparison_file: Path, name: str) -> None:
    """compare run rg or rmsd exits 1 with the polyzymd analyze command for this file."""
    result = CliRunner().invoke(compare, ["run", name, "-f", str(comparison_file)])

    assert result.exit_code == 1
    assert "no longer runs through 'polyzymd compare run'" in result.stderr
    assert "Unknown comparison type" not in result.stderr
    fix = next(line for line in result.stderr.splitlines() if line.startswith("fix: "))
    assert fix.startswith(f"fix: polyzymd analyze {name} -c ")
    assert "--label 'No Polymer'" in fix and "--label 'SBMA 50'" in fix
    assert fix.endswith("--replicates 1,2 --eq 200ns")


def test_compare_run_all_skips_retired(comparison_file, monkeypatch, tmp_path: Path) -> None:
    """compare run-all warns about rg and rmsd and runs only the plugins that remain."""
    seen = {}

    def _run_all(config, **kwargs):
        seen["plugins"] = config.plugins.get_enabled_plugins()
        return {"rmsf": {"comparison": {"ok": True}, "comparison_path": tmp_path / "r.json"}}

    monkeypatch.setattr("polyzymd.analyses.orchestrator.run_all_comparisons", _run_all)
    with pytest.warns(UserWarning, match="plugins.rmsd block, which is ignored"):
        result = CliRunner().invoke(compare, ["run-all", "-f", str(comparison_file)])

    assert result.exit_code == 0, result.output
    assert seen == {"plugins": ["rmsf"]}
