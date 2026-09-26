"""Tests for ``polyzymd analyze rg`` and the radius of gyration function it runs.

Each replicate is an OpenMM run directory with one DCD segment of four
unit-mass atoms on a cross, scaled on frame ``k`` so that its radius of
gyration is ``base + 0.01 * k``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from polyzymd.analyses.functions import radius_of_gyration
from polyzymd.analyses.protocols import ProtocolReport
from polyzymd.cli.analyze import EXIT_ANALYSIS_ERROR, analyze_command
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

EQUILIBRATION = "0.25ns"


@pytest.fixture()
def configs(tmp_path: Path) -> dict[str, Path]:
    """Two conditions of three replicates, B larger than A by one Å."""
    paths = {}
    for label, offset in (("A", 1.0), ("B", 2.0)):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            scales = [offset + 0.1 * replicate + 0.01 * k for k in range(10)]
            write_openmm_replicate(config, replicate, scales)
        paths[label] = config
    return paths


def test_radius_of_gyration_is_mass_weighted() -> None:
    """Masses 3 and 1 at x = 1 and x = -3 put the centre at 0 and give Rg sqrt(3)."""
    import MDAnalysis as mda

    universe = mda.Universe.empty(2, trajectory=True)
    universe.add_TopologyAttr("masses", [3.0, 1.0])
    universe.atoms.positions = np.array([[1.0, 0.0, 0.0], [-3.0, 0.0, 0.0]])
    assert radius_of_gyration(universe.atoms) == pytest.approx(np.sqrt(3.0))


def test_cli_rg_prints_the_report(configs: dict[str, Path], tmp_path: Path) -> None:
    """'polyzymd analyze rg' measures protein Rg by default and compares with the first config."""
    arguments = ["rg", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    arguments += ["--set", "selection=all", "--output-dir", str(tmp_path)]
    result = CliRunner().invoke(analyze_command, arguments)
    assert result.exit_code == 0, result.output
    lines = result.stdout.strip().split("\n")
    assert lines[0].startswith("# polyzymd analyze rg  metric mean_rg  unit A  eq 0.25ns")
    assert lines[1].startswith("A  n 3  mean 1.26")
    assert lines[2].startswith("B  n 3  mean 2.26")
    assert lines[3].startswith("A vs B  delta +1  ")
    assert lines[-1].startswith("verdict: B larger mean_rg than A")

    report = CliRunner().invoke(analyze_command, [*arguments, "--format", "json"])
    assert ProtocolReport.model_validate_json(report.stdout).frames_per_replicate == {
        "A": [7, 7, 7],
        "B": [7, 7, 7],
    }
    refused = CliRunner().invoke(analyze_command, [*arguments, "--run", "Protein"])
    assert refused.exit_code == EXIT_ANALYSIS_ERROR
    default = CliRunner().invoke(
        analyze_command, ["rg", "-c", str(configs["A"]), "--eq", EQUILIBRATION]
    )
    assert default.exit_code == EXIT_ANALYSIS_ERROR
    assert "matched no atoms" in default.stderr


def test_cli_rg_with_one_config_prints_a_summary(configs: dict[str, Path], tmp_path) -> None:
    """One -c gives the condition line and its verdict and no comparison."""
    arguments = ["rg", "-c", str(configs["A"]), "--eq", EQUILIBRATION, "--set", "selection=all"]
    result = CliRunner().invoke(analyze_command, [*arguments, "--output-dir", str(tmp_path)])
    assert result.exit_code == 0, result.output
    lines = result.stdout.strip().split("\n")
    # Seven frames per replicate are too few to judge the equilibrated start.
    (note,) = [line for line in lines if line.startswith("warning:")]
    assert "replicates 1, 2, 3 have fewer than 20 effective samples" in note
    quiet = CliRunner().invoke(
        analyze_command, [*arguments, "--output-dir", str(tmp_path), "--no-eq-check"]
    )
    assert "warning:" not in quiet.stdout and "eq_detected" not in quiet.stdout
    lines = [line for line in lines if not line.startswith("warning:")]
    assert len(lines) == 3
    assert lines[1].startswith("A  n 3  mean 1.26")
    assert lines[2].startswith("verdict: A mean_rg 1.26 A (95% CI")


def test_python_analyze_rg_gives_the_cli_report(configs: dict[str, Path], tmp_path) -> None:
    """protocols.analyze("rg", ...) runs the same path as the CLI and returns its report."""
    from polyzymd.analyses import analyze

    report = analyze(
        "rg",
        [configs["A"], configs["B"]],
        equilibration=EQUILIBRATION,
        labels=["control", "treated"],
        replicates=[1, 2],
        settings={"selection": "all"},
        output_dir=tmp_path,
    )

    assert report.analysis == "rg" and report.metric == "mean_rg" and report.unit == "A"
    assert [row.label for row in report.conditions] == ["control", "treated"]
    assert report.conditions[0].replicate_values == pytest.approx([1.16, 1.26], abs=1e-5)
    assert report.frames_per_replicate == {"control": [7, 7], "treated": [7, 7]}
    (row,) = report.pairwise
    assert (row.a, row.b, row.test) == ("control", "treated", "welch_t")
    assert row.delta == pytest.approx(1.0, abs=1e-5)
    assert (tmp_path / "polyzymd_results/rg/treated/replicate_2/record.json").is_file()

    arguments = ["rg", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    arguments += ["--label", "control", "--label", "treated", "--replicates", "1-2"]
    arguments += ["--set", "selection=all", "--output-dir", str(tmp_path), "--format", "json"]
    result = CliRunner().invoke(analyze_command, arguments)
    assert ProtocolReport.model_validate_json(result.stdout) == report


def test_python_analyze_rg_refuses_other_settings_and_names_rg(configs) -> None:
    """Unknown settings are refused, and an unknown analysis lists rg among the names."""
    from polyzymd.analyses import analyze
    from polyzymd.analyses.exceptions import ProtocolError

    with pytest.raises(ProtocolError, match="no setting other than selection"):
        analyze("rg", [configs["A"]], equilibration=EQUILIBRATION, settings={"runs": []})
    with pytest.raises(ProtocolError) as excinfo:
        analyze("not_an_analysis", [configs["A"]])
    assert "rg" in excinfo.value.hint.split("Use one of: ")[1].split(", ")
