"""study.yaml: the schema, Study("study.yaml"), analyze --study, study check and results().

The study folder holds two conditions of two synthetic OpenMM replicates
each (four unit-mass atoms on a cross): frame k of replicate r of a condition
with offset o has a radius of gyration of o + 0.1 r + 0.01 k, so after the
0.25 ns window (frames 3 to 9) its mean is o + 0.1 r + 0.06.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.results import read_results
from polyzymd.analyses.study_file import load_study_file
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

OFFSETS = {"No polymer": 1.0, "Polymer": 2.0}


def _write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


@pytest.fixture()
def study_dir(tmp_path: Path) -> Path:
    """A study folder: conditions/<label>/config.yaml, data under data/, and study.yaml."""
    root = tmp_path / "my_study"
    for folder, offset in (("no_polymer", 1.0), ("polymer", 2.0)):
        config = write_simulation_config(
            root / "conditions" / folder, scratch=tmp_path / "data" / folder
        )
        for replicate in (1, 2):
            write_openmm_replicate(
                config, replicate, [offset + 0.1 * replicate + 0.01 * k for k in range(10)]
            )
    _write(
        root / "study.yaml",
        "polyzymd: 1.3.0\n"
        "equilibration: 0.25ns\n"
        "conditions:\n"
        "  No polymer: conditions/no_polymer/config.yaml\n"
        "  Polymer: conditions/polymer/config.yaml\n"
        "analyses:\n"
        "  rg: {selection: all}\n"
        "  rg_first:\n"
        "    analysis: rg\n"
        "    selection: index 0\n",
    )
    return root


def _analyze(*arguments: str) -> object:
    return CliRunner().invoke(cli, ["analyze", *arguments, "--no-eq-check", "--no-plots"])


class TestSchema:
    def test_reads_conditions_relative_to_the_file(self, study_dir: Path) -> None:
        protocol = load_study_file(study_dir)
        assert list(protocol.conditions) == ["No polymer", "Polymer"]
        assert protocol.conditions["Polymer"] == study_dir / "conditions/polymer/config.yaml"
        assert protocol.analyses["rg_first"].analysis == "rg"
        assert protocol.analyses["rg_first"].settings == {"selection": "index 0"}
        assert protocol.results_dir("rg") == study_dir / "results" / "rg"

    @pytest.mark.parametrize(
        ("text", "message", "hint"),
        [
            (
                "equilibrium: 1ns\nconditions: {A: a.yaml}\n",
                "unknown key 'equilibrium'",
                "equilibration",
            ),
            ("conditions: {A: a.yaml}\n", "no 'equilibration'", None),
            ("equilibration: 1ns\n", "no 'conditions'", None),
            ("equilibration: soon\nconditions: {A: a.yaml}\n", "cannot read equilibration", None),
            (
                "equilibration: 1ns\nconditions: {A: a.yaml}\nanalyses: {rg: {selction: all}}\n",
                "unknown key 'selction'",
                "selection",
            ),
            (
                "equilibration: 1ns\nconditions: {A: a.yaml}\nanalyses: {rgg: {}}\n",
                "does not ship",
                "'rg'",
            ),
            ("equilibration: 1ns\nconditions: {A: a.yaml}\nstride: 0\n", "stride", None),
        ],
    )
    def test_refuses_with_a_hint(self, tmp_path: Path, text: str, message: str, hint) -> None:
        with pytest.raises(ProtocolError, match=message) as caught:
            load_study_file(_write(tmp_path / "study.yaml", text))
        if hint:
            assert hint in caught.value.hint

    def test_replicate_range(self, tmp_path: Path) -> None:
        path = _write(
            tmp_path / "study.yaml",
            "equilibration: 1ns\nconditions: {A: a.yaml}\nreplicates: 1-3\n",
        )
        assert load_study_file(path).replicates == [1, 2, 3]


class TestStudy:
    def test_conditions_from_the_file(self, study_dir: Path) -> None:
        study = pz.Study(study_dir / "study.yaml")
        assert study.labels == ["No polymer", "Polymer"]
        assert study.root == study_dir
        assert [r.index for r in study["Polymer"].replicates] == [1, 2]
        assert study["Polymer"].equilibration == "0.25ns"
        assert study.settings("rg_first") == {"selection": "index 0"}

    def test_settings_and_labels_need_no_trajectories(
        self, study_dir: Path, tmp_path: Path
    ) -> None:
        shutil.rmtree(tmp_path / "data")
        study = pz.Study(study_dir)
        assert study.labels == ["No polymer", "Polymer"]
        assert study.settings("rg") == {"selection": "all"}
        with pytest.raises(ProtocolError, match="no run directory"):
            study["Polymer"]

    def test_unknown_run(self, study_dir: Path) -> None:
        with pytest.raises(ProtocolError, match="no analysis run 'sasa'"):
            pz.Study(study_dir).settings("sasa")

    def test_from_configs_has_no_protocol(self, study_dir: Path) -> None:
        study = pz.Study.from_configs(
            {"A": study_dir / "conditions/no_polymer/config.yaml"}, equilibration="0ns"
        )
        with pytest.raises(ProtocolError, match="built from config paths"):
            study.results("rg")


class TestAnalyzeStudy:
    def test_runs_and_stores_beside_the_study(self, study_dir: Path) -> None:
        result = _analyze("rg", "--study", str(study_dir))
        assert result.exit_code == 0, result.output
        means = {
            item.label: item.mean
            for item in read_results(study_dir / "results" / "rg").report.conditions
        }
        assert means == pytest.approx({"No polymer": 1.21, "Polymer": 2.21})
        assert (study_dir / "results" / "rg" / "report.json").is_file()

    def test_same_values_as_explicit_options(self, study_dir: Path, tmp_path: Path) -> None:
        assert _analyze("rg_first", "--study", str(study_dir)).exit_code == 0
        explicit = _analyze(
            "rg",
            "-c", str(study_dir / "conditions/no_polymer/config.yaml"),
            "-c", str(study_dir / "conditions/polymer/config.yaml"),
            "--label", "No polymer", "--label", "Polymer",
            "--eq", "0.25ns", "--set", "selection=index 0",
            "--output-dir", str(tmp_path / "explicit"),
            "--format", "json", "-o", str(tmp_path / "explicit.json"),
        )  # fmt: skip
        assert explicit.exit_code == 0, explicit.output
        stored = read_results(study_dir / "results" / "rg_first").report
        direct = json.loads((tmp_path / "explicit.json").read_text())
        assert [c.mean for c in stored.conditions] == [c["mean"] for c in direct["conditions"]]

    def test_two_runs_of_one_analysis_keep_separate_results(self, study_dir: Path) -> None:
        assert _analyze("rg", "--study", str(study_dir)).exit_code == 0
        assert _analyze("rg_first", "--study", str(study_dir)).exit_code == 0
        whole = read_results(study_dir / "results" / "rg").table
        first = read_results(study_dir / "results" / "rg_first").table
        assert not np.allclose(whole["value"], first["value"])

    def test_command_line_overrides_the_file(self, study_dir: Path) -> None:
        result = _analyze("rg", "--study", str(study_dir), "--eq", "0ns", "--replicates", "1")
        assert result.exit_code == 0, result.output
        report = read_results(study_dir / "results" / "rg").report
        assert report.equilibration == "0ns"
        assert [c.replicates for c in report.conditions] == [[1], [1]]
        assert report.conditions[0].mean == pytest.approx(1.0 + 0.1 + 0.045)

    def test_refuses_configs_with_study(self, study_dir: Path) -> None:
        result = _analyze("rg", "--study", str(study_dir), "-c", "x.yaml")
        assert result.exit_code == 2
        assert "-c and --label cannot be given" in result.output

    def test_unlisted_shipped_analysis_runs_with_a_note(self, study_dir: Path) -> None:
        protocol = load_study_file(study_dir)
        assert "rmsd" not in protocol.analyses
        result = _analyze("rmsd", "--study", str(study_dir), "--set", "selection=all")
        assert "does not list rmsd" in result.output

    def test_submit_dry_run_reports_with_the_study(self, study_dir: Path) -> None:
        result = CliRunner().invoke(
            cli, ["analyze", "rg_first", "--study", str(study_dir), "--dry-run"]
        )
        assert result.exit_code == 0, result.output
        folder = next((study_dir / "results" / "rg_first" / "slurm").iterdir())
        report = (folder / "report.sbatch").read_text()
        tasks = (folder / "replicates.sbatch").read_text()
        assert "analyze rg_first --study" in report and " -c " not in report
        assert "analyze rg -c" in tasks and "selection=" in tasks
        assert len((folder / "tasks.tsv").read_text().splitlines()) == 4


class TestResultsWithoutTrajectories:
    def test_moved_study_reads_results(self, study_dir: Path, tmp_path: Path) -> None:
        assert _analyze("rg", "--study", str(study_dir)).exit_code == 0
        moved = shutil.copytree(study_dir, tmp_path / "elsewhere" / "my_study")
        shutil.rmtree(tmp_path / "data")
        results = pz.Study(moved).results("rg")
        assert results.names == ["rg"]
        assert set(results.table["condition"]) == {"No polymer", "Polymer"}
        assert len(results.table) == 2 * 2 * 7
        means = results.table.groupby(["condition", "replicate"])["value"].mean()
        assert means[("Polymer", 2)] == pytest.approx(2.0 + 0.2 + 0.06)
        assert results.report.conditions[1].mean == pytest.approx(2.21)

    def test_no_results_yet(self, study_dir: Path) -> None:
        with pytest.raises(ProtocolError, match="No stored results"):
            pz.Study(study_dir).results("rg")


class TestCheck:
    def test_reports_runs_results_and_citation(self, study_dir: Path) -> None:
        assert _analyze("rg", "--study", str(study_dir)).exit_code == 0
        result = CliRunner().invoke(cli, ["study", "check", str(study_dir)])
        assert result.exit_code == 0, result.output
        assert "control No polymer: runs [1, 2]" in result.output
        assert "analysis rg: selection=all; stored results" in result.output
        assert "analysis rg as rg_first: selection=index 0; no stored results" in result.output
        assert "cite: " in result.output and "PolyzyMD" in result.output

    def test_missing_trajectories_are_not_errors(self, study_dir: Path, tmp_path: Path) -> None:
        shutil.rmtree(tmp_path / "data")
        result = CliRunner().invoke(cli, ["study", "check", str(study_dir)])
        assert result.exit_code == 0, result.output
        assert "no runs found" in result.output

    def test_bad_file_exits_2(self, tmp_path: Path) -> None:
        _write(
            tmp_path / "study.yaml",
            "equilibration: 1ns\nconditions: {A: a.yaml}\nanalyses: {rg: {selction: x}}\n",
        )
        result = CliRunner().invoke(cli, ["study", "check", str(tmp_path)])
        assert result.exit_code == 2
        assert "Did you mean 'selection'" in result.output

    def test_unreadable_config_exits_2(self, tmp_path: Path) -> None:
        _write(tmp_path / "study.yaml", "equilibration: 1ns\nconditions: {A: missing.yaml}\n")
        result = CliRunner().invoke(cli, ["study", "check", str(tmp_path)])
        assert result.exit_code == 2
        assert "cannot read" in result.output


def test_citation_matches_citation_cff() -> None:
    import yaml

    from polyzymd import __version__
    from polyzymd.citation import REPOSITORY, TITLE

    cff = yaml.safe_load((Path(__file__).resolve().parents[2] / "CITATION.cff").read_text())
    assert cff["title"] == TITLE
    assert cff["repository-code"] == REPOSITORY
    assert str(cff["version"]) == __version__
    assert cff["authors"][0]["family-names"] == "Laforet"
