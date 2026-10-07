"""Study folders: study init, data.local.yaml, --data, study locate and git provenance.

The conditions are synthetic OpenMM runs (four unit-mass atoms on a cross),
as in test_study_file.py: frame k of replicate r of a condition with offset
o has a radius of gyration of o + 0.1 r + 0.01 k.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.results import read_results
from polyzymd.analyses.study_file import load_study_file
from polyzymd.analyses.study_git import git_state
from polyzymd.analyses.study_scaffold import condition_folder, create_study
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

GIT = shutil.which("git") is not None


@pytest.fixture()
def sources(tmp_path: Path) -> dict[str, Path]:
    """Two configs outside any study, each with two replicates of run data."""
    configs = {}
    for label, offset in (("No polymer", 1.0), ("Polymer", 2.0)):
        folder = tmp_path / "runs" / condition_folder(label)
        config = write_simulation_config(folder, scratch=tmp_path / "scratch" / folder.name)
        (folder / "test.pdb").write_text("REMARK input structure of the condition\nEND\n")
        for replicate in (1, 2):
            write_openmm_replicate(
                config, replicate, [offset + 0.1 * replicate + 0.01 * k for k in range(10)]
            )
        configs[label] = config
    return configs


@pytest.fixture(autouse=True)
def git_identity(monkeypatch: pytest.MonkeyPatch) -> None:
    """Let git commit in a test without the user's own identity."""
    for key, value in {
        "GIT_AUTHOR_NAME": "Test",
        "GIT_AUTHOR_EMAIL": "test@example.com",
        "GIT_COMMITTER_NAME": "Test",
        "GIT_COMMITTER_EMAIL": "test@example.com",
    }.items():
        monkeypatch.setenv(key, value)


def _study(tmp_path: Path, sources: dict[str, Path], **kwargs) -> Path:
    root = tmp_path / "my_study"
    create_study(root, conditions=sources, equilibration="0.25ns", **kwargs)
    text = (
        (root / "study.yaml")
        .read_text()
        .replace("analyses: {}", "analyses:\n  rg: {selection: all}")
    )
    (root / "study.yaml").write_text(text)
    return root


def _analyze(*arguments: str):
    return CliRunner().invoke(cli, ["analyze", *arguments, "--no-eq-check", "--no-plots"])


class TestInit:
    def test_layout(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = _study(tmp_path, sources, git=False)
        for name in (
            "study.yaml",
            "README.md",
            "data.example.yaml",
            ".gitignore",
            "LICENSE-data",
            "LICENSE-code",
            "environment/README.md",
            "conditions/no_polymer/config.yaml",
            "conditions/polymer/config.yaml",
        ):
            assert (root / name).is_file(), name
        for folder in ("structures", "analyses", "figures", "results"):
            assert (root / folder).is_dir()
        assert "data.local.yaml" in (root / ".gitignore").read_text()
        assert "CC-BY-4.0" in (root / "LICENSE-data").read_text()
        assert "MIT License" in (root / "LICENSE-code").read_text()
        assert "How to cite" in (root / "README.md").read_text()

    def test_conditions_are_listed_and_load(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        protocol = load_study_file(_study(tmp_path, sources, git=False))
        assert list(protocol.conditions) == ["No polymer", "Polymer"]
        assert protocol.equilibration == "0.25ns"
        study = pz.Study(protocol.path)
        assert [r.index for r in study["Polymer"].replicates] == [1, 2]

    def test_input_files_are_copied_with_relative_paths(
        self, tmp_path: Path, sources: dict[str, Path]
    ) -> None:
        root = _study(tmp_path, sources, git=False)
        text = (root / "conditions/polymer/config.yaml").read_text()
        assert "pdb_path: structures/" in text
        copied = list((root / "conditions/polymer/structures").iterdir())
        assert len(copied) == 1 and copied[0].suffix == ".pdb"
        # The copy must still load after the original inputs are gone.
        shutil.rmtree(tmp_path / "runs")
        assert pz.Study(root).labels == ["No polymer", "Polymer"]

    def test_refuses_an_existing_study(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = _study(tmp_path, sources, git=False)
        with pytest.raises(ProtocolError, match="already holds"):
            create_study(root)

    def test_refuses_clashing_folders(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        with pytest.raises(ProtocolError, match="own label and folder"):
            create_study(
                tmp_path / "s", conditions={"A b": sources["Polymer"], "a_b": sources["Polymer"]}
            )

    def test_new_condition_runs_polyzymd_init(self, tmp_path: Path) -> None:
        result = CliRunner().invoke(
            cli, ["study", "init", str(tmp_path / "s"), "--new-condition", "No polymer", "--no-git"]
        )
        assert result.exit_code == 0, result.output
        assert (tmp_path / "s/conditions/no_polymer/config.yaml").is_file()
        assert (tmp_path / "s/conditions/no_polymer/structures").is_dir()
        assert "set equilibration" in result.output

    def test_cli_with_conditions(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        result = CliRunner().invoke(
            cli,
            [
                "study", "init", str(tmp_path / "s"),
                "--condition", f"No polymer={sources['No polymer']}",
                "--equilibration", "1ns", "--holder", "Ada Lovelace", "--no-git",
            ],
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        assert "1 input files copied" in result.output
        assert "Ada Lovelace" in (tmp_path / "s/LICENSE-code").read_text()
        bad = CliRunner().invoke(
            cli, ["study", "init", str(tmp_path / "t"), "--condition", "nolabel"]
        )
        assert bad.exit_code == 2

    def test_missing_input_is_left_with_a_warning(
        self, tmp_path: Path, sources: dict[str, Path]
    ) -> None:
        (sources["Polymer"].parent / "test.pdb").unlink()
        result = CliRunner().invoke(
            cli,
            [
                "study",
                "init",
                str(tmp_path / "s"),
                "--condition",
                f"P={sources['Polymer']}",
                "--no-git",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "does not exist, so it was left as it was" in result.output

    @pytest.mark.skipif(not GIT, reason="git is not installed")
    def test_git_repository_with_one_commit(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = tmp_path / "my_study"
        created = create_study(root, conditions=sources, equilibration="0.25ns")
        assert created.commit
        state = git_state(root)
        assert state["commit"] == created.commit
        assert state["uncommitted"] == []


class TestData:
    def test_data_local_yaml_replaces_the_config_scratch(
        self, tmp_path: Path, sources: dict[str, Path]
    ) -> None:
        root = _study(tmp_path, sources, git=False)
        before = {c.label: c.config_hash for c in pz.Study(root)}
        shutil.move(tmp_path / "scratch", tmp_path / "moved")
        with pytest.raises(ProtocolError, match="no run directory"):
            pz.Study(root)["Polymer"]
        (root / "data.local.yaml").write_text(
            f"No polymer: {tmp_path / 'moved' / 'no_polymer'}\nPolymer: ../moved/polymer\n"
        )
        study = pz.Study(root)
        assert [r.index for r in study["Polymer"].replicates] == [1, 2]
        assert {c.label: c.config_hash for c in study} == before

    def test_unknown_label_in_data_file(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = _study(tmp_path, sources, git=False)
        (root / "data.local.yaml").write_text("Polymr: /x\n")
        with pytest.raises(ProtocolError, match="unknown key 'Polymr'") as caught:
            load_study_file(root)
        assert "Polymer" in caught.value.hint

    def test_data_option(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = _study(tmp_path, sources, git=False)
        shutil.move(tmp_path / "scratch" / "polymer", tmp_path / "elsewhere")
        import json

        result = _analyze(
            "rg",
            "--study",
            str(root),
            "--label",
            "Polymer",
            "--data",
            str(tmp_path / "elsewhere"),
            "--format",
            "json",
        )
        assert result.exit_code == 0, result.output
        report = json.loads(result.stdout[result.stdout.index("{") :])
        assert report["conditions"][0]["mean"] == pytest.approx(2.21)

    def test_locate_writes_data_local_yaml(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = _study(tmp_path, sources, git=False)
        download = tmp_path / "download" / "zenodo_123"
        download.mkdir(parents=True)
        for folder in ("no_polymer", "polymer"):
            shutil.move(tmp_path / "scratch" / folder, download / folder)
        result = CliRunner().invoke(
            cli, ["study", "locate", str(tmp_path / "download"), "--study", str(root)]
        )
        assert result.exit_code == 0, result.output
        assert "Polymer: replicates [1, 2]" in result.output
        assert load_study_file(root).data["Polymer"] == download / "polymer"
        check = CliRunner().invoke(cli, ["study", "check", str(root)])
        assert "(from data.local.yaml)" in check.output
        assert _analyze("rg", "--study", str(root)).exit_code == 0

    def test_locate_reports_what_it_cannot_find(
        self, tmp_path: Path, sources: dict[str, Path]
    ) -> None:
        root = _study(tmp_path, sources, git=False)
        (tmp_path / "empty").mkdir()
        result = CliRunner().invoke(
            cli, ["study", "locate", str(tmp_path / "empty"), "--study", str(root)]
        )
        assert result.exit_code == 2
        assert "no run directories named" in result.output

    def test_data_example_is_not_read(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = _study(tmp_path, sources, git=False)
        # study init records where the copied conditions' runs are; without
        # that file, data.example.yaml is still never read.
        (root / "data.local.yaml").unlink()
        assert load_study_file(root).data == {}


@pytest.mark.skipif(not GIT, reason="git is not installed")
class TestGitProvenance:
    def test_report_records_commit_and_uncommitted_inputs(
        self, tmp_path: Path, sources: dict[str, Path]
    ) -> None:
        root = tmp_path / "my_study"
        created = create_study(root, conditions=sources, equilibration="0.25ns")
        result = _analyze("rg", "--study", str(root), "--set", "selection=all")
        assert result.exit_code == 0, result.output
        study = json.loads((root / "results/rg/report.json").read_text())["provenance"]["study"]
        assert study["git"]["commit"] == created.commit
        assert study["git"]["inputs_uncommitted"] == []
        assert study["run"] == "rg" and len(study["sha256"]) == 64

        (root / "analyses" / "new.py").write_text("x = 1\n")
        result = _analyze("rg", "--study", str(root), "--set", "selection=all")
        assert "1 uncommitted input files (analyses/new.py)" in result.output
        study = json.loads((root / "results/rg/report.json").read_text())["provenance"]["study"]
        assert "analyses/new.py" in study["git"]["inputs_uncommitted"]

    def test_results_are_outputs_not_inputs(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = tmp_path / "my_study"
        create_study(root, conditions=sources, equilibration="0.25ns")
        (root / "results" / "x.json").write_text("{}")
        (root / "data.local.yaml").write_text("{}\n")
        state = git_state(root)
        assert state["inputs_uncommitted"] == []
        assert "results/x.json" in state["uncommitted"]

    def test_check_prints_the_commit(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = tmp_path / "my_study"
        create_study(root, conditions=sources, equilibration="0.25ns")
        result = CliRunner().invoke(cli, ["study", "check", str(root)])
        assert "git: commit" in result.output and "inputs committed" in result.output

    def test_outside_a_repository(self, tmp_path: Path, sources: dict[str, Path]) -> None:
        root = _study(tmp_path, sources, git=False)
        if (
            subprocess.run(["git", "-C", str(root), "rev-parse"], capture_output=True).returncode
            == 0
        ):
            pytest.skip("the test folder is inside a git repository")
        assert git_state(root) is None
        result = CliRunner().invoke(cli, ["study", "check", str(root)])
        assert "git: not a repository" in result.output
