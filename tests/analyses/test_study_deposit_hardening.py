# ruff: noqa: F811 - the tests take the study fixture imported from test_study_freeze
"""add-condition, check's production lengths, and a deposit free of machine paths and self-describing.

Uses the analysed, committed study of test_study_freeze.py.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from polyzymd.analyses.study_file import load_study_file
from polyzymd.analyses.study_freeze import freeze, without_machine_paths
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config
from tests.analyses.test_study_freeze import git_identity, study  # noqa: F401 - fixtures

pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
    pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed"),
]
SCHEMAS = Path(__file__).resolve().parents[2] / "src" / "polyzymd" / "analyses" / "schemas"


class TestAddCondition:
    def test_copies_and_lists_a_config(self, study: Path, tmp_path: Path) -> None:
        source = write_simulation_config(tmp_path / "extra", scratch=tmp_path / "extra_data")
        (source.parent / "test.pdb").write_text("REMARK\nEND\n")
        write_openmm_replicate(source, 1, [3.0 + 0.01 * k for k in range(10)])
        before = (study / "study.yaml").read_text()
        result = CliRunner().invoke(
            cli,
            ["study", "add-condition", "SBMA 100%", "--config", str(source), "--study", str(study)],
        )
        assert result.exit_code == 0, result.output
        protocol = load_study_file(study)
        assert list(protocol.conditions)[-1] == "SBMA 100%"
        assert protocol.conditions["SBMA 100%"] == study / "conditions" / "sbma_100" / "config.yaml"
        after = (study / "study.yaml").read_text()
        assert after.startswith(before.split("conditions:")[0])  # the rest of the file is kept
        assert "metadata:" in after

    def test_new_condition(self, study: Path) -> None:
        result = CliRunner().invoke(
            cli, ["study", "add-condition", "Draft", "--new", "--study", str(study)]
        )
        assert result.exit_code == 0, result.output
        assert (study / "conditions" / "draft" / "config.yaml").is_file()
        assert "Draft" in load_study_file(study).conditions

    def test_refusals(self, study: Path, tmp_path: Path) -> None:
        both = CliRunner().invoke(cli, ["study", "add-condition", "X", "--study", str(study)])
        assert both.exit_code == 2 and "either a config to copy or new" in both.output
        taken = CliRunner().invoke(
            cli, ["study", "add-condition", "Polymer", "--new", "--study", str(study)]
        )
        assert taken.exit_code == 2 and "already has a condition" in taken.output


class TestCheck:
    def test_production_lengths(self, study: Path) -> None:
        result = CliRunner().invoke(cli, ["study", "check", str(study)])
        assert "; production " not in result.output
        result = CliRunner().invoke(cli, ["study", "check", str(study), "--production"])
        assert "production 0.9 ns" in result.output

    def test_frozen_copy_says_how_to_reproduce(self, study: Path, tmp_path: Path) -> None:
        deposit = freeze(study).deposit
        copy = shutil.copytree(deposit / "study", tmp_path / "copy")
        result = CliRunner().invoke(cli, ["study", "check", str(copy)])
        assert "reproduce: this study was frozen" in result.output


class TestDeposit:
    def test_no_machine_paths(self, study: Path, tmp_path: Path) -> None:
        deposit = freeze(study).deposit
        for path in (deposit / "study").rglob("*"):
            if path.is_file() and path.suffix in {".json", ".yaml", ".md", ".cff"}:
                assert str(tmp_path) not in path.read_text(), path
        for name in ("manifest.json", "README.md", "CITATION.cff", "UPLOAD.md"):
            assert str(tmp_path) not in (deposit / name).read_text(), name
        config = yaml.safe_load(
            (deposit / "study" / "conditions" / "polymer" / "config.yaml").read_text()
        )
        assert config["output"]["projects_directory"] == "."
        assert config["output"]["scratch_directory"] == "data"

    def test_records_and_reports_use_relative_paths(self, study: Path) -> None:
        record = json.loads(
            next((study / "results" / "rg").glob("polyzymd_results/*/*/*/record.json")).read_text()
        )
        assert not Path(record["topology"]["path"]).is_absolute()
        assert all(not Path(t["path"]).is_absolute() for t in record["trajectories"])
        report = json.loads((study / "results" / "rg" / "report.json").read_text())
        assert all(not Path(v).is_absolute() for v in report["provenance"]["output_paths"].values())

    def test_without_machine_paths(self) -> None:
        text = "# Copied by polyzymd study init from /home/x/run/config.yaml\noutput:\n  projects_directory: /home/x\n  scratch_directory: /scratch/x\n  naming_template: r{replicate}\n"
        clean = without_machine_paths(text)
        assert "/home" not in clean and "/scratch" not in clean and "r{replicate}" in clean
        assert clean.startswith("# Copied by polyzymd from config.yaml\n")

    def test_readme_from_metadata(self, study: Path) -> None:
        readme = (freeze(study).deposit / "README.md").read_text()
        assert readme.startswith("# A test study") and "To test study freeze." in readme
        assert "## How to cite" in readme and "PolyzyMD:" in readme and "## Results" in readme
        assert "**rg:**" in readme and "CC-BY-4.0" in readme and "MIT" in readme

    def test_upload_lists_both_licences_and_the_schema(self, study: Path) -> None:
        result = freeze(study)
        assert "CC-BY-4.0 (data, results and figures) and MIT (code" in result.guide.read_text()
        assert (result.upload / "manifest-1.schema.json").is_file()

    def test_simulated_with_and_its_warning(self, study: Path) -> None:
        result = freeze(study)
        replicate = result.manifest["conditions"]["Polymer"]["replicates"]["1"]
        assert "simulated_with" in replicate
        assert any("record no OpenMM version" in w for w in result.warnings)


class TestSchemas:
    def test_schemas_are_json_with_ids(self) -> None:
        for name in ("study-1.schema.json", "manifest-1.schema.json"):
            schema = json.loads((SCHEMAS / name).read_text())
            assert schema["$id"].endswith(name)

    def test_manifest_and_study_validate(self, study: Path) -> None:
        jsonschema = pytest.importorskip("jsonschema")
        manifest = freeze(study).manifest
        jsonschema.Draft202012Validator(
            json.loads((SCHEMAS / "manifest-1.schema.json").read_text())
        ).validate(json.loads(json.dumps(manifest)))
        jsonschema.Draft202012Validator(
            json.loads((SCHEMAS / "study-1.schema.json").read_text())
        ).validate(yaml.safe_load((study / "study.yaml").read_text()))
