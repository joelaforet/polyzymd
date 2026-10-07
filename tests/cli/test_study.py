"""polyzymd study locate and add-condition: where a study's runs are and go."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from polyzymd.analyses.study_file import load_study_file
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_committed_study

pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore"), pytest.mark.usefixtures("git_identity")]


def test_alike_runs_are_told_apart_by_folder_name_or_refused(tmp_path: Path) -> None:
    """Two conditions whose runs share a name are never mapped to one folder."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    download = tmp_path / "download"
    shutil.copytree(tmp_path / "scratch", download / "by_condition")
    result = CliRunner().invoke(cli, ["study", "locate", str(download), "--study", str(root)])
    assert result.exit_code == 0, result.output
    data = load_study_file(root).data
    assert data["Polymer"].name == "polymer" and data["No polymer"].name == "no_polymer"
    mixed = tmp_path / "mixed"
    shutil.copytree(tmp_path / "scratch" / "polymer", mixed / "a")
    shutil.copytree(tmp_path / "scratch" / "no_polymer", mixed / "b")
    refused = CliRunner().invoke(cli, ["study", "locate", str(mixed), "--study", str(root)])
    assert refused.exit_code == 2 and "cannot be told apart" in refused.output


def _write_manifest(root: Path, folders: dict[str, Path]) -> None:
    """Write a manifest.json listing the size of every run file of each condition."""
    import json

    conditions = {
        label: {
            "replicates": {
                "1": {
                    "files": [
                        {"path": str(path.relative_to(folder)), "size": path.stat().st_size}
                        for path in sorted(folder.rglob("*"))
                        if path.is_file()
                    ]
                }
            }
        }
        for label, folder in folders.items()
    }
    (root / "manifest.json").write_text(json.dumps({"conditions": conditions}))


def test_runs_of_equal_size_go_to_the_folder_named_for_each_condition(tmp_path: Path) -> None:
    """Folders named for the conditions decide before file sizes; one folder never serves two."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    scratch = tmp_path / "scratch"
    _write_manifest(root, {"No polymer": scratch / "no_polymer", "Polymer": scratch / "polymer"})
    download = tmp_path / "download"
    shutil.copytree(scratch, download)
    result = CliRunner().invoke(cli, ["study", "locate", str(download), "--study", str(root)])
    assert result.exit_code == 0, result.output
    data = load_study_file(root).data
    assert data["Polymer"].name == "polymer" and data["No polymer"].name == "no_polymer"
    tied = tmp_path / "tied"
    shutil.copytree(scratch / "no_polymer", tied / "a")
    shutil.copytree(scratch / "polymer", tied / "b")
    refused = CliRunner().invoke(cli, ["study", "locate", str(tied), "--study", str(root)])
    assert refused.exit_code == 2, refused.output
    data = load_study_file(root).data
    assert data["Polymer"] != data["No polymer"]


def test_conditions_with_runs_named_apart_can_share_one_folder(tmp_path: Path) -> None:
    """One scratch folder holds the runs of two enzymes; locate writes it for both."""
    import yaml

    from polyzymd.analyses.study_scaffold import create_study
    from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

    configs = {}
    for label, enzyme in (("Wild type", "WT"), ("Mutant", "MUT")):
        folder = tmp_path / "runs" / enzyme
        config = write_simulation_config(folder, scratch=tmp_path / "scratch")
        data = yaml.safe_load(config.read_text())
        data["enzyme"]["name"] = enzyme
        config.write_text(yaml.safe_dump(data, sort_keys=False))
        (folder / "test.pdb").write_text("REMARK input\nEND\n")
        for replicate in (1, 2):
            write_openmm_replicate(config, replicate, [1.0 + 0.01 * k for k in range(10)])
        configs[label] = config
    root = tmp_path / "study"
    create_study(root, conditions=configs, equilibration="0.25ns")
    (root / "data.local.yaml").unlink(missing_ok=True)
    scratch = tmp_path / "scratch"
    result = CliRunner().invoke(cli, ["study", "locate", str(scratch), "--study", str(root)])
    assert result.exit_code == 0, result.output
    data = load_study_file(root).data
    assert data["Wild type"] == data["Mutant"] == scratch.resolve()


def test_locate_keeps_the_data_file_when_nothing_is_located(tmp_path: Path) -> None:
    """Entries written by hand survive a locate that finds no condition, and nothing is written."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    text = "# mine\nNo polymer: ../scratch/no_polymer\nPolymer: ../scratch/polymer\n"
    (root / "data.local.yaml").write_text(text)
    empty = tmp_path / "empty"
    empty.mkdir()
    result = CliRunner().invoke(cli, ["study", "locate", str(empty), "--study", str(root)])
    assert result.exit_code == 2
    assert "wrote" not in result.output
    assert (root / "data.local.yaml").read_text() == text


def test_add_condition_names_where_its_runs_are(tmp_path: Path) -> None:
    """add-condition records and prints the folder holding the config's runs, if it has any."""
    from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    config = write_simulation_config(tmp_path / "sims" / "new", scratch=tmp_path / "old_runs")
    (config.parent / "test.pdb").write_text("REMARK input\nEND\n")
    write_openmm_replicate(config, 1, [1.0, 1.1, 1.2])
    result = CliRunner().invoke(
        cli, ["study", "add-condition", "New", "--config", str(config), "--study", str(root)]
    )
    assert result.exit_code == 0, result.output
    where = (tmp_path / "old_runs").resolve()
    assert f"data New: {where} (from the config's scratch_directory)" in result.output
    help_text = CliRunner().invoke(cli, ["study", "add-condition", "--help"]).output
    assert "data.local.yaml" in help_text


def test_add_condition_warns_where_new_runs_go_and_their_disk_space(tmp_path: Path) -> None:
    """Each way of adding a condition prints where the runs go and how to send them to scratch."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    for label, how in (("Draft", ["--new"]), ("Draft 2", ["--from", "Draft"])):
        result = CliRunner().invoke(
            cli, ["study", "add-condition", label, *how, "--study", str(root)]
        )
        assert result.exit_code == 0, result.output
        warning = next(line for line in result.output.splitlines() if line.startswith("warning:"))
        assert "runs" in warning and "a lot of disk space" in warning
        assert "scratch_directory" in warning and "cluster" in warning


def test_the_new_condition_template_validates_once_its_pdb_is_set(tmp_path: Path) -> None:
    """The config add-condition --new writes sets the OpenMM platform and validates as written."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    result = CliRunner().invoke(
        cli, ["study", "add-condition", "Draft", "--new", "--study", str(root)]
    )
    assert result.exit_code == 0, result.output
    config = root / "conditions" / "draft" / "config.yaml"
    text = config.read_text(encoding="utf-8")
    assert "\nopenmm:\n  platform:" in text
    assert "{{" not in text and "{%" not in text
    assert text.count("PolyzyMD: Created by Joseph R. Laforet Jr.") == 1
    pdb = Path(__file__).resolve().parents[2] / "examples" / "quickstart" / "trpcage.pdb"
    config.write_text(text.replace("structures/protein_X.pdb", str(pdb)), encoding="utf-8")

    result = CliRunner().invoke(cli, ["validate", "-c", str(config)])

    assert result.exit_code == 0, result.output
    assert "Configuration is valid!" in result.output
    assert "Referenced file warnings" not in result.output


def test_check_names_the_control_of_each_stratum(tmp_path: Path) -> None:
    from tests.cli.test_analyze import _temperature_polymer_study

    root = _temperature_polymer_study(tmp_path, "comparison: {within: temperature_K}")
    output = CliRunner().invoke(cli, ["study", "check", str(root)]).output
    for kelvin in (300, 330, 360):
        assert f"control none_{kelvin}: replicates" in output
        assert f"condition SBMA_{kelvin}: replicates" in output
