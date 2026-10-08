"""study init: where create_study records the runs of each condition."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

from polyzymd.analyses.study_scaffold import create_study
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore"), pytest.mark.usefixtures("git_identity")]


def test_a_relative_scratch_is_relative_to_its_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """study init records where the runs are, not '.'."""
    folder = tmp_path / "water"
    config = write_simulation_config(folder, scratch=Path("."))
    (folder / "test.pdb").write_text("END\n")
    monkeypatch.chdir(folder)  # where a user runs it, so the runs land beside the config
    write_openmm_replicate(config, 1, [1.0, 1.1, 1.2])
    monkeypatch.chdir(tmp_path)
    root = tmp_path / "st"
    create_study(root, conditions={"Water": config}, equilibration="0ns")
    recorded = yaml.safe_load((root / "data.local.yaml").read_text())["Water"]
    assert Path(recorded) == folder.resolve()


def _output(config: Path) -> tuple[Path, object]:
    """Return where ``config`` writes its runs, resolved, and its scratch_directory."""
    output = yaml.safe_load(config.read_text())["output"]
    return (config.parent / output["projects_directory"]).resolve(), output.get("scratch_directory")


def test_a_new_condition_of_a_project_runs_into_the_project_runs_folder(tmp_path: Path) -> None:
    """add-condition --new writes a config whose runs go into <project>/runs/<study>/<condition>."""
    from polyzymd.analyses.project_scaffold import create_project
    from polyzymd.analyses.study_scaffold import add_condition

    project = create_project(tmp_path / "paper", ["lipa363"], git=False).root
    config = add_condition(project / "lipa363", "No polymer", new=True)

    assert _output(config) == ((project / "runs" / "lipa363" / "no_polymer").resolve(), None)
    text = config.read_text()
    assert "scratch_directory" in text and "a lot of disk space" in text
    assert "runs/" in (project / ".gitignore").read_text().splitlines()
    assert ".polymer_cache/" in (project / ".gitignore").read_text().splitlines()
    assert not (config.parent / "job_scripts").exists()


def test_a_condition_of_a_lone_study_runs_into_the_study_runs_folder(tmp_path: Path) -> None:
    """In a study with no project, the runs go into <study>/runs/<condition>, which git ignores."""
    created = create_study(tmp_path / "st", new_conditions=["Water"], git=False)

    assert _output(created.conditions["Water"]) == (
        (created.root / "runs" / "water").resolve(),
        None,
    )
    assert "runs/" in (created.root / ".gitignore").read_text().splitlines()
    assert ".polymer_cache/" in (created.root / ".gitignore").read_text().splitlines()


def test_from_copies_a_sibling_condition_with_its_inputs(tmp_path: Path) -> None:
    """add-condition --from copies the other condition's config and the files it names."""
    from polyzymd.analyses.study_file import load_study_file
    from polyzymd.analyses.study_scaffold import add_condition

    created = create_study(tmp_path / "st", new_conditions=["Water"], git=False)
    first = created.conditions["Water"]
    (first.parent / "structures" / "protein_X.pdb").write_text("REMARK protein\nEND\n")

    second = add_condition(created.root, "Water 350 K", source="Water")

    assert second == created.root / "conditions" / "water_350_k" / "config.yaml"
    copied = yaml.safe_load(second.read_text())
    assert copied["enzyme"]["pdb_path"] == "structures/protein_X.pdb"
    assert (second.parent / "structures" / "protein_X.pdb").read_text() == "REMARK protein\nEND\n"
    assert _output(second) == ((created.root / "runs" / "water_350_k").resolve(), None)
    assert list(load_study_file(created.root).conditions) == ["Water", "Water 350 K"]
    assert "# Enzyme Configuration (REQUIRED)" in second.read_text()  # comments are kept


def test_from_refuses_a_condition_the_study_does_not_have(tmp_path: Path) -> None:
    """--from names an existing condition; another label is refused and nothing is written."""
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study_scaffold import add_condition

    created = create_study(tmp_path / "st", new_conditions=["Water"], git=False)
    with pytest.raises(ProtocolError, match="no condition 'Polymer'"):
        add_condition(created.root, "Copy", source="Polymer")
    assert not (created.root / "conditions" / "copy").exists()


def test_a_copied_config_without_runs_records_no_data_location(tmp_path: Path) -> None:
    """--config of a config that has not run: no data.local.yaml, and its runs go into runs/."""
    from polyzymd.analyses.study_scaffold import add_condition

    config = write_simulation_config(tmp_path / "example", scratch=Path("."))
    (config.parent / "test.pdb").write_text("END\n")
    created = create_study(tmp_path / "st", git=False)

    copy = add_condition(created.root, "Water", config=config)

    assert not (created.root / "data.local.yaml").exists()
    assert _output(copy) == ((created.root / "runs" / "water").resolve(), None)


def test_add_condition_git_ignores_runs_in_a_study_made_before_runs_existed(
    tmp_path: Path,
) -> None:
    """A .gitignore without runs/ gets it when a condition whose runs go there is added."""
    from polyzymd.analyses.study_scaffold import add_condition

    created = create_study(tmp_path / "st", git=False)
    gitignore = created.root / ".gitignore"
    gitignore.write_text("data.local.yaml\n")

    add_condition(created.root, "Water", new=True)

    assert gitignore.read_text().splitlines() == ["data.local.yaml", "runs/"]


@pytest.mark.parametrize("label", ["runs", "Logs", "slurm_logs", "results", "deposit", "figures"])
def test_a_condition_may_not_take_a_folder_name_polyzymd_uses(tmp_path: Path, label: str) -> None:
    """conditions/runs/ would be skipped as machine files; such labels are refused."""
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study_scaffold import add_condition

    created = create_study(tmp_path / "st", git=False)
    with pytest.raises(ProtocolError, match="reserved"):
        add_condition(created.root, label, new=True)
    with pytest.raises(ProtocolError, match="reserved"):
        create_study(tmp_path / "other", new_conditions=[label], git=False)
    assert not (tmp_path / "other").exists()


def test_a_copy_drops_the_machine_path_comment_of_a_deposited_config(tmp_path: Path) -> None:
    """The output lines of a copy name runs/, so the deposit's 'machine path removed' note goes."""
    from polyzymd.analyses.study_freeze import without_machine_paths
    from polyzymd.analyses.study_scaffold import add_condition

    created = create_study(tmp_path / "st", new_conditions=["Water"], git=False)
    first = created.conditions["Water"]
    (first.parent / "structures" / "protein_X.pdb").write_text("REMARK protein\nEND\n")
    absolute = re.sub(
        r"(?m)^(\s*projects_directory:).*$", rf"\1 {tmp_path / 'jobs'}", first.read_text()
    )
    first.write_text(without_machine_paths(absolute))
    assert "machine path removed" in first.read_text()

    copy = add_condition(created.root, "Water 350 K", source="Water")

    assert "machine path removed" not in copy.read_text()
    assert "# Enzyme Configuration (REQUIRED)" in copy.read_text()  # comments are kept
    assert _output(copy) == ((created.root / "runs" / "water_350_k").resolve(), None)
