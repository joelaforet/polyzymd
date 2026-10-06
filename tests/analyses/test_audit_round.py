"""Regressions for the 1.3 audit round, wave A (analysis correctness and reproducibility).

Each test names its finding in the audit log
(``PAPERS/polyzymd_v1.3_refactor_handoff/audit_2026-10-06/AUDIT_LOG.md``) and
is built from the reviewer's reproduction.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.study_file import load_study_file
from polyzymd.analyses.study_freeze import freeze, stale_runs, without_machine_paths
from polyzymd.analyses.study_scaffold import condition_folder, create_study
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config
from tests.analyses.test_project import project  # noqa: F401  (fixture)

pytest.importorskip("MDAnalysis")
pytestmark = pytest.mark.filterwarnings("ignore")


@pytest.fixture(autouse=True)
def git_identity(monkeypatch: pytest.MonkeyPatch) -> None:
    for key, value in {
        "GIT_AUTHOR_NAME": "Test",
        "GIT_AUTHOR_EMAIL": "test@example.com",
        "GIT_COMMITTER_NAME": "Test",
        "GIT_COMMITTER_EMAIL": "test@example.com",
    }.items():
        monkeypatch.setenv(key, value)


def _git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", "-C", str(root), *arguments], capture_output=True, text=True
    ).stdout


def _study(tmp_path: Path, analyses: str, name: str = "my_study") -> Path:
    """A committed study of two conditions with three replicates of ten frames each."""
    configs = {}
    for label, offset in (("No polymer", 1.0), ("Polymer", 2.0)):
        folder = tmp_path / "runs" / condition_folder(label)
        config = write_simulation_config(folder, scratch=tmp_path / "scratch" / folder.name)
        (folder / "test.pdb").write_text("REMARK input\nEND\n")
        for replicate in (1, 2, 3):
            write_openmm_replicate(
                config, replicate, [offset + 0.1 * replicate + 0.01 * k for k in range(10)]
            )
        configs[label] = config
    root = tmp_path / name
    create_study(root, conditions=configs, equilibration="0.25ns")
    text = (root / "study.yaml").read_text().replace("analyses: {}", "analyses:\n" + analyses)
    (root / "study.yaml").write_text(text)
    _git(root, "commit", "-qam", "analyses")
    return root


def _analyze(*arguments: str, cwd: Path | None = None):
    old = Path.cwd()
    if cwd is not None:
        os.chdir(cwd)
    try:
        return CliRunner().invoke(cli, ["analyze", *arguments], catch_exceptions=False)
    finally:
        os.chdir(old)


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
    _git(root, "commit", "-qam", "analyses")
    return root


class TestStoredValues:
    def test_reordered_parts_recompute(self, tmp_path: Path) -> None:
        """ADV-1: the order of parts names the stored columns, so it keys reuse."""
        root = _study(
            tmp_path,
            "  pair: {function: analyses/two.py:two, kind: timeseries, universe: u, parts: [a, b]}\n",
        )
        (root / "analyses" / "two.py").write_text("def two(u):\n    return {'a': 1.0, 'b': 100.0}\n")
        options = ["--study", str(root), "--no-plots", "--no-eq-check", "--run", "a"]
        first = _json(_analyze("pair", *options, "--format", "json"))
        study_yaml = root / "study.yaml"
        study_yaml.write_text(study_yaml.read_text().replace("parts: [a, b]", "parts: [b, a]"))
        second = _json(_analyze("pair", *options, "--format", "json"))
        assert [c["mean"] for c in first["conditions"]] == [1.0, 1.0]
        assert [c["mean"] for c in second["conditions"]] == [1.0, 1.0]
        means = pz.Study(root).results("pair").table.groupby("part")["value"].mean()
        assert means["a"] == 1.0 and means["b"] == 100.0

    def test_replicate_table_keeps_each_quantity(self, tmp_path: Path) -> None:
        """ADV-2: two hydrogen-bond summaries in one run are two quantities, not two replicates."""
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

    def test_replicate_table_fills_missing_labels_as_the_report(self, tmp_path: Path) -> None:
        """ADV-3: a pair one replicate never formed counts as 0 in the table, as in the report."""
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

    def test_shipped_analyses_hash_their_modules(self, monkeypatch) -> None:
        """ADV-4: a fix in a module a shipped analysis imports changes its record."""
        from polyzymd.analyses import functions, timeseries

        record = timeseries._function_record(functions.rms_decomposition)
        assert record["hash_of"] == "polyzymd_modules"
        timeseries._shipped_code_hash.cache_clear()
        reference = Path(__import__("polyzymd.analyses.reference").analyses.reference.__file__)
        real = Path.read_bytes

        def edited(self):
            data = real(self)
            return data + b"# fixed\n" if self == reference else data

        monkeypatch.setattr(Path, "read_bytes", edited)
        try:
            assert timeseries._function_record(functions.rms_decomposition)["hash"] != record["hash"]
        finally:
            monkeypatch.undo()
            timeseries._shipped_code_hash.cache_clear()


class TestSettings:
    def test_a_reference_file_is_the_reference(self, tmp_path: Path) -> None:
        """ADV-5: rmsd and rmsf with a reference_file and no mode measure from that file."""
        from polyzymd.analyses.protocols import FUNCTION_ANALYSES
        from polyzymd.analyses.reference import reference

        for name in ("rmsd", "rmsf", "rmsd_per_residue"):
            assert FUNCTION_ANALYSES[name]["reference_mode"] is None, name
        with pytest.raises(ProtocolError, match="does not use the file"):
            reference("centroid", "all", file=__file__)

    def test_relative_files_of_shipped_analyses_follow_the_study(self, tmp_path: Path) -> None:
        """ADV-6: reference_file: structures/ref.pdb is relative to study.yaml, not the shell."""
        root = _study(tmp_path, "  rmsd: {selection: all, reference_file: structures/ref.pdb}\n")
        (root / "structures").mkdir(exist_ok=True)
        (root / "structures" / "ref.pdb").write_text("END\n")
        entry = load_study_file(root).analyses["rmsd"]
        assert Path(entry.settings["reference_file"]) == (root / "structures" / "ref.pdb").resolve()

    def test_set_values_round_trip(self) -> None:
        """ADV-16: a float such as 1e-05 stays a float through --set."""
        from polyzymd.cli.analyze import _set_value, _settings

        values = {"a": 1e-05, "b": "protein and name CA", "c": ["EGM", "SBM"], "d": None, "e": "x: y"}
        text = tuple(f"{key}={_set_value(value)}" for key, value in values.items())
        assert _settings(text) == values

    def test_repeated_pair_labels_are_refused(self, tmp_path: Path) -> None:
        """ADV-17: two pairs with one label would overwrite each other's values."""
        from polyzymd.analyses.protocols import _analyze_pairs

        pair = {"label": "d", "selection_a": "name C1", "selection_b": "name C2"}
        with pytest.raises(ProtocolError, match="repeat"):
            _analyze_pairs(
                "distances",
                None,
                {"pairs": [pair, pair]},
                None,
                recompute=False,
                output_dir=None,
                eq_check=False,
                plots=False,
            )


class TestFreezeAndLayout:
    def test_freeze_needs_a_git_identity_before_writing(self, tmp_path: Path, monkeypatch) -> None:
        """ADV-14: without a name and email nothing is written that names a tag."""
        root = _study(tmp_path, "  rg: {selection: all}\n")
        for key in ("GIT_AUTHOR_NAME", "GIT_AUTHOR_EMAIL", "GIT_COMMITTER_NAME", "GIT_COMMITTER_EMAIL"):
            monkeypatch.delenv(key, raising=False)
        monkeypatch.setenv("HOME", str(tmp_path / "home"))
        monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
        with pytest.raises(ProtocolError, match="no user name and email"):
            freeze(root)
        assert not (root / "manifest.json").exists()

    def test_a_flow_style_output_loses_its_machine_paths(self) -> None:
        """ADV-15: flow style and block scalars are rewritten through YAML."""
        flow = "output: {projects_directory: /home/u/p, scratch_directory: /scratch/u/r}\nx: 1\n"
        block = "output:\n  projects_directory: >-\n    /home/u/p\n  scratch_directory: /s/u\n"
        for text in (flow, block):
            output = yaml.safe_load(without_machine_paths(text))["output"]
            assert output == {"projects_directory": ".", "scratch_directory": "data"}

    def test_the_freeze_warning_names_no_machine_path(self, tmp_path: Path) -> None:
        """ADV-12: the manifest is published, so its warnings hold no scratch path."""
        root = _study(tmp_path, "  rg: {selection: all}\n")
        shutil.rmtree(tmp_path / "scratch")
        result = freeze(root)
        assert not any(str(tmp_path) in warning for warning in result.warnings)

    def test_a_subset_report_is_not_saved_and_the_report_job_gets_the_labels(
        self, tmp_path: Path
    ) -> None:
        """ADV-11: --label runs some conditions; their report never replaces the run's."""
        import shlex

        root = _study(tmp_path, "  rg: {selection: all}\n")
        result = _analyze("rg", "--study", str(root), "--label", "Polymer", "--no-plots")
        assert result.exit_code == 0 and "not saved" in result.output
        assert not (root / "results" / "rg" / "report.json").exists()
        dry = _analyze("rg", "--study", str(root), "--label", "Polymer", "--dry-run")
        assert dry.exit_code == 0, dry.output
        (report_job,) = (root / "results").rglob("report.sbatch")
        words = shlex.split(report_job.read_text().splitlines()[-1])
        assert words[words.index("--label") + 1] == "Polymer"

    def test_a_relative_scratch_is_relative_to_its_config(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """SIM-10: study init records where the runs are, not '.'."""
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


class TestProjects:
    def test_a_study_of_a_project_is_frozen_with_the_project(self, tmp_path: Path) -> None:
        """ADV-7: study freeze inside a project would leave out project.yaml and analyses/."""
        root = _study(tmp_path, "  rg: {selection: all}\n", name="lipa")
        (tmp_path / "project.yaml").write_text("studies: {lipa: lipa}\n")
        with pytest.raises(ProtocolError, match="is a study of the project") as info:
            freeze(root)
        assert "polyzymd project freeze" in info.value.hint

    @pytest.mark.parametrize("where", ["nested/lipa", "../outside/lipa"])
    def test_a_study_must_sit_directly_in_the_project(self, tmp_path: Path, where: str) -> None:
        """ADV-10: a study elsewhere never finds its paper, or is not published."""
        from polyzymd.analyses.project_file import load_project_file

        paper = tmp_path / "Paper"
        paper.mkdir()
        folder = (paper / where).resolve()
        folder.mkdir(parents=True)
        (folder / "study.yaml").write_text("conditions: {}\n")
        (paper / "project.yaml").write_text(f"studies: {{lipa: {where}}}\n")
        with pytest.raises(ProtocolError, match="directly inside the project"):
            load_project_file(paper)

    def test_an_analysis_entry_must_be_a_mapping(self, tmp_path: Path) -> None:
        """ADV-18: rg: notamapping gets a message, not a traceback."""
        from polyzymd.analyses.project_file import load_project_file

        paper = tmp_path / "Paper"
        (paper / "lipa").mkdir(parents=True)
        (paper / "lipa" / "study.yaml").write_text("conditions: {}\n")
        (paper / "project.yaml").write_text("studies: {lipa: lipa}\nanalyses: {rg: notamapping}\n")
        with pytest.raises(ProtocolError, match="must be a mapping"):
            load_project_file(paper)


def test_project_study_manifests_list_what_git_tracks(project: Path) -> None:  # noqa: F811
    """ADV-9: a study's manifest in a project lists the files the deposit holds."""
    from polyzymd.analyses.project_freeze import freeze_project
    from polyzymd.analyses.study_git import init_repository

    (project / ".gitignore").write_text("__pycache__/\ndata.local.yaml\n")
    init_repository(project, "start")
    (project / "lipa" / "__pycache__").mkdir()
    (project / "lipa" / "__pycache__" / "x.pyc").write_bytes(b"x")
    freeze_project(project)
    manifest = json.loads((project / "lipa" / "manifest.json").read_text())
    assert not any("__pycache__" in name for name in manifest["files"])
    assert manifest["git"]["commit"]


class TestLocate:
    def test_alike_runs_are_told_apart_by_folder_name_or_refused(self, tmp_path: Path) -> None:
        """NOV-3: two conditions whose runs share a name are never mapped to one folder."""
        root = _study(tmp_path, "  rg: {selection: all}\n")
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


def test_a_project_logs_once_in_its_own_folder(tmp_path: Path) -> None:
    """ADV-19 and WF-4: one log for the command, in the project's logs/."""
    from tests.analyses.test_project import _study as project_study

    paper = tmp_path / "Paper"
    project_study(paper, tmp_path / "data", "lipa", "{core: name C1 C2}")
    (paper / "project.yaml").write_text("studies: {lipa: lipa}\nanalyses: {rg: {selection: all}}\n")
    result = _analyze("--project", str(paper), "rg")
    assert result.exit_code == 0, result.output
    assert result.output.count("log: ") == 1
    assert list((paper / "logs").glob("polyzymd-analyze-*.log"))
    assert not (paper / "lipa" / "logs").exists() or not list((paper / "lipa" / "logs").iterdir())


def test_compiled_python_and_job_files_are_never_inputs() -> None:
    """D1 refuses uncommitted inputs, so bytecode and job logs must not count as inputs."""
    from polyzymd.analyses.study_git import is_output

    for path in ("lipa/analyses/__pycache__/helper.cpython-311.pyc", "x.pyc", "results/a/slurm/j.sh"):
        assert is_output(path), path
    for path in ("lipa/analyses/helper.py", "study.yaml", "stats/plan.py"):
        assert not is_output(path), path
