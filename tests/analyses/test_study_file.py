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
from tests._support.analysis_testkit import (
    write_committed_study,
    write_openmm_replicate,
    write_simulation_config,
)

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

    def test_reads_a_condition_config_outside_the_study(self, tmp_path: Path) -> None:
        """Analysis reads a config outside the study; only freeze refuses it."""
        (tmp_path / "study").mkdir()
        path = _write(
            tmp_path / "study" / "study.yaml",
            "equilibration: 1ns\nconditions: {Ext: ../external/cfg/config.yaml}\n",
        )
        config = load_study_file(path).conditions["Ext"]
        assert config == (tmp_path / "external" / "cfg" / "config.yaml").resolve()

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
        assert "-c cannot be given" in result.output

    def test_label_picks_conditions(self, study_dir: Path) -> None:
        import json

        result = _analyze("rg", "--study", str(study_dir), "--label", "Polymer", "--format", "json")
        assert result.exit_code == 0, result.output
        report = json.loads(result.stdout[result.stdout.index("{") :])
        assert [c["label"] for c in report["conditions"]] == ["Polymer"]
        # A run of some conditions never replaces the run's report.
        assert not (study_dir / "results" / "rg" / "report.json").exists()
        assert _analyze("rg", "--study", str(study_dir), "--label", "Nope").exit_code == 2

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
        assert "analyze rg_first --study" in tasks and '--label "$label"' in tasks
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
        assert "analysis rg: selection=all; window eq 0.25ns; stored results" in result.output
        assert (
            "analysis rg as rg_first: selection=index 0; window eq 0.25ns; no stored results"
            in result.output
        )
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


USER_MODULE = """
def _scale(value):
    return value * SCALE


SCALE = 1.0


def mean_rg(atoms, frames, offset=0.0):
    u = atoms.universe
    total = 0.0
    for _ in u.trajectory[frames]:
        total += atoms.radius_of_gyration()
    return _scale(total / len(frames)) + offset


def frame_rg(atoms):
    return _scale(atoms.radius_of_gyration())


def per_atom(atoms, frames):
    return list(atoms.indices + 1), [1.0] * len(atoms)
"""


@pytest.fixture()
def user_study(study_dir: Path) -> Path:
    (study_dir / "analyses").mkdir()
    (study_dir / "analyses" / "metrics.py").write_text(USER_MODULE)
    with (study_dir / "study.yaml").open("a") as handle:
        handle.write(
            "  my_rg:\n"
            "    function: analyses/metrics.py:mean_rg\n"
            "    kind: per_replicate\n"
            "    unit: A\n"
            "    selections: {atoms: all}\n"
            "    settings: {offset: 10.0}\n"
            "  my_series:\n"
            "    function: analyses/metrics.py:frame_rg\n"
            "    kind: timeseries\n"
            "    unit: A\n"
            "    selections: {atoms: all}\n"
            "  my_atoms:\n"
            "    function: analyses/metrics.py:per_atom\n"
            "    kind: per_replicate\n"
            "    labels: returned\n"
            "    selections: {atoms: all}\n"
        )
    return study_dir


class TestUserFunctions:
    def test_schema(self, user_study: Path) -> None:
        entry = load_study_file(user_study).analyses["my_rg"]
        assert entry.analysis is None
        assert entry.function.file == user_study / "analyses" / "metrics.py"
        assert entry.function.qualname == "mean_rg"
        assert entry.function.selections == {"atoms": "all"}

    @pytest.mark.parametrize(
        ("entry", "message"),
        [
            ("{function: analyses/nope.py:f, kind: per_replicate}", "no Python file"),
            ("{function: analyses/metrics.py, kind: per_replicate}", "file.py:function_name"),
            ("{function: analyses/metrics.py:mean_rg, kind: frames}", "kind must be"),
            ("{function: analyses/metrics.py:mean_rg, kind: per_replicate, unti: A}", "'unit'"),
        ],
    )
    def test_schema_refusals(self, user_study: Path, entry: str, message: str) -> None:
        with (user_study / "study.yaml").open("a") as handle:
            handle.write(f"  bad: {entry}\n")
        with pytest.raises(ProtocolError) as caught:
            load_study_file(user_study)
        assert message in str(caught.value) + str(caught.value.hint)

    def test_per_replicate(self, user_study: Path) -> None:
        result = _analyze("my_rg", "--study", str(user_study))
        assert result.exit_code == 0, result.output
        means = [c.mean for c in read_results(user_study / "results" / "my_rg").report.conditions]
        assert means == pytest.approx([11.21, 12.21])

    def test_timeseries(self, user_study: Path) -> None:
        result = _analyze("my_series", "--study", str(user_study))
        assert result.exit_code == 0, result.output
        results = read_results(user_study / "results" / "my_series")
        assert [c.mean for c in results.report.conditions] == pytest.approx([1.21, 2.21])
        assert len(results.table) == 2 * 2 * 7

    def test_returned_labels(self, user_study: Path) -> None:
        result = _analyze("my_atoms", "--study", str(user_study))
        assert result.exit_code == 0, result.output
        table = read_results(user_study / "results" / "my_atoms").table
        assert sorted(set(table["label"])) == [1, 2, 3, 4]

    def test_editing_a_helper_recomputes(self, user_study: Path) -> None:
        assert _analyze("my_rg", "--study", str(user_study)).exit_code == 0
        module = user_study / "analyses" / "metrics.py"
        module.write_text(module.read_text().replace("SCALE = 1.0", "SCALE = 2.0"))
        assert _analyze("my_rg", "--study", str(user_study)).exit_code == 0
        means = [c.mean for c in read_results(user_study / "results" / "my_rg").report.conditions]
        assert means == pytest.approx([2 * 1.21 + 10, 2 * 2.21 + 10])
        record = json.loads(
            next(
                (user_study / "results" / "my_rg").glob("polyzymd_results/*/*/*/record.json")
            ).read_text()
        )
        assert record["function"]["hash_of"] == "module_folder"
        assert record["function"]["module"] == "polyzymd_study.metrics"

    def test_check_imports_the_function(self, user_study: Path) -> None:
        result = CliRunner().invoke(cli, ["study", "check", str(user_study)])
        assert result.exit_code == 0, result.output
        assert "analysis my_rg (analyses/metrics.py:mean_rg, per_replicate)" in result.output
        (user_study / "analyses" / "metrics.py").write_text("raise RuntimeError('broken')\n")
        result = CliRunner().invoke(cli, ["study", "check", str(user_study)])
        assert result.exit_code == 2
        assert "RuntimeError: broken" in result.output

    def test_every_run(self, user_study: Path) -> None:
        result = CliRunner().invoke(
            cli, ["analyze", "--study", str(user_study), "--no-eq-check", "--no-plots"]
        )
        assert result.exit_code == 0, result.output
        for run in ("rg", "rg_first", "my_rg", "my_series", "my_atoms"):
            assert f"== {run}" in result.output
            assert (user_study / "results" / run / "report.json").is_file()

    def test_submit_tasks_run_the_study_function(self, user_study: Path) -> None:
        result = CliRunner().invoke(
            cli, ["analyze", "my_rg", "--study", str(user_study), "--dry-run"]
        )
        assert result.exit_code == 0, result.output
        folder = next((user_study / "results" / "my_rg" / "slurm").iterdir())
        assert "analyze my_rg --study" in (folder / "replicates.sbatch").read_text()


def test_no_name_needs_a_study() -> None:
    result = CliRunner().invoke(cli, ["analyze"])
    assert result.exit_code == 2
    assert "needs NAME" in result.output


class TestAnalysisWindow:
    """Issue 1: an analysis entry sets its own equilibration and until."""

    def _with_full(self, study_dir: Path) -> Path:
        text = (study_dir / "study.yaml").read_text()
        _write(
            study_dir / "study.yaml",
            text + "  rg_full:\n    analysis: rg\n    selection: all\n    equilibration: 0ns\n"
            "  rg_early:\n    analysis: rg\n    selection: all\n    until: 0.5ns\n",
        )
        return study_dir

    def test_entries_keep_their_own_window(self, study_dir: Path) -> None:
        import json

        from polyzymd.analyses.study_freeze import stale_runs

        root = self._with_full(study_dir)
        protocol = load_study_file(root)
        assert protocol.window("rg") == ("0.25ns", None)
        assert protocol.window("rg_full") == ("0ns", None)
        assert protocol.window("rg_early") == ("0.25ns", "0.5ns")
        for run in ("rg", "rg_full", "rg_early"):
            assert _analyze(run, "--study", str(root)).exit_code == 0
        study = pz.Study(root)
        frames = {
            run: len(study.results(run).table.query("condition == 'Polymer' and replicate == 1"))
            for run in ("rg", "rg_full", "rg_early")
        }
        assert frames == {"rg": 7, "rg_full": 10, "rg_early": 3}
        record = next((root / "results" / "rg_full").rglob("record.json"))
        assert json.loads(record.read_text())["equilibration"] == "0ns"
        stale = stale_runs(protocol)
        assert not {"rg", "rg_full", "rg_early"} & set(stale), stale
        check = CliRunner().invoke(cli, ["study", "check", str(root)])
        assert "window eq 0ns (its own)" in check.output
        assert "window eq 0.25ns until 0.5ns (its own)" in check.output

    def test_command_line_window_wins(self, study_dir: Path) -> None:
        root = self._with_full(study_dir)
        assert _analyze("rg_full", "--study", str(root), "--eq", "0.5ns").exit_code == 0
        table = pz.Study(root).results("rg_full").table
        assert len(table.query("condition == 'Polymer' and replicate == 1")) == 5

    def test_bad_window_is_refused(self, study_dir: Path) -> None:
        text = (study_dir / "study.yaml").read_text()
        _write(study_dir / "study.yaml", text + "  rg_bad: {analysis: rg, equilibration: soon}\n")
        with pytest.raises(ProtocolError, match="analyses.rg_bad: cannot read equilibration"):
            load_study_file(study_dir)


class TestCommonUntil:
    """Issue 3: until: common ends every replicate at the shortest one's last time."""

    def _unequal(self, tmp_path: Path, entry_until: str) -> Path:
        root = tmp_path / "s"
        for folder in ("a", "b"):
            config = write_simulation_config(
                root / "conditions" / folder, scratch=tmp_path / "data" / folder
            )
            for replicate, n_frames in ((1, 10), (2, 6 if folder == "b" else 8)):
                write_openmm_replicate(config, replicate, [1.0 + 0.01 * k for k in range(n_frames)])
        _write(
            root / "study.yaml",
            "equilibration: 0ns\n"
            "conditions: {A: conditions/a/config.yaml, B: conditions/b/config.yaml}\n"
            f"analyses:\n  rg: {{selection: all, until: {entry_until}}}\n",
        )
        return root

    def test_python_api_aligns_replicates(self, tmp_path: Path) -> None:
        root = self._unequal(tmp_path, "common")
        protocol = load_study_file(root)
        study = pz.Study.from_configs(
            dict(protocol.conditions), equilibration="0ns", until="common"
        )
        times = [r.times.tolist() for c in study for r in c.replicates]
        assert all(t == times[0] for t in times) and len(times[0]) == 6
        assert study["A"].until == "0.5ns"

    def test_command_line_records_the_common_end(self, tmp_path: Path) -> None:
        import json

        from polyzymd.analyses.study_freeze import stale_runs

        root = self._unequal(tmp_path, "common")
        assert _analyze("rg", "--study", str(root)).exit_code == 0
        table = pz.Study(root).results("rg").table
        assert set(table.groupby(["condition", "replicate"]).size()) == {6}
        record = next((root / "results" / "rg").rglob("record.json"))
        assert json.loads(record.read_text())["until_ns"] == pytest.approx(0.5)
        assert "rg" not in stale_runs(load_study_file(root))


def test_an_analysis_sets_its_own_stride(study_dir: Path) -> None:
    from polyzymd.analyses.study_freeze import stale_runs

    text = (study_dir / "study.yaml").read_text()
    _write(
        study_dir / "study.yaml", text + "  rg_sparse: {analysis: rg, selection: all, stride: 2}\n"
    )
    protocol = load_study_file(study_dir)
    assert protocol.stride_of("rg_sparse") == 2 and protocol.stride_of("rg") == 1
    assert _analyze("rg_sparse", "--study", str(study_dir)).exit_code == 0
    table = pz.Study(study_dir).results("rg_sparse").table
    assert len(table.query("condition == 'Polymer' and replicate == 1")) == 4
    assert "rg_sparse" not in stale_runs(protocol)
    check = CliRunner().invoke(cli, ["study", "check", str(study_dir)])
    assert "window eq 0.25ns stride 2 (its own)" in check.output


@pytest.mark.usefixtures("git_identity")
def test_relative_files_of_shipped_analyses_follow_the_study(tmp_path: Path) -> None:
    """reference_file: structures/ref.pdb is relative to study.yaml, not the shell."""
    root = write_committed_study(
        tmp_path, "  rmsd: {selection: all, reference_file: structures/ref.pdb}\n"
    )
    (root / "structures").mkdir(exist_ok=True)
    (root / "structures" / "ref.pdb").write_text("END\n")
    entry = load_study_file(root).analyses["rmsd"]
    assert Path(entry.settings["reference_file"]) == (root / "structures" / "ref.pdb").resolve()
