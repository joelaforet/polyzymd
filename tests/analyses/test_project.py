"""Projects: one study per protein, shared analyses, statistics, project init and freeze."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.study_file import load_study_file
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore::UserWarning")]

LID = "def lid_size(lid, reference):\n    assert reference.endswith('ref.pdb')\n    return float(len(lid))\n"


def _study(root: Path, data: Path, name: str, regions: str, analyses: str = "{}") -> Path:
    """Write one protein's study: a control and one polymer condition with a factor."""
    folder = root / name
    for condition, offset in (("none", 1.0), ("half", 2.0)):
        config = write_simulation_config(
            folder / "conditions" / condition, scratch=data / name / condition
        )
        for replicate in (1, 2):
            write_openmm_replicate(
                config, replicate, [offset + 0.1 * replicate + 0.01 * k for k in range(6)]
            )
    (folder / "structures").mkdir()
    (folder / "structures" / "ref.pdb").write_text("END\n")
    (folder / "study.yaml").write_text(
        f"description: protein {name}\n"
        "equilibration: 0ns\n"
        "structures: {reference: structures/ref.pdb}\n"
        f"regions: {regions}\n"
        "conditions:\n"
        "  No polymer: conditions/none\n"
        "  Half: {config: conditions/half, factors: {sbma_fraction: 0.5}}\n"
        f"analyses: {analyses}\n"
    )
    return folder


@pytest.fixture()
def project(tmp_path: Path) -> Path:
    """Two proteins: lipa has a lid, rml has none but an analysis of its own."""
    root = tmp_path / "Paper"
    (root / "analyses").mkdir(parents=True)
    (root / "analyses" / "lid.py").write_text(LID)
    _study(root, tmp_path / "data", "lipa", "{core: name C1 C2, lid: name C3}")
    _study(
        root,
        tmp_path / "data",
        "rml",
        "{core: name C1 C3 C4}",
        "{rg_all: {analysis: rg, selection: all}}",
    )
    (root / "project.yaml").write_text(
        "studies: {lipa: lipa, rml: rml}\n"
        "analyses:\n"
        "  rg: {selection: region core}\n"
        "  lid:\n"
        "    function: analyses/lid.py:lid_size\n"
        "    kind: timeseries\n"
        "    studies: [lipa]\n"
        "    selections: {lid: region lid and name C3}\n"
        "    settings: {reference: structure reference}\n"
        "metadata: {title: A paper}\n"
    )
    return root


def _analyze(*arguments: str):
    return CliRunner().invoke(cli, ["analyze", *arguments, "--no-eq-check", "--no-plots"])


class TestStudiesOfAProject:
    def test_project_analyses_use_each_proteins_names(self, project: Path) -> None:
        lipa, rml = load_study_file(project / "lipa"), load_study_file(project / "rml")
        assert lipa.project_label == "lipa" and lipa.description == "protein lipa"
        assert lipa.analyses["rg"].settings == {"selection": "(name C1 C2)"}
        assert rml.analyses["rg"].settings == {"selection": "(name C1 C3 C4)"}
        lid = lipa.analyses["lid"].function
        assert lid.selections == {"lid": "(name C3) and name C3"}
        assert lid.settings["reference"] == str((project / "lipa/structures/ref.pdb").resolve())
        assert lid.file == (project / "analyses" / "lid.py").resolve()
        assert set(lipa.analyses) == {"rg", "lid"} and set(rml.analyses) == {"rg", "rg_all"}
        assert lipa.factors == {"Half": {"sbma_fraction": 0.5}}
        assert lipa.conditions["Half"] == (project / "lipa/conditions/half/config.yaml").resolve()
        assert lipa.metadata == {"title": "A paper"}

    def test_a_region_a_study_lacks_is_an_error(self, project: Path) -> None:
        text = (project / "project.yaml").read_text().replace("    studies: [lipa]\n", "")
        (project / "project.yaml").write_text(text)
        with pytest.raises(ProtocolError, match="in study rml.* uses region lid"):
            load_study_file(project / "rml")

    def test_a_study_cannot_redefine_a_project_analysis(self, project: Path) -> None:
        text = (project / "rml" / "study.yaml").read_text().replace("rg_all:", "rg:")
        (project / "rml" / "study.yaml").write_text(text)
        with pytest.raises(ProtocolError, match="also an analysis of"):
            load_study_file(project / "rml")

    def test_unknown_study_in_studies_list(self, project: Path) -> None:
        text = (project / "project.yaml").read_text().replace("[lipa]", "[lipa, calb]")
        (project / "project.yaml").write_text(text)
        with pytest.raises(ProtocolError, match="lists studies \\['calb'\\]"):
            pz.Project(project)


class TestProjectResults:
    def test_analyze_project_and_read_one_table(self, project: Path) -> None:
        result = _analyze("rg", "--project", str(project))
        assert result.exit_code == 0, result.output
        assert "== study lipa" in result.output and "== study rml" in result.output
        assert (project / "lipa" / "results" / "rg" / "report.json").is_file()
        paper = pz.Project(project)
        assert paper.runs_in("lid") == ["lipa"]
        stored = paper.results("rg")
        assert list(stored.table.columns[:2]) == ["study", "name"]
        assert set(stored.table["study"]) == {"lipa", "rml"}
        assert set(stored.reports) == {"lipa", "rml"}
        half = stored.table.query("condition == 'Half'")["sbma_fraction"]
        assert (half == 0.5).all()
        assert stored.table.query("condition == 'No polymer'")["sbma_fraction"].isna().all()

    def test_a_study_without_results_is_named(self, project: Path) -> None:
        assert _analyze("rg", "--study", str(project / "lipa")).exit_code == 0
        with pytest.raises(ProtocolError, match="studies rml have no stored results"):
            pz.Project(project).results("rg")

    def test_only_listed_studies_run_an_analysis(self, project: Path) -> None:
        result = _analyze("lid", "--project", str(project))
        assert result.exit_code == 0, result.output
        assert "== study rml" not in result.output
        table = pz.Project(project).results("lid").table
        assert set(table["study"]) == {"lipa"} and (table["value"] == 1.0).all()

    def test_project_check(self, project: Path) -> None:
        result = CliRunner().invoke(cli, ["project", "check", str(project)])
        assert result.exit_code == 0, result.output
        assert "analysis lid: studies lipa" in result.output
        assert "analysis rg: studies lipa, rml" in result.output
        assert "region lid: name C3" in result.output
        assert "factors sbma_fraction=0.5" in result.output


def test_files_validate_against_the_schemas(project: Path) -> None:
    import json

    import yaml

    jsonschema = pytest.importorskip("jsonschema")
    schemas = Path(pz.__file__).parent / "analyses" / "schemas"
    for name, files in (
        ("project-1.schema.json", [project / "project.yaml"]),
        ("study-1.schema.json", [project / "lipa" / "study.yaml", project / "rml" / "study.yaml"]),
    ):
        validator = jsonschema.Draft202012Validator(json.loads((schemas / name).read_text()))
        for file in files:
            validator.validate(yaml.safe_load(file.read_text()))


@pytest.fixture()
def graded(tmp_path: Path) -> Path:
    """A project of one study whose Rg rises by 2 per unit of sbma_fraction, with a stats plan."""
    root = tmp_path / "Graded"
    folder = root / "prot"
    levels = {"none": None, "q1": 0.25, "q3": 0.75, "full": 1.0}
    for name, level in levels.items():
        config = write_simulation_config(
            folder / "conditions" / name, scratch=tmp_path / "data" / name
        )
        base = 1.0 if level is None else 1.0 + 2.0 * level
        for replicate in (1, 2, 3):
            write_openmm_replicate(
                config, replicate, [base + 0.01 * replicate + 0.001 * k for k in range(5)]
            )
    conditions = "".join(
        f"  {name}: conditions/{name}\n"
        if level is None
        else f"  {name}: {{config: conditions/{name}, factors: {{sbma_fraction: {level}}}}}\n"
        for name, level in levels.items()
    )
    (folder / "study.yaml").write_text(f"equilibration: 0ns\nconditions:\n{conditions}")
    (root / "stats").mkdir(parents=True)
    (root / "stats" / "plan.py").write_text(
        "def plan(project):\n"
        "    table = project.replicate_table('rg')\n"
        "    means = table.groupby('condition', sort=False)['value'].mean().reset_index()\n"
        "    return {'condition_means': means, 'n_replicates': int(len(table))}\n"
    )
    (root / "project.yaml").write_text(
        "studies: {prot: prot}\n"
        "analyses: {rg: {selection: all}}\n"
        "stats: {plan: stats/plan.py:plan}\n"
    )
    return root


class TestStatistics:
    def test_replicate_table(self, graded: Path) -> None:
        assert _analyze("rg", "--project", str(graded)).exit_code == 0
        table = pz.Project(graded).replicate_table("rg")
        assert len(table) == 4 * 3
        assert list(table.columns[:3]) == ["study", "condition", "replicate"]
        row = table.query("condition == 'q3' and replicate == 2").iloc[0]
        assert row["value"] == pytest.approx(1.0 + 1.5 + 0.02 + 0.002, abs=1e-6)
        assert row["sbma_fraction"] == 0.75
        assert table.query("condition == 'none'")["sbma_fraction"].isna().all()

    def test_trend_over_a_numeric_factor(self, graded: Path) -> None:
        assert _analyze("rg", "--study", str(graded / "prot")).exit_code == 0
        report = pz.Study(graded / "prot").results("rg").report
        (trend,) = report.trends
        assert trend.factor == "sbma_fraction" and trend.conditions == ["q1", "q3", "full"]
        assert trend.n_replicates == 9
        assert trend.slope == pytest.approx(2.0, abs=1e-6) and trend.significant
        assert any(v.startswith("mean_rg rises with sbma_fraction") for v in report.verdict)

    def test_stats_plan_runs_and_tracks_its_inputs(self, graded: Path) -> None:
        import pandas as pd

        assert _analyze("rg", "--project", str(graded)).exit_code == 0
        result = CliRunner().invoke(cli, ["stats", str(graded)])
        assert result.exit_code == 0, result.output
        folder = graded / "results" / "stats" / "plan"
        means = pd.read_csv(folder / "condition_means.csv")
        assert list(means["condition"]) == ["none", "q1", "q3", "full"]
        assert '"n_replicates": 12' in (folder / "values.json").read_text()
        check = CliRunner().invoke(cli, ["project", "check", str(graded)])
        assert "stats plan: up to date" in check.output
        text = (graded / "project.yaml").read_text()
        (graded / "project.yaml").write_text(text + "metadata: {title: Filled in later}\n")
        check = CliRunner().invoke(cli, ["project", "check", str(graded)])
        assert "stats plan: up to date" in check.output
        with (graded / "stats" / "plan.py").open("a") as handle:
            handle.write("# edited\n")
        check = CliRunner().invoke(cli, ["project", "check", str(graded)])
        assert "stats plan: stale: its code changed" in check.output

    def test_stats_without_a_plan(self, project: Path) -> None:
        result = CliRunner().invoke(cli, ["stats", str(project)])
        assert result.exit_code == 2 and "has no stats: plan" in result.output


def test_a_moved_project_reuses_its_results(project: Path, tmp_path: Path) -> None:
    import shutil

    from polyzymd.analyses.study_freeze import stale_runs

    assert _analyze("lid", "--project", str(project)).exit_code == 0
    moved = Path(shutil.copytree(project, tmp_path / "elsewhere" / "Paper"))
    assert "lid" not in stale_runs(load_study_file(moved / "lipa"))
    record = next((moved / "lipa" / "results" / "lid").rglob("record.json")).read_text()
    assert str(project) not in record


class TestProjectInit:
    def test_new_project_scaffold(self, tmp_path: Path) -> None:
        result = CliRunner().invoke(
            cli,
            ["project", "init", str(tmp_path / "P"), "--study", "lipa363", "--study", "rml333"]
            + ["--no-git"],
        )
        assert result.exit_code == 0, result.output
        text = (tmp_path / "P" / "project.yaml").read_text()
        assert "studies:" in text and "  lipa363: lipa363\n  rml333: rml333" in text
        assert "title: TODO" in text and "polyzymd:" not in text
        study = (tmp_path / "P" / "lipa363" / "study.yaml").read_text()
        assert "description: TODO" in study and "regions: {}" in study
        assert not (tmp_path / "P" / "lipa363" / "LICENSE-code").exists()
        assert (tmp_path / "P" / "LICENSE-code").is_file()

    def test_conditions_are_added_to_a_new_study(self, tmp_path: Path) -> None:
        config = write_simulation_config(tmp_path / "runs" / "a", scratch=tmp_path / "data")
        write_openmm_replicate(config, 1, [1.0, 1.1, 1.2])
        assert (
            CliRunner()
            .invoke(cli, ["project", "init", str(tmp_path / "P"), "--study", "lipa", "--no-git"])
            .exit_code
            == 0
        )
        result = CliRunner().invoke(
            cli,
            ["study", "add-condition", "No polymer", "--config", str(config)]
            + ["--study", str(tmp_path / "P" / "lipa")],
        )
        assert result.exit_code == 0, result.output
        protocol = load_study_file(tmp_path / "P" / "lipa")
        assert list(protocol.conditions) == ["No polymer"]

    def test_labels_must_be_folder_names(self, tmp_path: Path) -> None:
        result = CliRunner().invoke(
            cli, ["project", "init", str(tmp_path / "P"), "--study", "LipA 363", "--no-git"]
        )
        assert result.exit_code == 2 and "'lipa_363'" in result.output


class TestProjectFreeze:
    def test_one_manifest_citation_tag_and_deposit(self, project: Path, monkeypatch) -> None:
        import json
        import subprocess

        from polyzymd.analyses.project_freeze import freeze_project
        from polyzymd.analyses.study_git import init_repository

        for key in ("GIT_AUTHOR_NAME", "GIT_COMMITTER_NAME"):
            monkeypatch.setenv(key, "Test")
        for key in ("GIT_AUTHOR_EMAIL", "GIT_COMMITTER_EMAIL"):
            monkeypatch.setenv(key, "test@example.com")
        assert _analyze("--project", str(project)).exit_code == 0
        init_repository(project, "start")
        result = freeze_project(project)
        assert result.tag == "project-v1" and result.commit
        manifest = json.loads((project / "manifest.json").read_text())
        assert set(manifest["studies"]) == {"lipa", "rml"}
        assert "lipa / Half" in manifest["conditions"]
        assert manifest["analyses"]["lid"] == ["lipa"]
        assert (project / "lipa" / "manifest.json").is_file()
        assert (project / "CITATION.cff").is_file()
        assert (result.deposit / "study" / "rml" / "study.yaml").is_file()
        readme = (result.deposit / "README.md").read_text()
        assert "### lipa" in readme and "protein rml" in readme
        log = subprocess.run(
            ["git", "-C", str(project), "tag", "--list"], capture_output=True, text=True
        ).stdout
        assert "project-v1" in log and "study-v" not in log


def test_names_may_start_with_a_digit_and_unknown_ones_are_errors(project: Path) -> None:
    from polyzymd.analyses.project_file import resolve_names

    structures = {"4TGL_open": Path("/x/4tgl.pdb")}
    assert resolve_names("structure 4TGL_open", {}, structures, "here") == "/x/4tgl.pdb"
    assert resolve_names("region 1st", {"1st": "resid 1"}, {}, "here") == "(resid 1)"
    with pytest.raises(ProtocolError, match="uses structure 3TGL"):
        resolve_names("structure 3TGL", {}, structures, "here")
    with pytest.raises(ProtocolError, match="one name"):
        resolve_names("structure a b", {}, structures, "here")


class TestAuditFindings:
    """Regressions for the pre-release reproducibility audit."""

    def test_editing_a_helper_module_recomputes(self, project: Path) -> None:
        (project / "analyses" / "helper.py").write_text("K = 1.0\n")
        (project / "analyses" / "lid.py").write_text(
            "from helper import K\n\n"
            "def lid_size(lid, reference):\n    return K * float(len(lid))\n"
        )
        assert _analyze("lid", "--project", str(project)).exit_code == 0
        assert set(pz.Project(project).results("lid").table["value"]) == {1.0}
        (project / "analyses" / "helper.py").write_text("K = 2.0\n")
        from polyzymd.analyses.study_freeze import stale_runs

        assert any("code" in w for w in stale_runs(load_study_file(project / "lipa"))["lid"])
        assert _analyze("lid", "--project", str(project)).exit_code == 0
        assert set(pz.Project(project).results("lid").table["value"]) == {2.0}

    def test_freeze_sees_factor_selection_file_and_replicate_changes(
        self, project: Path, tmp_path: Path
    ) -> None:
        from polyzymd.analyses.study_freeze import stale_runs

        assert _analyze("--project", str(project)).exit_code == 0
        lipa = project / "lipa"
        assert stale_runs(load_study_file(lipa)) == {}
        text = (lipa / "study.yaml").read_text()
        (lipa / "study.yaml").write_text(
            text.replace("sbma_fraction: 0.5", "sbma_fraction: 0.6").replace(
                "lid: name C3", "lid: name C3 C4"
            )
        )
        (lipa / "structures" / "ref.pdb").write_text("REMARK changed\nEND\n")
        config = load_study_file(lipa).conditions["Half"]
        write_openmm_replicate(config, 3, [2.0 + 0.01 * k for k in range(6)])
        why = " | ".join(" | ".join(v) for v in stale_runs(load_study_file(lipa)).values())
        assert "factors changed" in why
        assert "selections differ" in why
        assert "content of ref.pdb changed" in why
        assert "replicates [1, 2, 3] on disk" in why

    def test_reports_and_manifests_hold_no_machine_paths(self, project: Path, monkeypatch) -> None:
        from polyzymd.analyses.project_freeze import freeze_project
        from polyzymd.analyses.study_git import init_repository

        for key in ("GIT_AUTHOR_NAME", "GIT_COMMITTER_NAME"):
            monkeypatch.setenv(key, "Test")
        for key in ("GIT_AUTHOR_EMAIL", "GIT_COMMITTER_EMAIL"):
            monkeypatch.setenv(key, "test@example.com")
        assert _analyze("--project", str(project)).exit_code == 0
        report = (project / "lipa" / "results" / "lid" / "report.json").read_text()
        assert str(project) not in report and '"structures/ref.pdb"' in report
        init_repository(project, "start")
        result = freeze_project(project)
        for name in ("manifest.json", "lipa/manifest.json"):
            assert str(project) not in (project / name).read_text(), name
        assert "analyses/lid.py" in result.manifest["files"]

    def test_freeze_names_partial_reports_and_stale_stats(self, graded: Path, monkeypatch) -> None:
        import json

        from polyzymd.analyses.project_freeze import freeze_project
        from polyzymd.analyses.study_git import init_repository

        for key in ("GIT_AUTHOR_NAME", "GIT_COMMITTER_NAME"):
            monkeypatch.setenv(key, "Test")
        for key in ("GIT_AUTHOR_EMAIL", "GIT_COMMITTER_EMAIL"):
            monkeypatch.setenv(key, "test@example.com")
        assert _analyze("rg", "--project", str(graded)).exit_code == 0
        report = graded / "prot" / "results" / "rg" / "report.json"
        data = json.loads(report.read_text())
        data.update(status="partial", problems=["condition q1 is left out: boom"])
        report.write_text(json.dumps(data))
        init_repository(graded, "start")
        warnings = " | ".join(freeze_project(graded).warnings)
        assert "partial report: condition q1 is left out: boom" in warnings
        assert "stats plan plan is not run" in warnings
        assert "PARTIAL REPORT" in (graded / "deposit" / "README.md").read_text()

    def test_trend_refuses_non_finite_values_and_two_levels(self) -> None:
        from polyzymd.analyses.protocols import ConditionReport
        from polyzymd.analyses.statistics_plan import trend_sentence, trend_tests

        def report(values):
            class R:
                conditions = [
                    ConditionReport(
                        label=f"c{i}", n_replicates=len(v), mean=0.0, replicate_values=v
                    )
                    for i, v in enumerate(values)
                ]

            return R()

        factors = {"c0": {"x": 0.1}, "c1": {"x": 0.5}, "c2": {"x": 0.9}}
        (nan,) = trend_tests(report([[1.0, 2.0], [1.5, float("nan")], [2.0, 3.0]]), factors)
        assert not nan.testable and "not finite" in nan.reason
        assert "not finite" in trend_sentence("m", None, nan)
        # Tight replicates around means 1.0, 1.0, 1.5: three points that are not a line.
        steps = report([[0.99, 1.0, 1.01], [0.99, 1.0, 1.01], [1.49, 1.5, 1.51]])
        (step,) = trend_tests(steps, {"c0": {"x": 0.0}, "c1": {"x": 0.5}, "c2": {"x": 1.0}})
        assert step.testable and not step.significant and step.p > 0.2
        assert step.slope == pytest.approx(0.5)
        two = {"c0": {"x": 0.1}, "c1": {"x": 0.5}}
        (pair,) = trend_tests(report([[1.0, 1.1], [2.0, 2.1]]), two)
        assert not pair.testable and "at least three" in pair.reason

    def test_report_records_uncommitted_project_files(self, project: Path, monkeypatch) -> None:
        from polyzymd.analyses.study_git import init_repository

        for key in ("GIT_AUTHOR_NAME", "GIT_COMMITTER_NAME"):
            monkeypatch.setenv(key, "Test")
        for key in ("GIT_AUTHOR_EMAIL", "GIT_COMMITTER_EMAIL"):
            monkeypatch.setenv(key, "test@example.com")
        init_repository(project, "start")
        with (project / "analyses" / "lid.py").open("a") as handle:
            handle.write("# edited\n")
        assert _analyze("lid", "--project", str(project)).exit_code == 0
        import json

        study = json.loads((project / "lipa" / "results" / "lid" / "report.json").read_text())[
            "provenance"
        ]["study"]
        assert "analyses/lid.py" in study["git"]["inputs_uncommitted"]
        assert study["project"]["label"] == "lipa" and len(study["project"]["sha256"]) == 64


class TestSecondAuditFindings:
    """Regressions for the second pre-release reproducibility audit."""

    def test_helper_packages_and_data_files_are_part_of_the_code(self, project: Path) -> None:
        (project / "analyses" / "util").mkdir()
        (project / "analyses" / "util" / "__init__.py").write_text("")
        (project / "analyses" / "util" / "k.py").write_text("K = 1.0\n")
        (project / "analyses" / "scale.json").write_text('{"s": 1.0}')
        (project / "analyses" / "lid.py").write_text(
            "import json, pathlib\nfrom util.k import K\n\n"
            "def lid_size(lid, reference):\n"
            "    s = json.loads((pathlib.Path(__file__).parent / 'scale.json').read_text())['s']\n"
            "    return K * s * float(len(lid))\n"
        )
        assert _analyze("lid", "--project", str(project)).exit_code == 0
        (project / "analyses" / "util" / "k.py").write_text("K = 2.0\n")
        assert _analyze("lid", "--project", str(project)).exit_code == 0
        assert set(pz.Project(project).results("lid").table["value"]) == {2.0}
        (project / "analyses" / "scale.json").write_text('{"s": 3.0}')
        assert _analyze("lid", "--project", str(project)).exit_code == 0
        assert set(pz.Project(project).results("lid").table["value"]) == {6.0}

    def test_freeze_sees_added_conditions_and_entry_edits(self, project: Path) -> None:
        from polyzymd.analyses.study_freeze import stale_runs

        assert _analyze("--project", str(project)).exit_code == 0
        rml = project / "rml"
        config = write_simulation_config(rml / "conditions" / "full", scratch=project.parent / "x")
        text = (
            (rml / "study.yaml")
            .read_text()
            .replace(
                "  Half: {config: conditions/half, factors: {sbma_fraction: 0.5}}\n",
                "  Half: {config: conditions/half, factors: {sbma_fraction: 0.5}}\n"
                f"  Full: {config.parent.relative_to(rml)}\n",
            )
        )
        (rml / "study.yaml").write_text(text)
        project_text = (project / "project.yaml").read_text()
        (project / "project.yaml").write_text(
            project_text.replace(
                "    kind: timeseries\n", "    kind: timeseries\n    reduce: std\n"
            )
        )
        stale = stale_runs(load_study_file(rml))
        assert any("not in the report ['Full']" in w for w in stale["rg"])
        lipa = stale_runs(load_study_file(project / "lipa"))
        assert any("the entry's reduce changed" in w for w in lipa["lid"])

    def test_submit_job_files_stay_out_of_the_deposit(self, project: Path, monkeypatch) -> None:
        from polyzymd.analyses.project_freeze import freeze_project
        from polyzymd.analyses.study_git import init_repository

        for key in ("GIT_AUTHOR_NAME", "GIT_COMMITTER_NAME"):
            monkeypatch.setenv(key, "Test")
        for key in ("GIT_AUTHOR_EMAIL", "GIT_COMMITTER_EMAIL"):
            monkeypatch.setenv(key, "test@example.com")
        assert _analyze("--project", str(project)).exit_code == 0
        jobs = project / "lipa" / "results" / "rg" / "slurm" / "rg_1"
        jobs.mkdir(parents=True)
        (jobs / "tasks.tsv").write_text(f"{project}/lipa/conditions/none/config.yaml\tA\t1\n")
        init_repository(project, "start")
        result = freeze_project(project)
        deposited = "\n".join(
            p.read_text(errors="ignore")
            for p in result.deposit.rglob("*")
            if p.is_file() and p.suffix in (".json", ".tsv", ".md", ".yaml", ".cff")
        )
        assert str(project) not in deposited
        assert not list((result.deposit / "study").rglob("slurm"))


class TestThirdAuditFindings:
    """Regressions for the third pre-release reproducibility audit."""

    def test_missing_nan_is_not_stale_after_a_fresh_analyze(self, tmp_path: Path) -> None:
        from polyzymd.analyses.study_freeze import stale_runs
        from tests.analyses.test_study_hardening import TestMissingLabels

        root = TestMissingLabels()._study(tmp_path, ".nan")
        assert _analyze("byf", "--study", str(root)).exit_code == 0
        assert "byf" not in stale_runs(load_study_file(root))

    def test_replicate_table_keeps_a_nan_frame(self, project: Path) -> None:
        (project / "analyses" / "lid.py").write_text(
            "import math\n\ndef lid_size(lid, reference):\n"
            "    return math.nan if lid.universe.trajectory.frame == 2 else 1.0\n"
        )
        assert _analyze("lid", "--project", str(project)).exit_code == 0
        table = pz.Project(project).replicate_table("lid")
        assert table["value"].isna().all()

    def test_root_level_code_hashes_only_python(self, tmp_path: Path) -> None:
        from polyzymd.analyses.timeseries import code_files

        root = tmp_path / "study"
        (root / "figures").mkdir(parents=True)
        (root / ".pixi" / "envs").mkdir(parents=True)
        for name in ("study.yaml", "count.py", "big.dat", "figures/fig.png", ".pixi/envs/x.py"):
            (root / name).write_text("x")
        assert [p.name for p in code_files(root)] == ["count.py"]
        analyses = root / "analyses"
        (analyses / "util").mkdir(parents=True)
        for name in ("f.py", "table.json", "util/k.py"):
            (analyses / name).write_text("x")
        assert [p.relative_to(analyses).as_posix() for p in code_files(analyses)] == [
            "f.py",
            "table.json",
            "util/k.py",
        ]

    def test_purge_spares_packages_installed_under_the_folder(self, tmp_path: Path) -> None:
        import types

        from polyzymd.analyses.user_functions import _imported_from

        folder = tmp_path / "study"
        helper = types.ModuleType("helper")
        helper.__file__ = str(folder / "helper.py")
        package = types.ModuleType("util.k")
        package.__file__ = str(folder / "util" / "k.py")
        installed = types.ModuleType("numpy")
        installed.__file__ = str(
            folder / ".pixi" / "envs" / "site-packages" / "numpy" / "__init__.py"
        )
        source = types.ModuleType("polyzymd")
        source.__file__ = str(folder / "src" / "polyzymd" / "__init__.py")
        resolved = folder.resolve()
        assert _imported_from(helper, "helper", resolved)
        assert _imported_from(package, "util.k", resolved)
        assert not _imported_from(installed, "numpy", resolved)
        assert not _imported_from(source, "polyzymd", resolved)

    def test_machine_and_hidden_folders_stay_out_of_the_deposit(
        self, project: Path, monkeypatch
    ) -> None:
        from polyzymd.analyses.project_freeze import freeze_project

        assert _analyze("--project", str(project)).exit_code == 0
        logs = project / "lipa" / "conditions" / "half" / "slurm_logs"
        logs.mkdir(parents=True)
        (logs / "run.out").write_text(f"{project}/somewhere\n")
        (project / ".venv").mkdir()
        (project / ".venv" / "pyvenv.cfg").write_text(f"home = {project}\n")
        result = freeze_project(project)  # no git: the copy is made file by file
        shipped = [p.relative_to(result.deposit).as_posix() for p in result.deposit.rglob("*")]
        assert not any("slurm_logs" in p or ".venv" in p for p in shipped)
        assert (
            "conditions/half/slurm_logs/run.out"
            not in (project / "lipa" / "manifest.json").read_text()
        )

    def test_project_init_refuses_a_folder_with_files(self, tmp_path: Path) -> None:
        (tmp_path / "P").mkdir()
        (tmp_path / "P" / "notes.txt").write_text("mine")
        result = CliRunner().invoke(
            cli, ["project", "init", str(tmp_path / "P"), "--study", "a", "--no-git"]
        )
        assert result.exit_code == 2 and "is not empty" in result.output


class TestFrictionLog:
    """Regressions for the Paper 1 case-study friction log."""

    def test_copied_configs_hold_no_machine_paths(self, tmp_path: Path) -> None:
        import yaml

        config = write_simulation_config(tmp_path / "runs" / "a", scratch=tmp_path / "scratch_a")
        (config.parent / "test.pdb").write_text("END\n")
        write_openmm_replicate(config, 1, [1.0, 1.1, 1.2])
        assert (
            CliRunner()
            .invoke(cli, ["project", "init", str(tmp_path / "P"), "--study", "lipa", "--no-git"])
            .exit_code
            == 0
        )
        study = tmp_path / "P" / "lipa"
        result = CliRunner().invoke(
            cli,
            [
                "study",
                "add-condition",
                "No polymer",
                "--config",
                str(config),
                "--study",
                str(study),
            ],
        )
        assert result.exit_code == 0, result.output
        copied = (study / "conditions" / "no_polymer" / "config.yaml").read_text()
        assert str(tmp_path) not in copied
        data = yaml.safe_load((study / "data.local.yaml").read_text())
        assert data == {"No polymer": str((tmp_path / "scratch_a").resolve())}
        assert load_study_file(study).conditions["No polymer"].is_file()
        assert pz.Study(study)["No polymer"].replicates[0].index == 1

    def test_results_warn_when_the_report_does_not_match(self, project: Path) -> None:
        lipa = str(project / "lipa")
        assert _analyze("rg", "--study", lipa, "--replicates", "1").exit_code == 0
        assert _analyze("rg", "--study", lipa, "--replicates", "2", "--task").exit_code == 0
        with pytest.warns(UserWarning, match="stored but not in the report"):
            stored = pz.Study(project / "lipa").results("rg")
        assert any("replicates 2" in note for note in stored.warnings)
        assert _analyze("rg", "--study", lipa).exit_code == 0
        assert pz.Study(project / "lipa").results("rg").warnings == []

    def test_study_lists_its_conditions(self, project: Path) -> None:
        study = pz.Study(project / "lipa")
        assert [c.label for c in study.conditions] == study.labels == ["No polymer", "Half"]

    def test_relative_file_settings_resolve_against_the_study(
        self, project: Path, monkeypatch
    ) -> None:
        text = (
            (project / "lipa" / "study.yaml")
            .read_text()
            .replace(
                "analyses: {}",
                "analyses:\n"
                "  own:\n"
                "    function: analyses/own.py:own\n"
                "    kind: per_replicate\n"
                "    selections: {atoms: all}\n"
                "    settings: {ref: structures/ref.pdb, other: structures/nothing.pdb}\n",
            )
        )
        (project / "lipa" / "study.yaml").write_text(text)
        (project / "lipa" / "analyses").mkdir()
        (project / "lipa" / "analyses" / "own.py").write_text(
            "def own(atoms, ref, other, frames):\n    open(ref).read()\n    return 1.0\n"
        )
        entry = load_study_file(project / "lipa").analyses["own"].function
        assert entry.settings["ref"] == str((project / "lipa/structures/ref.pdb").resolve())
        check = CliRunner().invoke(cli, ["study", "check", str(project / "lipa")])
        assert "setting structures/nothing.pdb names no file" in check.output
        monkeypatch.chdir(project.parent)  # not the study folder
        assert _analyze("own", "--study", str(project / "lipa")).exit_code == 0

    def test_per_replicate_functions_can_ask_for_times(self, tmp_path: Path) -> None:
        from tests.analyses.test_study_hardening import TestMissingLabels

        root = TestMissingLabels()._study(tmp_path, None)
        (root / "analyses" / "m.py").write_text(
            "def per_frame(atoms, frames, times):\n"
            "    return [f'{t:.1f}' for t in times], [1.0 for _ in times]\n"
        )
        text = (root / "study.yaml").read_text() + "until: common\n"
        (root / "study.yaml").write_text(text)
        result = _analyze("byf", "--study", str(root))
        assert result.exit_code == 0, result.output
        labels = set(pz.Study(root).results("byf").table["label"])
        assert labels == {"0.0", "0.1", "0.2", "0.3"}
