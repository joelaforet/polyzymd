"""Projects: one paper's studies, one per protein, with shared analyses (projects.md, slice P1)."""

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
        "polyzymd: 1.3.0\n"
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


def _old_study(tmp_path: Path, name: str) -> Path:
    """A study as Paper_1_REDO's: absolute condition paths and an absolute reference file."""
    outside = tmp_path / "configs" / name
    configs = {}
    for condition, offset in (("none", 1.0), ("half", 2.0)):
        config = write_simulation_config(
            outside / condition, scratch=tmp_path / "data" / name / condition
        )
        for replicate in (1, 2):
            write_openmm_replicate(config, replicate, [offset + 0.01 * k for k in range(5)])
        configs[condition] = config
    reference = tmp_path / "refs" / f"{name}_crystal.pdb"
    reference.parent.mkdir(parents=True, exist_ok=True)
    reference.write_text(f"REMARK {name}\nEND\n")
    root = tmp_path / "old" / name
    (root / "analyses").mkdir(parents=True)
    (root / "analyses" / "shared.py").write_text(
        "def size(atoms, reference):\n    return float(len(atoms))\n"
    )
    (root / "study.yaml").write_text(
        "equilibration: 0ns\n"
        f"conditions: {{No polymer: {configs['none']}, Half: {configs['half']}}}\n"
        "analyses:\n"
        "  rg: {selection: all}\n"
        "  size:\n"
        "    function: analyses/shared.py:size\n"
        "    kind: timeseries\n"
        "    selections: {atoms: all}\n"
        f"    settings: {{reference: {reference}}}\n"
        "metadata: {title: Old paper}\n"
    )
    return root


class TestProjectInit:
    def test_new_project_scaffold(self, tmp_path: Path) -> None:
        result = CliRunner().invoke(
            cli,
            ["project", "init", str(tmp_path / "P"), "--study", "lipa363", "--study", "rml333"]
            + ["--no-git"],
        )
        assert result.exit_code == 0, result.output
        text = (tmp_path / "P" / "project.yaml").read_text()
        assert "studies:\n  lipa363: lipa363\n  rml333: rml333" in text
        assert "title: TODO" in text
        study = (tmp_path / "P" / "lipa363" / "study.yaml").read_text()
        assert "description: TODO" in study and "regions: {}" in study

    def test_moving_studies_in_keeps_their_results(self, tmp_path: Path) -> None:
        import yaml

        old = {name: _old_study(tmp_path, name) for name in ("lipa", "rml")}
        for path in old.values():
            assert _analyze("--study", str(path)).exit_code == 0
        result = CliRunner().invoke(
            cli,
            ["project", "init", str(tmp_path / "P")]
            + [arg for name, path in old.items() for arg in ("--study", f"{name}={path}")]
            + ["--no-git"],
        )
        assert result.exit_code == 0, result.output
        project = tmp_path / "P"
        lipa = yaml.safe_load((project / "lipa" / "study.yaml").read_text())
        assert lipa["structures"] == {"reference": "structures/lipa_crystal.pdb"}
        assert lipa["conditions"]["Half"] == "conditions/half/config.yaml"
        shared = yaml.safe_load((project / "project.yaml").read_text())
        assert set(shared["analyses"]) == {"rg", "size"}
        assert shared["analyses"]["size"]["settings"] == {"reference": "structure reference"}
        assert shared["metadata"] == {"title": "Old paper"}
        assert (project / "analyses" / "shared.py").is_file()
        data = yaml.safe_load((project / "lipa" / "data.local.yaml").read_text())
        assert data["Half"] == str((tmp_path / "data" / "lipa" / "half").resolve())
        records = sorted((project / "lipa" / "results").rglob("record.json"))
        before = {p: p.stat().st_mtime_ns for p in records}
        assert records and _analyze("--project", str(project)).exit_code == 0
        # Every stored result was reused: nothing was measured again.
        assert {p: p.stat().st_mtime_ns for p in records} == before
        assert (old["lipa"] / "study.yaml").read_text().startswith("equilibration: 0ns")


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
