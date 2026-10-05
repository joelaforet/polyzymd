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
