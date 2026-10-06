"""polyzymd study freeze: metadata, citation files, manifest, checklist, git tag and deposit.

The study is built with study init from two synthetic OpenMM conditions of
two replicates (four unit-mass atoms on a cross), as in test_study_folder.py.
"""

from __future__ import annotations

import json
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
from polyzymd.analyses.study_metadata import (
    check_metadata,
    citation_cff,
    is_placeholder,
    zenodo_json,
)
from polyzymd.analyses.study_scaffold import condition_folder, create_study
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
    pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed"),
]

METADATA = """\
metadata:
  title: A test study
  description: Two synthetic conditions.
  purpose: To test study freeze.
  keywords: [molecular dynamics, test]
  system_type: [toy]
  authors:
    - {family-names: Lovelace, given-names: Ada, orcid: "0000-0002-1825-0097", affiliation: Analytical Engines}
  license: {data: CC-BY-4.0, code: MIT}
  related:
    paper: {title: A paper, status: in-preparation, doi: "10.XXXX/placeholder"}
    trajectories:
      - {doi: "10.5281/zenodo.1234567", conditions: [No polymer, Polymer]}
"""


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


@pytest.fixture()
def study(tmp_path: Path) -> Path:
    """An analysed study folder with metadata, all committed."""
    configs = {}
    for label, offset in (("No polymer", 1.0), ("Polymer", 2.0)):
        folder = tmp_path / "runs" / condition_folder(label)
        config = write_simulation_config(folder, scratch=tmp_path / "scratch" / folder.name)
        (folder / "test.pdb").write_text("REMARK input\nEND\n")
        for replicate in (1, 2):
            write_openmm_replicate(
                config, replicate, [offset + 0.1 * replicate + 0.01 * k for k in range(10)]
            )
        configs[label] = config
    root = tmp_path / "my_study"
    create_study(root, conditions=configs, equilibration="0.25ns")
    text = (
        (root / "study.yaml")
        .read_text()
        .replace("analyses: {}", "analyses:\n  rg: {selection: all}")
    )
    (root / "study.yaml").write_text(text + METADATA)
    _git(root, "commit", "-qam", "Add analyses and metadata")
    result = CliRunner().invoke(
        cli, ["analyze", "rg", "--study", str(root), "--no-eq-check", "--no-plots"]
    )
    assert result.exit_code == 0, result.output
    return root


class TestMetadata:
    def test_complete_block_has_no_gaps_but_the_placeholder_doi(self) -> None:
        meta, warnings = check_metadata(yaml.safe_load(METADATA)["metadata"])
        assert warnings == [
            "metadata.doi is not set: reserve a DOI for the study in Zenodo, add it here and "
            "refreeze (deposit/UPLOAD.md says how)",
            "the paper DOI is missing or a placeholder; refreeze once it is known",
        ]
        assert meta["license"] == {"data": "CC-BY-4.0", "code": "MIT"}

    def test_gaps_become_todos(self) -> None:
        meta, warnings = check_metadata({})
        assert meta["title"].startswith("TODO") and meta["authors"][0]["name"].startswith("TODO")
        assert meta["license"] == {"data": "CC-BY-4.0", "code": "MIT"}
        assert any("purpose" in w for w in warnings) and any("trajectories" in w for w in warnings)

    @pytest.mark.parametrize(
        ("raw", "message"),
        [
            ({"titel": "x"}, "unknown key 'titel'"),
            ({"authors": [{"given-names": "Ada"}]}, "no family-names"),
            ({"authors": [{"family-names": "L", "orcd": "x"}]}, "unknown key 'orcd'"),
            ({"related": {"paper": {"status": "published"}}}, "not one of"),
            ({"zenodo": {"access_right": "public"}}, "not one of"),
            ({"keywords": "md"}, "must be a list"),
        ],
    )
    def test_refusals(self, raw: dict, message: str) -> None:
        with pytest.raises(ProtocolError, match=message):
            check_metadata(raw)

    def test_placeholder(self) -> None:
        assert is_placeholder("10.XXXX/x") and is_placeholder(None) and is_placeholder("TODO")
        assert not is_placeholder("10.5281/zenodo.1234567")

    def test_citation_cff_cites_the_paper_and_polyzymd(self) -> None:
        meta, _ = check_metadata(yaml.safe_load(METADATA)["metadata"])
        cff = citation_cff(meta, version="study-v1", released="2026-01-01", commit="abc")
        assert cff["cff-version"] == "1.2.0" and cff["type"] == "dataset"
        assert cff["preferred-citation"]["title"] == "A paper"
        assert "doi" not in cff["preferred-citation"]
        assert cff["authors"][0]["orcid"] == "https://orcid.org/0000-0002-1825-0097"
        titles = [r["title"] for r in cff["references"]]
        assert any("PolyzyMD" in t for t in titles)
        assert any(r.get("doi") == "10.5281/zenodo.1234567" for r in cff["references"])

    def test_zenodo_json(self) -> None:
        meta, _ = check_metadata(yaml.safe_load(METADATA)["metadata"])
        zenodo = zenodo_json(meta, version="study-v1", released="2026-01-01", method="m")
        assert zenodo["upload_type"] == "dataset" and zenodo["license"] == "cc-by-4.0"
        assert zenodo["creators"] == [
            {
                "name": "Lovelace, Ada",
                "affiliation": "Analytical Engines",
                "orcid": "0000-0002-1825-0097",
            }
        ]
        relations = {(r["relation"], r["identifier"]) for r in zenodo["related_identifiers"]}
        assert ("requires", "https://github.com/joelaforet/polyzymd") in relations
        assert ("references", "10.5281/zenodo.1234567") in relations
        assert not any(r["relation"] == "isSupplementTo" for r in zenodo["related_identifiers"])


class TestStale:
    def test_fresh_results_are_not_stale(self, study: Path) -> None:
        assert stale_runs(load_study_file(study)) == {}

    def test_changed_window_and_missing_run(self, study: Path) -> None:
        text = (
            (study / "study.yaml")
            .read_text()
            .replace("equilibration: 0.25ns", "equilibration: 0ns")
        )
        (study / "study.yaml").write_text(
            text.replace("  rg: {selection: all}", "  rg: {selection: all}\n  rg2: {analysis: rg}")
        )
        reasons = stale_runs(load_study_file(study))
        assert any("equilibration" in r for r in reasons["rg"])
        assert reasons["rg2"] == ["no stored results; run polyzymd analyze rg2 --study"]

    def test_own_function_with_settings_is_not_stale(self, study: Path) -> None:
        (study / "analyses" / "metrics.py").write_text(
            "def scaled_rg(atoms, frames, factor=1.0):\n"
            "    total = 0.0\n"
            "    for _ in atoms.universe.trajectory[frames]:\n"
            "        total += atoms.radius_of_gyration()\n"
            "    return factor * total / len(frames)\n"
        )
        text = (
            (study / "study.yaml")
            .read_text()
            .replace(
                "  rg: {selection: all}",
                "  rg: {selection: all}\n  scaled:\n    function: analyses/metrics.py:scaled_rg\n"
                "    kind: per_replicate\n    selections: {atoms: all}\n    settings: {factor: 2.0}",
            )
        )
        (study / "study.yaml").write_text(text)
        result = CliRunner().invoke(
            cli, ["analyze", "scaled", "--study", str(study), "--no-eq-check", "--no-plots"]
        )
        assert result.exit_code == 0, result.output
        report = json.loads((study / "results" / "scaled" / "report.json").read_text())
        assert report["provenance"]["study"]["settings"] == {"factor": 2.0}
        assert report["provenance"]["study"]["path"] == "study.yaml"
        assert "scaled" not in stale_runs(load_study_file(study))

    def test_changed_setting(self, study: Path) -> None:
        text = (
            (study / "study.yaml")
            .read_text()
            .replace("rg: {selection: all}", "rg: {selection: index 0}")
        )
        (study / "study.yaml").write_text(text)
        assert any("setting selection" in r for r in stale_runs(load_study_file(study))["rg"])


class TestFreeze:
    def test_writes_commits_and_tags(self, study: Path) -> None:
        result = freeze(study)
        assert result.tag == "study-v1" and result.commit
        assert _git(study, "tag").split() == ["study-v1"]
        committed = _git(study, "show", "--name-only", "--format=", "HEAD").split()
        for name in (
            "manifest.json",
            "CITATION.cff",
            ".zenodo.json",
            "md_checklist.yaml",
            "system_summary.csv",
        ):
            assert name in committed
        assert any(p.startswith("results/rg/") for p in committed)
        assert "deposit/" in (study / ".gitignore").read_text()
        assert _git(study, "status", "--porcelain").strip() == ""

    def test_manifest(self, study: Path) -> None:
        manifest = freeze(study).manifest
        assert manifest["schema"] == "polyzymd-study-manifest/1"
        replicate = manifest["conditions"]["Polymer"]["replicates"]["2"]
        assert "frames_analysed" not in replicate and replicate["production_frames"] == 10
        assert all(len(f["sha256"]) == 64 for f in replicate["files"])
        assert not any(f["path"].startswith("/") for f in replicate["files"])
        assert manifest["conditions"]["Polymer"]["resolved_config"]["enzyme"][
            "pdb_path"
        ].startswith("conditions/")
        assert "study.yaml" in manifest["files"]
        assert str(study.parent) not in json.dumps(manifest)

    def test_deposit_layout(self, study: Path) -> None:
        result = freeze(study)
        deposit = result.deposit
        for name in (
            "manifest.json",
            "CITATION.cff",
            ".zenodo.json",
            "README.md",
            "UPLOAD.md",
            "trajectories.csv",
            "study/study.yaml",
            "study/results/rg/report.json",
        ):
            assert (deposit / name).is_file(), name
        assert list((deposit / "engine_inputs" / "polymer" / "replicate_1").glob("*.gz"))
        assert (deposit / "final_frames" / "polymer" / "replicate_1_final.pdb.gz").is_file()
        assert not (deposit / "study" / "data.local.yaml").exists()
        assert result.guide == deposit / "UPLOAD.md" and result.upload == deposit / "upload"

    def test_summary_table(self, study: Path) -> None:
        freeze(study)
        rows = (study / "system_summary.csv").read_text().splitlines()
        assert rows[0].startswith("condition,replicate,atoms")
        assert len(rows) == 1 + 4

    def test_checklist(self, study: Path) -> None:
        freeze(study)
        checklist = yaml.safe_load((study / "md_checklist.yaml").read_text())
        assert "Communications Biology" in checklist["source"]
        assert checklist["1c_replicates_and_statistics"]["evidence"][
            "replicates_per_condition"
        ] == {
            "No polymer": 2,
            "Polymer": 2,
        }

    def test_unrestrained_study_is_unbiased(self, study: Path) -> None:
        freeze(study)
        checklist = yaml.safe_load((study / "md_checklist.yaml").read_text())
        assert checklist["3c_enhanced_sampling"]["answer"].startswith("unbiased")

    def test_distance_restraints_are_reported_as_biased_sampling(self, study: Path) -> None:
        config = study / "conditions" / "polymer" / "config.yaml"
        data = yaml.safe_load(config.read_text())
        data["restraints"] = [
            {
                "type": "flat_bottom",
                "name": "ligand_in_pocket",
                "atom1": {"selection": "index 0"},
                "atom2": {"selection": "index 1"},
                "distance": 4.0,
                "force_constant": 5000.0,
            }
        ]
        config.write_text(yaml.safe_dump(data, sort_keys=False))
        _git(study, "commit", "-qam", "Restrain the ligand")
        manifest = freeze(study).manifest
        restraint = manifest["conditions"]["Polymer"]["restraints"][0]
        assert restraint["type"] == "flat_bottom" and restraint["distance_A"] == 4.0
        assert manifest["conditions"]["No polymer"]["restraints"] == []
        checklist = yaml.safe_load((study / "md_checklist.yaml").read_text())
        sampling = checklist["3c_enhanced_sampling"]
        assert sampling["answer"].startswith("restrained molecular dynamics")
        assert list(sampling["evidence"]) == ["Polymer"]
        assert sampling["evidence"]["Polymer"][0]["name"] == "ligand_in_pocket"
        assert checklist["4b_simulation_parameters"]["evidence"]["Polymer"]["restraints"]

    def test_refreeze_makes_the_next_tag(self, study: Path) -> None:
        freeze(study)
        assert freeze(study).tag == "study-v2"
        with pytest.raises(ProtocolError, match="already exists"):
            freeze(study, tag="study-v1")

    def test_uncommitted_inputs_are_refused(self, study: Path) -> None:
        """The tag and deposit hold committed files, so the manifest may describe only those."""
        (study / "analyses" / "draft.py").write_text("x = 1\n")
        with pytest.raises(ProtocolError, match="uncommitted input") as info:
            freeze(study)
        assert "analyses/draft.py" in str(info.value) and "git" in info.value.hint
        assert not _git(study, "tag")

    def test_without_trajectories(self, study: Path, tmp_path: Path) -> None:
        shutil.rmtree(tmp_path / "scratch")
        result = freeze(study)
        assert result.tag
        assert any("not on this machine" in w for w in result.warnings)

    def test_cli(self, study: Path) -> None:
        result = CliRunner().invoke(cli, ["study", "freeze", str(study), "--tag", "paper-v1"])
        assert result.exit_code == 0, result.output
        assert "as paper-v1" in result.output and "warning: the paper DOI" in result.output


class TestReproduce:
    def test_locate_verifies_against_the_manifest(self, study: Path, tmp_path: Path) -> None:
        result = freeze(study)
        copy = shutil.copytree(result.deposit / "study", tmp_path / "reproducer" / "my_study")
        download = tmp_path / "download"
        shutil.copytree(tmp_path / "scratch", download)
        located = CliRunner().invoke(
            cli, ["study", "locate", str(download), "--study", str(copy), "--verify"]
        )
        assert located.exit_code == 0, located.output
        assert "files match manifest.json (SHA-256)" in located.output
        assert pz.Study(copy).results("rg").report.conditions[1].mean == pytest.approx(2.21)

    def test_locate_reports_a_changed_file(self, study: Path, tmp_path: Path) -> None:
        result = freeze(study)
        copy = shutil.copytree(result.deposit / "study", tmp_path / "reproducer" / "my_study")
        download = shutil.copytree(tmp_path / "scratch", tmp_path / "download")
        dcd = next(download.rglob("*.dcd"))
        dcd.write_bytes(dcd.read_bytes()[:-8] + b"\0" * 8)
        located = CliRunner().invoke(
            cli, ["study", "locate", str(download), "--study", str(copy), "--verify"]
        )
        assert located.exit_code == 2
        assert "has another SHA-256" in located.output

    def test_check_reports_metadata_gaps_and_the_next_step(self, study: Path) -> None:
        result = CliRunner().invoke(cli, ["study", "check", str(study)])
        assert "metadata (study.yaml): 2 gaps" in result.output
        assert "publish: when the analyses are final, run polyzymd study freeze" in result.output
        freeze(study)
        result = CliRunner().invoke(cli, ["study", "check", str(study)])
        assert "publish: follow" in result.output and "UPLOAD.md" in result.output

    def test_study_doi_reaches_the_citation_files(self, study: Path) -> None:
        text = (
            (study / "study.yaml")
            .read_text()
            .replace("metadata:\n", 'metadata:\n  doi: "10.5281/zenodo.7654321"\n')
        )
        (study / "study.yaml").write_text(text)
        _git(study, "commit", "-qam", "Add the reserved DOI")
        result = freeze(study)
        assert not any("metadata.doi" in w for w in result.warnings)
        assert (
            yaml.safe_load((study / "CITATION.cff").read_text())["doi"] == "10.5281/zenodo.7654321"
        )
        assert json.loads((study / ".zenodo.json").read_text())["doi"] == "10.5281/zenodo.7654321"
        assert "already set: `10.5281/zenodo.7654321`" in result.guide.read_text()


def test_freeze_needs_a_git_identity_before_writing(tmp_path: Path, monkeypatch) -> None:
    """Without a git name and email nothing is written that names a tag."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    for key in ("GIT_AUTHOR_NAME", "GIT_AUTHOR_EMAIL", "GIT_COMMITTER_NAME", "GIT_COMMITTER_EMAIL"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    with pytest.raises(ProtocolError, match="no user name and email"):
        freeze(root)
    assert not (root / "manifest.json").exists()


def test_a_flow_style_output_loses_its_machine_paths() -> None:
    """Flow-style and block-scalar output paths are rewritten through YAML."""
    flow = "output: {projects_directory: /home/u/p, scratch_directory: /scratch/u/r}\nx: 1\n"
    block = "output:\n  projects_directory: >-\n    /home/u/p\n  scratch_directory: /s/u\n"
    for text in (flow, block):
        output = yaml.safe_load(without_machine_paths(text))["output"]
        assert output == {"projects_directory": ".", "scratch_directory": "data"}


def test_the_freeze_warning_names_no_machine_path(tmp_path: Path) -> None:
    """The manifest is published, so its warnings hold no scratch path."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    shutil.rmtree(tmp_path / "scratch")
    result = freeze(root)
    assert not any(str(tmp_path) in warning for warning in result.warnings)


def test_a_study_of_a_project_is_frozen_with_the_project(tmp_path: Path) -> None:
    """study freeze inside a project refuses, since it would leave out project.yaml and analyses/."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n", name="lipa")
    (tmp_path / "project.yaml").write_text("studies: {lipa: lipa}\n")
    with pytest.raises(ProtocolError, match="is a study of the project") as info:
        freeze(root)
    assert "polyzymd project freeze" in info.value.hint


def test_freeze_names_cosolvents_and_missing_build_files(tmp_path: Path) -> None:
    """Freeze recognises co-solvents in the composition check and names a listed build file that is missing."""
    from types import SimpleNamespace

    from polyzymd.analyses.study_freeze import _missing_build_files, composition_warnings

    config = SimpleNamespace(
        substrate=None,
        polymers=None,
        solvent=SimpleNamespace(co_solvents=[SimpleNamespace(name="sds", residue_name="SDS")]),
    )

    class Residues:
        resnames = ["SDS", "SDS"]

    universe = SimpleNamespace(select_atoms=lambda selection: SimpleNamespace(residues=Residues()))
    assert composition_warnings("SDS", config, universe) == []
    (tmp_path / "build_manifest.json").write_text(
        json.dumps({"artifacts": {"system.prmtop": {}, "system.xml": {}}})
    )
    (tmp_path / "system.xml").write_text("<x/>")
    assert _missing_build_files(tmp_path) == ["system.prmtop"]
