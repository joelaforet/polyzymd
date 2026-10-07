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
            "refreeze",
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

    def test_an_extended_trajectory_makes_its_run_stale(self, study: Path) -> None:
        config = study.parent / "runs" / condition_folder("Polymer") / "config.yaml"
        write_openmm_replicate(config, 2, [2.2 + 0.01 * k for k in range(12)])
        (warning,) = [w for w in freeze(study).warnings if w.startswith("run rg may be stale")]
        assert "Polymer replicate 2" in warning

    def test_record_without_trajectory_hashes_is_stale_not_an_error(self, study: Path) -> None:
        record_path = next((study / "results").rglob("record.json"))
        record = json.loads(record_path.read_text())
        segment = {k: v for k, v in record["trajectories"][0].items() if k != "sha256"}
        record["trajectories"] = [segment, dict(segment)]
        record_path.write_text(json.dumps(record))
        stale = [w for w in freeze(study).warnings if w.startswith("run rg may be stale")]
        assert stale and "no trajectory hash recorded" in stale[0]

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

    def test_locate_verify_reads_a_file_replaced_with_the_same_size_and_time(
        self, study: Path, tmp_path: Path
    ) -> None:
        import os

        result = freeze(study)
        copy = shutil.copytree(result.deposit / "study", tmp_path / "reproducer" / "my_study")
        download = shutil.copytree(tmp_path / "scratch", tmp_path / "download")
        command = ["study", "locate", str(download), "--study", str(copy), "--verify"]
        assert CliRunner().invoke(cli, command).exit_code == 0
        dcd = next(download.rglob("*.dcd"))
        stat = dcd.stat()
        dcd.write_bytes(dcd.read_bytes()[:-8] + b"\0" * 8)
        os.utime(dcd, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        located = CliRunner().invoke(cli, command)
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


def test_identical_warnings_for_several_conditions_are_one_line() -> None:
    """Conditions that share one warning give one line naming them all."""
    from polyzymd.analyses.study_freeze import group_warnings

    warnings = [
        "A: no hashes; run x",
        "B: no hashes; run x",
        "metadata.doi is not set",
        "A replicate 1: odd",
    ]
    assert group_warnings(warnings, ["A", "B"]) == [
        "A, B: no hashes; run x",
        "metadata.doi is not set",
        "A replicate 1: odd",
    ]


def test_freeze_deposits_only_names_polyzymd_chooses(tmp_path: Path) -> None:
    """Stray files in a study are neither hashed nor published, and freeze says so."""
    from polyzymd.analyses.study_freeze import _listed_files, left_out_files

    root = tmp_path / "study"
    for name in (
        "study.yaml",
        "README.md",
        "analyses/f.py",
        "analyses/data/t.csv",
        "conditions/A/config.yaml",
        "conditions/A/enzyme.pdb",
        "structures/crystal.pdb",
        "results/rg/A/replicate_1/record.json",
        "notes.txt",
        "scratch/copy.xtc",
        "scratch/more.xtc",
    ):
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_text("x")
    assert _listed_files(root, None) == [
        "README.md",
        "analyses/data/t.csv",
        "analyses/f.py",
        "conditions/A/config.yaml",
        "conditions/A/enzyme.pdb",
        "results/rg/A/replicate_1/record.json",
        "structures/crystal.pdb",
        "study.yaml",
    ]
    message = left_out_files(root, None)
    assert message.startswith("not deposited: notes.txt, scratch/.")


def test_a_project_applies_the_deposit_rule_inside_each_study(tmp_path: Path) -> None:
    """In a project, each study's stray files are left out as in a lone study."""
    from polyzymd.analyses.study_freeze import _listed_files, left_out_files

    root = tmp_path / "paper"
    for name in ("project.yaml", "stats/plan.py", "lipa/study.yaml", "lipa/notes.docx", "todo.md"):
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_text("x")
    assert _listed_files(root, None) == ["lipa/study.yaml", "project.yaml", "stats/plan.py"]
    assert left_out_files(root, None).startswith("not deposited: lipa/notes.docx, todo.md.")


def test_the_runs_folder_is_neither_listed_nor_named_as_left_out(tmp_path: Path) -> None:
    """runs/ holds the simulations: freeze neither deposits it nor warns that it does not."""
    from polyzymd.analyses.study_freeze import EXCLUDE_MACHINE_FILES, _listed_files, left_out_files

    root = tmp_path / "paper"
    for name in ("project.yaml", "lipa/study.yaml", "runs/lipa/water/w_run1/traj.dcd"):
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_text("x")
    assert _listed_files(root, None) == ["lipa/study.yaml", "project.yaml"]
    assert left_out_files(root, None) is None
    assert ":(exclude)**/runs/**" in EXCLUDE_MACHINE_FILES


@pytest.fixture()
def gromacs_study(tmp_path: Path) -> Path:
    """A committed study of one GROMACS condition of two replicates, without polymers."""
    import MDAnalysis as mda
    import numpy as np

    from polyzymd.config.schema import SimulationConfig

    folder = tmp_path / "runs" / "water"
    config = write_simulation_config(folder, scratch=tmp_path / "scratch")
    data = yaml.safe_load(config.read_text())
    data["engine"] = "gromacs"
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    (folder / "test.pdb").write_text("REMARK input\nEND\n")
    for replicate in (1, 2):
        working = SimulationConfig.from_yaml(config).get_working_directory(replicate) / "gromacs"
        working.mkdir(parents=True)
        universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
        universe.add_TopologyAttr("names", ["C1", "C2", "C3", "C4"])
        universe.add_TopologyAttr("resnames", ["MOL"])
        universe.add_TopologyAttr("masses", [1.0] * 4)
        universe.dimensions = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]
        cross = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], dtype=np.float32)
        universe.atoms.positions = cross + 10.0
        universe.atoms.write(str(working / "system.gro"))
        with mda.Writer(str(working / "prod.xtc"), n_atoms=4, dt=100.0) as writer:
            for k in range(10):
                universe.atoms.positions = cross * (1.0 + 0.01 * k) + 10.0
                universe.trajectory.ts.frame = k  # the writer stamps time dt * frame
                writer.write(universe.atoms)
        (working / "prod.log").write_text(
            "                      :-) GROMACS - gmx mdrun, 2025.2 (-:\n"
            "GROMACS version:     2025.2\nPrecision:           mixed\n"
        )
    root = tmp_path / "my_study"
    create_study(root, conditions={"Water": config}, equilibration="0.25ns")
    text = (root / "study.yaml").read_text() + METADATA.replace("[No polymer, Polymer]", "[Water]")
    (root / "study.yaml").write_text(text)
    _git(root, "commit", "-qam", "Add metadata")
    return root


class TestGromacsFreeze:
    def test_a_gromacs_run_is_analysed_without_gromacs(
        self, gromacs_study: Path, monkeypatch
    ) -> None:
        """Freeze reads a GROMACS run on a machine with no gmx on PATH and no GMX_BIN."""
        monkeypatch.delenv("GMX_BIN", raising=False)
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
        result = freeze(gromacs_study)
        assert result.manifest["conditions"]["Water"]["replicates"]

    def test_manifest_records_the_gromacs_version(self, gromacs_study: Path) -> None:
        """The version comes from each replicate's production log."""
        result = freeze(gromacs_study)
        replicates = result.manifest["conditions"]["Water"]["replicates"]
        assert {r["simulated_with"]["gromacs_version"] for r in replicates.values()} == {"2025.2"}
        assert not any("gromacs version" in w.lower() for w in result.warnings)

    def test_checklist_gives_what_gromacs_ran(self, gromacs_study: Path) -> None:
        """4b lists the GROMACS integrator and barostat; 1d has no polymers for a study without them."""
        freeze(gromacs_study)
        checklist = yaml.safe_load((gromacs_study / "md_checklist.yaml").read_text())
        ran = checklist["4b_simulation_parameters"]["evidence"]["Water"]["gromacs_production"]
        assert ran["integrator"] == "sd" and ran["pcoupl"] == "c-rescale"
        assert "polymer" not in checklist["1d_independent_starting_configurations"]["answer"]


def _zip_names(deposit: Path, part: str) -> set[str]:
    """Return the file names in the upload zip of ``part`` (study, engine_inputs or final_frames)."""
    import zipfile

    pattern = "*-study-v*.zip" if part == "study" else f"{part}.zip"
    (archive,) = (deposit / "upload").glob(pattern)
    with zipfile.ZipFile(archive) as opened:
        return {
            name.removeprefix(f"{part}/") for name in opened.namelist() if not name.endswith("/")
        }


def test_tracked_stray_files_stay_out_of_the_deposit(study: Path) -> None:
    """The deposited study holds the files the manifest lists and the files freeze writes."""
    from polyzymd.analyses.study_freeze import GENERATED

    (study / "notes.txt").write_text("private\n")
    (study / "copied.dcd").write_bytes(b"\0" * 16)
    _git(study, "add", "notes.txt", "copied.dcd")
    _git(study, "commit", "-qm", "Notes")
    result = freeze(study)
    copied = result.deposit / "study"
    assert not (copied / "notes.txt").exists() and not (copied / "copied.dcd").exists()
    written = {name for name in GENERATED if (study / name).is_file()}
    assert _zip_names(result.deposit, "study") == set(result.manifest["files"]) | written


def test_a_dropped_replicate_leaves_the_deposit(tmp_path: Path) -> None:
    """Engine inputs and final frames of a replicate no longer in the study are not deposited."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    freeze(root)
    (root / "study.yaml").write_text((root / "study.yaml").read_text() + "replicates: [1, 2]\n")
    _git(root, "commit", "-qam", "Two replicates")
    result = freeze(root)
    for part in ("engine_inputs", "final_frames"):
        assert not list((result.deposit / part).rglob("*replicate_3*")), part
        assert not any("replicate_3" in name for name in _zip_names(result.deposit, part)), part


def test_an_absolute_input_path_inside_the_study_is_deposited_relative(
    study: Path, tmp_path: Path
) -> None:
    """A reproducer's copy of the deposit finds the input and its results are not stale."""
    config = study / "conditions" / "polymer" / "config.yaml"
    data = yaml.safe_load(config.read_text())
    data["enzyme"]["pdb_path"] = str((config.parent / data["enzyme"]["pdb_path"]).resolve())
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    _git(study, "commit", "-qam", "Absolute path")
    deposit = freeze(study).deposit
    elsewhere = tmp_path / "downloaded"
    shutil.copytree(deposit / "study", elsewhere)
    shutil.rmtree(study)
    deposited = yaml.safe_load((elsewhere / "conditions" / "polymer" / "config.yaml").read_text())
    assert not Path(deposited["enzyme"]["pdb_path"]).is_absolute()
    assert stale_runs(load_study_file(elsewhere)) == {}


def test_an_input_outside_the_study_is_named_in_a_warning(study: Path, tmp_path: Path) -> None:
    """Freeze warns that an input outside the study is not deposited, naming no machine path."""
    config = study / "conditions" / "polymer" / "config.yaml"
    data = yaml.safe_load(config.read_text())
    data["enzyme"]["pdb_path"] = str(tmp_path / "runs" / "polymer" / "test.pdb")
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    _git(study, "commit", "-qam", "Outside path")
    warnings = freeze(study).warnings
    assert any("test.pdb" in w and "outside the study" in w for w in warnings), warnings
    assert not any(str(tmp_path) in w for w in warnings)


def test_the_manifest_names_the_parent_of_the_tagged_commit(study: Path) -> None:
    """The manifest is inside the tagged commit, so it records that commit's parent."""
    freeze(study)
    manifest = json.loads((study / "manifest.json").read_text())
    assert manifest["git"]["parent_commit"] == _git(study, "rev-parse", "study-v1^").strip()
    assert "commit" not in manifest["git"]


def test_a_failed_commit_leaves_no_tag_in_the_deposit(study: Path) -> None:
    """When git cannot commit, freeze exits with an error and no file names the tag."""
    hook = study / ".git" / "hooks" / "pre-commit"
    hook.write_text("#!/bin/sh\nexit 1\n")
    hook.chmod(0o755)
    result = CliRunner().invoke(cli, ["study", "freeze", str(study)])
    assert result.exit_code != 0, result.output
    assert not _git(study, "tag").strip()
    for path in (study / "deposit").rglob("*"):
        if path.is_file() and path.suffix not in (".zip", ".gz"):
            assert "study-v1" not in path.read_text(errors="ignore"), path


def _write_segment(run_dir: Path, index: int, istart: int, n_frames: int = 5) -> None:
    """Write production_<index>: ``n_frames`` frames 100 ps apart, the first at step ``istart``."""
    import MDAnalysis as mda
    import numpy as np

    from tests._support.analysis_testkit import CROSS

    segment = run_dir / f"production_{index}"
    segment.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("names", ["C1", "C2", "C3", "C4"])
    universe.add_TopologyAttr("resnames", ["MOL"])
    universe.add_TopologyAttr("masses", [1.0] * 4)
    universe.atoms.positions = np.asarray(CROSS, dtype=np.float32)
    universe.atoms.write(str(run_dir / "solvated_system.pdb"))
    path = segment / f"production_{index}_trajectory.dcd"
    with mda.Writer(str(path), n_atoms=4, dt=100.0, istart=istart, nsavc=1) as writer:
        for k in range(n_frames):
            universe.atoms.positions = np.asarray(CROSS, dtype=np.float32) * (1.0 + 0.01 * k)
            writer.write(universe.atoms)


def test_production_length_skips_a_repeated_boundary_frame(tmp_path: Path) -> None:
    """The manifest's production length is the replicate's, with a repeated frame left out."""
    from polyzymd.analyses.study import Study
    from polyzymd.config.schema import SimulationConfig
    from tests._support.analysis_testkit import write_simulation_config

    config = write_simulation_config(tmp_path / "runs" / "toy", scratch=tmp_path / "scratch")
    (config.parent / "test.pdb").write_text("REMARK input\nEND\n")
    run_dir = SimulationConfig.from_yaml(config).get_working_directory(1)
    _write_segment(run_dir, 0, 0)
    _write_segment(run_dir, 1, 4)
    root = tmp_path / "study"
    create_study(root, conditions={"Toy": config}, equilibration="0ns")
    expected = Study.from_configs({"Toy": config}, equilibration="0ns")["Toy"].replicates[0]
    manifest = freeze(root).manifest
    recorded = manifest["conditions"]["Toy"]["replicates"]["1"]["production_ns"]
    assert recorded == pytest.approx(expected.production_ns)


def test_a_new_study_has_no_left_out_warning(study: Path) -> None:
    """The files study init writes are all deposited."""
    assert not any(w.startswith("not deposited") for w in freeze(study).warnings)


def test_an_unknown_package_version_is_not_recorded_and_the_lock_file_is_hashed(
    study: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A package reporting version 0.0.0 is recorded as unknown; the deposited pixi.lock pins it."""
    import hashlib
    import sys

    import numpy

    monkeypatch.setattr(sys, "prefix", str(study.parent / "env"))
    (study / "environment" / "pixi.lock").write_text("version: 6\n")
    _git(study, "add", "environment/pixi.lock")
    _git(study, "commit", "-qm", "Lock file")
    monkeypatch.setattr(numpy, "__version__", "0.0.0")
    versions = freeze(study).manifest["versions"]
    assert versions["numpy"] is None
    assert versions["pixi.lock"] == hashlib.sha256(b"version: 6\n").hexdigest()


def test_versions_come_from_the_pixi_environment_that_runs_freeze(
    study: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A 0.0.0 package has its conda record's version, OpenMM its full version, and the
    pixi.lock of the running workspace pins the environment when the study has none."""
    import hashlib
    import sys

    import numpy

    from polyzymd.utils.version import get_openmm_version

    workspace = study.parent / "workspace"
    prefix = workspace / ".pixi" / "envs" / "analysis"
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "conda-meta" / "numpy-9.9.1-py312_0.json").write_text(
        '{"name": "numpy", "version": "9.9.1"}'
    )
    (prefix / "conda-meta" / "numpy-base-1.0-py312_0.json").write_text(
        '{"name": "numpy-base", "version": "1.0"}'
    )
    (workspace / "pixi.lock").write_text("version: 6\n")
    monkeypatch.setattr(sys, "prefix", str(prefix))
    monkeypatch.setattr(numpy, "__version__", "0.0.0")
    versions = freeze(study).manifest["versions"]
    assert versions["numpy"] == "9.9.1"
    assert versions["openmm"] == get_openmm_version()
    assert versions["pixi.lock"] == hashlib.sha256(b"version: 6\n").hexdigest()


def test_a_non_ascii_file_name_is_deposited(study: Path) -> None:
    """git quotes non-ASCII names in its listings; freeze still lists and deposits the file."""
    (study / "results" / "résumé.csv").write_text("x\n")
    result = freeze(study)
    assert "results/résumé.csv" in result.manifest["files"]
    assert "results/résumé.csv" in _zip_names(result.deposit, "study")


def test_a_file_name_with_brackets_deposits_only_that_file(tmp_path: Path) -> None:
    """A listed name such as a[1].csv is a file name, not a pattern that also matches a1.csv."""
    from polyzymd.analyses.study_freeze import _copy_frozen_folder

    root = tmp_path / "study"
    (root / "results").mkdir(parents=True)
    for name in ("a[1].csv", "a1.csv"):
        (root / "results" / name).write_text("x\n")
    _git(root, "init", "-q")
    _git(root, "add", ".")
    _git(root, "commit", "-qm", "Results")
    _git(root, "tag", "v1")
    commit = _git(root, "rev-parse", "HEAD").strip()
    _copy_frozen_folder(root, tmp_path / "deposit", "v1", commit, ["results/a[1].csv"], [])
    copied = tmp_path / "deposit" / "study" / "results"
    assert sorted(p.name for p in copied.iterdir()) == ["a[1].csv"]


def test_the_gitignore_written_by_a_first_freeze_is_deposited(study: Path) -> None:
    """A study without .gitignore gets one from freeze, committed and deposited."""
    _git(study, "rm", "-q", ".gitignore")
    _git(study, "commit", "-qm", "No .gitignore")
    result = freeze(study)
    assert ".gitignore" in _git(study, "ls-tree", "--name-only", "study-v1")
    assert ".gitignore" in _zip_names(result.deposit, "study")


def test_the_manifest_lists_the_gitignore_written_by_a_first_freeze(study: Path) -> None:
    """The study zip holds the files manifest.json lists and the files freeze writes."""
    from polyzymd.analyses.study_freeze import GENERATED

    _git(study, "rm", "-q", ".gitignore")
    _git(study, "commit", "-qm", "No .gitignore")
    result = freeze(study)
    assert ".gitignore" in result.manifest["files"]
    assert _zip_names(result.deposit, "study") - {"manifest.json"} == set(result.manifest["files"])


def test_freeze_refuses_a_condition_config_outside_the_study(tmp_path: Path) -> None:
    """Analysis reads a config outside the study, but freeze could not deposit it."""
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    outside = tmp_path / "runs" / "polymer" / "config.yaml"
    text = (root / "study.yaml").read_text().replace("conditions/polymer/config.yaml", str(outside))
    (root / "study.yaml").write_text(text)
    _git(root, "commit", "-qam", "Polymer config outside")
    result = CliRunner().invoke(
        cli, ["analyze", "rg", "--study", str(root), "--no-eq-check", "--no-plots"]
    )
    assert result.exit_code == 0, result.output
    check = CliRunner().invoke(cli, ["study", "check", str(root)])
    assert "freeze will refuse" in check.output
    with pytest.raises(ProtocolError, match="outside the study") as caught:
        freeze(root)
    assert "conditions:" in caught.value.hint
    assert 'polyzymd study add-condition "Polymer" --config' in caught.value.hint


def test_a_study_whose_gitignore_lacks_runs_freezes_with_a_run_in_it(study: Path) -> None:
    """A study made before runs/ was git-ignored: a trajectory in runs/ is not an input."""
    gitignore = study / ".gitignore"
    gitignore.write_text(
        "".join(line for line in gitignore.read_text().splitlines(True) if "runs" not in line)
    )
    _git(study, "commit", "-qam", "Old .gitignore")
    (study / "runs" / "water" / "w_run1").mkdir(parents=True)
    (study / "runs" / "water" / "w_run1" / "traj.dcd").write_text("x")

    assert freeze(study).tag == "study-v1"
    assert "runs/" in gitignore.read_text().splitlines()
    assert "runs/water/w_run1/traj.dcd" not in _git(study, "ls-tree", "-r", "--name-only", "HEAD")


def test_freeze_ignores_the_polymer_cache_of_a_build(study: Path) -> None:
    """A dynamic polymer build writes .polymer_cache/ into the folder it runs in."""
    gitignore = study / ".gitignore"
    gitignore.write_text(
        "".join(line for line in gitignore.read_text().splitlines(True) if "polymer" not in line)
    )
    _git(study, "commit", "-qam", "Old .gitignore")
    (study / ".polymer_cache").mkdir()
    (study / ".polymer_cache" / "chain.sdf").write_text("x")

    assert freeze(study).tag == "study-v1"
    assert ".polymer_cache/" in gitignore.read_text().splitlines()
    assert ".polymer_cache/chain.sdf" not in _git(study, "ls-tree", "-r", "--name-only", "HEAD")


def test_a_deposited_config_keeps_relative_directories() -> None:
    """Only absolute directories and the path of an old copy header leave a deposited config."""
    body = (
        "#   the runs go into ../../runs/water (relative to this file) unless you set\n"
        "#   scratch_directory in config.yaml.\n"
        "name: water\n"
        "output:\n"
        '  projects_directory: "../../runs/water"\n'
        "  scratch_directory: null\n"
    )
    copied = "# Copied by polyzymd from config.yaml\n" + body
    assert without_machine_paths(copied) == copied
    old = "# Copied by polyzymd study init from /home/u/runs/config.yaml\n" + body
    assert without_machine_paths(old) == copied
    absolute = copied.replace("null", "/scratch/u/water")
    assert yaml.safe_load(without_machine_paths(absolute))["output"] == {
        "projects_directory": "../../runs/water",
        "scratch_directory": "data",
    }


def test_every_deposited_file_matches_its_manifest_entry(study: Path, tmp_path: Path) -> None:
    """The deposit verifies against its manifest, configs with absolute directories included.

    The polymer config gets absolute directories, which the deposit removes;
    the other keeps the relative ones study init wrote, which it keeps.
    """
    import hashlib
    import zipfile

    config = study / "conditions" / "polymer" / "config.yaml"
    data = yaml.safe_load(config.read_text())
    data["output"] = {
        "projects_directory": str(tmp_path / "jobs"),
        "scratch_directory": str(tmp_path / "scratch" / "polymer"),
    }
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    kept = (study / "conditions" / "no_polymer" / "config.yaml").read_text()
    _git(study, "commit", "-qam", "Absolute directories")
    result = freeze(study)
    manifest = result.manifest

    def entry(content: bytes) -> dict:
        return {"size": len(content), "sha256": hashlib.sha256(content).hexdigest()}

    (archive,) = result.upload.glob("*-study-v*.zip")
    with zipfile.ZipFile(archive) as opened:
        members = {
            name.removeprefix("study/"): opened.read(name)
            for name in opened.namelist()
            if not name.endswith("/")
        }
    assert members["conditions/no_polymer/config.yaml"].decode() == kept
    assert str(tmp_path).encode() not in members["conditions/polymer/config.yaml"]
    for name in ("md_checklist.yaml", "system_summary.csv", "CITATION.cff", ".zenodo.json"):
        assert name in manifest["files"], name
    for name, content in members.items():
        if name != "manifest.json":
            assert manifest["files"].get(name) == entry(content), name
    replicates = [r for c in manifest["conditions"].values() for r in c["replicates"].values()]
    recorded = [
        {"size": f["size"], "sha256": f["sha256"]}
        for r in replicates
        for f in [*r.get("engine_inputs", []), r.get("final_frame")]
        if f
    ]
    for part in ("engine_inputs", "final_frames"):
        with zipfile.ZipFile(result.upload / f"{part}.zip") as opened:
            for name in opened.namelist():
                if not name.endswith("/"):
                    assert entry(opened.read(name)) in recorded, name
    assert (result.upload / "CITATION.cff").read_bytes() == members["CITATION.cff"]
