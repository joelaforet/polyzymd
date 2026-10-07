"""project freeze: what each study's manifest in a project lists."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.analyses.test_project import project  # noqa: F401  (fixture)

pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore"), pytest.mark.usefixtures("git_identity")]


def test_project_study_manifests_list_what_git_tracks(project: Path) -> None:  # noqa: F811
    """A study's manifest in a project lists the files the deposit holds."""
    from polyzymd.analyses.project_freeze import freeze_project
    from polyzymd.analyses.study_git import init_repository

    (project / ".gitignore").write_text("__pycache__/\ndata.local.yaml\n")
    init_repository(project, "start")
    (project / "lipa" / "__pycache__").mkdir()
    (project / "lipa" / "__pycache__" / "x.pyc").write_bytes(b"x")
    freeze_project(project)
    manifest = json.loads((project / "lipa" / "manifest.json").read_text())
    assert not any("__pycache__" in name for name in manifest["files"])
    assert manifest["git"]["parent_commit"]


def test_a_stray_file_is_named_once(project: Path) -> None:  # noqa: F811
    """Project freeze names a file it does not deposit in one warning line."""
    from polyzymd.analyses.project_freeze import freeze_project
    from polyzymd.analyses.study_git import init_repository

    (project / ".gitignore").write_text("data.local.yaml\n")
    (project / "lipa" / "notes.txt").write_text("private\n")
    init_repository(project, "start")
    warnings = freeze_project(project).warnings
    assert len([w for w in warnings if "not deposited" in w]) == 1, warnings


def test_every_deposited_project_file_matches_a_manifest_entry(project: Path) -> None:  # noqa: F811
    """Each file in the project zip has the size and SHA-256 a manifest gives for it.

    One condition config keeps absolute directories, which the deposit removes,
    and one has relative directories, which it keeps.
    """
    import hashlib
    import zipfile

    import yaml

    from polyzymd.analyses.project_freeze import freeze_project
    from polyzymd.analyses.study_git import init_repository

    relative = project / "rml" / "conditions" / "half" / "config.yaml"
    data = yaml.safe_load(relative.read_text())
    data["output"] = {"projects_directory": "../../runs/half", "scratch_directory": None}
    relative.write_text(yaml.safe_dump(data, sort_keys=False))
    (project / ".gitignore").write_text("data.local.yaml\n")
    init_repository(project, "start")
    result = freeze_project(project)
    (archive,) = result.upload.glob("*-project-v*.zip")
    manifests = {"": json.loads((project / "manifest.json").read_text())}
    for label in ("lipa", "rml"):
        manifests[f"{label}/"] = json.loads((project / label / "manifest.json").read_text())
    with zipfile.ZipFile(archive) as opened:
        members = {
            name.removeprefix("study/"): opened.read(name)
            for name in opened.namelist()
            if not name.endswith("/")
        }
    assert "rml/conditions/half/config.yaml" in members
    assert b"../../runs/half" in members["rml/conditions/half/config.yaml"]
    assert str(project.parent).encode() not in members["lipa/conditions/half/config.yaml"]
    for name, content in members.items():
        if Path(name).name == "manifest.json":
            continue
        prefix = f"{Path(name).parts[0]}/" if f"{Path(name).parts[0]}/" in manifests else ""
        entry = manifests[prefix]["files"].get(name.removeprefix(prefix))
        assert entry, f"{name} is in no manifest"
        assert entry == {"size": len(content), "sha256": hashlib.sha256(content).hexdigest()}, name
    for path in result.upload.iterdir():
        if path.suffix != ".zip" and path.name not in ("manifest.json", "README.md"):
            if path.name != "manifest-1.schema.json":
                assert path.read_bytes() == members[path.name], path.name
