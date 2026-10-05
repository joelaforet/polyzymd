"""Freeze a project for publication: every study, one manifest, one citation, one deposit.

Slice P3 of the "Projects and studies" design. :func:`freeze_project`

1. freezes each study as :func:`~polyzymd.analyses.study_freeze.freeze` does
   with ``publish=False``: its ``manifest.json``, ``md_checklist.yaml``,
   ``system_summary.csv``, engine inputs and final frames, and its warnings,
   each prefixed with the study's label;
2. writes the project's ``manifest.json``, which lists each study's manifest
   by SHA-256 and every condition as ``<study> / <condition>``, and its
   ``CITATION.cff`` and ``.zenodo.json`` from ``project.yaml``'s metadata;
3. commits the generated files and the results of every study, and tags the
   project ``project-v<n>``;
4. lays out ``deposit/`` at the project root: the tagged project (configs
   without machine paths), the engine inputs and final frames of every study
   under its label, and ``upload/`` with ``UPLOAD.md``, as a study's deposit.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.study_freeze import (
    CHECKLIST,
    CITATION,
    DEPOSIT,
    MANIFEST,
    MANIFEST_SCHEMA_FILE,
    SUMMARY,
    ZENODO,
    FreezeResult,
    _git,
    _versions,
    freeze,
    without_machine_paths,
)

PROJECT_MANIFEST_SCHEMA = "polyzymd-project-manifest/1"
#: Files freeze writes in the project folder and commits.
PROJECT_GENERATED = (MANIFEST, CITATION, ZENODO)
#: Files freeze writes in each study folder of a project.
STUDY_GENERATED = (MANIFEST, CHECKLIST, SUMMARY)


def _sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()


def _next_project_tag(root: Path) -> str:
    existing = (_git(root, "tag", "--list", "project-v*") or "").split()
    numbers = [int(t.split("v")[-1]) for t in existing if t.split("v")[-1].isdigit()]
    return f"project-v{max(numbers, default=0) + 1}"


def freeze_project(root: str | Path, *, tag: str | None = None) -> FreezeResult:
    """Freeze the project in ``root`` and every study it lists; see the module docstring.

    Raises
    ------
    ProtocolError
        When the project or a study file cannot be read, or the tag exists.
        Everything else is a warning in the result.
    """
    import polyzymd
    from polyzymd.analyses.project import Project
    from polyzymd.analyses.study_git import git_state
    from polyzymd.analyses.study_metadata import check_metadata, citation_cff, dump_cff, zenodo_json
    from polyzymd.analyses.study_upload_guide import deposit_readme, prepare_upload
    from polyzymd.citation import citation_line

    project = Project(root)
    root = project.root
    meta, warnings = check_metadata(project.protocol.metadata)
    state = git_state(root)
    if state is None:
        warnings.append("the project is not a git repository, so freeze cannot commit or tag it")
    elif state["inputs_uncommitted"]:
        warnings.append(
            "uncommitted inputs are not part of the tagged project: "
            + ", ".join(state["inputs_uncommitted"])
        )
    tag = tag or (_next_project_tag(root) if state else None)
    if state and tag and _git(root, "rev-parse", "--verify", "--quiet", f"refs/tags/{tag}"):
        raise ProtocolError(f"The tag {tag} already exists.", hint="Give another with --tag.")

    studies: dict[str, Any] = {}
    conditions: dict[str, Any] = {}
    for label in project.labels:
        study = project[label]
        result = freeze(study.root, publish=False)
        warnings.extend(f"{label}: {text}" for text in result.warnings)
        folder = study.root.relative_to(root)
        studies[label] = {
            "folder": str(folder),
            "manifest": {"path": f"{folder}/{MANIFEST}", "sha256": _sha256(study.root / MANIFEST)},
            "description": study.protocol.description,
        }
        for name, condition in result.manifest["conditions"].items():
            conditions[f"{label} / {name}"] = condition

    released = date.today().isoformat()
    version = tag or "unversioned"
    manifest: dict[str, Any] = {
        "schema": PROJECT_MANIFEST_SCHEMA,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tag": tag,
        "git": {
            "commit": state["commit"] if state else None,
            "inputs_uncommitted": state["inputs_uncommitted"] if state else None,
        },
        "project_file": {
            "path": project.protocol.path.name,
            "sha256": _sha256(project.protocol.path),
        },
        "versions": _versions(),
        "metadata": meta,
        "studies": studies,
        "conditions": conditions,
        "analyses": {run: project.runs_in(run) for run in project.protocol.analyses},
        "trajectory_deposits": meta["related"]["trajectories"],
        "cite": {"polyzymd": citation_line()},
        "warnings": warnings,
    }
    (root / MANIFEST).write_text(json.dumps(manifest, indent=1) + "\n")
    method = (
        f"Analysed with PolyzyMD {polyzymd.__version__}: {len(studies)} studies, one per "
        f"protein ({', '.join(studies)}), each against its own control; the replicate is "
        "the sampling unit."
    )
    (root / CITATION).write_text(
        dump_cff(citation_cff(meta, version=version, released=released, commit=None))
    )
    (root / ZENODO).write_text(
        json.dumps(zenodo_json(meta, version=version, released=released, method=method), indent=2)
        + "\n"
    )
    gitignore = root / ".gitignore"
    lines = gitignore.read_text().splitlines() if gitignore.exists() else []
    for entry in (f"{DEPOSIT}/", "logs/", "data.local.yaml"):
        if entry not in lines:
            lines.append(entry)
    gitignore.write_text("\n".join([*lines, ""]))

    commit = None
    if state and tag:
        paths = [p for p in (*PROJECT_GENERATED, ".gitignore") if (root / p).exists()]
        for label in project.labels:
            folder = project[label].root.relative_to(root)
            paths += [
                str(folder / p)
                for p in (*STUDY_GENERATED, "results")
                if (root / folder / p).exists()
            ]
        results = root / "results"
        if results.exists():
            paths.append("results")
        _git(root, "add", "--", *paths)
        if _git(root, "commit", "--quiet", "-m", f"Freeze project as {tag}", "--", *paths) is None:
            warnings.append("git could not commit the frozen files (is user.name set?)")
        elif (
            _git(root, "tag", "-a", tag, "-m", f"Project frozen by PolyzyMD {polyzymd.__version__}")
            is None
        ):
            warnings.append(f"git could not create the tag {tag}")
        else:
            commit = (_git(root, "rev-parse", "HEAD") or "").strip() or None

    deposit = root / DEPOSIT
    if deposit.exists():
        shutil.rmtree(deposit)
    copy = deposit / "study"
    copy.mkdir(parents=True)
    if commit:
        archive = subprocess.run(
            ["git", "-C", str(root), "archive", "--format=tar", tag], capture_output=True
        )
        subprocess.run(["tar", "-x", "-C", str(copy)], input=archive.stdout, check=False)
    else:
        for path in root.rglob("*"):
            relative = path.relative_to(root)
            if (
                path.is_file()
                and relative.parts[0] not in (DEPOSIT, ".git")
                and DEPOSIT not in relative.parts
                and path.name != "data.local.yaml"
                and "logs" not in relative.parts
            ):
                target = copy / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)
    for label in project.labels:
        study = project[label]
        for config in study.protocol.conditions.values():
            if config.is_relative_to(root) and (copy / config.relative_to(root)).is_file():
                deposited = copy / config.relative_to(root)
                deposited.write_text(without_machine_paths(deposited.read_text()))
        for part in ("engine_inputs", "final_frames"):
            source = study.root / DEPOSIT / part
            if source.is_dir():
                shutil.copytree(source, deposit / part / label, dirs_exist_ok=True)
    for name in (MANIFEST, CITATION, ZENODO):
        shutil.copy2(root / name, deposit / name)
    shutil.copy2(
        Path(__file__).parent / "schemas" / MANIFEST_SCHEMA_FILE, deposit / MANIFEST_SCHEMA_FILE
    )
    readme = deposit_readme(
        study_name=root.name, tag=tag if commit else None, meta=meta, analyses={}, root=root
    )
    (deposit / "README.md").write_text(readme + _studies_section(project))
    prepared = prepare_upload(
        deposit,
        study_name=root.name,
        tag=tag if commit else None,
        manifest=manifest,
        zenodo=json.loads((root / ZENODO).read_text()),
        warnings=warnings,
    )
    return FreezeResult(
        root,
        tag if commit else None,
        commit,
        deposit,
        manifest,
        warnings,
        guide=prepared["guide"],
        upload=prepared["upload"],
    )


def _studies_section(project: Any) -> str:
    """The deposit README's account of each study: its protein and what each analysis found."""
    lines = [
        "## Studies",
        "",
        "This deposit is a project: `project.yaml` lists one study per protein, and the "
        "analyses each runs. Reproduce every study's analyses with "
        "`polyzymd analyze --project .`, and read every result with "
        '`pz.Project(".").results(run)`.',
        "",
    ]
    for label in project.labels:
        study = project[label]
        lines += [f"### {label}", "", str(study.protocol.description or ""), ""]
        for run in study.protocol.analyses:
            report = study.protocol.results_dir(run) / "report.json"
            try:
                verdicts = json.loads(report.read_text()).get("verdict", [])
            except (OSError, ValueError):
                verdicts = ["no stored report"]
            lines.append(f"- **{run}:** " + " ".join(verdicts))
        lines.append("")
    return "\n".join(lines)
