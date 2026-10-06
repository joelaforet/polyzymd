"""Freeze a project for publication: every study, one manifest, one citation, one deposit.

The module's public function is :func:`freeze_project`. It freezes each
study of a project with :func:`~polyzymd.analyses.study_freeze.freeze`
(``publish=False``), writes the project's ``manifest.json``,
``CITATION.cff`` and ``.zenodo.json``, commits and tags the project in git
when it is a repository, and lays out one ``deposit/`` folder for the whole
project, in the same form as a study's deposit.
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
    EXCLUDE_MACHINE_FILES,
    MANIFEST,
    MANIFEST_SCHEMA_FILE,
    SUMMARY,
    ZENODO,
    FreezeResult,
    _git,
    _versions,
    committable,
    drop_machine_files,
    freeze,
    is_machine_file,
    stats_warnings,
    without_machine_paths,
)

#: Value of the ``schema`` key of a project's ``manifest.json``.
PROJECT_MANIFEST_SCHEMA = "polyzymd-project-manifest/1"
#: Files freeze writes in the project folder and commits.
PROJECT_GENERATED = (MANIFEST, CITATION, ZENODO)
#: Files freeze writes in each study folder of a project.
STUDY_GENERATED = (MANIFEST, CHECKLIST, SUMMARY)


def _sha256(path: Path) -> str:
    """Return the SHA-256 hex digest of the bytes of ``path``."""
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()


def _next_project_tag(root: Path) -> str:
    """Return ``project-v<n>``, with ``n`` one more than the highest existing such tag."""
    existing = (_git(root, "tag", "--list", "project-v*") or "").split()
    numbers = [int(t.split("v")[-1]) for t in existing if t.split("v")[-1].isdigit()]
    return f"project-v{max(numbers, default=0) + 1}"


def freeze_project(root: str | Path, *, tag: str | None = None) -> FreezeResult:
    """Freeze the project in ``root`` and every study it lists, for one publication.

    The function does the following, in order:

    - Freezes each study with :func:`~polyzymd.analyses.study_freeze.freeze`
      and ``publish=False``, which writes the study's ``manifest.json``,
      ``md_checklist.yaml`` and ``system_summary.csv`` and copies its engine
      inputs and final frames to its ``deposit/``. Each study's warnings are
      added to the result, prefixed with its label.
    - Writes the project's ``manifest.json``, which lists each study's
      manifest with its SHA-256, every condition as ``<study> / <condition>``,
      the studies that run each analysis, and the software versions; and
      writes ``CITATION.cff`` and ``.zenodo.json`` from the ``metadata:`` of
      ``project.yaml``. Adds ``deposit/``, ``logs/`` and ``data.local.yaml``
      to ``.gitignore``.
    - When the project is a git repository, commits the generated files,
      ``.gitignore`` and the ``results/`` of the project and of every study,
      and tags the commit (``project-v<n>`` unless ``tag`` is given).
    - Rebuilds ``deposit/`` at the project root: under ``deposit/study`` the
      tagged commit (from ``git archive``) or, without a commit, a copy of
      the project folder without ``.git``, ``deposit/``, ``logs/`` and
      ``data.local.yaml``, with machine-specific directories removed from
      the condition configs inside the project; the engine inputs and final frames of each
      study under ``deposit/<part>/<label>``; the project's manifest,
      citation files and manifest schema; a ``README.md`` with a section per
      study giving each analysis's stored verdicts; and the upload folder
      and ``UPLOAD.md`` guide from
      :func:`~polyzymd.analyses.study_upload_guide.prepare_upload`.

    Nothing is uploaded.

    Parameters
    ----------
    root : str or Path
        The project folder or its ``project.yaml``.
    tag : str, optional
        Git tag for the frozen project. Defaults to the next ``project-v<n>``.

    Returns
    -------
    FreezeResult
        The project root, the tag and commit (both ``None`` when nothing was
        committed), the deposit folder, the project manifest, the warnings,
        and the paths of the upload guide and upload folder.

    Raises
    ------
    ProtocolError
        If ``project.yaml`` or a study's files cannot be read, or if the tag
        already exists. Other problems (missing metadata, uncommitted inputs,
        a failed commit or tag) are added to the result's warnings.
    """
    import polyzymd
    from polyzymd.analyses.project import Project
    from polyzymd.analyses.study_git import git_state
    from polyzymd.analyses.study_metadata import check_metadata, citation_cff, dump_cff, zenodo_json
    from polyzymd.analyses.study_upload_guide import deposit_readme, prepare_upload
    from polyzymd.citation import citation_line

    project = Project(root)
    root = project.root
    meta, warnings = check_metadata(project.protocol.metadata, what="project")
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

    if project.protocol.stats is not None:
        warnings.extend(stats_warnings(root, project.protocol.stats, project=True))

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
        # The project's own files (project.yaml, shared analyses/ and stats/
        # code, figures, the stats plan's output); each study's are in its manifest.
        "files": _project_files(project, state),
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
        paths += EXCLUDE_MACHINE_FILES
        _git(root, "add", "--", *paths)
        paths = committable(root, paths)
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
                and not is_machine_file(str(relative))
            ):
                target = copy / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)
    drop_machine_files(copy)
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
        study_name=root.name,
        tag=tag if commit else None,
        meta=meta,
        analyses={},
        root=root,
        project=True,
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


def _project_files(project: Any, state: dict | None) -> dict[str, dict[str, Any]]:
    """Return the size and SHA-256 of every project file outside the study folders.

    Tracked files and untracked ``results/`` files in a git repository,
    otherwise every file; ``deposit/``, ``logs/``, ``data.local.yaml`` and the
    files freeze writes are left out.
    """
    root = project.root
    studies = {project[label].root.relative_to(root).parts[0] for label in project.labels}
    if state is not None:
        listed = (_git(root, "ls-files") or "").splitlines() + (
            _git(root, "ls-files", "--others", "--exclude-standard", "--", "results") or ""
        ).splitlines()
    else:
        listed = [
            str(p.relative_to(root))
            for p in root.rglob("*")
            if p.is_file() and not is_machine_file(str(p.relative_to(root)))
        ]
    files = {}
    for name in sorted(set(listed)):
        parts = Path(name).parts
        if (
            not parts
            or parts[0] in studies
            or parts[0] in (DEPOSIT, ".git", "logs")
            or is_machine_file(name)
            or name in (*PROJECT_GENERATED, "data.local.yaml")
            or not (root / name).is_file()
        ):
            continue
        path = root / name
        files[name] = {"size": path.stat().st_size, "sha256": _sha256(path)}
    return files


def _studies_section(project: Any) -> str:
    """Return the deposit README's ``## Studies`` section: each study's description and verdicts.

    For every study it gives the study's description and, for each of its
    analyses, the ``verdict`` lines of the stored ``report.json``, led by
    ``PARTIAL REPORT (...)`` when the report is partial, or ``no stored
    report`` when that file is missing or unreadable.
    """
    from polyzymd.analyses.study_upload_guide import report_summary

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
            lines.append(
                f"- **{run}:** " + report_summary(study.protocol.results_dir(run) / "report.json")
            )
        lines.append("")
    return "\n".join(lines)
