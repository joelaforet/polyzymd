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
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

from polyzymd.analyses.shared.file_hashes import file_sha256
from polyzymd.analyses.study_freeze import (
    CHECKLIST,
    CITATION,
    DEPOSIT,
    MANIFEST,
    SUMMARY,
    ZENODO,
    FreezeResult,
    _commit_and_tag,
    _copy_frozen_folder,
    _drop_tag,
    _finish_deposit,
    _git_preflight,
    _listed_files,
    _versions,
    _write_citation,
    _write_gitignore,
    deposited_entry,
    freeze,
    group_warnings,
    left_out_files,
)

#: Value of the ``schema`` key of a project's ``manifest.json``.
PROJECT_MANIFEST_SCHEMA = "polyzymd-project-manifest/1"
#: Files freeze writes in the project folder and commits.
PROJECT_GENERATED = (MANIFEST, CITATION, ZENODO)
#: Files freeze writes in each study folder of a project.
STUDY_GENERATED = (MANIFEST, CHECKLIST, SUMMARY)


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
    from polyzymd.analyses.study_metadata import check_metadata
    from polyzymd.analyses.study_upload_guide import deposit_readme
    from polyzymd.citation import citation_line

    project = Project(root)
    root = project.root
    meta, warnings = check_metadata(project.protocol.metadata, what="project")
    state, tag, found = _git_preflight(root, tag, "project")
    warnings += found

    studies: dict[str, Any] = {}
    conditions: dict[str, Any] = {}
    for label in project.labels:
        study = project[label]
        result = freeze(study.root, publish=False)
        warnings.extend(f"{label}: {text}" for text in result.warnings)
        folder = study.root.relative_to(root)
        studies[label] = {
            "folder": str(folder),
            "manifest": {"path": f"{folder}/{MANIFEST}", "sha256": file_sha256(study.root / MANIFEST)},
            "description": study.protocol.description,
        }
        for name, condition in result.manifest["conditions"].items():
            conditions[f"{label} / {name}"] = condition

    # The same warning for several studies or conditions is one line naming them.
    names = [*project.labels, *(c.replace(" / ", ": ") for c in conditions)]
    warnings[:] = group_warnings(warnings, names)
    released = date.today().isoformat()
    version = tag or "unversioned"
    _write_gitignore(root)
    method = (
        f"Analysed with PolyzyMD {polyzymd.__version__}: {len(studies)} studies "
        f"({', '.join(studies)}), each against its own control; the replicate is "
        "the sampling unit."
    )
    _write_citation(root, meta, version=version, released=released, commit=None, method=method)
    manifest: dict[str, Any] = {
        "schema": PROJECT_MANIFEST_SCHEMA,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tag": tag,
        "git": {
            "parent_commit": state["commit"] if state else None,
            "inputs_uncommitted": state["inputs_uncommitted"] if state else None,
        },
        "project_file": {
            "path": project.protocol.path.name,
            "sha256": file_sha256(project.protocol.path),
        },
        "versions": _versions(root),
        "metadata": meta,
        "studies": studies,
        "conditions": conditions,
        "analyses": {run: project.runs_in(run) for run in project.protocol.analyses},
        # The project's own files (project.yaml, shared analyses/ and stats/
        # code, figures and what they wrote) and its citation files; each
        # study's are in its manifest.
        "files": _project_files(project, state),
        "trajectory_deposits": meta["related"]["trajectories"],
        "cite": {"polyzymd": citation_line()},
        "warnings": warnings,
    }
    (root / MANIFEST).write_text(json.dumps(manifest, indent=1) + "\n")

    commit = None
    if state and tag:
        paths = [p for p in (*PROJECT_GENERATED, ".gitignore", "results") if (root / p).exists()]
        for label in project.labels:
            folder = project[label].root.relative_to(root)
            paths += [
                str(folder / p)
                for p in (*STUDY_GENERATED, "results")
                if (root / folder / p).exists()
            ]
        commit = _commit_and_tag(root, paths, tag, "project", warnings)
        if commit is None:
            _drop_tag(root, manifest, meta, released, method)

    deposit = root / DEPOSIT
    if deposit.exists():
        shutil.rmtree(deposit)
    configs = [
        str(config.relative_to(root))
        for label in project.labels
        for config in project[label].protocol.conditions.values()
        if config.is_relative_to(root)
    ]
    left_out = left_out_files(root, None)
    if left_out:
        warnings.append(left_out)
    _copy_frozen_folder(root, deposit, tag, commit, _listed_files(root, None), configs)
    for label in project.labels:
        for part in ("engine_inputs", "final_frames"):
            source = project[label].root / DEPOSIT / part
            if source.is_dir():
                shutil.copytree(source, deposit / part / label, dirs_exist_ok=True)
    readme = deposit_readme(
        study_name=root.name,
        tag=tag if commit else None,
        meta=meta,
        analyses={},
        root=root,
        project=True,
    )
    result = _finish_deposit(
        root, deposit, tag, commit, manifest, readme + _studies_section(project), warnings
    )
    result.git_failed = bool(state and tag and commit is None)
    return result


def _project_files(project: Any, state: dict | None) -> dict[str, dict[str, Any]]:
    """Return the size and SHA-256 of every project file outside the study folders, as deposited.

    Tracked files and untracked ``results/`` files in a git repository,
    otherwise every file, and the citation files freeze writes;
    ``deposit/``, ``logs/``, ``data.local.yaml`` and ``manifest.json`` are
    left out. A condition config is described without its machine paths
    (:func:`~polyzymd.analyses.study_freeze.deposited_entry`).
    """
    root = project.root
    studies = {project[label].root.relative_to(root).parts[0] for label in project.labels}
    configs = {
        config.resolve()
        for label in project.labels
        for config in project[label].protocol.conditions.values()
    }
    names = [*_listed_files(root, state), CITATION, ZENODO]
    return {
        name: deposited_entry(root / name, root, configs)
        for name in dict.fromkeys(names)
        if Path(name).parts[0] not in studies and name != MANIFEST and (root / name).is_file()
    }


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
        "This deposit is a project: `project.yaml` lists its studies, and the "
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
