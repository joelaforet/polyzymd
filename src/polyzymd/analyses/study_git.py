"""The git state of a study folder, recorded with every analysis run from it.

A study folder can be a git repository (``polyzymd study init`` makes one).
:func:`git_state` reads its commit and the files that differ from it, so a
report says which version of the study produced it. Stored results are
reused by content, never by commit, so committing changes nothing that is
stored, and uncommitted changes give an analysis run a warning, never a
refusal. ``freeze`` refuses them (:mod:`polyzymd.analyses.study_freeze`).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

#: Folders and files of a study whose changes are outputs, not inputs, of an analysis.
OUTPUTS = ("results/", "data.local.yaml", "logs/", "deposit/")


def is_output(path: str) -> bool:
    """Return whether ``path`` is an analysis output or machine file, at any depth.

    ``results/``, ``logs/`` and ``deposit/`` folders and ``data.local.yaml``
    count wherever they are, so a project's ``<study>/results/`` is an output too.
    So do compiled Python (``__pycache__/``, ``*.pyc``) and one machine's job
    and log folders (``slurm/``, ``slurm_logs/``), which are never inputs.
    """
    parts = Path(path).parts
    return bool(parts) and (
        parts[-1] == "data.local.yaml"
        or parts[-1].endswith(".pyc")
        or any(f"{part}/" in OUTPUTS for part in parts[:-1])
        or any(part in ("__pycache__", "slurm", "slurm_logs") for part in parts[:-1])
    )


def _git(root: Path, *arguments: str) -> str | None:
    """Run git in ``root`` and return its output, or ``None`` if it fails."""
    try:
        result = subprocess.run(
            ["git", "-C", str(root), *arguments],
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout if result.returncode == 0 else None


def git_state(root: str | Path) -> dict[str, Any] | None:
    """Return the git commit of the study folder ``root`` and its uncommitted files.

    Returns
    -------
    dict or None
        ``None`` when git is not installed or ``root`` is not inside a git
        repository. Otherwise ``commit`` (the full SHA of ``HEAD``, or
        ``None`` before the first commit), ``uncommitted`` (paths relative to
        ``root`` that are modified, staged or untracked, ignored files left
        out) and ``inputs_uncommitted``, those of ``uncommitted`` that are not
        analysis outputs (``results/``, ``logs/``, ``deposit/``) or the machine's
        ``data.local.yaml``.
    """
    root = Path(root).resolve()
    if shutil.which("git") is None or _git(root, "rev-parse", "--show-toplevel") is None:
        return None
    head = _git(root, "rev-parse", "--verify", "--quiet", "HEAD")
    status = _git(root, "status", "--porcelain", "--untracked-files=all", "--", ".") or ""
    top = Path((_git(root, "rev-parse", "--show-toplevel") or str(root)).strip())
    uncommitted = []
    for line in status.splitlines():
        path = line[3:].split(" -> ")[-1].strip().strip('"')
        try:
            uncommitted.append(str((top / path).resolve().relative_to(root)))
        except ValueError:
            uncommitted.append(path)
    inputs = [path for path in uncommitted if not is_output(path)]
    return {
        "commit": head.strip() if head else None,
        "uncommitted": sorted(uncommitted),
        "inputs_uncommitted": sorted(inputs),
    }


def describe(state: dict[str, Any] | None) -> str:
    """Return one line describing ``state`` for ``polyzymd study check``."""
    if state is None:
        return "git: not a repository (polyzymd study init makes one)"
    commit = (state["commit"] or "no commit yet")[:12]
    paths = state["inputs_uncommitted"]
    if paths:
        more = f" and {len(paths) - 5} more" if len(paths) > 5 else ""
        return (
            f"git: commit {commit}; {len(paths)} uncommitted inputs: {', '.join(paths[:5])}{more}"
        )
    return f"git: commit {commit}; inputs committed"


def init_repository(root: Path, message: str) -> str | None:
    """Make ``root`` a git repository and commit everything in it.

    Returns the commit's SHA, or ``None`` when git is missing or the commit
    fails (for example with no ``user.name`` configured); the folder is then
    left as it is, with a repository if ``git init`` worked.
    """
    if shutil.which("git") is None or _git(root, "init", "--quiet") is None:
        return None
    if _git(root, "add", "--all") is None or _git(root, "commit", "--quiet", "-m", message) is None:
        return None
    head = _git(root, "rev-parse", "HEAD")
    return head.strip() if head else None
