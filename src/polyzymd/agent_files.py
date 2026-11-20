"""The guidance for AI coding agents that PolyzyMD writes into a project.

``polyzymd project init`` and ``polyzymd study init`` write ``AGENTS.md``,
``CLAUDE.md`` and the user skills (``.claude/skills/<name>/SKILL.md`` for
Claude Code, ``.agents/skills/<name>/SKILL.md`` for Codex and other agents)
into the new folder, and ``polyzymd agent sync`` rewrites them to the
installed version.

The skills have one source, ``.claude/skills/`` of the PolyzyMD repository,
which Claude Code also reads when PolyzyMD itself is developed. The wheel
carries a copy at ``polyzymd/templates/agents/skills/`` (``force-include`` in
``pyproject.toml``); an editable install reads the repository's copy. Every
skill there is shipped except :data:`DEVELOPER_SKILLS`.

Each written file holds a marker line with the PolyzyMD version that wrote it
and the SHA-256 of the rest of the file. A file whose content no longer
matches its marker was edited, and sync leaves it alone unless forced.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from polyzymd.analyses.exceptions import ProtocolError

#: Skills under ``.claude/skills/`` that are for developing PolyzyMD, not using it.
DEVELOPER_SKILLS = frozenset({"livecoms-check"})

#: Where the skills are written in a project: Claude Code, then Codex and other agents.
SKILL_FOLDERS = (".claude/skills", ".agents/skills")

#: The guides written at the top of a project or study. Like the skills, they
#: are neither an input of an analysis nor deposited by freeze.
GUIDES = ("AGENTS.md", "CLAUDE.md")

#: Where the wheel holds the skills, inside the ``polyzymd`` package.
PACKAGED_SKILLS = "templates/agents/skills"

MARKER = re.compile(
    r"^<!-- polyzymd-managed: version=(?P<version>\S+) sha256=(?P<sha>[0-9a-f]{64}) .*-->\n",
    re.MULTILINE,
)
MARKER_NOTE = "polyzymd agent sync rewrites this file unless you edit it"


def skills_source() -> Path:
    """Return the folder holding one folder per skill, as shipped or in the repository.

    The installed wheel's copy comes first; an editable or pixi development
    install has none and reads ``.claude/skills/`` of the checkout it runs from.

    Raises
    ------
    ProtocolError
        If neither exists.
    """
    packaged = files("polyzymd").joinpath(PACKAGED_SKILLS)
    if packaged.is_dir():
        return Path(str(packaged))
    import polyzymd

    checkout = Path(polyzymd.__file__).resolve().parents[2] / ".claude" / "skills"
    if checkout.is_dir():
        return checkout
    raise ProtocolError(
        "This PolyzyMD installation has no agent skills.",
        hint="Reinstall PolyzyMD from a wheel or from a full checkout of its repository.",
    )


def user_skills() -> dict[str, str]:
    """Return each user skill's ``SKILL.md`` text by skill name, without the developer skills."""
    source = skills_source()
    return {
        folder.name: (folder / "SKILL.md").read_text()
        for folder in sorted(source.iterdir())
        if folder.name not in DEVELOPER_SKILLS and (folder / "SKILL.md").is_file()
    }


def managed_files() -> dict[str, str]:
    """Return the text of every file PolyzyMD manages in a project, without markers, by path."""
    shipped = files("polyzymd.templates").joinpath("agents")
    managed = {name: shipped.joinpath(name).read_text() for name in GUIDES}
    for name, text in user_skills().items():
        for folder in SKILL_FOLDERS:
            managed[f"{folder}/{name}/SKILL.md"] = text
    return managed


def _sha(body: str) -> str:
    return hashlib.sha256(body.encode()).hexdigest()


def with_marker(body: str, version: str) -> str:
    """Return ``body`` with its marker line: after the YAML front matter, or first."""
    line = f"<!-- polyzymd-managed: version={version} sha256={_sha(body)} {MARKER_NOTE} -->\n"
    if body.startswith("---\n"):
        end = body.find("\n---\n", 4)
        if end != -1:
            cut = end + len("\n---\n")
            return body[:cut] + line + body[cut:]
    return line + body


def read_marker(text: str) -> tuple[str, str, str] | None:
    """Return the version, recorded hash and body of a managed file, or ``None`` without a marker."""
    match = MARKER.search(text)
    if match is None:
        return None
    return match["version"], match["sha"], text[: match.start()] + text[match.end() :]


@dataclass(frozen=True)
class SyncAction:
    """What sync does, or would do, to one managed file.

    Attributes
    ----------
    path : str
        The file, relative to the project or study folder.
    action : str
        ``created``, ``updated``, ``unchanged`` or ``skipped``.
    """

    path: str
    action: str

    def line(self) -> str:
        """Return the one-line agent output for this file."""
        if self.action == "skipped":
            return f"skipped: {self.path} (edited; --force to overwrite)"
        return f"{self.action}: {self.path}"


def sync(root: str | Path, *, force: bool = False, dry_run: bool = False) -> list[SyncAction]:
    """Write the installed version's agent files into the folder ``root``.

    A missing file is created. A file whose body still matches the hash in
    its marker is rewritten when the installed version's body differs, and
    left as it is when the body is the same. A file without a marker, or
    whose body no longer matches its hash, was edited and is skipped unless
    ``force``. With ``dry_run`` nothing is written. Nothing is committed.
    """
    from polyzymd import __version__

    root = Path(root)
    actions = []
    for relative, body in managed_files().items():
        target = root / relative
        if not target.exists():
            action = "created"
        else:
            marker = read_marker(target.read_text())
            edited = marker is None or marker[1] != _sha(marker[2])
            if edited and not force:
                action = "skipped"
            elif not edited and marker[2] == body:
                action = "unchanged"
            else:
                action = "updated"
        if action in ("created", "updated") and not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(with_marker(body, __version__))
        actions.append(SyncAction(relative, action))
    return actions


def find_root(path: str | Path) -> Path:
    """Return the project folder holding ``path``, or else the study folder holding it.

    ``path`` and each of its parents are searched for ``project.yaml``
    first, so a study inside a project gives the project.

    Raises
    ------
    ProtocolError
        If ``path`` is in neither a project nor a study.
    """
    from polyzymd.analyses.project_file import PROJECT_FILE
    from polyzymd.analyses.study_file import STUDY_FILE

    start = Path(path).expanduser().resolve()
    folders = [start, *start.parents]
    for name in (PROJECT_FILE, STUDY_FILE):
        for folder in folders:
            if (folder / name).is_file():
                return folder
    raise ProtocolError(
        f"{start} is not inside a PolyzyMD project or study.",
        hint="Run this inside a folder made by polyzymd project init or study init, or give its path.",
    )
