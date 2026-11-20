"""The per-user list of PolyzyMD projects and lone studies.

``polyzymd project init`` and ``study init`` add the new folder to
``projects.yaml`` in ``$POLYZYMD_CONFIG_DIR`` (default
``$XDG_CONFIG_HOME/polyzymd``, else ``~/.config/polyzymd``), so an AI agent
started in another folder can find the project with ``polyzymd project list``.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from datetime import date
from pathlib import Path

import yaml

from polyzymd.analyses.exceptions import ProtocolError

REGISTRY_FILE = "projects.yaml"


def config_dir() -> Path:
    """Return the folder of PolyzyMD's per-user settings."""
    override = os.environ.get("POLYZYMD_CONFIG_DIR")
    if override:
        return Path(override).expanduser()
    xdg = os.environ.get("XDG_CONFIG_HOME", "")
    base = Path(xdg) if os.path.isabs(xdg) else Path.home() / ".config"
    return base / "polyzymd"


def registry_path() -> Path:
    """Return the path of the list of projects."""
    return config_dir() / REGISTRY_FILE


@dataclass(frozen=True)
class Entry:
    """One registered project or lone study.

    Attributes
    ----------
    name : str
        A unique name, the folder name unless another entry has it.
    path : str
        The absolute path of the folder.
    kind : str
        ``project`` or ``study``.
    studies : list of str
        The study labels of a project; empty for a study.
    created : str
        The date it was registered, ISO format.
    """

    name: str
    path: str
    kind: str
    studies: list[str] = field(default_factory=list)
    created: str = ""

    def exists(self) -> bool:
        """Whether the folder is there now: one stat, so a mount that is down only shows as missing."""
        return os.path.isdir(self.path)


def load() -> list[Entry]:
    """Return every registered entry, or none when there is no list yet."""
    path = registry_path()
    if not path.is_file():
        return []
    try:
        data = yaml.safe_load(path.read_text()) or {}
        return [Entry(**item) for item in data.get("projects") or []]
    except (OSError, yaml.YAMLError, TypeError, AttributeError) as exc:
        raise ProtocolError(
            f"Cannot read the list of projects {path}: {exc}",
            hint="Fix or delete the file, then register your projects again.",
        ) from exc


def _save(entries: list[Entry]) -> None:
    path = registry_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(".yaml.tmp")
    partial.write_text(yaml.safe_dump({"projects": [asdict(e) for e in entries]}, sort_keys=False))
    os.replace(partial, path)


def _describe(root: Path) -> tuple[str, list[str]]:
    """Return the kind of the folder ``root`` and its study labels."""
    from polyzymd.analyses.project_file import PROJECT_FILE

    if (root / PROJECT_FILE).is_file():
        data = yaml.safe_load((root / PROJECT_FILE).read_text()) or {}
        return "project", [str(label) for label in data.get("studies") or {}]
    return "study", []


def register(root: str | Path) -> tuple[Entry, bool]:
    """Add the project or study folder ``root`` to the list, or refresh its entry.

    The name is the folder name, or ``<name>-2``, ``<name>-3`` and so on when
    another folder has it. Returns the entry and whether it is new.
    """
    root = Path(root).expanduser().resolve()
    kind, studies = _describe(root)
    entries = load()
    for index, entry in enumerate(entries):
        if entry.path == str(root):
            entries[index] = Entry(entry.name, entry.path, kind, studies, entry.created)
            _save(entries)
            return entries[index], False
    taken = {entry.name for entry in entries}
    name, number = root.name, 2
    while name in taken:
        name, number = f"{root.name}-{number}", number + 1
    entry = Entry(name, str(root), kind, studies, date.today().isoformat())
    _save([*entries, entry])
    return entry, True


def unregister(name: str) -> Entry:
    """Remove the entry ``name`` from the list; the folder is not touched.

    Raises
    ------
    ProtocolError
        If no entry has that name.
    """
    entries = load()
    for entry in entries:
        if entry.name == name:
            _save([other for other in entries if other is not entry])
            return entry
    raise ProtocolError(
        f"No registered project is named {name!r}.",
        hint="polyzymd project list gives the names.",
    )


def register_quietly(root: str | Path) -> str | None:
    """Register ``root`` and return a warning instead of raising when that fails."""
    try:
        register(root)
    except Exception as exc:  # the list is a convenience; it must never fail an init
        return (
            f"could not add {root} to {registry_path()}: {exc}; "
            "run polyzymd project register later"
        )
    return None
