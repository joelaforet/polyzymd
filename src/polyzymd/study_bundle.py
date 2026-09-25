"""Package a study for publication, and check a package someone received.

:func:`export_study` writes one zip of a study folder with a
``bundle_manifest.json`` inside it. The study is what its comparisons
reference: a condition folder no ``comparison.yaml`` lists, such as a test or a
rerun, is left out. Trajectories and checkpoints are left out too, since they
are archived separately; the manifest lists every trajectory the published
results were computed from, with its size and content fingerprint, so a reader
can check the files they download.

:func:`verify_study` checks an unpacked study against its manifest: every
packaged file by SHA-256, and every listed trajectory that is present by size
and fingerprint.
"""

from __future__ import annotations

import hashlib
import json
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import yaml

from polyzymd.config.study import STUDY_FILE, StudyConfig

MANIFEST_NAME = "bundle_manifest.json"
MANIFEST_VERSION = 1

#: Files never packaged: trajectories and checkpoints are archived elsewhere.
EXCLUDED_SUFFIXES = frozenset({".dcd", ".xtc", ".trr", ".nc", ".chk", ".cpt"})

#: Folders never packaged: caches, environments and scheduler logs.
EXCLUDED_DIRS = frozenset(
    {"__pycache__", ".pytest_cache", ".git", ".pixi", ".polymer_cache", "slurm_logs"}
)


@dataclass
class ExportPlan:
    """What a study export contains and what it leaves out."""

    root: Path
    files: list[Path] = field(default_factory=list)
    conditions: list[Path] = field(default_factory=list)
    unreferenced: list[Path] = field(default_factory=list)
    outside: list[str] = field(default_factory=list)
    trajectories: list[dict[str, Any]] = field(default_factory=list)

    @property
    def size_bytes(self) -> int:
        """Total size of the packaged files."""
        return sum(path.stat().st_size for path in self.files)


@dataclass
class VerifyReport:
    """Result of checking an unpacked study against its manifest."""

    checked_files: int = 0
    missing_files: list[str] = field(default_factory=list)
    changed_files: list[str] = field(default_factory=list)
    trajectories_ok: list[str] = field(default_factory=list)
    trajectories_absent: list[str] = field(default_factory=list)
    trajectories_changed: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """No packaged file is missing or changed and no trajectory changed."""
        return not (self.missing_files or self.changed_files or self.trajectories_changed)


def plan_export(root: Path | str) -> ExportPlan:
    """Decide which files of a study go into its package.

    Parameters
    ----------
    root : Path or str
        Study root, the folder holding ``study.yaml``.

    Returns
    -------
    ExportPlan
        Files to package, the condition folders they include, condition
        folders left out because no comparison lists them, conditions that
        live outside the study, and the trajectories the results were computed
        from.
    """
    root = Path(root).resolve()
    if not (root / STUDY_FILE).is_file():
        raise FileNotFoundError(f"{root} is not a study: it holds no {STUDY_FILE}")
    plan = ExportPlan(root=root)
    labels_by_comparison: dict[Path, dict[str, Path]] = {}
    for comparison_yaml in sorted(_walk(root, "comparison.yaml")):
        labels = _conditions_of(comparison_yaml)
        labels_by_comparison[comparison_yaml.parent] = labels
        for label, config in labels.items():
            if _inside(config, root):
                plan.conditions.append(config.parent)
            else:
                plan.outside.append(f"{label}: {config}")
    plan.conditions = sorted(set(plan.conditions))
    plan.unreferenced = sorted(
        {config.parent for config in _walk(root, "config.yaml")} - set(plan.conditions)
    )
    plan.files = [
        path
        for path in sorted(_walk(root))
        if path.name != MANIFEST_NAME
        and not any(_inside(path, left_out) for left_out in plan.unreferenced)
    ]
    for directory, labels in labels_by_comparison.items():
        plan.trajectories += _trajectories(root, directory, labels)
    return plan


def export_study(root: Path | str, output: Path | str) -> ExportPlan:
    """Write a study's package: a zip of its files plus ``bundle_manifest.json``.

    Every path in the zip starts with the study folder's name, so unpacking
    it recreates that folder with the manifest at its root.
    """
    from polyzymd import __version__
    from polyzymd.analyses.identity import _git_commit

    plan = plan_export(root)
    output = Path(output).resolve()
    plan.files = [path for path in plan.files if path.resolve() != output]
    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "study": StudyConfig.from_yaml(plan.root / STUDY_FILE).name,
        "polyzymd_version": __version__,
        "polyzymd_git_commit": _git_commit(),
        "files": [
            {
                "path": path.relative_to(plan.root).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in plan.files
        ],
        "trajectories": plan.trajectories,
        "left_out_conditions": [
            path.relative_to(plan.root).as_posix() for path in plan.unreferenced
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    prefix = plan.root.name
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in plan.files:
            archive.write(path, f"{prefix}/{path.relative_to(plan.root).as_posix()}")
        archive.writestr(f"{prefix}/{MANIFEST_NAME}", json.dumps(manifest, indent=2) + "\n")
    return plan


def verify_study(
    root: Path | str,
    *,
    working_dir: Callable[[Path, int], Path | None] | None = None,
) -> VerifyReport:
    """Check an unpacked study against the manifest at its root.

    Parameters
    ----------
    root : Path or str
        Folder holding ``bundle_manifest.json``.
    working_dir : callable, optional
        ``(config_path, replicate) -> directory`` giving where a replicate's
        trajectories live. Defaults to reading each condition's config.

    Returns
    -------
    VerifyReport
        Missing and changed files, and which trajectories are present and
        match, are absent, or changed.
    """
    root = Path(root).resolve()
    manifest = json.loads((root / MANIFEST_NAME).read_text())
    report = VerifyReport()
    for entry in manifest["files"]:
        path = root / entry["path"]
        report.checked_files += 1
        if not path.is_file():
            report.missing_files.append(entry["path"])
        elif _sha256(path) != entry["sha256"]:
            report.changed_files.append(entry["path"])
    locate = working_dir or _config_working_dir
    for entry in manifest.get("trajectories", []):
        name = f"{entry['condition']} replicate {entry['replicate']}: {entry['relative_path']}"
        directory = locate(root / entry["config"], int(entry["replicate"]))
        path = directory / entry["relative_path"] if directory is not None else None
        if path is None or not path.is_file():
            report.trajectories_absent.append(name)
            continue
        if path.stat().st_size != entry.get("size_bytes") or (
            entry.get("fingerprint") and _fingerprint(path) != entry["fingerprint"]
        ):
            report.trajectories_changed.append(name)
        else:
            report.trajectories_ok.append(name)
    return report


def _walk(root: Path, name: str | None = None) -> list[Path]:
    """Files under ``root`` outside excluded folders, optionally of one name."""
    found = []
    for path in root.rglob(name or "*"):
        if not path.is_file() or path.suffix.lower() in EXCLUDED_SUFFIXES:
            continue
        if any(part in EXCLUDED_DIRS for part in path.relative_to(root).parts[:-1]):
            continue
        found.append(path)
    return found


def _conditions_of(comparison_yaml: Path) -> dict[str, Path]:
    """Condition label to config path, as a comparison file lists them."""
    data = yaml.safe_load(comparison_yaml.read_text()) or {}
    conditions: dict[str, Path] = {}
    for item in data.get("conditions") or []:
        if not isinstance(item, dict) or "config" not in item:
            continue
        config = Path(str(item["config"])).expanduser()
        if not config.is_absolute():
            config = comparison_yaml.parent / config
        conditions[str(item.get("label", config.parent.name))] = config.resolve()
    return conditions


def _trajectories(
    root: Path, comparison_dir: Path, labels: dict[str, Path]
) -> list[dict[str, Any]]:
    """Trajectories the replicate results of one comparison were computed from."""
    seen: set[tuple[str, int, str]] = set()
    listed: list[dict[str, Any]] = []
    for result in sorted(comparison_dir.glob("analysis/*/*/run_*/result.json")):
        try:
            artifact = json.loads(result.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        label = artifact.get("condition_label")
        config = labels.get(label)
        if config is None or not _inside(config, root):
            continue
        replicate = int(artifact.get("replicate", 0))
        identity = (artifact.get("provenance") or {}).get("identity") or {}
        for item in identity.get("inputs") or []:
            relative = item.get("relative_path")
            if not relative or Path(relative).suffix.lower() not in EXCLUDED_SUFFIXES:
                continue
            key = (config.relative_to(root).as_posix(), replicate, relative)
            if key in seen:
                continue
            seen.add(key)
            listed.append(
                {
                    "condition": label,
                    "config": key[0],
                    "replicate": replicate,
                    "relative_path": relative,
                    "size_bytes": item.get("size_bytes"),
                    "fingerprint": item.get("fingerprint"),
                }
            )
    return listed


def _config_working_dir(config: Path, replicate: int) -> Path | None:
    """Working directory of one replicate, read from its condition's config."""
    from polyzymd.config.loader import load_config

    try:
        return Path(load_config(config).get_working_directory(replicate))
    except Exception:
        return None


def _inside(path: Path, directory: Path) -> bool:
    try:
        path.resolve().relative_to(directory.resolve())
    except ValueError:
        return False
    return True


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def _fingerprint(path: Path) -> str:
    from polyzymd.analyses.identity import file_fingerprint

    return file_fingerprint(path)
