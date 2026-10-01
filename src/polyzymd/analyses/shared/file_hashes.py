"""SHA-256 of input files, computed once per file and location and kept in a cache.

Stored results and frozen studies identify trajectory and topology files by
their content, so a copied or downloaded study keeps its results. Hashing a
large trajectory takes about a second per gigabyte, so :func:`file_sha256`
does it once: it uses the hash the simulation recorded in ``progress.json``
when the size matches, and otherwise the cache, one small JSON file per
hashed path in ``$POLYZYMD_CACHE_DIR/hashes`` (default
``~/.cache/polyzymd/hashes``), valid while the file's size and modification
time are unchanged. Each entry is written atomically, so parallel jobs can
share the cache.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


def cache_dir() -> Path:
    """Return the folder of cached file hashes."""
    base = os.environ.get("POLYZYMD_CACHE_DIR")
    root = Path(base).expanduser() if base else Path.home() / ".cache" / "polyzymd"
    return root / "hashes"


def _entry(path: Path) -> Path:
    return cache_dir() / f"{hashlib.sha1(str(path).encode()).hexdigest()}.json"


def file_sha256(path: str | Path, known: tuple[str, int] | None = None) -> str:
    """Return the SHA-256 of the file at ``path``.

    ``known`` is a ``(sha256, size)`` the run recorded for this file, used
    when its size equals the file's; otherwise a cached hash valid for the
    file's current size and modification time is used, and failing that the
    file is read and the hash cached.
    """
    file = Path(path).resolve()
    stat = file.stat()
    if known is not None and known[0] and known[1] == stat.st_size:
        return known[0]
    entry = _entry(file)
    try:
        cached = json.loads(entry.read_text())
        if cached["size"] == stat.st_size and cached["mtime_ns"] == stat.st_mtime_ns:
            return cached["sha256"]
    except (OSError, ValueError, KeyError):
        pass
    digest = hashlib.sha256()
    with file.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    value = digest.hexdigest()
    try:
        entry.parent.mkdir(parents=True, exist_ok=True)
        temporary = entry.with_suffix(f".{os.getpid()}.tmp")
        temporary.write_text(
            json.dumps(
                {
                    "path": str(file),
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                    "sha256": value,
                }
            )
        )
        temporary.replace(entry)
    except OSError:
        pass  # an unwritable cache only costs hashing again
    return value


def recorded_segment_hashes(working_dir: str | Path) -> dict[str, tuple[str, int]]:
    """Return the trajectory hashes ``progress.json`` recorded, by trajectory file name."""
    from polyzymd.simulation.progress import load_progress

    try:
        progress = load_progress(working_dir)
    except Exception:  # noqa: BLE001 - an unreadable progress file records nothing
        return {}
    if progress is None:
        return {}
    return {
        f"production_{s.index}_trajectory.dcd": (s.trajectory_sha256, s.trajectory_bytes)
        for s in progress.segments
        if s.trajectory_sha256 and s.trajectory_bytes is not None
    }
