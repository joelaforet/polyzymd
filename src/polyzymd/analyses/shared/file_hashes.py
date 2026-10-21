"""SHA-256 of input files, computed once per file and location and kept in a cache.

Stored results and frozen studies identify trajectory and topology files by
their content, so a copied or downloaded study keeps its results. Hashing a
large trajectory takes about a second per gigabyte, so :func:`file_sha256`
does it once: it uses the hash the run recorded in ``progress.json`` (read
through the simulation engine,
:meth:`~polyzymd.engines.base.SimulationEngine.recorded_trajectory_hashes`)
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
import tempfile
from pathlib import Path


def cache_dir() -> Path:
    """Return the folder of cached file hashes."""
    base = os.environ.get("POLYZYMD_CACHE_DIR")
    root = Path(base).expanduser() if base else Path.home() / ".cache" / "polyzymd"
    return root / "hashes"


def _entry(path: Path) -> Path:
    return cache_dir() / f"{hashlib.sha1(str(path).encode()).hexdigest()}.json"


def file_sha256(
    path: str | Path, known: tuple[str, int] | None = None, use_cache: bool = True
) -> str:
    """Return the SHA-256 of the file at ``path``.

    ``known`` is a ``(sha256, size)`` the run recorded for this file, used
    when its size equals the file's; otherwise a cached hash valid for the
    file's current size and modification time is used, and failing that the
    file is read and the hash cached. With ``use_cache=False`` the file is
    always read and nothing is cached: integrity checks use this, because a
    file replaced with the same size and modification time keeps its cache
    entry.
    """
    file = Path(path).resolve()
    stat = file.stat()
    if known is not None and known[0] and known[1] == stat.st_size:
        return known[0]
    entry = _entry(file)
    if use_cache:
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
    if not use_cache:
        return value
    try:
        entry.parent.mkdir(parents=True, exist_ok=True)
        # A unique temporary name, so nodes sharing the cache never collide.
        with tempfile.NamedTemporaryFile(
            "w", dir=entry.parent, prefix=f".{entry.name}.", suffix=".tmp", delete=False
        ) as handle:
            json.dump(
                {
                    "path": str(file),
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                    "sha256": value,
                },
                handle,
            )
        os.replace(handle.name, entry)
    except OSError:
        pass  # an unwritable cache only costs hashing again
    return value
