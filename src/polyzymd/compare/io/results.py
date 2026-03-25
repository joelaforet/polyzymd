from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable, Sequence, TypeVar

logger = logging.getLogger(__name__)

T = TypeVar("T")


def find_comparison_result(
    data: dict[str, Any],
    labels: Sequence[str],
    glob_patterns: list[str],
    loader: Callable[[Path], T],
    *,
    fallback_subdir: str = "comparison",
    fallback_filenames: list[str] | None = None,
    log: logging.Logger | None = None,
) -> T | None:
    """Locate and load a saved comparison result JSON.

    Implements the standard two-phase discovery strategy used by all plotters:

    1. Primary search in ``data["__meta__"]["results_dir"]`` for files matching
       ``glob_patterns`` and load the most recently modified one
    2. Fallback search via per-condition path navigation into ``fallback_subdir``

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of ``condition_label -> condition data``, plus ``"__meta__"`` entry.
    labels : Sequence[str]
        Condition labels.
    glob_patterns : list[str]
        Glob patterns to search for result JSON files.
    loader : Callable[[Path], T]
        Loader function that reads and returns a result object.
    fallback_subdir : str, optional
        Subdirectory name to search in fallback mode, by default ``"comparison"``.
    fallback_filenames : list[str] | None, optional
        Exact fallback filenames to try before globs, by default None.
    log : logging.Logger | None, optional
        Logger to use, by default module logger.

    Returns
    -------
    T | None
        Loaded result object, or None when no loadable result is found.
    """
    if log is None:
        log = logger

    meta = data.get("__meta__")
    if meta is not None:
        results_dir = meta.get("results_dir")
        if results_dir is not None:
            loaded = _try_load_from_dir(Path(results_dir), glob_patterns, loader, log)
            if loaded is not None:
                return loaded
            log.debug(f"No matching result JSON in {results_dir} - trying fallback")

    searched: set[Path] = set()

    for label in labels:
        condition_data = data.get(label)
        if condition_data is None:
            continue

        analysis_dir = condition_data.get("analysis_dir")
        if analysis_dir is not None:
            project_root = Path(analysis_dir).parent.parent
            candidate = project_root / fallback_subdir
            if candidate not in searched and candidate.is_dir():
                searched.add(candidate)

                if fallback_filenames is not None:
                    for filename in fallback_filenames:
                        file_path = candidate / filename
                        if file_path.exists():
                            try:
                                result = loader(file_path)
                                log.debug(f"Loaded result from {file_path}")
                                return result
                            except Exception as exc:  # noqa: BLE001
                                log.debug(f"Could not load {file_path}: {exc}")

                loaded = _try_load_from_dir(candidate, glob_patterns, loader, log)
                if loaded is not None:
                    return loaded

        condition = condition_data.get("condition")
        if condition is not None:
            config_path = getattr(condition, "config", None)
            if config_path is not None:
                config_path = Path(config_path)
                for parent in [config_path.parent, config_path.parent.parent]:
                    candidate = parent / fallback_subdir
                    if candidate not in searched and candidate.is_dir():
                        searched.add(candidate)
                        loaded = _try_load_from_dir(candidate, glob_patterns, loader, log)
                        if loaded is not None:
                            return loaded

    return None


def _try_load_from_dir(
    directory: Path,
    glob_patterns: list[str],
    loader: Callable[[Path], T],
    log: logging.Logger,
) -> T | None:
    """Try loading the newest matching result file from a directory.

    Parameters
    ----------
    directory : Path
        Directory to search.
    glob_patterns : list[str]
        Glob patterns used to collect candidate files.
    loader : Callable[[Path], T]
        Loader function for candidate files.
    log : logging.Logger
        Logger used for debug messages.

    Returns
    -------
    T | None
        Loaded result object, or None if no loadable candidate exists.
    """
    if not directory.is_dir():
        return None

    files: list[Path] = []
    for pattern in glob_patterns:
        files.extend(directory.glob(pattern))

    if not files:
        return None

    result_file = max(files, key=lambda path: path.stat().st_mtime)
    try:
        result = loader(result_file)
        log.debug(f"Loaded result from {result_file}")
        return result
    except Exception as exc:  # noqa: BLE001
        log.debug(f"Could not load {result_file}: {exc}")

    return None
