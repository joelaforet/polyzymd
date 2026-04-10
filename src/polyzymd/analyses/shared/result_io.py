"""Result discovery helpers for analysis plugin comparison outputs."""

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
    analysis_type: str | None = None,
    fallback_subdir: str = "comparison",
    fallback_filenames: list[str] | None = None,
    log: logging.Logger | None = None,
) -> T | None:
    """Locate and load a saved comparison result JSON.

    Implements the standard two-phase discovery strategy used by all plugin ``plot()`` methods:

    1. Primary search in ``data["__meta__"]`` using canonical comparison metadata
       (exact ``comparison_result_path`` first, then ``comparison_dir``)
    2. Legacy search in ``data["__meta__"]["results_dir"]``
    3. Fallback search via per-condition path navigation into ``fallback_subdir``

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
    analysis_type : str | None, optional
        Analysis type used to check the canonical
        ``comparison/<analysis_type>/result.json`` location first.
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
        comparison_result_path = meta.get("comparison_result_path")
        if comparison_result_path is not None:
            loaded = _try_load_exact(Path(comparison_result_path), loader, log)
            if loaded is not None:
                return loaded
            log.debug(
                f"Could not load comparison result from {comparison_result_path} - trying fallback"
            )

        comparison_dir = meta.get("comparison_dir")
        if comparison_dir is not None:
            loaded = _try_load_from_comparison_root(
                Path(comparison_dir),
                analysis_type,
                glob_patterns,
                loader,
                log,
                fallback_filenames=fallback_filenames,
            )
            if loaded is not None:
                return loaded
            log.debug(f"No matching result JSON in {comparison_dir} - trying legacy paths")

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

                loaded = _try_load_from_comparison_root(
                    candidate,
                    analysis_type,
                    glob_patterns,
                    loader,
                    log,
                    fallback_filenames=fallback_filenames,
                )
                if loaded is not None:
                    return loaded

        condition = condition_data.get("condition")
        if condition is not None:
            # Support both Condition.config_path (plugin system) and
            # ConditionConfig.config (comparison config) for backward compat.
            config_path = getattr(condition, "config_path", None) or getattr(
                condition, "config", None
            )
            if config_path is not None:
                config_path = Path(config_path)
                for parent in [config_path.parent, config_path.parent.parent]:
                    candidate = parent / fallback_subdir
                    if candidate not in searched and candidate.is_dir():
                        searched.add(candidate)
                        loaded = _try_load_from_comparison_root(
                            candidate,
                            analysis_type,
                            glob_patterns,
                            loader,
                            log,
                            fallback_filenames=fallback_filenames,
                        )
                        if loaded is not None:
                            return loaded

    return None


def canonical_comparison_result_path(results_dir: Path, analysis_type: str) -> Path:
    """Return the canonical comparison result path for one analysis."""
    return results_dir / analysis_type / "result.json"


def _try_load_exact(path: Path, loader: Callable[[Path], T], log: logging.Logger) -> T | None:
    """Try loading one exact file path."""
    if not path.exists():
        return None

    try:
        result = loader(path)
        log.debug(f"Loaded result from {path}")
        return result
    except Exception as exc:  # noqa: BLE001
        log.debug(f"Could not load {path}: {exc}")
        return None


def _try_load_from_comparison_root(
    comparison_dir: Path,
    analysis_type: str | None,
    glob_patterns: list[str],
    loader: Callable[[Path], T],
    log: logging.Logger,
    *,
    fallback_filenames: list[str] | None = None,
) -> T | None:
    """Try canonical comparison paths before legacy filename searches."""
    if not comparison_dir.is_dir():
        return None

    if analysis_type is not None:
        loaded = _try_load_exact(
            canonical_comparison_result_path(comparison_dir, analysis_type),
            loader,
            log,
        )
        if loaded is not None:
            return loaded

    if fallback_filenames is not None:
        for filename in fallback_filenames:
            loaded = _try_load_exact(comparison_dir / filename, loader, log)
            if loaded is not None:
                return loaded

    return _try_load_from_dir(comparison_dir, glob_patterns, loader, log)


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
