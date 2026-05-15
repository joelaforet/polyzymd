"""Result I/O and path helpers for analysis plugins."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from pydantic import BaseModel

if TYPE_CHECKING:
    from polyzymd.analyses._contexts import PlotContext
    from polyzymd.config.schema import SimulationConfig

logger = logging.getLogger("polyzymd.analyses")


def load_aggregated_result(analysis: Any, aggregated_dir: Path) -> Any:
    """Load the aggregated result from disk.

    Parameters
    ----------
    analysis : Any
        Analysis instance requesting the result.
    aggregated_dir : Path
        Directory containing aggregated result files.

    Returns
    -------
    Any
        Loaded result, or ``None`` if no file exists.
    """
    if not aggregated_dir.exists():
        return None
    canonical = analysis.aggregate_result_path(aggregated_dir)
    if canonical.exists():
        return analysis._deserialize_result(canonical)
    json_files = sorted(aggregated_dir.glob("*.json"), key=lambda path: path.stat().st_mtime)
    if not json_files:
        return None
    chosen = json_files[-1]
    logger.warning(
        "%s: canonical result.json not found in %s — falling back to %s (%d JSON file(s) present)",
        analysis.name,
        aggregated_dir,
        chosen.name,
        len(json_files),
    )
    return analysis._deserialize_result(chosen)


def deserialize_result(analysis: Any, path: Path) -> Any:
    """Deserialize an aggregated JSON result.

    Parameters
    ----------
    analysis : Any
        Analysis instance that owns the result class configuration.
    path : Path
        Path to the JSON result file.

    Returns
    -------
    Any
        Deserialized result.
    """
    result_cls = type(analysis).AggregatedResultClass
    if result_cls is not None:
        if hasattr(result_cls, "load"):
            return result_cls.load(path)
        if hasattr(result_cls, "model_validate_json"):
            return result_cls.model_validate_json(path.read_text())
    return json.loads(path.read_text())


def deserialize_replicate_result(analysis: Any, path: Path) -> Any:
    """Deserialize a per-replicate JSON result.

    Parameters
    ----------
    analysis : Any
        Analysis instance that owns the result class configuration.
    path : Path
        Path to the JSON result file.

    Returns
    -------
    Any
        Deserialized result.
    """
    result_cls = type(analysis).ReplicateResultClass
    if result_cls is not None:
        if hasattr(result_cls, "load"):
            return result_cls.load(path)
        if hasattr(result_cls, "model_validate_json"):
            return result_cls.model_validate_json(path.read_text())
    loaded = json.loads(path.read_text())
    if isinstance(loaded, dict) and loaded.get("artifact_type") == "replicate":
        from polyzymd.analyses.mda.artifacts import ReplicateArtifact

        return ReplicateArtifact.model_validate(loaded)
    return loaded


def load_replicate_result(analysis: Any, run_dir: Path) -> Any | None:
    """Load a per-replicate result from a run directory.

    Parameters
    ----------
    analysis : Any
        Analysis instance requesting the result.
    run_dir : Path
        Replicate run directory.

    Returns
    -------
    Any | None
        Deserialized result, or ``None`` if no canonical result exists.
    """
    if not run_dir.exists():
        return None
    result_path = analysis.replicate_result_path(run_dir)
    if not result_path.exists():
        return None
    return analysis._deserialize_replicate_result(result_path)


def check_cache(
    analysis: Any,
    result_cls: type,
    cache_path: Path,
    *,
    recompute: bool,
    sim_config: SimulationConfig | None = None,
    settings: BaseModel | None = None,
) -> Any | None:
    """Load a cached result from disk when it is valid.

    Parameters
    ----------
    analysis : Any
        Analysis instance requesting the cache lookup.
    result_cls : type
        Result class with a ``load()`` method.
    cache_path : Path
        Path to the cached JSON result file.
    recompute : bool
        If ``True``, skip cache loading.
    sim_config : SimulationConfig, optional
        Simulation configuration used for hash validation.
    settings : BaseModel, optional
        Settings model used for fingerprint validation.

    Returns
    -------
    Any | None
        Cached result on hit, otherwise ``None``.
    """
    if recompute or not cache_path.exists():
        return None

    if not hasattr(result_cls, "load"):
        raise TypeError(
            f"_check_cache requires result_cls to have a .load() method, got {result_cls!r}."
        )

    logger.info("Loading cached %s result from %s", analysis.name, cache_path)
    result = result_cls.load(cache_path)

    if sim_config is not None and hasattr(result, "config_hash"):
        from polyzymd.analyses.shared.config_hash import validate_config_hash

        stored_hash = getattr(result, "config_hash", None)
        if stored_hash not in (None, "", "unknown") and not validate_config_hash(
            str(stored_hash),
            sim_config,
        ):
            return None

    if settings is not None:
        from polyzymd.analyses.shared.config_hash import (
            extract_settings_fingerprint_from_path,
            validate_settings_fingerprint,
        )

        stored_fingerprint = getattr(result, "settings_fingerprint", None)
        if stored_fingerprint is None:
            stored_fingerprint = getattr(result, "settings_fp", None)
        if stored_fingerprint is None:
            metadata = getattr(result, "metadata", None)
            if isinstance(metadata, dict):
                stored_fingerprint = metadata.get("settings_fingerprint")
                if stored_fingerprint is None:
                    stored_fingerprint = metadata.get("settings_fp")
        if stored_fingerprint is None:
            stored_fingerprint = extract_settings_fingerprint_from_path(cache_path)

        if not validate_settings_fingerprint(
            stored_fingerprint,
            settings,
            source=cache_path,
        ):
            return None

    return result


def replicate_result_path(output_dir: Path) -> Path:
    """Return the canonical per-replicate cache path."""
    return output_dir / "result.json"


def aggregate_result_path(output_dir: Path) -> Path:
    """Return the canonical aggregated cache path."""
    return output_dir / "result.json"


def comparison_result_path(results_dir: Path) -> Path:
    """Return the canonical comparison cache path."""
    return results_dir / "result.json"


def figures_output_dir(analysis: Any, figures_root: Path) -> Path:
    """Return the analysis-specific figure directory."""
    return figures_root / analysis.name


def save_result(result: Any, path: Path) -> Path:
    """Save a result object to disk using the framework contract.

    Parameters
    ----------
    result : Any
        Result object or JSON-serializable payload.
    path : Path
        Destination path.

    Returns
    -------
    Path
        Saved path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(result, "save"):
        return result.save(path)
    if hasattr(result, "model_dump_json"):
        path.write_text(result.model_dump_json(indent=2))
        return path
    path.write_text(json.dumps(result, indent=2))
    return path


def resolve_output_dir(analysis: Any, analysis_root: Path, condition_label: str) -> Path:
    """Build the analysis output directory for a condition.

    Parameters
    ----------
    analysis : Any
        Analysis instance.
    analysis_root : Path
        Root analysis directory.
    condition_label : str
        Condition label to sanitize for filesystem use.

    Returns
    -------
    Path
        Analysis output directory.
    """
    from polyzymd.analyses.shared.paths import sanitize_label

    return analysis_root / sanitize_label(condition_label) / analysis.name


def format_replicate_range(replicates: Sequence[int]) -> str:
    """Format replicate numbers as a compact string.

    Parameters
    ----------
    replicates : Sequence[int]
        Replicate numbers.

    Returns
    -------
    str
        Compact replicate range string.
    """
    if not replicates:
        return "no_replicates"
    reps = sorted(set(replicates))
    if reps == list(range(reps[0], reps[-1] + 1)):
        return f"reps{reps[0]}-{reps[-1]}"
    return "reps" + "_".join(map(str, reps))


def build_plot_data(
    ctx: PlotContext,
    *,
    include_replicates: bool = False,
) -> tuple[dict[str, Any], list[str]]:
    """Build the data and label structures consumed by plotters.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.
    include_replicates : bool, optional
        If ``True``, include replicate numbers for each condition.

    Returns
    -------
    tuple[dict[str, Any], list[str]]
        Plot data mapping and ordered condition labels.
    """
    data: dict[str, Any] = {}
    labels: list[str] = []

    for condition in ctx.conditions:
        label = condition.label
        if label == "__meta__":
            logger.warning("Condition label '__meta__' conflicts with reserved key — skipping.")
            continue
        analysis_dir = ctx.analysis_dirs.get(label)
        if analysis_dir is not None:
            entry: dict[str, Any] = {
                "analysis_dir": analysis_dir,
                "aggregated_dir": analysis_dir / "aggregated",
            }
            if include_replicates:
                entry["replicates"] = list(condition.replicates)
            data[label] = entry
            labels.append(label)

    data["__meta__"] = {
        "results_dir": ctx.results_dir,
        "comparison_result_path": ctx.results_dir / "result.json",
        "comparison_dir": ctx.results_dir,
        "settings": ctx.settings,
        "control_label": ctx.control_label,
    }
    return data, labels
