"""Result I/O and path helpers for analysis plugins."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

if TYPE_CHECKING:
    from polyzymd.analyses._framework.contexts import PlotContext

logger = logging.getLogger("polyzymd.analyses")

_STRICT_CANONICAL_AGGREGATE_ANALYSES = frozenset(
    {
        "catalytic_triad",
        "contacts",
        "distances",
        "hydrogen_bonds",
        "rg",
        "rmsd",
        "rmsf",
        "sasa",
        "secondary_structure",
    }
)


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
    canonical = analysis.aggregate_result_path(aggregated_dir)
    if canonical.exists():
        return analysis._deserialize_result(canonical)
    return None


def deserialize_result(analysis: Any, path: Path) -> Any:
    """Deserialize an aggregated artifact or JSON result.

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
    text: str | None = None
    if path.exists():
        text = path.read_text()
        artifact = _try_load_mda_condition_artifact(analysis, path.parent, path, text)
        if artifact is not None:
            return artifact
    if getattr(analysis, "name", None) in _STRICT_CANONICAL_AGGREGATE_ANALYSES:
        raise ValueError(
            f"{analysis.name}: aggregate result at {path} is not a canonical MDA "
            "ConditionArtifact. Remove stale legacy JSON caches and recompute this analysis."
        )
    result_cls = type(analysis).AggregatedResultClass
    if result_cls is not None:
        if hasattr(result_cls, "load"):
            return result_cls.load(path)
        if hasattr(result_cls, "model_validate_json") and text is not None:
            return result_cls.model_validate_json(text)
    if text is None:
        text = path.read_text()
    return json.loads(text)


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
    artifact = _try_load_mda_replicate_artifact(analysis, path.parent, path)
    if artifact is not None:
        return artifact
    return json.loads(path.read_text())


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
    return _load_replicate_result_from_path(analysis, run_dir, result_path)


def _load_replicate_result_from_path(analysis: Any, run_dir: Path, result_path: Path) -> Any:
    """Load a replicate result with MDA artifact-store routing.

    Parameters
    ----------
    analysis : Any
        Analysis instance requesting the result.
    run_dir : Path
        Replicate output directory containing the result.
    result_path : Path
        Canonical replicate result path.

    Returns
    -------
    Any
        Loaded plugin result or MDA replicate artifact.
    """

    result_cls = type(analysis).ReplicateResultClass
    if result_cls is not None:
        return analysis._deserialize_replicate_result(result_path)
    artifact = _try_load_mda_replicate_artifact(analysis, run_dir, result_path)
    if artifact is not None:
        return artifact
    return analysis._deserialize_replicate_result(result_path)


def _try_load_mda_replicate_artifact(analysis: Any, run_dir: Path, result_path: Path) -> Any | None:
    """Try loading an MDA replicate artifact through ``ArtifactStore``.

    Parameters
    ----------
    analysis : Any
        Analysis instance requesting the artifact.
    run_dir : Path
        Replicate output directory used as the artifact-store root.
    result_path : Path
        Candidate artifact path.

    Returns
    -------
    Any or None
        Loaded ``ReplicateArtifact`` when the file is an MDA artifact, otherwise
        ``None`` so standard result deserialization can continue.
    """

    from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError

    store = ArtifactStore(run_dir)
    relative_path = _artifact_relative_path(run_dir, result_path)
    try:
        return store.read_replicate_result(relative_path)
    except ArtifactStoreError as exc:
        try:
            loaded = json.loads(result_path.read_text())
        except (OSError, json.JSONDecodeError) as parse_exc:
            raise _artifact_load_error(analysis, run_dir, result_path, exc) from parse_exc
        if isinstance(loaded, dict) and loaded.get("artifact_type") == "replicate":
            raise _artifact_load_error(analysis, run_dir, result_path, exc) from exc
    return None


def _try_load_mda_condition_artifact(
    analysis: Any,
    aggregated_dir: Path,
    result_path: Path,
    text: str | None = None,
) -> Any | None:
    """Try loading an MDA condition artifact through ``ArtifactStore``.

    Parameters
    ----------
    analysis : Any
        Analysis instance requesting the artifact.
    aggregated_dir : Path
        Aggregate output directory used as the artifact-store root.
    result_path : Path
        Candidate artifact path.
    text : str or None, optional
        Already-read JSON text for artifact-type probing.

    Returns
    -------
    Any or None
        Loaded ``ConditionArtifact`` when the file is an MDA condition artifact,
        otherwise ``None`` so standard result deserialization can continue.
    """

    try:
        loaded = json.loads(text if text is not None else result_path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(loaded, dict) or loaded.get("artifact_type") != "condition":
        return None

    from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError

    store = ArtifactStore(aggregated_dir)
    relative_path = _artifact_relative_path(aggregated_dir, result_path)
    try:
        return store.read_condition_result(relative_path)
    except ArtifactStoreError as exc:
        raise ArtifactStoreError(
            f"{analysis.name}: failed to load MDA condition artifact from {result_path}: {exc}"
        ) from exc


def _artifact_relative_path(run_dir: Path, result_path: Path) -> Path:
    """Return the store-relative artifact path when possible.

    Parameters
    ----------
    run_dir : Path
        Artifact-store root.
    result_path : Path
        Candidate artifact path.

    Returns
    -------
    Path
        Path relative to ``run_dir``.
    """

    try:
        return result_path.relative_to(run_dir)
    except ValueError:
        return Path(result_path.name)


def _artifact_load_error(
    analysis: Any, run_dir: Path, result_path: Path, exc: Exception
) -> Exception:
    """Build a contextual MDA artifact loading error.

    Parameters
    ----------
    analysis : Any
        Analysis instance requesting the artifact.
    run_dir : Path
        Replicate output directory.
    result_path : Path
        Artifact path that failed to load.
    exc : Exception
        Original artifact-store failure.

    Returns
    -------
    Exception
        Contextual artifact-store error.
    """

    from polyzymd.analyses.mda.store import ArtifactStoreError

    replicate = _replicate_id_from_run_dir(run_dir)
    return ArtifactStoreError(
        f"{analysis.name}: failed to load MDA replicate artifact for replicate={replicate} "
        f"from {result_path}: {exc}"
    )


def _replicate_id_from_run_dir(run_dir: Path) -> str:
    """Extract a replicate ID from a canonical ``run_N`` directory name.

    Parameters
    ----------
    run_dir : Path
        Replicate output directory.

    Returns
    -------
    str
        Replicate ID text or ``"unknown"`` when the directory is non-standard.
    """

    name = run_dir.name
    if name.startswith("run_") and name.removeprefix("run_").isdigit():
        return name.removeprefix("run_")
    return "unknown"


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
    from polyzymd.analyses.mda.artifacts import ComparisonArtifact, ConditionArtifact

    if isinstance(result, ConditionArtifact):
        from polyzymd.analyses.mda.store import ArtifactStore

        return ArtifactStore(path.parent).write_condition_result(
            result, path.relative_to(path.parent)
        )
    if isinstance(result, ComparisonArtifact):
        from polyzymd.analyses.mda.store import ArtifactStore

        return ArtifactStore(path.parent).write_comparison_result(
            result, path.relative_to(path.parent)
        )
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
