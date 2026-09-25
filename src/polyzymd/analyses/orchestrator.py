"""Run an analysis over the replicates and conditions of a comparison.

One module owns the whole pipeline:

1. :func:`run_replicate_once` loads a replicate through
   :func:`polyzymd.analyses.loading.open_replicate`, calls the plugin's
   ``compute()``, reduces each observable to its replicate value and writes
   ``run_<N>/result.json``. A cached result is reused only when its identity
   block still matches; see :mod:`polyzymd.analyses.identity`.
2. :func:`run_analysis` does that for every replicate of a condition and
   aggregates them into ``aggregated/result.json``.
   :func:`aggregate_condition_from_disk` aggregates results a worker wrote.
3. :func:`finalize_comparison_from_disk` checks every aggregate, tests the
   observables across conditions under one Benjamini-Hochberg family, writes
   ``comparison/<name>/result.json`` and draws the figures.
4. :func:`run_comparison` chains the three for one analysis, and
   :func:`run_all_comparisons` for several.

The replicate is the sampling unit throughout; see
:mod:`polyzymd.analyses.contract` for the statistics.

On disk, for a comparison file in ``<dir>``::

    <dir>/analysis/<condition>/<name>/run_<N>/result.json
    <dir>/analysis/<condition>/<name>/aggregated/result.json
    <dir>/comparison/<name>/result.json
    <plot_settings.output_dir>/<name>/*.png
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from pydantic import BaseModel

from polyzymd.analyses import identity, loading
from polyzymd.analyses.base import (
    Analysis,
    ComparisonContext,
    Condition,
    PlotContext,
    _measurement_warnings,
    _unpack,
)
from polyzymd.analyses.contract import ObservableEstimate, aggregate_observables, reduce_replicate
from polyzymd.analyses.exceptions import (
    AggregateValidationError,
    AggregationError,
    AnalysisError,
    ComparisonError,
    PlotError,
    PluginContractError,
    ReplicateError,
    ReplicateSkippedError,
    StaleCacheError,
)
from polyzymd.analyses.mda.artifacts import ComparisonArtifact, ConditionArtifact, ReplicateArtifact
from polyzymd.analyses.mda.frame_selection import frame_selection_payload
from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError

if TYPE_CHECKING:
    from polyzymd.config.comparison import ComparisonConfig

logger = logging.getLogger("polyzymd.analyses")

RESULT_FILE = "result.json"
_MANY_TASKS_THRESHOLD = 10
_PLOT_ONLY_EXPECTED_FAILURES = (PlotError, ArtifactStoreError, OSError, ValueError)


# ---------------------------------------------------------------------------
# Replicates
# ---------------------------------------------------------------------------


def run_replicate_once(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    output_dir: Path,
    replicate: int,
    recompute: bool,
) -> ReplicateArtifact:
    """Compute one replicate, or reuse its cached result, and write it.

    Parameters
    ----------
    analysis : Analysis
        Analysis to run.
    condition : Condition
        Condition the replicate belongs to.
    settings : BaseModel
        Plugin settings.
    equilibration : str
        Window discarded from the start of the trajectory, for example
        ``"10ns"``.
    output_dir : Path
        Replicate directory, ``.../<name>/run_<N>``.
    replicate : int
        One-indexed replicate number.
    recompute : bool
        Ignore a cached result.

    Returns
    -------
    ReplicateArtifact
        The reduced observables with their identity block.

    Raises
    ------
    FileNotFoundError, OSError, ReplicateSkippedError
        When the replicate cannot be read; :func:`run_analysis` skips it.
    PluginContractError
        When the plugin breaks the contract. Never wrapped.
    ReplicateError
        For any other failure while computing.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    if not recompute:
        cached = _reusable(analysis, condition, settings, equilibration, output_dir, replicate)
        if cached is not None:
            return cached
    try:
        artifact, arrays = _compute_replicate(
            analysis, condition, settings, equilibration, replicate
        )
    except (FileNotFoundError, OSError, ReplicateSkippedError, PluginContractError):
        raise
    except Exception as exc:
        raise ReplicateError(
            f"{analysis.name}: compute stage failed for condition='{condition.label}' "
            f"replicate={replicate}: {type(exc).__name__}: {exc}"
        ) from exc
    store = ArtifactStore(output_dir)
    artifact.sidecars = [store.write_npz_sidecar(path, **data) for path, data in arrays.items()]
    store.write_replicate_result(artifact, RESULT_FILE)
    return artifact


def _compute_replicate(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    replicate: int,
) -> tuple[ReplicateArtifact, dict[str, dict[str, Any]]]:
    """Load the replicate and run ``compute()``.

    Returns the artifact without sidecars, and the arrays to write as sidecars
    keyed by their path in the replicate directory: the full per-frame series of
    every observable, and any extra table the plugin returned.
    """
    universe, frames, provenance = loading.open_replicate(
        condition.sim_config, replicate, equilibration
    )
    if frames.warning_message:
        logger.warning(
            "%s: %s [condition=%s, replicate=%d]",
            analysis.name,
            frames.warning_message,
            condition.label,
            replicate,
        )
    observables, extras = _unpack(analysis.plugin.compute(universe, frames, settings))
    estimates = reduce_replicate(observables)
    warnings = [frames.warning_message] if frames.warning_message else []
    warnings += [str(message) for message in provenance.get("warnings") or []]
    warnings += _measurement_warnings(observables)
    artifact = ReplicateArtifact(
        analysis_name=analysis.name,
        condition_label=condition.label,
        replicate=replicate,
        payload={"observables": [estimate.model_dump(mode="json") for estimate in estimates]},
        provenance={
            "source": "observable_contract",
            "identity": identity.replicate_identity(
                analysis.plugin,
                condition.sim_config,
                settings,
                equilibration,
                identity.describe_inputs(
                    identity.input_files(provenance),
                    _working_dir(condition, replicate),
                    fingerprint=True,
                ),
            ),
            "frame_selection": frame_selection_payload(frames),
            "universe_policy": {
                "condition_label": condition.label,
                "replicate": replicate,
                "provenance": provenance,
                "metadata": {"equilibration": equilibration},
            },
        },
        metadata={
            "result_kind": "observables",
            "settings_fingerprint": identity.settings_fingerprint(settings),
            "equilibration": str(equilibration),
        },
        warnings=list(dict.fromkeys(warnings)),
    )
    arrays = {"observables.npz": {item.name: item.values for item in observables}}
    arrays.update({f"sidecars/{stem}.npz": {stem: array} for stem, array in extras.items()})
    return artifact, arrays


def _reusable(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    output_dir: Path,
    replicate: int,
) -> ReplicateArtifact | None:
    """Return the cached replicate when its identity still matches, else ``None``."""
    result_path = output_dir / RESULT_FILE
    if not result_path.exists():
        return None
    try:
        cached = ArtifactStore(output_dir).read_replicate_result(RESULT_FILE)
    except (ArtifactStoreError, OSError, ValueError) as exc:
        logger.info(
            "%s: recomputing replicate %d because %s is unreadable: %s",
            analysis.name,
            replicate,
            result_path,
            exc,
        )
        return None
    current = identity.replicate_identity(
        analysis.plugin,
        condition.sim_config,
        settings,
        equilibration,
        identity.describe_inputs(
            identity.input_files(loading.replicate_provenance(condition.sim_config, replicate)),
            _working_dir(condition, replicate),
        ),
    )
    reason = identity.identity_mismatch(cached.provenance.get("identity"), current)
    if reason is not None:
        logger.info(
            "%s: recomputing replicate %d for '%s' because %s",
            analysis.name,
            replicate,
            condition.label,
            reason,
        )
        return None
    identity.warn_on_version_mismatch(cached.polyzymd_version, result_path)
    logger.info(
        "%s: reusing cached replicate %d for '%s' from %s",
        analysis.name,
        replicate,
        condition.label,
        result_path,
    )
    return cached


# ---------------------------------------------------------------------------
# Conditions
# ---------------------------------------------------------------------------


def run_analysis(
    analysis: Analysis,
    condition: Condition,
    settings: Any,
    equilibration: str = "0ns",
    output_dir: Path | None = None,
    recompute: bool = False,
) -> ConditionArtifact:
    """Compute every replicate of one condition and aggregate them.

    A replicate that cannot be read (a missing file, or a
    :class:`~polyzymd.analyses.exceptions.ReplicateSkippedError`) is skipped
    with a warning. Any other failure stops the condition.

    Returns
    -------
    ConditionArtifact
        The condition's aggregate, also written to ``aggregated/result.json``.

    Raises
    ------
    ValueError
        If no replicate succeeds.
    """
    if output_dir is None:
        output_dir = condition.config_path.parent / "analysis" / analysis.name
    logger.info(
        f"Running {analysis.name} for '{condition.label}' "
        f"(replicates {list(condition.replicates)})"
    )
    results: list[ReplicateArtifact] = []
    successful: list[int] = []
    failed: list[int] = []
    reasons: list[str] = []
    for replicate in condition.replicates:
        try:
            results.append(
                run_replicate_once(
                    analysis,
                    condition,
                    settings,
                    equilibration,
                    output_dir / f"run_{replicate}",
                    replicate,
                    recompute,
                )
            )
            successful.append(replicate)
        except (FileNotFoundError, OSError) as exc:
            logger.warning(
                "  Skipping %s rep %d: data not found — %s", condition.label, replicate, exc
            )
            failed.append(replicate)
            reasons.append(f"replicate {replicate}: {type(exc).__name__}: {exc}")
        except ReplicateSkippedError as exc:
            logger.warning("  Skipping %s rep %d: %s", condition.label, replicate, exc)
            failed.append(replicate)
            reasons.append(f"replicate {replicate}: {exc}")
    if not results:
        raise ValueError(
            f"{analysis.name}: condition '{condition.label}' has 0 successful replicates, "
            f"need at least 1.  Failed: {failed}{_format_reasons(reasons)}"
        )
    if failed:
        logger.warning(
            "  %s: %d replicate(s) failed %s, using %d of %d",
            condition.label,
            len(failed),
            failed,
            len(results),
            len(condition.replicates),
        )
    aggregated = _aggregate(
        analysis, condition, settings, equilibration, output_dir, results, successful, recompute
    )
    logger.info(f"  Aggregated {len(results)} replicates for '{condition.label}'")
    return aggregated


def aggregate_condition_from_disk(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    output_dir: Path,
    replicates: Sequence[int],
    recompute: bool = False,
) -> ConditionArtifact:
    """Aggregate replicate results already on disk.

    Nothing is recomputed here, so a replicate whose identity no longer
    matches is an error rather than a cache miss.

    Raises
    ------
    StaleCacheError
        If a replicate result no longer describes its inputs or this run.
    ValueError
        If no replicate result is on disk.
    """
    loaded: list[ReplicateArtifact] = []
    successful: list[int] = []
    missing: list[Path] = []
    warnings: list[str] = []
    for replicate in replicates:
        path = output_dir / f"run_{replicate}" / RESULT_FILE
        if not path.exists():
            missing.append(path)
            logger.warning(
                "%s: missing replicate result for '%s' rep %d at %s",
                analysis.name,
                condition.label,
                replicate,
                path,
            )
            continue
        try:
            artifact = ArtifactStore(path.parent).read_replicate_result(RESULT_FILE)
        except ArtifactStoreError as exc:
            raise ArtifactStoreError(
                f"{analysis.name}: failed to load replicate result for "
                f"condition='{condition.label}' replicate={replicate} from {path}: {exc}"
            ) from exc
        warnings += _require_fresh(
            analysis, condition, settings, equilibration, artifact, path, replicate
        )
        identity.warn_on_version_mismatch(artifact.polyzymd_version, path)
        loaded.append(artifact)
        successful.append(replicate)
    if not loaded:
        shown = ", ".join(str(path) for path in missing[:5])
        more = f", ... {len(missing) - 5} more" if len(missing) > 5 else ""
        raise ValueError(
            f"{analysis.name}: condition '{condition.label}' has 0 replicate result(s) on "
            f"disk, need at least 1. Expected missing replicate output path(s): {shown}{more}"
        )
    return _aggregate(
        analysis,
        condition,
        settings,
        equilibration,
        output_dir,
        loaded,
        successful,
        recompute,
        warnings=warnings,
    )


def _require_fresh(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    artifact: ReplicateArtifact,
    path: Path,
    replicate: int,
) -> list[str]:
    """Refuse a replicate result that no longer matches this run.

    An input that is gone cannot be checked, which is the normal state of a
    published study whose trajectories are archived elsewhere. Its recorded
    value is aggregated as it stands, and the returned warning says so. An
    input that is present but changed is refused.

    Returns
    -------
    list of str
        Warnings about inputs that could not be checked.

    Raises
    ------
    StaleCacheError
        If the result records no identity, or an input or any other compared
        field changed.
    """
    stored = artifact.provenance.get("identity") or {}
    recorded = stored.get("inputs") or []
    current = identity.restat(recorded, _working_dir(condition, replicate))
    unverified = [entry for entry in current if entry.get("missing")]
    comparable = [
        original if now.get("missing") else now
        for original, now in zip(recorded, current, strict=True)
    ]
    reason = identity.identity_mismatch(
        stored,
        identity.replicate_identity(
            analysis.plugin, condition.sim_config, settings, equilibration, comparable
        ),
    )
    if reason is not None:
        raise StaleCacheError(
            f"{analysis.name}: cached result {path} for condition='{condition.label}' "
            f"replicate={replicate} no longer matches this run: {reason}. Rerun with "
            "--recompute to recompute this replicate."
        )
    if not unverified:
        return []
    message = (
        f"replicate {replicate}: {len(unverified)} input file(s) are not on disk, so the "
        f"recorded values were aggregated without checking them "
        f"({', '.join(str(entry.get('relative_path') or entry['path']) for entry in unverified[:3])})"
    )
    logger.warning("%s: %s [condition=%s]", analysis.name, message, condition.label)
    return [message]


def _working_dir(condition: Condition, replicate: int) -> Path | None:
    """The replicate's working directory, or ``None`` when the config cannot say."""
    try:
        return Path(condition.sim_config.get_working_directory(replicate))
    except Exception:
        return None


def _aggregate(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    output_dir: Path,
    results: Sequence[ReplicateArtifact],
    replicates: Sequence[int],
    recompute: bool,
    *,
    warnings: Sequence[str] = (),
) -> ConditionArtifact:
    """Aggregate replicate artifacts by observable kind and write the result."""
    aggregated_dir = output_dir / "aggregated"
    if recompute and aggregated_dir.exists():
        shutil.rmtree(aggregated_dir)
    aggregated_dir.mkdir(parents=True, exist_ok=True)
    try:
        estimates = [
            [ObservableEstimate.model_validate(item) for item in result.payload["observables"]]
            for result in results
        ]
        aggregates = aggregate_observables(estimates)
    except PluginContractError:
        raise
    except Exception as exc:
        raise AggregationError(
            f"{analysis.name}: aggregate failed for condition='{condition.label}': "
            f"{type(exc).__name__}: {exc}"
        ) from exc
    first = dict(results[0].provenance.get("identity") or {})
    first.pop("inputs", None)
    sources = [
        {"replicate": int(replicate), "fingerprint": identity.file_fingerprint(path)}
        for replicate in replicates
        if (path := output_dir / f"run_{replicate}" / RESULT_FILE).is_file()
    ]
    artifact = ConditionArtifact(
        analysis_name=analysis.name,
        condition_label=condition.label,
        replicates=[int(replicate) for replicate in replicates],
        source_replicates=sources,
        payload={"observables": [item.model_dump(mode="json") for item in aggregates]},
        provenance={"source": "observable_contract", "identity": first},
        metadata={
            "settings_fingerprint": identity.settings_fingerprint(settings),
            "equilibration": str(equilibration),
            "config_hash": identity.compute_config_hash(condition.sim_config),
            "n_replicates": len(results),
        },
        warnings=list(warnings),
    )
    try:
        ArtifactStore(aggregated_dir).write_condition_result(artifact, RESULT_FILE)
    except OSError as exc:
        raise AggregationError(
            f"{analysis.name}: failed to save aggregated result for "
            f"condition='{condition.label}': {exc}"
        ) from exc
    return artifact


def _check_aggregate(
    artifact: ConditionArtifact,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
) -> None:
    """Refuse an aggregate computed for other settings, window or config.

    Raises
    ------
    AggregateValidationError
        If the settings fingerprint is missing or differs, the equilibration
        window or config hash differs, or the aggregate covers a replicate the
        condition does not list.
    """
    if artifact.metadata.get("settings_fingerprint") is None:
        raise AggregateValidationError("aggregate is missing settings fingerprint metadata")
    expected = {
        "settings_fingerprint": identity.settings_fingerprint(settings),
        "equilibration": str(equilibration),
        "config_hash": identity.compute_config_hash(condition.sim_config),
    }
    for key, value in expected.items():
        stored = artifact.metadata.get(key)
        if stored is not None and stored != value:
            raise AggregateValidationError(f"{key} mismatch: stored {stored}, current {value}")
    extra = sorted(set(artifact.replicates) - set(condition.replicates))
    if not artifact.replicates or extra:
        raise AggregateValidationError(
            f"aggregate covers replicates {artifact.replicates}, not a subset of "
            f"{list(condition.replicates)}"
        )


# ---------------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------------


def prepare_comparison_run(
    analysis: Analysis,
    config: ComparisonConfig,
    equilibration: str | None,
) -> dict[str, Any]:
    """Resolve the conditions, settings, window and output root of a comparison.

    Returns
    -------
    dict
        ``all_conditions``, ``valid_conditions``, ``excluded_conditions``
        (always empty), ``condition_by_label``, ``settings``,
        ``equilibration`` and ``analysis_root``.
    """
    conditions = [Condition.from_condition_config(item) for item in config.conditions]
    source_path = getattr(config, "source_path", None)
    return {
        "all_conditions": conditions,
        "valid_conditions": list(conditions),
        "excluded_conditions": [],
        "condition_by_label": {condition.label: condition for condition in conditions},
        "settings": _resolve_settings(analysis, config),
        "equilibration": equilibration or config.defaults.equilibration_time,
        "analysis_root": source_path.parent / "analysis" if source_path else Path("analysis"),
    }


def run_comparison(
    analysis: Analysis,
    config: ComparisonConfig,
    recompute: bool = False,
    equilibration: str | None = None,
) -> dict[str, Any]:
    """Run one analysis over every condition of a comparison.

    A condition that fails is logged and left out, and finalizing then refuses
    the incomplete comparison. A plugin contract error stops everything.

    Returns
    -------
    dict
        ``aggregated`` (label to aggregate), ``comparison``,
        ``comparison_path`` and ``plots``.
    """
    from polyzymd.analyses.shared.paths import sanitize_label

    prepared = prepare_comparison_run(analysis, config, equilibration)
    conditions = prepared["valid_conditions"]
    settings = prepared["settings"]
    equilibration = prepared["equilibration"]
    analysis_root = prepared["analysis_root"]
    _print_execution_summary(analysis, conditions, settings, equilibration)

    analysis_dirs: dict[str, Path] = {}
    aggregated: dict[str, Any] = {}
    failures: list[tuple[str, Exception]] = []
    for condition in conditions:
        directory = analysis_root / sanitize_label(condition.label) / analysis.name
        try:
            aggregated[condition.label] = run_analysis(
                analysis, condition, settings, equilibration, directory, recompute
            )
            analysis_dirs[condition.label] = directory
        except PluginContractError:
            raise
        except (AnalysisError, ValueError, FileNotFoundError, OSError) as exc:
            logger.error(f"  {condition.label}: {type(exc).__name__} — {exc}")
            failures.append((condition.label, exc))
    if not aggregated:
        raise ValueError(_no_conditions_message(analysis.name, failures))

    control = config.control if config.control in aggregated else None
    final_config = config.model_copy(deep=True)
    final_config.defaults.equilibration_time = equilibration
    result = finalize_comparison_from_disk(
        analysis=analysis,
        config=final_config,
        analysis_dirs=analysis_dirs,
        aggregated_results=aggregated,
        results_dir=analysis_root.parent / "comparison" / analysis.name,
        figures_dir=analysis.figures_output_dir(_figures_root(config, analysis_root)),
        settings=settings,
        effective_control=control,
        prepared_state=prepared,
        recompute=recompute,
    )
    return {"aggregated": aggregated, **result}


def finalize_comparison_from_disk(
    analysis: Analysis,
    config: ComparisonConfig,
    analysis_dirs: dict[str, Path],
    aggregated_results: dict[str, Any],
    results_dir: Path,
    figures_dir: Path,
    settings: BaseModel,
    effective_control: str | None,
    prepared_state: dict[str, Any] | None = None,
    allow_partial: bool = False,
    recompute: bool = False,
) -> dict[str, Any]:
    """Compare the conditions' aggregates and draw the figures.

    An aggregate not passed in is read from ``<condition dir>/aggregated``.
    One that is missing, unreadable, or computed for other settings counts as
    missing. Without ``allow_partial`` any missing aggregate is an error.

    Returns
    -------
    dict
        ``comparison``, ``comparison_path`` and ``plots``.
    """
    if prepared_state is None:
        prepared = prepare_comparison_run(analysis, config, config.defaults.equilibration_time)
        prepared["settings"] = settings
    else:
        prepared = prepared_state
    by_label = prepared["condition_by_label"]
    equilibration = prepared["equilibration"]
    settings = prepared["settings"]

    valid: dict[str, ConditionArtifact] = {}
    for label, directory in analysis_dirs.items():
        condition = by_label.get(label)
        if condition is None:
            continue
        artifact = aggregated_results.get(label)
        try:
            if artifact is None:
                artifact = analysis._load_aggregated_result(directory / "aggregated")
            if artifact is None:
                logger.warning(
                    "%s: missing aggregated result for '%s' at %s. Use --allow-partial (CLI) / "
                    "allow_partial=True (API) to continue with available conditions.",
                    analysis.name,
                    label,
                    directory / "aggregated" / RESULT_FILE,
                )
                continue
            _check_aggregate(artifact, condition, settings, equilibration)
        except (AggregateValidationError, ArtifactStoreError) as exc:
            logger.warning("%s: invalid aggregated result for '%s': %s", analysis.name, label, exc)
            continue
        valid[label] = artifact

    dropped = [c for c in prepared["valid_conditions"] if c.label not in valid]
    if dropped:
        expected = _format_expected_paths(
            [(c.label, _expected_aggregate(analysis, c, analysis_dirs, prepared)) for c in dropped]
        )
        if not allow_partial:
            raise ValueError(
                f"{analysis.name}: missing aggregated results for condition(s): "
                f"{[c.label for c in dropped]}.{expected}\n"
                "Re-run failed condition aggregate jobs or use --allow-partial (CLI) / "
                "allow_partial=True (API) to continue with available conditions."
            )
        logger.warning(
            "%s: dropped condition(s) during finalize: %s",
            analysis.name,
            [c.label for c in dropped],
        )
    if not valid:
        raise ValueError(
            f"{analysis.name}: no successful conditions remain after finalize; no aggregate "
            "files were found for the finalized condition set. Re-run aggregate jobs before "
            "finalizing."
        )

    control = config.control if config.control is not None else effective_control
    if control is not None and control not in valid:
        if not allow_partial:
            raise ValueError(
                f"{analysis.name}: control condition '{control}' is missing from "
                "successful finalized conditions."
            )
        logger.warning(
            "%s: configured control '%s' was dropped during partial finalization; "
            "comparison will proceed without a designated control (all-vs-all).",
            analysis.name,
            control,
        )
        control = None
    conditions = [c for c in prepared["valid_conditions"] if c.label in valid]
    logger.info(
        "%s: finalizing comparison with conditions=%s effective_control=%s",
        analysis.name,
        [c.label for c in conditions],
        control,
    )

    results_dir.mkdir(parents=True, exist_ok=True)
    comparison_path = results_dir / RESULT_FILE
    if recompute and comparison_path.exists():
        comparison_path.unlink()
    comparison_ctx = ComparisonContext(
        name=config.name,
        conditions=conditions,
        excluded_conditions=[c for c in prepared["all_conditions"] if c.label not in valid],
        control_label=control,
        analysis_dirs={label: analysis_dirs[label] for label in valid},
        results_dir=results_dir,
        equilibration=equilibration,
        settings=settings,
        fdr_alpha=getattr(config.defaults, "fdr_alpha", 0.05),
        ttest_method=getattr(config.defaults, "ttest_method", "student"),
        posthoc_method=getattr(config.defaults, "posthoc_method", "ttest_bh"),
        result_path=comparison_path,
        failed_conditions=dropped,
        aggregated_results=valid,
        recompute=recompute,
    )
    try:
        comparison = analysis.compare(comparison_ctx)
    except PluginContractError:
        raise
    except Exception as exc:
        raise ComparisonError(
            f"{analysis.name}: compare failed for comparison='{config.name}': "
            f"{type(exc).__name__}: {exc}"
        ) from exc
    try:
        ArtifactStore(results_dir).write_comparison_result(comparison, RESULT_FILE)
    except OSError as exc:
        raise ComparisonError(f"{analysis.name}: failed to save comparison result: {exc}") from exc

    if recompute and figures_dir.exists():
        shutil.rmtree(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)
    plot_ctx = PlotContext(
        conditions=conditions,
        analysis_dirs={label: analysis_dirs[label] for label in valid},
        results_dir=results_dir,
        output_dir=figures_dir,
        settings=settings,
        plot_settings=_plot_settings(config),
        recompute=recompute,
        comparison_path=comparison_path,
        control_label=control,
        equilibration=equilibration,
    )
    try:
        plots = analysis.plot(plot_ctx)
    except PluginContractError:
        raise
    except Exception as exc:
        raise PlotError(
            f"{analysis.name}: plot failed for comparison='{config.name}': "
            f"{type(exc).__name__}: {exc}"
        ) from exc
    return {"comparison": comparison, "comparison_path": comparison_path, "plots": list(plots)}


# ---------------------------------------------------------------------------
# Plots and several analyses
# ---------------------------------------------------------------------------


def run_plot_only(
    analysis: Analysis,
    config: ComparisonConfig,
    equilibration: str | None = None,
) -> tuple[list[Path], list[tuple[str, str]]]:
    """Redraw one analysis's figures from results already on disk.

    Returns
    -------
    tuple
        ``(paths, failures)``. An expected plotting failure is returned as
        ``[(name, message)]`` instead of raised.
    """
    from polyzymd.analyses.shared.paths import sanitize_label

    prepared = prepare_comparison_run(analysis, config, equilibration)
    analysis_root = prepared["analysis_root"]
    conditions = prepared["valid_conditions"]
    results_dir = analysis_root.parent / "comparison" / analysis.name
    figures_dir = analysis.figures_output_dir(_figures_root(config, analysis_root).resolve())
    figures_dir.mkdir(parents=True, exist_ok=True)
    control = config.control if any(c.label == config.control for c in conditions) else None
    plot_ctx = PlotContext(
        conditions=conditions,
        analysis_dirs={
            c.label: analysis_root / sanitize_label(c.label) / analysis.name for c in conditions
        },
        results_dir=results_dir,
        output_dir=figures_dir,
        settings=prepared["settings"],
        plot_settings=_plot_settings(config),
        comparison_path=results_dir / RESULT_FILE,
        control_label=control,
        equilibration=prepared["equilibration"],
    )
    try:
        return list(analysis.plot(plot_ctx)), []
    except PluginContractError:
        raise
    except _PLOT_ONLY_EXPECTED_FAILURES as exc:
        error = (
            exc
            if isinstance(exc, PlotError)
            else PlotError(
                f"{analysis.name}: plot failed for comparison='{config.name}': "
                f"{type(exc).__name__}: {exc}"
            )
        )
        logger.error("Failed to generate plots for %s: %s", analysis.name, error)
        return [], [(analysis.name, str(error))]


def run_all_plots(
    config: ComparisonConfig,
    analysis_names: list[str] | None = None,
    equilibration: str | None = None,
) -> tuple[list[Path], list[tuple[str, str]]]:
    """Redraw the figures of several analyses, collecting failures."""
    from polyzymd.analyses.discovery import get_analysis

    generated: list[Path] = []
    failures: list[tuple[str, str]] = []
    for name in analysis_names if analysis_names is not None else _enabled(config):
        try:
            analysis_cls = get_analysis(name)
        except KeyError:
            failures.append((name, f"Unknown analysis type {name!r}"))
            continue
        paths, errors = run_plot_only(analysis_cls(), config, equilibration=equilibration)
        generated.extend(paths)
        failures.extend(errors)
    return generated, failures


def run_all_comparisons(
    config: ComparisonConfig,
    analysis_names: list[str] | None = None,
    recompute: bool = False,
    equilibration: str | None = None,
) -> dict[str, dict[str, Any]]:
    """Run several analyses, recording a failure as ``{"error": message}``."""
    from polyzymd.analyses.discovery import get_analysis

    results: dict[str, dict[str, Any]] = {}
    for name in analysis_names if analysis_names is not None else _enabled(config):
        try:
            analysis = get_analysis(name)()
        except KeyError:
            logger.warning(f"Unknown analysis type {name!r} — skipping.")
            continue
        logger.info(f"{'=' * 60}\nRunning {analysis.name} comparison\n{'=' * 60}")
        try:
            results[analysis.name] = run_comparison(
                analysis, config, recompute, equilibration=equilibration
            )
        except PluginContractError:
            raise
        except (AnalysisError, ValueError, FileNotFoundError, OSError) as exc:
            logger.error(f"{analysis.name} comparison failed: {exc}")
            results[analysis.name] = {"error": str(exc)}
    return results


def order_analyses_for_execution(analysis_names: Sequence[str]) -> list[str]:
    """Canonical analysis names in the order given, without duplicates.

    Analyses do not depend on each other, so the order is the requested one.

    Raises
    ------
    KeyError
        If a name is not a known analysis.
    """
    from polyzymd.analyses.discovery import get_analysis

    ordered: list[str] = []
    for name in analysis_names:
        canonical = get_analysis(name).name
        if canonical not in ordered:
            ordered.append(canonical)
    return ordered


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_settings(analysis: Analysis, config: ComparisonConfig) -> Any:
    """The analysis's settings from the comparison config, or its defaults."""
    plugins = getattr(config, "plugins", None)
    value = plugins.get(analysis.name) if plugins is not None else None
    if value is None:
        return analysis.Settings()
    if isinstance(value, analysis.Settings):
        return value
    if hasattr(value, "model_dump"):
        value = value.model_dump()
    return analysis.Settings.model_validate(value)


def _plot_settings(config: ComparisonConfig) -> Any:
    """Plot settings of a comparison, defaults when it has none."""
    from polyzymd.config.comparison import PlotSettings

    settings = getattr(config, "plot_settings", None)
    return settings if isinstance(settings, PlotSettings) else PlotSettings()


def _figures_root(config: ComparisonConfig, analysis_root: Path) -> Path:
    """Figure directory of a comparison, a relative path resolved against its file."""
    output_dir = Path(_plot_settings(config).output_dir)
    if output_dir.is_absolute():
        return output_dir
    source_path = getattr(config, "source_path", None)
    base = source_path.parent if source_path is not None else analysis_root.parent
    return base / output_dir


def _enabled(config: ComparisonConfig) -> list[str]:
    """Names of the analyses a comparison config configures."""
    plugins = getattr(config, "plugins", None)
    return plugins.get_enabled_plugins() if hasattr(plugins, "get_enabled_plugins") else []


def _expected_aggregate(
    analysis: Analysis,
    condition: Condition,
    analysis_dirs: dict[str, Path],
    prepared: dict[str, Any],
) -> Path | None:
    """Where a condition's aggregate was expected, for error messages."""
    from polyzymd.analyses.shared.paths import sanitize_label

    directory = analysis_dirs.get(condition.label)
    root = prepared.get("analysis_root")
    if directory is None and root is not None:
        directory = root / sanitize_label(condition.label) / analysis.name
    return None if directory is None else directory / "aggregated" / RESULT_FILE


def _format_expected_paths(entries: Sequence[tuple[str, Path | None]], *, limit: int = 5) -> str:
    """List where each missing condition's aggregate was expected."""
    if not entries:
        return ""
    shown = list(entries[:limit])
    lines = ["Expected paths:"]
    lines += [f"  - {label}: {path if path is not None else 'unknown'}" for label, path in shown]
    if len(entries) > limit:
        lines.append(f"  - ... {len(entries) - limit} more path(s) omitted")
    return "\n" + "\n".join(lines)


def _format_reasons(reasons: Sequence[str], *, limit: int = 3) -> str:
    """List why replicates failed, for the error that says none succeeded."""
    if not reasons:
        return ""
    lines = [" Failure reasons:"] + [f"  - {reason}" for reason in reasons[:limit]]
    if len(reasons) > limit:
        lines.append(f"  - ... {len(reasons) - limit} more failure(s) omitted")
    return "\n" + "\n".join(lines)


def _no_conditions_message(name: str, failures: Sequence[tuple[str, BaseException]]) -> str:
    """Quote the shared cause when every condition failed the same way.

    Typed analysis errors carry their own remedy, so a user sees what to do
    rather than only the summary line.
    """
    base = f"{name}: no conditions succeeded analysis."
    if not failures or len({type(error) for _, error in failures}) > 1:
        return base
    label, error = failures[0]
    return f"{base} Every condition failed with {type(error).__name__}. {label}: {error}"


def _print_execution_summary(
    analysis: Analysis,
    conditions: list[Condition],
    settings: BaseModel,
    equilibration: str,
) -> None:
    """Log what is about to run and suggest SLURM for a long local run."""
    del settings
    total = sum(len(condition.replicates) for condition in conditions)
    logger.info(
        "%s\n  %s — %d conditions × %d total replicate tasks\n  Mode: sequential (local)\n"
        "  Equilibration: %s\n%s",
        "=" * 60,
        analysis.name,
        len(conditions),
        total,
        equilibration,
        "=" * 60,
    )
    expensive = getattr(analysis, "execution_cost_hint", "medium") == "high"
    if not expensive and total <= _MANY_TASKS_THRESHOLD:
        return
    if shutil.which("sbatch") is not None:
        logger.warning(
            "This analysis may take a long time to run locally\n"
            "Consider submitting to SLURM for parallel execution:\n"
            "  polyzymd compare submit %s",
            analysis.name,
        )
    else:
        logger.warning(
            "This analysis may take a long time to run locally\n"
            "If you have access to an HPC cluster with SLURM, consider:\n"
            "  polyzymd compare submit %s",
            analysis.name,
        )
