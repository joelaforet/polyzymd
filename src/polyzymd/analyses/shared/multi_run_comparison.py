"""Shared helpers for multi-run comparison orchestration.

These helpers keep run-wise comparison logic concise across plugins that
compare multiple named runs (RMSD, Rg, SASA).
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping, Sequence
from typing import Any


def filter_summaries_with_run(
    summaries: dict[str, Any],
    run_label: str,
    get_run_fn: Callable[[Any, str], Any],
    logger: logging.Logger | None = None,
) -> dict[str, Any]:
    """Filter condition summaries to those containing a specific run.

    Parameters
    ----------
    summaries : dict[str, Any]
        Mapping from condition label to condition summary.
    run_label : str
        Run label to keep.
    get_run_fn : Callable[[Any, str], Any]
        Callback that returns run summary for ``(summary, run_label)`` and
        raises ``KeyError`` when the run is missing.
    logger : logging.Logger | None, optional
        Optional logger for missing-run warnings.

    Returns
    -------
    dict[str, Any]
        Subset of ``summaries`` with run data available.
    """
    filtered: dict[str, Any] = {}
    for label, summary in summaries.items():
        try:
            get_run_fn(summary, run_label)
        except KeyError:
            if logger is not None:
                logger.warning(
                    "Run '%s' missing for condition '%s'; excluding from run-level comparison",
                    run_label,
                    label,
                )
            continue
        filtered[label] = summary
    return filtered


def build_condition_pairs(
    condition_labels: list[str],
    control_label: str | None,
    on_control_missing: str = "all_pairs",
    logger: logging.Logger | None = None,
) -> list[tuple[str, str]]:
    """Build pairwise condition pairs for comparison.

    Parameters
    ----------
    condition_labels : list[str]
        Ordered condition labels to compare.
    control_label : str | None
        Preferred control label for control-vs-treatment comparisons.
    on_control_missing : str, optional
        Behavior when ``control_label`` is requested but unavailable.

        Supported values:

        - ``"all_pairs"``: fall back to all-vs-all
        - ``"skip"``: return no pairs
    logger : logging.Logger | None, optional
        Optional logger for fallback/skip messages.

    Returns
    -------
    list[tuple[str, str]]
        Pair list as ``(condition_a, condition_b)`` tuples.

    Raises
    ------
    ValueError
        Raised when ``on_control_missing`` is not ``"all_pairs"`` or ``"skip"``.
    """
    if on_control_missing not in ("all_pairs", "skip"):
        raise ValueError(
            f"on_control_missing must be 'all_pairs' or 'skip', got {on_control_missing!r}"
        )

    if len(condition_labels) < 2:
        return []

    if control_label is not None:
        if control_label in condition_labels:
            return [(control_label, label) for label in condition_labels if label != control_label]

        if on_control_missing == "skip":
            if logger is not None:
                logger.warning(
                    "Control condition '%s' unavailable; skipping pairwise comparisons",
                    control_label,
                )
            return []

        if logger is not None:
            logger.warning(
                "Control condition '%s' unavailable; falling back to all-vs-all pairwise comparisons",
                control_label,
            )

    return [
        (condition_labels[i], condition_labels[j])
        for i in range(len(condition_labels))
        for j in range(i + 1, len(condition_labels))
    ]


def apply_fdr_correction(
    pairwise_results: list[Any],
    anova_by_run: dict[Any, Any] | list[Any] | None = None,
    fdr_alpha: float = 0.05,
    get_p_value: Callable[[Any], float | None] | None = None,
    set_corrected: Callable[[Any, Any], None] | None = None,
) -> None:
    """Apply Benjamini-Hochberg FDR correction across statistical result families.

    Parameters
    ----------
    pairwise_results : list[Any]
        Pairwise comparison result objects.
    anova_by_run : dict[Any, Any] | list[Any] | None, optional
        ANOVA result objects, as either list-like or dict-like container.
    fdr_alpha : float, optional
        FDR threshold.
    get_p_value : Callable[[Any], float | None] | None, optional
        Callback extracting raw p-value from a result object. Defaults to
        reading ``.p_value``.
    set_corrected : Callable[[Any, Any], None] | None, optional
        Callback applying BH output to each result object. Defaults to setting
        ``.p_value_adjusted`` (when available) and ``.significant``.
    """
    from polyzymd.compare.statistics import benjamini_hochberg

    def _default_get_p_value(result: Any) -> float | None:
        return getattr(result, "p_value", None)

    def _default_set_corrected(result: Any, bh_result: Any) -> None:
        if hasattr(result, "p_value_adjusted"):
            result.p_value_adjusted = bh_result.adjusted_p_value
        result.significant = bh_result.significant

    p_getter = get_p_value or _default_get_p_value
    corrected_setter = set_corrected or _default_set_corrected

    if pairwise_results:
        pairwise_p_values = [p_getter(result) for result in pairwise_results]
        pairwise_bh = benjamini_hochberg(pairwise_p_values, alpha=fdr_alpha)
        for result, bh_result in zip(pairwise_results, pairwise_bh, strict=False):
            corrected_setter(result, bh_result)

    anova_items = _coerce_result_sequence(anova_by_run)
    if anova_items:
        anova_p_values = [p_getter(result) for result in anova_items]
        anova_bh = benjamini_hochberg(anova_p_values, alpha=fdr_alpha)
        for result, bh_result in zip(anova_items, anova_bh, strict=False):
            corrected_setter(result, bh_result)


def _coerce_result_sequence(results: Mapping[Any, Any] | Sequence[Any] | None) -> list[Any]:
    """Normalize mapping or sequence result containers to a list."""
    if results is None:
        return []
    if isinstance(results, Mapping):
        return list(results.values())
    return list(results)
