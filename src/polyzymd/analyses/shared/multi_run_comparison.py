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
    """Correct a run-wise comparison using the package-wide family policy.

    This is a thin wrapper over
    :func:`polyzymd.analyses.shared.inferential_statistics.apply_family_correction`,
    which owns the definition of the correction family. All pairwise tests
    form one Benjamini-Hochberg family; ANOVA results are reported as
    uncorrected omnibus tests.

    Parameters
    ----------
    pairwise_results : list[Any]
        Pairwise comparison result objects, mutated in place.
    anova_by_run : dict[Any, Any] | list[Any] | None, optional
        Omnibus ANOVA result objects, as either list-like or dict-like
        container.
    fdr_alpha : float, optional
        FDR threshold for the pairwise family and plain alpha for the
        omnibus ANOVA.
    get_p_value : Callable[[Any], float | None] | None, optional
        Callback extracting raw p-value from a result object. Defaults to
        reading ``.p_value``.
    set_corrected : Callable[[Any, Any], None] | None, optional
        Callback applying BH output to each result object. Defaults to
        setting ``.p_value_adjusted`` (when available) and ``.significant``.
    """
    from polyzymd.analyses.shared.inferential_statistics import apply_family_correction

    anova_items = _coerce_result_sequence(anova_by_run)
    if set_corrected is None:
        _validate_default_setter_targets(pairwise_results, "pairwise_results")
    _validate_anova_targets(anova_items)

    apply_family_correction(
        pairwise_results,
        fdr_alpha=fdr_alpha,
        get_p_value=get_p_value,
        set_corrected=set_corrected,
        anova_results=anova_items,
    )


def _validate_default_setter_targets(results: Sequence[Any], label: str) -> None:
    """Check that the default setter can write to every result object.

    Parameters
    ----------
    results : Sequence[Any]
        Result objects the default setter will mutate.
    label : str
        Name of the container, used in the error message.

    Raises
    ------
    TypeError
        Raised when a result lacks a writable ``significant`` or
        ``p_value_adjusted`` attribute.
    """
    for idx, result in enumerate(results):
        result_type = type(result).__name__
        for attribute in ("significant", "p_value_adjusted"):
            if not hasattr(result, attribute):
                raise TypeError(
                    "apply_fdr_correction() default setter requires results with a "
                    f"{attribute!r} attribute. {label}[{idx}] has type {result_type}. "
                    "Provide set_corrected=... for custom result objects."
                )
            try:
                setattr(result, attribute, getattr(result, attribute))
            except (AttributeError, TypeError) as exc:
                raise TypeError(
                    "apply_fdr_correction() default setter requires a mutable "
                    f"{attribute!r} attribute. {label}[{idx}] has type {result_type}. "
                    "Provide set_corrected=... for custom result objects."
                ) from exc


def _validate_anova_targets(results: Sequence[Any]) -> None:
    """Check that omnibus ANOVA results can record the policy's verdict.

    Parameters
    ----------
    results : Sequence[Any]
        ANOVA result objects.

    Raises
    ------
    TypeError
        Raised when an ANOVA result lacks a writable ``significant`` or
        ``p_value_adjusted`` attribute.
    """
    _validate_default_setter_targets(results, "anova_by_run")


def _coerce_result_sequence(results: Mapping[Any, Any] | Sequence[Any] | None) -> list[Any]:
    """Normalize mapping or sequence result containers to a list.

    Parameters
    ----------
    results : Mapping[Any, Any] | Sequence[Any] | None
        Container of result objects.

    Returns
    -------
    list[Any]
        The contained results, in order.
    """
    if results is None:
        return []
    if isinstance(results, Mapping):
        return list(results.values())
    return list(results)
