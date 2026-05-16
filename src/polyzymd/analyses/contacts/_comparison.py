"""Comparison and statistics helpers for the contacts analysis facade."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from datetime import datetime
from typing import Any, Sequence

import numpy as np

from polyzymd.analyses.base import ComparisonContext, Condition
from polyzymd.analyses.contacts._aggregator import CONTACTS_LEGACY_RECOMPUTE_GUIDANCE
from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact
from polyzymd.analyses.mda.store import ArtifactStoreError

logger = logging.getLogger("polyzymd.analyses.contacts")

NOT_TESTABLE_SINGLETON_NOTE = (
    "Inferential statistics require at least two replicates per condition."
)
CONTACTS_COMPARED_METRICS = ("coverage", "mean_contact_fraction")
METRIC_STAT_TOLERANCE = 1e-12


def compare(analysis: Any, ctx: ComparisonContext) -> Any:
    """Compare contacts metrics across conditions.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for cache validation and helper wrappers.
    ctx : ComparisonContext
        Framework-provided comparison context.

    Returns
    -------
    Any
        Contacts comparison result, or ``None`` when no condition data is available.
    """

    from polyzymd import __version__
    from polyzymd.analyses.contacts._comparison_results import (
        ContactsANOVASummary,
        ContactsComparisonResult,
        ContactsConditionSummary,
        ContactsPairwiseComparison,
    )

    settings = ctx.settings

    logger.info(f"Starting contacts comparison: {ctx.name}")
    logger.info(f"Conditions: {len(ctx.conditions)}")
    logger.info(f"Equilibration: {ctx.equilibration}")

    if not ctx.conditions:
        logger.warning("contacts: no conditions provided — skipping comparison.")
        return None

    condition_data: list[tuple[Condition, dict[str, Any]]] = []
    for cond in ctx.conditions:
        agg_result = _load_condition_artifact(analysis, ctx, cond)
        if agg_result is None:
            logger.warning(f"No aggregated result for '{cond.label}' — skipping.")
            continue

        coverage_per_rep = analysis._compute_coverage_per_replicate(agg_result)
        contact_fraction_per_rep = analysis._compute_contact_fraction_per_replicate(agg_result)

        condition_data.append(
            (
                cond,
                {
                    "agg_result": agg_result,
                    "coverage_per_replicate": coverage_per_rep,
                    "contact_fraction_per_replicate": contact_fraction_per_rep,
                },
            )
        )

    if not condition_data:
        logger.warning("contacts: no conditions have valid results — skipping.")
        return None

    analysis._validate_residue_sets(condition_data)

    summaries: list[ContactsConditionSummary] = []
    for cond, data in condition_data:
        agg_result = data["agg_result"]
        coverage_metric = agg_result.payload["metrics"]["coverage"]
        contact_metric = agg_result.payload["metrics"]["mean_contact_fraction"]
        residence_summary = (
            agg_result.payload.get("residence_time_by_polymer_type", {})
            if settings.compute_residence_times
            else {}
        )
        summary = ContactsConditionSummary(
            label=cond.label,
            config_path=str(cond.config_path),
            n_replicates=len(agg_result.replicates),
            n_residues=int(agg_result.payload.get("n_residues", 0)),
            coverage_mean=float(coverage_metric["mean"]),
            coverage_sem=float(coverage_metric["sem"]),
            mean_contact_fraction=float(contact_metric["mean"]),
            mean_contact_fraction_sem=float(contact_metric["sem"]),
            residence_time_by_polymer_type=residence_summary,
            compute_residence_times=settings.compute_residence_times,
        )
        summaries.append(summary)

    effective_control = analysis._resolve_effective_control(ctx.effective_control, summaries)

    comparisons: list[ContactsPairwiseComparison] = []
    if len(summaries) >= 2:
        comparisons = analysis._compute_contacts_pairwise(
            summaries, condition_data, effective_control
        )

    anova_results: list[ContactsANOVASummary] = []
    if len(summaries) >= 3:
        anova_results = analysis._compute_contacts_anova(condition_data)

    analysis._apply_fdr_correction(comparisons, anova_results, settings.fdr_alpha)
    analysis._apply_effect_size_threshold(comparisons, settings.min_effect_size)

    ranked_coverage = sorted(summaries, key=lambda s: s.coverage_mean, reverse=True)
    ranked_contact = sorted(summaries, key=lambda s: s.mean_contact_fraction, reverse=True)
    top_residues_data = analysis._compute_top_contacted_residues(
        condition_data, settings.top_residues
    )
    return ContactsComparisonResult(
        name=ctx.name,
        contacts_name="polymer_contacts",
        contacts_description=None,
        polymer_selection=settings.polymer_selection,
        protein_selection=settings.protein_selection,
        cutoff=settings.cutoff,
        contact_criteria="distance",
        compute_residence_times=settings.compute_residence_times,
        fdr_alpha=settings.fdr_alpha,
        min_effect_size=settings.min_effect_size,
        top_residues=settings.top_residues,
        control_label=effective_control,
        conditions=summaries,
        pairwise_comparisons=comparisons,
        anova=anova_results,
        ranking_by_coverage=[s.label for s in ranked_coverage],
        ranking_by_contact_fraction=[s.label for s in ranked_contact],
        excluded_conditions=[c.label for c in ctx.excluded_conditions],
        top_contacted_residues=top_residues_data,
        equilibration_time=ctx.equilibration,
        created_at=datetime.now(),
        polyzymd_version=__version__,
    )


def _load_condition_artifact(
    analysis: Any, ctx: ComparisonContext, cond: Condition
) -> ConditionArtifact | None:
    """Load and validate one contacts condition artifact for comparison."""

    in_memory = ctx.aggregated_results.get(cond.label)
    if in_memory is not None:
        if not isinstance(in_memory, ConditionArtifact):
            raise ValueError(CONTACTS_LEGACY_RECOMPUTE_GUIDANCE)
        _validate_condition_artifact(in_memory, ctx, cond)
        return in_memory

    aggregate_path = ctx.analysis_dirs[cond.label] / "aggregated" / "result.json"
    if not aggregate_path.exists():
        return None
    try:
        artifact = ArtifactStore(aggregate_path.parent).read_condition_result("result.json")
    except ArtifactStoreError as exc:
        raise ValueError(
            f"contacts: aggregate at {aggregate_path} is not a canonical ConditionArtifact. "
            f"{CONTACTS_LEGACY_RECOMPUTE_GUIDANCE}"
        ) from exc
    _validate_condition_artifact(artifact, ctx, cond)
    del analysis
    return artifact


def _validate_condition_artifact(
    artifact: ConditionArtifact, ctx: ComparisonContext, cond: Condition
) -> None:
    """Validate contacts condition artifact identity for comparison."""

    if artifact.analysis_name != "contacts":
        raise ValueError(
            f"contacts: condition artifact for {cond.label!r} has analysis "
            f"{artifact.analysis_name!r}; expected 'contacts'"
        )
    if artifact.condition_label != cond.label:
        raise ValueError(
            f"contacts: condition artifact label mismatch: stored {artifact.condition_label!r}, "
            f"expected {cond.label!r}"
        )
    active_replicates = {int(replicate) for replicate in cond.replicates}
    stored_replicates = {int(replicate) for replicate in artifact.replicates}
    unexpected = sorted(stored_replicates - active_replicates)
    if unexpected:
        raise ValueError(
            f"contacts: condition {cond.label!r} contains unexpected replicate IDs {unexpected}; "
            "recompute contacts or clear stale aggregate result.json files"
        )
    if not stored_replicates:
        raise ValueError(f"contacts: condition {cond.label!r} has no successful replicates")
    _validate_condition_frame_selection(artifact, ctx, cond)
    _validate_condition_residence_setting(artifact, ctx, cond)
    _validate_condition_metric_integrity(artifact, cond)
    from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint

    expected_fingerprint = contacts_detection_fingerprint(ctx.settings)
    stored_fingerprint = artifact.metadata.get("contacts_detection_fingerprint")
    if stored_fingerprint != expected_fingerprint:
        raise ValueError(
            f"contacts: condition {cond.label!r} detection fingerprint mismatch; expected "
            f"{expected_fingerprint!r}, got {stored_fingerprint!r}. Recompute contacts."
        )


def _validate_condition_residence_setting(
    artifact: ConditionArtifact, ctx: ComparisonContext, cond: Condition
) -> None:
    """Validate residence-time aggregation identity for comparison."""

    expected = bool(getattr(ctx.settings, "compute_residence_times", True))
    stored = artifact.metadata.get("compute_residence_times")
    if stored is None:
        raise ValueError(
            f"contacts: condition {cond.label!r} aggregate lacks compute_residence_times "
            "metadata. Recompute contacts or clear stale aggregate result.json files."
        )
    if bool(stored) != expected:
        raise ValueError(
            f"contacts: condition {cond.label!r} compute_residence_times mismatch: stored "
            f"{bool(stored)!r}, expected {expected!r}. Recompute contacts or clear stale "
            "aggregate result.json files."
        )


def _validate_condition_metric_integrity(artifact: ConditionArtifact, cond: Condition) -> None:
    """Validate compared metric vectors and replicate-level identity."""

    expected_replicates = [int(replicate) for replicate in artifact.replicates]
    expected_keys = {str(replicate) for replicate in expected_replicates}
    replicate_metrics = artifact.payload.get("replicate_metrics")
    if not isinstance(replicate_metrics, Mapping):
        raise ValueError(
            f"contacts: condition {cond.label!r} lacks replicate_metrics. Recompute contacts."
        )
    stored_keys = {str(key) for key in replicate_metrics}
    if stored_keys != expected_keys:
        raise ValueError(
            f"contacts: condition {cond.label!r} replicate_metrics keys {sorted(stored_keys)} "
            f"do not match artifact replicates {sorted(expected_keys)}. Recompute contacts."
        )

    metrics = artifact.payload.get("metrics")
    if not isinstance(metrics, Mapping):
        raise ValueError(f"contacts: condition {cond.label!r} lacks metrics. Recompute contacts.")
    for metric_name in CONTACTS_COMPARED_METRICS:
        metric = metrics.get(metric_name)
        if not isinstance(metric, Mapping):
            raise ValueError(
                f"contacts: condition {cond.label!r} lacks metric {metric_name!r}. "
                "Recompute contacts."
            )
        _validate_metric_summary(metric, metric_name, expected_replicates, cond)


def _validate_metric_summary(
    metric: Mapping[str, Any], metric_name: str, replicates: Sequence[int], cond: Condition
) -> None:
    """Validate one compared metric summary against replicate values."""

    values = _metric_values(metric, metric_name, cond)
    stored_n = _metric_int(metric.get("n"), metric_name, "n", cond)
    expected_n = len(replicates)
    if len(values) != expected_n or stored_n != expected_n:
        raise ValueError(
            f"contacts: condition {cond.label!r} metric {metric_name!r} length mismatch: "
            f"len(values)={len(values)}, n={stored_n}, replicates={expected_n}. Recompute contacts."
        )
    mean = float(np.mean(values)) if values else 0.0
    std = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
    sem = float(std / np.sqrt(len(values))) if len(values) > 1 else 0.0
    expected_stats = {"mean": mean, "std": std, "sem": sem}
    for stat_name, expected_value in expected_stats.items():
        stored_value = _metric_float(metric.get(stat_name), metric_name, stat_name, cond)
        if not np.isclose(
            stored_value,
            expected_value,
            rtol=METRIC_STAT_TOLERANCE,
            atol=METRIC_STAT_TOLERANCE,
        ):
            raise ValueError(
                f"contacts: condition {cond.label!r} metric {metric_name!r} {stat_name} "
                f"mismatch: stored {stored_value!r}, recomputed {expected_value!r}. "
                "Recompute contacts."
            )


def _metric_values(metric: Mapping[str, Any], metric_name: str, cond: Condition) -> list[float]:
    """Return finite replicate values from one metric summary."""

    raw_values = metric.get("values")
    if not isinstance(raw_values, Sequence) or isinstance(raw_values, (str, bytes, bytearray)):
        raise ValueError(
            f"contacts: condition {cond.label!r} metric {metric_name!r} lacks a values vector. "
            "Recompute contacts."
        )
    values = [_metric_float(value, metric_name, "values", cond) for value in raw_values]
    if not all(np.isfinite(values)):
        raise ValueError(
            f"contacts: condition {cond.label!r} metric {metric_name!r} contains non-finite "
            "replicate values. Recompute contacts."
        )
    return values


def _metric_float(value: Any, metric_name: str, field: str, cond: Condition) -> float:
    """Coerce and validate one finite metric float."""

    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"contacts: condition {cond.label!r} metric {metric_name!r} has invalid {field}. "
            "Recompute contacts."
        ) from exc
    if not np.isfinite(result):
        raise ValueError(
            f"contacts: condition {cond.label!r} metric {metric_name!r} has non-finite "
            f"{field}. Recompute contacts."
        )
    return result


def _metric_int(value: Any, metric_name: str, field: str, cond: Condition) -> int:
    """Coerce and validate one metric integer field."""

    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"contacts: condition {cond.label!r} metric {metric_name!r} has invalid {field}. "
            "Recompute contacts."
        ) from exc


def _validate_condition_frame_selection(
    artifact: ConditionArtifact, ctx: ComparisonContext, cond: Condition
) -> None:
    """Validate condition artifact frame-selection metadata for comparison."""

    frame_selection = artifact.provenance.get("frame_selection")
    if not isinstance(frame_selection, Mapping):
        raise ValueError(
            f"contacts: condition {cond.label!r} lacks frame-selection provenance; "
            "recompute contacts or clear stale aggregate result.json files"
        )
    _validate_equilibration_provenance(
        stored=frame_selection.get("equilibration"),
        expected=ctx.equilibration,
        label=f"condition {cond.label!r} frame selection",
    )
    stored_metadata_equilibration = artifact.metadata.get("equilibration")
    if stored_metadata_equilibration is not None:
        _validate_equilibration_provenance(
            stored=stored_metadata_equilibration,
            expected=ctx.equilibration,
            label=f"condition {cond.label!r} metadata",
        )


def _validate_equilibration_provenance(*, stored: Any, expected: str, label: str) -> None:
    """Validate stored equilibration provenance against the active context."""

    if stored is None:
        raise ValueError(f"contacts: {label} lacks equilibration provenance; recompute contacts")
    if _equilibration_to_ps(stored) != _equilibration_to_ps(expected):
        raise ValueError(
            f"contacts: {label} equilibration mismatch: stored {stored!r}, expected "
            f"{expected!r}. Recompute contacts or clear stale aggregate result.json files."
        )


def _equilibration_to_ps(value: Any) -> float:
    """Normalize an equilibration string to picoseconds for comparisons."""

    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    try:
        numeric_value, unit = parse_time_string(str(value))
    except ValueError as exc:
        raise ValueError(f"contacts: invalid equilibration provenance {value!r}") from exc
    return float(convert_time(numeric_value, unit, "ps"))


def compute_coverage_per_replicate(result: Any) -> list[float]:
    """Compute coverage per replicate from residue statistics.

    Parameters
    ----------
    result : Any
        Aggregated contacts result with per-residue replicate fractions.

    Returns
    -------
    list[float]
        Fraction of residues contacted in each replicate.
    """

    if isinstance(result, ConditionArtifact):
        return [float(value) for value in result.payload["metrics"]["coverage"]["values"]]

    n_replicates = result.n_replicates
    n_residues = result.n_residues

    coverages = []
    for rep_idx in range(n_replicates):
        contacted = sum(
            1 for rs in result.residue_stats if rs.contact_fraction_per_replicate[rep_idx] > 0
        )
        coverages.append(contacted / n_residues if n_residues > 0 else 0.0)

    return coverages


def compute_contact_fraction_per_replicate(result: Any) -> list[float]:
    """Compute mean contact fraction per replicate.

    Parameters
    ----------
    result : Any
        Aggregated contacts result with per-residue replicate fractions.

    Returns
    -------
    list[float]
        Mean contact fraction for each replicate.
    """

    if isinstance(result, ConditionArtifact):
        return [
            float(value) for value in result.payload["metrics"]["mean_contact_fraction"]["values"]
        ]

    n_replicates = result.n_replicates

    fractions = []
    for rep_idx in range(n_replicates):
        rep_fractions = [rs.contact_fraction_per_replicate[rep_idx] for rs in result.residue_stats]
        mean_frac = float(np.mean(rep_fractions)) if rep_fractions else 0.0
        fractions.append(mean_frac)

    return fractions


def validate_residue_sets(condition_data: list[tuple[Condition, dict[str, Any]]]) -> None:
    """Validate that all conditions have identical residue sets.

    Parameters
    ----------
    condition_data : list[tuple[Condition, dict[str, Any]]]
        Condition data to validate.

    Raises
    ------
    ValueError
        If residue sets differ between conditions.
    """

    if len(condition_data) < 2:
        return

    first_cond, first_data = condition_data[0]
    first_result = first_data["agg_result"]
    first_resids = _residue_identity_set(first_result)

    for cond, data in condition_data[1:]:
        result = data["agg_result"]
        resids = _residue_identity_set(result)
        if resids != first_resids:
            missing_in_first = resids - first_resids
            missing_in_other = first_resids - resids
            raise ValueError(
                f"Residue set mismatch between '{first_cond.label}' and '{cond.label}'. "
                f"Missing in {first_cond.label}: {sorted(missing_in_first)}, "
                f"Missing in {cond.label}: {sorted(missing_in_other)}."
            )


def _residue_identity_set(result: Any) -> set[tuple[int, str, str]]:
    """Return chain-safe residue identities for legacy or condition artifacts."""

    if isinstance(result, ConditionArtifact):
        return {
            (
                int(row.get("protein_resid", 0)),
                str(row.get("protein_resname", "")),
                str(row.get("protein_chain_id", "")),
            )
            for row in result.payload.get("residue_stats", [])
        }
    return {
        (
            int(getattr(row, "protein_resid", 0)),
            _identity_text(getattr(row, "protein_resname", "")),
            _identity_text(getattr(row, "protein_chain_id", "")),
        )
        for row in getattr(result, "residue_stats", [])
    }


def _identity_text(value: Any) -> str:
    """Return stable identity text for optional legacy residue attributes."""

    return value if isinstance(value, str) else ""


def compute_contacts_pairwise(
    analysis: Any,
    summaries: list[Any],
    condition_data: list[tuple[Condition, dict[str, Any]]],
    effective_control: str | None,
) -> list[Any]:
    """Compute pairwise statistical comparisons for contacts.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for pairwise helper dispatch.
    summaries : list[Any]
        Condition summaries.
    condition_data : list[tuple[Condition, dict[str, Any]]]
        Raw condition data with per-replicate values.
    effective_control : str or None
        Control condition label.

    Returns
    -------
    list[Any]
        Pairwise comparison results.
    """

    comparisons = []
    label_to_data = {cond.label: data for cond, data in condition_data}
    label_to_summary = {s.label: s for s in summaries}

    if effective_control:
        control_data = label_to_data[effective_control]
        control_summary = label_to_summary[effective_control]

        for summary in summaries:
            if summary.label == effective_control:
                continue
            data = label_to_data[summary.label]
            comp = analysis._compare_contacts_pair(
                effective_control,
                control_summary,
                control_data,
                summary.label,
                summary,
                data,
            )
            comparisons.append(comp)
    else:
        labels = [s.label for s in summaries]
        for i, label_a in enumerate(labels):
            for label_b in labels[i + 1 :]:
                comp = analysis._compare_contacts_pair(
                    label_a,
                    label_to_summary[label_a],
                    label_to_data[label_a],
                    label_b,
                    label_to_summary[label_b],
                    label_to_data[label_b],
                )
                comparisons.append(comp)

    return comparisons


def resolve_effective_control(
    requested_control: str | None, summaries: Sequence[Any]
) -> str | None:
    """Return a control label that exists in validated summaries.

    Parameters
    ----------
    requested_control : str or None
        Control label from the comparison context.
    summaries : Sequence[Any]
        Validated condition summaries.

    Returns
    -------
    str or None
        The requested control when present, otherwise ``None``.
    """

    if requested_control is None:
        return None
    labels = {summary.label for summary in summaries}
    if requested_control in labels:
        return requested_control
    logger.warning(
        "contacts: configured control '%s' has no validated aggregate; falling back to all-pairs",
        requested_control,
    )
    return None


def compare_contacts_pair(
    label_a: str,
    summary_a: Any,
    data_a: dict[str, Any],
    label_b: str,
    summary_b: Any,
    data_b: dict[str, Any],
) -> Any:
    """Compare two conditions for coverage and contact fraction.

    Parameters
    ----------
    label_a, label_b : str
        Condition labels.
    summary_a, summary_b : Any
        Condition summaries.
    data_a, data_b : dict[str, Any]
        Raw data with per-replicate values.

    Returns
    -------
    Any
        Contacts pairwise comparison result.
    """

    from polyzymd.analyses.contacts._comparison_results import (
        AggregateComparisonResult,
        ContactsPairwiseComparison,
    )
    from polyzymd.analyses.shared.inferential_statistics import (
        cohens_d,
        independent_ttest,
        percent_change,
    )
    from polyzymd.analyses.stats import interpret_direction

    aggregate_comps = []

    coverage_a = data_a["coverage_per_replicate"]
    coverage_b = data_b["coverage_per_replicate"]
    coverage_testable = len(coverage_a) >= 2 and len(coverage_b) >= 2
    ttest = independent_ttest(coverage_a, coverage_b)
    effect = cohens_d(coverage_a, coverage_b)
    pct = percent_change(summary_a.coverage_mean, summary_b.coverage_mean)
    direction = interpret_direction(
        pct,
        direction_labels=("decreased", "unchanged", "increased"),
        threshold=1.0,
    )

    aggregate_comps.append(
        AggregateComparisonResult(
            metric="coverage",
            condition_a=label_a,
            condition_b=label_b,
            condition_a_mean=summary_a.coverage_mean,
            condition_a_sem=summary_a.coverage_sem,
            condition_b_mean=summary_b.coverage_mean,
            condition_b_sem=summary_b.coverage_sem,
            t_statistic=ttest.t_statistic,
            p_value=ttest.p_value,
            cohens_d=effect.cohens_d,
            effect_size_interpretation=effect.interpretation,
            significant=ttest.significant if coverage_testable else False,
            percent_change=pct,
            direction=direction,
            testable=coverage_testable,
            note=None if coverage_testable else NOT_TESTABLE_SINGLETON_NOTE,
        )
    )

    contact_a = data_a["contact_fraction_per_replicate"]
    contact_b = data_b["contact_fraction_per_replicate"]
    contact_testable = len(contact_a) >= 2 and len(contact_b) >= 2
    ttest = independent_ttest(contact_a, contact_b)
    effect = cohens_d(contact_a, contact_b)
    pct = percent_change(summary_a.mean_contact_fraction, summary_b.mean_contact_fraction)
    direction = interpret_direction(
        pct,
        direction_labels=("decreased", "unchanged", "increased"),
        threshold=1.0,
    )

    aggregate_comps.append(
        AggregateComparisonResult(
            metric="mean_contact_fraction",
            condition_a=label_a,
            condition_b=label_b,
            condition_a_mean=summary_a.mean_contact_fraction,
            condition_a_sem=summary_a.mean_contact_fraction_sem,
            condition_b_mean=summary_b.mean_contact_fraction,
            condition_b_sem=summary_b.mean_contact_fraction_sem,
            t_statistic=ttest.t_statistic,
            p_value=ttest.p_value,
            cohens_d=effect.cohens_d,
            effect_size_interpretation=effect.interpretation,
            significant=ttest.significant if contact_testable else False,
            percent_change=pct,
            direction=direction,
            testable=contact_testable,
            note=None if contact_testable else NOT_TESTABLE_SINGLETON_NOTE,
        )
    )

    return ContactsPairwiseComparison(
        condition_a=label_a,
        condition_b=label_b,
        aggregate_comparisons=aggregate_comps,
    )


def compute_contacts_anova(condition_data: list[tuple[Condition, dict[str, Any]]]) -> list[Any]:
    """Compute one-way ANOVA for both aggregate metrics.

    Parameters
    ----------
    condition_data : list[tuple[Condition, dict[str, Any]]]
        Condition data.

    Returns
    -------
    list[Any]
        ANOVA summaries for coverage and mean contact fraction.
    """

    from polyzymd.analyses.contacts._comparison_results import ContactsANOVASummary
    from polyzymd.analyses.shared.inferential_statistics import one_way_anova

    results = []

    coverage_groups = [data["coverage_per_replicate"] for _, data in condition_data]
    coverage_testable = all(len(group) >= 2 for group in coverage_groups)
    anova_coverage = one_way_anova(*coverage_groups)
    results.append(
        ContactsANOVASummary(
            metric="coverage",
            f_statistic=anova_coverage.f_statistic,
            p_value=anova_coverage.p_value,
            significant=anova_coverage.significant if coverage_testable else False,
            testable=coverage_testable,
            note=None if coverage_testable else NOT_TESTABLE_SINGLETON_NOTE,
        )
    )

    contact_groups = [data["contact_fraction_per_replicate"] for _, data in condition_data]
    contact_testable = all(len(group) >= 2 for group in contact_groups)
    anova_contact = one_way_anova(*contact_groups)
    results.append(
        ContactsANOVASummary(
            metric="mean_contact_fraction",
            f_statistic=anova_contact.f_statistic,
            p_value=anova_contact.p_value,
            significant=anova_contact.significant if contact_testable else False,
            testable=contact_testable,
            note=None if contact_testable else NOT_TESTABLE_SINGLETON_NOTE,
        )
    )

    return results


def apply_fdr_correction(
    comparisons: list[Any], anova_results: list[Any], fdr_alpha: float
) -> None:
    """Apply Benjamini-Hochberg FDR correction to comparison p-values.

    Parameters
    ----------
    comparisons : list[Any]
        Pairwise comparison results.
    anova_results : list[Any]
        ANOVA summaries.
    fdr_alpha : float
        False discovery rate alpha.
    """

    from polyzymd.analyses.shared.inferential_statistics import benjamini_hochberg

    all_pairwise_agg = []
    for comp in comparisons:
        all_pairwise_agg.extend(comp.aggregate_comparisons)

    if all_pairwise_agg:
        logger.debug(
            "Starting BH correction for contacts pairwise family: size=%d, alpha=%.4f",
            len(all_pairwise_agg),
            fdr_alpha,
        )
        raw_p = [agg.p_value if agg.testable else None for agg in all_pairwise_agg]
        bh_results = benjamini_hochberg(raw_p, alpha=fdr_alpha)
        changed_significance = 0
        for agg, bh in zip(all_pairwise_agg, bh_results, strict=False):
            if agg.significant != bh.significant:
                changed_significance += 1
            agg.p_value_adjusted = bh.adjusted_p_value
            agg.significant = bh.significant if agg.testable else False
        n_significant = sum(1 for agg in all_pairwise_agg if agg.significant)
        logger.info(
            "Applied BH correction to %d contacts pairwise tests at α=%.3f: "
            "%d remain significant, %d changed significance",
            len(all_pairwise_agg),
            fdr_alpha,
            n_significant,
            changed_significance,
        )

    if anova_results:
        logger.debug(
            "Starting BH correction for contacts ANOVA family: size=%d, alpha=%.4f",
            len(anova_results),
            fdr_alpha,
        )
        raw_p = [a.p_value if a.testable else None for a in anova_results]
        bh_results = benjamini_hochberg(raw_p, alpha=fdr_alpha)
        changed_significance = 0
        for anova, bh in zip(anova_results, bh_results, strict=False):
            if anova.significant != bh.significant:
                changed_significance += 1
            anova.p_value_adjusted = bh.adjusted_p_value
            anova.significant = bh.significant if anova.testable else False
        n_significant = sum(1 for anova in anova_results if anova.significant)
        logger.info(
            "Applied BH correction to %d contacts ANOVA tests at α=%.3f: "
            "%d remain significant, %d changed significance",
            len(anova_results),
            fdr_alpha,
            n_significant,
            changed_significance,
        )


def apply_effect_size_threshold(comparisons: list[Any], min_effect_size: float) -> None:
    """Tag aggregate comparisons by effect-size threshold.

    Parameters
    ----------
    comparisons : list[Any]
        Pairwise comparison results.
    min_effect_size : float
        Minimum absolute Cohen's d threshold.
    """

    met_threshold = 0
    failed_threshold = 0
    for comp in comparisons:
        for agg in comp.aggregate_comparisons:
            agg.meets_effect_size_threshold = (
                bool(agg.testable) and abs(agg.cohens_d) >= min_effect_size
            )
            if agg.meets_effect_size_threshold:
                met_threshold += 1
            else:
                failed_threshold += 1
    logger.info(
        "Applied contacts effect-size threshold |d| >= %.3f: %d meet, %d fail",
        min_effect_size,
        met_threshold,
        failed_threshold,
    )


def compute_top_contacted_residues(
    condition_data: list[tuple[Any, dict[str, Any]]], top_n: int
) -> dict[str, list[tuple[int, str, float]]] | None:
    """Compute top contacted residues per condition.

    Parameters
    ----------
    condition_data : list[tuple[Any, dict[str, Any]]]
        Condition data with aggregated contacts results.
    top_n : int
        Number of top residues to include. Values less than one disable output.

    Returns
    -------
    dict[str, list[tuple[int, str, float]]] or None
        Top residue tuples by condition, or ``None`` when disabled.
    """

    if top_n <= 0:
        return None

    logger.debug(
        "Computing top contacted residues: top_n=%d across %d conditions",
        top_n,
        len(condition_data),
    )

    def _as_contact_fraction(residue_stat: Any) -> float:
        """Convert residue contact fraction to float for robust sorting."""

        if isinstance(residue_stat, Mapping):
            value = residue_stat.get("contact_fraction_mean", 0.0)
        else:
            value = getattr(residue_stat, "contact_fraction_mean", 0.0)
        try:
            return float(value)
        except (TypeError, ValueError):
            return 0.0

    result: dict[str, list[tuple[int, str, float]]] = {}
    conditions_with_residue_data: list[str] = []
    for cond, data in condition_data:
        agg_result = data["agg_result"]
        if isinstance(agg_result, ConditionArtifact):
            residue_stats = list(agg_result.payload.get("residue_stats", []))
        else:
            residue_stats = getattr(agg_result, "residue_stats", [])
        if residue_stats:
            conditions_with_residue_data.append(cond.label)
        sorted_residues = sorted(residue_stats, key=_as_contact_fraction, reverse=True)
        result[cond.label] = [
            (
                int(
                    rs.get("protein_resid", 0)
                    if isinstance(rs, Mapping)
                    else getattr(rs, "protein_resid", 0)
                ),
                str(
                    rs.get("protein_resname", "UNK")
                    if isinstance(rs, Mapping)
                    else getattr(rs, "protein_resname", "UNK")
                ),
                _as_contact_fraction(rs),
            )
            for rs in sorted_residues[:top_n]
        ]

    logger.info(
        "Computed top %d contacted residues for %d conditions with residue data: %s",
        top_n,
        len(conditions_with_residue_data),
        ", ".join(conditions_with_residue_data) if conditions_with_residue_data else "none",
    )

    return result
