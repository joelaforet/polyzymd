"""Comparison and statistics helpers for the contacts analysis facade."""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Any, Sequence

import numpy as np

from polyzymd.analyses.base import ComparisonContext, Condition

logger = logging.getLogger("polyzymd.analyses.contacts")


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
        agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
        agg_result = ctx.aggregated_results.get(cond.label)
        if agg_result is not None and not analysis._cache_matches_context(
            agg_result,
            settings=settings,
            equilibration=ctx.equilibration,
            sim_config=cond.sim_config,
            replicates=cond.replicates,
            allow_replicate_subset=True,
        ):
            logger.warning(f"Invalid in-memory aggregate for '{cond.label}' — reloading.")
            agg_result = None

        if agg_result is None:
            agg_result = analysis._load_validated_aggregated_result(
                agg_dir,
                settings=settings,
                equilibration=ctx.equilibration,
                replicates=cond.replicates,
                sim_config=cond.sim_config,
                recompute=False,
                allow_replicate_subset=True,
            )
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
        summary = ContactsConditionSummary(
            label=cond.label,
            config_path=str(cond.config_path),
            n_replicates=agg_result.n_replicates,
            n_residues=agg_result.n_residues,
            coverage_mean=agg_result.coverage_mean,
            coverage_sem=agg_result.coverage_sem,
            mean_contact_fraction=agg_result.mean_contact_fraction,
            mean_contact_fraction_sem=agg_result.mean_contact_fraction_sem,
            residence_time_by_polymer_type=agg_result.residence_time_by_polymer_type,
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
    first_resids = {rs.protein_resid for rs in first_result.residue_stats}

    for cond, data in condition_data[1:]:
        result = data["agg_result"]
        resids = {rs.protein_resid for rs in result.residue_stats}
        if resids != first_resids:
            missing_in_first = resids - first_resids
            missing_in_other = first_resids - resids
            raise ValueError(
                f"Residue set mismatch between '{first_cond.label}' and '{cond.label}'. "
                f"Missing in {first_cond.label}: {sorted(missing_in_first)}, "
                f"Missing in {cond.label}: {sorted(missing_in_other)}."
            )


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
            significant=ttest.significant,
            percent_change=pct,
            direction=direction,
        )
    )

    contact_a = data_a["contact_fraction_per_replicate"]
    contact_b = data_b["contact_fraction_per_replicate"]
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
            significant=ttest.significant,
            percent_change=pct,
            direction=direction,
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
    anova_coverage = one_way_anova(*coverage_groups)
    results.append(
        ContactsANOVASummary(
            metric="coverage",
            f_statistic=anova_coverage.f_statistic,
            p_value=anova_coverage.p_value,
            significant=anova_coverage.significant,
        )
    )

    contact_groups = [data["contact_fraction_per_replicate"] for _, data in condition_data]
    anova_contact = one_way_anova(*contact_groups)
    results.append(
        ContactsANOVASummary(
            metric="mean_contact_fraction",
            f_statistic=anova_contact.f_statistic,
            p_value=anova_contact.p_value,
            significant=anova_contact.significant,
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
        raw_p = [agg.p_value for agg in all_pairwise_agg]
        bh_results = benjamini_hochberg(raw_p, alpha=fdr_alpha)
        changed_significance = 0
        for agg, bh in zip(all_pairwise_agg, bh_results, strict=False):
            if agg.significant != bh.significant:
                changed_significance += 1
            agg.p_value_adjusted = bh.adjusted_p_value
            agg.significant = bh.significant
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
        raw_p = [a.p_value for a in anova_results]
        bh_results = benjamini_hochberg(raw_p, alpha=fdr_alpha)
        changed_significance = 0
        for anova, bh in zip(anova_results, bh_results, strict=False):
            if anova.significant != bh.significant:
                changed_significance += 1
            anova.p_value_adjusted = bh.adjusted_p_value
            anova.significant = bh.significant
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
            agg.meets_effect_size_threshold = abs(agg.cohens_d) >= min_effect_size
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

        value = getattr(residue_stat, "contact_fraction_mean", 0.0)
        try:
            return float(value)
        except (TypeError, ValueError):
            return 0.0

    result: dict[str, list[tuple[int, str, float]]] = {}
    conditions_with_residue_data: list[str] = []
    for cond, data in condition_data:
        agg_result = data["agg_result"]
        residue_stats = getattr(agg_result, "residue_stats", [])
        if residue_stats:
            conditions_with_residue_data.append(cond.label)
        sorted_residues = sorted(residue_stats, key=_as_contact_fraction, reverse=True)
        result[cond.label] = [
            (
                int(getattr(rs, "protein_resid", 0)),
                str(getattr(rs, "protein_resname", "UNK")),
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
