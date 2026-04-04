"""Polymer affinity score analysis plugin — composite interaction scoring.

Quantifies total polymer-protein interaction strength by summing per-contact
free energy contributions weighted by the number of simultaneous contacts:

    S_{p,g} = N_{p,g} × ΔG_sel(p,g)

where:
    N_{p,g}   = mean_contact_fraction × n_exposed_in_group
    ΔG_sel(p,g) = -ln(contact_share / expected_share)  [kT]

The total score for a condition is:

    S = Σ_{p,g} S_{p,g}

Independence assumption
-----------------------
This formulation assumes contacts are thermodynamically independent — each
contact contributes the same free energy regardless of what other contacts
exist simultaneously.  The absolute values are NOT rigorous binding free
energies.  Only the *relative differences* between polymer compositions are
meaningful as a comparative scoring metric.

Sign convention
---------------
    S < 0  →  net favorable polymer-protein interaction
    S > 0  →  net unfavorable (avoidance dominates)
    S = 0  →  contacts match the surface-availability reference

This is a **comparator-only** analysis: ``compute_replicate()`` and
``aggregate()`` return ``None``.  All computation is orchestrated within
``compare()`` which loads cached binding preference data from the contacts
analysis layer.

Plugin contract
---------------
- ``name = "polymer_affinity"``
- ``aliases = ("pa",)``
- ``dependencies = ("contacts",)``
- ``min_replicates = 1``
- ``compare()`` is fully overridden (temperature-aware, pairwise on total scores)
- ``filter_conditions()`` excludes no-polymer conditions
- ``plot()`` calls ``_plot_affinity_stacked_bars()`` and ``_plot_affinity_group_bars()``
"""

from __future__ import annotations

import logging
import math
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator

import polyzymd.analyses.polymer_affinity._plot_settings as _plot_settings  # register plot settings  # noqa: F401
from polyzymd.analyses.base import (
    Analysis,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
)
from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    get_colors,
    get_output_path,
    grouped_bars,
    save_figure,
)

if TYPE_CHECKING:
    from polyzymd.analyses.polymer_affinity._comparison_results import (
        AffinityScoreConditionSummary,
        AffinityScoreEntry,
        AffinityScorePairwiseEntry,
        PolymerAffinityScoreResult,
        PolymerTypeScore,
    )
    from polyzymd.analyses.shared.binding_preference import (
        AggregatedBindingPreferenceEntry,
        AggregatedBindingPreferenceResult,
        BindingPreferenceResult,
    )

logger = logging.getLogger("polyzymd.analyses.polymer_affinity")


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class PolymerAffinitySettings(BaseModel):
    """Unified settings for polymer affinity score plugin.

    Merges fields from ``PolymerAffinityScoreSettings`` and
    ``PolymerAffinityScoreComparisonSettings`` into a single Pydantic model.

    Parameters
    ----------
    compute_binding_preference : bool
        When True, compute binding preference from contacts data if no cached
        result is found.
    surface_exposure_threshold : float
        Relative SASA threshold for surface exposure (0.0–1.0).
    enzyme_pdb_for_sasa : str | None
        Path to enzyme PDB for SASA calculation.
    include_default_aa_groups : bool
        Include default amino acid class groupings.
    protein_groups : dict | None
        Custom protein groups ``{name: [resid1, resid2, ...]}``.
    protein_partitions : dict | None
        Custom partitions for mutually-exclusive group comparison.
    polymer_type_selections : dict | None
        Custom polymer type MDAnalysis selections.
    polymer_chain : str
        Chain ID for polymer auto-detection when *polymer_type_selections*
        is None. Defaults to ``"C"`` (PolyzyMD chain convention).
    fdr_alpha : float
        FDR alpha for Benjamini-Hochberg correction.
    """

    compute_binding_preference: bool = True
    surface_exposure_threshold: float = 0.2
    enzyme_pdb_for_sasa: str | None = None
    include_default_aa_groups: bool = True
    protein_groups: dict[str, list[int]] | None = None
    protein_partitions: dict[str, list[str]] | None = None
    polymer_type_selections: dict[str, str] | None = None
    polymer_chain: str = "C"
    fdr_alpha: float = 0.05

    @field_validator("surface_exposure_threshold")
    @classmethod
    def _validate_threshold(cls, v: float) -> float:
        if not 0.0 <= v <= 1.0:
            raise ValueError("surface_exposure_threshold must be between 0 and 1")
        return v


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class PolymerAffinityAnalysis(Analysis):
    """Polymer affinity score analysis plugin.

    Computes a composite interaction score for each (polymer_type,
    protein_group) pair by multiplying the mean number of simultaneous
    contacts by the per-contact selectivity free energy:

        S = N × ΔG_sel   [kT]

    The total score is summed across all polymer types and protein groups.
    More negative = stronger net polymer-protein affinity.

    This is a **comparator-only** plugin: ``compute_replicate()`` and
    ``aggregate()`` return ``None``.  The full pipeline is orchestrated
    within ``compare()``.

    Statistical comparisons use per-replicate total scores and are only
    computed between conditions at the same simulation temperature.

    Class Attributes
    ----------------
    name : str
        ``"polymer_affinity"``
    Settings : type
        :class:`PolymerAffinitySettings`
    aliases : tuple[str, ...]
        ``("pa",)``
    dependencies : tuple[str, ...]
        ``("contacts",)``
    min_replicates : int
        ``1``
    """

    name: ClassVar[str] = "polymer_affinity"
    execution_cost_hint: ClassVar[str] = "high"
    Settings: ClassVar[type] = PolymerAffinitySettings
    aliases: ClassVar[tuple[str, ...]] = ("pa",)
    dependencies: ClassVar[tuple[str, ...]] = ("contacts",)
    min_replicates: ClassVar[int] = 1
    has_compute_stage: ClassVar[bool] = False
    has_aggregate_stage: ClassVar[bool] = False

    # === Compare (full override) ===

    def compare(self, ctx: ComparisonContext) -> Any | None:
        """Run polymer affinity score comparison across all conditions.

        Steps:
        1. Load binding preference data per condition (cached or computed)
        2. Compute affinity score entries from enrichments
        3. Detect temperature groups
        4. Pairwise total-score comparisons (stats suppressed cross-temperature)
        5. Build and return result

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided context.

        Returns
        -------
        PolymerAffinityScoreResult | None
            Comparison result, or ``None`` if no conditions have data.
        """
        from polyzymd import __version__
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScoreConditionSummary,
            PolymerAffinityScoreResult,
        )

        settings: PolymerAffinitySettings = ctx.settings
        logger.info(f"Starting polymer affinity score comparison: {ctx.name}")
        logger.info(f"  Conditions: {len(ctx.conditions)}")

        # Step 1: Build condition summaries
        condition_summaries: list[AffinityScoreConditionSummary] = []
        for cond in ctx.conditions:
            try:
                summary = self._build_condition_summary(cond, ctx, settings)
                condition_summaries.append(summary)
            except Exception as e:
                logger.warning(f"  Skipping condition '{cond.label}': {e}")

        if not condition_summaries:
            logger.warning("No binding preference data found for any condition")
            return None

        # Step 2: Collect metadata
        all_polymer_types: list[str] = sorted(
            {e.polymer_type for s in condition_summaries for e in s.entries}
        )
        all_protein_groups: list[str] = sorted(
            {e.protein_group for s in condition_summaries for e in s.entries}
        )

        # Step 3: Temperature grouping
        temp_groups: dict[float, list[str]] = {}
        for s in condition_summaries:
            temp_groups.setdefault(s.temperature_K, []).append(s.label)
        mixed_temperatures = len(temp_groups) > 1

        # Step 4: Pairwise comparisons on total scores
        pairwise = self._compute_pairwise(condition_summaries, ctx.effective_control)

        # Step 5: Build result
        surface_threshold = settings.surface_exposure_threshold

        return PolymerAffinityScoreResult(
            name=ctx.name,
            mixed_temperatures=mixed_temperatures,
            temperature_groups={str(k): v for k, v in temp_groups.items()},
            conditions=condition_summaries,
            pairwise_comparisons=pairwise,
            polymer_types=all_polymer_types,
            protein_groups=all_protein_groups,
            surface_exposure_threshold=surface_threshold,
            equilibration_time=ctx.equilibration or "",
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    # === Filter conditions ===

    def filter_conditions(
        self,
        conditions: list[Condition],
        settings: "BaseModel | None" = None,
    ) -> list[Condition]:
        """Exclude conditions without polymer (no affinity to score).

        No-polymer conditions (e.g. bare enzyme controls) have no polymer
        contacts and therefore no affinity data.  They are excluded from
        comparison but noted in logs.
        """
        filtered = [c for c in conditions if _condition_has_polymer(c)]
        excluded = len(conditions) - len(filtered)
        if excluded:
            logger.info(
                f"Excluded {excluded} no-polymer condition(s) from polymer affinity comparison"
            )
        return filtered

    # === Plot ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate polymer affinity score plots.

        Calls module-level private functions:
        - :func:`_plot_affinity_stacked_bars`
        - :func:`_plot_affinity_group_bars`
        """
        plots: list[Path] = []

        data: dict[str, Any] = {}
        labels: list[str] = []
        for cond in ctx.conditions:
            label = cond.label
            labels.append(label)
            analysis_dir = ctx.analysis_dirs.get(label)
            if analysis_dir is not None:
                data[label] = {
                    "analysis_dir": analysis_dir,
                    "aggregated_dir": analysis_dir / "aggregated",
                    "replicates": list(cond.replicates),
                }

        if not labels:
            return plots

        data["__meta__"] = {"results_dir": ctx.results_dir}
        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        plot_settings = ctx.plot_settings

        try:
            result = _plot_affinity_stacked_bars(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Affinity stacked bars plot failed: {exc}")

        try:
            result = _plot_affinity_group_bars(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Affinity group bars plot failed: {exc}")

        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format polymer affinity results without legacy dispatch."""
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        return format_affinity_result(result, format=self._normalize_output_format(output_format))

    # === extract_metrics (empty — full compare() override) ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Return empty dict — polymer affinity uses full compare() override."""
        return {}

    @staticmethod
    def _normalize_output_format(output_format: str) -> str:
        return "table" if output_format == "text" else output_format

    # === Private helpers ===

    def _build_condition_summary(
        self,
        cond: Condition,
        ctx: ComparisonContext,
        settings: PolymerAffinitySettings,
    ) -> Any:
        """Build affinity score entries for one condition.

        Data flow:
        1. Load binding preference (cached or computed)
        2. Use ``bp_result.entries`` (``AggregatedBindingPreferenceEntry`` list)
           for ``mean_contact_fraction`` and ``n_exposed_in_group``
        3. Also try to load per-replicate files for exact per-replicate scores
        4. Aggregate into per-polymer-type scores and total condition score

        Parameters
        ----------
        cond : Condition
            Condition being processed.
        ctx : ComparisonContext
            Comparison context.
        settings : PolymerAffinitySettings
            Plugin settings.

        Returns
        -------
        AffinityScoreConditionSummary
        """
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScoreConditionSummary,
        )

        logger.info(f"  Processing condition: {cond.label}")

        # Get temperature from sim config
        temperature_K = self._get_temperature(cond)

        # Load binding preference
        bp_result = self._load_binding_preference(cond, ctx, settings)

        if bp_result is None:
            return AffinityScoreConditionSummary(
                label=cond.label,
                config_path=str(cond.config_path),
                temperature_K=temperature_K,
                n_replicates=0,
                entries=[],
                polymer_types=[],
                protein_groups=[],
            )

        # Try to load per-replicate binding preference files
        per_rep_data = self._load_per_replicate_entries(cond, ctx)

        # Compute affinity score entries
        entries = self._compute_affinity_entries(bp_result, temperature_K, per_rep_data)

        # Aggregate into per-polymer-type scores
        polymer_type_scores = self._aggregate_polymer_type_scores(entries)

        # Compute total condition score
        total_score = sum(pts.total_score for pts in polymer_type_scores)
        total_n_contacts = sum(pts.total_n_contacts for pts in polymer_type_scores)

        # Total per-replicate scores (sum across all polymer types)
        total_per_rep = self._compute_total_per_replicate_scores(polymer_type_scores)
        total_unc: float | None = None
        if len(total_per_rep) >= 2:
            total_unc = float(np.std(total_per_rep, ddof=1) / np.sqrt(len(total_per_rep)))

        polymer_types = sorted({e.polymer_type for e in entries})
        protein_groups = sorted({e.protein_group for e in entries})
        n_replicates = max((e.n_replicates for e in entries), default=0)

        return AffinityScoreConditionSummary(
            label=cond.label,
            config_path=str(cond.config_path),
            temperature_K=temperature_K,
            n_replicates=n_replicates,
            total_score=total_score,
            total_score_uncertainty=total_unc,
            total_score_per_replicate=total_per_rep,
            total_n_contacts=total_n_contacts,
            polymer_type_scores=polymer_type_scores,
            entries=entries,
            polymer_types=polymer_types,
            protein_groups=protein_groups,
        )

    def _get_temperature(self, cond: Condition) -> float:
        """Extract simulation temperature from condition's sim_config.

        Parameters
        ----------
        cond : Condition
            Condition with sim_config.

        Returns
        -------
        float
            Temperature in Kelvin.
        """
        sim_config = cond.sim_config
        return float(sim_config.thermodynamics.temperature)

    def _load_binding_preference(
        self,
        cond: Condition,
        ctx: ComparisonContext,
        settings: PolymerAffinitySettings,
    ) -> Any:
        """Load or compute binding preference for a condition.

        Attempts to load cached aggregated or single-replicate binding
        preference results from the contacts analysis directory.

        Parameters
        ----------
        cond : Condition
            Condition being processed.
        ctx : ComparisonContext
            Comparison context with analysis_dirs.
        settings : PolymerAffinitySettings
            Plugin settings.

        Returns
        -------
        AggregatedBindingPreferenceResult | BindingPreferenceResult | None
        """
        from polyzymd.analyses.shared.binding_preference_helpers import (
            CondProxy,
            compute_condition_binding_preference,
            resolve_enzyme_pdb,
            try_load_cached_binding_preference,
        )

        # Resolve contacts analysis dir (sibling of polymer_affinity analysis dir)
        pa_analysis_dir = ctx.analysis_dirs.get(cond.label)
        contacts_analysis_dir: Path | None = None
        if pa_analysis_dir is not None:
            contacts_analysis_dir = pa_analysis_dir.parent / "contacts"

        # Build a minimal cond-like object for the helpers
        cond_proxy = CondProxy(label=cond.label, config=str(cond.config_path))

        # Try cached first
        if not ctx.recompute and contacts_analysis_dir is not None:
            bp = try_load_cached_binding_preference(cond_proxy, contacts_analysis_dir)
            if bp is not None:
                logger.info(f"    Loaded cached binding preference for {cond.label}")
                return bp

        # Compute if enabled
        if not settings.compute_binding_preference:
            logger.warning(
                f"    No cached binding preference for '{cond.label}' and "
                "compute_binding_preference=False"
            )
            return None

        logger.info(f"    Computing binding preference for {cond.label}...")
        sim_config = cond.sim_config

        # Resolve enzyme PDB
        enzyme_pdb = resolve_enzyme_pdb(
            enzyme_pdb_setting=settings.enzyme_pdb_for_sasa,
            source_path=cond.config_path,
            sim_config=sim_config,
        )

        if enzyme_pdb is None or not enzyme_pdb.exists():
            logger.warning(
                f"    Cannot compute binding preference for '{cond.label}': "
                "enzyme PDB not found. Set enzyme_pdb_for_sasa in settings."
            )
            return None

        if contacts_analysis_dir is None:
            logger.warning(
                f"    Cannot compute binding preference for '{cond.label}': "
                "no contacts analysis directory"
            )
            return None

        bp = compute_condition_binding_preference(
            cond=cond_proxy,
            sim_config=sim_config,
            analysis_dir=contacts_analysis_dir,
            enzyme_pdb=enzyme_pdb,
            threshold=settings.surface_exposure_threshold,
            include_default_aa_groups=settings.include_default_aa_groups,
            custom_protein_groups=settings.protein_groups,
            protein_partitions=settings.protein_partitions,
            polymer_type_selections=settings.polymer_type_selections,
            polymer_chain=settings.polymer_chain,
        )

        if bp is not None:
            logger.info(f"    Computed binding preference for {cond.label}")
        else:
            logger.warning(f"    Failed to compute binding preference for '{cond.label}'")

        return bp

    def _load_per_replicate_entries(
        self,
        cond: Condition,
        ctx: ComparisonContext,
    ) -> dict[int, list[Any]] | None:
        """Load per-replicate BindingPreferenceEntry objects.

        Attempts to load ``binding_preference_rep{N}.json`` files from the
        contacts analysis directory.  Returns None if files are not available.

        Parameters
        ----------
        cond : Condition
            Condition being processed.
        ctx : ComparisonContext
            Comparison context with analysis_dirs.

        Returns
        -------
        dict[int, list] or None
            Mapping of replicate number to list of ``BindingPreferenceEntry``
            objects, or None if files unavailable.
        """
        from polyzymd.analyses.shared.binding_preference import (
            BindingPreferenceResult,
        )

        pa_analysis_dir = ctx.analysis_dirs.get(cond.label)
        if pa_analysis_dir is None:
            return None

        contacts_analysis_dir = pa_analysis_dir.parent / "contacts"
        if not contacts_analysis_dir.exists():
            return None

        per_rep: dict[int, list[Any]] = {}

        for rep in cond.replicates:
            rep_path = contacts_analysis_dir / f"binding_preference_rep{rep}.json"
            if rep_path.exists():
                try:
                    rep_result = BindingPreferenceResult.load(rep_path)
                    per_rep[rep] = rep_result.entries
                except Exception as exc:
                    logger.debug(f"Failed to load per-replicate BP for rep{rep}: {exc}")

        return per_rep if per_rep else None

    def _compute_affinity_entries(
        self,
        bp_result: Any,
        temperature_K: float,
        per_rep_data: dict[int, list[Any]] | None,
    ) -> list[Any]:
        """Compute affinity score entries from binding preference data.

        Uses ``AggregatedBindingPreferenceEntry`` objects from
        ``bp_result.entries`` which have the ``mean_contact_fraction``
        and ``n_exposed_in_group`` fields needed for N_contacts.

        Parameters
        ----------
        bp_result : AggregatedBindingPreferenceResult | BindingPreferenceResult
            Binding preference data for one condition.
        temperature_K : float
            Simulation temperature in Kelvin.
        per_rep_data : dict or None
            Mapping of replicate number to list of per-replicate
            ``BindingPreferenceEntry`` objects, if loaded.

        Returns
        -------
        list[AffinityScoreEntry]
            One entry per (polymer_type, protein_group) pair.
        """
        from polyzymd.analyses.shared.binding_preference import (
            AggregatedBindingPreferenceResult,
        )

        if isinstance(bp_result, AggregatedBindingPreferenceResult):
            return self._entries_from_aggregated(bp_result, temperature_K, per_rep_data)
        else:
            return self._entries_from_single(bp_result, temperature_K)

    def _entries_from_aggregated(
        self,
        result: Any,
        temperature_K: float,
        per_rep_data: dict[int, list[Any]] | None,
    ) -> list[Any]:
        """Compute affinity entries from aggregated binding preference.

        Uses the group-level ``entries`` list which has
        ``mean_contact_fraction`` and ``n_exposed_in_group``.

        Parameters
        ----------
        result : AggregatedBindingPreferenceResult
            Multi-replicate aggregated result.
        temperature_K : float
            Simulation temperature.
        per_rep_data : dict or None
            Per-replicate BindingPreferenceEntry objects by replicate number.

        Returns
        -------
        list[AffinityScoreEntry]
        """
        entries: list[Any] = []

        for agg_entry in result.entries:
            entry = self._entry_from_agg_bp_entry(agg_entry, temperature_K, per_rep_data)
            if entry is not None:
                entries.append(entry)

        return entries

    def _entry_from_agg_bp_entry(
        self,
        agg_entry: Any,
        temperature_K: float,
        per_rep_data: dict[int, list[Any]] | None,
    ) -> Any | None:
        """Convert one AggregatedBindingPreferenceEntry to an AffinityScoreEntry.

        Parameters
        ----------
        agg_entry : AggregatedBindingPreferenceEntry
            Aggregated binding preference entry for one
            (polymer_type, protein_group) pair.
        temperature_K : float
            Simulation temperature.
        per_rep_data : dict or None
            Per-replicate data if available.

        Returns
        -------
        AffinityScoreEntry or None
            Affinity score entry, or None if data is insufficient.
        """
        from polyzymd.analyses.polymer_affinity._comparison_results import AffinityScoreEntry

        mcf = agg_entry.mean_contact_fraction
        n_exposed = agg_entry.n_exposed_in_group
        cs = agg_entry.mean_contact_share
        es = agg_entry.expected_share
        polymer_type = agg_entry.polymer_type
        protein_group = agg_entry.protein_group

        # N_contacts = mean_contact_fraction × n_exposed_in_group
        n_contacts = mcf * n_exposed

        # ΔG_sel per contact = -ln(contact_share / expected_share) [kT]
        delta_g: float | None = None
        if cs > 0 and es > 0:
            delta_g = -math.log(cs / es)

        # Affinity score = N × ΔG_sel
        score: float | None = None
        if delta_g is not None:
            score = n_contacts * delta_g

        # Per-replicate scores
        score_per_rep: list[float] = []

        if per_rep_data is not None:
            # Use exact per-replicate N_rep × ΔG_sel,rep
            score_per_rep = self._per_rep_scores_from_files(
                per_rep_data, polymer_type, protein_group
            )
        elif agg_entry.per_replicate_enrichments:
            # Approximate: use mean N with per-replicate ΔG_sel
            for enrichment_rep in agg_entry.per_replicate_enrichments:
                ratio_rep = enrichment_rep + 1.0
                if ratio_rep > 0:
                    dg_rep = -math.log(ratio_rep)
                    score_per_rep.append(n_contacts * dg_rep)
                else:
                    score_per_rep.append(float("nan"))

        # Clean NaN values
        valid_reps = [v for v in score_per_rep if not math.isnan(v)]

        # Uncertainty
        score_unc: float | None = None
        if len(valid_reps) >= 2:
            score_unc = float(np.std(valid_reps, ddof=1) / np.sqrt(len(valid_reps)))
        elif score is not None and agg_entry.sem_contact_fraction is not None:
            # Analytical error propagation: σ(S) = √[(N·σ_ΔG)² + (ΔG·σ_N)²]
            sem_cs = getattr(agg_entry, "sem_contact_share", None)
            if delta_g is not None and sem_cs and cs > 0:
                sigma_n = agg_entry.sem_contact_fraction * n_exposed
                sigma_dg = sem_cs / cs  # σ(ln(ratio)) ≈ σ_cs/cs
                score_unc = math.sqrt((n_contacts * sigma_dg) ** 2 + (delta_g * sigma_n) ** 2)

        return AffinityScoreEntry(
            polymer_type=polymer_type,
            protein_group=protein_group,
            partition_name="aa_class",
            n_contacts=n_contacts,
            delta_G_per_contact=delta_g,
            affinity_score=score,
            affinity_score_uncertainty=score_unc,
            affinity_score_per_replicate=valid_reps,
            mean_contact_fraction=mcf,
            n_exposed_in_group=n_exposed,
            contact_share=cs,
            expected_share=es,
            temperature_K=temperature_K,
            n_replicates=len(valid_reps) if valid_reps else agg_entry.n_replicates,
        )

    def _entries_from_single(
        self,
        result: Any,
        temperature_K: float,
    ) -> list[Any]:
        """Compute affinity entries from a single-replicate result.

        Parameters
        ----------
        result : BindingPreferenceResult
            Single-replicate binding preference result.
        temperature_K : float
            Simulation temperature.

        Returns
        -------
        list[AffinityScoreEntry]
        """
        from polyzymd.analyses.polymer_affinity._comparison_results import AffinityScoreEntry

        entries: list[Any] = []

        for bp_entry in result.entries:
            mcf = bp_entry.mean_contact_fraction
            n_exposed = bp_entry.n_exposed_in_group
            cs = bp_entry.contact_share
            es = bp_entry.expected_share

            n_contacts = mcf * n_exposed

            delta_g: float | None = None
            if cs > 0 and es > 0:
                delta_g = -math.log(cs / es)

            score: float | None = None
            if delta_g is not None:
                score = n_contacts * delta_g

            entries.append(
                AffinityScoreEntry(
                    polymer_type=bp_entry.polymer_type,
                    protein_group=bp_entry.protein_group,
                    partition_name="aa_class",
                    n_contacts=n_contacts,
                    delta_G_per_contact=delta_g,
                    affinity_score=score,
                    affinity_score_uncertainty=None,
                    affinity_score_per_replicate=[score] if score is not None else [],
                    mean_contact_fraction=mcf,
                    n_exposed_in_group=n_exposed,
                    contact_share=cs,
                    expected_share=es,
                    temperature_K=temperature_K,
                    n_replicates=1,
                )
            )

        return entries

    # =========================================================================
    # Per-replicate data loading
    # =========================================================================

    def _per_rep_scores_from_files(
        self,
        per_rep_data: dict[int, list[Any]],
        polymer_type: str,
        protein_group: str,
    ) -> list[float]:
        """Compute per-replicate N×ΔG_sel from loaded per-replicate entries.

        Parameters
        ----------
        per_rep_data : dict
            Replicate number → list of BindingPreferenceEntry.
        polymer_type : str
            Polymer type to match.
        protein_group : str
            Protein group to match.

        Returns
        -------
        list[float]
            Per-replicate affinity scores.
        """
        scores: list[float] = []

        for _rep_num, rep_entries in sorted(per_rep_data.items()):
            matched = False
            for entry in rep_entries:
                if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                    mcf = entry.mean_contact_fraction
                    n_exposed = entry.n_exposed_in_group
                    n_rep = mcf * n_exposed

                    cs = entry.contact_share
                    es = entry.expected_share
                    if cs > 0 and es > 0:
                        dg_rep = -math.log(cs / es)
                        scores.append(n_rep * dg_rep)
                    else:
                        scores.append(0.0)
                    matched = True
                    break

            if not matched:
                # This (polymer_type, group) pair doesn't exist in this replicate
                scores.append(0.0)

        return scores

    # =========================================================================
    # Aggregation helpers
    # =========================================================================

    def _aggregate_polymer_type_scores(
        self,
        entries: list[Any],
    ) -> list[Any]:
        """Group entries by polymer type and compute per-type totals.

        Parameters
        ----------
        entries : list[AffinityScoreEntry]
            All (polymer_type, protein_group) entries for one condition.

        Returns
        -------
        list[PolymerTypeScore]
            One score per polymer type, sorted by total score ascending.
        """
        from polyzymd.analyses.polymer_affinity._comparison_results import PolymerTypeScore

        # Group entries by polymer type
        by_polymer: dict[str, list[Any]] = {}
        for e in entries:
            by_polymer.setdefault(e.polymer_type, []).append(e)

        scores: list[Any] = []

        for polymer_type, group_entries in sorted(by_polymer.items()):
            total = sum(e.affinity_score for e in group_entries if e.affinity_score is not None)
            total_n = sum(e.n_contacts for e in group_entries)

            # Per-replicate total scores for this polymer type
            per_rep_totals = self._sum_per_replicate_across_groups(group_entries)
            total_unc: float | None = None
            if len(per_rep_totals) >= 2:
                total_unc = float(np.std(per_rep_totals, ddof=1) / np.sqrt(len(per_rep_totals)))

            scores.append(
                PolymerTypeScore(
                    polymer_type=polymer_type,
                    total_score=total,
                    total_score_uncertainty=total_unc,
                    total_score_per_replicate=per_rep_totals,
                    total_n_contacts=total_n,
                    group_contributions=group_entries,
                )
            )

        return sorted(scores, key=lambda s: s.total_score)

    def _sum_per_replicate_across_groups(
        self,
        entries: list[Any],
    ) -> list[float]:
        """Sum per-replicate scores across protein groups for one polymer type.

        Parameters
        ----------
        entries : list[AffinityScoreEntry]
            Entries for a single polymer type.

        Returns
        -------
        list[float]
            Per-replicate total scores (sum across groups).  Empty if any
            entry lacks per-replicate data.
        """
        if not entries:
            return []

        # Check all entries have the same number of replicates
        rep_counts = [len(e.affinity_score_per_replicate) for e in entries]
        if not rep_counts or min(rep_counts) == 0:
            return []

        n_reps = min(rep_counts)
        totals: list[float] = []

        for rep_idx in range(n_reps):
            rep_total = sum(
                e.affinity_score_per_replicate[rep_idx]
                for e in entries
                if rep_idx < len(e.affinity_score_per_replicate)
            )
            totals.append(rep_total)

        return totals

    def _compute_total_per_replicate_scores(
        self,
        polymer_type_scores: list[Any],
    ) -> list[float]:
        """Sum per-replicate scores across all polymer types.

        Parameters
        ----------
        polymer_type_scores : list[PolymerTypeScore]
            Per-polymer-type scores.

        Returns
        -------
        list[float]
            Per-replicate grand total scores.
        """
        if not polymer_type_scores:
            return []

        rep_counts = [len(pts.total_score_per_replicate) for pts in polymer_type_scores]
        if not rep_counts or min(rep_counts) == 0:
            return []

        n_reps = min(rep_counts)
        totals: list[float] = []

        for rep_idx in range(n_reps):
            rep_total = sum(
                pts.total_score_per_replicate[rep_idx]
                for pts in polymer_type_scores
                if rep_idx < len(pts.total_score_per_replicate)
            )
            totals.append(rep_total)

        return totals

    # =========================================================================
    # Pairwise comparison
    # =========================================================================

    def _compute_pairwise(
        self,
        summaries: list[Any],
        effective_control: str | None,
    ) -> list[Any]:
        """Compute pairwise total score comparisons with temperature grouping.

        Unlike the BFE plugin which compares per-(polymer, group) pairs,
        this plugin compares total scores between conditions.  This is
        the natural unit for the polymer affinity score since it sums across
        all interactions.

        Parameters
        ----------
        summaries : list[AffinityScoreConditionSummary]
            Condition summaries.
        effective_control : str | None
            Control label (or None).

        Returns
        -------
        list[AffinityScorePairwiseEntry]
        """
        label_to_summary = {s.label: s for s in summaries}
        labels = [s.label for s in summaries]
        comparisons: list[Any] = []

        # Check if control has usable data
        control_has_data = (
            effective_control is not None
            and effective_control in label_to_summary
            and len(label_to_summary[effective_control].entries) > 0
        )

        if control_has_data:
            summary_a = label_to_summary[effective_control]
            for label_b in labels:
                if label_b == effective_control:
                    continue
                summary_b = label_to_summary[label_b]
                if not summary_b.entries:
                    continue
                pw = self._compare_total_scores(summary_a, summary_b)
                comparisons.append(pw)
        else:
            # Fall back to all-pairs among conditions with data
            if effective_control is not None and effective_control in label_to_summary:
                logger.info(
                    f"Control '{effective_control}' has no affinity data (e.g. no polymer "
                    "contacts). Falling back to all-pairs comparison."
                )
            valid_labels = [lbl for lbl in labels if label_to_summary[lbl].entries]
            for i in range(len(valid_labels)):
                for j in range(i + 1, len(valid_labels)):
                    sa = label_to_summary[valid_labels[i]]
                    sb = label_to_summary[valid_labels[j]]
                    comparisons.append(self._compare_total_scores(sa, sb))

        return comparisons

    def _compare_total_scores(
        self,
        summary_a: Any,
        summary_b: Any,
    ) -> Any:
        """Compare total affinity scores between two conditions.

        Parameters
        ----------
        summary_a : AffinityScoreConditionSummary
            First condition (typically control/reference).
        summary_b : AffinityScoreConditionSummary
            Second condition.

        Returns
        -------
        AffinityScorePairwiseEntry
        """
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScorePairwiseEntry,
        )
        from polyzymd.compare.statistics import independent_ttest

        cross_temperature = not math.isclose(
            summary_a.temperature_K, summary_b.temperature_K, rel_tol=1e-3
        )

        score_a = summary_a.total_score
        score_b = summary_b.total_score
        delta = score_b - score_a

        # Stats: use per-replicate totals when enough replicates are available
        t_stat: float | None = None
        p_val: float | None = None

        if not cross_temperature:
            reps_a = [v for v in summary_a.total_score_per_replicate if not math.isnan(v)]
            reps_b = [v for v in summary_b.total_score_per_replicate if not math.isnan(v)]
            if len(reps_a) >= 2 and len(reps_b) >= 2:
                try:
                    ttest = independent_ttest(reps_a, reps_b)
                    t_stat = ttest.t_statistic
                    p_val = ttest.p_value
                except Exception as exc:
                    logger.debug(f"T-test failed for {summary_a.label} vs {summary_b.label}: {exc}")

        return AffinityScorePairwiseEntry(
            condition_a=summary_a.label,
            condition_b=summary_b.label,
            temperature_a_K=summary_a.temperature_K,
            temperature_b_K=summary_b.temperature_K,
            cross_temperature=cross_temperature,
            score_a=score_a,
            score_b=score_b,
            delta_score=delta,
            t_statistic=t_stat,
            p_value=p_val,
        )


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def _condition_has_polymer(cond: Condition, polymer_selection: str = "chainID C") -> bool:
    """Check whether a condition's simulation includes polymer chains.

    Uses three detection strategies in order:

    1. ``sim_config.polymers`` — fast check for polymer builder section.
    2. ``sim_config.topology.chains`` — check for chain C in topology model.
    3. MDAnalysis topology inspection — uses ``TrajectoryLoader`` to find
       the topology file and checks for atoms matching *polymer_selection*.

    Parameters
    ----------
    cond : Condition
        Condition to check.
    polymer_selection : str
        MDAnalysis selection string for polymer atoms.  Defaults to
        ``"chainID C"`` (PolyzyMD chain convention).

    Returns
    -------
    bool
        True if the simulation config mentions polymer chains.
    """
    sim_config = cond.sim_config

    # Check 1: polymers builder section exists, is non-None, and enabled
    polymers = getattr(sim_config, "polymers", None)
    if polymers is not None:
        enabled = getattr(polymers, "enabled", True)
        if enabled:
            return True

    # Check 2: chain C in topology model
    if hasattr(sim_config, "topology"):
        topo = sim_config.topology
        if hasattr(topo, "chains") and topo.chains:
            chain_ids = [c.chain_id if hasattr(c, "chain_id") else c for c in topo.chains]
            if "C" in chain_ids:
                return True

    # Check 3: MDAnalysis topology inspection (same approach as contacts plugin)
    try:
        for rep in cond.replicates:
            run_dir = sim_config.get_working_directory(rep)

            # Use TrajectoryLoader's robust topology search rather than
            # hardcoding "solvated_system.pdb".
            try:
                from polyzymd.analyses.shared.loader import TrajectoryLoader

                loader = TrajectoryLoader(sim_config)
                topology_path = loader.find_topology(run_dir)
            except (FileNotFoundError, ImportError):
                continue

            try:
                import MDAnalysis as mda

                universe = mda.Universe(str(topology_path))
                polymer_atoms = universe.select_atoms(polymer_selection)
                if len(polymer_atoms) > 0:
                    logger.debug(f"  {cond.label} rep {rep}: {len(polymer_atoms)} polymer atoms")
                    return True
                else:
                    logger.debug(f"  {cond.label} rep {rep}: 0 polymer atoms")
            except Exception as e:
                logger.warning(f"  Error checking {cond.label} rep {rep}: {e}")
                continue
    except (AttributeError, TypeError):
        # sim_config may not have get_working_directory (e.g. in tests)
        pass

    return False


# ===========================================================================
# Inlined plotting functions (from compare/plotters/polymer_affinity.py)
# ===========================================================================


def _find_affinity_result(data: dict[str, Any], labels: Sequence[str]) -> Any | None:
    """Find and load PolymerAffinityScoreResult from the comparison cache.

    Parameters
    ----------
    data : dict
        Condition data dict with optional ``"__meta__"`` key.
    labels : sequence of str
        Condition labels in display order.

    Returns
    -------
    PolymerAffinityScoreResult or None
    """
    from polyzymd.analyses.polymer_affinity._comparison_results import PolymerAffinityScoreResult
    from polyzymd.compare.io.results import find_comparison_result

    return find_comparison_result(
        data,
        labels,
        glob_patterns=[
            "polymer_affinity_comparison_*.json",
            "affinity_comparison_*.json",
        ],
        loader=PolymerAffinityScoreResult.load,
        analysis_type="polymer_affinity",
        log=logger,
    )


def _plot_affinity_stacked_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Stacked bar chart of total affinity score per condition.

    Each bar represents one condition's total affinity score, with segments
    colored by polymer type contribution.

    Parameters
    ----------
    data : dict
        Condition data dict from orchestrator.
    labels : sequence of str
        Condition labels in display order.
    output_dir : Path
        Directory to save plot files.
    plot_settings : PlotSettings
        Plot configuration.

    Returns
    -------
    list[Path]
        Paths to generated plot files.
    """
    import matplotlib.pyplot as plt

    t = plot_settings.theme

    result = _find_affinity_result(data, labels)
    if result is None:
        return []

    affinity_settings = plot_settings.polymer_affinity

    # Order conditions by labels
    cond_labels = [c.label for c in result.conditions]
    display_labels = [lbl for lbl in labels if lbl in cond_labels]
    if not display_labels:
        display_labels = cond_labels

    # Collect polymer types across all conditions
    all_polymer_types = result.polymer_types
    if not all_polymer_types:
        logger.info("No polymer types in affinity result — skipping stacked bars")
        return []

    n_conds = len(display_labels)
    n_poly = len(all_polymer_types)
    colors = get_colors(n_poly, plot_settings)

    fig, ax = plt.subplots(figsize=affinity_settings.figsize_stacked, dpi=plot_settings.dpi)

    x = np.arange(n_conds)
    bottoms_neg = np.zeros(n_conds)
    bottoms_pos = np.zeros(n_conds)

    for poly_idx, poly_type in enumerate(all_polymer_types):
        values = []
        for cond_label in display_labels:
            cond = result.get_condition(cond_label)
            if cond is None:
                values.append(0.0)
                continue
            # Find this polymer type's score
            pts = [s for s in cond.polymer_type_scores if s.polymer_type == poly_type]
            if pts:
                values.append(pts[0].total_score)
            else:
                values.append(0.0)

        vals = np.array(values)

        # Stack negative bars downward, positive upward
        neg_vals = np.where(vals < 0, vals, 0)
        pos_vals = np.where(vals >= 0, vals, 0)

        if np.any(neg_vals != 0):
            ax.bar(
                x,
                neg_vals,
                bottom=bottoms_neg,
                color=colors[poly_idx],
                label=poly_type,
                alpha=t.bar_alpha,
                edgecolor="white",
                linewidth=t.bar_linewidth,
            )
            bottoms_neg += neg_vals

        if np.any(pos_vals != 0):
            ax.bar(
                x,
                pos_vals,
                bottom=bottoms_pos,
                color=colors[poly_idx],
                label=poly_type if np.all(neg_vals == 0) else None,
                alpha=t.bar_alpha,
                edgecolor="white",
                linewidth=t.bar_linewidth,
            )
            bottoms_pos += pos_vals

    # Add total score markers with uncertainty
    if affinity_settings.show_error_bars:
        totals = []
        errors = []
        for cond_label in display_labels:
            cond = result.get_condition(cond_label)
            if cond is not None:
                totals.append(cond.total_score)
                errors.append(cond.total_score_uncertainty if cond.total_score_uncertainty else 0.0)
            else:
                totals.append(0.0)
                errors.append(0.0)
        ax.errorbar(
            x,
            totals,
            yerr=errors,
            fmt="k_",
            capsize=t.bar_capsize,
            capthick=1.5,
            linewidth=0,
            elinewidth=1.5,
            label="Total ± SEM",
            zorder=10,
        )

    ax.axhline(y=0, color="black", linewidth=1.0, linestyle="-")
    ax.set_xticks(x)
    ax.set_xticklabels(display_labels, rotation=35, ha="right", fontsize=t.tick_fontsize)
    # Temperature string
    temp_str = ""
    if result.conditions:
        temps = {c.temperature_K for c in result.conditions}
        if len(temps) == 1:
            temp_str = f" ({next(iter(temps)):.0f} K)"

    apply_axis_style(
        ax,
        plot_settings,
        title=f"Polymer Affinity Score by Condition{temp_str}",
        ylabel=r"Affinity Score ($k_\mathrm{b}T$)",
    )

    # De-duplicate legend entries
    handles, legend_labels = ax.get_legend_handles_labels()
    seen: dict[str, Any] = {}
    unique_handles = []
    unique_labels = []
    for handle, lbl in zip(handles, legend_labels):
        if lbl not in seen:
            seen[lbl] = True
            unique_handles.append(handle)
            unique_labels.append(lbl)
    apply_legend(
        ax,
        plot_settings,
        fontsize=t.small_fontsize,
        handles=unique_handles,
        labels=unique_labels,
        framealpha=0.7,
    )

    plt.tight_layout()

    output_path = get_output_path(output_dir, "affinity_stacked_bars", plot_settings)
    return [
        save_figure(
            fig,
            output_path,
            plot_settings,
            experimental_features=("polymer_affinity",),
        )
    ]


def _plot_affinity_group_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Grouped bar chart of per-group affinity score contributions.

    Creates one figure per polymer type with groups on x-axis (AA classes),
    bars per condition, error bars, and reference line at score = 0.

    Parameters
    ----------
    data : dict
        Condition data dict from orchestrator.
    labels : sequence of str
        Condition labels in display order.
    output_dir : Path
        Directory to save plot files.
    plot_settings : PlotSettings
        Plot configuration.

    Returns
    -------
    list[Path]
        Paths to generated plot files.
    """
    import matplotlib.pyplot as plt

    from polyzymd.analyses.shared.aa_classification import CANONICAL_AA_CLASS_ORDER

    t = plot_settings.theme

    result = _find_affinity_result(data, labels)
    if result is None:
        return []

    affinity_settings = plot_settings.polymer_affinity

    cond_labels = [c.label for c in result.conditions]
    display_labels = [lbl for lbl in labels if lbl in cond_labels]
    if not display_labels:
        display_labels = cond_labels

    # Filter to conditions with data
    valid_labels = [
        lbl
        for lbl in display_labels
        if lbl in cond_labels
        and any(e.affinity_score is not None for e in result.get_condition(lbl).entries)
    ]
    if not valid_labels:
        logger.info("No conditions with affinity score data — skipping group bars")
        return []

    polymer_types = result.polymer_types
    protein_groups = result.protein_groups

    if not polymer_types or not protein_groups:
        return []

    # Sort protein groups canonically
    ordered_groups = [g for g in CANONICAL_AA_CLASS_ORDER if g in protein_groups]
    for g in sorted(protein_groups):
        if g not in ordered_groups:
            ordered_groups.append(g)

    n_conds = len(valid_labels)
    n_groups = len(ordered_groups)
    colors = get_colors(n_conds, plot_settings)
    n_poly = len(polymer_types)

    # Temperature string
    temp_str = ""
    if result.conditions:
        temps = {c.temperature_K for c in result.conditions}
        if len(temps) == 1:
            temp_str = f" ({next(iter(temps)):.0f} K)"

    output_paths: list[Path] = []

    for poly_type in polymer_types:
        fig, ax = plt.subplots(figsize=affinity_settings.figsize_group_bars, dpi=plot_settings.dpi)

        x = np.arange(n_groups)

        series: list[tuple[str, list[float], list[float]]] = []
        for cond_label in valid_labels:
            cond = result.get_condition(cond_label)
            means: list[float] = []
            sems: list[float] = []

            for group in ordered_groups:
                # Find matching entry
                entry = None
                if cond is not None:
                    for e in cond.entries:
                        if e.polymer_type == poly_type and e.protein_group == group:
                            entry = e
                            break

                if entry is not None and entry.affinity_score is not None:
                    means.append(entry.affinity_score)
                    # Prefer replicate-based SEM
                    per_rep = entry.affinity_score_per_replicate
                    if len(per_rep) >= 2:
                        sem = float(np.std(per_rep, ddof=1) / np.sqrt(len(per_rep)))
                    elif entry.affinity_score_uncertainty is not None:
                        sem = entry.affinity_score_uncertainty
                    else:
                        sem = 0.0
                    sems.append(sem)
                else:
                    means.append(0.0)
                    sems.append(0.0)

            series.append((cond_label, means, sems))

        grouped_bars(
            ax,
            x,
            series,
            colors,
            plot_settings,
            show_error=affinity_settings.show_error_bars,
            reference_label="Score = 0 (neutral)",
            bar_edgecolor="none",
        )

        poly_label = f": {poly_type}" if n_poly > 1 else ""
        apply_axis_style(
            ax,
            plot_settings,
            title=f"Per-Group Affinity Score{poly_label}{temp_str}",
            xlabel="Amino Acid Group",
            ylabel=r"Affinity Score ($k_\mathrm{b}T$)",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(ordered_groups, rotation=35, ha="right", fontsize=t.tick_fontsize)
        apply_legend(
            ax,
            plot_settings,
            fontsize=t.small_fontsize,
            framealpha=0.7,
        )

        # Guide lines at ±1 kT
        ax.axhline(y=1.0, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
        ax.axhline(y=-1.0, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)

        plt.tight_layout()

        stem = f"affinity_group_bars_{poly_type.lower()}" if n_poly > 1 else "affinity_group_bars"
        output_path = get_output_path(output_dir, stem, plot_settings)
        output_paths.append(
            save_figure(
                fig,
                output_path,
                plot_settings,
                experimental_features=("polymer_affinity",),
            )
        )

    return output_paths
