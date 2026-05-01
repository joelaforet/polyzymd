"""Contacts analysis plugin.

Computes polymer-protein contacts from MD trajectories using parallel
neighbour searching, aggregates per-residue contact fractions across
replicates, and performs cross-condition comparison with dual metrics
(coverage and mean contact fraction).

Contact computation uses :class:`ParallelContactAnalyzer` which delegates to
MDAnalysis ``capped_distance`` for O(N) neighbour searching — typically
10–100× faster than naïve pairwise distance calculations.

Unlike single-scalar analyses (RMSF, catalytic_triad), contacts has **two**
primary metrics — coverage (fraction of residues contacted) and mean
contact fraction (average per-residue contact fraction).  Therefore
``compare()`` is overridden entirely and ``extract_metrics()`` is not used.

Additional sub-pipeline: binding preference comparison is optionally
computed when ``compute_binding_preference=True`` in settings.

Condition filtering
-------------------
No-polymer conditions (e.g. "No Polymer" controls) are automatically
excluded via :meth:`filter_conditions`. Detection checks topology with the
active MDAnalysis polymer selection and does not treat stale contacts caches
as polymer evidence.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, ValidationError, field_validator, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.contacts import _cache as _contacts_cache
from polyzymd.analyses.contacts import _lifecycle as _contacts_lifecycle
from polyzymd.analyses.contacts._aggregator import AggregatedContactResult
from polyzymd.analyses.contacts._identity import (
    contacts_settings_fingerprint,
    contacts_settings_fingerprint_candidates,
)
from polyzymd.analyses.contacts._plot_settings import ContactsPlotSettings
from polyzymd.analyses.contacts._plotters import (
    _plot_binding_preference_bars,
    _plot_binding_preference_heatmap,
    _plot_cf_by_aa_class_bars,
    _plot_cf_by_partition_bars,
    _plot_contact_fraction_profile,
    _plot_residence_time_profile,
    _plot_rt_by_aa_class_bars,
    _plot_rt_by_partition_bars,
    _plot_system_coverage_bars,
    _plot_system_coverage_heatmap,
    _plot_user_partition_bars,
)
from polyzymd.analyses.contacts._results import ContactResult
from polyzymd.analyses.contacts._runner import (
    ContactsReplicateRunner,
    ParallelContactAnalyzer,
    build_contact_grouping,
)

if TYPE_CHECKING:
    from polyzymd.analyses.contacts._comparison_results import ContactsComparisonResult

logger = logging.getLogger("polyzymd.analyses.contacts")


# Default cutoff matching the existing settings module
DEFAULT_CONTACT_CUTOFF = 4.5


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class ContactsSettings(BaseModel):
    """Settings for contacts analysis.

    Unified settings model for the contacts analysis plugin.

    Attributes
    ----------
    polymer_selection : str
        MDAnalysis selection for polymer atoms.
    protein_selection : str
        MDAnalysis selection for protein atoms.
    cutoff : float
        Contact distance cutoff in Angstroms.
    polymer_types : list[str] | None
        Filter contacts by polymer residue names.
    grouping : str
        How to group protein residues: ``aa_class``, ``secondary_structure``,
        or ``none``.
    compute_residence_times : bool
        If ``True``, compute residence time statistics.
    compute_binding_preference : bool
        If ``True``, compute binding preference enrichment analysis.
    surface_exposure_threshold : float
        Relative SASA threshold for surface-exposed residues.
    enzyme_pdb_for_sasa : str | None
        Path to enzyme PDB for SASA calculation.
    include_default_aa_groups : bool
        Include default AA class groupings in binding preference.
    protein_groups : dict[str, list[int]] | None
        Custom protein groups as ``{name: [resid, ...]}``.
    protein_partitions : dict[str, list[str]] | None
        Custom partitions as ``{partition_name: [group1, ...]}``.
    polymer_type_selections : dict[str, str] | None
        Custom polymer type selections as ``{name: "MDAnalysis selection"}``.
    polymer_chain : str
        Chain ID for polymer auto-detection when *polymer_type_selections*
        is None. Defaults to ``"C"`` (PolyzyMD chain convention).
    fdr_alpha : float
        False discovery rate alpha for Benjamini-Hochberg correction.
    min_effect_size : float
        Minimum Cohen's d effect size to highlight.
    top_residues : int
        Number of top residues to display in console.
    enrichment_normalization : str
        **DEPRECATED** — kept for backward compatibility. Ignored.
    """

    # --- Analysis settings ---
    polymer_selection: str = Field(
        default="chainID C", description="MDAnalysis selection for polymer atoms"
    )
    protein_selection: str = Field(
        default="protein", description="MDAnalysis selection for protein atoms"
    )
    cutoff: float = Field(
        default=DEFAULT_CONTACT_CUTOFF,
        description="Contact distance cutoff in Angstroms",
    )
    polymer_types: list[str] | None = Field(
        default=None, description="Filter by polymer residue names"
    )
    grouping: str = Field(
        default="aa_class",
        description="Group by: aa_class, secondary_structure, or none",
    )
    compute_residence_times: bool = Field(
        default=True, description="Compute residence time statistics"
    )

    # --- Binding preference settings ---
    compute_binding_preference: bool = Field(
        default=False,
        description="Compute binding preference enrichment analysis",
    )
    surface_exposure_threshold: float = Field(
        default=0.2,
        description="Relative SASA threshold for surface-exposed residues",
    )
    enzyme_pdb_for_sasa: str | None = Field(
        default=None,
        description="Path to enzyme PDB for SASA calculation",
    )
    include_default_aa_groups: bool = Field(
        default=True,
        description="Include default AA class groupings",
    )
    protein_groups: dict[str, list[int]] | None = Field(
        default=None,
        description="Custom protein groups as {name: [resid, ...]}",
    )
    protein_partitions: dict[str, list[str]] | None = Field(
        default=None,
        description="Custom partitions as {partition_name: [group1, ...]}",
    )
    polymer_type_selections: dict[str, str] | None = Field(
        default=None,
        description="Custom polymer type selections as {name: 'MDAnalysis selection'}",
    )
    polymer_chain: str = Field(
        default="C",
        description=(
            "Chain ID for polymer auto-detection when polymer_type_selections "
            "is None. Defaults to 'C' (PolyzyMD chain convention)."
        ),
    )
    enrichment_normalization: str = Field(
        default="residue",
        description="DEPRECATED: ignored. Kept for backward compat.",
    )

    # --- Comparison settings ---
    fdr_alpha: float = Field(
        default=0.05,
        description="FDR alpha for Benjamini-Hochberg correction",
    )
    min_effect_size: float = Field(
        default=0.5,
        description="Minimum Cohen's d to highlight",
    )
    top_residues: int = Field(
        default=10,
        description="Number of top residues to display in console",
    )

    @field_validator("grouping", mode="after")
    @classmethod
    def validate_grouping(cls, v: str) -> str:
        """Validate grouping mode."""
        valid = {"aa_class", "secondary_structure", "none"}
        if v not in valid:
            raise ValueError(f"grouping must be one of {valid}, got '{v}'")
        return v

    @field_validator("fdr_alpha", mode="after")
    @classmethod
    def validate_fdr_alpha(cls, v: float) -> float:
        """Validate FDR alpha is in valid range."""
        if not 0 < v < 1:
            raise ValueError(f"fdr_alpha must be between 0 and 1, got {v}")
        return v

    @model_validator(mode="after")
    def validate_protein_partitions(self) -> "ContactsSettings":
        """Validate protein_partitions references and mutual exclusivity."""
        if not self.protein_partitions:
            return self
        if not self.protein_groups:
            raise ValueError("protein_partitions requires protein_groups to be defined.")
        for partition_name, group_names in self.protein_partitions.items():
            if not group_names:
                raise ValueError(f"Partition '{partition_name}' is empty.")
            for group_name in group_names:
                if group_name not in self.protein_groups:
                    available = ", ".join(sorted(self.protein_groups.keys()))
                    raise ValueError(
                        f"Partition '{partition_name}' references undefined "
                        f"group '{group_name}'. Available: {available}"
                    )
            # Detect overlapping groups within this partition
            seen_resids: dict[int, str] = {}
            for group_name in group_names:
                for resid in self.protein_groups[group_name]:
                    if resid in seen_resids:
                        raise ValueError(
                            f"Partition '{partition_name}' has overlapping groups: "
                            f"residue {resid} is in both '{seen_resids[resid]}' "
                            f"and '{group_name}'."
                        )
                    seen_resids[resid] = group_name
        return self


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class ContactsAnalysis(Analysis):
    """Contacts analysis: polymer-protein contacts from MD trajectories.

    This plugin delegates trajectory-native contact detection to the contacts
    runner seam, while keeping aggregation, comparison, plotting, and
    downstream orchestration in the plugin class.

    The ``compare()`` method is **fully overridden** because:

    - Two primary metrics: coverage and mean_contact_fraction.
    - Auto-exclusion of no-polymer conditions.
    - Optional binding preference sub-pipeline.
    - Residue set validation across conditions.

    Plots
    -----
    Generates 11 plot types via private module-level functions:

    - Contact fraction / residence time profiles
    - Grouped bar charts (by AA class, by partition)
    - System coverage bar / heatmap
    - Binding preference bar / heatmap
    """

    name: ClassVar[str] = "contacts"
    execution_cost_hint: ClassVar[str] = "high"
    Settings: ClassVar[type] = ContactsSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = ContactsPlotSettings
    AggregatedResultClass: ClassVar[type] = AggregatedContactResult
    ReplicateResultClass: ClassVar[type | None] = ContactResult
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    @staticmethod
    def _replicate_sidecar_path(
        output_dir: Path,
        settings: BaseModel,
        equilibration: str,
        replicate: int,
    ) -> Path:
        """Build the legacy per-replicate contacts side-output path."""

        return _contacts_cache.replicate_sidecar_path(
            output_dir, settings, equilibration, replicate
        )

    @staticmethod
    def _replicate_sidecar_candidates(
        output_dir: Path,
        settings: BaseModel,
        equilibration: str,
        replicate: int,
    ) -> tuple[Path, ...]:
        """Build accepted per-replicate contacts side-output paths."""

        return _contacts_cache.replicate_sidecar_candidates(
            output_dir, settings, equilibration, replicate
        )

    @staticmethod
    def _replicate_cache_candidates(
        output_dir: Path,
        settings: BaseModel,
        equilibration: str,
        replicate: int,
        *,
        result_path: Path | None = None,
    ) -> tuple[Path, ...]:
        """Build accepted per-replicate contacts cache candidates."""

        return _contacts_cache.replicate_cache_candidates(
            output_dir,
            settings,
            equilibration,
            replicate,
            result_path=result_path,
        )

    @staticmethod
    def _aggregated_sidecar_path(
        output_dir: Path,
        settings: BaseModel,
        equilibration: str,
        replicates: Sequence[int],
    ) -> Path:
        """Build the legacy aggregated contacts side-output path."""

        return _contacts_cache.aggregated_sidecar_path(
            output_dir, settings, equilibration, replicates
        )

    @staticmethod
    def _aggregated_sidecar_candidates(
        output_dir: Path,
        settings: BaseModel,
        equilibration: str,
        replicates: Sequence[int],
    ) -> tuple[Path, ...]:
        """Build accepted aggregated contacts side-output paths."""

        return _contacts_cache.aggregated_sidecar_candidates(
            output_dir, settings, equilibration, replicates
        )

    @staticmethod
    def _aggregated_cache_candidates(
        output_dir: Path,
        settings: BaseModel,
        equilibration: str,
        replicates: Sequence[int],
        *,
        result_path: Path | None = None,
    ) -> tuple[Path, ...]:
        """Build accepted aggregated contacts cache candidates."""

        return _contacts_cache.aggregated_cache_candidates(
            output_dir,
            settings,
            equilibration,
            replicates,
            result_path=result_path,
        )

    @staticmethod
    def _effective_polymer_selection(settings: "ContactsSettings") -> str:
        """Return the polymer selection constrained by polymer type filters."""

        return _contacts_cache.effective_polymer_selection(settings)

    @staticmethod
    def _cache_matches_window(result: Any, equilibration: str) -> bool:
        """Return whether a cached result matches the requested analysis window."""

        return _contacts_cache.cache_matches_window(result, equilibration)

    def _load_cache_candidate(
        self,
        result_cls: type,
        candidate: Path,
        *,
        recompute: bool,
        sim_config: Any | None = None,
        settings: BaseModel | None = None,
    ) -> Any | None:
        """Load one cache candidate while tolerating invalid legacy JSON."""

        return _contacts_cache.load_cache_candidate(
            self,
            result_cls,
            candidate,
            recompute=recompute,
            sim_config=sim_config,
            settings=settings,
        )

    @staticmethod
    def _result_window_string(result: Any) -> str | None:
        """Return the equilibration window encoded on a result."""

        return _contacts_cache.result_window_string(result)

    @staticmethod
    def _cache_matches_replicate_id(
        result: Any,
        replicate: int,
        *,
        source: Path | None = None,
    ) -> bool:
        """Return whether a cached per-replicate result matches the request."""

        return _contacts_cache.cache_matches_replicate_id(result, replicate, source=source)

    @staticmethod
    def _embedded_contacts_settings_fingerprint(result: Any) -> str | None:
        """Return contacts settings fingerprint embedded in cache content."""

        return _contacts_cache.embedded_contacts_settings_fingerprint(result)

    @classmethod
    def _cache_has_contacts_identity_proof(cls, result: Any, settings: BaseModel) -> bool:
        """Return whether cache content fully proves contacts settings identity."""

        del cls
        return _contacts_cache.cache_has_contacts_identity_proof(result, settings)

    @staticmethod
    def _coerce_optional_bool(value: Any) -> bool | None:
        """Coerce common metadata values to booleans."""

        return _contacts_cache.coerce_optional_bool(value)

    @classmethod
    def _cache_matches_residence_time_setting(
        cls,
        result: Any,
        settings: BaseModel,
        *,
        allow_missing: bool,
        source: Path | None = None,
    ) -> bool:
        """Return whether a cache proves the active residence-time setting."""

        del cls
        return _contacts_cache.cache_matches_residence_time_setting(
            result,
            settings,
            allow_missing=allow_missing,
            source=source,
        )

    @staticmethod
    def _cache_has_residence_time_summaries(result: Any) -> bool:
        """Return whether an aggregate artifact contains RT summary maps."""

        return _contacts_cache.cache_has_residence_time_summaries(result)

    @classmethod
    def _attach_contacts_identity_metadata(cls, result: Any, settings: BaseModel) -> None:
        """Attach current contacts identity metadata to a cache result."""

        del cls
        _contacts_cache.attach_contacts_identity_metadata(result, settings)

    @staticmethod
    def _cache_matches_contacts_settings(
        result: Any,
        settings: BaseModel,
        *,
        source: Path | None = None,
    ) -> bool:
        """Return whether a cached contacts artifact matches active settings."""

        return _contacts_cache.cache_matches_contacts_settings(result, settings, source=source)

    @staticmethod
    def _align_replicate_results(
        results: Sequence[Any],
        replicates: Sequence[int],
    ) -> list[Any]:
        """Return replicate results ordered to match the requested IDs."""

        return _contacts_lifecycle.align_replicate_results(results, replicates)

    def _resolve_plot_equilibration(
        self,
        ctx: PlotContext,
        aggregated_dir: Path | None = None,
    ) -> str:
        """Resolve the equilibration window to use when validating plot inputs.

        Plot contexts created by the orchestrator carry the resolved window.
        Direct plot calls may omit it, so finalized comparison metadata is used
        first and canonical aggregate metadata is used only as a last-resort
        fallback for direct calls with the dataclass default.

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided plot context.
        aggregated_dir : Path or None, optional
            Aggregated result directory for direct-call fallback.

        Returns
        -------
        str
            Requested equilibration window for finalized aggregate validation.
        """

        comparison_path = ctx.comparison_path or self.comparison_result_path(ctx.results_dir)
        if comparison_path.exists():
            try:
                from polyzymd.analyses.contacts._comparison_results import ContactsComparisonResult

                comparison = ContactsComparisonResult.load(comparison_path)
            except (OSError, ValueError, ValidationError):
                logger.warning(
                    "contacts: could not load comparison metadata for plot validation at %s",
                    comparison_path,
                )
            else:
                return comparison.equilibration_time

        if ctx.equilibration != "0ns" or aggregated_dir is None:
            return ctx.equilibration

        aggregate_path = self.aggregate_result_path(aggregated_dir)
        if not aggregate_path.exists():
            return ctx.equilibration
        try:
            aggregate = AggregatedContactResult.load(aggregate_path)
        except (OSError, ValueError, ValidationError):
            return ctx.equilibration

        return self._result_window_string(aggregate) or ctx.equilibration

    def _cache_matches_context(
        self,
        result: Any,
        *,
        settings: BaseModel,
        equilibration: str,
        sim_config: Any,
        replicates: Sequence[int] | None = None,
        allow_replicate_subset: bool = False,
        source: Path | None = None,
    ) -> bool:
        """Return whether an aggregate result matches the active context."""

        return _contacts_cache.cache_matches_context(
            self,
            result,
            settings=settings,
            equilibration=equilibration,
            sim_config=sim_config,
            replicates=replicates,
            allow_replicate_subset=allow_replicate_subset,
            source=source,
        )

    def _cache_matches_replicates(
        self,
        result: Any,
        replicates: Sequence[int],
        *,
        allow_subset: bool = False,
        source: Path | None = None,
    ) -> bool:
        """Return whether cached aggregate replicate identity matches the request."""

        return _contacts_cache.cache_matches_replicates(
            self,
            result,
            replicates,
            allow_subset=allow_subset,
            source=source,
        )

    @staticmethod
    def _replicate_vector_length_mismatches(
        result: Any,
        expected_replicates: Sequence[int],
    ) -> list[str]:
        """Return per-replicate vector fields whose lengths are inconsistent."""

        return _contacts_cache.replicate_vector_length_mismatches(result, expected_replicates)

    def _load_validated_aggregated_result(
        self,
        aggregated_dir: Path,
        *,
        settings: BaseModel,
        equilibration: str,
        replicates: Sequence[int],
        sim_config: Any,
        recompute: bool,
        allow_replicate_subset: bool = False,
    ) -> Any | None:
        """Load an aggregated contacts result after cache identity validation."""

        return _contacts_cache.load_validated_aggregated_result(
            self,
            aggregated_dir,
            settings=settings,
            equilibration=equilibration,
            replicates=replicates,
            sim_config=sim_config,
            recompute=recompute,
            allow_replicate_subset=allow_replicate_subset,
        )

    def run_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Run contacts for a single replicate."""

        return _contacts_lifecycle.run_replicate(self, ctx, replicate)

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the contacts analysis window with a dt-aware fallback."""

        return _contacts_lifecycle.get_trajectory_window(ctx, replicate, loader, universe)

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any:
        """Build the runner-backed contacts execution object."""

        return _contacts_lifecycle.build_runner(self, ctx, replicate, universe, window)

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Attach framework metadata to the runner-produced contact result."""

        return _contacts_lifecycle.summarize_replicate(self, ctx, replicate, runner, window)

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate contact results across replicates for one condition."""

        return _contacts_lifecycle.aggregate(self, ctx, results)

    # === filter_conditions() — exclude no-polymer conditions ===

    def filter_conditions(
        self,
        conditions: list[Condition],
        settings: "BaseModel | None" = None,
    ) -> list[Condition]:
        """Filter conditions to only those with polymer atoms.

        Conditions without polymer atoms (e.g. "No Polymer" controls) are
        excluded since there are no polymer-protein contacts to analyse.

        Detection uses MDAnalysis topology inspection with the active
        ``polymer_selection``. Stale contacts caches are not trusted as
        polymer evidence.

        Parameters
        ----------
        conditions : list[Condition]
            All conditions from the comparison config.
        settings : BaseModel or None
            Resolved plugin settings from the orchestrator.

        Returns
        -------
        list[Condition]
            Conditions to include in analysis.
        """
        resolved = settings if isinstance(settings, self.Settings) else self.Settings()
        valid: list[Condition] = []

        for cond in conditions:
            try:
                if self._condition_has_polymer(cond, resolved):
                    valid.append(cond)
                else:
                    logger.info(
                        f"  Excluding '{cond.label}': no polymer atoms found "
                        f"with selection "
                        f"'{self._effective_polymer_selection(resolved)}'"
                    )
            except (AttributeError, ValueError, KeyError, OSError) as e:
                logger.warning(f"  Error checking condition '{cond.label}': {e} — including anyway")
                valid.append(cond)

        return valid

    # === compare() — fully overridden ===

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare contacts metrics across conditions.

        Dual primary metrics:

        - **coverage**: fraction of residues contacted (higher = better).
        - **mean_contact_fraction**: average per-residue contact fraction
          across all residues (higher = better).

        Additionally computes optional binding preference comparison
        when ``compute_binding_preference=True``.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided context.

        Returns
        -------
        ContactsComparisonResult | None
            Comparison result, or ``None`` when no conditions are available.
        """
        from polyzymd import __version__
        from polyzymd.analyses.contacts._comparison_results import (
            BindingPreferenceComparisonEntry,
            BindingPreferenceComparisonSummary,
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

        # Load aggregated results and build condition data
        condition_data: list[tuple[Condition, dict[str, Any]]] = []
        for cond in ctx.conditions:
            agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
            agg_result = ctx.aggregated_results.get(cond.label)
            if agg_result is not None and not self._cache_matches_context(
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
                # Comparison consumes finalize outputs even when recompute was requested
                agg_result = self._load_validated_aggregated_result(
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

            # Compute per-replicate values for statistical tests
            coverage_per_rep = self._compute_coverage_per_replicate(agg_result)
            contact_fraction_per_rep = self._compute_contact_fraction_per_replicate(agg_result)

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

        # Validate identical residue sets
        self._validate_residue_sets(condition_data)

        # Build condition summaries
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

        # Resolve effective control
        effective_control = self._resolve_effective_control(ctx.effective_control, summaries)

        # Compute pairwise comparisons for both aggregate metrics
        comparisons: list[ContactsPairwiseComparison] = []
        if len(summaries) >= 2:
            comparisons = self._compute_contacts_pairwise(
                summaries, condition_data, effective_control
            )

        # Compute ANOVA when at least three conditions are available
        anova_results: list[ContactsANOVASummary] = []
        if len(summaries) >= 3:
            anova_results = self._compute_contacts_anova(condition_data)

        # Apply BH FDR correction
        self._apply_fdr_correction(comparisons, anova_results, settings.fdr_alpha)

        # Apply effect size threshold
        self._apply_effect_size_threshold(comparisons, settings.min_effect_size)

        # Rank conditions by primary contact metrics
        ranked_coverage = sorted(summaries, key=lambda s: s.coverage_mean, reverse=True)
        ranked_contact = sorted(summaries, key=lambda s: s.mean_contact_fraction, reverse=True)

        # Collect top contacted residues for display
        top_residues_data = self._compute_top_contacted_residues(
            condition_data, settings.top_residues
        )

        # Add binding preference comparison when enabled
        binding_pref_summary = self._load_or_compute_binding_preference(ctx, condition_data)

        # Build comparison result
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
            binding_preference=binding_pref_summary,
            top_contacted_residues=top_residues_data,
            equilibration_time=ctx.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    # === plot() ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate contacts comparison plots.

        Calls 11 private module-level plotting functions covering:

        - Contact fraction / residence time profiles
        - Grouped bars by AA class and by partition
        - System coverage bar / heatmap
        - Binding preference bar / heatmap (if data exists)

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided context.

        Returns
        -------
        list[Path]
            Paths to generated figure files.
        """
        plots: list[Path] = []

        data, labels = self._build_plot_data(ctx, include_replicates=True)
        if not labels:
            return plots

        valid_labels: list[str] = []
        for cond in ctx.conditions:
            entry = data.get(cond.label)
            if entry is None:
                continue
            aggregated_dir = entry["aggregated_dir"]
            if not (aggregated_dir.exists() and any(aggregated_dir.glob("*.json"))):
                valid_labels.append(cond.label)
                continue
            equilibration = self._resolve_plot_equilibration(ctx, aggregated_dir)
            summary = self._load_validated_aggregated_result(
                aggregated_dir,
                settings=ctx.settings,
                equilibration=equilibration,
                replicates=cond.replicates,
                sim_config=cond.sim_config,
                recompute=False,
                allow_replicate_subset=True,
            )
            if summary is not None:
                valid_labels.append(cond.label)

        labels = [label for label in labels if label in valid_labels]
        data = {label: data[label] for label in labels} | {"__meta__": data["__meta__"]}
        if not labels:
            return plots

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve plot settings
        plot_settings = ctx.plot_settings

        # Residence-time plots are skipped when summaries are disabled
        plot_functions = [
            _plot_contact_fraction_profile,
            _plot_cf_by_aa_class_bars,
            _plot_cf_by_partition_bars,
            _plot_user_partition_bars,
            _plot_system_coverage_bars,
            _plot_system_coverage_heatmap,
            _plot_binding_preference_bars,
            _plot_binding_preference_heatmap,
        ]
        if ctx.settings.compute_residence_times:
            plot_functions[1:1] = [
                _plot_residence_time_profile,
                _plot_rt_by_aa_class_bars,
                _plot_rt_by_partition_bars,
            ]

        for plot_fn in plot_functions:
            try:
                result = plot_fn(data, labels, ctx.output_dir, plot_settings)
                if result:
                    plots.extend(result)
            except (ValueError, RuntimeError, OSError) as exc:
                fn_name = getattr(plot_fn, "__name__", repr(plot_fn))
                logger.warning(f"{fn_name} plot failed: {exc}")

        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format contacts comparison results without legacy dispatch."""
        from polyzymd.analyses.contacts._formatters import format_contacts_result

        return format_contacts_result(result, format=self._normalize_output_format(output_format))

    @staticmethod
    def _normalize_output_format(output_format: str) -> str:
        return "table" if output_format == "text" else output_format

    # === Private helpers: condition filtering ===

    def _condition_has_polymer(
        self, cond: Condition, settings: "ContactsAnalysis.Settings"
    ) -> bool:
        """Check whether a condition has polymer atoms.

        Detection uses topology inspection via MDAnalysis with the effective
        polymer selection.

        Parameters
        ----------
        cond : Condition
            Condition to check.
        settings : ContactsAnalysis.Settings
            Resolved settings with polymer query filters.

        Returns
        -------
        bool
            True if polymer atoms detected.
        """
        # Topology-based polymer detection
        for rep in cond.replicates:
            run_dir = cond.sim_config.get_working_directory(rep)

            # Use TrajectoryLoader's robust topology search rather than
            # hardcoding "solvated_system.pdb".
            try:
                from polyzymd.analyses.shared.loader import TrajectoryLoader

                loader = TrajectoryLoader(cond.sim_config)
                topo_path = loader.find_topology(run_dir)
            except (FileNotFoundError, ImportError):
                continue

            try:
                import MDAnalysis as mda

                universe = mda.Universe(str(topo_path))
                polymer_atoms = universe.select_atoms(self._effective_polymer_selection(settings))
                if len(polymer_atoms) > 0:
                    logger.debug(f"  {cond.label} rep {rep}: {len(polymer_atoms)} polymer atoms")
                    return True
                else:
                    logger.debug(f"  {cond.label} rep {rep}: 0 polymer atoms")
            except (AttributeError, ValueError, KeyError, OSError) as e:
                logger.warning(f"  Error checking {cond.label} rep {rep}: {e}")
                continue

        return False

    # === Private helpers: per-replicate metric computation ===

    @staticmethod
    def _compute_coverage_per_replicate(
        result: "AggregatedContactResult",
    ) -> list[float]:
        """Compute coverage per replicate from residue stats.

        Coverage = fraction of residues with any contact in a replicate.

        Parameters
        ----------
        result : AggregatedContactResult
            Aggregated result.

        Returns
        -------
        list[float]
            Coverage for each replicate.
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

    @staticmethod
    def _compute_contact_fraction_per_replicate(
        result: "AggregatedContactResult",
    ) -> list[float]:
        """Compute mean contact fraction per replicate.

        Parameters
        ----------
        result : AggregatedContactResult
            Aggregated result.

        Returns
        -------
        list[float]
            Mean contact fraction for each replicate.
        """
        n_replicates = result.n_replicates

        fractions = []
        for rep_idx in range(n_replicates):
            rep_fractions = [
                rs.contact_fraction_per_replicate[rep_idx] for rs in result.residue_stats
            ]
            mean_frac = float(np.mean(rep_fractions)) if rep_fractions else 0.0
            fractions.append(mean_frac)

        return fractions

    # === Private helpers: residue set validation ===

    @staticmethod
    def _validate_residue_sets(
        condition_data: list[tuple[Condition, dict[str, Any]]],
    ) -> None:
        """Validate that all conditions have identical residue sets.

        Parameters
        ----------
        condition_data : list[tuple[Condition, dict]]
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
                    f"Residue set mismatch between '{first_cond.label}' and "
                    f"'{cond.label}'. "
                    f"Missing in {first_cond.label}: {sorted(missing_in_first)}, "
                    f"Missing in {cond.label}: {sorted(missing_in_other)}."
                )

    # === Private helpers: pairwise comparisons ===

    def _compute_contacts_pairwise(
        self,
        summaries: list[Any],
        condition_data: list[tuple[Condition, dict[str, Any]]],
        effective_control: str | None,
    ) -> list[Any]:
        """Compute pairwise statistical comparisons for contacts.

        Unlike single-metric analyses, contacts compares TWO metrics
        (coverage and mean_contact_fraction) in each pairwise comparison.

        Parameters
        ----------
        summaries : list[ContactsConditionSummary]
            Condition summaries.
        condition_data : list[tuple[Condition, dict]]
            Raw condition data with per-replicate values.
        effective_control : str or None
            Control condition label.

        Returns
        -------
        list[ContactsPairwiseComparison]
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
                comp = self._compare_contacts_pair(
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
            for i in range(len(labels)):
                for j in range(i + 1, len(labels)):
                    label_a, label_b = labels[i], labels[j]
                    comp = self._compare_contacts_pair(
                        label_a,
                        label_to_summary[label_a],
                        label_to_data[label_a],
                        label_b,
                        label_to_summary[label_b],
                        label_to_data[label_b],
                    )
                    comparisons.append(comp)

        return comparisons

    @staticmethod
    def _resolve_effective_control(
        requested_control: str | None,
        summaries: Sequence[Any],
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
            The requested control when present, otherwise ``None`` to trigger
            all-pairs comparison.
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

    @staticmethod
    def _compare_contacts_pair(
        label_a: str,
        summary_a: Any,
        data_a: dict[str, Any],
        label_b: str,
        summary_b: Any,
        data_b: dict[str, Any],
    ) -> Any:
        """Compare two conditions for both coverage and contact fraction.

        Parameters
        ----------
        label_a, label_b : str
            Condition labels.
        summary_a, summary_b : ContactsConditionSummary
            Condition summaries.
        data_a, data_b : dict
            Raw data with per-replicate values.

        Returns
        -------
        ContactsPairwiseComparison
            Comparison with both metrics.
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

        # Coverage comparison
        coverage_a = data_a["coverage_per_replicate"]
        coverage_b = data_b["coverage_per_replicate"]

        ttest = independent_ttest(coverage_a, coverage_b)
        effect = cohens_d(coverage_a, coverage_b)
        pct = percent_change(summary_a.coverage_mean, summary_b.coverage_mean)

        # Direction: higher contact = increased
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

        # Mean contact fraction comparison
        contact_a = data_a["contact_fraction_per_replicate"]
        contact_b = data_b["contact_fraction_per_replicate"]

        ttest = independent_ttest(contact_a, contact_b)
        effect = cohens_d(contact_a, contact_b)
        pct = percent_change(
            summary_a.mean_contact_fraction,
            summary_b.mean_contact_fraction,
        )

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

    # === Private helpers: ANOVA ===

    @staticmethod
    def _compute_contacts_anova(
        condition_data: list[tuple[Condition, dict[str, Any]]],
    ) -> list[Any]:
        """Compute one-way ANOVA for both aggregate metrics.

        Parameters
        ----------
        condition_data : list[tuple[Condition, dict]]
            Condition data.

        Returns
        -------
        list[ContactsANOVASummary]
            ANOVA results for coverage and mean_contact_fraction.
        """
        from polyzymd.analyses.contacts._comparison_results import ContactsANOVASummary
        from polyzymd.analyses.shared.inferential_statistics import one_way_anova

        results = []

        # Coverage ANOVA
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

        # Mean contact fraction ANOVA
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

    @staticmethod
    def _apply_fdr_correction(
        comparisons: list[Any],
        anova_results: list[Any],
        fdr_alpha: float,
    ) -> None:
        """Apply Benjamini-Hochberg FDR correction to pairwise and ANOVA p-values.

        Treats all pairwise aggregate comparisons as one family and ANOVA p-values
        as a separate family
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

    @staticmethod
    def _apply_effect_size_threshold(
        comparisons: list[Any],
        min_effect_size: float,
    ) -> None:
        """Tag aggregate comparisons that don't meet the effect size threshold."""
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

    @staticmethod
    def _compute_top_contacted_residues(
        condition_data: list[tuple[Any, dict[str, Any]]],
        top_n: int,
    ) -> dict[str, list[tuple[int, str, float]]] | None:
        """Compute top contacted residues per condition by contact_fraction_mean.

        Parameters
        ----------
        condition_data : list
            (Condition, data_dict) pairs where data_dict["agg_result"] has residue_stats
        top_n : int
            Number of top residues to include. 0 means disabled

        Returns
        -------
        dict or None
            {condition_label: [(resid, resname, contact_fraction_mean), ...]}
            sorted descending by contact_fraction_mean. None if top_n <= 0
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
            sorted_residues = sorted(
                residue_stats,
                key=_as_contact_fraction,
                reverse=True,
            )
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

    # === Private helpers: binding preference ===

    def _load_or_compute_binding_preference(
        self,
        ctx: ComparisonContext,
        condition_data: list[tuple[Condition, dict[str, Any]]],
    ) -> Any | None:
        """Load or compute binding preference results across conditions.

        Parameters
        ----------
        ctx : ComparisonContext
            Comparison context.
        condition_data : list[tuple[Condition, dict]]
            Already-loaded contacts data for each condition.

        Returns
        -------
        BindingPreferenceComparisonSummary or None
            Cross-condition comparison summary, or None if unavailable.
        """
        from polyzymd.analyses.contacts._paths import find_contact_results_for_replicates
        from polyzymd.analyses.contacts._results import ContactResult
        from polyzymd.analyses.shared.binding_preference import compute_condition_binding_preference
        from polyzymd.analyses.shared.binding_preference_helpers import (
            binding_preference_settings_fingerprint,
            resolve_enzyme_pdb,
            try_load_cached_binding_preference,
        )

        settings = ctx.settings
        bp_settings_fp = binding_preference_settings_fingerprint(settings)
        compute_enabled = getattr(settings, "compute_binding_preference", False)
        recompute = ctx.recompute

        logger.info(f"Binding preference: compute_enabled={compute_enabled}, recompute={recompute}")

        condition_results: dict[str, Any] = {}
        surface_threshold: float | None = None

        for cond, data in condition_data:
            try:
                analysis_dir = ctx.analysis_dirs[cond.label]
                contact_results_by_replicate: dict[int, Path] = {}
                contacts_settings_fp = contacts_settings_fingerprint(settings)
                for candidate_fp in contacts_settings_fingerprint_candidates(settings):
                    contact_results_by_replicate = find_contact_results_for_replicates(
                        analysis_dir,
                        cond.replicates,
                        settings_fp=candidate_fp,
                        equilibration=ctx.equilibration,
                    )
                    if contact_results_by_replicate:
                        contacts_settings_fp = candidate_fp
                        break

                # Try cached first
                if not recompute:
                    cached = try_load_cached_binding_preference(
                        cond,
                        analysis_dir,
                        settings_fp=bp_settings_fp,
                        contact_settings_fp=contacts_settings_fp,
                        equilibration=ctx.equilibration,
                        successful_replicates=tuple(contact_results_by_replicate) or None,
                    )
                    if cached is not None:
                        condition_results[cond.label] = cached
                        if surface_threshold is None:
                            surface_threshold = cached.surface_exposure_threshold
                        logger.debug(f"  Loaded cached binding preference for {cond.label}")
                        continue

                # Compute if enabled
                if compute_enabled:
                    logger.info(f"  Computing binding preference for {cond.label}...")
                    enzyme_pdb = resolve_enzyme_pdb(
                        enzyme_pdb_setting=getattr(settings, "enzyme_pdb_for_sasa", None),
                        source_path=cond.config_path,
                        sim_config=cond.sim_config,
                    )
                    if enzyme_pdb is None or not enzyme_pdb.exists():
                        logger.warning(
                            f"Cannot compute binding preference for {cond.label}: "
                            f"enzyme PDB not found."
                        )
                        continue

                    computed = compute_condition_binding_preference(
                        cond=cond,
                        sim_config=cond.sim_config,
                        analysis_dir=analysis_dir,
                        enzyme_pdb=enzyme_pdb,
                        contact_results_by_replicate=contact_results_by_replicate,
                        load_contact_result=ContactResult.load,
                        threshold=getattr(settings, "surface_exposure_threshold", 0.2),
                        include_default_aa_groups=getattr(
                            settings, "include_default_aa_groups", True
                        ),
                        custom_protein_groups=getattr(settings, "protein_groups", None),
                        protein_partitions=getattr(settings, "protein_partitions", None),
                        polymer_type_selections=getattr(settings, "polymer_type_selections", None),
                        polymer_chain=getattr(settings, "polymer_chain", "C"),
                        settings_fp=bp_settings_fp,
                        contact_settings_fp=contacts_settings_fp,
                        equilibration=ctx.equilibration,
                    )
                    if computed is not None:
                        condition_results[cond.label] = computed
                        if surface_threshold is None:
                            surface_threshold = computed.surface_exposure_threshold
                        logger.info(f"  Computed binding preference for {cond.label}")
                        continue

            except (FileNotFoundError, ValueError, OSError, ValidationError, KeyError) as e:
                logger.warning(f"Could not load/compute binding preference for {cond.label}: {e}")
                continue

        if not condition_results:
            if compute_enabled:
                logger.warning(
                    "compute_binding_preference is enabled but no results "
                    "could be loaded or computed for any condition"
                )
            return None

        # Build comparison summary
        return self._build_binding_preference_summary(condition_results, surface_threshold)

    def _build_binding_preference_summary(
        self,
        condition_results: dict[str, Any],
        surface_threshold: float | None,
    ) -> Any:
        """Build binding preference comparison summary from per-condition results.

        Parameters
        ----------
        condition_results : dict
            Mapping of condition_label to binding preference result.
        surface_threshold : float or None
            SASA threshold used for surface filtering.

        Returns
        -------
        BindingPreferenceComparisonSummary
        """
        from polyzymd.analyses.contacts._comparison_results import (
            BindingPreferenceComparisonEntry,
            BindingPreferenceComparisonSummary,
        )
        from polyzymd.analyses.shared.binding_preference import (
            AggregatedBindingPreferenceResult,
            AggregatedPartitionBindingResult,
            BindingPreferenceResult,
            PartitionBindingResult,
        )
        from polyzymd.analyses.shared.inferential_statistics import independent_ttest

        # Collect all polymer types and AA classes
        all_polymer_types: set[str] = set()
        all_aa_classes: set[str] = set()

        for result in condition_results.values():
            bp = None
            if isinstance(result, AggregatedBindingPreferenceResult):
                bp = result.binding_preference
            elif isinstance(result, BindingPreferenceResult):
                bp = result.binding_preference
            if bp is not None:
                all_polymer_types.update(bp.polymer_types)
                all_aa_classes.update(bp.aa_class_names())

        polymer_types = sorted(all_polymer_types)
        canonical_order = ["aromatic", "polar", "nonpolar", "charged_positive", "charged_negative"]
        protein_groups = [aa for aa in canonical_order if aa in all_aa_classes]
        condition_labels = sorted(condition_results.keys())

        # Build comparison entries for each (polymer_type, aa_class) pair
        entries = []
        for poly_type in polymer_types:
            for aa_class in protein_groups:
                condition_values: dict[str, tuple[float, float]] = {}
                enrichments_for_ranking: list[tuple[str, float]] = []
                per_replicate_data: dict[str, list[float]] = {}

                for cond_label, result in condition_results.items():
                    bp = None
                    if isinstance(result, AggregatedBindingPreferenceResult):
                        bp = result.binding_preference
                    elif isinstance(result, BindingPreferenceResult):
                        bp = result.binding_preference

                    if bp is None:
                        continue

                    aa_binding = bp.aa_class_binding.get(poly_type)
                    if aa_binding is None:
                        continue

                    if isinstance(aa_binding, AggregatedPartitionBindingResult):
                        entry = None
                        for e in aa_binding.entries:
                            if e.partition_element == aa_class:
                                entry = e
                                break
                        if entry is not None:
                            mean_val = entry.mean_enrichment
                            sem_val = entry.sem_enrichment
                            if mean_val is not None:
                                condition_values[cond_label] = (
                                    mean_val,
                                    sem_val or 0.0,
                                )
                                enrichments_for_ranking.append((cond_label, mean_val))
                            if entry.per_replicate_enrichments:
                                per_replicate_data[cond_label] = entry.per_replicate_enrichments

                    elif isinstance(aa_binding, PartitionBindingResult):
                        entry = aa_binding.get_entry(aa_class)
                        if entry is not None and entry.enrichment is not None:
                            condition_values[cond_label] = (entry.enrichment, 0.0)
                            enrichments_for_ranking.append((cond_label, entry.enrichment))

                if not condition_values:
                    continue

                highest_cond = None
                lowest_cond = None
                if enrichments_for_ranking:
                    sorted_by_enrichment = sorted(
                        enrichments_for_ranking, key=lambda x: x[1], reverse=True
                    )
                    highest_cond = sorted_by_enrichment[0][0]
                    lowest_cond = sorted_by_enrichment[-1][0]

                pairwise_p_values = self._compute_bp_pairwise_pvalues(per_replicate_data)

                entries.append(
                    BindingPreferenceComparisonEntry(
                        polymer_type=poly_type,
                        protein_group=aa_class,
                        condition_values=condition_values,
                        pairwise_p_values=pairwise_p_values,
                        highest_enrichment_condition=highest_cond,
                        lowest_enrichment_condition=lowest_cond,
                    )
                )

        return BindingPreferenceComparisonSummary(
            entries=entries,
            polymer_types=polymer_types,
            protein_groups=protein_groups,
            n_conditions=len(condition_results),
            condition_labels=condition_labels,
            surface_exposure_threshold=surface_threshold,
        )

    @staticmethod
    def _compute_bp_pairwise_pvalues(
        per_replicate_data: dict[str, list[float]],
    ) -> dict[str, float]:
        """Compute pairwise t-test p-values from per-replicate enrichment data.

        Parameters
        ----------
        per_replicate_data : dict[str, list[float]]
            Mapping of condition_label to list of enrichment values.

        Returns
        -------
        dict[str, float]
            Mapping of "condA_vs_condB" to p-value.
        """
        from polyzymd.analyses.shared.inferential_statistics import independent_ttest

        if len(per_replicate_data) < 2:
            return {}

        pairwise_p_values: dict[str, float] = {}
        cond_labels = sorted(per_replicate_data.keys())

        for i, cond_a in enumerate(cond_labels):
            for cond_b in cond_labels[i + 1 :]:
                values_a = per_replicate_data[cond_a]
                values_b = per_replicate_data[cond_b]

                if len(values_a) < 2 or len(values_b) < 2:
                    continue

                try:
                    ttest_result = independent_ttest(values_a, values_b)
                    key = f"{cond_a}_vs_{cond_b}"
                    pairwise_p_values[key] = ttest_result.p_value
                except (ValueError, RuntimeError) as e:
                    logger.warning(f"T-test failed for {cond_a} vs {cond_b}: {e}")

        return pairwise_p_values
