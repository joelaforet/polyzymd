"""Binding free energy analysis plugin — selectivity ΔG via Boltzmann inversion.

Converts cached binding preference (enrichment) data into selectivity free
energies via:

    ΔG_sel = -k_B·T · ln(contact_share / expected_share)

This is a **compare-only** analysis: ``compute_replicate()`` and
``aggregate()`` return ``None``.  All computation is orchestrated within
``compare()`` which delegates to private helpers for enrichment → ΔG
conversion, pairwise comparison, and result assembly.

Plugin contract
---------------
- ``name = "binding_free_energy"``
- ``aliases = ("bfe",)``
- ``dependencies = ("contacts",)``
- ``min_replicates = 1``
- ``compare()`` is fully overridden (temperature-aware, multi-partition)
- ``filter_conditions()`` keeps all conditions (no-polymer conditions get
  empty entries but are still valid for comparison)
- ``plot()`` calls ``_plot_bfe_heatmap()`` and ``_plot_bfe_bars()``
"""

from __future__ import annotations

import logging
import math
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, ValidationError, field_validator

from polyzymd.analyses.base import (
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
)
from polyzymd.analyses.binding_free_energy._plot_settings import BFEPlotSettings
from polyzymd.analyses.shared.plotting import (
    annotate_cells,
    apply_axis_style,
    apply_legend,
    get_colors,
    get_output_path,
    grouped_bars,
    save_figure,
    symmetric_clim,
)

if TYPE_CHECKING:
    from polyzymd.analyses.binding_free_energy._comparison_results import (
        BindingFreeEnergyResult,
        FreeEnergyConditionSummary,
        FreeEnergyEntry,
        FreeEnergyPairwiseEntry,
    )
    from polyzymd.analyses.shared.binding_preference import (
        AggregatedBindingPreferenceResult,
        AggregatedPartitionBindingEntry,
        BindingPreferenceResult,
    )

logger = logging.getLogger("polyzymd.analyses.binding_free_energy")


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

# Boltzmann constants in various unit systems
_K_B_KCAL = 0.0019872041  # kcal/(mol·K)
_K_B_KJ = 0.0083144626  # kJ/(mol·K)

_VALID_UNITS = {"kT", "kcal/mol", "kJ/mol"}


class BFESettings(BaseModel):
    """Unified settings for binding free energy plugin.

    Merges fields from ``BindingFreeEnergyAnalysisSettings`` and
    ``BindingFreeEnergyComparisonSettings`` into a single Pydantic model.

    Parameters
    ----------
    units : str
        Energy units: ``"kT"`` (default, dimensionless), ``"kcal/mol"``, or
        ``"kJ/mol"``.
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

    units: str = "kT"
    compute_binding_preference: bool = True
    surface_exposure_threshold: float = 0.2
    enzyme_pdb_for_sasa: str | None = None
    include_default_aa_groups: bool = True
    protein_groups: dict[str, list[int]] | None = None
    protein_partitions: dict[str, list[str]] | None = None
    polymer_type_selections: dict[str, str] | None = None
    polymer_chain: str = "C"
    fdr_alpha: float = 0.05

    @field_validator("units")
    @classmethod
    def _validate_units(cls, v: str) -> str:
        if v not in _VALID_UNITS:
            raise ValueError(f"units must be one of {sorted(_VALID_UNITS)}, got '{v}'")
        return v

    @field_validator("surface_exposure_threshold")
    @classmethod
    def _validate_threshold(cls, v: float) -> float:
        if not 0.0 <= v <= 1.0:
            raise ValueError("surface_exposure_threshold must be between 0 and 1")
        return v

    def k_b(self) -> float:
        """Return Boltzmann constant in selected energy units.

        Returns 0.0 for ``"kT"`` — callers should use ``kT=1.0`` directly.
        """
        if self.units == "kT":
            return 0.0
        if self.units == "kJ/mol":
            return _K_B_KJ
        return _K_B_KCAL


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class BindingFreeEnergyAnalysis(Analysis):
    """Binding free energy (selectivity ΔG) analysis plugin.

    Converts binding preference enrichment data into physically meaningful
    free energy values via Boltzmann inversion:

        ΔG_sel = -k_B·T · ln(contact_share / expected_share)

    This is a **comparator-only** plugin: ``compute_replicate()`` and
    ``aggregate()`` return ``None``.  The full pipeline is orchestrated
    within ``compare()`` using the mixin's ``_load_or_compute()`` pattern.

    Temperature-aware: cross-temperature pairwise statistics are suppressed.

    Class Attributes
    ----------------
    name : str
        ``"binding_free_energy"``
    Settings : type
        :class:`BFESettings`
    aliases : tuple[str, ...]
        ``("bfe",)``
    dependencies : tuple[str, ...]
        ``("contacts",)``
    min_replicates : int
        ``1``
    """

    name: ClassVar[str] = "binding_free_energy"
    execution_cost_hint: ClassVar[str] = "high"
    Settings: ClassVar[type] = BFESettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = BFEPlotSettings
    aliases: ClassVar[tuple[str, ...]] = ("bfe",)
    dependencies: ClassVar[tuple[str, ...]] = ("contacts",)
    min_replicates: ClassVar[int] = 1
    has_compute_stage: ClassVar[bool] = False
    has_aggregate_stage: ClassVar[bool] = False
    settings_path_fields: ClassVar[tuple[str, ...]] = ("enzyme_pdb_for_sasa",)

    # === Compare (full override) ===

    def compare(self, ctx: ComparisonContext) -> Any | None:
        """Run binding free energy comparison across all conditions.

        Workflow:

        - Load binding preference data per condition (cached or computed).
        - Convert enrichments to ΔG_sel via Boltzmann inversion.
        - Detect temperature groups.
        - Compute pairwise ΔΔG comparisons with cross-temperature statistics
          suppressed.
        - Build and return the comparison result.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided context.

        Returns
        -------
        BindingFreeEnergyResult | None
            Comparison result, or ``None`` if no conditions have data.
        """
        from polyzymd import __version__
        from polyzymd.analyses.binding_free_energy._comparison_results import (
            BindingFreeEnergyResult,
            FreeEnergyConditionSummary,
            FreeEnergyPairwiseEntry,
        )

        settings: BFESettings = ctx.settings
        logger.info("Starting binding free energy comparison")
        logger.info(f"  Units: {settings.units}")
        logger.info(f"  Conditions: {len(ctx.conditions)}")

        # Build condition summaries from binding preference data
        condition_summaries: list[FreeEnergyConditionSummary] = []
        for cond in ctx.conditions:
            try:
                summary = self._build_condition_summary(cond, ctx, settings)
                condition_summaries.append(summary)
            except (
                FileNotFoundError,
                OSError,
                ValueError,
                TypeError,
                KeyError,
                ValidationError,
            ) as e:
                logger.warning(f"  Skipping condition '{cond.label}': {e}")

        if not condition_summaries:
            logger.warning("No binding preference data found for any condition")
            return None

        if all(not summary.entries for summary in condition_summaries):
            logger.warning("All condition summaries are empty - no binding free energy entries")
            return None

        # Collect shared metadata across condition summaries
        all_polymer_types: list[str] = sorted(
            {e.polymer_type for s in condition_summaries for e in s.entries}
        )
        all_protein_groups: list[str] = sorted(
            {e.protein_group for s in condition_summaries for e in s.entries}
        )

        # Group conditions by simulation temperature
        temp_groups: dict[float, list[str]] = {}
        for s in condition_summaries:
            temp_groups.setdefault(s.temperature_K, []).append(s.label)
        mixed_temperatures = len(temp_groups) > 1

        # Compute pairwise free-energy differences
        pairwise = self._compute_pairwise(condition_summaries, ctx.effective_control)

        # Apply BH FDR correction within comparable temperature groups
        self._apply_fdr_correction(pairwise, settings.fdr_alpha)

        # Report the formula that matches the requested output units
        if settings.units == "kT":
            formula = "ΔG_sel = -ln(contact_share / expected_share)  [units: k_bT]"
        else:
            formula = "ΔG_sel = -k_B·T · ln(contact_share / expected_share)"

        return BindingFreeEnergyResult(
            name=ctx.name,
            units=settings.units,
            formula=formula,
            mixed_temperatures=mixed_temperatures,
            temperature_groups={str(k): v for k, v in temp_groups.items()},
            conditions=condition_summaries,
            pairwise_comparisons=pairwise,
            polymer_types=all_polymer_types,
            protein_groups=all_protein_groups,
            surface_exposure_threshold=settings.surface_exposure_threshold,
            fdr_alpha=settings.fdr_alpha,
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
        """Keep all conditions — no-polymer conditions get empty entries.

        Unlike exposure/contacts, BFE includes all conditions. No-polymer
        conditions simply have no entries and serve as reference points.
        """
        return list(conditions)

    # === Plot ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate BFE comparison plots.

        Calls module-level private functions:

        - :func:`_plot_bfe_heatmap`
        - :func:`_plot_bfe_bars`
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

        data["__meta__"] = {
            "results_dir": ctx.results_dir,
            "comparison_result_path": ctx.results_dir / "result.json",
            "comparison_dir": ctx.results_dir,
        }
        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        plot_settings = ctx.plot_settings

        result = _plot_bfe_heatmap(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        result = _plot_bfe_bars(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format binding free energy results without legacy dispatch."""
        from polyzymd.analyses.binding_free_energy._formatters import format_bfe_result

        return format_bfe_result(result, format=self._normalize_output_format(output_format))

    # === extract_metrics (empty — full compare() override) ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Return empty dict — BFE uses full compare() override."""
        return {}

    @staticmethod
    def _normalize_output_format(output_format: str) -> str:
        return "table" if output_format == "text" else output_format

    # === Private helpers ===

    def _build_condition_summary(
        self,
        cond: Condition,
        ctx: ComparisonContext,
        settings: BFESettings,
    ) -> Any:
        """Build ΔG_sel entries for one condition.

        Loads or computes binding preference data, then converts enrichments
        to free energies via Boltzmann inversion.

        Parameters
        ----------
        cond : Condition
            Condition to process.
        ctx : ComparisonContext
            Comparison context.
        settings : BFESettings
            Plugin settings.

        Returns
        -------
        FreeEnergyConditionSummary
        """
        from polyzymd.analyses.binding_free_energy._comparison_results import (
            FreeEnergyConditionSummary,
        )
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"  Processing condition: {cond.label}")

        # Read temperature from sim config
        temperature_K = self._get_temperature(cond)

        # Convert temperature to thermal energy in requested units
        if settings.units == "kT":
            kT = 1.0
        else:
            kT = settings.k_b() * temperature_K

        # Load binding preference from contacts artifacts
        bp_result = self._load_binding_preference(cond, ctx, settings)

        if bp_result is None:
            return FreeEnergyConditionSummary(
                label=cond.label,
                config_path=str(cond.config_path),
                temperature_K=temperature_K,
                n_replicates=0,
                units=settings.units,
                entries=[],
                polymer_types=[],
                protein_groups=[],
            )

        entries = self._compute_dg_entries(bp_result, kT, settings.units, temperature_K)
        polymer_types = sorted({e.polymer_type for e in entries})
        protein_groups = sorted({e.protein_group for e in entries})
        n_replicates = max((e.n_replicates for e in entries), default=0)

        return FreeEnergyConditionSummary(
            label=cond.label,
            config_path=str(cond.config_path),
            temperature_K=temperature_K,
            n_replicates=n_replicates,
            units=settings.units,
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
        settings: BFESettings,
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
        settings : BFESettings
            Plugin settings.

        Returns
        -------
        AggregatedBindingPreferenceResult | BindingPreferenceResult | None
        """
        from polyzymd.analyses.contacts._paths import (
            find_contact_results_for_replicates,
            has_unproven_fingerprinted_contact_artifacts,
            infer_contact_results_settings_fingerprint,
            infer_contacts_artifact_settings_fingerprint,
        )
        from polyzymd.analyses.contacts._results import ContactResult
        from polyzymd.analyses.shared.binding_preference import compute_condition_binding_preference
        from polyzymd.analyses.shared.binding_preference_helpers import (
            binding_preference_settings_fingerprint,
            resolve_enzyme_pdb,
            try_load_cached_binding_preference,
        )

        # Resolve contacts analysis dir from the sibling dependency output
        bfe_analysis_dir = ctx.analysis_dirs.get(cond.label)
        contacts_analysis_dir: Path | None = None
        if bfe_analysis_dir is not None:
            contacts_analysis_dir = bfe_analysis_dir.parent / "contacts"

        contact_results_by_replicate: dict[int, Path] = {}
        contacts_settings_fp: str | None = None
        bp_settings_fp = binding_preference_settings_fingerprint(settings)
        if contacts_analysis_dir is not None:
            if has_unproven_fingerprinted_contact_artifacts(
                contacts_analysis_dir,
                cond.replicates,
                equilibration=ctx.equilibration,
            ):
                logger.warning(
                    "    Cannot use binding preference for '%s': contacts identity is unproven",
                    cond.label,
                )
                return None
            contact_results_by_replicate = find_contact_results_for_replicates(
                contacts_analysis_dir,
                cond.replicates,
                equilibration=ctx.equilibration,
                strict_identity=True,
            )
            inferred_artifact_fp = infer_contacts_artifact_settings_fingerprint(
                contacts_analysis_dir,
                cond.replicates,
                equilibration=ctx.equilibration,
            )
            if inferred_artifact_fp is not None:
                contacts_settings_fp = inferred_artifact_fp
            else:
                inferred_fp = infer_contact_results_settings_fingerprint(
                    contact_results_by_replicate
                )
                if inferred_fp is not None:
                    contacts_settings_fp = inferred_fp
            if contacts_settings_fp is not None:
                contact_results_by_replicate = find_contact_results_for_replicates(
                    contacts_analysis_dir,
                    cond.replicates,
                    settings_fp=contacts_settings_fp,
                    equilibration=ctx.equilibration,
                    strict_identity=True,
                )

        # Load cached binding preference before recomputing from contacts
        if not ctx.recompute and contacts_analysis_dir is not None:
            bp = try_load_cached_binding_preference(
                cond,
                contacts_analysis_dir,
                settings_fp=bp_settings_fp,
                contact_settings_fp=contacts_settings_fp,
                equilibration=ctx.equilibration,
                successful_replicates=tuple(contact_results_by_replicate) or None,
            )
            if bp is not None:
                logger.info(f"    Loaded cached binding preference for {cond.label}")
                return bp

        # Compute only when the caller opted into the fallback path
        if not settings.compute_binding_preference:
            logger.warning(
                f"    No cached binding preference for '{cond.label}' and "
                "compute_binding_preference=False"
            )
            return None

        logger.info(f"    Computing binding preference for {cond.label}...")
        sim_config = cond.sim_config

        # Resolve enzyme coordinates for SASA grouping
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
            cond=cond,
            sim_config=sim_config,
            analysis_dir=contacts_analysis_dir,
            enzyme_pdb=enzyme_pdb,
            contact_results_by_replicate=contact_results_by_replicate,
            load_contact_result=ContactResult.load,
            threshold=settings.surface_exposure_threshold,
            include_default_aa_groups=settings.include_default_aa_groups,
            custom_protein_groups=settings.protein_groups,
            protein_partitions=settings.protein_partitions,
            polymer_type_selections=settings.polymer_type_selections,
            polymer_chain=settings.polymer_chain,
            settings_fp=bp_settings_fp,
            contact_settings_fp=contacts_settings_fp,
            equilibration=ctx.equilibration,
        )

        if bp is not None:
            logger.info(f"    Computed binding preference for {cond.label}")
        else:
            logger.warning(f"    Failed to compute binding preference for '{cond.label}'")

        return bp

    def _compute_dg_entries(
        self,
        bp_result: Any,
        kT: float,
        units: str,
        temperature_K: float,
    ) -> list[Any]:
        """Compute ΔG_sel entries from a binding preference result.

        Handles both ``AggregatedBindingPreferenceResult`` (multi-replicate) and
        ``BindingPreferenceResult`` (single-replicate).

        Parameters
        ----------
        bp_result : AggregatedBindingPreferenceResult | BindingPreferenceResult
            Binding preference result.
        kT : float
            Thermal energy in chosen units.
        units : str
            Energy units label.
        temperature_K : float
            Simulation temperature.

        Returns
        -------
        list[FreeEnergyEntry]
        """
        from polyzymd.analyses.binding_free_energy._comparison_results import FreeEnergyEntry
        from polyzymd.analyses.shared.binding_preference import (
            AggregatedBindingPreferenceResult,
        )

        entries: list[FreeEnergyEntry] = []

        if isinstance(bp_result, AggregatedBindingPreferenceResult):
            entries = self._entries_from_aggregated(bp_result, kT, units, temperature_K)
        else:
            entries = self._entries_from_single(bp_result, kT, units, temperature_K)

        return entries

    def _entries_from_aggregated(
        self,
        result: Any,
        kT: float,
        units: str,
        temperature_K: float,
    ) -> list[Any]:
        """Extract ΔG_sel entries from an AggregatedBindingPreferenceResult."""
        from polyzymd.analyses.binding_free_energy._comparison_results import FreeEnergyEntry

        entries: list[FreeEnergyEntry] = []
        bp = result.binding_preference

        # AA-class partition
        for polymer_type, partition_result in bp.aa_class_binding.items():
            for entry in partition_result.entries:
                fe = self._entry_from_agg_bp_entry(
                    entry, polymer_type, "aa_class", kT, units, temperature_K
                )
                if fe is not None:
                    entries.append(fe)

        # User-defined partitions
        for partition_name, partition_dict in bp.user_defined_partitions.items():
            for polymer_type, partition_result in partition_dict.items():
                for entry in partition_result.entries:
                    fe = self._entry_from_agg_bp_entry(
                        entry, polymer_type, partition_name, kT, units, temperature_K
                    )
                    if fe is not None:
                        entries.append(fe)

        return entries

    def _entry_from_agg_bp_entry(
        self,
        entry: Any,
        polymer_type: str,
        partition_name: str,
        kT: float,
        units: str,
        temperature_K: float,
    ) -> Any | None:
        """Convert one AggregatedPartitionBindingEntry to a FreeEnergyEntry.

        Returns None if data is insufficient (zero shares).
        """
        from polyzymd.analyses.binding_free_energy._comparison_results import FreeEnergyEntry

        cs = entry.mean_contact_share
        es = entry.expected_share
        sem_cs = entry.sem_contact_share if entry.sem_contact_share is not None else 0.0

        if cs <= 0.0 or es <= 0.0:
            return None

        enrichment_ratio = cs / es
        delta_G = -kT * math.log(enrichment_ratio)

        # Delta-method uncertainty: σ(ΔG) = kT * σ_cs / cs
        delta_G_unc: float | None = None
        if sem_cs > 0:
            sigma_es = 0.0  # expected_share from single PDB — no variance
            delta_G_unc = kT * math.sqrt((sem_cs / cs) ** 2 + (sigma_es / max(es, 1e-12)) ** 2)

        # Per-replicate ΔG from per_replicate_enrichments
        dg_per_rep: list[float] = []
        for enrich_rep in entry.per_replicate_enrichments:
            ratio_rep = enrich_rep + 1.0
            if ratio_rep > 0:
                dg_per_rep.append(-kT * math.log(ratio_rep))
            else:
                dg_per_rep.append(float("nan"))

        return FreeEnergyEntry(
            polymer_type=polymer_type,
            protein_group=entry.partition_element,
            partition_name=partition_name,
            contact_share=cs,
            expected_share=es,
            enrichment_ratio=enrichment_ratio,
            delta_G=delta_G,
            delta_G_uncertainty=delta_G_unc,
            delta_G_per_replicate=dg_per_rep,
            units=units,
            temperature_K=temperature_K,
            n_replicates=len(dg_per_rep),
            n_exposed_in_group=getattr(entry, "n_exposed_in_group", 0),
        )

    def _entries_from_single(
        self,
        result: Any,
        kT: float,
        units: str,
        temperature_K: float,
    ) -> list[Any]:
        """Extract ΔG_sel entries from a single-replicate BindingPreferenceResult."""
        from polyzymd.analyses.binding_free_energy._comparison_results import FreeEnergyEntry

        entries: list[FreeEnergyEntry] = []
        bp = result.binding_preference

        # AA-class partition
        for polymer_type, partition_result in bp.aa_class_binding.items():
            for entry in partition_result.entries:
                cs = entry.contact_share
                es = entry.expected_share
                if cs <= 0.0 or es <= 0.0:
                    continue
                enrichment_ratio = cs / es
                delta_G = -kT * math.log(enrichment_ratio)

                entries.append(
                    FreeEnergyEntry(
                        polymer_type=polymer_type,
                        protein_group=entry.partition_element,
                        partition_name="aa_class",
                        contact_share=cs,
                        expected_share=es,
                        enrichment_ratio=enrichment_ratio,
                        delta_G=delta_G,
                        delta_G_uncertainty=None,
                        delta_G_per_replicate=[delta_G],
                        units=units,
                        temperature_K=temperature_K,
                        n_replicates=1,
                        n_exposed_in_group=getattr(entry, "n_exposed_in_group", 0),
                    )
                )

        # User-defined partitions
        for partition_name, partition_dict in bp.user_defined_partitions.items():
            for polymer_type, partition_result in partition_dict.items():
                for entry in partition_result.entries:
                    cs = entry.contact_share
                    es = entry.expected_share
                    if cs <= 0.0 or es <= 0.0:
                        continue
                    enrichment_ratio = cs / es
                    delta_G = -kT * math.log(enrichment_ratio)

                    entries.append(
                        FreeEnergyEntry(
                            polymer_type=polymer_type,
                            protein_group=entry.partition_element,
                            partition_name=partition_name,
                            contact_share=cs,
                            expected_share=es,
                            enrichment_ratio=enrichment_ratio,
                            delta_G=delta_G,
                            delta_G_uncertainty=None,
                            delta_G_per_replicate=[delta_G],
                            units=units,
                            temperature_K=temperature_K,
                            n_replicates=1,
                            n_exposed_in_group=getattr(entry, "n_exposed_in_group", 0),
                        )
                    )

        return entries

    # === Pairwise comparison ===

    def _compute_pairwise(
        self,
        summaries: list[Any],
        effective_control: str | None,
    ) -> list[Any]:
        """Compute pairwise ΔΔG comparisons respecting temperature grouping.

        Parameters
        ----------
        summaries : list[FreeEnergyConditionSummary]
            Condition summaries.
        effective_control : str | None
            Control label (or None).

        Returns
        -------
        list[FreeEnergyPairwiseEntry]
        """
        label_to_summary = {s.label: s for s in summaries}
        labels = [s.label for s in summaries]
        comparisons: list[Any] = []

        # Control comparisons require usable reference data
        control_has_data = (
            effective_control is not None
            and effective_control in label_to_summary
            and any(e.delta_G is not None for e in label_to_summary[effective_control].entries)
        )

        if control_has_data:
            summary_a = label_to_summary[effective_control]
            for label_b in labels:
                if label_b == effective_control:
                    continue
                summary_b = label_to_summary[label_b]
                if not any(e.delta_G is not None for e in summary_b.entries):
                    continue
                comparisons.extend(self._compare_condition_pair(summary_a, summary_b))
        else:
            # All-pairs among conditions with valid data
            if effective_control is not None and effective_control in label_to_summary:
                logger.info(
                    f"Control '{effective_control}' has no ΔG_sel data. "
                    "Falling back to all-pairs comparison."
                )
            valid_labels = [
                lbl
                for lbl in labels
                if any(e.delta_G is not None for e in label_to_summary[lbl].entries)
            ]
            for i in range(len(valid_labels)):
                for j in range(i + 1, len(valid_labels)):
                    sa = label_to_summary[valid_labels[i]]
                    sb = label_to_summary[valid_labels[j]]
                    comparisons.extend(self._compare_condition_pair(sa, sb))

        return comparisons

    def _compare_condition_pair(
        self,
        summary_a: Any,
        summary_b: Any,
    ) -> list[Any]:
        """Compare two conditions for all shared (polymer, partition, group) entries.

        Parameters
        ----------
        summary_a : FreeEnergyConditionSummary
            First condition (typically control).
        summary_b : FreeEnergyConditionSummary
            Second condition.

        Returns
        -------
        list[FreeEnergyPairwiseEntry]
        """
        from polyzymd.analyses.binding_free_energy._comparison_results import (
            FreeEnergyPairwiseEntry,
        )
        from polyzymd.analyses.shared.inferential_statistics import independent_ttest

        cross_temperature = not math.isclose(
            summary_a.temperature_K, summary_b.temperature_K, rel_tol=1e-3
        )

        pairs_a = {
            (e.polymer_type, e.partition_name, e.protein_group): e for e in summary_a.entries
        }
        pairs_b = {
            (e.polymer_type, e.partition_name, e.protein_group): e for e in summary_b.entries
        }
        shared_pairs = sorted(set(pairs_a.keys()) & set(pairs_b.keys()))

        pairwise: list[FreeEnergyPairwiseEntry] = []

        for polymer_type, partition_name, protein_group in shared_pairs:
            entry_a = pairs_a[(polymer_type, partition_name, protein_group)]
            entry_b = pairs_b[(polymer_type, partition_name, protein_group)]

            dg_a = entry_a.delta_G
            dg_b = entry_b.delta_G

            delta_delta_G: float | None = None
            if dg_a is not None and dg_b is not None:
                delta_delta_G = dg_b - dg_a

            t_stat: float | None = None
            p_val: float | None = None

            if not cross_temperature:
                reps_a = [v for v in entry_a.delta_G_per_replicate if not math.isnan(v)]
                reps_b = [v for v in entry_b.delta_G_per_replicate if not math.isnan(v)]
                if len(reps_a) >= 2 and len(reps_b) >= 2:
                    try:
                        ttest = independent_ttest(reps_a, reps_b)
                        t_stat = ttest.t_statistic
                        p_val = ttest.p_value
                    except (ValueError, TypeError, RuntimeError) as exc:
                        logger.debug(
                            "T-test failed for (%s, %s, %s): %s",
                            polymer_type,
                            partition_name,
                            protein_group,
                            exc,
                        )

            pairwise.append(
                FreeEnergyPairwiseEntry(
                    polymer_type=polymer_type,
                    protein_group=protein_group,
                    partition_name=partition_name,
                    condition_a=summary_a.label,
                    condition_b=summary_b.label,
                    temperature_a_K=summary_a.temperature_K,
                    temperature_b_K=summary_b.temperature_K,
                    cross_temperature=cross_temperature,
                    delta_G_a=dg_a,
                    delta_G_b=dg_b,
                    delta_delta_G=delta_delta_G,
                    t_statistic=t_stat,
                    p_value=p_val,
                )
            )

        return pairwise

    @staticmethod
    def _apply_fdr_correction(
        pairwise: list[Any],
        fdr_alpha: float,
    ) -> None:
        """Apply BH correction to pairwise p-values per temperature group.

        Same-temperature pairs form one hypothesis family per temperature
        Cross-temperature pairs are skipped
        """
        from polyzymd.analyses.shared.inferential_statistics import benjamini_hochberg

        same_temp = [p for p in pairwise if not p.cross_temperature and p.p_value is not None]
        if not same_temp:
            return

        temp_groups: dict[float, list[Any]] = {}
        for p in same_temp:
            temp_groups.setdefault(p.temperature_a_K, []).append(p)

        for temp, group in temp_groups.items():
            logger.debug(
                "Starting BH correction for BFE temperature group %.2f K: size=%d, alpha=%.4f",
                temp,
                len(group),
                fdr_alpha,
            )
            raw_p = [e.p_value for e in group]
            bh_results = benjamini_hochberg(raw_p, alpha=fdr_alpha)
            changed_significance = 0
            for entry, bh in zip(group, bh_results, strict=False):
                if entry.significant != bh.significant:
                    changed_significance += 1
                entry.p_value_adjusted = bh.adjusted_p_value
                entry.significant = bh.significant
            n_significant = sum(1 for entry in group if entry.significant)
            logger.info(
                "Applied BH correction to %d BFE tests at %.2f K (α=%.3f): "
                "%d remain significant, %d changed significance",
                len(group),
                temp,
                fdr_alpha,
                n_significant,
                changed_significance,
            )


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------


def _unit_label_mpl(units: str) -> str:
    """Return matplotlib-compatible unit label with subscript for kT."""
    if units == "kT":
        return r"$k_\mathrm{b}T$"
    return units


def _find_bfe_result(data: dict[str, Any], labels: Sequence[str]) -> Any | None:
    """Find and load BindingFreeEnergyResult from the comparison cache."""
    from polyzymd.analyses.binding_free_energy._comparison_results import BindingFreeEnergyResult
    from polyzymd.analyses.shared.result_io import find_comparison_result

    return find_comparison_result(
        data,
        labels,
        glob_patterns=[
            "binding_free_energy_comparison_*.json",
            "bfe_comparison_*.json",
        ],
        loader=BindingFreeEnergyResult.load,
        analysis_type="binding_free_energy",
        log=logger,
    )


def _sorted_groups(groups: list[str]) -> list[str]:
    """Sort AA groups in canonical order, with non-canonical groups appended."""
    from polyzymd.analyses.shared.aa_classification import CANONICAL_AA_CLASS_ORDER

    ordered = [g for g in CANONICAL_AA_CLASS_ORDER if g in groups]
    for g in sorted(groups):
        if g not in ordered:
            ordered.append(g)
    return ordered


def _get_partitions(result: Any) -> dict[str, list[str]]:
    """Build a mapping of partition_name -> sorted list of protein groups."""
    partition_groups: dict[str, set[str]] = {}
    for cond in result.conditions:
        for entry in cond.entries:
            partition_groups.setdefault(entry.partition_name, set()).add(entry.protein_group)
    ordered_partitions: dict[str, list[str]] = {}
    partition_names = sorted(partition_groups.keys())
    if "aa_class" in partition_names:
        partition_names.remove("aa_class")
        partition_names.insert(0, "aa_class")
    for pname in partition_names:
        ordered_partitions[pname] = _sorted_groups(list(partition_groups[pname]))
    return ordered_partitions


def _partition_display_name(partition_name: str) -> str:
    """Convert a partition name to a human-readable display string."""
    return partition_name.replace("_", " ").title()


def _build_stem(
    prefix: str,
    partition_name: str,
    poly_type: str,
    multi_partition: bool,
    multi_poly: bool,
) -> str:
    """Build output filename stem from partition and polymer type."""
    parts = [prefix]
    if multi_partition:
        parts.append(partition_name.lower())
    if multi_poly:
        parts.append(poly_type.lower())
    return "_".join(parts)


def _plot_bfe_heatmap(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate ΔG_sel heatmap comparing binding free energy across conditions.

    Creates one figure per (partition, polymer_type) combination.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    result = _find_bfe_result(data, labels)
    if result is None:
        return []

    bfe_settings = plot_settings.binding_free_energy
    if not bfe_settings.generate_heatmap:
        return []

    t = plot_settings.theme
    units = result.units

    cond_labels = [c.label for c in result.conditions]
    display_labels = [lbl for lbl in labels if lbl in cond_labels]
    if not display_labels:
        display_labels = cond_labels

    polymer_types = result.polymer_types
    partitions = _get_partitions(result)

    if not polymer_types or not partitions:
        logger.warning("BFE result has no polymer types or protein groups - skipping heatmap")
        return []

    n_conds = len(display_labels)
    n_poly = len(polymer_types)
    n_partitions = len(partitions)
    multi_partition = n_partitions > 1

    temp_str = ""
    if result.conditions:
        temps = {c.temperature_K for c in result.conditions}
        if len(temps) == 1:
            temp_str = f" at {next(iter(temps)):.0f} K"

    output_paths: list[Path] = []

    for partition_name, protein_groups in partitions.items():
        n_groups = len(protein_groups)

        partition_vals: list[float] = []
        for cond_summary in result.conditions:
            for entry in cond_summary.entries:
                if entry.partition_name == partition_name and entry.delta_G is not None:
                    partition_vals.append(entry.delta_G)

        if not partition_vals:
            logger.debug(f"No ΔG_sel values for partition '{partition_name}' - skipping")
            continue

        vmin, vmax = symmetric_clim(partition_vals, pad=0.05)
        max_abs = vmax - 0.05

        for poly_type in polymer_types:
            if bfe_settings.figsize_heatmap is not None:
                figsize = bfe_settings.figsize_heatmap
            else:
                figsize = (
                    max(6, 1.5 * n_conds + 1.5),
                    max(4, 0.9 * n_groups + 1.5),
                )

            fig, ax = plt.subplots(figsize=figsize, dpi=plot_settings.dpi)

            matrix = np.full((n_groups, n_conds), np.nan)
            sem_matrix = np.full((n_groups, n_conds), np.nan)

            for col_idx, cond_label in enumerate(display_labels):
                try:
                    cond_summary = result.get_condition(cond_label)
                except KeyError:
                    continue
                for row_idx, group in enumerate(protein_groups):
                    entry = cond_summary.get_entry(poly_type, group, partition_name=partition_name)
                    if entry is not None and entry.delta_G is not None:
                        matrix[row_idx, col_idx] = entry.delta_G
                        if entry.delta_G_uncertainty is not None:
                            sem_matrix[row_idx, col_idx] = entry.delta_G_uncertainty

            valid = matrix[~np.isnan(matrix)]
            if len(valid) == 0:
                logger.debug(
                    f"No ΔG_sel data for partition '{partition_name}', "
                    f"polymer '{poly_type}' - skipping"
                )
                plt.close(fig)
                continue

            im = ax.imshow(
                matrix,
                cmap=bfe_settings.colormap,
                vmin=vmin,
                vmax=vmax,
                aspect="auto",
            )

            if bfe_settings.annotate_heatmap:
                annotate_cells(
                    ax,
                    matrix,
                    plot_settings,
                    fontsize=t.small_fontsize,
                    threshold=0.35 * max_abs,
                    sem_matrix=sem_matrix,
                    linespacing=1.2,
                )

            ax.set_xticks(range(n_conds))
            ax.set_xticklabels(display_labels, rotation=35, ha="right")
            ax.set_yticks(range(n_groups))
            ax.set_yticklabels(protein_groups)

            if multi_partition:
                ylabel = f"Protein Group ({_partition_display_name(partition_name)})"
            else:
                ylabel = "Amino Acid Group"

            poly_label = poly_type if n_poly > 1 else ""
            if multi_partition:
                part_label = _partition_display_name(partition_name)
                title_parts = [r"$\Delta G_{\mathrm{sel}}$", part_label]
                if poly_label:
                    title_parts.append(poly_label)
                if temp_str:
                    title_parts.append(temp_str.strip())
                title = " — ".join(title_parts[:2])
                if poly_label:
                    title += f" ({poly_label})"
                if temp_str:
                    title += temp_str
            else:
                parts = [r"$\Delta G_{\mathrm{sel}}$"]
                if poly_label:
                    parts.append(poly_label)
                if temp_str:
                    parts.append(temp_str.strip())
                title = (
                    " ".join(parts)
                    if len(parts) > 1
                    else r"$\Delta G_{\mathrm{sel}}$ Binding Selectivity"
                )
            apply_axis_style(ax, plot_settings, title=title, xlabel="Condition", ylabel=ylabel)

            cbar = fig.colorbar(im, ax=ax, shrink=0.85)
            unit_lbl = _unit_label_mpl(units)
            cbar.set_label(
                r"$\Delta G_{\mathrm{sel}}$" + f" ({unit_lbl})",
                rotation=270,
                labelpad=14,
                fontsize=t.legend_fontsize,
            )
            cbar.ax.axhline(
                y=0.0,
                color=t.reference_line_color,
                linewidth=t.reference_line_width,
                linestyle=t.reference_line_style,
            )

            plt.tight_layout()

            stem = _build_stem(
                "bfe_heatmap",
                partition_name,
                poly_type,
                multi_partition,
                n_poly > 1,
            )
            output_path = get_output_path(output_dir, stem, plot_settings)
            output_paths.append(
                save_figure(
                    fig,
                    output_path,
                    plot_settings,
                    experimental_features=("binding_free_energy",),
                )
            )

    return output_paths


def _plot_bfe_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate ΔG_sel grouped bar charts comparing binding free energy across conditions.

    Creates one figure per (partition, polymer_type) combination.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    result = _find_bfe_result(data, labels)
    if result is None:
        return []

    bfe_settings = plot_settings.binding_free_energy
    if not bfe_settings.generate_bars:
        return []

    t = plot_settings.theme
    units = result.units

    cond_labels = [c.label for c in result.conditions]
    display_labels = [lbl for lbl in labels if lbl in cond_labels]
    if not display_labels:
        display_labels = cond_labels

    valid_labels = [
        lbl
        for lbl in display_labels
        if any(e.delta_G is not None for e in result.get_condition(lbl).entries)
        if lbl in cond_labels
    ]
    if not valid_labels:
        logger.info("No conditions with ΔG_sel values - skipping bar charts")
        return []

    polymer_types = result.polymer_types
    partitions = _get_partitions(result)

    if not polymer_types or not partitions:
        return []

    n_conds = len(valid_labels)
    colors = get_colors(n_conds, plot_settings)
    n_poly = len(polymer_types)
    n_partitions = len(partitions)
    multi_partition = n_partitions > 1

    temp_str = ""
    if result.conditions:
        temps = {c.temperature_K for c in result.conditions}
        if len(temps) == 1:
            temp_str = f" ({next(iter(temps)):.0f} K)"

    if units == "kT":
        kt: float | None = 1.0
    else:
        temps_list = [c.temperature_K for c in result.conditions]
        kt = None
        if temps_list:
            t_med = float(np.median(temps_list))
            tmp_settings = BFESettings(units=units)
            kt = tmp_settings.k_b() * t_med

    output_paths: list[Path] = []

    for partition_name, protein_groups in partitions.items():
        n_groups = len(protein_groups)

        for poly_type in polymer_types:
            figsize = bfe_settings.figsize_bars
            fig, ax = plt.subplots(figsize=figsize, dpi=plot_settings.dpi)

            x = np.arange(n_groups)

            series: list[tuple[str, list[float], list[float]]] = []
            for cond_label in valid_labels:
                cond_summary = result.get_condition(cond_label)
                means: list[float] = []
                sems: list[float] = []

                for group in protein_groups:
                    entry = cond_summary.get_entry(poly_type, group, partition_name=partition_name)
                    if entry is not None and entry.delta_G is not None:
                        means.append(entry.delta_G)
                        per_rep = entry.delta_G_per_replicate
                        if len(per_rep) >= 2:
                            sem = float(np.std(per_rep, ddof=1) / np.sqrt(len(per_rep)))
                        elif entry.delta_G_uncertainty is not None:
                            sem = entry.delta_G_uncertainty
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
                show_error=bfe_settings.show_error_bars,
                reference_label=r"$\Delta G_{\mathrm{sel}}$ = 0 (neutral)",
                bar_edgecolor="none",
            )

            poly_label = f": {poly_type}" if n_poly > 1 else ""
            if multi_partition:
                part_label = _partition_display_name(partition_name)
                title = r"$\Delta G_{\mathrm{sel}}$" + f" — {part_label}{poly_label}{temp_str}"
            else:
                title = r"$\Delta G_{\mathrm{sel}}$" + f"{poly_label}{temp_str}"

            if multi_partition:
                xlabel = f"Protein Group ({_partition_display_name(partition_name)})"
            else:
                xlabel = "Amino Acid Group"
            unit_lbl = _unit_label_mpl(units)
            ylabel = r"$\Delta G_{\mathrm{sel}}$" + f" ({unit_lbl})"

            apply_axis_style(ax, plot_settings, title=title, xlabel=xlabel, ylabel=ylabel)
            ax.set_xticks(x)
            ax.set_xticklabels(protein_groups, rotation=35, ha="right")
            apply_legend(
                ax,
                plot_settings,
                fontsize=t.small_fontsize,
                framealpha=0.7,
            )

            if kt is not None:
                ax.axhline(y=kt, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
                ax.axhline(y=-kt, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
                kt_label = r"$k_\mathrm{b}T$"
                ax.text(
                    n_groups - 0.5,
                    kt,
                    f"+{kt_label}",
                    color="gray",
                    fontsize=t.tiny_fontsize,
                    va="bottom",
                    ha="right",
                )
                ax.text(
                    n_groups - 0.5,
                    -kt,
                    f"\u2212{kt_label}",
                    color="gray",
                    fontsize=t.tiny_fontsize,
                    va="top",
                    ha="right",
                )

            plt.tight_layout()

            stem = _build_stem(
                "bfe_bars",
                partition_name,
                poly_type,
                multi_partition,
                n_poly > 1,
            )
            output_path = get_output_path(output_dir, stem, plot_settings)
            output_paths.append(
                save_figure(
                    fig,
                    output_path,
                    plot_settings,
                    experimental_features=("binding_free_energy",),
                )
            )

    return output_paths
