"""Binding free energy analysis plugin — selectivity ΔG via Boltzmann inversion.

Converts cached binding preference (enrichment) data into selectivity free
energies via:

    ΔG_sel = -k_B·T · ln(contact_share / expected_share)

This is a **comparator-only** analysis: ``compute_replicate()`` and
``aggregate()`` return ``None``.  All computation is orchestrated within
``compare()`` which delegates to the existing
:class:`~polyzymd.compare.comparators.binding_free_energy.BindingFreeEnergyComparator`
via its private helpers.

Plugin contract
---------------
- ``name = "binding_free_energy"``
- ``aliases = ("bfe",)``
- ``dependencies = ("contacts",)``
- ``min_replicates = 1``
- ``compare()`` is fully overridden (temperature-aware, multi-partition)
- ``filter_conditions()`` keeps all conditions (no-polymer conditions get
  empty entries but are still valid for comparison)
- ``plot()`` delegates to ``BFEHeatmapPlotter`` and ``BFEBarPlotter``
"""

from __future__ import annotations

import logging
import math
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
    ReplicateContext,
)

if TYPE_CHECKING:
    from polyzymd.analyses._contacts_binding_preference import (
        AggregatedBindingPreferenceResult,
        AggregatedPartitionBindingEntry,
        BindingPreferenceResult,
    )
    from polyzymd.compare.results.binding_free_energy import (
        BindingFreeEnergyResult,
        FreeEnergyConditionSummary,
        FreeEnergyEntry,
        FreeEnergyPairwiseEntry,
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
    Settings: ClassVar[type] = BFESettings
    aliases: ClassVar[tuple[str, ...]] = ("bfe",)
    dependencies: ClassVar[tuple[str, ...]] = ("contacts",)
    min_replicates: ClassVar[int] = 1
    has_compute_stage: ClassVar[bool] = False
    has_aggregate_stage: ClassVar[bool] = False

    # === Compare (full override) ===

    def compare(self, ctx: ComparisonContext) -> Any | None:
        """Run binding free energy comparison across all conditions.

        Steps:
        1. Load binding preference data per condition (cached or computed)
        2. Convert enrichments to ΔG_sel via Boltzmann inversion
        3. Detect temperature groups
        4. Pairwise ΔΔG comparisons (stats suppressed for cross-temperature)
        5. Build and return result

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
        from polyzymd.compare.results.binding_free_energy import (
            BindingFreeEnergyResult,
            FreeEnergyConditionSummary,
            FreeEnergyPairwiseEntry,
        )
        from polyzymd.compare.statistics import independent_ttest

        settings: BFESettings = ctx.settings
        logger.info("Starting binding free energy comparison")
        logger.info(f"  Units: {settings.units}")
        logger.info(f"  Conditions: {len(ctx.conditions)}")

        # Step 1: Build condition summaries
        condition_summaries: list[FreeEnergyConditionSummary] = []
        for cond in ctx.conditions:
            try:
                summary = self._build_condition_summary(cond, ctx, settings)
                condition_summaries.append(summary)
            except Exception as e:
                logger.warning(f"  Skipping condition '{cond.label}': {e}")

        if not condition_summaries:
            logger.warning("No binding preference data found for any condition")
            return None

        # Step 2: Metadata
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

        # Step 4: Pairwise comparisons
        pairwise = self._compute_pairwise(condition_summaries, ctx.effective_control)

        # Step 5: Build result
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
            equilibration_time=ctx.equilibration or "",
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    # === Filter conditions ===

    def filter_conditions(self, conditions: list[Condition]) -> list[Condition]:
        """Keep all conditions — no-polymer conditions get empty entries.

        Unlike exposure/contacts, BFE includes all conditions. No-polymer
        conditions simply have no entries and serve as reference points.
        """
        return list(conditions)

    # === Plot ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate BFE comparison plots.

        Delegates to:
        - :class:`~polyzymd.compare.plotters.binding_free_energy.BFEHeatmapPlotter`
        - :class:`~polyzymd.compare.plotters.binding_free_energy.BFEBarPlotter`
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
        if plot_settings is None:
            from polyzymd.compare.config import PlotSettings

            plot_settings = PlotSettings()

        plotter_specs: list[tuple[str, str]] = [
            ("polyzymd.compare.plotters.binding_free_energy", "BFEHeatmapPlotter"),
            ("polyzymd.compare.plotters.binding_free_energy", "BFEBarPlotter"),
        ]

        for module_path, class_name in plotter_specs:
            try:
                import importlib

                mod = importlib.import_module(module_path)
                plotter_cls = getattr(mod, class_name)
                plotter = plotter_cls(settings=plot_settings)
                result = plotter.plot(data, labels, ctx.output_dir)
                if result:
                    plots.extend(result)
            except Exception as exc:
                logger.warning(f"{class_name} plot failed: {exc}")

        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format binding free energy results without legacy dispatch."""
        from polyzymd.compare.binding_free_energy_formatters import format_bfe_result

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
        from polyzymd.compare.results.binding_free_energy import (
            FreeEnergyConditionSummary,
        )
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"  Processing condition: {cond.label}")

        # Get temperature from sim config
        temperature_K = self._get_temperature(cond)

        # Compute kT
        if settings.units == "kT":
            kT = 1.0
        else:
            kT = settings.k_b() * temperature_K

        # Load binding preference
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
        from polyzymd.analyses._binding_preference_helpers import (
            CondProxy,
            compute_condition_binding_preference,
            resolve_enzyme_pdb,
            try_load_cached_binding_preference,
        )

        # Resolve contacts analysis dir (sibling of BFE analysis dir)
        bfe_analysis_dir = ctx.analysis_dirs.get(cond.label)
        contacts_analysis_dir: Path | None = None
        if bfe_analysis_dir is not None:
            contacts_analysis_dir = bfe_analysis_dir.parent / "contacts"

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
        from polyzymd.analyses._contacts_binding_preference import (
            AggregatedBindingPreferenceResult,
        )
        from polyzymd.compare.results.binding_free_energy import FreeEnergyEntry

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
        from polyzymd.compare.results.binding_free_energy import FreeEnergyEntry

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
        from polyzymd.compare.results.binding_free_energy import FreeEnergyEntry

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
        from polyzymd.compare.results.binding_free_energy import FreeEnergyEntry

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

        # Check if control has usable data
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
        """Compare two conditions for all shared (polymer_type, protein_group) pairs.

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
        from polyzymd.compare.results.binding_free_energy import FreeEnergyPairwiseEntry
        from polyzymd.compare.statistics import independent_ttest

        cross_temperature = not math.isclose(
            summary_a.temperature_K, summary_b.temperature_K, rel_tol=1e-3
        )

        pairs_a = {(e.polymer_type, e.protein_group): e for e in summary_a.entries}
        pairs_b = {(e.polymer_type, e.protein_group): e for e in summary_b.entries}
        shared_pairs = sorted(set(pairs_a.keys()) & set(pairs_b.keys()))

        pairwise: list[FreeEnergyPairwiseEntry] = []

        for polymer_type, protein_group in shared_pairs:
            entry_a = pairs_a[(polymer_type, protein_group)]
            entry_b = pairs_b[(polymer_type, protein_group)]

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
                    except Exception as exc:
                        logger.debug(f"T-test failed for ({polymer_type}, {protein_group}): {exc}")

            pairwise.append(
                FreeEnergyPairwiseEntry(
                    polymer_type=polymer_type,
                    protein_group=protein_group,
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


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Module-level polymer detection
# ---------------------------------------------------------------------------
