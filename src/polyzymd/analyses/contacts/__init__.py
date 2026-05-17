"""Contacts analysis plugin.

Computes polymer-protein contacts from MD trajectories using parallel
neighbour searching, aggregates per-residue contact fractions across
replicates, and performs cross-condition comparison with dual metrics
(coverage and mean contact fraction).

Contact computation uses the MDAnalysis-native sparse contact job seam, which
delegates to MDAnalysis ``capped_distance`` for O(N) neighbour searching.

Unlike single-scalar analyses (RMSF, catalytic_triad), contacts has **two**
primary metrics — coverage (fraction of residues contacted) and mean
contact fraction (average per-residue contact fraction).  Therefore
``compare()`` is overridden entirely and ``extract_metrics()`` is not used.

Condition filtering
-------------------
No-polymer conditions (e.g. "No Polymer" controls) are automatically
excluded via :meth:`filter_conditions`. Detection checks topology with the
active MDAnalysis polymer selection and does not treat stale contact artifacts
as polymer evidence.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence, cast

from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.contacts import _comparison as _contacts_comparison
from polyzymd.analyses.contacts import _filters as _contacts_filters
from polyzymd.analyses.contacts import _lifecycle as _contacts_lifecycle
from polyzymd.analyses.contacts._aggregator import AggregatedContactResult
from polyzymd.analyses.contacts._identity import (
    contacts_detection_fingerprint,
    normalize_polymer_types,
)
from polyzymd.analyses.contacts._mda import (
    ContactsArtifactCollector,
    build_contacts_jobs,
)
from polyzymd.analyses.contacts._plot_settings import ContactsPlotSettings
from polyzymd.analyses.mda import MDAAnalysisJob

logger = logging.getLogger("polyzymd.analyses.contacts")

if TYPE_CHECKING:
    from polyzymd.analyses.mda import MDACollectorContext, MDAReplicateJobContext


# Default cutoff matching the existing settings module
DEFAULT_CONTACT_CUTOFF = 4.5

ARCHIVED_BINDING_PREFERENCE_SETTINGS: frozenset[str] = frozenset(
    {
        "compute_binding_preference",
        "surface_exposure_threshold",
        "enzyme_pdb_for_sasa",
        "include_default_aa_groups",
        "polymer_type_selections",
        "polymer_chain",
        "enrichment_normalization",
    }
)

ARCHIVE_DIAGNOSTIC = (
    "Contacts binding-preference support has been archived and is no longer shipped as "
    "an active PolyzyMD contacts subpipeline. To recover the archived implementation, "
    "use git tag 'archive_experimental_analysis' from branch 'feature/mda-analysis-migration'."
)


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
    protein_groups : dict[str, list[int]] | None
        Custom protein groups as ``{name: [resid, ...]}``.
    protein_partitions : dict[str, list[str]] | None
        Custom partitions as ``{partition_name: [group1, ...]}``.
    fdr_alpha : float
        False discovery rate alpha for Benjamini-Hochberg correction.
    min_effect_size : float
        Minimum Cohen's d effect size to highlight.
    top_residues : int
        Number of top residues to display in console.
    """

    # --- Analysis settings ---
    polymer_selection: str = Field(
        default="chainid C", description="MDAnalysis selection for polymer atoms"
    )
    protein_selection: str = Field(
        default="chainid A", description="MDAnalysis selection for protein atoms"
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

    protein_groups: dict[str, list[int]] | None = Field(
        default=None,
        description="Custom protein groups as {name: [resid, ...]}",
    )
    protein_partitions: dict[str, list[str]] | None = Field(
        default=None,
        description="Custom partitions as {partition_name: [group1, ...]}",
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

    @model_validator(mode="before")
    @classmethod
    def reject_archived_binding_preference_settings(cls, data: Any) -> Any:
        """Reject archived contacts binding-preference settings."""

        if not isinstance(data, dict):
            return data
        archived_keys = sorted(ARCHIVED_BINDING_PREFERENCE_SETTINGS.intersection(data))
        if archived_keys:
            joined = ", ".join(archived_keys)
            raise ValueError(
                f"Archived contacts binding-preference setting(s): {joined}. {ARCHIVE_DIAGNOSTIC}"
            )
        return data

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
    MDAnalysis job seam, while keeping aggregation, comparison, plotting, and
    downstream orchestration in the plugin class.

    The ``compare()`` method is **fully overridden** because:

    - Two primary metrics: coverage and mean_contact_fraction.
    - Auto-exclusion of no-polymer conditions.
    - Residue set validation across conditions.

    Plots
    -----
    Generates stable contact plot types via private module-level functions:

    - Contact fraction / residence time profiles
    - Grouped bar charts (by AA class, by partition)
    """

    name: ClassVar[str] = "contacts"
    execution_cost_hint: ClassVar[str] = "high"
    Settings: ClassVar[type] = ContactsSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = ContactsPlotSettings
    AggregatedResultClass: ClassVar[type] = AggregatedContactResult
    ReplicateResultClass: ClassVar[type | None] = None
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 1

    # === Required methods ===

    def aggregate_settings_fingerprint(self, settings: BaseModel | None) -> str | None:
        """Return the sparse contact-detection settings fingerprint.

        Parameters
        ----------
        settings : BaseModel or None
            Active contacts settings.

        Returns
        -------
        str or None
            Contact-detection fingerprint, or ``None`` when settings are absent.
        """

        if settings is None:
            return None
        return contacts_detection_fingerprint(settings)

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[MDAAnalysisJob] | None:
        """Build the sparse contacts MDAnalysis job for one replicate."""

        return build_contacts_jobs(ctx, cast(ContactsSettings, ctx.settings))

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the contacts artifact collector."""

        del ctx
        return ContactsArtifactCollector()

    @staticmethod
    def _effective_polymer_selection(settings: "ContactsSettings") -> str:
        """Return the polymer selection constrained by polymer type filters."""

        polymer_types = normalize_polymer_types(settings.polymer_types)
        if not polymer_types:
            return settings.polymer_selection
        return f"({settings.polymer_selection}) and (resname {' '.join(polymer_types)})"

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the contacts analysis window with a dt-aware fallback."""

        return _contacts_lifecycle.get_trajectory_window(ctx, replicate, loader, universe)

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
        """Filter conditions to only those with polymer atoms."""

        return _contacts_filters.filter_conditions(self, conditions, settings)

    # === compare() — fully overridden ===

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare contacts metrics across conditions."""

        return _contacts_comparison.compare(self, ctx)

    # === plot() ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate contacts comparison plots from canonical artifacts only."""

        from polyzymd.analyses.contacts._plotters import (
            _plot_cf_by_aa_class_bars,
            _plot_cf_by_partition_bars,
            _plot_contact_fraction_profile,
            _plot_residence_time_profile,
            _plot_rt_by_aa_class_bars,
            _plot_rt_by_partition_bars,
            load_contacts_plot_data,
        )

        contacts_plot_settings = (ctx.plot_settings.model_extra or {}).get("contacts")
        if contacts_plot_settings is None:
            contacts_plot_settings = ctx.plot_settings.contacts
        plot_functions = []
        if contacts_plot_settings.generate_contact_fraction_profile:
            plot_functions.append(_plot_contact_fraction_profile)
        if ctx.settings.compute_residence_times:
            if contacts_plot_settings.generate_residence_time_profile:
                plot_functions.append(_plot_residence_time_profile)
            if contacts_plot_settings.generate_rt_by_aa_class_bars:
                plot_functions.append(_plot_rt_by_aa_class_bars)
            if contacts_plot_settings.generate_rt_by_partition_bars:
                plot_functions.append(_plot_rt_by_partition_bars)
        if contacts_plot_settings.generate_cf_by_aa_class_bars:
            plot_functions.append(_plot_cf_by_aa_class_bars)
        if contacts_plot_settings.generate_cf_by_partition_bars:
            plot_functions.append(_plot_cf_by_partition_bars)
        if not plot_functions:
            return []

        plot_data = load_contacts_plot_data(ctx)
        if not plot_data.labels:
            return []

        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        plots: list[Path] = []
        for plot_function in plot_functions:
            plots.extend(plot_function(plot_data, ctx.output_dir, ctx.plot_settings))
        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format contacts comparison results without legacy dispatch."""
        from polyzymd.analyses.contacts._formatters import format_contacts_result

        return format_contacts_result(result, format=self._normalize_output_format(output_format))

    @staticmethod
    def _normalize_output_format(output_format: str) -> str:
        return "table" if output_format == "text" else output_format

    # === Private helper wrappers: condition filtering ===

    def _condition_has_polymer(
        self, cond: Condition, settings: "ContactsAnalysis.Settings"
    ) -> bool:
        """Check whether a condition has polymer atoms."""

        return _contacts_filters.condition_has_polymer(self, cond, settings)

    # === Private helper wrappers: contacts comparison ===

    @staticmethod
    def _compute_coverage_per_replicate(result: "AggregatedContactResult") -> list[float]:
        """Compute coverage per replicate from residue statistics."""

        return _contacts_comparison.compute_coverage_per_replicate(result)

    @staticmethod
    def _compute_contact_fraction_per_replicate(result: "AggregatedContactResult") -> list[float]:
        """Compute mean contact fraction per replicate."""

        return _contacts_comparison.compute_contact_fraction_per_replicate(result)

    @staticmethod
    def _validate_residue_sets(condition_data: list[tuple[Condition, dict[str, Any]]]) -> None:
        """Validate that all conditions have identical residue sets."""

        _contacts_comparison.validate_residue_sets(condition_data)

    def _compute_contacts_pairwise(
        self,
        summaries: list[Any],
        condition_data: list[tuple[Condition, dict[str, Any]]],
        effective_control: str | None,
    ) -> list[Any]:
        """Compute pairwise statistical comparisons for contacts."""

        return _contacts_comparison.compute_contacts_pairwise(
            self, summaries, condition_data, effective_control
        )

    @staticmethod
    def _resolve_effective_control(
        requested_control: str | None,
        summaries: Sequence[Any],
    ) -> str | None:
        """Return a control label that exists in validated summaries."""

        return _contacts_comparison.resolve_effective_control(requested_control, summaries)

    @staticmethod
    def _compare_contacts_pair(
        label_a: str,
        summary_a: Any,
        data_a: dict[str, Any],
        label_b: str,
        summary_b: Any,
        data_b: dict[str, Any],
    ) -> Any:
        """Compare two conditions for coverage and contact fraction."""

        return _contacts_comparison.compare_contacts_pair(
            label_a, summary_a, data_a, label_b, summary_b, data_b
        )

    @staticmethod
    def _compute_contacts_anova(
        condition_data: list[tuple[Condition, dict[str, Any]]],
    ) -> list[Any]:
        """Compute one-way ANOVA for both aggregate metrics."""

        return _contacts_comparison.compute_contacts_anova(condition_data)

    @staticmethod
    def _apply_fdr_correction(
        comparisons: list[Any],
        anova_results: list[Any],
        fdr_alpha: float,
    ) -> None:
        """Apply Benjamini-Hochberg FDR correction to comparison p-values."""

        _contacts_comparison.apply_fdr_correction(comparisons, anova_results, fdr_alpha)

    @staticmethod
    def _apply_effect_size_threshold(comparisons: list[Any], min_effect_size: float) -> None:
        """Tag aggregate comparisons by effect-size threshold."""

        _contacts_comparison.apply_effect_size_threshold(comparisons, min_effect_size)

    @staticmethod
    def _compute_top_contacted_residues(
        condition_data: list[tuple[Any, dict[str, Any]]],
        top_n: int,
    ) -> dict[str, list[tuple[int, str, float]]] | None:
        """Compute top contacted residues per condition."""

        return _contacts_comparison.compute_top_contacted_residues(condition_data, top_n)
