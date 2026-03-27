"""Contacts analysis plugin.

Computes polymer-protein contacts from MD trajectories using parallel
neighbour searching, aggregates per-residue contact fractions across
replicates, and performs cross-condition comparison with dual metrics
(coverage and mean contact fraction).

Contact computation uses :class:`ParallelContactAnalyzer` (inlined from the
former ``analysis.contacts.calculator_parallel`` module) which delegates to
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
excluded via :meth:`filter_conditions`.  Detection checks for:

1. Cached per-replicate contact results (if found, condition had polymer).
2. Topology inspection (MDAnalysis atom selection) as fallback.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup, Residue
    from MDAnalysis.core.universe import Universe
    from numpy.typing import NDArray

    from polyzymd.analyses._contacts_aggregator import (
        AggregatedContactResult,
    )
    from polyzymd.analyses._contacts_results import ContactResult
    from polyzymd.compare.results.contacts import ContactsComparisonResult

logger = logging.getLogger("polyzymd.analyses.contacts")


# Default cutoff matching the existing settings module
DEFAULT_CONTACT_CUTOFF = 4.5


# ---------------------------------------------------------------------------
# Inlined computation classes (formerly in analysis.contacts.calculator_parallel
# and analysis.contacts._utils)
# ---------------------------------------------------------------------------


def identify_polymer_chains(
    query_residues: "AtomGroup",
) -> list[list[tuple[int, str, "Residue"]]]:
    """Identify polymer chains and their segments from an atom group.

    Groups residues by connected fragments (molecular graph components)
    and returns structured chain information.

    Parameters
    ----------
    query_residues : AtomGroup
        MDAnalysis AtomGroup containing the polymer residues to classify.

    Returns
    -------
    list[list[tuple[int, str, Residue]]]
        List of chains, where each chain is a list of
        (chain_idx, resname, residue) tuples.
    """
    all_atoms = query_residues.atoms
    fragments = all_atoms.fragments if all_atoms.fragments else [all_atoms]

    chains = []
    for chain_idx, frag in enumerate(fragments):
        chain_residues = []
        for res in frag.residues:
            chain_residues.append((chain_idx, res.resname, res))
        if chain_residues:
            chains.append(chain_residues)

    return chains


def _get_contact_analysis_base_cls() -> type:
    """Return a ``_ContactAnalysisBase`` class that extends AnalysisBase.

    MDAnalysis is a heavy dependency imported lazily.  Because
    ``AnalysisBase`` is a C-extension-style class whose ``__bases__`` cannot
    be reassigned post-hoc, we define the subclass *inside* this factory on
    the first call and cache it for subsequent use.
    """
    cached = getattr(_get_contact_analysis_base_cls, "_cls", None)
    if cached is not None:
        return cached

    from MDAnalysis.analysis.base import AnalysisBase

    class _ContactAnalysisBase(AnalysisBase):  # type: ignore[misc]
        """MDAnalysis AnalysisBase for residue-level contact detection.

        Iterates trajectory frames using ``capped_distance`` for O(N)
        neighbour searching (serial execution, ~26 s / 1 900 frames).
        """

        def __init__(
            self,
            target_atoms: "AtomGroup",
            query_atoms: "AtomGroup",
            target_residue_indices: "NDArray[np.int64]",
            query_residue_indices: "NDArray[np.int64]",
            cutoff: float,
            n_target_residues: int,
            n_query_residues: int,
            **kwargs: Any,
        ):
            super().__init__(target_atoms.universe.trajectory, **kwargs)

            self.target_atoms = target_atoms
            self.query_atoms = query_atoms
            self.target_residue_indices = target_residue_indices
            self.query_residue_indices = query_residue_indices
            self.cutoff = cutoff
            self.n_target_residues = n_target_residues
            self.n_query_residues = n_query_residues

        def _prepare(self) -> None:
            """Initialize results array."""
            self._contact_matrix = np.zeros(
                (self.n_frames, self.n_target_residues, self.n_query_residues),
                dtype=np.uint8,
            )

        def _single_frame(self) -> None:
            """Compute contacts for a single frame using capped_distance."""
            from MDAnalysis.lib.distances import capped_distance

            target_pos = self.target_atoms.positions
            query_pos = self.query_atoms.positions
            box = self._ts.dimensions

            pairs, distances = capped_distance(
                query_pos,
                target_pos,
                max_cutoff=self.cutoff,
                box=box,
                return_distances=True,
            )

            if len(pairs) == 0:
                return

            query_atom_indices = pairs[:, 0]
            target_atom_indices = pairs[:, 1]

            query_res_indices = self.query_residue_indices[query_atom_indices]
            target_res_indices = self.target_residue_indices[target_atom_indices]

            self._contact_matrix[self._frame_index, target_res_indices, query_res_indices] = 1

        def _conclude(self) -> None:
            """Store results."""
            self.results.contact_matrix = self._contact_matrix

    _get_contact_analysis_base_cls._cls = _ContactAnalysisBase  # type: ignore[attr-defined]
    return _ContactAnalysisBase


class ParallelContactAnalyzer:
    """Optimized contact analyzer using cell-list based neighbour searching.

    Uses MDAnalysis ``capped_distance`` for O(N) neighbour searching,
    providing ~10–100× speedup over naïve distance calculations.

    Parameters
    ----------
    target_selector : MolecularSelector
        Selector for target residues (typically protein).
    query_selector : MolecularSelector
        Selector for query residues (typically polymer).
    cutoff : float
        Contact distance cutoff in Ångströms.  Default 4.0.
    grouping : ResidueGrouping | None
        Classification scheme for target residues.
        Default: ``ProteinAAClassification()``.
    """

    def __init__(
        self,
        target_selector: Any,
        query_selector: Any,
        cutoff: float = 4.0,
        grouping: Any | None = None,
    ) -> None:
        from polyzymd.analyses.shared.groupings import ProteinAAClassification

        self.target_selector = target_selector
        self.query_selector = query_selector
        self.cutoff = cutoff
        self.grouping = grouping or ProteinAAClassification()

    def run(
        self,
        universe: "Universe",
        start: int = 0,
        stop: int | None = None,
        step: int = 1,
        verbose: bool = False,
    ) -> "ContactResult":
        """Run contact analysis.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe with loaded trajectory.
        start : int
            First frame to analyse (0-indexed).
        stop : int | None
            Last frame to analyse (exclusive).  ``None`` → all frames.
        step : int
            Frame stride.
        verbose : bool
            Print progress information.

        Returns
        -------
        ContactResult
            Complete contact analysis results.
        """
        from polyzymd.analyses._contacts_results import (
            ContactResult,
            PolymerSegmentContacts,
            ResidueContactData,
            compress_contact_array,
        )

        target_result = self.target_selector.select(universe)
        query_result = self.query_selector.select(universe)

        target_atoms = target_result.atoms
        query_atoms = query_result.atoms
        target_residues = target_result.residues
        query_residues = query_result.residues

        if verbose:
            logger.info(
                f"Analyzing contacts: {len(query_residues)} query residues "
                f"({len(query_atoms)} atoms) -> {len(target_residues)} target residues "
                f"({len(target_atoms)} atoms)"
            )
            logger.info(f"Cutoff: {self.cutoff} A")

        target_atom_to_res = self._build_atom_to_residue_map(target_atoms, target_residues)
        query_atom_to_res = self._build_atom_to_residue_map(query_atoms, query_residues)

        try:
            timestep_ps = universe.trajectory.dt
        except (AttributeError, ValueError):
            timestep_ps = 1.0

        analysis = _get_contact_analysis_base_cls()(
            target_atoms=target_atoms,
            query_atoms=query_atoms,
            target_residue_indices=target_atom_to_res,
            query_residue_indices=query_atom_to_res,
            cutoff=self.cutoff,
            n_target_residues=len(target_residues),
            n_query_residues=len(query_residues),
        )
        analysis.run(start=start, stop=stop, step=step, verbose=verbose)

        contact_matrix = analysis.results.contact_matrix
        n_frames = contact_matrix.shape[0]

        if verbose:
            logger.info(f"Processing {n_frames} frames of contact data...")

        polymer_chains = identify_polymer_chains(query_residues)

        # Build residue index → chain/segment mapping
        res_to_chain_seg: dict[int, tuple[int, int, str, int]] = {}
        for chain_idx, chain in enumerate(polymer_chains):
            for seg_idx, (_, _, res) in enumerate(chain):
                for i, qr in enumerate(query_residues):
                    if qr.resid == res.resid and qr.resname == res.resname:
                        res_to_chain_seg[i] = (chain_idx, seg_idx, res.resname, res.resid)
                        break

        # Convert contact matrix → ContactResult
        residue_contacts: list[ResidueContactData] = []
        for target_idx, target_res in enumerate(target_residues):
            group = self.grouping.classify(target_res.resname)

            segment_contacts: list[PolymerSegmentContacts] = []
            target_contacts = contact_matrix[:, target_idx, :]

            for query_idx in range(target_contacts.shape[1]):
                contact_array = target_contacts[:, query_idx].astype(bool)
                if not np.any(contact_array):
                    continue

                if query_idx in res_to_chain_seg:
                    chain_idx, seg_idx, resname, resid = res_to_chain_seg[query_idx]
                else:
                    qres = query_residues[query_idx]
                    chain_idx, seg_idx = 0, query_idx
                    resname, resid = qres.resname, qres.resid

                events = compress_contact_array(contact_array)
                if events:
                    segment_contacts.append(
                        PolymerSegmentContacts(
                            polymer_resname=resname,
                            polymer_resid=resid,
                            polymer_chain_idx=chain_idx,
                            events=events,
                        )
                    )

            residue_contacts.append(
                ResidueContactData(
                    protein_resid=target_res.resid,
                    protein_resname=target_res.resname,
                    protein_group=group,
                    segment_contacts=segment_contacts,
                )
            )

        result = ContactResult(
            residue_contacts=residue_contacts,
            n_frames=n_frames,
            timestep_ps=timestep_ps * step,
            criteria_label=f"any_atom_{self.cutoff:.1f}A",
            criteria_cutoff=self.cutoff,
            start_frame=start,
            metadata={
                "target_selector": self.target_selector.label,
                "query_selector": self.query_selector.label,
                "n_polymer_chains": len(polymer_chains),
                "n_polymer_segments": sum(len(c) for c in polymer_chains),
                "optimized": True,
                "algorithm": "capped_distance",
            },
        )
        result.compute_per_residue_statistics()

        if verbose:
            logger.info(
                f"Analysis complete: {result.n_contacted_residues}/{result.n_protein_residues} "
                f"residues contacted ({result.coverage_fraction():.1%})"
            )

        return result

    # -- helpers ----------------------------------------------------------

    @staticmethod
    def _build_atom_to_residue_map(
        atoms: "AtomGroup",
        residues: "AtomGroup",
    ) -> "NDArray[np.int64]":
        """Build mapping from atom index to residue index.

        Returns an array where ``arr[atom_local_idx] = residue_idx``.
        """
        resid_to_idx = {res.resid: i for i, res in enumerate(residues)}
        atom_to_res = np.zeros(len(atoms), dtype=np.int64)
        for i, atom in enumerate(atoms):
            atom_to_res[i] = resid_to_idx.get(atom.residue.resid, 0)
        return atom_to_res


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class ContactsSettings(BaseModel):
    """Settings for contacts analysis.

    Mirrors :class:`~polyzymd.compare.settings.ContactsAnalysisSettings`
    and :class:`~polyzymd.compare.settings.ContactsComparisonSettings`
    in a single flat model.

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
            # Check for overlapping groups within this partition
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

    This plugin uses :class:`ParallelContactAnalyzer` (defined in this
    module) for per-replicate computation, aggregates across replicates,
    and performs cross-condition comparison with dual metrics (coverage and
    mean contact fraction).

    The ``compare()`` method is **fully overridden** because:

    - Two primary metrics: coverage and mean_contact_fraction.
    - Auto-exclusion of no-polymer conditions.
    - Optional binding preference sub-pipeline.
    - Residue set validation across conditions.

    Plots
    -----
    Delegates to 11 existing plotter classes across 5 modules:

    - Contact fraction / residence time profiles
    - Grouped bar charts (by AA class, by partition)
    - System coverage bar / heatmap
    - Binding preference bar / heatmap
    """

    name: ClassVar[str] = "contacts"
    Settings: ClassVar[type] = ContactsSettings
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute contacts for a single replicate.

        Uses :class:`ParallelContactAnalyzer` for optimised neighbour-search
        based contact detection.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        ContactResult
            Per-replicate contact result.
        """
        import MDAnalysis as mda

        from polyzymd.analyses.shared.loader import (
            TrajectoryLoader,
            convert_time,
            parse_time_string,
        )
        from polyzymd.analyses.shared.selectors import MDAnalysisSelector

        settings = ctx.settings

        logger.info(f"  Computing contacts for replicate {replicate}...")

        loader = TrajectoryLoader(ctx.sim_config)
        try:
            traj_info = loader.get_trajectory_info(replicate)
        except FileNotFoundError as e:
            logger.warning(f"  Skipping replicate {replicate}: trajectory data not found. {e}")
            return None

        try:
            traj_files = [str(p) for p in traj_info.trajectory_files]
            universe = mda.Universe(str(traj_info.topology_file), traj_files)

            # Convert equilibration time to start frame
            eq_value, eq_unit = parse_time_string(ctx.equilibration)
            eq_time_ps = convert_time(eq_value, eq_unit, "ps")

            try:
                timestep_ps = universe.trajectory.dt
            except (AttributeError, ValueError):
                timestep_ps = 1.0

            start_frame = int(eq_time_ps / timestep_ps) if timestep_ps > 0 else 0
            logger.info(f"    Equilibration: {ctx.equilibration} -> frame {start_frame}")

            # Create selectors
            target_selector = MDAnalysisSelector(settings.protein_selection)
            query_selector = MDAnalysisSelector(settings.polymer_selection)

            # Create analyzer and run
            analyzer = ParallelContactAnalyzer(
                target_selector=target_selector,
                query_selector=query_selector,
                cutoff=settings.cutoff,
            )
            result = analyzer.run(universe, start=start_frame)

            # Save result
            ctx.output_dir.mkdir(parents=True, exist_ok=True)
            output_file = ctx.output_dir / f"contacts_rep{replicate}.json"
            result.save(output_file)
            if ctx.result_path is not None:
                self.save_result(result, ctx.result_path)
            logger.info(f"  Saved: {output_file}")

            return result
        except Exception as e:
            logger.warning(f"  Skipping replicate {replicate}: analysis failed with error: {e}")
            return None

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate contact results across replicates for one condition.

        Uses :func:`~polyzymd.analysis.contacts.aggregator.aggregate_contact_results`
        (to be migrated in Tier 2).

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[ContactResult]
            Per-replicate contact results.

        Returns
        -------
        AggregatedContactResult
            Aggregated result with per-residue statistics.
        """
        from polyzymd.analyses._contacts_aggregator import aggregate_contact_results

        logger.info(f"  Aggregating {len(results)} replicates...")
        agg_result = aggregate_contact_results(list(results))

        # Save
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        reps = sorted(ctx.replicates)
        rep_range = f"{reps[0]}-{reps[-1]}"
        agg_file = ctx.output_dir / f"contacts_aggregated_reps{rep_range}.json"
        agg_result.save(agg_file)
        if ctx.result_path is not None:
            self.save_result(agg_result, ctx.result_path)
        logger.info(f"  Saved aggregated contacts: {agg_file}")

        return agg_result

    # === filter_conditions() — exclude no-polymer conditions ===

    def filter_conditions(self, conditions: list[Condition]) -> list[Condition]:
        """Filter conditions to only those with polymer atoms.

        Conditions without polymer atoms (e.g. "No Polymer" controls) are
        excluded since there are no polymer-protein contacts to analyse.

        Detection strategy:

        1. Check for cached per-replicate contact results — if found, the
           condition must have had polymer.
        2. Fall back to MDAnalysis topology inspection.
        3. If neither works, include the condition (let analysis fail later).

        Parameters
        ----------
        conditions : list[Condition]
            All conditions from the comparison config.

        Returns
        -------
        list[Condition]
            Conditions to include in analysis.
        """
        valid: list[Condition] = []

        for cond in conditions:
            try:
                if self._condition_has_polymer(cond):
                    valid.append(cond)
                else:
                    logger.info(
                        f"  Excluding '{cond.label}': no polymer atoms found "
                        f"with selection "
                        f"'{self._get_polymer_selection(cond)}'"
                    )
            except Exception as e:
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
            Comparison result, or ``None`` if fewer than 2 conditions.
        """
        from polyzymd import __version__
        from polyzymd.compare.results.contacts import (
            AggregateComparisonResult,
            BindingPreferenceComparisonEntry,
            BindingPreferenceComparisonSummary,
            ContactsANOVASummary,
            ContactsComparisonResult,
            ContactsConditionSummary,
            ContactsPairwiseComparison,
        )
        from polyzymd.compare.statistics import (
            cohens_d,
            independent_ttest,
            one_way_anova,
            percent_change,
        )

        settings = ctx.settings

        logger.info(f"Starting contacts comparison: {ctx.name}")
        logger.info(f"Conditions: {len(ctx.conditions)}")
        logger.info(f"Equilibration: {ctx.equilibration}")

        if len(ctx.conditions) < 2:
            logger.warning("contacts: fewer than 2 conditions — skipping comparison.")
            return None

        # Step 1: Load aggregated results and build condition data
        condition_data: list[tuple[Condition, dict[str, Any]]] = []
        for cond in ctx.conditions:
            agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
            agg_result = self._load_aggregated_result(agg_dir)
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

        if len(condition_data) < 2:
            logger.warning("contacts: fewer than 2 conditions have results — skipping.")
            return None

        # Step 2: Validate identical residue sets
        self._validate_residue_sets(condition_data)

        # Step 3: Build condition summaries
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
            )
            summaries.append(summary)

        # Step 4: Effective control
        effective_control = ctx.effective_control

        # Step 5: Pairwise comparisons (dual metrics)
        comparisons = self._compute_contacts_pairwise(summaries, condition_data, effective_control)

        # Step 6: ANOVA if 3+ conditions
        anova_results: list[ContactsANOVASummary] = []
        if len(summaries) >= 3:
            anova_results = self._compute_contacts_anova(condition_data)

        # Step 7: Rankings
        ranked_coverage = sorted(summaries, key=lambda s: s.coverage_mean, reverse=True)
        ranked_contact = sorted(summaries, key=lambda s: s.mean_contact_fraction, reverse=True)

        # Step 8: Binding preference comparison (optional)
        binding_pref_summary = self._load_or_compute_binding_preference(ctx, condition_data)

        # Step 9: Build result
        return ContactsComparisonResult(
            name=ctx.name,
            contacts_name="polymer_contacts",
            contacts_description=None,
            polymer_selection=settings.polymer_selection,
            protein_selection=settings.protein_selection,
            cutoff=settings.cutoff,
            contact_criteria="distance",
            fdr_alpha=settings.fdr_alpha,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova_results,
            ranking_by_coverage=[s.label for s in ranked_coverage],
            ranking_by_contact_fraction=[s.label for s in ranked_contact],
            excluded_conditions=[c.label for c in ctx.excluded_conditions],
            binding_preference=binding_pref_summary,
            equilibration_time=ctx.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    # === plot() ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate contacts comparison plots.

        Delegates to the existing 11 plotter classes:

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

        # Build data dict expected by the old plotter API
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

        # Add __meta__ with results_dir for comparison result lookup
        data["__meta__"] = {"results_dir": ctx.results_dir}

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve plot settings
        plot_settings = ctx.plot_settings
        if plot_settings is None:
            from polyzymd.compare.config import PlotSettings

            plot_settings = PlotSettings()

        # All 11 plotter classes
        plotter_specs: list[tuple[str, str]] = [
            ("polyzymd.compare.plotters.contacts_profiles", "ContactFractionProfilePlotter"),
            ("polyzymd.compare.plotters.contacts_profiles", "ResidenceTimeProfilePlotter"),
            (
                "polyzymd.compare.plotters.contacts_grouped_bars",
                "ContactFractionByAAClassBarPlotter",
            ),
            (
                "polyzymd.compare.plotters.contacts_grouped_bars",
                "ContactFractionByPartitionBarPlotter",
            ),
            ("polyzymd.compare.plotters.contacts_grouped_bars", "ResidenceTimeByAAClassBarPlotter"),
            (
                "polyzymd.compare.plotters.contacts_grouped_bars",
                "ResidenceTimeByPartitionBarPlotter",
            ),
            ("polyzymd.compare.plotters.contacts_grouped_bars", "UserPartitionBarPlotter"),
            ("polyzymd.compare.plotters.contacts_coverage", "SystemCoverageBarPlotter"),
            ("polyzymd.compare.plotters.contacts_coverage", "SystemCoverageHeatmapPlotter"),
            (
                "polyzymd.compare.plotters.contacts_binding_preference",
                "BindingPreferenceBarPlotter",
            ),
            (
                "polyzymd.compare.plotters.contacts_binding_preference",
                "BindingPreferenceHeatmapPlotter",
            ),
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
        """Format contacts comparison results without legacy dispatch."""
        from polyzymd.compare.contacts_formatters import format_contacts_result

        return format_contacts_result(result, format=self._normalize_output_format(output_format))

    # === Framework hooks ===

    def _deserialize_result(self, path: Path) -> Any:
        """Load an aggregated contact result from JSON.

        Parameters
        ----------
        path : Path
            Path to the JSON result file.

        Returns
        -------
        AggregatedContactResult
            Loaded result.
        """
        from polyzymd.analyses._contacts_aggregator import AggregatedContactResult

        return AggregatedContactResult.load(path)

    @staticmethod
    def _normalize_output_format(output_format: str) -> str:
        return "table" if output_format == "text" else output_format

    # === Private helpers: condition filtering ===

    def _get_polymer_selection(self, cond: Condition) -> str:
        """Get the polymer selection string, trying settings then default."""
        # Settings may not be resolved yet during filter_conditions;
        # use a safe fallback.
        return "chainID C"

    def _condition_has_polymer(self, cond: Condition) -> bool:
        """Check whether a condition has polymer atoms.

        Checks:
        1. Cached per-replicate contact results in the standard output
           directories.
        2. Topology inspection via MDAnalysis.

        Parameters
        ----------
        cond : Condition
            Condition to check.

        Returns
        -------
        bool
            True if polymer atoms detected.
        """
        # Check 1: Cached results
        for rep in cond.replicates:
            # Check multiple standard locations
            run_dir = cond.sim_config.get_working_directory(rep)
            projects_dir = cond.sim_config.output.projects_directory

            result_paths = [
                projects_dir / "analysis" / "contacts" / f"contacts_rep{rep}.json",
            ]
            for rp in result_paths:
                if rp.exists():
                    logger.debug(f"  {cond.label} rep {rep}: found cached result at {rp}")
                    return True

        # Check 2: Topology inspection
        for rep in cond.replicates:
            run_dir = cond.sim_config.get_working_directory(rep)
            topology_path = run_dir / "solvated_system.pdb"

            if not topology_path.exists():
                continue

            try:
                import MDAnalysis as mda

                universe = mda.Universe(str(topology_path))
                polymer_atoms = universe.select_atoms("chainID C")
                if len(polymer_atoms) > 0:
                    logger.debug(f"  {cond.label} rep {rep}: {len(polymer_atoms)} polymer atoms")
                    return True
                else:
                    logger.debug(f"  {cond.label} rep {rep}: 0 polymer atoms")
            except Exception as e:
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
        from polyzymd.compare.results.contacts import (
            AggregateComparisonResult,
            ContactsPairwiseComparison,
        )
        from polyzymd.compare.statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )

        aggregate_comps = []

        # Coverage comparison
        coverage_a = data_a["coverage_per_replicate"]
        coverage_b = data_b["coverage_per_replicate"]

        ttest = independent_ttest(coverage_a, coverage_b)
        effect = cohens_d(coverage_a, coverage_b, rmsf_mode=False)
        pct = percent_change(summary_a.coverage_mean, summary_b.coverage_mean)

        # Direction: higher contact = increased
        if abs(pct) < 1:
            direction = "unchanged"
        elif pct < 0:
            direction = "decreased"
        else:
            direction = "increased"

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
        effect = cohens_d(contact_a, contact_b, rmsf_mode=False)
        pct = percent_change(
            summary_a.mean_contact_fraction,
            summary_b.mean_contact_fraction,
        )

        if abs(pct) < 1:
            direction = "unchanged"
        elif pct < 0:
            direction = "decreased"
        else:
            direction = "increased"

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
        from polyzymd.compare.results.contacts import ContactsANOVASummary
        from polyzymd.compare.statistics import one_way_anova

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
        from polyzymd.analyses._binding_preference_helpers import (
            compute_condition_binding_preference,
            resolve_enzyme_pdb,
            try_load_cached_binding_preference,
        )
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
            BindingPreferenceResult,
        )

        settings = ctx.settings
        compute_enabled = getattr(settings, "compute_binding_preference", False)
        recompute = ctx.recompute

        logger.info(f"Binding preference: compute_enabled={compute_enabled}, recompute={recompute}")

        condition_results: dict[str, Any] = {}
        surface_threshold: float | None = None

        for cond, data in condition_data:
            try:
                analysis_dir = ctx.analysis_dirs[cond.label]

                # Try cached first
                if not recompute:
                    # Build a minimal ConditionConfig-like object for the helper
                    cached = self._try_load_cached_bp(cond, analysis_dir)
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
                        source_path=None,  # Plugin doesn't have source_path
                        sim_config=cond.sim_config,
                    )
                    if enzyme_pdb is None or not enzyme_pdb.exists():
                        logger.warning(
                            f"Cannot compute binding preference for {cond.label}: "
                            f"enzyme PDB not found."
                        )
                        continue

                    # Build a minimal ConditionConfig-like object
                    computed = self._compute_bp_for_condition(
                        cond, analysis_dir, enzyme_pdb, settings
                    )
                    if computed is not None:
                        condition_results[cond.label] = computed
                        if surface_threshold is None:
                            surface_threshold = computed.surface_exposure_threshold
                        logger.info(f"  Computed binding preference for {cond.label}")
                        continue

            except Exception as e:
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

    def _try_load_cached_bp(
        self,
        cond: Condition,
        analysis_dir: Path,
    ) -> Any | None:
        """Try to load cached binding preference results.

        Parameters
        ----------
        cond : Condition
            Condition to load.
        analysis_dir : Path
            Analysis directory for this condition.

        Returns
        -------
        AggregatedBindingPreferenceResult | BindingPreferenceResult | None
        """
        import glob as glob_module

        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
            BindingPreferenceResult,
            aggregate_binding_preference,
        )

        # Try aggregated result first
        agg_path = analysis_dir / "binding_preference_aggregated.json"
        if agg_path.exists():
            return AggregatedBindingPreferenceResult.load(agg_path)

        # Try with rep range in name
        agg_pattern = str(analysis_dir / "binding_preference_aggregated_reps*.json")
        agg_matches = sorted(glob_module.glob(agg_pattern))
        if agg_matches:
            return AggregatedBindingPreferenceResult.load(agg_matches[-1])

        # Try single replicate result
        single_path = analysis_dir / "binding_preference.json"
        if single_path.exists():
            return BindingPreferenceResult.load(single_path)

        # Try per-replicate results and aggregate them
        rep_results = []
        for rep in cond.replicates:
            rep_path = analysis_dir / f"binding_preference_rep{rep}.json"
            if rep_path.exists():
                rep_results.append(BindingPreferenceResult.load(rep_path))

        if rep_results:
            return aggregate_binding_preference(rep_results)

        return None

    def _compute_bp_for_condition(
        self,
        cond: Condition,
        analysis_dir: Path,
        enzyme_pdb: Path,
        settings: Any,
    ) -> Any | None:
        """Compute binding preference for a single condition.

        Parameters
        ----------
        cond : Condition
            Condition to compute.
        analysis_dir : Path
            Analysis directory containing contacts results.
        enzyme_pdb : Path
            Path to enzyme PDB for SASA calculation.
        settings : ContactsSettings
            Plugin settings.

        Returns
        -------
        AggregatedBindingPreferenceResult or None
        """
        import MDAnalysis as mda

        from polyzymd.analysis.contacts.binding_preference import (
            PolymerComposition,
            aggregate_binding_preference,
            compute_binding_preference,
            extract_polymer_composition,
            resolve_protein_groups_from_surface_exposure,
        )
        from polyzymd.analyses._contacts_results import ContactResult
        from polyzymd.analyses._contacts_surface_exposure import SurfaceExposureFilter

        threshold = getattr(settings, "surface_exposure_threshold", 0.2)
        include_defaults = getattr(settings, "include_default_aa_groups", True)
        custom_groups = getattr(settings, "protein_groups", None)
        protein_partitions = getattr(settings, "protein_partitions", None)
        polymer_type_selections = getattr(settings, "polymer_type_selections", None)

        # Compute surface exposure
        try:
            exposure_filter = SurfaceExposureFilter(threshold=threshold)
            surface_exposure = exposure_filter.calculate(str(enzyme_pdb))
        except Exception as e:
            logger.warning(f"Failed to compute surface exposure for {cond.label}: {e}")
            return None

        # Resolve protein groups
        protein_groups = resolve_protein_groups_from_surface_exposure(
            surface_exposure,
            include_default_aa_groups=include_defaults,
            custom_protein_groups=custom_groups,
        )
        if not protein_groups:
            logger.warning(f"No protein groups resolved for {cond.label}")
            return None

        # Extract polymer composition
        polymer_composition = None
        first_rep = cond.replicates[0] if cond.replicates else 1
        run_dir = cond.sim_config.get_working_directory(first_rep)
        topology_path = run_dir / "solvated_system.pdb"

        if topology_path.exists():
            try:
                universe = mda.Universe(str(topology_path))
                polymer_composition = extract_polymer_composition(universe, polymer_type_selections)
            except Exception as e:
                logger.warning(f"Failed to extract polymer composition for {cond.label}: {e}")

        if polymer_composition is None:
            polymer_composition = PolymerComposition()

        # Compute per-replicate
        rep_results = []
        for rep in cond.replicates:
            contact_path = analysis_dir / f"contacts_rep{rep}.json"
            if not contact_path.exists():
                # Also check run_N subdirectory
                contact_path = analysis_dir / f"run_{rep}" / f"contacts_rep{rep}.json"
                if not contact_path.exists():
                    logger.warning(f"Contacts file not found: {contact_path}")
                    continue

            try:
                contact_result = ContactResult.load(contact_path)
                bp_result = compute_binding_preference(
                    contact_result=contact_result,
                    surface_exposure=surface_exposure,
                    protein_groups=protein_groups,
                    polymer_composition=polymer_composition,
                    protein_partitions=protein_partitions,
                )
                rep_results.append(bp_result)

                rep_bp_path = analysis_dir / f"binding_preference_rep{rep}.json"
                bp_result.save(rep_bp_path)
            except Exception as e:
                logger.warning(
                    f"Failed to compute binding preference for {cond.label} rep{rep}: {e}"
                )
                continue

        if not rep_results:
            return None

        agg_result = aggregate_binding_preference(rep_results)
        rep_range = f"{min(cond.replicates)}-{max(cond.replicates)}"
        agg_path = analysis_dir / f"binding_preference_aggregated_reps{rep_range}.json"
        agg_result.save(agg_path)
        logger.info(
            f"Computed binding preference for {cond.label}: "
            f"{len(rep_results)} replicates, {len(protein_groups)} protein groups"
        )

        return agg_result

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
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
            AggregatedPartitionBindingResult,
            BindingPreferenceResult,
            PartitionBindingResult,
        )
        from polyzymd.compare.results.contacts import (
            BindingPreferenceComparisonEntry,
            BindingPreferenceComparisonSummary,
        )
        from polyzymd.compare.statistics import independent_ttest

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
        from polyzymd.compare.statistics import independent_ttest

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
                except Exception as e:
                    logger.warning(f"T-test failed for {cond_a} vs {cond_b}: {e}")

        return pairwise_p_values
