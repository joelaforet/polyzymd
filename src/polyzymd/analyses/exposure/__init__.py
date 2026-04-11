"""Exposure dynamics analysis plugin — chaperone-like polymer activity detection.

Combines per-frame SASA (MDTraj shrake_rupley) with polymer-protein contact data
to classify residues, detect chaperone events, compute enrichment, and statistically
compare conditions.

This is a **comparator-only** analysis: the real work happens in ``compare()``.
``compute_replicate()`` and ``aggregate()`` return ``None`` because the exposure
pipeline needs both SASA and contacts simultaneously — the per-replicate computation
is orchestrated within ``compare()`` directly.

Plugin contract
---------------
- ``name = "exposure"``
- ``aliases = ()``
- ``dependencies = ("contacts",)`` — contacts must be pre-computed
- ``min_replicates = 2``
- ``compare()`` is fully overridden (custom multi-metric flow)
- ``filter_conditions()`` excludes no-polymer conditions
- ``plot()`` generates 2 figure types via private module-level functions
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    Analysis,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
)
from polyzymd.analyses.shared.plotting import (
    annotate_cells,
    apply_axis_style,
    get_colors,
    get_output_path,
    save_figure,
)

if TYPE_CHECKING:
    from polyzymd.analyses.contacts._results import ContactResult
    from polyzymd.analyses.exposure._dynamics import ExposureDynamicsResult
    from polyzymd.analyses.exposure._enrichment import ChaperoneEnrichmentResult
    from polyzymd.analyses.exposure._sasa_trajectory import SASATrajectoryResult

logger = logging.getLogger("polyzymd.analyses.exposure")

# ---------------------------------------------------------------------------
# Default values (shared with ExposureAnalysisSettings in settings.py)
# ---------------------------------------------------------------------------
DEFAULT_SURFACE_EXPOSURE_THRESHOLD = 0.2


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class ExposureSettings(BaseModel):
    """Unified settings for exposure dynamics analysis plugin.

    Merges fields from ``ExposureAnalysisSettings`` and
    ``ExposureComparisonSettings`` into a single flat Pydantic model.

    Parameters
    ----------
    protein_selection : str
        MDAnalysis selection for the protein (default ``"protein"``).
    polymer_selection : str
        MDAnalysis selection for polymer (default ``"chainID C"``).
    exposure_threshold : float
        Relative SASA threshold distinguishing exposed/buried (default 0.2).
    transient_lower : float
        Fraction below which a residue is "stably buried" (default 0.2).
    transient_upper : float
        Fraction above which a residue is "stably exposed" (default 0.8).
    min_event_length : int
        Minimum contiguous exposed frames per event (default 1).
    probe_radius_nm : float
        Probe radius for SASA calculation in nm (default 0.14).
    n_sphere_points : int
        Sphere points for shrake_rupley (default 960).
    protein_chain : str
        Chain letter for protein (default ``"A"``).
    polymer_resnames : list[str] | None
        Subset of polymer residue names to consider (``None`` = all).
    """

    protein_selection: str = "protein"
    polymer_selection: str = "chainID C"
    exposure_threshold: float = DEFAULT_SURFACE_EXPOSURE_THRESHOLD
    transient_lower: float = 0.2
    transient_upper: float = 0.8
    min_event_length: int = 1
    probe_radius_nm: float = 0.14
    n_sphere_points: int = 960
    protein_chain: str = "A"
    polymer_resnames: list[str] | None = None

    @field_validator("exposure_threshold")
    @classmethod
    def _validate_threshold(cls, v: float) -> float:
        if not 0.0 < v < 1.0:
            raise ValueError("exposure_threshold must be between 0 and 1 (exclusive)")
        return v

    @field_validator("transient_upper")
    @classmethod
    def _validate_transient_bounds(cls, v: float, info: Any) -> float:
        lower = info.data.get("transient_lower", 0.2)
        if v <= lower:
            raise ValueError("transient_upper must be > transient_lower")
        return v


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class ExposureAnalysis(Analysis):
    """Exposure dynamics / chaperone activity analysis plugin.

    Detects polymer-assisted "chaperone events" by combining:
    - Per-frame SASA (protein surface accessibility trajectory)
    - Polymer-protein contact data (from contacts analysis)
    - Residue classification (stably buried / transient / stably exposed)

    This is a **comparator-only** plugin: ``compute_replicate()`` and
    ``aggregate()`` return ``None``. All computation is orchestrated
    within ``compare()``.

    Class Attributes
    ----------------
    name : str
        ``"exposure"``
    Settings : type
        :class:`ExposureSettings`
    aliases : tuple[str, ...]
        ``()``
    dependencies : tuple[str, ...]
        ``("contacts",)`` — contacts must be pre-computed.
    min_replicates : int
        ``2``
    """

    name: ClassVar[str] = "exposure"
    execution_cost_hint: ClassVar[str] = "high"
    Settings: ClassVar[type] = ExposureSettings
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ("contacts",)
    min_replicates: ClassVar[int] = 2
    has_compute_stage: ClassVar[bool] = False
    has_aggregate_stage: ClassVar[bool] = False

    # === Compare (full override) ===

    def compare(self, ctx: ComparisonContext) -> Any | None:
        """Run exposure dynamics comparison across all conditions.

        Orchestrates the full pipeline per replicate:
        1. Load pre-computed ContactResult
        2. Compute/load cached SASATrajectoryResult
        3. Compute/load cached ExposureDynamicsResult
        4. Compute ChaperoneEnrichmentResult
        5. Build per-condition summaries
        6. Pairwise t-tests on chaperone_fraction
        7. ANOVA if 3+ conditions
        8. Dual rankings: chaperone_fraction and transient_fraction

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided context.

        Returns
        -------
        ExposureComparisonResult or None
            Comparison result, or ``None`` when no valid conditions are available.
        """
        from polyzymd import __version__
        from polyzymd.analyses.base import ANOVAResult, PairwiseResult
        from polyzymd.analyses.exposure._comparison_results import (
            ExposureComparisonResult,
            ExposureConditionSummary,
        )
        from polyzymd.analyses.shared.inferential_statistics import one_way_anova

        settings = ctx.settings
        fdr_alpha = ctx.fdr_alpha
        ttest_method = ctx.ttest_method
        posthoc_method = ctx.posthoc_method

        logger.info("Starting exposure dynamics comparison")
        logger.info(f"  Conditions: {len(ctx.conditions)}")
        logger.info(f"  Equilibration: {ctx.equilibration}")

        conditions = ctx.conditions
        if not conditions:
            logger.warning("No conditions provided — skipping exposure comparison")
            return None

        # Step 1: Load/compute per condition
        summaries: list[ExposureConditionSummary] = []
        for cond in conditions:
            try:
                summary = self._build_condition_summary(cond, ctx, settings)
                summaries.append(summary)
            except (FileNotFoundError, ValueError, OSError, IndexError) as e:
                logger.debug("Full traceback:", exc_info=True)
                logger.warning(f"  Skipping condition '{cond.label}': {e}")

        if not summaries:
            logger.warning("No conditions with valid data — skipping exposure comparison")
            return None

        # Step 2: Determine effective control
        effective_control = ctx.effective_control

        # Step 3: Pairwise comparisons on chaperone_fraction
        comparisons: list[PairwiseResult] = []
        if len(summaries) >= 2:
            if effective_control:
                control_summary = next((s for s in summaries if s.label == effective_control), None)
                if control_summary:
                    treatments = [s for s in summaries if s.label != effective_control]
                    for treatment in treatments:
                        comp = self._compare_pair(
                            control_summary,
                            treatment,
                            fdr_alpha=fdr_alpha,
                            ttest_method=ttest_method,
                            posthoc_method=posthoc_method,
                        )
                        comparisons.append(comp)
            else:
                for i, cond_a in enumerate(summaries):
                    for cond_b in summaries[i + 1 :]:
                        comp = self._compare_pair(
                            cond_a,
                            cond_b,
                            fdr_alpha=fdr_alpha,
                            ttest_method=ttest_method,
                            posthoc_method=posthoc_method,
                        )
                        comparisons.append(comp)

        # Step 4: ANOVA if 3+ conditions
        anova: ANOVAResult | None = None
        if len(summaries) >= 3:
            all_groups = [s.replicate_values for s in summaries]
            result = one_way_anova(*all_groups)
            anova = ANOVAResult(
                metric="chaperone_fraction",
                f_statistic=result.f_statistic,
                p_value=result.p_value,
                significant=result.p_value <= fdr_alpha,
            )

        # Step 5: Dual rankings
        ranked_chaperone = sorted(summaries, key=lambda s: s.mean_chaperone_fraction, reverse=True)
        ranked_transient = sorted(summaries, key=lambda s: s.mean_transient_fraction, reverse=True)

        # Build excluded list from context
        excluded_labels = [c.label for c in ctx.excluded_conditions]

        return ExposureComparisonResult(
            name="exposure_comparison",
            metric="chaperone_fraction",
            control_label=effective_control,
            fdr_alpha=fdr_alpha,
            ttest_method=ttest_method,
            posthoc_method=posthoc_method,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova,
            ranking=[s.label for s in ranked_chaperone],
            ranking_by_transient_fraction=[s.label for s in ranked_transient],
            excluded_conditions=excluded_labels,
            equilibration_time=ctx.equilibration or "0ns",
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    # === Filter conditions (exclude no-polymer) ===

    def filter_conditions(
        self,
        conditions: list[Condition],
        settings: "BaseModel | None" = None,
    ) -> list[Condition]:
        """Exclude conditions without polymer (no chaperone events possible).

        A condition is excluded if:
        1. Its ``sim_config.polymers`` is ``None`` or disabled — no polymer.
        2. No pre-computed contact JSON exists for any replicate.

        Parameters
        ----------
        conditions : list[Condition]
            All candidate conditions.

        Returns
        -------
        list[Condition]
            Conditions that have polymer and contacts data.
        """
        valid: list[Condition] = []

        for cond in conditions:
            try:
                if not self._condition_has_polymer(cond):
                    logger.info(f"  Excluding '{cond.label}': no polymer configured")
                    continue
                valid.append(cond)
            except Exception as e:
                logger.warning(f"  Error checking condition '{cond.label}': {e}")
                valid.append(cond)  # fail-open

        return valid

    # === Plot ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate exposure comparison plots."""
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
            result = _plot_chaperone_fraction(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Chaperone fraction plot failed: {exc}")

        try:
            result = _plot_enrichment_heatmap(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Enrichment heatmap plot failed: {exc}")

        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format exposure comparison results without legacy dispatch."""
        from polyzymd.analyses.exposure._formatters import format_exposure_result

        return format_exposure_result(result, format=self._normalize_output_format(output_format))

    # === extract_metrics (empty — full compare() override) ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Return empty dict — exposure uses full compare() override."""
        return {}

    @staticmethod
    def _normalize_output_format(output_format: str) -> str:
        return "table" if output_format == "text" else output_format

    # === Private helpers ===

    def _build_condition_summary(
        self,
        cond: Condition,
        ctx: ComparisonContext,
        settings: ExposureSettings,
    ) -> Any:
        """Build condition summary from per-replicate exposure dynamics results.

        For each replicate, this method:
        1. Loads cached ContactResult
        2. Computes/loads SASA trajectory
        3. Computes/loads exposure dynamics
        4. Computes chaperone enrichment

        Then aggregates across replicates.

        Parameters
        ----------
        cond : Condition
            Condition to process.
        ctx : ComparisonContext
            Comparison context with analysis dirs.
        settings : ExposureSettings
            Plugin settings.

        Returns
        -------
        ExposureConditionSummary
        """
        from polyzymd.analyses.contacts._results import ContactResult
        from polyzymd.analyses.exposure._comparison_results import ExposureConditionSummary
        from polyzymd.analyses.exposure._config import ExposureConfig
        from polyzymd.analyses.exposure._dynamics import (
            ExposureDynamicsResult,
            analyze_exposure_dynamics,
        )
        from polyzymd.analyses.exposure._enrichment import compute_chaperone_enrichment
        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import compute_trajectory_sasa
        from polyzymd.analyses.shared.loader import TrajectoryLoader
        from polyzymd.analyses.shared.statistics import compute_sem

        logger.info(f"  Processing condition: {cond.label}")

        # Resolve analysis directories
        exposure_analysis_dir = ctx.analysis_dirs.get(cond.label)
        contacts_analysis_dir: Path | None = None
        if exposure_analysis_dir is not None:
            contacts_analysis_dir = exposure_analysis_dir.parent / "contacts"

        dynamics_per_rep: list[ExposureDynamicsResult] = []
        enrichment_per_rep: list[ChaperoneEnrichmentResult] = []

        for rep in cond.replicates:
            result = self._load_or_compute_replicate(
                cond=cond,
                replicate=rep,
                settings=settings,
                exposure_analysis_dir=exposure_analysis_dir,
                contacts_analysis_dir=contacts_analysis_dir,
                recompute=ctx.recompute,
            )
            if result is not None:
                dynamics, enrichment = result
                dynamics_per_rep.append(dynamics)
                enrichment_per_rep.append(enrichment)

        if not dynamics_per_rep:
            raise ValueError(f"No successful replicates for condition: {cond.label}")

        n_replicates = len(dynamics_per_rep)

        # Per-replicate primary metric: chaperone_fraction over transient residues
        chap_fractions: list[float] = []
        transient_fractions: list[float] = []
        n_transient_per_rep: list[float] = []
        total_chap_events: list[float] = []
        total_unassisted: list[float] = []

        for dyn in dynamics_per_rep:
            n_total = dyn.n_residues
            n_transient = dyn.n_transient()
            n_transient_per_rep.append(float(n_transient))
            transient_fractions.append(float(n_transient / n_total) if n_total > 0 else 0.0)

            transient_residues = dyn.transient_residues()
            if transient_residues:
                chap_frac = float(np.mean([r.chaperone_fraction for r in transient_residues]))
            else:
                chap_frac = 0.0
            chap_fractions.append(chap_frac)
            total_chap_events.append(float(dyn.total_chaperone_events()))
            total_unassisted.append(float(dyn.total_unassisted_events()))

        mean_chap = float(np.mean(chap_fractions))
        _chap_stats = compute_sem(chap_fractions) if n_replicates > 1 else None
        sem_chap = _chap_stats.sem if _chap_stats else 0.0
        mean_transient = float(np.mean(transient_fractions))
        _trans_stats = compute_sem(transient_fractions) if n_replicates > 1 else None
        sem_transient = _trans_stats.sem if _trans_stats else 0.0

        # Aggregate enrichment: mean residue-based enrichment per (polymer_type, aa_group)
        enrichment_by_ptype: dict[str, dict[str, float]] = {}
        polymer_types_set: set[str] = set()
        aa_groups_set: set[str] = set()

        for enr in enrichment_per_rep:
            for e in enr.entries:
                polymer_types_set.add(e.polymer_type)
                aa_groups_set.add(e.aa_group)

        polymer_types = sorted(polymer_types_set)
        aa_groups = sorted(aa_groups_set)

        for ptype in polymer_types:
            enrichment_by_ptype[ptype] = {}
            for ag in aa_groups:
                vals = []
                for enr in enrichment_per_rep:
                    entry = enr.get(ptype, ag)
                    if entry is not None:
                        vals.append(entry.enrichment_residue)
                enrichment_by_ptype[ptype][ag] = float(np.mean(vals)) if vals else float("nan")

        return ExposureConditionSummary(
            label=cond.label,
            config_path=str(cond.config_path),
            n_replicates=n_replicates,
            replicate_values=chap_fractions,
            mean_chaperone_fraction=mean_chap,
            sem_chaperone_fraction=sem_chap,
            mean_transient_fraction=mean_transient,
            sem_transient_fraction=sem_transient,
            mean_n_transient=float(np.mean(n_transient_per_rep)),
            mean_total_chaperone_events=float(np.mean(total_chap_events)),
            mean_total_unassisted_events=float(np.mean(total_unassisted)),
            enrichment_by_polymer_type=enrichment_by_ptype,
            polymer_types=polymer_types,
            aa_groups=aa_groups,
        )

    def _load_or_compute_replicate(
        self,
        cond: Condition,
        replicate: int,
        settings: ExposureSettings,
        exposure_analysis_dir: Path | None,
        contacts_analysis_dir: Path | None,
        recompute: bool,
    ) -> tuple[Any, Any] | None:
        """Load or compute exposure dynamics for a single replicate.

        Parameters
        ----------
        cond : Condition
            Condition being processed.
        replicate : int
            1-indexed replicate number.
        settings : ExposureSettings
            Plugin settings.
        exposure_analysis_dir : Path or None
            Base exposure analysis directory.
        contacts_analysis_dir : Path or None
            Contacts sibling directory.
        recompute : bool
            Force recompute.

        Returns
        -------
        tuple[ExposureDynamicsResult, ChaperoneEnrichmentResult] or None
        """
        from polyzymd.analyses.contacts._results import ContactResult
        from polyzymd.analyses.exposure._config import ExposureConfig
        from polyzymd.analyses.exposure._dynamics import (
            ExposureDynamicsResult,
            analyze_exposure_dynamics,
        )
        from polyzymd.analyses.exposure._enrichment import compute_chaperone_enrichment
        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import compute_trajectory_sasa
        from polyzymd.analyses.shared.loader import TrajectoryLoader

        # Find cached ContactResult
        contact_result_path = self._find_contact_result(contacts_analysis_dir, replicate)
        if contact_result_path is None or not contact_result_path.exists():
            logger.warning(
                f"    Contacts not found for rep {replicate}. Run contacts analysis first."
            )
            return None

        contact_result = ContactResult.load(contact_result_path)
        logger.info(f"    Loaded contacts for rep {replicate}: {contact_result_path}")

        # Load trajectory
        loader = TrajectoryLoader(cond.sim_config)
        try:
            traj_info = loader.get_trajectory_info(replicate)
        except FileNotFoundError as e:
            logger.warning(f"    Skipping rep {replicate}: trajectory not found. {e}")
            return None

        try:
            topology_path = traj_info.topology_file
            trajectory_paths = traj_info.trajectory_files

            # Analysis dir for this replicate — use run_N convention
            # (matches orchestrator's output_dir layout used by all plugins)
            if exposure_analysis_dir is None:
                raise ValueError(
                    "exposure_analysis_dir is required — ensure the exposure plugin is "
                    "configured with an output directory in comparison.yaml"
                )
            analysis_dir = exposure_analysis_dir / f"run_{replicate}"

            # SASA config
            sasa_config = SASAConfig(
                exposure_threshold=settings.exposure_threshold,
                probe_radius_nm=settings.probe_radius_nm,
                n_sphere_points=settings.n_sphere_points,
                chain_id=settings.protein_chain,
                cache_sasa=True,
            )

            # Check cached dynamics
            dynamics_cache_path = ExposureDynamicsResult.cache_path(analysis_dir)
            if not recompute and dynamics_cache_path.exists():
                logger.info(f"    Loading cached dynamics: {dynamics_cache_path}")
                dynamics = ExposureDynamicsResult.load(dynamics_cache_path)
                # Still need enrichment (always recomputed)
                sasa_result = compute_trajectory_sasa(
                    topology_path=topology_path,
                    trajectory_path=trajectory_paths,
                    config=sasa_config,
                    analysis_dir=analysis_dir,
                    recompute=False,
                )
                enrichment = compute_chaperone_enrichment(
                    sasa_result=sasa_result,
                    contact_result=contact_result,
                    polymer_resnames=settings.polymer_resnames,
                )
                return dynamics, enrichment

            # Compute SASA
            logger.info(f"    Computing SASA for rep {replicate}...")
            sasa_result = compute_trajectory_sasa(
                topology_path=topology_path,
                trajectory_path=trajectory_paths,
                config=sasa_config,
                analysis_dir=analysis_dir,
                recompute=recompute,
            )

            # Compute exposure dynamics
            exposure_config = ExposureConfig(
                transient_lower=settings.transient_lower,
                transient_upper=settings.transient_upper,
                min_event_length=settings.min_event_length,
            )

            logger.info(f"    Analyzing exposure dynamics for rep {replicate}...")
            dynamics = analyze_exposure_dynamics(
                sasa_result=sasa_result,
                contact_result=contact_result,
                config=exposure_config,
                analysis_dir=analysis_dir,
                recompute=recompute,
            )

            # Compute enrichment
            enrichment = compute_chaperone_enrichment(
                sasa_result=sasa_result,
                contact_result=contact_result,
                polymer_resnames=settings.polymer_resnames,
            )

            return dynamics, enrichment
        except (FileNotFoundError, ValueError, OSError, IndexError) as e:
            logger.debug("Full traceback:", exc_info=True)
            logger.warning(f"    Skipping rep {replicate}: analysis failed with error: {e}")
            return None

    @staticmethod
    def _find_contact_result(
        contacts_dir: Path | None,
        replicate: int,
    ) -> Path | None:
        """Find cached contact result JSON for a replicate.

        Only checks the condition-specific contacts directory — does NOT
        fall back to ``projects_directory`` to prevent cross-contamination
        between conditions.

        Parameters
        ----------
        contacts_dir : Path or None
            Condition-specific contacts directory.
        replicate : int
            Replicate number.

        Returns
        -------
        Path or None
        """

        if contacts_dir is None:
            return None

        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        resolved = find_contact_result_for_replicate(contacts_dir, replicate)
        if resolved is not None:
            return resolved

        run_dir = contacts_dir / f"run_{replicate}"
        fallback_patterns = (
            (contacts_dir, f"contacts*rep{replicate}.json"),
            (run_dir, f"contacts*rep{replicate}.json"),
        )
        for search_dir, pattern in fallback_patterns:
            matches = sorted(p for p in search_dir.glob(pattern) if p.is_file())
            if len(matches) == 1:
                return matches[0]
            if len(matches) > 1:
                raise ValueError(
                    f"Ambiguous contacts cache for replicate {replicate}: "
                    f"found {len(matches)} files matching '{pattern}' in {search_dir}: "
                    + ", ".join(str(match.name) for match in matches)
                )

        return None

    def _condition_has_polymer(self, cond: Condition) -> bool:
        """Check if a condition has polymer configured.

        Parameters
        ----------
        cond : Condition
            Condition to check.

        Returns
        -------
        bool
            True if polymer is configured and enabled.
        """
        sim_config = cond.sim_config
        if sim_config.polymers is None:
            return False
        if hasattr(sim_config.polymers, "enabled") and not sim_config.polymers.enabled:
            return False
        return True

    @staticmethod
    def _compare_pair(
        cond_a: Any,
        cond_b: Any,
        *,
        fdr_alpha: float = 0.05,
        ttest_method: str = "student",
        posthoc_method: str = "ttest_bh",
    ) -> Any:
        """Compare two condition summaries statistically.

        Parameters
        ----------
        cond_a : ExposureConditionSummary
            First condition (typically control).
        cond_b : ExposureConditionSummary
            Second condition (typically treatment).

        Returns
        -------
        PairwiseResult
        """
        from polyzymd.analyses.base import PairwiseResult
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )
        from polyzymd.analyses.stats import interpret_direction

        values_a = cond_a.replicate_values
        values_b = cond_b.replicate_values
        mean_a = cond_a.mean_chaperone_fraction
        mean_b = cond_b.mean_chaperone_fraction

        ttest = independent_ttest(values_a, values_b, method=ttest_method)
        effect = cohens_d(values_a, values_b)
        pct = percent_change(mean_a, mean_b)

        # Direction: higher chaperone fraction = more polymer-assisted events
        direction = interpret_direction(
            pct,
            direction_labels=("decreased", "unchanged", "increased"),
            threshold=5.0,
        )

        return PairwiseResult(
            condition_a=cond_a.label,
            condition_b=cond_b.label,
            metric="chaperone_fraction",
            t_statistic=ttest.t_statistic,
            p_value=ttest.p_value,
            posthoc_method=posthoc_method,
            cohens_d=effect.cohens_d,
            effect_size_interpretation=effect.interpretation,
            direction=direction,
            significant=ttest.p_value <= fdr_alpha,
            percent_change=pct,
        )


def _find_exposure_comparison_result(
    data: dict[str, Any],
    labels: Sequence[str],
) -> Any | None:
    """Locate a saved ExposureComparisonResult JSON.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition labels to plot input data
    labels : Sequence[str]
        Condition labels to search for

    Returns
    -------
    Any | None
        Loaded comparison result if found
    """
    from polyzymd.analyses.exposure._comparison_results import ExposureComparisonResult
    from polyzymd.analyses.shared.result_io import find_comparison_result

    return find_comparison_result(
        data,
        labels,
        glob_patterns=["exposure_comparison*.json"],
        loader=ExposureComparisonResult.load,
        analysis_type="exposure",
        fallback_filenames=["exposure_comparison.json", "comparison_result.json"],
        log=logger,
    )


def _plot_chaperone_fraction(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate chaperone fraction bar chart.

    Parameters
    ----------
    data : dict[str, Any]
        Plot input data keyed by condition label
    labels : Sequence[str]
        Ordered condition labels
    output_dir : Path
        Directory to save figure into
    plot_settings : Any
        Plot settings model with theme and output fields

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt
    import numpy as np

    result = _find_exposure_comparison_result(data, labels)
    if result is None:
        logger.warning(
            "No ExposureComparisonResult found; skipping chaperone fraction plot. "
            "Run exposure comparison first"
        )
        return []

    t = plot_settings.theme
    conditions = result.conditions
    n = len(conditions)

    cond_labels = [c.label for c in conditions]
    means = [c.mean_chaperone_fraction for c in conditions]
    sems = [c.sem_chaperone_fraction for c in conditions]
    colors = get_colors(n, plot_settings)

    fig, ax = plt.subplots(figsize=(max(6, n * 1.4), 5))
    x = np.arange(n)
    ax.bar(
        x,
        means,
        yerr=sems,
        capsize=t.bar_capsize,
        color=colors,
        alpha=t.bar_alpha,
        edgecolor=t.bar_edgecolor,
        linewidth=t.bar_linewidth,
    )

    rng = np.random.default_rng(seed=42)
    bar_width = 0.8
    for i, cond in enumerate(conditions):
        rep_vals = getattr(cond, "replicate_values", None)
        if rep_vals:
            rep_arr = np.asarray(rep_vals, dtype=float)
            jitter = rng.uniform(-bar_width * 0.25, bar_width * 0.25, size=len(rep_arr))
            ax.scatter(
                np.full_like(rep_arr, float(x[i])) + jitter,
                rep_arr,
                color=t.dot_color,
                s=t.dot_size,
                zorder=5,
                alpha=t.dot_alpha,
                edgecolors="none",
            )

    ax.set_xticks(x)
    ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
    apply_axis_style(
        ax,
        plot_settings,
        title="Chaperone fraction across conditions\n(transient residues only)",
        ylabel="Mean chaperone fraction",
    )
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    output_path = get_output_path(output_dir, "exposure_chaperone_fraction", plot_settings)
    return [
        save_figure(
            fig,
            output_path,
            plot_settings,
            experimental_features=("exposure",),
        )
    ]


def _plot_enrichment_heatmap(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate enrichment heatmap.

    Parameters
    ----------
    data : dict[str, Any]
        Plot input data keyed by condition label
    labels : Sequence[str]
        Ordered condition labels
    output_dir : Path
        Directory to save figure into
    plot_settings : Any
        Plot settings model with theme and output fields

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt
    import numpy as np

    result = _find_exposure_comparison_result(data, labels)
    if result is None:
        logger.warning(
            "No ExposureComparisonResult found; skipping enrichment heatmap. "
            "Run exposure comparison first"
        )
        return []

    t = plot_settings.theme
    conditions = result.conditions

    all_ptypes: list[str] = sorted({pt for c in conditions for pt in c.polymer_types})
    all_groups: list[str] = sorted({ag for c in conditions for ag in c.aa_groups})

    if not all_ptypes or not all_groups:
        logger.warning("No enrichment data to plot")
        return []

    n_conds = len(conditions)
    n_ptypes = len(all_ptypes)
    n_groups = len(all_groups)

    matrices = np.full((n_conds, n_ptypes, n_groups), np.nan)
    for ci, cond in enumerate(conditions):
        for pi, pt in enumerate(all_ptypes):
            for gi, ag in enumerate(all_groups):
                val = cond.enrichment_by_polymer_type.get(pt, {}).get(ag, float("nan"))
                matrices[ci, pi, gi] = val

    finite_vals = matrices[np.isfinite(matrices)]
    if len(finite_vals) == 0:
        logger.warning("All enrichment values are NaN; skipping heatmap")
        return []

    floor = 0.1
    vmax_raw = max(abs(finite_vals.min()), abs(finite_vals.max()), floor)
    vmin, vmax = -vmax_raw, vmax_raw

    fig_width = max(8, n_groups * 1.2 + 2)
    fig_height = max(4, n_ptypes * 0.8 * n_conds + 1)
    fig, axes = plt.subplots(
        1, n_conds, figsize=(fig_width, fig_height), sharey=True, squeeze=False
    )

    im = None
    for ci, (cond, ax) in enumerate(zip(conditions, axes[0])):
        mat = matrices[ci]
        im = ax.imshow(mat, vmin=vmin, vmax=vmax, cmap="RdBu_r", aspect="auto")
        ax.set_xticks(range(n_groups))
        ax.set_xticklabels(all_groups, rotation=45, ha="right", fontsize=t.tick_fontsize)
        ax.set_title(cond.label, fontsize=t.title_fontsize, fontweight=t.title_fontweight)
        if ci == 0:
            ax.set_yticks(range(n_ptypes))
            ax.set_yticklabels(all_ptypes, fontsize=t.tick_fontsize)
        else:
            ax.set_yticks([])

        annotate_cells(
            ax,
            mat,
            plot_settings,
            fmt="+.2f",
            fontsize=t.small_fontsize,
            threshold=vmax * 0.6,
            show_sign=False,
        )

    if im is not None:
        cbar = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.04, pad=0.04)
        cbar.set_label("Chaperone enrichment (residue-based)", fontsize=t.legend_fontsize)

    fig.suptitle(
        "Dynamic chaperone enrichment by AA group",
        fontsize=t.suptitle_fontsize,
        y=1.01,
    )
    fig.tight_layout()

    output_path = get_output_path(output_dir, "exposure_enrichment_heatmap", plot_settings)
    return [
        save_figure(
            fig,
            output_path,
            plot_settings,
            experimental_features=("exposure",),
        )
    ]
