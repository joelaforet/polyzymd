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
- ``plot()`` delegates to 2 existing plotters
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
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
    from polyzymd.analysis.contacts.results import ContactResult
    from polyzymd.analysis.exposure.dynamics import ExposureDynamicsResult
    from polyzymd.analysis.exposure.enrichment import ChaperoneEnrichmentResult
    from polyzymd.analysis.sasa.trajectory import SASATrajectoryResult

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
    Settings: ClassVar[type] = ExposureSettings
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ("contacts",)
    min_replicates: ClassVar[int] = 2

    # === Compute (no-op — comparator-only) ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> None:
        """No-op — exposure dynamics is comparator-only.

        All per-replicate computation (SASA, dynamics, enrichment) is
        orchestrated within ``compare()`` because the pipeline requires
        both SASA and contacts data simultaneously.

        Returns
        -------
        None
        """
        return None

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> None:
        """No-op — exposure dynamics is comparator-only.

        Returns
        -------
        None
        """
        return None

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
            Comparison result, or ``None`` if fewer than 2 valid conditions.
        """
        from polyzymd import __version__
        from polyzymd.analysis.core.statistics import compute_sem
        from polyzymd.compare.core.base import ANOVASummary, PairwiseComparison
        from polyzymd.compare.results.exposure import (
            ExposureComparisonResult,
            ExposureConditionSummary,
        )
        from polyzymd.compare.statistics import (
            cohens_d,
            independent_ttest,
            one_way_anova,
            percent_change,
        )

        settings = ctx.settings

        logger.info("Starting exposure dynamics comparison")
        logger.info(f"  Conditions: {len(ctx.conditions)}")
        logger.info(f"  Equilibration: {ctx.equilibration}")

        conditions = ctx.conditions
        if len(conditions) < 2:
            logger.warning("Fewer than 2 conditions — skipping exposure comparison")
            return None

        # Step 1: Load/compute per condition
        summaries: list[ExposureConditionSummary] = []
        for cond in conditions:
            try:
                summary = self._build_condition_summary(cond, ctx, settings)
                summaries.append(summary)
            except Exception as e:
                logger.warning(f"  Skipping condition '{cond.label}': {e}")

        if len(summaries) < 2:
            logger.warning(
                f"Only {len(summaries)} condition(s) with valid data — "
                "need at least 2 for comparison"
            )
            return None

        # Step 2: Determine effective control
        effective_control = ctx.effective_control

        # Step 3: Pairwise comparisons on chaperone_fraction
        comparisons: list[PairwiseComparison] = []
        if effective_control:
            control_summary = next((s for s in summaries if s.label == effective_control), None)
            if control_summary:
                treatments = [s for s in summaries if s.label != effective_control]
                for treatment in treatments:
                    comp = self._compare_pair(control_summary, treatment)
                    comparisons.append(comp)
        else:
            for i, cond_a in enumerate(summaries):
                for cond_b in summaries[i + 1 :]:
                    comp = self._compare_pair(cond_a, cond_b)
                    comparisons.append(comp)

        # Step 4: ANOVA if 3+ conditions
        anova: ANOVASummary | None = None
        if len(summaries) >= 3:
            all_groups = [s.replicate_values for s in summaries]
            result = one_way_anova(*all_groups)
            anova = ANOVASummary(
                metric="chaperone_fraction",
                f_statistic=result.f_statistic,
                p_value=result.p_value,
                significant=result.significant,
            )

        # Step 5: Dual rankings
        ranked_chaperone = sorted(summaries, key=lambda s: s.mean_chaperone_fraction, reverse=True)
        ranked_transient = sorted(summaries, key=lambda s: s.mean_transient_fraction, reverse=True)

        # Build excluded list from context
        excluded_labels = (
            [
                c.label
                for c in getattr(ctx, "_all_conditions", [])
                if c.label not in {s.label for s in summaries}
            ]
            if hasattr(ctx, "_all_conditions")
            else []
        )

        return ExposureComparisonResult(
            name="exposure_comparison",
            metric="chaperone_fraction",
            control_label=effective_control,
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

    def filter_conditions(self, conditions: list[Condition]) -> list[Condition]:
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
        """Generate exposure comparison plots.

        Delegates to:
        - :class:`~polyzymd.compare.plotters.exposure.ExposureChaperoneFractionPlotter`
        - :class:`~polyzymd.compare.plotters.exposure.ExposureEnrichmentHeatmapPlotter`

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

        # Add __meta__ for results_dir
        data["__meta__"] = {"results_dir": ctx.results_dir}

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve plot settings
        plot_settings = ctx.plot_settings
        if plot_settings is None:
            from polyzymd.compare.config import PlotSettings

            plot_settings = PlotSettings()

        plotter_specs: list[tuple[str, str]] = [
            ("polyzymd.compare.plotters.exposure", "ExposureChaperoneFractionPlotter"),
            ("polyzymd.compare.plotters.exposure", "ExposureEnrichmentHeatmapPlotter"),
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

    # === extract_metrics (empty — full compare() override) ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Return empty dict — exposure uses full compare() override."""
        return {}

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
        from polyzymd.analysis.contacts.results import ContactResult
        from polyzymd.analysis.core.loader import TrajectoryLoader
        from polyzymd.analysis.core.statistics import compute_sem
        from polyzymd.analysis.exposure.config import ExposureConfig
        from polyzymd.analysis.exposure.dynamics import (
            ExposureDynamicsResult,
            analyze_exposure_dynamics,
        )
        from polyzymd.analysis.exposure.enrichment import compute_chaperone_enrichment
        from polyzymd.analysis.sasa.config import SASAConfig
        from polyzymd.analysis.sasa.trajectory import compute_trajectory_sasa
        from polyzymd.compare.results.exposure import ExposureConditionSummary

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
        from polyzymd.analysis.contacts.results import ContactResult
        from polyzymd.analysis.core.loader import TrajectoryLoader
        from polyzymd.analysis.exposure.config import ExposureConfig
        from polyzymd.analysis.exposure.dynamics import (
            ExposureDynamicsResult,
            analyze_exposure_dynamics,
        )
        from polyzymd.analysis.exposure.enrichment import compute_chaperone_enrichment
        from polyzymd.analysis.sasa.config import SASAConfig
        from polyzymd.analysis.sasa.trajectory import compute_trajectory_sasa

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
            analysis_dir: Path
            if exposure_analysis_dir is not None:
                analysis_dir = exposure_analysis_dir / f"run_{replicate}"
            else:
                analysis_dir = Path(f"/tmp/exposure_run_{replicate}")

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
        except Exception as e:
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
        result_filename = f"contacts_rep{replicate}.json"

        if contacts_dir is not None:
            # Flat path (legacy layout)
            cond_path = contacts_dir / result_filename
            if cond_path.exists():
                return cond_path

            # Orchestrator layout: run_<N>/ subdirectory per replicate
            run_path = contacts_dir / f"run_{replicate}" / result_filename
            if run_path.exists():
                return run_path

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
        PairwiseComparison
        """
        from polyzymd.compare.core.base import PairwiseComparison
        from polyzymd.compare.statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )

        values_a = cond_a.replicate_values
        values_b = cond_b.replicate_values
        mean_a = cond_a.mean_chaperone_fraction
        mean_b = cond_b.mean_chaperone_fraction

        ttest = independent_ttest(values_a, values_b)
        effect = cohens_d(values_a, values_b)
        pct = percent_change(mean_a, mean_b)

        # Direction: higher chaperone fraction = more polymer-assisted events
        if pct > 5.0:
            direction = "increased"
        elif pct < -5.0:
            direction = "decreased"
        else:
            direction = "unchanged"

        return PairwiseComparison(
            condition_a=cond_a.label,
            condition_b=cond_b.label,
            metric="chaperone_fraction",
            t_statistic=ttest.t_statistic,
            p_value=ttest.p_value,
            cohens_d=effect.cohens_d,
            effect_size_interpretation=effect.interpretation,
            direction=direction,
            significant=ttest.significant,
            percent_change=pct,
        )

    def _deserialize_result(self, path: Path) -> Any:
        """Load an exposure comparison result from JSON.

        Parameters
        ----------
        path : Path
            Path to JSON file.

        Returns
        -------
        ExposureComparisonResult
        """
        from polyzymd.compare.results.exposure import ExposureComparisonResult

        return ExposureComparisonResult.load(path)
