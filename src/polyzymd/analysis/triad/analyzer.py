"""Catalytic triad analysis module.

This module provides the CatalyticTriadAnalyzer class for analyzing
catalytic triad/active site geometry from MD trajectories.

Key Features
------------
- Config-based trajectory loading via CatalyticTriadConfig
- Reuses DistanceCalculator for per-pair distance computation
- Computes simultaneous contact fraction (all pairs below threshold)
- Autocorrelation-corrected uncertainty for contact fraction
- Multi-replicate aggregation with SEM

The key metric is the "simultaneous contact fraction" - the percentage
of frames where ALL distance pairs are below the contact threshold at
the same time. This indicates triad integrity and catalytic competence.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

import numpy as np
from numpy.typing import NDArray

from polyzymd.analysis.core.autocorrelation import estimate_correlation_time
from polyzymd.analysis.core.config_hash import compute_config_hash, validate_config_hash
from polyzymd.analysis.core.loader import (
    TrajectoryLoader,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analysis.core.statistics import compute_sem
from polyzymd.analysis.distances.calculator import DistanceCalculator
from polyzymd.analysis.results.base import get_polyzymd_version
from polyzymd.analysis.results.distances import (
    DistancePairAggregatedResult,
    DistancePairResult,
)
from polyzymd.analysis.results.triad import TriadAggregatedResult, TriadResult
from polyzymd.compare.config import CatalyticTriadConfig

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig

LOGGER = logging.getLogger(__name__)


class CatalyticTriadAnalyzer:
    """Analyzer for catalytic triad/active site geometry.

    This class handles the complete triad analysis workflow:
    1. Load trajectories from SimulationConfig
    2. Compute per-pair distances (delegated to DistanceCalculator)
    3. Compute simultaneous contact fraction
    4. Apply autocorrelation correction for uncertainty
    5. Aggregate across replicates with SEM

    Parameters
    ----------
    config : SimulationConfig
        PolyzyMD simulation configuration
    triad_config : CatalyticTriadConfig
        Catalytic triad configuration with pairs and threshold
    equilibration : str, optional
        Equilibration time to skip. Default is "0ns".

    Examples
    --------
    >>> from polyzymd.config.schema import SimulationConfig
    >>> from polyzymd.compare.config import CatalyticTriadConfig, TriadPairConfig
    >>>
    >>> config = SimulationConfig.load("config.yaml")
    >>>
    >>> triad = CatalyticTriadConfig(
    ...     name="LipA_triad",
    ...     description="Ser-His-Asp catalytic triad",
    ...     threshold=3.5,
    ...     pairs=[
    ...         TriadPairConfig(
    ...             label="Asp133-His156",
    ...             selection_a="midpoint(resid 133 and name OD1 OD2)",
    ...             selection_b="resid 156 and name ND1",
    ...         ),
    ...         TriadPairConfig(
    ...             label="His156-Ser77",
    ...             selection_a="resid 156 and name NE2",
    ...             selection_b="resid 77 and name OG",
    ...         ),
    ...     ],
    ... )
    >>>
    >>> analyzer = CatalyticTriadAnalyzer(
    ...     config,
    ...     triad_config=triad,
    ...     equilibration="100ns",
    ... )
    >>>
    >>> result = analyzer.compute(replicate=1)
    >>> print(f"Simultaneous contact: {result.simultaneous_contact_fraction * 100:.1f}%")
    """

    def __init__(
        self,
        config: "SimulationConfig",
        triad_config: CatalyticTriadConfig,
        equilibration: str = "0ns",
    ) -> None:
        self.config = config
        self.triad_config = triad_config

        # Parse equilibration time
        eq_value, eq_unit = parse_time_string(equilibration)
        self.equilibration_time = eq_value
        self.equilibration_unit = eq_unit

        # Initialize loader and config hash
        self._loader = TrajectoryLoader(config)
        self._config_hash = compute_config_hash(config)

        # Create internal DistanceCalculator with triad pairs
        self._distance_calc = DistanceCalculator(
            config,
            pairs=triad_config.get_pair_selections(),
            equilibration=equilibration,
            threshold=triad_config.threshold,
        )

    def compute(
        self,
        replicate: int,
        save: bool = True,
        output_dir: Path | None = None,
        recompute: bool = False,
        store_timeseries: bool = False,
    ) -> TriadResult:
        """Compute triad analysis for a single replicate.

        Parameters
        ----------
        replicate : int
            Replicate number (1-indexed)
        save : bool, optional
            If True (default), save result to JSON
        output_dir : Path, optional
            Directory to save results
        recompute : bool, optional
            If True, recompute even if cached
        store_timeseries : bool, optional
            If True, store full simultaneous contact timeseries

        Returns
        -------
        TriadResult
            Triad analysis results with per-pair distances and
            simultaneous contact fraction
        """
        if output_dir is None:
            output_dir = (
                self.config.output.projects_directory / "analysis" / "triad" / f"run_{replicate}"
            )

        result_file = output_dir / self._make_result_filename()

        # Check cache
        if not recompute and result_file.exists():
            LOGGER.info(f"Loading cached triad result from {result_file}")
            result = TriadResult.load(result_file)
            validate_config_hash(result.config_hash, self.config)
            return result

        LOGGER.info(
            f"Computing triad analysis '{self.triad_config.name}' for replicate {replicate}"
        )

        # Compute per-pair distances using DistanceCalculator
        # Don't save individual pair results - we'll include them in the triad result
        distance_result = self._distance_calc.compute(
            replicate=replicate,
            save=False,
            recompute=recompute,
            store_distributions=True,  # Need full arrays for simultaneous contact
        )

        # Get pair results and update labels from triad config
        pair_results = []
        pair_distances_arrays = []

        for i, pr in enumerate(distance_result.pair_results):
            # Update pair label from triad config
            pr_updated = pr.model_copy(
                update={
                    "pair_label": self.triad_config.pairs[i].label,
                    "replicate": replicate,
                }
            )
            pair_results.append(pr_updated)

            # Collect distance arrays for simultaneous contact calculation
            if pr.distances is not None:
                pair_distances_arrays.append(np.array(pr.distances))
            else:
                raise ValueError(
                    f"Distance array not available for pair {i}. "
                    "This shouldn't happen - store_distributions=True was set."
                )

        # Compute simultaneous contact fraction
        n_frames_used = distance_result.n_frames_used
        threshold = self.triad_config.threshold

        # Frame-by-frame AND of all pairs below threshold
        all_below = np.ones(n_frames_used, dtype=bool)
        for dist_arr in pair_distances_arrays:
            all_below &= dist_arr < threshold

        simultaneous_fraction = float(all_below.mean())
        n_frames_simultaneous = int(all_below.sum())

        LOGGER.info(
            f"  Simultaneous contact: {simultaneous_fraction * 100:.1f}% "
            f"({n_frames_simultaneous}/{n_frames_used} frames)"
        )

        # Autocorrelation analysis for simultaneous contact timeseries
        # Convert boolean to float for ACF computation
        contact_timeseries = all_below.astype(np.float64)

        sim_contact_sem = None
        sim_contact_tau = None
        sim_contact_tau_unit = None
        sim_contact_n_ind = None
        sim_contact_warning = None

        if n_frames_used >= 20:
            try:
                # Get timestep for correlation time estimation
                timestep = self._loader.get_timestep(replicate, unit="ps")

                tau_result = estimate_correlation_time(
                    contact_timeseries,
                    timestep=timestep,
                    timestep_unit="ps",
                    method="integration",
                    n_frames=n_frames_used,
                )

                sim_contact_tau = tau_result.tau
                sim_contact_tau_unit = tau_result.tau_unit
                sim_contact_n_ind = tau_result.n_independent
                sim_contact_warning = tau_result.warning

                # Compute autocorrelation-corrected SEM for the fraction
                # For a proportion p with n_independent samples:
                # SEM = sqrt(p * (1-p) / n_independent)
                if sim_contact_n_ind > 0:
                    p = simultaneous_fraction
                    sim_contact_sem = float(np.sqrt(p * (1 - p) / sim_contact_n_ind))

                LOGGER.debug(
                    f"  Contact autocorrelation: τ={sim_contact_tau:.1f} {sim_contact_tau_unit}, "
                    f"n_ind={sim_contact_n_ind}, SEM={sim_contact_sem:.3f}"
                )

            except Exception as e:
                LOGGER.warning(f"Autocorrelation analysis for contact timeseries failed: {e}")
                # Fall back to naive SEM
                p = simultaneous_fraction
                sim_contact_sem = float(np.sqrt(p * (1 - p) / n_frames_used))

        # Create result
        result = TriadResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            selection_string="; ".join(
                f"({p.selection_a} : {p.selection_b})" for p in self.triad_config.pairs
            ),
            triad_name=self.triad_config.name,
            triad_description=self.triad_config.description,
            pair_results=pair_results,
            threshold=threshold,
            simultaneous_contact_fraction=simultaneous_fraction,
            n_frames_simultaneous=n_frames_simultaneous,
            simultaneous_contact_timeseries=(all_below.tolist() if store_timeseries else None),
            sim_contact_sem=sim_contact_sem,
            sim_contact_correlation_time=sim_contact_tau,
            sim_contact_correlation_time_unit=sim_contact_tau_unit,
            sim_contact_n_independent=sim_contact_n_ind,
            sim_contact_warning=sim_contact_warning,
            n_frames_total=distance_result.n_frames_total,
            n_frames_used=n_frames_used,
        )

        if save:
            result.save(result_file)
            LOGGER.info(f"Saved triad result to {result_file}")

        return result

    def compute_aggregated(
        self,
        replicates: Sequence[int],
        save: bool = True,
        output_dir: Path | None = None,
        recompute: bool = False,
    ) -> TriadAggregatedResult:
        """Compute aggregated triad analysis across replicates.

        Parameters
        ----------
        replicates : sequence of int
            Replicate numbers to aggregate
        save : bool, optional
            If True (default), save result
        output_dir : Path, optional
            Output directory
        recompute : bool, optional
            If True, recompute

        Returns
        -------
        TriadAggregatedResult
            Aggregated results with SEM across replicates
        """
        replicates = list(replicates)

        if output_dir is None:
            output_dir = self.config.output.projects_directory / "analysis" / "triad" / "aggregated"

        # Compute individual replicates
        individual_results = []
        for rep in replicates:
            result = self.compute(
                replicate=rep,
                save=save,
                recompute=recompute,
                store_timeseries=False,
            )
            individual_results.append(result)

        # Aggregate per-pair statistics
        n_pairs = len(self.triad_config.pairs)
        aggregated_pairs = []

        for pair_idx in range(n_pairs):
            pair_config = self.triad_config.pairs[pair_idx]

            # Collect per-replicate statistics
            per_rep_means = []
            per_rep_stds = []
            per_rep_medians = []
            per_rep_fractions = []
            per_rep_kde_peaks = []

            for result in individual_results:
                pr = result.pair_results[pair_idx]
                per_rep_means.append(pr.mean_distance)
                per_rep_stds.append(pr.std_distance)
                per_rep_medians.append(pr.median_distance)
                if pr.fraction_below_threshold is not None:
                    per_rep_fractions.append(pr.fraction_below_threshold)
                if pr.kde_peak is not None:
                    per_rep_kde_peaks.append(pr.kde_peak)

            # Compute aggregated statistics
            mean_stats = compute_sem(per_rep_means)
            median_stats = compute_sem(per_rep_medians)

            fraction_stats = None
            if per_rep_fractions:
                fraction_stats = compute_sem(per_rep_fractions)

            kde_peak_stats = None
            if per_rep_kde_peaks:
                kde_peak_stats = compute_sem(per_rep_kde_peaks)

            agg_pair = DistancePairAggregatedResult(
                config_hash=self._config_hash,
                polyzymd_version=get_polyzymd_version(),
                replicate=None,
                equilibration_time=self.equilibration_time,
                equilibration_unit=self.equilibration_unit,
                selection_string=f"{pair_config.selection_a} : {pair_config.selection_b}",
                replicates=replicates,
                n_replicates=len(replicates),
                pair_label=pair_config.label,
                selection1=pair_config.selection_a,
                selection2=pair_config.selection_b,
                overall_mean=mean_stats.mean,
                overall_sem=mean_stats.sem,
                overall_median=median_stats.mean,
                per_replicate_means=per_rep_means,
                per_replicate_stds=per_rep_stds,
                per_replicate_medians=per_rep_medians,
                threshold=self.triad_config.threshold,
                overall_fraction_below=(fraction_stats.mean if fraction_stats else None),
                sem_fraction_below=(fraction_stats.sem if fraction_stats else None),
                per_replicate_fractions_below=(per_rep_fractions if per_rep_fractions else None),
                overall_kde_peak=(kde_peak_stats.mean if kde_peak_stats else None),
                sem_kde_peak=(kde_peak_stats.sem if kde_peak_stats else None),
                per_replicate_kde_peaks=(per_rep_kde_peaks if per_rep_kde_peaks else None),
            )
            aggregated_pairs.append(agg_pair)

        # Aggregate simultaneous contact fraction
        per_rep_simultaneous = [r.simultaneous_contact_fraction for r in individual_results]
        sim_stats = compute_sem(per_rep_simultaneous)

        # Create aggregated result
        agg_result = TriadAggregatedResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            selection_string="; ".join(
                f"({p.selection_a} : {p.selection_b})" for p in self.triad_config.pairs
            ),
            replicates=replicates,
            n_replicates=len(replicates),
            triad_name=self.triad_config.name,
            triad_description=self.triad_config.description,
            pair_results=aggregated_pairs,
            threshold=self.triad_config.threshold,
            overall_simultaneous_contact=sim_stats.mean,
            sem_simultaneous_contact=sim_stats.sem,
            per_replicate_simultaneous=per_rep_simultaneous,
            source_result_files=[
                str(
                    self.config.output.projects_directory
                    / "analysis"
                    / "triad"
                    / f"run_{r.replicate}"
                    / self._make_result_filename()
                )
                for r in individual_results
            ],
        )

        if save:
            result_file = output_dir / self._make_aggregated_filename(replicates)
            agg_result.save(result_file)
            LOGGER.info(f"Saved aggregated triad result to {result_file}")

        return agg_result

    def _make_result_filename(self) -> str:
        """Generate filename for result JSON."""
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"
        # Sanitize triad name for filename
        name_safe = self.triad_config.name.replace(" ", "_").replace("/", "-")
        return f"triad_{name_safe}_{eq_str}.json"

    def _make_aggregated_filename(self, replicates: Sequence[int]) -> str:
        """Generate filename for aggregated result."""
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
        name_safe = self.triad_config.name.replace(" ", "_").replace("/", "-")
        return f"triad_{name_safe}_{rep_str}_{eq_str}.json"
