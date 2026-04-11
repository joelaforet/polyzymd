"""Secondary structure analysis plugin.

Computes per-residue DSSP secondary structure from MD trajectories (via mdtraj),
aggregates persistence fractions across replicates, compares conditions via the
default scalar comparison pipeline (primary metric: helix fraction), and produces
comparison plots (timeline heatmaps, content bars, persistence diff heatmaps).

All heavy computation is self-contained — no delegation to legacy
``analysis.secondary_structure`` calculator classes.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.secondary_structure._plot_settings import SSPlotSettings
from polyzymd.analyses.secondary_structure._plotters import (
    _plot_ss_content_bars,
    _plot_ss_individual_bars,
    _plot_ss_persistence_diff_heatmap,
    _plot_ss_timeline_heatmap,
)
from polyzymd.analyses.secondary_structure._results import (
    SecondaryStructureAggregatedResult,
)
from polyzymd.analyses.shared.config_hash import compute_config_hash
from polyzymd.analyses.shared.loader import TrajectoryLoader, convert_time, parse_time_string

if TYPE_CHECKING:
    from polyzymd.analyses.secondary_structure._results import (
        SecondaryStructureResult,
    )

logger = logging.getLogger("polyzymd.analyses.secondary_structure")

# Canonical mapping between DSSP letters and integer codes
SS_CHAR_TO_INT: dict[str, int] = {"C": 0, "H": 1, "E": 2}


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class SecondaryStructureSettings(BaseModel):
    """Settings for secondary structure analysis.

    Attributes
    ----------
    chain_id : str
        Chain letter for the protein chain.  Default is ``"A"``
        (PolyzyMD convention: chain A = protein).
    """

    chain_id: str = Field(
        default="A",
        description="Chain letter for the protein chain",
    )


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class SecondaryStructureAnalysis(Analysis):
    """Secondary structure analysis: DSSP from MD trajectories.

    Performs the complete secondary structure workflow inline:

    1. Load trajectory from config (via mdtraj)
    2. Select protein chain and skip equilibration frames
    3. Compute DSSP assignments using ``mdtraj.compute_dssp(simplified=True)``
    4. Encode assignments as integer matrix and compute persistence fractions
    5. Save per-replicate result (JSON + NPZ sidecar)
    6. Aggregate across replicates with mean/SEM of persistence

    The ``compare()`` method is NOT overridden — it uses the default
    implementation which calls ``extract_metrics()`` to get
    ``helix_fraction`` as a single scalar metric with
    ``higher_is_better=True`` (more helix = more structured = stabilising).

    Plots
    -----
    - ``ss_timeline_<condition>.png``: Per-condition residue x time heatmap.
    - ``ss_content_bars.png``: Grouped bar chart of helix/strand/coil fractions.
    - ``ss_helix_bars.png``, ``ss_strand_bars.png``, ``ss_coil_bars.png``:
      Individual bar charts per SS type.
    - ``ss_persistence_diff_heatmap.png``: Delta(helix persistence) vs control.
    """

    name: ClassVar[str] = "secondary_structure"
    Settings: ClassVar[type] = SecondaryStructureSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = SSPlotSettings
    AggregatedResultClass: ClassVar[type] = SecondaryStructureAggregatedResult
    aliases: ClassVar[tuple[str, ...]] = ("ss",)
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute DSSP secondary structure for a single replicate.

        Performs trajectory loading (via mdtraj), chain selection,
        equilibration skipping, DSSP calculation, and persistence
        fraction computation inline.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        SecondaryStructureResult
            Per-replicate secondary structure result.
        """
        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.secondary_structure._results import SecondaryStructureResult

        settings = ctx.settings
        sim_config = ctx.sim_config

        chain_id = settings.chain_id

        # Parse equilibration time
        eq_value, eq_unit = parse_time_string(ctx.equilibration)

        # Initialize loader and config hash
        loader = TrajectoryLoader(sim_config)
        config_hash = compute_config_hash(sim_config)

        # Determine output path and check cache
        output_dir = ctx.output_dir
        eq_str = f"eq{eq_value:.2f}{eq_unit}"
        result_prefix = f"secondary_structure_{eq_str}"
        result_file = output_dir / f"{result_prefix}.json"

        cached = self._check_cache(
            SecondaryStructureResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

        logger.info(f"Computing secondary structure for replicate {replicate}")

        # =================================================================
        # Load trajectory with mdtraj
        # =================================================================
        import mdtraj as md

        traj_info = loader.get_trajectory_info(replicate)
        traj_files_str = [str(f) for f in traj_info.trajectory_files]
        topology_path = str(traj_info.topology_file)

        traj = md.load(traj_files_str, top=topology_path)
        logger.info(f"Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")

        # =================================================================
        # Select protein chain
        # =================================================================
        chain_idx = _chain_letter_to_index(chain_id)
        protein_indices = traj.topology.select(f"chainid {chain_idx}")
        if len(protein_indices) == 0:
            raise ValueError(
                f"Chain '{chain_id}' (index {chain_idx}) not found in topology. "
                f"Available chains: {[c.index for c in traj.topology.chains]}. "
                "Check your Settings.chain_id configuration."
            )

        protein_traj = traj.atom_slice(protein_indices)
        logger.info(
            f"Protein sub-trajectory: {protein_traj.n_atoms} atoms, "
            f"{protein_traj.n_residues} residues"
        )

        # =================================================================
        # Skip equilibration frames
        # =================================================================
        n_frames_total = protein_traj.n_frames
        # Use TrajectoryLoader for reliable timestep (mdtraj often has
        # incorrect time metadata for legacy DCD files)
        timestep_ps = loader.get_timestep(replicate)
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")
        start_frame = int(eq_time_ps / timestep_ps) if timestep_ps > 0 else 0
        start_frame = min(start_frame, n_frames_total - 1)

        if start_frame > 0:
            protein_traj = protein_traj[start_frame:]
            logger.info(
                f"Skipped {start_frame} equilibration frames "
                f"({eq_time_ps / 1000:.1f} ns), "
                f"{protein_traj.n_frames} frames remaining"
            )

        n_frames = protein_traj.n_frames

        # =================================================================
        # Compute DSSP
        # =================================================================
        dssp_raw = md.compute_dssp(protein_traj, simplified=True)
        # dssp_raw is (n_frames, n_residues) of single-char strings: H, E, C

        logger.info(f"DSSP complete: {dssp_raw.shape[0]} frames x {dssp_raw.shape[1]} residues")

        # =================================================================
        # Integer-encode the matrix
        # =================================================================
        ss_matrix = _encode_dssp_matrix(dssp_raw)
        n_residues = ss_matrix.shape[1]

        # =================================================================
        # Compute persistence fractions (per-residue)
        # =================================================================
        persistence_coil = (ss_matrix == 0).mean(axis=0).tolist()
        persistence_helix = (ss_matrix == 1).mean(axis=0).tolist()
        persistence_strand = (ss_matrix == 2).mean(axis=0).tolist()

        # =================================================================
        # Compute overall content fractions
        # =================================================================
        total_entries = ss_matrix.size
        overall_coil = float(np.sum(ss_matrix == 0)) / total_entries
        overall_helix = float(np.sum(ss_matrix == 1)) / total_entries
        overall_strand = float(np.sum(ss_matrix == 2)) / total_entries

        # =================================================================
        # Collect residue metadata
        # =================================================================
        protein_top = protein_traj.topology
        residue_ids = [res.resSeq for res in protein_top.residues]
        residue_names = [res.name.upper() for res in protein_top.residues]

        # =================================================================
        # Build result
        # =================================================================
        result = SecondaryStructureResult(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string=f"chainid {chain_idx} (chain {chain_id})",
            n_frames=n_frames,
            n_residues=n_residues,
            residue_ids=residue_ids,
            residue_names=residue_names,
            persistence_coil=persistence_coil,
            persistence_helix=persistence_helix,
            persistence_strand=persistence_strand,
            overall_helix_fraction=overall_helix,
            overall_strand_fraction=overall_strand,
            overall_coil_fraction=overall_coil,
        )

        # Save JSON + NPZ sidecar
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        json_path = result.save_with_matrix(
            directory=output_dir,
            matrix=ss_matrix,
            filename_prefix=result_prefix,
        )
        logger.info(f"Saved SS result to {json_path}")

        return result

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate secondary structure across replicates for one condition.

        Computes mean and SEM of persistence fractions and overall content
        fractions from the already-computed per-replicate results.  Does NOT
        re-run the calculator.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[SecondaryStructureResult]
            Per-replicate secondary structure results.

        Returns
        -------
        SecondaryStructureAggregatedResult
            Aggregated result with mean/SEM across replicates.
        """
        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.secondary_structure._results import (
            SecondaryStructureAggregatedResult,
        )

        first = results[0]
        n_reps = len(results)

        # Stack persistence arrays: shape (n_reps, n_residues)
        coil_stack = np.array([r.persistence_coil for r in results])
        helix_stack = np.array([r.persistence_helix for r in results])
        strand_stack = np.array([r.persistence_strand for r in results])

        mean_coil = coil_stack.mean(axis=0).tolist()
        mean_helix = helix_stack.mean(axis=0).tolist()
        mean_strand = strand_stack.mean(axis=0).tolist()

        sem_coil = (coil_stack.std(axis=0, ddof=1) / np.sqrt(n_reps)).tolist()
        sem_helix = (helix_stack.std(axis=0, ddof=1) / np.sqrt(n_reps)).tolist()
        sem_strand = (strand_stack.std(axis=0, ddof=1) / np.sqrt(n_reps)).tolist()

        # Overall content fractions per replicate
        per_rep_helix = [r.overall_helix_fraction for r in results]
        per_rep_strand = [r.overall_strand_fraction for r in results]
        per_rep_coil = [r.overall_coil_fraction for r in results]

        mean_overall_helix = float(np.mean(per_rep_helix))
        mean_overall_strand = float(np.mean(per_rep_strand))
        mean_overall_coil = float(np.mean(per_rep_coil))

        sem_overall_helix = float(np.std(per_rep_helix, ddof=1) / np.sqrt(n_reps))
        sem_overall_strand = float(np.std(per_rep_strand, ddof=1) / np.sqrt(n_reps))
        sem_overall_coil = float(np.std(per_rep_coil, ddof=1) / np.sqrt(n_reps))

        agg_result = SecondaryStructureAggregatedResult(
            config_hash=first.config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=first.equilibration_time,
            equilibration_unit=first.equilibration_unit,
            selection_string=first.selection_string,
            n_replicates=n_reps,
            replicates=list(ctx.replicates),
            residue_ids=first.residue_ids,
            residue_names=first.residue_names,
            mean_persistence_coil=mean_coil,
            mean_persistence_helix=mean_helix,
            mean_persistence_strand=mean_strand,
            sem_persistence_coil=sem_coil,
            sem_persistence_helix=sem_helix,
            sem_persistence_strand=sem_strand,
            mean_overall_helix=mean_overall_helix,
            mean_overall_strand=mean_overall_strand,
            mean_overall_coil=mean_overall_coil,
            sem_overall_helix=sem_overall_helix,
            sem_overall_strand=sem_overall_strand,
            sem_overall_coil=sem_overall_coil,
            per_replicate_helix=per_rep_helix,
            per_replicate_strand=per_rep_strand,
            per_replicate_coil=per_rep_coil,
        )

        target_path = ctx.result_path
        if target_path is None:
            filename = self._make_aggregated_filename(ctx.replicates, first)
            target_path = ctx.output_dir / filename
        self.save_result(agg_result, target_path)
        logger.info(f"Saved aggregated SS result to {target_path}")

        return agg_result

    # === Optional methods ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract helix fraction for scalar comparison.

        The primary metric is helix fraction because helix loss is the
        most common signature of thermal unfolding in globular enzymes.

        Parameters
        ----------
        summary : SecondaryStructureAggregatedResult
            Aggregated secondary structure result.

        Returns
        -------
        dict[str, MetricValue]
            Single entry ``"helix_fraction"`` with
            ``higher_is_better=True`` (more helix = more stable).
        """
        return {
            "helix_fraction": MetricValue(
                name="helix_fraction",
                mean=summary.mean_overall_helix,
                sem=summary.sem_overall_helix,
                replicate_values=summary.per_replicate_helix,
                higher_is_better=True,
                direction_labels=("destabilizing", "unchanged", "stabilizing"),
            ),
        }

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format secondary structure comparison result for CLI display.

        Parameters
        ----------
        result : ComparisonResult or BaseModel
            Comparison result to format.
        output_format : str
            ``"text"``, ``"markdown"``, or ``"json"``.

        Returns
        -------
        str
            Formatted output.
        """
        from polyzymd.analyses.base import ComparisonResult
        from polyzymd.analyses.stats import format_scalar_comparison

        if isinstance(result, ComparisonResult):
            return format_scalar_comparison(
                result,
                title="Secondary Structure Comparison",
                metric_label="Helix Fraction",
                metric_unit="",
                metric_key="helix_fraction",
                output_format=output_format,
                higher_is_better=True,
            )
        return super().format(result, output_format)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate secondary structure comparison plots.

        Produces:

        - ``ss_timeline_<condition>.png``: Per-condition residue x time heatmap.
        - ``ss_content_bars.png``: Grouped bar chart of helix/strand/coil fractions.
        - ``ss_helix_bars.png``, ``ss_strand_bars.png``, ``ss_coil_bars.png``:
          Individual bar charts per SS type.
        - ``ss_persistence_diff_heatmap.png``: Delta(helix persistence) vs control.

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

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        plot_settings = ctx.plot_settings

        result = _plot_ss_timeline_heatmap(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        result = _plot_ss_content_bars(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        result = _plot_ss_individual_bars(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        result = _plot_ss_persistence_diff_heatmap(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        return plots

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Backward-compatible filename helper retained for tests."""
        eq_str = f"eq{first_result.equilibration_time:.2f}{first_result.equilibration_unit}"
        rep_str = Analysis._format_replicate_range(replicates)
        return f"secondary_structure_{rep_str}_{eq_str}.json"


# ---------------------------------------------------------------------------
# Private helper functions (extracted from legacy calculator)
# ---------------------------------------------------------------------------


def _chain_letter_to_index(chain_id: str) -> int:
    """Convert chain letter (A, B, C...) to 0-indexed integer for mdtraj."""
    return ord(chain_id.upper()) - ord("A")


def _encode_dssp_matrix(dssp_raw: np.ndarray) -> np.ndarray:
    """Convert character DSSP matrix to integer encoding.

    Parameters
    ----------
    dssp_raw : np.ndarray
        Character array from ``mdtraj.compute_dssp()``, shape
        ``(n_frames, n_residues)`` with values ``'H'``, ``'E'``,
        ``'C'`` (or ``'NA'`` for non-protein residues).

    Returns
    -------
    np.ndarray
        Integer-encoded matrix, shape ``(n_frames, n_residues)``,
        dtype ``int8``.  Encoding: 0=C, 1=H, 2=E.
    """
    n_frames, n_residues = dssp_raw.shape
    matrix = np.zeros((n_frames, n_residues), dtype=np.int8)

    for char, code in SS_CHAR_TO_INT.items():
        matrix[dssp_raw == char] = code

    return matrix
