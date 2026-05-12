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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    MetricValue,
    PlotContext,
    ReplicateContext,
    SlurmResourceHint,
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
    SecondaryStructureResult,
)
from polyzymd.analyses.secondary_structure._runner import SecondaryStructureReplicateRunner
from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
from polyzymd.analyses.shared.loader import TrajectoryLoader, convert_time, parse_time_string

logger = logging.getLogger("polyzymd.analyses.secondary_structure")

# Canonical mapping between DSSP letters and integer codes
SS_CHAR_TO_INT: dict[str, int] = {"C": 0, "H": 1, "E": 2}


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class SecondaryStructureSettings(BaseModel):
    """Settings for secondary structure analysis.

    Chain ``A`` is the default protein chain under the PolyzyMD chain
    convention.
    """

    chain_id: str = Field(
        default="A",
        description="Chain letter for the protein chain",
    )


@dataclass(frozen=True)
class _SecondaryStructureTrajectoryWindow:
    """Trajectory window augmented with mdtraj file metadata."""

    base_window: Any
    trajectory_files: tuple[Path, ...]
    topology_file: Path

    @property
    def start(self) -> int:
        """Return the inclusive start frame."""

        return self.base_window.start

    @property
    def stop(self) -> int:
        """Return the exclusive stop frame."""

        return self.base_window.stop

    @property
    def step(self) -> int:
        """Return the frame stride."""

        return self.base_window.step

    @property
    def timestep_ps(self) -> float:
        """Return the trajectory timestep in picoseconds."""

        return self.base_window.timestep_ps

    @property
    def warning_message(self) -> str | None:
        """Return a non-fatal window warning, when present."""

        return self.base_window.warning_message

    def run_kwargs(self) -> dict[str, int]:
        """Return runner-compatible frame-window keyword arguments."""

        return self.base_window.run_kwargs()


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class SecondaryStructureAnalysis(Analysis):
    """Secondary structure analysis: DSSP from MD trajectories.

    Performs the complete secondary structure workflow through a dedicated
    runner-backed compute path:

    - Load trajectory metadata through the shared loader
    - Select protein chain and resolve the equilibration frame window
    - Compute DSSP assignments using ``mdtraj.compute_dssp(simplified=True)``
    - Encode assignments as integer matrix and compute persistence fractions
    - Save per-replicate result (JSON + NPZ sidecar)
    - Aggregate across replicates with mean/SEM of persistence

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
    ReplicateResultClass: ClassVar[type | None] = SecondaryStructureResult
    aliases: ClassVar[tuple[str, ...]] = ("ss",)
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 1
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = SlurmResourceHint(mem="16G")

    # === Required methods ===

    def run_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Run DSSP secondary structure for a single replicate.

        This thin wrapper preserves the legacy cache filename while delegating
        trajectory analysis to the runner-backed base implementation.

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
        from polyzymd.analyses.secondary_structure._results import SecondaryStructureResult

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        output_dir = ctx.output_dir
        eq_str = f"eq{eq_value:g}{eq_unit}"
        result_file = output_dir / f"secondary_structure_{eq_str}.json"

        cached = self._check_cache(
            SecondaryStructureResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=ctx.sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

        logger.info(f"Computing secondary structure for replicate {replicate}")
        return super().run_replicate(ctx, replicate)

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return the secondary-structure loader class for the runner seam.

        Returns
        -------
        type[Any]
            Loader class patched by secondary-structure unit tests.
        """

        return TrajectoryLoader

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve a validated frame window and attach mdtraj paths.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        loader : Any
            Trajectory loader for this replicate.
        universe : Any
            Loaded universe used for frame-count metadata.

        Returns
        -------
        _SecondaryStructureTrajectoryWindow
            Shared frame window plus trajectory files and topology path.
        """

        from polyzymd.analyses.shared.window import resolve_trajectory_window

        min_frames = 10
        n_frames_total = len(universe.trajectory)
        timestep_ps = float(loader.get_timestep(replicate, unit="ps"))
        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")
        start_frame = int(eq_time_ps / timestep_ps) if timestep_ps > 0 else 0
        if start_frame >= n_frames_total - min_frames:
            raise ValueError(
                f"Equilibration of {eq_time_ps / 1000:.1f} ns "
                f"({start_frame} frames) leaves fewer than {min_frames} "
                f"production frames from {n_frames_total} total. "
                "Reduce equilibration time or use a longer trajectory."
            )

        base_window = resolve_trajectory_window(
            equilibration=ctx.equilibration,
            n_frames_total=n_frames_total,
            timestep_ps=timestep_ps,
            min_frames=min_frames,
        )
        traj_info = loader.get_trajectory_info(replicate)
        return _SecondaryStructureTrajectoryWindow(
            base_window=base_window,
            trajectory_files=tuple(Path(path) for path in traj_info.trajectory_files),
            topology_file=Path(traj_info.topology_file),
        )

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any:
        """Build the runner-backed secondary-structure execution object.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        universe : Any
            Loaded universe for the replicate.
        window : Any
            Resolved trajectory window with mdtraj file metadata.

        Returns
        -------
        Any
            Runner object compatible with the shared trajectory seam.
        """

        del replicate, universe
        return SecondaryStructureReplicateRunner(
            trajectory_files=window.trajectory_files,
            topology_file=window.topology_file,
            chain_id=ctx.settings.chain_id,
            chain_index_func=_chain_letter_to_index,
            encode_dssp_func=_encode_dssp_matrix,
        )

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Serialize runner output into the legacy result schema.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        runner : Any
            Executed secondary-structure runner.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        SecondaryStructureResult
            Cache-compatible per-replicate secondary-structure result.
        """

        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.secondary_structure._results import SecondaryStructureResult

        del window
        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        payload = runner.results
        if payload is None:
            raise ValueError("Secondary-structure runner did not produce results")

        result = SecondaryStructureResult(
            config_hash=compute_config_hash(ctx.sim_config),
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string=payload.selection_string,
            n_frames=payload.n_frames,
            n_residues=payload.n_residues,
            residue_ids=payload.residue_ids,
            residue_names=payload.residue_names,
            persistence_coil=payload.persistence_coil,
            persistence_helix=payload.persistence_helix,
            persistence_strand=payload.persistence_strand,
            overall_helix_fraction=payload.overall_helix_fraction,
            overall_strand_fraction=payload.overall_strand_fraction,
            overall_coil_fraction=payload.overall_coil_fraction,
            settings_fingerprint=settings_fingerprint(ctx.settings),
        )

        output_dir = Path(ctx.output_dir)
        eq_str = f"eq{eq_value:g}{eq_unit}"
        result_prefix = f"secondary_structure_{eq_str}"
        json_path = result.save_with_matrix(
            directory=output_dir,
            matrix=payload.ss_matrix,
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

        if n_reps == 1:
            sem_coil = np.zeros(coil_stack.shape[1], dtype=float).tolist()
            sem_helix = np.zeros(helix_stack.shape[1], dtype=float).tolist()
            sem_strand = np.zeros(strand_stack.shape[1], dtype=float).tolist()
        else:
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

        if n_reps == 1:
            sem_overall_helix = 0.0
            sem_overall_strand = 0.0
            sem_overall_coil = 0.0
        else:
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
            settings_fingerprint=settings_fingerprint(ctx.settings),
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
        eq_str = f"eq{first_result.equilibration_time:g}{first_result.equilibration_unit}"
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
