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

from polyzymd.analyses.secondary_structure._results import (
    SecondaryStructureAggregatedResult,
)
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.shared.config_hash import compute_config_hash, validate_config_hash
from polyzymd.analyses.shared.loader import TrajectoryLoader, convert_time, parse_time_string
from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    find_json,
    get_output_path,
    get_theme,
    grouped_bars,
    save_figure,
    symmetric_clim,
)

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

        chain_id = getattr(settings, "chain_id", "A")

        # Parse equilibration time
        eq_value, eq_unit = parse_time_string(ctx.equilibration)

        # Initialize loader and config hash
        loader = TrajectoryLoader(sim_config)
        config_hash = compute_config_hash(sim_config)

        # Determine output path and check cache
        output_dir = ctx.output_dir
        eq_str = f"eq{eq_value:.0f}{eq_unit}"
        result_prefix = f"secondary_structure_{eq_str}"
        result_file = output_dir / f"{result_prefix}.json"

        if not ctx.recompute and result_file.exists():
            logger.info(f"Loading cached SS result from {result_file}")
            result = SecondaryStructureResult.load(result_file)
            validate_config_hash(result.config_hash, sim_config)
            return result

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
            # Fallback: try selecting all protein atoms
            protein_indices = traj.topology.select("protein")
            logger.warning(
                f"Chain '{chain_id}' not found; falling back to "
                f"'protein' selection ({len(protein_indices)} atoms)"
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

        # Build the data dict expected by the plotting helpers
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

        data["__meta__"] = {"results_dir": ctx.results_dir}

        if not labels:
            return plots

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        plot_settings = ctx.plot_settings
        if plot_settings is None:
            from polyzymd.compare.config import PlotSettings

            plot_settings = PlotSettings()

        try:
            result = _plot_ss_timeline_heatmap(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"SS timeline heatmap plot failed: {exc}")

        try:
            result = _plot_ss_content_bars(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"SS content bars plot failed: {exc}")

        try:
            result = _plot_ss_individual_bars(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"SS individual bars plot failed: {exc}")

        try:
            result = _plot_ss_persistence_diff_heatmap(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"SS persistence diff heatmap plot failed: {exc}")

        return plots

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Backward-compatible filename helper retained for tests."""
        eq_str = f"eq{first_result.equilibration_time:.0f}{first_result.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
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


# SS integer encoding colors (0=coil/grey, 1=helix/red, 2=strand/blue)
_SS_COLORS = {0: "#CCCCCC", 1: "#E74C3C", 2: "#3498DB"}
_SS_NAMES = {0: "No SS", 1: "Helix", 2: "\u03b2-Sheet"}


# Per-SS-type metadata: (internal_key, display_name, bar_color)
_SS_INDIVIDUAL_SPECS: list[tuple[str, str, str]] = [
    ("helix", "Helix", "#E74C3C"),
    ("strand", "\u03b2-Sheet", "#3498DB"),
    ("coil", "No SS", "#95A5A6"),
]


def _find_ss_comparison_result(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> Any | None:
    """Try to locate a saved SSComparisonResult JSON."""
    from polyzymd.compare.io.results import find_comparison_result
    from polyzymd.compare.results.secondary_structure import SSComparisonResult

    return find_comparison_result(
        data,
        labels,
        glob_patterns=["secondary_structure_comparison*.json"],
        loader=SSComparisonResult.load,
        analysis_type="secondary_structure",
        fallback_filenames=["secondary_structure_comparison.json"],
        log=log,
    )


def _plot_ss_timeline_heatmap(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate SS timeline heatmaps for each condition."""
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    t = get_theme(plot_settings)
    generated: list[Path] = []

    cmap = mcolors.ListedColormap([_SS_COLORS[0], _SS_COLORS[1], _SS_COLORS[2]])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        matrix, residue_ids = _load_ss_timeline_matrix(cond_data)
        if matrix is None:
            logger.debug(f"No SS timeline data for {label}")
            continue

        n_frames, n_residues = matrix.shape
        fig_width = max(14, n_residues * 0.05)
        fig_height = max(4, n_frames * 0.008)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        ax.imshow(
            matrix.T,
            aspect="auto",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
            origin="lower",
        )

        tick_stride = max(1, n_residues // 30)
        yticks = list(range(0, n_residues, tick_stride))
        yticklabels = [str(residue_ids[i]) for i in yticks]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=t.tiny_fontsize)

        apply_axis_style(
            ax,
            plot_settings,
            title=f"Secondary Structure Timeline — {label}",
            xlabel="Frame",
            ylabel="Residue",
        )

        legend_patches = [Patch(facecolor=_SS_COLORS[i], label=_SS_NAMES[i]) for i in [1, 2, 0]]
        apply_legend(
            ax,
            plot_settings,
            loc="upper right",
            bbox_to_anchor=None,
            fontsize=t.small_fontsize,
            handles=legend_patches,
            framealpha=0.8,
        )

        plt.tight_layout()

        from polyzymd.compare.io.paths import sanitize_label

        safe_label = sanitize_label(label)
        output_path = get_output_path(output_dir, f"ss_timeline_{safe_label}", plot_settings)
        generated.append(save_figure(fig, output_path, plot_settings))

    return generated


def _load_ss_timeline_matrix(cond_data: dict[str, Any]) -> tuple[np.ndarray | None, list[int]]:
    """Load per-frame SS matrix from replicate NPZ data."""
    import json as json_mod

    analysis_dir = cond_data.get("analysis_dir")
    if analysis_dir is None:
        return None, []

    analysis_dir = Path(analysis_dir)
    replicates = cond_data.get("replicates", [1])

    for rep in replicates:
        rep_dir = analysis_dir / f"run_{rep}"
        if not rep_dir.is_dir():
            continue

        npz_files = sorted(rep_dir.glob("*_matrix.npz"))
        if not npz_files:
            continue

        json_files = sorted(rep_dir.glob("secondary_structure*.json"))
        if not json_files:
            with np.load(str(npz_files[0])) as npz_data:
                matrix = np.asarray(npz_data["ss_matrix"])
            residue_ids = list(range(1, matrix.shape[1] + 1))
            return matrix, residue_ids

        try:
            with json_files[0].open() as handle:
                result_data = json_mod.load(handle)
            residue_ids = result_data.get("residue_ids", [])
        except Exception:
            residue_ids = []

        try:
            with np.load(str(npz_files[0])) as npz_data:
                matrix = np.asarray(npz_data["ss_matrix"])
            if not residue_ids:
                residue_ids = list(range(1, matrix.shape[1] + 1))
            return matrix, residue_ids
        except Exception as exc:
            logger.debug(f"Failed to load NPZ from {npz_files[0]}: {exc}")

    return None, []


def _plot_ss_content_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped SS-content bars from comparison or aggregated results."""
    import json as json_mod

    import matplotlib.pyplot as plt

    t = get_theme(plot_settings)
    result = _find_ss_comparison_result(data, labels)
    if result is not None:
        conditions = result.conditions
        n = len(conditions)

        cond_labels = [c.label for c in conditions]
        helix_means = [c.mean_helix for c in conditions]
        helix_sems = [c.sem_helix for c in conditions]
        strand_means = [c.mean_strand for c in conditions]
        strand_sems = [c.sem_strand for c in conditions]
        coil_means = [c.mean_coil for c in conditions]
        coil_sems = [c.sem_coil for c in conditions]

        helix_reps = [c.per_replicate_helix for c in conditions]
        strand_reps = [c.per_replicate_strand for c in conditions]
        coil_reps = [c.per_replicate_coil for c in conditions]

        x = np.arange(n)
        series = [
            ("Helix", helix_means, helix_sems),
            ("\u03b2-Sheet", strand_means, strand_sems),
            ("No SS", coil_means, coil_sems),
        ]
        ss_bar_colors = ["#E74C3C", "#3498DB", "#95A5A6"]
        replicate_values = [
            [helix_reps[i] for i in range(n)],
            [strand_reps[i] for i in range(n)],
            [coil_reps[i] for i in range(n)],
        ]

        fig, ax = plt.subplots(figsize=(max(8, n * 1.6), 5))
        grouped_bars(
            ax,
            x,
            series,
            ss_bar_colors,
            plot_settings,
            reference_line=None,
            replicate_values=replicate_values,
        )
        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
        apply_legend(ax, plot_settings)
        ax.set_ylim(bottom=0)
        apply_axis_style(
            ax,
            plot_settings,
            title="Secondary Structure Content Comparison",
            xlabel=None,
            ylabel="Fraction of (residue, frame) entries",
        )
        plt.tight_layout()

        output_path = get_output_path(output_dir, "ss_content_bars", plot_settings)
        return [save_figure(fig, output_path, plot_settings)]

    cond_labels: list[str] = []
    helix_means: list[float] = []
    helix_sems: list[float] = []
    strand_means: list[float] = []
    strand_sems: list[float] = []
    coil_means: list[float] = []
    coil_sems: list[float] = []
    helix_reps: list[list[float]] = []
    strand_reps: list[list[float]] = []
    coil_reps: list[list[float]] = []

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue
        aggregated_dir = Path(aggregated_dir)

        result_file = find_json(aggregated_dir, "secondary_structure_aggregated.json")
        if result_file is None:
            result_file = find_json(aggregated_dir, "secondary_structure.json")
        if result_file is None:
            continue

        try:
            with result_file.open() as handle:
                agg = json_mod.load(handle)

            cond_labels.append(label)
            helix_means.append(agg.get("mean_overall_helix", 0.0))
            helix_sems.append(agg.get("sem_overall_helix", 0.0))
            strand_means.append(agg.get("mean_overall_strand", 0.0))
            strand_sems.append(agg.get("sem_overall_strand", 0.0))
            coil_means.append(agg.get("mean_overall_coil", 0.0))
            coil_sems.append(agg.get("sem_overall_coil", 0.0))
            helix_reps.append(agg.get("per_replicate_helix", []))
            strand_reps.append(agg.get("per_replicate_strand", []))
            coil_reps.append(agg.get("per_replicate_coil", []))
        except Exception as exc:
            logger.warning(f"Failed to load aggregated SS for {label}: {exc}")

    if not cond_labels:
        logger.warning("No aggregated SS data found for content bars")
        return []

    n = len(cond_labels)
    x = np.arange(n)
    series = [
        ("Helix", helix_means, helix_sems),
        ("\u03b2-Sheet", strand_means, strand_sems),
        ("No SS", coil_means, coil_sems),
    ]
    ss_bar_colors = ["#E74C3C", "#3498DB", "#95A5A6"]
    replicate_values = [helix_reps, strand_reps, coil_reps]
    has_reps = any(any(r for r in reps) for reps in replicate_values)

    fig, ax = plt.subplots(figsize=(max(8, n * 1.6), 5))
    grouped_bars(
        ax,
        x,
        series,
        ss_bar_colors,
        plot_settings,
        reference_line=None,
        replicate_values=replicate_values if has_reps else None,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
    apply_legend(ax, plot_settings)
    ax.set_ylim(bottom=0)
    apply_axis_style(
        ax,
        plot_settings,
        title="Secondary Structure Content Comparison",
        xlabel=None,
        ylabel="Fraction of (residue, frame) entries",
    )

    plt.tight_layout()
    output_path = get_output_path(output_dir, "ss_content_bars", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


def _plot_ss_individual_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate per-SS-type bar charts."""
    import json as json_mod

    result = _find_ss_comparison_result(data, labels)
    if result is not None:
        conditions = result.conditions
        cond_labels = [c.label for c in conditions]
        ss_data = {
            "helix": {
                "means": [c.mean_helix for c in conditions],
                "sems": [c.sem_helix for c in conditions],
                "reps": [c.per_replicate_helix for c in conditions],
            },
            "strand": {
                "means": [c.mean_strand for c in conditions],
                "sems": [c.sem_strand for c in conditions],
                "reps": [c.per_replicate_strand for c in conditions],
            },
            "coil": {
                "means": [c.mean_coil for c in conditions],
                "sems": [c.sem_coil for c in conditions],
                "reps": [c.per_replicate_coil for c in conditions],
            },
        }
        return _render_ss_individual_plots(
            cond_labels,
            len(cond_labels),
            ss_data,
            output_dir,
            plot_settings,
            has_reps=True,
        )

    cond_labels: list[str] = []
    ss_data: dict[str, dict[str, list[Any]]] = {
        "helix": {"means": [], "sems": [], "reps": []},
        "strand": {"means": [], "sems": [], "reps": []},
        "coil": {"means": [], "sems": [], "reps": []},
    }

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue
        aggregated_dir = Path(aggregated_dir)

        result_file = find_json(aggregated_dir, "secondary_structure_aggregated.json")
        if result_file is None:
            result_file = find_json(aggregated_dir, "secondary_structure.json")
        if result_file is None:
            continue

        try:
            with result_file.open() as handle:
                agg = json_mod.load(handle)

            cond_labels.append(label)
            ss_data["helix"]["means"].append(agg.get("mean_overall_helix", 0.0))
            ss_data["helix"]["sems"].append(agg.get("sem_overall_helix", 0.0))
            ss_data["helix"]["reps"].append(agg.get("per_replicate_helix", []))
            ss_data["strand"]["means"].append(agg.get("mean_overall_strand", 0.0))
            ss_data["strand"]["sems"].append(agg.get("sem_overall_strand", 0.0))
            ss_data["strand"]["reps"].append(agg.get("per_replicate_strand", []))
            ss_data["coil"]["means"].append(agg.get("mean_overall_coil", 0.0))
            ss_data["coil"]["sems"].append(agg.get("sem_overall_coil", 0.0))
            ss_data["coil"]["reps"].append(agg.get("per_replicate_coil", []))
        except Exception as exc:
            logger.warning(f"Failed to load aggregated SS for {label}: {exc}")

    if not cond_labels:
        logger.warning("No aggregated SS data found for individual bars")
        return []

    has_reps = any(any(r for r in ss_data[key]["reps"]) for key in ("helix", "strand", "coil"))
    return _render_ss_individual_plots(
        cond_labels,
        len(cond_labels),
        ss_data,
        output_dir,
        plot_settings,
        has_reps=has_reps,
    )


def _render_ss_individual_plots(
    cond_labels: list[str],
    n: int,
    ss_data: dict[str, dict[str, list[Any]]],
    output_dir: Path,
    plot_settings: Any,
    *,
    has_reps: bool,
) -> list[Path]:
    """Render one bar chart per SS type."""
    import matplotlib.pyplot as plt

    t = get_theme(plot_settings)
    generated: list[Path] = []
    x = np.arange(n)

    tab10 = plt.cm.get_cmap("tab10")
    condition_colors = [tab10(i % 10) for i in range(n)]

    for internal_key, display_name, _bar_color in _SS_INDIVIDUAL_SPECS:
        means = ss_data[internal_key]["means"]
        sems = ss_data[internal_key]["sems"]

        fig, ax = plt.subplots(figsize=(max(8, n * 1.2), 5))
        ax.bar(
            x,
            means,
            yerr=sems,
            color=condition_colors,
            alpha=t.bar_alpha,
            edgecolor=t.bar_edgecolor,
            linewidth=t.bar_linewidth,
            capsize=t.bar_capsize,
        )

        if has_reps and "reps" in ss_data[internal_key]:
            rng = np.random.default_rng(seed=42)
            reps = ss_data[internal_key]["reps"]
            for j in range(n):
                if j < len(reps) and reps[j]:
                    rep_vals = np.asarray(reps[j], dtype=float)
                    jitter = rng.uniform(-0.2, 0.2, size=len(rep_vals))
                    ax.scatter(
                        np.full_like(rep_vals, float(j)) + jitter,
                        rep_vals,
                        color=t.dot_color,
                        s=t.dot_size,
                        zorder=5,
                        alpha=t.dot_alpha,
                        edgecolors="none",
                    )

        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
        ax.set_ylim(bottom=0)
        apply_axis_style(
            ax,
            plot_settings,
            title=f"{display_name} Content by Condition",
            xlabel=None,
            ylabel="Fraction of (residue, frame) entries",
        )

        plt.tight_layout()
        output_path = get_output_path(output_dir, f"ss_{internal_key}_bars", plot_settings)
        generated.append(save_figure(fig, output_path, plot_settings))

    return generated


def _plot_ss_persistence_diff_heatmap(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate condition x residue Delta(helix persistence) heatmap."""
    import json as json_mod

    import matplotlib.pyplot as plt

    t = get_theme(plot_settings)
    persistence_data: dict[str, dict[str, Any]] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue
        aggregated_dir = Path(aggregated_dir)

        result_file = find_json(aggregated_dir, "secondary_structure_aggregated.json")
        if result_file is None:
            result_file = find_json(aggregated_dir, "secondary_structure.json")
        if result_file is None:
            continue

        try:
            with result_file.open() as handle:
                agg = json_mod.load(handle)

            helix_persist = agg.get("mean_persistence_helix")
            residue_ids = agg.get("residue_ids")
            if helix_persist is not None and residue_ids is not None:
                persistence_data[label] = {
                    "helix": np.array(helix_persist),
                    "residue_ids": residue_ids,
                }
        except Exception as exc:
            logger.warning(f"Failed to load SS persistence for {label}: {exc}")

    if len(persistence_data) < 2:
        logger.warning("Need at least 2 conditions for persistence difference heatmap")
        return []

    meta = data.get("__meta__", {})
    source_path = meta.get("comparison_source_path")
    control_label: str | None = None
    if source_path is not None:
        try:
            from polyzymd.compare.config import ComparisonConfig

            comp_config = ComparisonConfig.from_yaml(source_path)
            control_label = comp_config.control_label
        except Exception:
            pass

    available_labels = [lbl for lbl in labels if lbl in persistence_data]
    if not available_labels:
        return []

    if control_label is None or control_label not in persistence_data:
        control_label = available_labels[0]

    control_helix = persistence_data[control_label]["helix"]
    residue_ids = persistence_data[control_label]["residue_ids"]
    n_residues = len(residue_ids)

    diff_labels: list[str] = []
    diff_rows: list[np.ndarray] = []
    for label in available_labels:
        if label == control_label:
            continue
        cond_helix = persistence_data[label]["helix"]
        if len(cond_helix) != n_residues:
            logger.warning(
                f"Residue count mismatch for {label}: {len(cond_helix)} vs {n_residues} (control)"
            )
            continue
        diff_rows.append(cond_helix - control_helix)
        diff_labels.append(label)

    if not diff_rows:
        logger.warning("No valid conditions for persistence diff heatmap")
        return []

    diff_matrix = np.array(diff_rows)
    n_conds = len(diff_labels)
    fig_width = max(14, n_residues * 0.05)
    fig_height = max(3, n_conds * 0.6 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    vmin, vmax = symmetric_clim(diff_matrix.ravel(), pad=0.01)
    im = ax.imshow(
        diff_matrix,
        aspect="auto",
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )

    ax.set_yticks(range(n_conds))
    ax.set_yticklabels(diff_labels, fontsize=t.tick_fontsize)

    tick_stride = max(1, n_residues // 40)
    xticks = list(range(0, n_residues, tick_stride))
    xticklabels = [str(residue_ids[i]) for i in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=t.tiny_fontsize, rotation=90)

    apply_axis_style(
        ax,
        plot_settings,
        title=f"Per-Residue Helix Persistence Change vs {control_label}",
        xlabel="Residue",
        ylabel=None,
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label(
        r"$\Delta$ Helix Persistence (condition $-$ control)",
        fontsize=t.tick_fontsize,
    )

    plt.tight_layout()
    output_path = get_output_path(output_dir, "ss_persistence_diff_heatmap", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]
