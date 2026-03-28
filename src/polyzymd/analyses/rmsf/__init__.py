"""RMSF analysis plugin.

Computes per-residue Root Mean Square Fluctuation from MD trajectories,
aggregates across replicates, compares conditions via the default scalar
comparison pipeline, and produces comparison/profile plots.

All heavy computation is self-contained — no delegation to legacy
``analysis.rmsf`` calculator classes.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory
from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    get_colors,
    get_output_path,
    save_figure,
)
from polyzymd.analyses.shared.autocorrelation import (
    compute_acf,
    estimate_correlation_time,
    get_independent_indices,
)
from polyzymd.analyses.shared.centroid import ReferenceMode
from polyzymd.analyses.shared.config_hash import compute_config_hash, validate_config_hash
from polyzymd.analyses.shared.diagnostics import (
    get_selection_diagnostics,
    validate_equilibration_time,
)
from polyzymd.analyses.shared.loader import (
    TrajectoryLoader,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analyses.shared.statistics import aggregate_per_residue_stats, compute_sem

from polyzymd.analyses.rmsf._results import RMSFAggregatedResult

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from polyzymd.analyses.rmsf._results import RMSFResult

logger = logging.getLogger("polyzymd.analyses.rmsf")


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class RMSFSettings(BaseModel):
    """Settings for RMSF analysis.

    Attributes
    ----------
    selection : str
        MDAnalysis selection string for RMSF calculation.
    reference_mode : str
        Reference structure mode: centroid, average, frame, or external.
    reference_frame : int | None
        Frame number if reference_mode is 'frame' (1-indexed).
    reference_file : str | None
        Path to external PDB file if reference_mode is 'external'.
    alignment_selection : str
        MDAnalysis selection for trajectory alignment.
    centroid_selection : str
        MDAnalysis selection for centroid finding (K-Means clustering).
    """

    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection string for RMSF calculation",
    )
    reference_mode: str = Field(
        default="centroid",
        description="Reference structure mode: centroid, average, frame, or external",
    )
    reference_frame: int | None = Field(
        default=None,
        description="Frame number if reference_mode is 'frame' (1-indexed)",
    )
    reference_file: str | None = Field(
        default=None,
        description="Path to external PDB file if reference_mode is 'external'",
    )
    alignment_selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for trajectory alignment",
    )
    centroid_selection: str = Field(
        default="protein",
        description="MDAnalysis selection for centroid finding (K-Means clustering)",
    )

    @field_validator("reference_mode", mode="after")
    @classmethod
    def validate_reference_mode(cls, v: str) -> str:
        valid = {"centroid", "average", "frame", "external"}
        if v not in valid:
            raise ValueError(f"reference_mode must be one of {valid}, got {v!r}")
        return v


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class RMSFAnalysis(Analysis):
    """RMSF analysis: per-residue flexibility from MD trajectories.

    Performs the complete RMSF workflow inline:

    1. Load trajectories from config
    2. Apply equilibration offset
    3. Find reference frame and align trajectory
    4. Compute autocorrelation and select independent frames
    5. Calculate per-residue RMSF
    6. Aggregate across replicates with SEM

    The ``compare()`` method is NOT overridden — it uses the default
    implementation which calls ``extract_metrics()`` to get ``mean_rmsf``
    as a single scalar metric with ``higher_is_better=False`` (lower RMSF
    = more stable = better ranking).

    Plots
    -----
    - ``rmsf_comparison.png``: Horizontal bar chart of mean RMSF per condition.
    - ``rmsf_profile.png``: Per-residue RMSF overlay with optional SS annotation.
    """

    name: ClassVar[str] = "rmsf"
    Settings: ClassVar[type] = RMSFSettings
    AggregatedResultClass: ClassVar[type] = RMSFAggregatedResult
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute RMSF for a single replicate.

        Performs trajectory loading, alignment, autocorrelation-based
        subsampling, and per-residue RMSF calculation inline.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        RMSFResult
            Per-replicate RMSF result.
        """
        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rmsf._results import RMSFResult

        settings = ctx.settings
        sim_config = ctx.sim_config

        selection = getattr(settings, "selection", "protein and name CA")
        reference_mode: ReferenceMode = getattr(settings, "reference_mode", "centroid")
        reference_frame = getattr(settings, "reference_frame", None)
        reference_file = getattr(settings, "reference_file", None)
        alignment_selection = getattr(settings, "alignment_selection", "protein and name CA")
        centroid_selection = getattr(settings, "centroid_selection", "protein")

        # Validate reference_mode + reference_frame/file combinations
        if reference_mode == "frame" and reference_frame is None:
            raise ValueError("reference_frame is required when reference_mode='frame'")
        if reference_mode == "external" and reference_file is None:
            raise ValueError(
                "reference_file is required when reference_mode='external'. "
                "Provide a path to the external PDB reference structure."
            )
        if reference_mode == "external" and reference_file is not None:
            ref_path = Path(reference_file)
            if not ref_path.exists():
                raise ValueError(
                    f"reference_file does not exist: {ref_path}. "
                    "Provide a valid path to the external PDB reference structure."
                )

        # Parse equilibration time
        eq_value, eq_unit = parse_time_string(ctx.equilibration)

        # Initialize loader and config hash
        loader = TrajectoryLoader(sim_config)
        config_hash = compute_config_hash(sim_config)

        # Determine output path and check cache
        output_dir = ctx.output_dir
        eq_str = f"eq{eq_value:.0f}{eq_unit}"
        result_filename = f"rmsf_{eq_str}.json"
        result_file = output_dir / result_filename

        if not ctx.recompute and result_file.exists():
            logger.info(f"Loading cached result from {result_file}")
            result = RMSFResult.load(result_file)
            validate_config_hash(result.config_hash, sim_config)
            return result

        logger.info(f"Computing RMSF for replicate {replicate}")

        # Load universe
        u = loader.load_universe(replicate)
        traj_info = loader.get_trajectory_info(replicate)

        # Get atom selection for RMSF
        atoms = u.select_atoms(selection)
        if len(atoms) == 0:
            diag = get_selection_diagnostics(u, selection)
            raise ValueError(f"Selection '{selection}' matched no atoms.\n\n{diag}")

        logger.info(f"Selected {len(atoms)} atoms with '{selection}'")

        # Get timestep
        timestep = loader.get_timestep(replicate, unit="ps")

        # Determine start frame after equilibration
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")
        start_frame = time_to_frame(eq_time_ps, "ps", timestep, "ps")

        n_frames_total = len(u.trajectory)
        n_frames_after_eq = n_frames_total - start_frame

        # Validate equilibration time against trajectory length
        eq_time_ns = convert_time(eq_value, eq_unit, "ns")
        traj_time_ns = (n_frames_total * timestep) / 1000.0
        is_valid, eq_message = validate_equilibration_time(eq_time_ns, traj_time_ns)
        if not is_valid:
            raise ValueError(eq_message)
        if eq_message:
            logger.warning(eq_message)

        logger.info(
            f"Trajectory: {n_frames_total} frames, skipping first {start_frame} for equilibration"
        )

        # ===== ALIGNMENT STEP =====
        alignment_config = AlignmentConfig(
            enabled=True,
            reference_mode=reference_mode,
            reference_frame=reference_frame,
            selection=alignment_selection,
            centroid_selection=centroid_selection,
            reference_file=(Path(reference_file) if reference_file is not None else None),
        )
        ref_frame_idx = align_trajectory(
            u, alignment_config, start_frame=start_frame, stop_frame=n_frames_total
        )
        ref_frame_1indexed = ref_frame_idx + 1 if ref_frame_idx is not None else None

        logger.info(
            f"Alignment: mode='{reference_mode}', "
            f"reference_frame={ref_frame_1indexed}, "
            f"selection='{alignment_selection}'"
        )

        # ===== AUTOCORRELATION & FRAME SELECTION =====
        correlation_time: float | None = None
        correlation_time_unit: str | None = None
        n_independent: int | None = None
        frame_indices: NDArray[np.int64]

        if n_frames_after_eq > 100:
            rmsd_timeseries = _compute_rmsd_timeseries(u, atoms, start_frame)
            acf_result = compute_acf(rmsd_timeseries, timestep=timestep, timestep_unit="ps")
            tau_result = estimate_correlation_time(acf_result, n_frames=n_frames_after_eq)

            correlation_time = tau_result.tau
            correlation_time_unit = tau_result.tau_unit
            n_independent = tau_result.n_independent

            logger.info(
                f"Correlation time: {correlation_time:.2f} {correlation_time_unit}, "
                f"~{n_independent} independent frames"
            )

            frame_indices = get_independent_indices(
                n_frames=n_frames_total,
                correlation_time=correlation_time,
                timestep=timestep,
                start_frame=start_frame,
            )
        else:
            frame_indices = np.arange(start_frame, n_frames_total, dtype=np.int64)

        n_frames_used = len(frame_indices)
        logger.info(f"Using {n_frames_used} frames for RMSF calculation")

        # ===== RMSF CALCULATION =====
        # Load external reference positions if needed
        external_ref_positions: NDArray[np.float64] | None = None
        if reference_mode == "external" and reference_file is not None:
            import MDAnalysis as mda_ext

            ref_path = Path(reference_file)
            logger.info(f"Loading external reference positions from: {ref_path}")
            ref_universe = mda_ext.Universe(str(ref_path))
            ref_atoms = ref_universe.select_atoms(selection)

            if len(ref_atoms) != len(atoms):
                raise ValueError(
                    f"External PDB atom count ({len(ref_atoms)}) does not match "
                    f"trajectory selection ({len(atoms)}) for '{selection}'. "
                    f"Cannot use external PDB positions as RMSF reference."
                )
            external_ref_positions = ref_atoms.positions.copy().astype(np.float64)
            logger.info(
                f"Using external PDB positions as RMSF reference "
                f"({len(ref_atoms)} atoms from '{selection}')"
            )

        rmsf_values = _compute_rmsf(u, atoms, frame_indices, external_ref_positions)

        # Get residue information
        residue_ids = [int(r.resid) for r in atoms.residues]
        residue_names = [r.resname for r in atoms.residues]

        if "NAME CA" in selection.upper():
            per_residue_rmsf = rmsf_values
        else:
            per_residue_rmsf = _aggregate_per_residue(atoms, rmsf_values)
            unique_residues = atoms.residues
            residue_ids = [int(r.resid) for r in unique_residues]
            residue_names = [r.resname for r in unique_residues]

        # Summary statistics
        mean_rmsf = float(np.mean(per_residue_rmsf))
        std_rmsf = float(np.std(per_residue_rmsf))
        min_rmsf = float(np.min(per_residue_rmsf))
        max_rmsf = float(np.max(per_residue_rmsf))

        result = RMSFResult(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string=selection,
            correlation_time=correlation_time,
            correlation_time_unit=correlation_time_unit,
            n_independent_frames=n_independent,
            residue_ids=residue_ids,
            residue_names=residue_names,
            rmsf_values=per_residue_rmsf.tolist(),
            mean_rmsf=mean_rmsf,
            std_rmsf=std_rmsf,
            min_rmsf=min_rmsf,
            max_rmsf=max_rmsf,
            reference_mode=reference_mode,
            reference_frame=ref_frame_1indexed,
            alignment_selection=alignment_selection,
            reference_file=(str(reference_file) if reference_file is not None else None),
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
            trajectory_files=[str(f) for f in traj_info.trajectory_files],
        )

        result.save(result_file)
        logger.info(f"Saved result to {result_file}")

        return result

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate RMSF across replicates for one condition.

        Computes per-residue mean +/- SEM and overall statistics from
        the already-computed per-replicate results.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[RMSFResult]
            Per-replicate RMSF results.

        Returns
        -------
        RMSFAggregatedResult
            Aggregated result with per-residue and overall statistics.
        """
        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.rmsf._results import RMSFAggregatedResult

        settings = ctx.settings

        # Collect per-residue RMSF arrays from each replicate
        per_replicate_rmsf = [np.array(r.rmsf_values) for r in results]

        # Aggregate per-residue statistics
        per_residue_stats = aggregate_per_residue_stats(
            per_replicate_rmsf,
            residue_ids=np.array(results[0].residue_ids),
        )

        # Aggregate whole-protein statistics
        per_replicate_means = [r.mean_rmsf for r in results]
        overall_stats = compute_sem(per_replicate_means)

        config_hash = results[0].config_hash

        agg_result = RMSFAggregatedResult(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=results[0].equilibration_time,
            equilibration_unit=results[0].equilibration_unit,
            selection_string=getattr(settings, "selection", results[0].selection_string),
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            residue_ids=results[0].residue_ids,
            residue_names=results[0].residue_names,
            mean_rmsf_per_residue=per_residue_stats.means.tolist(),
            sem_rmsf_per_residue=per_residue_stats.sems.tolist(),
            per_replicate_mean_rmsf=per_replicate_means,
            overall_mean_rmsf=overall_stats.mean,
            overall_sem_rmsf=overall_stats.sem,
            overall_min_rmsf=float(np.min(per_residue_stats.means)),
            overall_max_rmsf=float(np.max(per_residue_stats.means)),
            source_result_files=[],
        )

        target_path = ctx.result_path
        if target_path is None:
            filename = self._make_aggregated_filename(ctx.replicates, results[0])
            target_path = ctx.output_dir / filename
        self.save_result(agg_result, target_path)
        logger.info(f"Saved aggregated RMSF to {target_path}")

        return agg_result

    # === Optional methods ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract the mean RMSF metric for scalar comparison.

        Parameters
        ----------
        summary : RMSFAggregatedResult
            Aggregated RMSF result.

        Returns
        -------
        dict[str, MetricValue]
            Single entry ``"mean_rmsf"`` with ``higher_is_better=False``
            (lower RMSF = more stable).
        """
        return {
            "mean_rmsf": MetricValue(
                name="mean_rmsf",
                mean=summary.overall_mean_rmsf,
                sem=summary.overall_sem_rmsf,
                replicate_values=summary.per_replicate_mean_rmsf,
                higher_is_better=False,
                direction_labels=("stabilizing", "unchanged", "destabilizing"),
            ),
        }

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format RMSF comparison result for CLI display.

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
                title="RMSF Comparison",
                metric_label="Mean RMSF",
                metric_unit="A",
                metric_key="mean_rmsf",
                output_format=output_format,
                higher_is_better=False,
            )
        return super().format(result, output_format)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate RMSF comparison and profile plots.

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

        # Build the data dict for plotting helpers
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
            result = _plot_rmsf_comparison(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"RMSF comparison plot failed: {exc}")

        try:
            result = _plot_rmsf_profile(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"RMSF profile plot failed: {exc}")

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
        return f"rmsf_{rep_str}_{eq_str}.json"


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Private plotting helpers (inlined from compare/plotters/rmsf.py)
# ---------------------------------------------------------------------------


def _plot_rmsf_comparison(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate RMSF comparison bar chart.

    Looks for a pre-computed comparison result JSON first. If not found,
    falls back to loading aggregated RMSF results.
    """
    comparison_result = _find_rmsf_comparison_result(data, labels)

    if comparison_result is not None:
        return _plot_rmsf_comparison_from_result(comparison_result, output_dir, plot_settings)
    else:
        return _plot_rmsf_comparison_from_aggregated(data, labels, output_dir, plot_settings)


def _find_rmsf_comparison_result(
    data: dict[str, Any],
    labels: Sequence[str],
) -> Any | None:
    """Try to find a pre-computed RMSF comparison result."""
    from polyzymd.compare.io.results import find_comparison_result
    from polyzymd.compare.results.rmsf import RMSFComparisonResult
    from polyzymd.compare.results.rmsf_legacy import ComparisonResult

    def _try_load(path: Path) -> Any | None:
        try:
            return RMSFComparisonResult.load(path)
        except Exception:
            try:
                return ComparisonResult.load(path)
            except Exception as e:
                logger.debug(f"Could not load {path}: {e}")
        return None

    return find_comparison_result(
        data,
        labels,
        glob_patterns=["rmsf_comparison*.json"],
        loader=_try_load,
        analysis_type="rmsf",
        fallback_filenames=["rmsf_comparison.json", "comparison_result.json"],
        log=logger,
    )


def _plot_rmsf_comparison_from_result(
    result: Any,
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate horizontal bar chart from comparison result."""
    import matplotlib.pyplot as plt
    import numpy as np

    t = plot_settings.theme

    # Get conditions sorted by RMSF (lowest first)
    labels_sorted = (
        result.ranking if hasattr(result, "ranking") else [c.label for c in result.conditions]
    )

    means = []
    sems = []
    replicate_data: list[list[float]] = []

    for label in labels_sorted:
        cond = result.get_condition(label)
        means.append(cond.mean_rmsf)
        sems.append(cond.sem_rmsf)
        replicate_data.append(getattr(cond, "replicate_values", None) or [])

    n = len(labels_sorted)
    means_arr = np.array(means)
    sems_arr = np.array(sems)
    positions = np.arange(n)
    colors = get_colors(n, plot_settings)

    fig, ax = plt.subplots(figsize=plot_settings.rmsf.figsize_comparison)

    bar_height = 0.7
    ax.barh(
        positions,
        means_arr,
        xerr=sems_arr,
        color=colors,
        edgecolor=t.bar_edgecolor,
        linewidth=t.bar_linewidth,
        capsize=t.bar_capsize,
        height=bar_height,
    )

    # Overlay jittered replicate dots
    rng = np.random.default_rng(seed=42)
    for i, rep_vals in enumerate(replicate_data):
        if rep_vals:
            rep_arr = np.asarray(rep_vals, dtype=float)
            jitter = rng.uniform(-bar_height * 0.25, bar_height * 0.25, size=len(rep_arr))
            ax.scatter(
                rep_arr,
                np.full_like(rep_arr, float(positions[i])) + jitter,
                color=t.dot_color,
                s=t.dot_size,
                zorder=5,
                alpha=t.dot_alpha,
                edgecolors="none",
            )

    ax.set_yticks(positions)
    ax.set_yticklabels(labels_sorted)
    apply_axis_style(ax, plot_settings, title="RMSF Comparison", xlabel="Mean RMSF (Å)")
    ax.invert_yaxis()

    plt.tight_layout()

    output_path = get_output_path(output_dir, "rmsf_comparison", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


def _plot_rmsf_comparison_from_aggregated(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate simple bar chart from aggregated RMSF data."""
    import json

    import matplotlib.pyplot as plt
    import numpy as np

    # Collect mean RMSF and SEM for each condition
    plot_labels = []
    means = []
    sems = []

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue

        aggregated_dir = Path(aggregated_dir)

        # Look for aggregated RMSF result
        result_file = aggregated_dir / "rmsf_aggregated.json"
        if not result_file.exists():
            json_files = list(aggregated_dir.glob("*.json"))
            if json_files:
                result_file = json_files[0]
            else:
                continue

        try:
            with open(result_file) as f:
                agg_data = json.load(f)

            # Support multiple key naming conventions
            mean_val = (
                agg_data.get("overall_mean_rmsf")
                or agg_data.get("overall_mean")
                or agg_data.get("mean_rmsf")
            )
            sem_val = (
                agg_data.get("overall_sem_rmsf")
                or agg_data.get("overall_sem")
                or agg_data.get("sem_rmsf", 0)
            )

            if mean_val is not None:
                plot_labels.append(label)
                means.append(mean_val)
                sems.append(sem_val)

        except Exception as e:
            logger.warning(f"Failed to load aggregated RMSF for {label}: {e}")

    if not plot_labels:
        logger.warning("No aggregated RMSF data found")
        return []

    # Create simple bar chart
    fig, ax = plt.subplots(figsize=plot_settings.rmsf.figsize_comparison)

    t = plot_settings.theme
    positions = np.arange(len(plot_labels))
    colors = get_colors(len(plot_labels), plot_settings)

    ax.barh(
        positions,
        means,
        xerr=sems,
        color=colors,
        edgecolor=t.bar_edgecolor,
        linewidth=t.bar_linewidth,
        capsize=t.bar_capsize,
        height=0.7,
    )

    ax.set_yticks(positions)
    ax.set_yticklabels(plot_labels)
    apply_axis_style(ax, plot_settings, title="RMSF Comparison", xlabel="Mean RMSF (Å)")
    ax.invert_yaxis()

    plt.tight_layout()

    output_path = get_output_path(output_dir, "rmsf_comparison", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


def _plot_rmsf_profile(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate per-residue RMSF profile plot with optional SS annotation."""
    import matplotlib.pyplot as plt
    import numpy as np

    t = plot_settings.theme
    colors = get_colors(len(labels), plot_settings)

    # Load per-residue RMSF data for each condition
    profiles: dict[str, dict] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue

        profile_data = _load_rmsf_profile(Path(aggregated_dir))
        if profile_data:
            profiles[label] = profile_data

    if not profiles:
        logger.warning("No per-residue RMSF data found for profile plot")
        return []

    # Try to load SS annotation from reference PDB
    ss_annotation = _load_reference_ss(data)

    # Create figure: 2-row if SS available, 1-row otherwise
    figsize = plot_settings.rmsf.figsize_profile
    if ss_annotation is not None:
        fig, (ax_rmsf, ax_ss) = plt.subplots(
            2,
            1,
            figsize=figsize,
            gridspec_kw={"height_ratios": [4, 1]},
            sharex=True,
        )
    else:
        fig, ax_rmsf = plt.subplots(figsize=figsize)
        ax_ss = None

    for idx, label in enumerate(labels):
        if label not in profiles:
            continue

        profile = profiles[label]
        residues = np.array(profile["residues"])
        rmsf = np.array(profile["rmsf"])

        color = colors[idx] if idx < len(colors) else f"C{idx}"

        if plot_settings.rmsf.show_error and "sem" in profile:
            sem = np.array(profile["sem"])
            ax_rmsf.fill_between(
                residues,
                rmsf - sem,
                rmsf + sem,
                alpha=t.fill_alpha,
                color=color,
            )

        ax_rmsf.plot(residues, rmsf, label=label, color=color, linewidth=1.5)

    # Highlight residues if configured
    for resid in plot_settings.rmsf.highlight_residues:
        ax_rmsf.axvline(
            resid, color="red", linestyle="--", alpha=t.highlight_line_alpha, linewidth=1
        )

    apply_axis_style(ax_rmsf, plot_settings, title="Per-Residue RMSF Comparison", ylabel="RMSF (Å)")
    apply_legend(ax_rmsf, plot_settings)

    # Draw SS annotation bar if available
    if ax_ss is not None and ss_annotation is not None:
        _draw_ss_bar(ax_ss, ss_annotation, plot_settings)
    else:
        ax_rmsf.set_xlabel("Residue Number", fontsize=t.label_fontsize)

    plt.tight_layout()

    output_path = get_output_path(output_dir, "rmsf_profile", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


def _load_reference_ss(data: dict[str, Any]) -> dict | None:
    """Load reference SS assignment from the crystal/input PDB.

    Reads ``plugins.rmsf.reference_file`` from the comparison config and
    runs mdtraj DSSP on it to get per-residue SS assignments.

    Returns
    -------
    dict or None
        ``{"residue_ids": [...], "ss_codes": [...]}`` where ss_codes
        are integers (0=coil, 1=helix, 2=strand), or None on failure.
    """
    meta = data.get("__meta__", {})
    source_path = meta.get("comparison_source_path")
    if source_path is None:
        return None

    try:
        from polyzymd.compare.config import ComparisonConfig

        comp_config = ComparisonConfig.from_yaml(source_path)
        rmsf_settings = comp_config.plugins.get("rmsf")
        if rmsf_settings is None:
            return None
        reference_file = getattr(rmsf_settings, "reference_file", None)
        if reference_file is None:
            return None
    except Exception as exc:
        logger.debug(f"Could not load comparison config for SS bar: {exc}")
        return None

    ref_path = Path(reference_file)
    if not ref_path.is_absolute():
        ref_path = Path(source_path).parent / ref_path
    if not ref_path.exists():
        logger.debug(f"Reference PDB not found: {ref_path}")
        return None

    try:
        import mdtraj as md

        traj = md.load(str(ref_path))

        # Select protein atoms only
        protein_indices = traj.topology.select("protein")
        if len(protein_indices) == 0:
            return None
        traj_protein = traj.atom_slice(protein_indices)

        dssp = md.compute_dssp(traj_protein, simplified=True)
        ss_string = dssp[0]  # Single frame -> 1D array of chars

        # Map chars to integers
        char_to_int = {"C": 0, "H": 1, "E": 2, "NA": 0}
        ss_codes = [char_to_int.get(c, 0) for c in ss_string]

        # Get residue IDs
        residue_ids = [r.resSeq for r in traj_protein.topology.residues]

        return {"residue_ids": residue_ids, "ss_codes": ss_codes}

    except ImportError:
        logger.debug("mdtraj not available; skipping SS annotation bar")
        return None
    except Exception as exc:
        logger.debug(f"Failed to compute reference SS: {exc}")
        return None


def _draw_ss_bar(ax: Any, ss_annotation: dict, plot_settings: Any) -> None:
    """Draw a colored SS annotation bar on the given axes."""
    import matplotlib.colors as mcolors
    import numpy as np
    from matplotlib.patches import Patch

    t = plot_settings.theme

    residue_ids = np.array(ss_annotation["residue_ids"])
    ss_codes = np.array(ss_annotation["ss_codes"])

    # SS colors: 0=coil(grey), 1=helix(red), 2=strand(blue)
    ss_colors = {0: "#CCCCCC", 1: "#E74C3C", 2: "#3498DB"}
    ss_names = {0: "No SS", 1: "Helix", 2: "\u03b2-Sheet"}

    cmap = mcolors.ListedColormap([ss_colors[0], ss_colors[1], ss_colors[2]])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Plot as a 1-row heatmap: reshape to (1, n_residues)
    ss_row = ss_codes.reshape(1, -1)

    ax.imshow(
        ss_row,
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        extent=[
            residue_ids[0] - 0.5,
            residue_ids[-1] + 0.5,
            0,
            1,
        ],
    )

    ax.set_yticks([])
    ax.set_ylabel(
        "Ref.\nSS",
        fontsize=t.small_fontsize,
        rotation=0,
        ha="right",
        va="center",
        fontstyle="italic",
    )
    ax.set_xlabel("Residue Number", fontsize=t.label_fontsize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Place SS legend outside the bar, to the right
    legend_patches = [Patch(facecolor=ss_colors[i], label=ss_names[i]) for i in [1, 2, 0]]
    apply_legend(
        ax,
        plot_settings,
        fontsize=t.small_fontsize,
        handles=legend_patches,
        ncol=1,
        framealpha=0.8,
        borderpad=0.4,
        handlelength=1.0,
        title="Reference SS",
        title_fontsize=t.small_fontsize,
    )


def _load_rmsf_profile(aggregated_dir: Path) -> dict | None:
    """Load per-residue RMSF data from aggregated directory.

    Returns
    -------
    dict or None
        {"residues": [...], "rmsf": [...], "sem": [...]}
    """
    import json

    result_file = aggregated_dir / "rmsf_aggregated.json"
    if not result_file.exists():
        json_files = list(aggregated_dir.glob("*.json"))
        if json_files:
            result_file = json_files[0]
        else:
            return None

    try:
        with open(result_file) as f:
            data = json.load(f)

        # Check for per-residue data (support multiple key naming conventions)
        if "mean_rmsf_per_residue" in data:
            per_res = data["mean_rmsf_per_residue"]
            return {
                "residues": data.get("residue_ids", list(range(1, len(per_res) + 1))),
                "rmsf": per_res,
                "sem": data.get("sem_rmsf_per_residue", []),
            }
        elif "per_residue_rmsf" in data:
            per_res = data["per_residue_rmsf"]
            return {
                "residues": data.get("residue_ids", list(range(1, len(per_res) + 1))),
                "rmsf": per_res,
                "sem": data.get("per_residue_sem", []),
            }
        elif "residue_rmsf" in data:
            return {
                "residues": data.get("residue_ids", list(range(len(data["residue_rmsf"])))),
                "rmsf": data["residue_rmsf"],
                "sem": data.get("residue_sem", []),
            }

        return None

    except Exception as e:
        logger.debug(f"Failed to load RMSF profile from {result_file}: {e}")
        return None


# ---------------------------------------------------------------------------
# Private helper functions (extracted from legacy RMSFCalculator)
# ---------------------------------------------------------------------------


def _compute_rmsd_timeseries(
    u: "Universe",
    atoms: Any,
    start_frame: int,
) -> NDArray[np.float64]:
    """Compute RMSD timeseries for autocorrelation analysis.

    Parameters
    ----------
    u : Universe
        MDAnalysis universe (should be aligned).
    atoms : AtomGroup
        Atom selection for RMSD calculation.
    start_frame : int
        First frame to include (0-indexed).

    Returns
    -------
    NDArray[np.float64]
        RMSD values for each frame from *start_frame* onward.
    """
    u.trajectory[start_frame]
    ref_pos = atoms.positions.copy()

    rmsd_values = []
    for _ts in u.trajectory[start_frame:]:
        diff = atoms.positions - ref_pos
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        rmsd_values.append(rmsd)

    return np.array(rmsd_values, dtype=np.float64)


def _compute_rmsf(
    u: "Universe",
    atoms: Any,
    frame_indices: NDArray[np.int64],
    reference_positions: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Compute RMSF using selected frames.

    Parameters
    ----------
    u : Universe
        Aligned MDAnalysis Universe.
    atoms : AtomGroup
        Atom selection for RMSF calculation.
    frame_indices : NDArray[np.int64]
        Frame indices to use for the calculation.
    reference_positions : NDArray[np.float64] or None
        External reference positions (n_atoms, 3) to use instead of the
        trajectory average.

    Returns
    -------
    NDArray[np.float64]
        Per-atom RMSF values in Angstroms.
    """
    n_frames = len(frame_indices)

    if reference_positions is not None:
        avg_positions = reference_positions
    else:
        positions_sum = np.zeros_like(atoms.positions)
        for idx in frame_indices:
            u.trajectory[int(idx)]
            positions_sum += atoms.positions
        avg_positions = positions_sum / n_frames

    sq_diff_sum = np.zeros(len(atoms), dtype=np.float64)
    for idx in frame_indices:
        u.trajectory[int(idx)]
        diff = atoms.positions - avg_positions
        sq_diff_sum += np.sum(diff**2, axis=1)

    return np.sqrt(sq_diff_sum / n_frames)


def _aggregate_per_residue(
    atoms: Any,
    atom_rmsf: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Aggregate per-atom RMSF to per-residue (mean within residue).

    Parameters
    ----------
    atoms : AtomGroup
        MDAnalysis atom selection.
    atom_rmsf : NDArray[np.float64]
        Per-atom RMSF values.

    Returns
    -------
    NDArray[np.float64]
        Per-residue mean RMSF values.
    """
    residues = atoms.residues
    n_residues = len(residues)
    per_residue = np.zeros(n_residues, dtype=np.float64)

    for i, res in enumerate(residues):
        res_atoms = atoms.select_atoms(f"resid {res.resid}")
        mask = np.isin(atoms.indices, res_atoms.indices)
        per_residue[i] = np.mean(atom_rmsf[mask])

    return per_residue
