"""RMSF analysis plugin.

Computes per-residue Root Mean Square Fluctuation from MD trajectories,
aggregates across replicates, compares conditions via the default scalar
comparison pipeline, and produces comparison/profile plots.

All heavy computation is self-contained within this plugin package.
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
    BasePlotSettings,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.rmsf._plot_settings import RMSFPlotSettings
from polyzymd.analyses.rmsf._plotters import _plot_rmsf_comparison, _plot_rmsf_profile
from polyzymd.analyses.rmsf._results import RMSFAggregatedResult, RMSFResult
from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory
from polyzymd.analyses.shared.autocorrelation import (
    compute_acf,
    estimate_correlation_time,
    get_independent_indices,
)
from polyzymd.analyses.shared.centroid import ReferenceMode
from polyzymd.analyses.shared.config_hash import compute_config_hash
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

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

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
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = RMSFPlotSettings
    AggregatedResultClass: ClassVar[type] = RMSFAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = RMSFResult
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

        selection = settings.selection
        reference_mode: ReferenceMode = settings.reference_mode
        reference_frame = settings.reference_frame
        reference_file = settings.reference_file
        alignment_selection = settings.alignment_selection
        centroid_selection = settings.centroid_selection

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
        eq_str = f"eq{eq_value:g}{eq_unit}"
        result_filename = f"rmsf_{eq_str}.json"
        result_file = output_dir / result_filename

        cached = self._check_cache(
            RMSFResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

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
            selection_string=settings.selection,
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

        data, labels = self._build_plot_data(ctx)
        if not labels:
            return plots

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        plot_settings = ctx.plot_settings

        result = _plot_rmsf_comparison(data, labels, ctx.output_dir, plot_settings)
        plots.extend(result)

        result = _plot_rmsf_profile(data, labels, ctx.output_dir, plot_settings)
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
        return f"rmsf_{rep_str}_{eq_str}.json"


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

    atom_indices_set = set(atoms.indices)
    for i, res in enumerate(residues):
        # Keep residue grouping chain-aware when resid values repeat across chains
        res_atom_indices = [idx for idx in res.atoms.indices if idx in atom_indices_set]
        if res_atom_indices:
            mask = np.isin(atoms.indices, res_atom_indices)
            per_residue[i] = np.mean(atom_rmsf[mask])

    return per_residue
