"""Secondary structure calculator using mdtraj DSSP.

Computes per-frame, per-residue secondary structure assignments using
``mdtraj.compute_dssp(simplified=True)`` and aggregates persistence
fractions across replicates.

The calculator follows the same patterns as :class:`RMSFCalculator`:

- Config-based trajectory loading via :class:`TrajectoryLoader`
- Equilibration time skipping
- Per-replicate computation with caching
- Multi-replicate aggregation with SEM

Design decisions
----------------
- **mdtraj backend**: mdtraj's DSSP is already a dependency (used for SASA)
  and runs in compiled C, making it significantly faster than pure-Python
  alternatives.
- **Simplified 3-state**: ``simplified=True`` maps all helix types to ``H``,
  all strand types to ``E``, and everything else to ``C``.  This is
  appropriate for tracking overall secondary structure content.
- **NPZ sidecar**: Per-frame assignment matrices (shape ``(n_frames,
  n_residues)``, dtype ``int8``) are stored as compressed NPZ files
  alongside the JSON result.  This keeps JSON results small while
  preserving full temporal resolution for timeline heatmaps.
- **Integer encoding**: 0=coil, 1=helix, 2=strand (see ``results.py``).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

import numpy as np

from polyzymd.analysis.core.aggregation import collect_replicate_results
from polyzymd.analysis.core.config_hash import compute_config_hash, validate_config_hash
from polyzymd.analysis.core.loader import TrajectoryLoader, convert_time, parse_time_string
from polyzymd.analysis.core.registry import AnalyzerRegistry, BaseAnalysisSettings, BaseAnalyzer
from polyzymd.analysis.results.base import get_polyzymd_version
from polyzymd.analysis.secondary_structure.results import (
    SS_CHAR_TO_INT,
    SecondaryStructureAggregatedResult,
    SecondaryStructureResult,
)

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig

LOGGER = logging.getLogger(__name__)


def _chain_letter_to_index(chain_id: str) -> int:
    """Convert chain letter (A, B, C...) to 0-indexed integer for mdtraj."""
    return ord(chain_id.upper()) - ord("A")


@AnalyzerRegistry.register("secondary_structure")
class SecondaryStructureCalculator(BaseAnalyzer):
    """Calculator for DSSP secondary structure analysis.

    This class handles the complete secondary structure analysis workflow:

    1. Load trajectory from config (via mdtraj)
    2. Skip equilibration frames
    3. Compute DSSP assignments using ``mdtraj.compute_dssp(simplified=True)``
    4. Encode assignments as integer matrix and compute persistence fractions
    5. Save per-replicate result (JSON + NPZ sidecar)
    6. Aggregate across replicates with mean/SEM of persistence

    Parameters
    ----------
    config : SimulationConfig
        PolyzyMD simulation configuration.
    chain_id : str, optional
        Chain letter for the protein chain. Default is ``"A"`` (PolyzyMD
        convention: chain A = protein).
    equilibration : str, optional
        Equilibration time to skip (e.g., ``"200ns"``, ``"5000ps"``).
        Default is ``"0ns"`` (no equilibration skip).

    Examples
    --------
    >>> from polyzymd.config import load_config
    >>> config = load_config("config.yaml")
    >>> calc = SecondaryStructureCalculator(
    ...     config, chain_id="A", equilibration="200ns"
    ... )
    >>> result = calc.compute(replicate=1)
    >>> print(result.summary())
    SS (4800 frames, 181 residues): helix=22.3%, strand=18.1%, coil=59.6%

    >>> agg = calc.compute_aggregated(replicates=[1, 2, 3, 4, 5])
    >>> print(agg.summary())
    """

    def __init__(
        self,
        config: "SimulationConfig",
        chain_id: str = "A",
        equilibration: str = "0ns",
    ) -> None:
        self.config = config
        self.chain_id = chain_id

        # Parse equilibration time
        eq_value, eq_unit = parse_time_string(equilibration)
        self.equilibration_time = eq_value
        self.equilibration_unit = eq_unit

        # Initialize loader (for trajectory path discovery)
        self._loader = TrajectoryLoader(config)
        self._config_hash = compute_config_hash(config)

    @classmethod
    def analysis_type(cls) -> str:
        """Return the unique identifier for this analyzer.

        Returns
        -------
        str
            Analysis type identifier.
        """
        return "secondary_structure"

    @classmethod
    def from_config(
        cls,
        analysis_settings: BaseAnalysisSettings,
        sim_config: "SimulationConfig",
        equilibration: str = "0ns",
    ) -> "SecondaryStructureCalculator":
        """Create a secondary structure calculator from analysis settings.

        Parameters
        ----------
        analysis_settings : BaseAnalysisSettings
            Secondary-structure-compatible settings object.
        sim_config : SimulationConfig
            Simulation configuration for trajectory loading.
        equilibration : str, optional
            Equilibration time to skip.

        Returns
        -------
        SecondaryStructureCalculator
            Configured secondary structure calculator.
        """
        return cls(
            config=sim_config,
            chain_id=getattr(analysis_settings, "chain_id", "A"),
            equilibration=equilibration,
        )

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def compute(
        self,
        replicate: int,
        save: bool = True,
        output_dir: Path | None = None,
        recompute: bool = False,
    ) -> SecondaryStructureResult:
        """Compute DSSP secondary structure for a single replicate.

        Parameters
        ----------
        replicate : int
            Replicate number (1-indexed).
        save : bool, optional
            If True (default), save result JSON + NPZ to disk.
        output_dir : Path, optional
            Directory to save results.  Default is
            ``{projects_dir}/analysis/secondary_structure/run_{replicate}/``.
        recompute : bool, optional
            If True, recompute even if a cached result exists.

        Returns
        -------
        SecondaryStructureResult
            Single-replicate secondary structure result.
        """
        # Determine output path
        if output_dir is None:
            output_dir = (
                self.config.output.projects_directory
                / "analysis"
                / "secondary_structure"
                / f"run_{replicate}"
            )

        result_file = output_dir / self._make_result_filename()

        # Check for cached result
        if not recompute and result_file.exists():
            LOGGER.info(f"Loading cached SS result from {result_file}")
            result = SecondaryStructureResult.load(result_file)
            validate_config_hash(result.config_hash, self.config)
            return result

        LOGGER.info(f"Computing secondary structure for replicate {replicate}")

        # ----------------------------------------------------------
        # Load trajectory with mdtraj
        # ----------------------------------------------------------
        import mdtraj as md

        traj_info = self._loader.get_trajectory_info(replicate)
        traj_files_str = [str(f) for f in traj_info.trajectory_files]
        topology_path = str(traj_info.topology_file)

        traj = md.load(traj_files_str, top=topology_path)
        LOGGER.info(f"Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")

        # ----------------------------------------------------------
        # Select protein chain
        # ----------------------------------------------------------
        chain_idx = _chain_letter_to_index(self.chain_id)
        protein_indices = traj.topology.select(f"chainid {chain_idx}")
        if len(protein_indices) == 0:
            # Fallback: try selecting all protein atoms
            protein_indices = traj.topology.select("protein")
            LOGGER.warning(
                f"Chain '{self.chain_id}' not found; falling back to "
                f"'protein' selection ({len(protein_indices)} atoms)"
            )

        protein_traj = traj.atom_slice(protein_indices)
        LOGGER.info(
            f"Protein sub-trajectory: {protein_traj.n_atoms} atoms, "
            f"{protein_traj.n_residues} residues"
        )

        # ----------------------------------------------------------
        # Skip equilibration frames
        # ----------------------------------------------------------
        n_frames_total = protein_traj.n_frames
        # Use TrajectoryLoader for reliable timestep (mdtraj often has
        # incorrect time metadata for legacy DCD files)
        timestep_ps = self._loader.get_timestep(replicate)
        eq_time_ps = convert_time(self.equilibration_time, self.equilibration_unit, "ps")
        start_frame = int(eq_time_ps / timestep_ps) if timestep_ps > 0 else 0
        start_frame = min(start_frame, n_frames_total - 1)

        if start_frame > 0:
            protein_traj = protein_traj[start_frame:]
            LOGGER.info(
                f"Skipped {start_frame} equilibration frames "
                f"({eq_time_ps / 1000:.1f} ns), "
                f"{protein_traj.n_frames} frames remaining"
            )

        n_frames = protein_traj.n_frames

        # ----------------------------------------------------------
        # Compute DSSP
        # ----------------------------------------------------------
        dssp_raw = md.compute_dssp(protein_traj, simplified=True)
        # dssp_raw is (n_frames, n_residues) of single-char strings: H, E, C

        LOGGER.info(f"DSSP complete: {dssp_raw.shape[0]} frames × {dssp_raw.shape[1]} residues")

        # ----------------------------------------------------------
        # Integer-encode the matrix
        # ----------------------------------------------------------
        ss_matrix = self._encode_dssp_matrix(dssp_raw)
        n_residues = ss_matrix.shape[1]

        # ----------------------------------------------------------
        # Compute persistence fractions (per-residue)
        # ----------------------------------------------------------
        persistence_coil = (ss_matrix == 0).mean(axis=0).tolist()
        persistence_helix = (ss_matrix == 1).mean(axis=0).tolist()
        persistence_strand = (ss_matrix == 2).mean(axis=0).tolist()

        # ----------------------------------------------------------
        # Compute overall content fractions
        # ----------------------------------------------------------
        total_entries = ss_matrix.size
        overall_coil = float(np.sum(ss_matrix == 0)) / total_entries
        overall_helix = float(np.sum(ss_matrix == 1)) / total_entries
        overall_strand = float(np.sum(ss_matrix == 2)) / total_entries

        # ----------------------------------------------------------
        # Collect residue metadata
        # ----------------------------------------------------------
        protein_top = protein_traj.topology
        residue_ids = [res.resSeq for res in protein_top.residues]
        residue_names = [res.name.upper() for res in protein_top.residues]

        # ----------------------------------------------------------
        # Build result
        # ----------------------------------------------------------
        result = SecondaryStructureResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            selection_string=f"chainid {chain_idx} (chain {self.chain_id})",
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

        if save:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            json_path = result.save_with_matrix(
                directory=output_dir,
                matrix=ss_matrix,
                filename_prefix=self._make_result_prefix(),
            )
            LOGGER.info(f"Saved SS result to {json_path}")

        return result

    def compute_aggregated(
        self,
        replicates: Sequence[int],
        save: bool = True,
        output_dir: Path | None = None,
        recompute: bool = False,
    ) -> SecondaryStructureAggregatedResult:
        """Aggregate secondary structure across multiple replicates.

        Computes mean and SEM of persistence fractions and overall content
        fractions across the supplied replicates.

        Parameters
        ----------
        replicates : sequence of int
            Replicate numbers to aggregate (1-indexed).
        save : bool, optional
            If True (default), save result to JSON.
        output_dir : Path, optional
            Directory to save aggregated results.  Default is
            ``{projects_dir}/analysis/secondary_structure/aggregated/``.
        recompute : bool, optional
            If True, recompute per-replicate results even if cached.

        Returns
        -------
        SecondaryStructureAggregatedResult
            Aggregated result with mean/SEM across replicates.

        Notes
        -----
        Missing or problematic replicates are skipped with a warning.
        At least 2 successful replicates are required for aggregation.
        """
        if output_dir is None:
            output_dir = (
                self.config.output.projects_directory
                / "analysis"
                / "secondary_structure"
                / "aggregated"
            )

        # Derive per-replicate output base from aggregated output_dir
        per_rep_base = output_dir.parent

        collection = collect_replicate_results(
            self.compute,
            replicates,
            output_dir_base=per_rep_base,
            save=save,
            recompute=recompute,
        )
        individual_results: list[SecondaryStructureResult] = collection.results
        successful_replicates = collection.successful_replicates

        n_reps = len(individual_results)

        # ----------------------------------------------------------
        # Stack persistence arrays: shape (n_reps, n_residues)
        # ----------------------------------------------------------
        coil_stack = np.array([r.persistence_coil for r in individual_results])
        helix_stack = np.array([r.persistence_helix for r in individual_results])
        strand_stack = np.array([r.persistence_strand for r in individual_results])

        mean_coil = coil_stack.mean(axis=0).tolist()
        mean_helix = helix_stack.mean(axis=0).tolist()
        mean_strand = strand_stack.mean(axis=0).tolist()

        sem_coil = (coil_stack.std(axis=0, ddof=1) / np.sqrt(n_reps)).tolist()
        sem_helix = (helix_stack.std(axis=0, ddof=1) / np.sqrt(n_reps)).tolist()
        sem_strand = (strand_stack.std(axis=0, ddof=1) / np.sqrt(n_reps)).tolist()

        # ----------------------------------------------------------
        # Overall content fractions per replicate
        # ----------------------------------------------------------
        per_rep_helix = [r.overall_helix_fraction for r in individual_results]
        per_rep_strand = [r.overall_strand_fraction for r in individual_results]
        per_rep_coil = [r.overall_coil_fraction for r in individual_results]

        mean_overall_helix = float(np.mean(per_rep_helix))
        mean_overall_strand = float(np.mean(per_rep_strand))
        mean_overall_coil = float(np.mean(per_rep_coil))

        sem_overall_helix = float(np.std(per_rep_helix, ddof=1) / np.sqrt(n_reps))
        sem_overall_strand = float(np.std(per_rep_strand, ddof=1) / np.sqrt(n_reps))
        sem_overall_coil = float(np.std(per_rep_coil, ddof=1) / np.sqrt(n_reps))

        # ----------------------------------------------------------
        # Build aggregated result
        # ----------------------------------------------------------
        agg_result = SecondaryStructureAggregatedResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            selection_string=individual_results[0].selection_string,
            n_replicates=n_reps,
            replicates=successful_replicates,
            residue_ids=individual_results[0].residue_ids,
            residue_names=individual_results[0].residue_names,
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

        if save:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            result_file = output_dir / self._make_aggregated_filename(successful_replicates)
            agg_result.save(result_file)
            LOGGER.info(f"Saved aggregated SS result to {result_file}")

        return agg_result

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
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

    @staticmethod
    def _get_mdtraj_timestep_ps(traj) -> float:
        """Extract the timestep in picoseconds from an mdtraj trajectory.

        mdtraj stores time in picoseconds.  If the trajectory has fewer
        than 2 frames, returns 0.0.

        Parameters
        ----------
        traj : mdtraj.Trajectory
            Loaded trajectory.

        Returns
        -------
        float
            Timestep in picoseconds.
        """
        if traj.n_frames < 2:
            return 0.0
        return float(traj.time[1] - traj.time[0])

    def _make_result_prefix(self) -> str:
        """Generate filename prefix for result files (JSON + NPZ)."""
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"
        return f"secondary_structure_{eq_str}"

    def _make_result_filename(self) -> str:
        """Generate filename for result JSON."""
        return f"{self._make_result_prefix()}.json"

    def _make_aggregated_filename(self, replicates: Sequence[int]) -> str:
        """Generate filename for aggregated result."""
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
        return f"secondary_structure_{rep_str}_{eq_str}.json"
