"""OpenMM engine adapter over existing PolyzyMD simulation workflow."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, ClassVar

from polyzymd.engines.base import EngineSubmitRequest, SimulationEngine, TrajectoryLayout
from polyzymd.simulation.progress import (
    SegmentStatus,
    SimulationProgress,
    calculate_report_interval,
    load_progress,
)

LOGGER = logging.getLogger(__name__)

# A segment in one of these states has no finished trajectory: ``running`` is
# still being appended to and ``failed`` stopped at an arbitrary step. Reading
# either one gives a trajectory that is short for a reason the analysis cannot
# see. ``interrupted`` is excluded from this set on purpose, because the
# continuation chain resumes from an interrupted segment's saved state, so its
# frames are part of the time line.
_INCOMPLETE_STATUSES = frozenset({SegmentStatus.RUNNING, SegmentStatus.FAILED})


def _topology_format(path: Path | None) -> str:
    """Return the layout format name for a topology path."""
    if path is not None and path.suffix.lower() == ".prmtop":
        return "prmtop"
    return "pdb"


class OpenMMEngine(SimulationEngine):
    """Thin adapter for the OpenMM execution path."""

    name: ClassVar[str] = "openmm"
    # OpenMM stores outputs directly in the replicate root

    def __init__(self, config: object):
        """Initialize the engine adapter.

        Parameters
        ----------
        config : object
            Simulation configuration object.
        """
        self._config = config

    @classmethod
    def from_config(cls, config: object) -> OpenMMEngine:
        """Create an OpenMM engine instance from configuration.

        Parameters
        ----------
        config : object
            Simulation configuration object.

        Returns
        -------
        OpenMMEngine
            Configured OpenMM engine.
        """
        return cls(config=config)

    def run_local(self, replicate: int, working_dir: Path, skip_build: bool = False) -> None:
        """Run the OpenMM simulation locally.

        Parameters
        ----------
        replicate : int
            Replicate index.
        working_dir : Path
            Working directory for simulation output.
        skip_build : bool, optional
            Whether to skip build and reuse existing artifacts.
        """
        from polyzymd.cli.main import _run_initial_segment

        prod = self._config.simulation_phases.production
        total_steps = int(prod.duration * 1e6 / prod.time_step)
        report_interval = calculate_report_interval(total_steps, prod.samples)
        _run_initial_segment(
            sim_config=self._config,
            working_dir=working_dir,
            replicate=replicate,
            skip_build=skip_build,
            duration_ns=prod.duration,
            num_samples=prod.samples,
            timestep_fs=prod.time_step,
            report_interval=report_interval,
            checkpoint_interval_s=prod.checkpoint_interval,
        )

    def prepare_submission(self, request: EngineSubmitRequest) -> Path:
        """Prepare a SLURM submission script for OpenMM.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Path
            Path to the generated SLURM script.
        """
        from polyzymd.workflow.slurm import SlurmScriptGenerator

        if request.slurm_config is None:
            raise ValueError("OpenMM submission requires slurm_config")

        script_dir = request.working_dir / "daisy_chain_scripts"
        script_dir.mkdir(parents=True, exist_ok=True)
        script_path = script_dir / f"run_rep{request.replicate}.sh"

        generator = SlurmScriptGenerator(config=request.slurm_config)
        script = generator.generate_job_script(
            config_path=str(request.config_path),
            replicate=request.replicate,
            working_dir=str(request.working_dir),
            job_name=request.job_name,
        )
        generator.save_script(script, script_path)
        return script_path

    def submit(self, request: EngineSubmitRequest) -> Any:
        """Submit an OpenMM replicate through the daisy-chain workflow.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Any
            Submission result object from the existing workflow.
        """
        from polyzymd.workflow.daisy_chain import DaisyChainConfig, DaisyChainSubmitter

        if request.slurm_config is None:
            raise ValueError("OpenMM submission requires slurm_config")

        dc_config = DaisyChainConfig.from_simulation_config(
            sim_config=self._config,
            slurm_config=request.slurm_config,
            replicates=[request.replicate],
            output_script_dir=request.working_dir / "daisy_chain_scripts",
            config_path=str(request.config_path),
        )
        submitter = DaisyChainSubmitter(self._config, dc_config)
        return submitter.submit_replicate(request.replicate)

    def load_or_scan_progress(self, working_dir: Path, replicate: int) -> SimulationProgress:
        """Load OpenMM progress from progress.json or filesystem scan.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        SimulationProgress
            Current progress state for the replicate.
        """
        from polyzymd.simulation.progress import load_or_scan_progress

        prod = self._config.simulation_phases.production
        total_steps = int(prod.duration * 1e6 / prod.time_step)

        return load_or_scan_progress(
            working_dir=working_dir,
            config_path="",
            total_steps=total_steps,
            total_samples=prod.samples,
            timestep_fs=prod.time_step,
            replicate=replicate,
        )

    def resolve_trajectory_layout(
        self,
        working_dir: Path,
        replicate: int,
        *,
        require_complete: bool = True,
    ) -> TrajectoryLayout:
        """Resolve OpenMM trajectory and topology paths.

        Uses the canonical OpenMM output layout rooted at ``working_dir``. The
        status recorded for each segment in ``progress.json`` decides whether
        the segment is read. A segment marked ``running`` or ``failed`` is
        still being written or was abandoned mid-write, so its DCD ends at an
        arbitrary frame; including it would silently shorten the analysis
        window. Such segments are left out unless ``require_complete`` is
        False, and either way the status of every segment is reported on the
        layout.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index (unused, kept for interface parity).
        require_complete : bool, optional
            Leave out segments recorded as running or failed, by default True.

        Returns
        -------
        TrajectoryLayout
            Resolved DCD/PDB layout for downstream analyses.
        """
        _ = replicate

        topology_path = self._find_openmm_topology(working_dir)
        trajectory_paths, segment_status, skipped = self._find_openmm_trajectories(
            working_dir, require_complete=require_complete
        )

        return TrajectoryLayout(
            topology_path=topology_path,
            trajectory_paths=trajectory_paths,
            trajectory_format="dcd",
            topology_format=_topology_format(topology_path),
            segment_status=segment_status,
            excluded_segments=skipped if require_complete else [],
            incomplete_segments=[] if require_complete else skipped,
        )

    @staticmethod
    def _find_openmm_topology(working_dir: Path) -> Path | None:
        """Find the topology using the canonical OpenMM search order.

        ``system.prmtop`` is preferred when the build wrote it, because it
        carries every bond and has no atom limit. ``solvated_system.pdb`` is
        the fallback for runs built before it existed; above 99,999 atoms its
        CONECT records are unreadable, so run ``polyzymd analysis-topology``
        on such a run first.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.

        Returns
        -------
        Path or None
            Topology path, or None if not found.
        """
        candidate = working_dir / "system.prmtop"
        if candidate.exists():
            return candidate

        candidate = working_dir / "solvated_system.pdb"
        if candidate.exists():
            return candidate

        candidate = working_dir / "production_0" / "production_0_topology.pdb"
        if candidate.exists():
            return candidate

        # Keep exact read-only support for expensive JRL 2025 LipA pre-PolyzyMD data
        candidate = working_dir / "production" / "production_topology.pdb"
        if candidate.exists():
            return candidate

        # Arbitrary topology discovery is intentionally disallowed to avoid
        # analyzing unrelated PDBs in mixed scratch or archival directories
        return None

    @staticmethod
    def _find_openmm_trajectories(
        working_dir: Path,
        *,
        require_complete: bool = True,
    ) -> tuple[list[Path], dict[int, str], list[int]]:
        """Find trajectory DCD files using the canonical OpenMM search order.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        require_complete : bool, optional
            Leave out segments recorded as running or failed, by default True.

        Returns
        -------
        tuple
            Ordered trajectory files, the status recorded for each segment
            index, and the indices of the segments that are not complete.
        """
        prod_re = re.compile(r"production_(\d+)$")
        segment_dirs = {
            int(match.group(1)): path
            for path in working_dir.iterdir()
            if path.is_dir() and (match := prod_re.fullmatch(path.name)) is not None
        }

        if segment_dirs:
            max_index = max(segment_dirs)
            trajectory_paths = []
            progress = load_progress(working_dir)
            segments = progress.segments if progress is not None else []
            segment_status = {segment.index: segment.status.value for segment in segments}
            zero_frame_segments = {
                segment.index
                for segment in segments
                if segment.status == SegmentStatus.COMPLETED and segment.samples_written == 0
            }
            incomplete_segments = [
                segment.index for segment in segments if segment.status in _INCOMPLETE_STATUSES
            ]
            for index in range(max_index + 1):
                segment_dir = working_dir / f"production_{index}"
                file_path = segment_dir / f"production_{index}_trajectory.dcd"
                if not segment_dir.is_dir():
                    raise ValueError(f"Missing OpenMM production segment directory: {segment_dir}")
                if index in zero_frame_segments and (
                    not file_path.is_file() or file_path.stat().st_size == 0
                ):
                    continue
                if require_complete and index in incomplete_segments:
                    LOGGER.warning(
                        "Excluding OpenMM production segment %d from analysis: "
                        "progress.json records it as %s, so %s may be mid-write. "
                        "Pass require_complete=False to read it anyway.",
                        index,
                        segment_status.get(index, "incomplete"),
                        file_path,
                    )
                    continue
                if not file_path.is_file():
                    raise ValueError(f"Missing OpenMM trajectory segment: {file_path}")
                if file_path.stat().st_size == 0:
                    raise ValueError(f"Empty OpenMM trajectory segment: {file_path}")
                trajectory_paths.append(file_path)
            included = {
                int(path.parent.name.removeprefix("production_")) for path in trajectory_paths
            }
            if require_complete:
                reported = [index for index in incomplete_segments if index not in included]
            else:
                reported = [index for index in incomplete_segments if index in included]
            return trajectory_paths, segment_status, reported

        # Keep exact read-only support for expensive JRL 2025 LipA pre-PolyzyMD data
        single_production = working_dir / "production" / "production_trajectory.dcd"
        if single_production.exists():
            if not single_production.is_file():
                raise ValueError(f"OpenMM trajectory path is not a file: {single_production}")
            if single_production.stat().st_size == 0:
                raise ValueError(f"Empty OpenMM trajectory: {single_production}")
            return [single_production], {}, []

        # Broad recursive globs are intentionally disallowed so old datasets
        # must match approved legacy names instead of accidental local files
        return [], {}, []
