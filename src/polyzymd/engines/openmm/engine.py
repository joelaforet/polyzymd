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


def _dcd_has_no_frames(path: Path) -> bool:
    """Return True when a DCD file is empty or holds only its header.

    A DCD file starts with an 84-byte record, a title record and a one-integer
    atom-count record, as written by OpenMM's ``DCDFile``. A file whose size
    is at most the end of those records holds no frame. A file whose first
    record marker is not 84 in either byte order is not read here and counts
    as holding frames, so the trajectory reader reports what is wrong with it.
    """
    import struct

    size = path.stat().st_size
    if size == 0:
        return True
    try:
        with path.open("rb") as handle:
            head = handle.read(96)
            for order in ("<", ">"):
                if len(head) >= 96 and struct.unpack(order + "i", head[:4])[0] == 84:
                    (title,) = struct.unpack(order + "i", head[92:96])
                    return size <= 92 + 4 + title + 4 + 12
    except OSError:
        return False
    return False


def _topology_format(path: Path | None) -> str:
    """Return the layout format name for a topology path."""
    if path is not None and path.suffix.lower() == ".prmtop":
        return "prmtop"
    return "pdb"


def _missing_completed_segments(working_dir: Path) -> list[int]:
    """Return the production segments ``progress.json`` records as run that are not on disk.

    A segment recorded as completed or interrupted wrote frames (an
    interrupted one stopped early, and its sample count is not updated), so
    its missing folder means missing production. A completed segment that
    wrote no frame is not counted.
    """
    progress = load_progress(working_dir)
    if progress is None:
        return []
    ran = {SegmentStatus.COMPLETED, SegmentStatus.INTERRUPTED}
    return sorted(
        segment.index
        for segment in progress.segments
        if segment.status in ran
        and not (segment.status == SegmentStatus.COMPLETED and segment.samples_written == 0)
        and not (working_dir / f"production_{segment.index}").is_dir()
    )


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

    def trajectory_files(
        self, working_dir: Path, progress: SimulationProgress | None
    ) -> list[Path]:
        """Return the DCD of each production segment that analyses read.

        With ``progress.json``, a segment recorded as running or failed is left
        out, as :meth:`resolve_trajectory_layout` leaves it out. Without one,
        as in a downsampled copy, every segment's DCD is returned.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        progress : SimulationProgress or None
            The run's progress, or None when it has no ``progress.json``.

        Returns
        -------
        list of Path
            ``production_<n>/production_<n>_trajectory.dcd`` files, by segment index.
        """
        left_out = {
            s.index
            for s in (progress.segments if progress else [])
            if s.status in _INCOMPLETE_STATUSES
        }
        found = []
        for folder in Path(working_dir).glob("production_*"):
            match = re.fullmatch(r"production_(\d+)", folder.name)
            if match is None or int(match.group(1)) in left_out:
                continue
            dcd = folder / f"{folder.name}_trajectory.dcd"
            if dcd.is_file():
                found.append((int(match.group(1)), dcd))
        return [path for _, path in sorted(found)]

    def recorded_trajectory_hashes(
        self, working_dir: Path, progress: SimulationProgress | None = None
    ) -> dict[Path, tuple[str, int]]:
        """Return recorded hashes, the runner's hash of each segment taking precedence.

        The runner records a segment's hash when the segment completes, so it
        is newer than an entry ``hash-trajectories`` wrote before the segment
        was resumed.
        """
        from polyzymd.simulation.progress import load_progress

        working_dir = Path(working_dir)
        progress = progress if progress is not None else load_progress(working_dir)
        recorded = super().recorded_trajectory_hashes(working_dir, progress)
        for segment in progress.segments if progress else []:
            if segment.trajectory_sha256 and segment.trajectory_bytes is not None:
                path = (
                    working_dir
                    / f"production_{segment.index}"
                    / f"production_{segment.index}_trajectory.dcd"
                ).resolve()
                recorded[path] = (segment.trajectory_sha256, segment.trajectory_bytes)
        return recorded

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
        layout. A segment whose DCD file exists but holds no frame, such as
        one interrupted at start-up before its first report and resumed by
        the next segment, is left out with a warning and listed in
        ``empty_segments``; the loader still checks that the remaining
        segments form one contiguous time line.

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
        trajectory_paths, segment_status, skipped, empty = self._find_openmm_trajectories(
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
            empty_segments=empty,
            missing_segments=_missing_completed_segments(working_dir),
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
    ) -> tuple[list[Path], dict[int, str], list[int], list[int]]:
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
            index, the indices of the segments that are not complete, and the
            indices of the segments left out because their DCD file holds no
            frame.
        """
        prod_re = re.compile(r"production_(\d+)$")
        segment_dirs = {
            int(match.group(1)): path
            for path in working_dir.iterdir()
            if path.is_dir() and (match := prod_re.fullmatch(path.name)) is not None
        }

        if segment_dirs:
            max_index = max(segment_dirs)
            trajectory_paths, empty = [], []
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
                if _dcd_has_no_frames(file_path):
                    LOGGER.warning(
                        "Skipping OpenMM production segment %d: %s holds no frames. The "
                        "remaining segments must still form one contiguous time line.",
                        index,
                        file_path,
                    )
                    empty.append(index)
                    continue
                trajectory_paths.append(file_path)
            included = {
                int(path.parent.name.removeprefix("production_")) for path in trajectory_paths
            }
            if require_complete:
                reported = [index for index in incomplete_segments if index not in included]
            else:
                reported = [index for index in incomplete_segments if index in included]
            return trajectory_paths, segment_status, reported, empty

        # Keep exact read-only support for expensive JRL 2025 LipA pre-PolyzyMD data
        single_production = working_dir / "production" / "production_trajectory.dcd"
        if single_production.exists():
            if not single_production.is_file():
                raise ValueError(f"OpenMM trajectory path is not a file: {single_production}")
            if single_production.stat().st_size == 0:
                raise ValueError(f"Empty OpenMM trajectory: {single_production}")
            return [single_production], {}, [], []

        # Broad recursive globs are intentionally disallowed so old datasets
        # must match approved legacy names instead of accidental local files
        return [], {}, [], []
