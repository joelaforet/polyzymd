"""Engine abstraction layer for simulation backends."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, ClassVar

from pydantic import BaseModel, Field

from polyzymd.simulation.progress import SimulationProgress


class TrajectoryLayout(BaseModel):
    """Canonical trajectory/topology layout emitted by an engine.

    Parameters
    ----------
    topology_path : Path | None
        Path to topology file used for analysis loading.
    trajectory_paths : list[Path]
        Ordered trajectory files to read.
    trajectory_format : str
        Trajectory format identifier, for example ``"dcd"`` or ``"xtc"``.
    topology_format : str
        Topology format identifier, for example ``"pdb"`` or ``"gro"``.
    segment_status : dict[int, str]
        Status recorded for each production segment index by the engine's own
        progress tracking. Empty when the engine does not track per-segment
        status or when the run has no progress file.
    excluded_segments : list[int]
        Segment indices left out of ``trajectory_paths`` because they are not
        complete. Empty when every discovered segment was included.
    incomplete_segments : list[int]
        Segment indices that were included even though they are not complete.
        Non-empty only when the caller asked for incomplete segments.
    empty_segments : list[int]
        Segment indices left out of ``trajectory_paths`` because their
        trajectory file holds no frame.
    missing_segments : list[int]
        Segment indices the engine's progress records as completed whose
        files are not on disk, for example in a copy of a run that kept only
        some segments.
    """

    topology_path: Path | None = None
    #: The GROMACS ``.top`` that replaces an unreadable ``prod.tpr``, by the name the run uses.
    gromacs_topology_path: Path | None = None
    trajectory_paths: list[Path] = Field(default_factory=list)
    trajectory_format: str
    topology_format: str
    segment_status: dict[int, str] = Field(default_factory=dict)
    excluded_segments: list[int] = Field(default_factory=list)
    incomplete_segments: list[int] = Field(default_factory=list)
    empty_segments: list[int] = Field(default_factory=list)
    missing_segments: list[int] = Field(default_factory=list)


class EngineSubmitRequest(BaseModel):
    """Submission request information for engine HPC execution.

    Parameters
    ----------
    replicate : int
        Replicate index.
    config_path : Path
        Path to simulation configuration file.
    working_dir : Path
        Engine working directory for this replicate.
    job_name : str | None, optional
        Scheduler job name override.
    slurm_config : Any, optional
        Scheduler-specific config object.
    extra : dict[str, Any], optional
        Engine-specific metadata for submission handling.
    """

    replicate: int
    config_path: Path
    working_dir: Path
    job_name: str | None = None
    slurm_config: Any = None
    extra: dict[str, Any] = Field(default_factory=dict)


class SimulationEngine(ABC):
    """Abstract interface for simulation execution engines."""

    name: ClassVar[str]
    engine_subdir: ClassVar[str | None] = None

    def get_engine_working_directory(self, sim_config: object, replicate: int) -> Path:
        """Resolve the engine-specific working directory for a replicate.

        Combines the shared scratch-based replicate root from
        ``sim_config.get_working_directory(replicate)`` with this engine's
        ``engine_subdir`` (if any).

        Parameters
        ----------
        sim_config : object
            Simulation configuration with ``get_working_directory`` method.
        replicate : int
            Replicate index.

        Returns
        -------
        Path
            Engine-specific working directory.
        """
        root = sim_config.get_working_directory(replicate)
        if self.engine_subdir:
            return root / self.engine_subdir
        return root

    def resolve_engine_working_directory(self, replicate_root: Path) -> Path:
        """Append the engine subdirectory to a discovered replicate root.

        Used by ``status`` and ``recover`` to map a discovered replicate
        directory to the engine-specific working directory.

        Parameters
        ----------
        replicate_root : Path
            Replicate root directory (e.g. from ``discover_replicate_dirs``).

        Returns
        -------
        Path
            Engine working directory (root / engine_subdir, or root itself).
        """
        if self.engine_subdir:
            return replicate_root / self.engine_subdir
        return replicate_root

    @classmethod
    @abstractmethod
    def from_config(cls, config: object) -> SimulationEngine:
        """Create an engine instance from a simulation config.

        Parameters
        ----------
        config : object
            Simulation configuration object.

        Returns
        -------
        SimulationEngine
            Configured engine instance.
        """

    @abstractmethod
    def run_local(self, replicate: int, working_dir: Path, skip_build: bool = False) -> None:
        """Run the simulation locally for a replicate.

        Parameters
        ----------
        replicate : int
            Replicate index.
        working_dir : Path
            Working directory for simulation outputs.
        skip_build : bool, optional
            Whether to skip system construction if cached files exist.
        """

    @abstractmethod
    def prepare_submission(self, request: EngineSubmitRequest) -> Path:
        """Prepare scheduler artifacts for a submission request.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Path
            Path to generated submission script.
        """

    @abstractmethod
    def submit(self, request: EngineSubmitRequest) -> Any:
        """Submit a simulation job to a scheduler.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Any
            Scheduler-specific submission result.
        """

    @abstractmethod
    def load_or_scan_progress(self, working_dir: Path, replicate: int) -> SimulationProgress:
        """Load persisted progress or reconstruct it from output files.

        Parameters
        ----------
        working_dir : Path
            Working directory for this replicate.
        replicate : int
            Replicate index.

        Returns
        -------
        SimulationProgress
            Current progress model for the replicate.
        """

    def trajectory_files(
        self, working_dir: Path, progress: SimulationProgress | None
    ) -> list[Path]:
        """Return the finished trajectory files of a run, whose content identifies it.

        ``polyzymd hash-trajectories`` hashes these, so an engine returns only
        files that will not be written again: the ones analyses read. The
        default returns none, meaning the engine records no trajectory hashes.

        Parameters
        ----------
        working_dir : Path
            Engine working directory of one replicate.
        progress : SimulationProgress or None
            The run's ``progress.json``, or None when it has none.

        Returns
        -------
        list of Path
            Existing finished trajectory files, in the engine's order.
        """
        return []

    def recorded_trajectory_hashes(
        self, working_dir: Path, progress: SimulationProgress | None = None
    ) -> dict[Path, tuple[str, int]]:
        """Return the SHA-256 and size recorded for a run's trajectory files, by resolved path.

        The default reads ``trajectory_hashes.json``, which
        ``polyzymd hash-trajectories`` writes. An engine whose runner also
        records hashes, such as OpenMM per segment, overrides this so the
        runner's hash takes precedence.
        """
        from polyzymd.simulation.progress import load_trajectory_hashes

        working_dir = Path(working_dir)
        return {
            (working_dir / key).resolve(): (item.sha256, item.bytes)
            for key, item in load_trajectory_hashes(working_dir).items()
        }

    def record_trajectory_hashes(
        self,
        working_dir: Path,
        *,
        verify: bool = False,
        dry_run: bool = False,
        rehash_changed: bool = False,
    ) -> dict[str, Any]:
        """Record the SHA-256 of each finished trajectory file of a run in ``trajectory_hashes.json``.

        The same for every engine; each engine says which files
        (:meth:`trajectory_files`) and what is already recorded
        (:meth:`recorded_trajectory_hashes`). ``progress.json`` is read, never
        written, so the run's segment statuses stay as the runner left them.
        It is idempotent: a file whose recorded hash has the file's size is
        left as it is without reading it, so a second call changes nothing,
        and ``trajectory_hashes.json`` is written once, atomically, only when
        a hash was added. A recorded hash is never overwritten: a recorded
        size that differs from the file's, or with ``verify`` a recomputed
        hash that differs, is reported as a conflict. With
        ``rehash_changed``, an entry of ``trajectory_hashes.json`` whose size
        differs, such as a GROMACS run extended after it was hashed, is
        hashed again and replaced.

        Returns
        -------
        dict
            ``hashed``, ``recorded`` (already present), ``verified`` and
            ``rehashed`` as lists of paths relative to ``working_dir``,
            ``conflicts`` as messages, and ``skipped`` with the reason when
            the run was left alone.
        """
        from polyzymd.simulation.progress import (
            TrajectoryHash,
            load_progress,
            load_trajectory_hashes,
            save_trajectory_hashes,
            trajectory_digest,
        )

        working_dir = Path(working_dir)
        report: dict[str, Any] = {
            "hashed": [],
            "recorded": [],
            "verified": [],
            "rehashed": [],
            "conflicts": [],
            "skipped": None,
        }
        if not working_dir.is_dir():
            report["skipped"] = f"no engine working directory {working_dir.name}"
            return report
        progress = load_progress(working_dir)
        stored = load_trajectory_hashes(working_dir)
        recorded = self.recorded_trajectory_hashes(working_dir, progress)
        changed = False
        for path in self.trajectory_files(working_dir, progress):
            name = path.relative_to(working_dir).as_posix()
            size = path.stat().st_size
            known = recorded.get(path.resolve())
            if known is not None and known[1] != size:
                own = name in stored and (stored[name].sha256, stored[name].bytes) == known
                if not (rehash_changed and own):
                    where = "trajectory_hashes.json" if own else "progress.json"
                    report["conflicts"].append(
                        f"{name}: {where} records {known[1]} bytes, the file has {size}; "
                        "the trajectory changed after it was recorded"
                        + (" (--rehash-changed records it again)" if own else "")
                    )
                    continue
                known = None
                report["rehashed"].append(name)
            elif known is not None:
                if not verify:
                    report["recorded"].append(name)
                elif trajectory_digest(path)["trajectory_sha256"] == known[0]:
                    report["verified"].append(name)
                else:
                    report["conflicts"].append(f"{name}: the SHA-256 differs from the one recorded")
                continue
            if name not in report["rehashed"]:
                report["hashed"].append(name)
            if not dry_run:
                digest = trajectory_digest(path)
                stored[name] = TrajectoryHash(
                    sha256=digest["trajectory_sha256"], bytes=digest["trajectory_bytes"]
                )
                changed = True
        if changed:
            save_trajectory_hashes(working_dir, stored)
        return report

    @abstractmethod
    def resolve_trajectory_layout(
        self,
        working_dir: Path,
        replicate: int,
        *,
        require_complete: bool = True,
    ) -> TrajectoryLayout:
        """Resolve trajectory files and topology for downstream analysis.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.
        require_complete : bool, optional
            Leave out segments the engine records as still running or failed,
            by default True. Engines without per-segment status ignore it.

        Returns
        -------
        TrajectoryLayout
            Engine-specific file layout normalized for analysis.
        """
