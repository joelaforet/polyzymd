"""Engine abstraction layer for simulation backends."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

from pydantic import BaseModel, Field

from polyzymd.simulation.progress import SimulationProgress

if TYPE_CHECKING:
    from polyzymd.engines.bias import EnhancedSamplingProtocol, ExternalBiasSpec


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

    @abstractmethod
    def trajectory_files(
        self, working_dir: Path, progress: SimulationProgress | None
    ) -> list[Path]:
        """Return the finished trajectory files of a run, whose content identifies it.

        Parameters
        ----------
        working_dir : Path
            Engine working directory of one replicate.
        progress : SimulationProgress or None
            The run's progress, when known, to leave out files still being
            written.

        Returns
        -------
        list of Path
            Existing trajectory files, in the engine's order.
        """

    def recorded_trajectory_hashes(
        self, working_dir: Path, progress: SimulationProgress | None = None
    ) -> dict[Path, tuple[str, int]]:
        """Return the SHA-256 and size the run recorded for its trajectory files, by path.

        The default reads the engine-neutral ``trajectory_hashes`` of
        ``progress.json``; an engine that also records hashes elsewhere, such
        as per segment, adds them by overriding this method.
        """
        from polyzymd.simulation.progress import load_progress

        working_dir = Path(working_dir)
        progress = progress if progress is not None else load_progress(working_dir)
        if progress is None:
            return {}
        return {
            (working_dir / key).resolve(): (item.sha256, item.bytes)
            for key, item in progress.trajectory_hashes.items()
        }

    def store_trajectory_hash(
        self, progress: SimulationProgress, working_dir: Path, path: Path, sha256: str, size: int
    ) -> None:
        """Record one file's hash in ``progress``: in the engine-neutral ``trajectory_hashes``."""
        from polyzymd.simulation.progress import TrajectoryHash

        key = str(Path(path).resolve().relative_to(Path(working_dir).resolve()))
        progress.trajectory_hashes[key] = TrajectoryHash(sha256=sha256, bytes=size)

    def record_trajectory_hashes(
        self,
        working_dir: Path,
        replicate: int,
        *,
        verify: bool = False,
        dry_run: bool = False,
        force: bool = False,
    ) -> dict[str, Any]:
        """Record the SHA-256 of each finished trajectory file of a run in its ``progress.json``.

        The same for every engine; each engine says which files
        (:meth:`trajectory_files`) and what was recorded
        (:meth:`recorded_trajectory_hashes`, :meth:`store_trajectory_hash`).
        It is idempotent: a file whose recorded hash has the file's size is
        left as it is without reading it, so a second call changes nothing.
        A recorded hash is never overwritten: a recorded size that differs
        from the file's, or with ``verify`` a recomputed hash that differs, is
        reported as a conflict. A run recorded as running, whose job may still
        write ``progress.json``, is left alone unless ``force``.
        ``progress.json`` is written once, atomically, only when a hash was
        added and not ``dry_run``, and nothing else in it changes; a run
        without one gets one from the engine's scan of its files.

        Returns
        -------
        dict
            ``hashed``, ``recorded`` (already present), ``verified`` and
            ``conflicts`` as lists of paths relative to ``working_dir`` or
            messages, ``created`` (whether ``progress.json`` was written for the
            first time), and ``skipped`` with the reason when the run was left
            alone.
        """
        from polyzymd.simulation.progress import (
            SegmentStatus,
            SimulationStatus,
            load_progress,
            save_progress,
            trajectory_digest,
        )

        working_dir = Path(working_dir)
        report: dict[str, Any] = {
            "hashed": [],
            "recorded": [],
            "verified": [],
            "conflicts": [],
            "created": False,
            "skipped": None,
        }
        if not working_dir.is_dir():
            report["skipped"] = f"no engine working directory {working_dir.name}"
            return report
        # The stored progress is used as it is, so only hashes are added to it;
        # the engine scans the files only for a run that has none.
        progress = load_progress(working_dir)
        existed = progress is not None
        if progress is None:
            progress = self.load_or_scan_progress(working_dir, replicate)
        running = progress.status == SimulationStatus.RUNNING or any(
            s.status == SegmentStatus.RUNNING for s in progress.segments
        )
        if running and not force:
            report["skipped"] = (
                "recorded as running; its job may still write progress.json "
                "(use --force once it has stopped)"
            )
            return report
        recorded = self.recorded_trajectory_hashes(working_dir, progress)
        changed = False
        for path in self.trajectory_files(working_dir, progress):
            name = str(path.relative_to(working_dir))
            size = path.stat().st_size
            known = recorded.get(path.resolve())
            if known is not None:
                if known[1] != size:
                    report["conflicts"].append(
                        f"{name}: progress.json records {known[1]} bytes, the file has {size}; "
                        "the trajectory changed after it was recorded"
                    )
                elif verify:
                    if trajectory_digest(path)["trajectory_sha256"] == known[0]:
                        report["verified"].append(name)
                    else:
                        report["conflicts"].append(
                            f"{name}: the SHA-256 differs from the one recorded"
                        )
                else:
                    report["recorded"].append(name)
                continue
            if not dry_run:
                digest = trajectory_digest(path)
                self.store_trajectory_hash(
                    progress,
                    working_dir,
                    path,
                    digest["trajectory_sha256"],
                    digest["trajectory_bytes"],
                )
                changed = True
            report["hashed"].append(name)
        if changed:
            save_progress(working_dir, progress)
            report["created"] = not existed
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

    # ------------------------------------------------------------------
    # Optional extension hooks — enhanced sampling & state persistence
    # ------------------------------------------------------------------
    # These hooks establish the contract for future PLUMED, metadynamics,
    # and checkpoint integrations.  Each raises ``NotImplementedError``
    # with an engine-specific message by default.  Engine subclasses
    # override them when support is added — no abstract decorator needed.
    # ------------------------------------------------------------------

    def _raise_unsupported(self, feature: str) -> None:
        """Raise a standardised ``NotImplementedError`` for unimplemented hooks.

        Parameters
        ----------
        feature : str
            Human-readable description of the unsupported feature.
        """
        raise NotImplementedError(
            f"{feature} is not yet supported for {self.name}. "
            "This extension point will be implemented in a future release."
        )

    def attach_external_bias(
        self,
        request: EngineSubmitRequest,
        bias_spec: ExternalBiasSpec,
    ) -> None:
        """Attach external biasing inputs to the simulation workflow.

        This hook is called during submission to inject bias scripts
        (e.g. PLUMED input files) into the engine's working directory
        and command-line arguments.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.
        bias_spec : ExternalBiasSpec
            Engine-neutral bias definition.

        Raises
        ------
        NotImplementedError
            Always, until a concrete engine implements this hook.
        """
        self._raise_unsupported("External bias attachment")

    def configure_enhanced_sampling(
        self,
        protocol: EnhancedSamplingProtocol,
    ) -> None:
        """Configure an enhanced-sampling protocol for the engine.

        This hook is called before submission to apply protocol-level
        settings (e.g. replica exchange temperature ladder, metadynamics
        hill parameters) to the engine's simulation inputs.

        Parameters
        ----------
        protocol : EnhancedSamplingProtocol
            Engine-neutral enhanced-sampling protocol definition.

        Raises
        ------
        NotImplementedError
            Always, until a concrete engine implements this hook.
        """
        self._raise_unsupported("Enhanced sampling")

    def save_engine_state(self, working_dir: Path, replicate: int) -> Path:
        """Persist an engine-specific restart state artifact.

        Used by the continuation framework to checkpoint the engine's
        internal state (e.g. OpenMM ``state.xml``, GROMACS ``state.cpt``).

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        Path
            Path to the saved state artifact.

        Raises
        ------
        NotImplementedError
            Always, until a concrete engine implements this hook.
        """
        self._raise_unsupported("Engine state persistence")

    def save_bias_state(self, working_dir: Path, replicate: int) -> Path:
        """Persist an external-bias restart state artifact.

        Used by the continuation framework to checkpoint bias-specific
        state (e.g. PLUMED HILLS file, metadynamics collective-variable
        history).

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        Path
            Path to the saved bias state artifact.

        Raises
        ------
        NotImplementedError
            Always, until a concrete engine implements this hook.
        """
        self._raise_unsupported("Bias state persistence")
