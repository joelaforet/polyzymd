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
    """

    topology_path: Path | None = None
    trajectory_paths: list[Path] = Field(default_factory=list)
    trajectory_format: str
    topology_format: str


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
    def resolve_trajectory_layout(self, working_dir: Path, replicate: int) -> TrajectoryLayout:
        """Resolve trajectory files and topology for downstream analysis.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        TrajectoryLayout
            Engine-specific file layout normalized for analysis.
        """

    def attach_external_bias(self, request: EngineSubmitRequest, bias_spec: Any) -> None:
        """Attach external biasing inputs to the simulation workflow.

        Planned for future release. Not yet implemented.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.
        bias_spec : Any
            Bias definition payload.
        """
        raise NotImplementedError("External bias integration is not implemented")

    def configure_enhanced_sampling(self, protocol: Any) -> None:
        """Configure an enhanced-sampling protocol.

        Planned for future release. Not yet implemented.

        Parameters
        ----------
        protocol : Any
            Protocol definition payload.
        """
        raise NotImplementedError("Enhanced sampling is not implemented")

    def save_engine_state(self, working_dir: Path, replicate: int) -> Path:
        """Persist an engine-specific restart state artifact.

        Planned for future release. Not yet implemented.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        Path
            Path to the state artifact.
        """
        raise NotImplementedError("Engine state persistence is not implemented")

    def save_bias_state(self, working_dir: Path, replicate: int) -> Path:
        """Persist an external-bias restart state artifact.

        Planned for future release. Not yet implemented.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        Path
            Path to the bias state artifact.
        """
        raise NotImplementedError("Bias state persistence is not implemented")
