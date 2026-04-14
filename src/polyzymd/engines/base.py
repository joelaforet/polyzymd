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
