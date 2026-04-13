"""GROMACS engine implementation for PolyzyMD."""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar

from polyzymd.engines.base import EngineSubmitRequest, SimulationEngine, TrajectoryLayout
from polyzymd.simulation.progress import SimulationProgress

from .binary import resolve_gromacs_binary


class GromacsEngine(SimulationEngine):
    """Phase-1 GROMACS execution adapter."""

    name: ClassVar[str] = "gromacs"

    def __init__(self, config: object, gmx_binary: str = "gmx"):
        """Initialize a GROMACS engine adapter.

        Parameters
        ----------
        config : object
            Simulation configuration object.
        gmx_binary : str, optional
            Resolved GROMACS executable.
        """
        self._config = config
        self._gmx_binary = gmx_binary

    @classmethod
    def from_config(cls, config: object) -> GromacsEngine:
        """Create a GROMACS engine from simulation config.

        Parameters
        ----------
        config : object
            Simulation configuration object.

        Returns
        -------
        GromacsEngine
            Configured GROMACS engine instance.
        """
        gmx_binary = resolve_gromacs_binary(config=config)
        return cls(config=config, gmx_binary=gmx_binary)

    def run_local(self, replicate: int, working_dir: Path, skip_build: bool = False) -> None:
        """Run local GROMACS workflow with exported files.

        Parameters
        ----------
        replicate : int
            Replicate index.
        working_dir : Path
            Working directory with GROMACS input files.
        skip_build : bool, optional
            Build-skip placeholder for interface parity.
        """
        from polyzymd.exporters.gromacs import GromacsRunner

        _ = replicate
        _ = skip_build

        prefix = self._config.generate_system_name()
        eq_mdps = sorted(path.name for path in working_dir.glob("eq_*.mdp"))

        runner = GromacsRunner(
            working_dir=working_dir,
            prefix=prefix,
            equilibration_mdps=eq_mdps,
            gmx_command=self._gmx_binary,
        )
        runner.run_full_workflow()

    def prepare_submission(self, request: EngineSubmitRequest) -> Path:
        """Prepare scheduler artifacts for a GROMACS job.

        Planned for Phase 3.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Path
            Placeholder return type for planned SLURM integration.
        """
        _ = request
        raise NotImplementedError("GROMACS SLURM submission will be implemented in Phase 3")

    def submit(self, request: EngineSubmitRequest) -> Any:
        """Submit GROMACS jobs to scheduler.

        Planned for Phase 3.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Any
            Placeholder scheduler result for planned implementation.
        """
        _ = request
        raise NotImplementedError("GROMACS SLURM submission will be implemented in Phase 3")

    def load_or_scan_progress(self, working_dir: Path, replicate: int) -> SimulationProgress:
        """Load or reconstruct GROMACS progress state.

        Planned for Phase 3.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        SimulationProgress
            Planned engine progress model.
        """
        _ = working_dir
        _ = replicate
        raise NotImplementedError("GROMACS progress tracking will be implemented in Phase 3")

    def resolve_trajectory_layout(self, working_dir: Path, replicate: int) -> TrajectoryLayout:
        """Resolve GROMACS trajectory layout for downstream analyses.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        TrajectoryLayout
            Resolved XTC layout with PDB topology fallback to GRO.
        """
        _ = replicate

        trajectory_paths = sorted(working_dir.glob("*.xtc"))

        pdb_candidates = sorted(working_dir.glob("*.pdb"))
        gro_candidates = sorted(working_dir.glob("*.gro"))

        topology_path: Path | None = None
        topology_format = "pdb"
        if pdb_candidates:
            topology_path = pdb_candidates[0]
        elif gro_candidates:
            topology_path = gro_candidates[0]
            topology_format = "gro"

        return TrajectoryLayout(
            topology_path=topology_path,
            trajectory_paths=trajectory_paths,
            trajectory_format="xtc",
            topology_format=topology_format,
        )
