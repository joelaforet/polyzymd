"""Trajectory loading utilities for PolyzyMD analysis.

This module provides config-aware trajectory loading that understands
PolyzyMD's directory structure and daisy-chain continuation patterns.

Key Features
------------
- Config-based path resolution (config.yaml is single source of truth)
- Automatic detection of daisy-chain trajectory segments
- Support for both scratch and projects directories
- Lazy loading and memory-efficient iteration
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Iterator, Sequence

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig

LOGGER = logging.getLogger(__name__)

# MDAnalysis is optional - only import when needed
try:
    import MDAnalysis as mda
    from MDAnalysis.core.universe import Universe

    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
    Universe = None  # type: ignore


def _require_mdanalysis() -> None:
    """Raise ImportError if MDAnalysis is not available."""
    if not HAS_MDANALYSIS:
        raise ImportError(
            "MDAnalysis is required for trajectory analysis.\n"
            "Install with: pip install polyzymd[analysis]"
        )


@dataclass
class TrajectoryInfo:
    """Information about discovered trajectory files.

    Attributes
    ----------
    topology_file : Path
        Path to topology file (PDB)
    trajectory_files : list[Path]
        List of trajectory files (DCD) in order
    n_segments : int
        Number of daisy-chain segments
    working_directory : Path
        Base working directory for this replicate
    replicate : int
        Replicate number
    """

    topology_file: Path
    trajectory_files: list[Path] = field(default_factory=list)
    n_segments: int = 0
    working_directory: Path = field(default_factory=Path)
    replicate: int = 1

    @property
    def n_trajectory_files(self) -> int:
        """Number of trajectory files found."""
        return len(self.trajectory_files)

    def validate(self) -> None:
        """Validate that all files exist."""
        if not self.topology_file.exists():
            raise FileNotFoundError(f"Topology not found: {self.topology_file}")

        missing = [f for f in self.trajectory_files if not f.exists()]
        if missing:
            raise FileNotFoundError(f"Missing trajectory files: {missing}")


class TrajectoryLoader:
    """Config-aware trajectory loader for PolyzyMD simulations.

    This class handles the complexity of finding and loading trajectories
    from PolyzyMD's output structure, including:
    - Daisy-chain continuation segments
    - Scratch vs projects directory resolution
    - Multiple replicates

    Parameters
    ----------
    config : SimulationConfig
        PolyzyMD simulation configuration

    Examples
    --------
    >>> from polyzymd.config import load_config
    >>> config = load_config("config.yaml")
    >>> loader = TrajectoryLoader(config)
    >>>
    >>> # Load single replicate
    >>> u = loader.load_universe(replicate=1)
    >>> print(f"Loaded {len(u.trajectory)} frames")
    >>>
    >>> # Get trajectory info without loading
    >>> info = loader.get_trajectory_info(replicate=1)
    >>> print(f"Found {info.n_segments} segments")
    >>>
    >>> # Load multiple replicates
    >>> for rep in range(1, 6):
    ...     u = loader.load_universe(replicate=rep)
    ...     # ... analyze

    Notes
    -----
    Frame indices in MDAnalysis are 0-indexed. For user-facing output,
    add 1 to follow PyMOL convention (1-indexed frames).
    """

    def __init__(self, config: "SimulationConfig") -> None:
        _require_mdanalysis()
        self.config = config
        self._universe_cache: dict[int, "Universe"] = {}

    def get_trajectory_info(self, replicate: int) -> TrajectoryInfo:
        """Get trajectory file information for a replicate.

        Parameters
        ----------
        replicate : int
            Replicate number (1-indexed)

        Returns
        -------
        TrajectoryInfo
            Information about discovered trajectory files

        Raises
        ------
        FileNotFoundError
            If working directory or required files don't exist
        """
        # Get working directory from config
        working_dir = self.config.get_working_directory(replicate)

        if not working_dir.exists():
            available = self._find_available_replicates()
            available_str = ", ".join(str(r) for r in available) if available else "none found"
            raise FileNotFoundError(
                f"Working directory not found: {working_dir}\n"
                f"Has replicate {replicate} been simulated?\n"
                f"Available replicates: {available_str}"
            )

        # Find topology and trajectories
        # Methods handle both new (production_N/) and legacy (production/) structures
        topology_file = self._find_topology(working_dir)
        trajectory_files = self._find_trajectories(working_dir)

        return TrajectoryInfo(
            topology_file=topology_file,
            trajectory_files=trajectory_files,
            n_segments=len(trajectory_files),
            working_directory=working_dir,
            replicate=replicate,
        )

    def load_universe(
        self,
        replicate: int,
        cache: bool = True,
    ) -> "Universe":
        """Load MDAnalysis Universe for a replicate.

        Parameters
        ----------
        replicate : int
            Replicate number (1-indexed)
        cache : bool, optional
            If True (default), cache the Universe for reuse

        Returns
        -------
        Universe
            MDAnalysis Universe with trajectory loaded

        Notes
        -----
        For daisy-chain trajectories, all segments are loaded as a
        continuous trajectory using MDAnalysis's ChainReader.
        """
        _require_mdanalysis()

        if cache and replicate in self._universe_cache:
            return self._universe_cache[replicate]

        info = self.get_trajectory_info(replicate)
        info.validate()

        # Load universe - MDAnalysis handles multiple trajectory files
        if len(info.trajectory_files) == 1:
            u = mda.Universe(
                str(info.topology_file),
                str(info.trajectory_files[0]),
            )
        else:
            # Multiple segments - use ChainReader
            u = mda.Universe(
                str(info.topology_file),
                [str(f) for f in info.trajectory_files],
            )

        if cache:
            self._universe_cache[replicate] = u

        return u

    def iter_replicates(
        self,
        replicates: Sequence[int],
    ) -> Iterator[tuple[int, "Universe"]]:
        """Iterate over multiple replicates.

        Parameters
        ----------
        replicates : sequence of int
            Replicate numbers to load

        Yields
        ------
        tuple of (int, Universe)
            Replicate number and loaded Universe

        Examples
        --------
        >>> for rep, u in loader.iter_replicates([1, 2, 3, 4, 5]):
        ...     rmsf = compute_rmsf(u)
        ...     results[rep] = rmsf
        """
        for rep in replicates:
            yield rep, self.load_universe(rep)

    def get_frame_times(
        self,
        replicate: int,
        unit: str = "ns",
    ) -> NDArray[np.float64]:
        """Get time values for each frame.

        Parameters
        ----------
        replicate : int
            Replicate number
        unit : str, optional
            Time unit for output. Options: "ps", "ns". Default is "ns".

        Returns
        -------
        NDArray[np.float64]
            Array of time values for each frame
        """
        u = self.load_universe(replicate)

        # Get times from trajectory
        times = np.array([ts.time for ts in u.trajectory], dtype=np.float64)

        # Convert units (MDAnalysis uses ps internally)
        if unit == "ns":
            times = times / 1000.0
        elif unit != "ps":
            raise ValueError(f"Unknown time unit: {unit}. Use 'ps' or 'ns'.")

        return times

    def get_timestep(self, replicate: int, unit: str = "ps") -> float:
        """Get the trajectory timestep (time between frames).

        Parameters
        ----------
        replicate : int
            Replicate number
        unit : str, optional
            Time unit. Options: "ps", "ns". Default is "ps".

        Returns
        -------
        float
            Time between consecutive frames
        """
        u = self.load_universe(replicate)

        # Get timestep from trajectory
        if len(u.trajectory) < 2:
            raise ValueError("Need at least 2 frames to determine timestep")

        u.trajectory[0]
        t0 = u.trajectory.time
        u.trajectory[1]
        t1 = u.trajectory.time

        dt = t1 - t0  # in ps (MDAnalysis default)

        if unit == "ns":
            dt = dt / 1000.0
        elif unit != "ps":
            raise ValueError(f"Unknown time unit: {unit}")

        return float(dt)

    def clear_cache(self) -> None:
        """Clear the Universe cache to free memory."""
        self._universe_cache.clear()

    def _find_available_replicates(self) -> list[int]:
        """Find available replicate numbers from existing run_* directories.

        Returns
        -------
        list[int]
            Sorted list of replicate numbers that have simulation directories
        """
        scratch_dir = self.config.output.effective_scratch_directory
        if not scratch_dir.exists():
            return []

        replicates = []
        for d in scratch_dir.iterdir():
            if d.is_dir() and d.name.startswith("run_"):
                try:
                    # Extract replicate number from directory name (e.g., "run_1" -> 1)
                    rep_num = int(d.name.split("_")[1])
                    replicates.append(rep_num)
                except (IndexError, ValueError):
                    continue
        return sorted(replicates)

    def _find_topology(self, working_dir: Path) -> Path:
        """Find topology file in working directory.

        Search order:
        1. solvated_system.pdb (primary - always exists after build)
        2. production_0/production_0_topology.pdb (daisy-chain structure)
        3. production/production_topology.pdb (legacy pre-daisy-chain)
        4. Glob fallback for any topology PDB
        """
        # Primary: solvated_system.pdb in working_dir root
        topology = working_dir / "solvated_system.pdb"
        if topology.exists():
            return topology

        # Daisy-chain structure: production_0/production_0_topology.pdb
        topology = working_dir / "production_0" / "production_0_topology.pdb"
        if topology.exists():
            return topology

        # Legacy structure: production/production_topology.pdb
        topology = working_dir / "production" / "production_topology.pdb"
        if topology.exists():
            return topology

        # Glob fallback: search for any topology PDB
        for pattern in [
            "production_*/*_topology.pdb",  # daisy-chain
            "production/*_topology.pdb",  # legacy with segments (unlikely)
            "*.pdb",  # any PDB in root
        ]:
            pdbs = sorted(working_dir.glob(pattern))
            if pdbs:
                return pdbs[0]

        raise FileNotFoundError(f"No topology file found in {working_dir}")

    def _find_trajectories(self, working_dir: Path) -> list[Path]:
        """Find trajectory files, handling daisy-chain segments.

        Search order:
        1. production_N/production_N_trajectory.dcd (daisy-chain, multiple segments)
        2. production/production_trajectory.dcd (legacy single file, deprecated)
        3. Glob fallback for any production DCD files
        """
        # New daisy-chain structure: production_N/production_N_trajectory.dcd
        # Use OS-agnostic pattern matching
        pattern = re.compile(r"production_(\d+)[/\\]production_\d+_trajectory\.dcd$")

        segments: dict[int, Path] = {}
        for f in working_dir.glob("production_*/production_*_trajectory.dcd"):
            match = pattern.search(str(f))
            if match:
                idx = int(match.group(1))
                segments[idx] = f

        if segments:
            # Return in segment order
            return [segments[i] for i in sorted(segments.keys())]

        # Legacy structure: production/production_trajectory.dcd (single file, no segments)
        legacy_traj = working_dir / "production" / "production_trajectory.dcd"
        if legacy_traj.exists():
            LOGGER.warning(
                f"Using deprecated trajectory structure: {legacy_traj}. "
                "This format (production/production_trajectory.dcd) is deprecated. "
                "Re-run simulation with current PolyzyMD to use production_N/ structure."
            )
            return [legacy_traj]

        # Last resort: any production DCD files (excluding equilibration)
        dcds = sorted(working_dir.glob("**/production*trajectory.dcd"))
        if dcds:
            return dcds

        raise FileNotFoundError(f"No production trajectory files found in {working_dir}")


def parse_time_string(time_str: str) -> tuple[float, str]:
    """Parse a time string with units into value and unit.

    Parameters
    ----------
    time_str : str
        Time string like "100ns", "5000ps", "100 ns", etc.

    Returns
    -------
    tuple of (float, str)
        Numeric value and unit string

    Examples
    --------
    >>> parse_time_string("100ns")
    (100.0, "ns")
    >>> parse_time_string("5000 ps")
    (5000.0, "ps")
    >>> parse_time_string("100")  # Default to ns
    (100.0, "ns")
    """
    time_str = time_str.strip()

    # Try to extract number and unit
    match = re.match(r"^([\d.]+)\s*([a-zA-Z]*)$", time_str)
    if not match:
        raise ValueError(f"Cannot parse time string: {time_str}")

    value = float(match.group(1))
    unit = match.group(2).lower() if match.group(2) else "ns"

    if unit not in ("ns", "ps", "fs"):
        raise ValueError(f"Unknown time unit: {unit}. Use 'ns', 'ps', or 'fs'.")

    return value, unit


def convert_time(value: float, from_unit: str, to_unit: str) -> float:
    """Convert time between units.

    Parameters
    ----------
    value : float
        Time value
    from_unit : str
        Source unit ("fs", "ps", "ns")
    to_unit : str
        Target unit ("fs", "ps", "ns")

    Returns
    -------
    float
        Converted time value
    """
    # Convert to picoseconds first
    to_ps = {"fs": 0.001, "ps": 1.0, "ns": 1000.0}
    from_ps = {"fs": 1000.0, "ps": 1.0, "ns": 0.001}

    if from_unit not in to_ps or to_unit not in from_ps:
        raise ValueError(f"Unknown unit: {from_unit} or {to_unit}")

    ps_value = value * to_ps[from_unit]
    return ps_value * from_ps[to_unit]


def time_to_frame(
    time: float,
    time_unit: str,
    timestep: float,
    timestep_unit: str = "ps",
) -> int:
    """Convert time to frame index.

    Parameters
    ----------
    time : float
        Time value
    time_unit : str
        Unit of time value
    timestep : float
        Time between frames
    timestep_unit : str
        Unit of timestep (default: "ps")

    Returns
    -------
    int
        Frame index (0-indexed)
    """
    # Convert both to same units
    time_ps = convert_time(time, time_unit, "ps")
    dt_ps = convert_time(timestep, timestep_unit, "ps")

    return int(time_ps / dt_ps)
