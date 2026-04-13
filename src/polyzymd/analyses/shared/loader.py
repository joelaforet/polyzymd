"""Trajectory loading utilities for PolyzyMD analysis.

This module provides config-aware trajectory loading that understands
PolyzyMD's directory structure and daisy-chain continuation patterns.
File discovery is delegated to the active simulation engine so that
both OpenMM and GROMACS directory layouts are handled transparently.

Key Features
------------
- Config-based path resolution (config.yaml is single source of truth)
- Engine-aware file discovery (OpenMM daisy-chain, GROMACS flat layout)
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
    from MDAnalysis.core.universe import Universe

    from polyzymd.config.schema import SimulationConfig
    from polyzymd.engines.base import SimulationEngine, TrajectoryLayout

LOGGER = logging.getLogger(__name__)


def _require_mdanalysis(feature_name: str = "trajectory analysis") -> None:
    """Raise ImportError if MDAnalysis is not available."""
    try:
        import MDAnalysis  # noqa: F401
    except ImportError:
        raise ImportError(
            f"MDAnalysis is required for {feature_name}.\n"
            "Ensure MDAnalysis is available in the PolyzyMD pixi environment "
            '(for example: pixi run -e build python -c "import MDAnalysis")'
        ) from None


def _require_matplotlib(feature_name: str = "plotting") -> None:
    """Raise ImportError if matplotlib is not available."""
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        raise ImportError(
            f"matplotlib is required for {feature_name}.\n"
            "Ensure matplotlib is available in the PolyzyMD pixi environment "
            '(for example: pixi run -e build python -c "import matplotlib")'
        ) from None


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

    - Daisy-chain continuation segments (OpenMM)
    - Flat production directories (GROMACS)
    - Scratch vs projects directory resolution
    - Multiple replicates

    File discovery is delegated to the simulation engine resolved from the
    config's ``engine`` field.  The engine is created lazily on the first
    call that needs it, so construction remains cheap.

    Parameters
    ----------
    config : SimulationConfig
        PolyzyMD simulation configuration.
    engine_override : str or None, optional
        Force a specific engine name (``"openmm"`` or ``"gromacs"``)
        instead of reading ``config.engine``.

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
    >>>
    >>> # Explicit engine override for GROMACS directories
    >>> loader = TrajectoryLoader(config, engine_override="gromacs")

    Notes
    -----
    Frame indices in MDAnalysis are 0-indexed. For user-facing output,
    add 1 to follow PyMOL convention (1-indexed frames).
    """

    def __init__(
        self,
        config: "SimulationConfig",
        engine_override: str | None = None,
    ) -> None:
        _require_mdanalysis()
        self.config = config
        self._engine_override = engine_override
        self._engine: SimulationEngine | None = None
        self._universe_cache: dict[int, "Universe"] = {}
        self._warned_gro_topologies: set[Path] = set()

    # ------------------------------------------------------------------
    # Engine delegation helpers
    # ------------------------------------------------------------------

    def _get_engine(self) -> "SimulationEngine":
        """Lazily create and cache the simulation engine.

        Falls back to OpenMM when the config's ``engine`` field is
        unrecognised (e.g. a mock object in tests).

        Returns
        -------
        SimulationEngine
            Engine instance resolved from config.
        """
        if self._engine is None:
            from polyzymd.engines import create_engine

            try:
                self._engine = create_engine(self.config, override=self._engine_override)
            except (ValueError, TypeError):
                # Unrecognised engine name (e.g. MagicMock in tests) —
                # fall back to OpenMM which works with any directory layout.
                LOGGER.debug(
                    "Could not resolve engine from config (%s); "
                    "falling back to OpenMM layout resolver.",
                    getattr(self.config, "engine", "<no engine attr>"),
                )
                from polyzymd.engines.openmm import OpenMMEngine

                self._engine = OpenMMEngine.from_config(self.config)
        return self._engine

    def _resolve_layout(
        self,
        working_dir: Path,
        replicate: int | None = None,
    ) -> "TrajectoryLayout":
        """Resolve trajectory layout via the engine.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int or None, optional
            Replicate index.  When ``None`` (e.g. from
            ``find_topology(working_dir)``), the replicate is inferred
            from the directory name (``run_<N>``) with a fallback to 1.

        Returns
        -------
        TrajectoryLayout
            Engine-resolved file layout.

        Raises
        ------
        FileNotFoundError
            If the engine cannot resolve the layout (e.g. invalid paths).
        """
        if replicate is None:
            replicate = self._infer_replicate(working_dir)
        engine = self._get_engine()
        try:
            layout = engine.resolve_trajectory_layout(working_dir, replicate)
        except Exception as exc:
            # Pydantic ValidationError, TypeError, etc. when paths are
            # invalid (e.g. MagicMock in tests).  Translate to
            # FileNotFoundError so callers' existing handlers work.
            if isinstance(exc, FileNotFoundError):
                raise
            raise FileNotFoundError(
                f"Engine could not resolve trajectory layout in {working_dir}: {exc}"
            ) from exc
        self._warn_gro_topology(layout)
        return layout

    @staticmethod
    def _infer_replicate(working_dir: Path) -> int:
        """Best-effort replicate number from a ``run_<N>`` directory name.

        Parameters
        ----------
        working_dir : Path
            Directory whose name may encode the replicate index.

        Returns
        -------
        int
            Parsed replicate number or ``1`` as a safe fallback.
        """
        try:
            dir_name = str(Path(working_dir).name)
        except (TypeError, ValueError):
            return 1
        match = re.match(r"run_(\d+)", dir_name)
        if match:
            return int(match.group(1))
        return 1

    def _warn_gro_topology(self, layout: "TrajectoryLayout") -> None:
        """Emit a one-time warning when the resolved topology is a GRO file.

        GRO files do not reliably preserve chain identifiers, which can
        break chain-based selections (``chainid A/B/C``) used by many
        analysis plugins.

        Parameters
        ----------
        layout : TrajectoryLayout
            Resolved layout from the engine.
        """
        if (
            layout.topology_path is not None
            and layout.topology_format.lower() == "gro"
            and layout.topology_path not in self._warned_gro_topologies
        ):
            self._warned_gro_topologies.add(layout.topology_path)
            LOGGER.warning(
                "Using GRO topology %s — GRO files may not preserve chain "
                "identifiers.  Chain-based selections (chainid A/B/C) used "
                "by analysis plugins may be unreliable.  Prefer a PDB "
                "topology when available.",
                layout.topology_path,
            )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

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

        # Delegate file discovery to the simulation engine
        layout = self._resolve_layout(working_dir, replicate=replicate)

        if layout.topology_path is None:
            raise FileNotFoundError(f"No topology file found in {working_dir}")
        if not layout.trajectory_paths:
            raise FileNotFoundError(f"No production trajectory files found in {working_dir}")

        return TrajectoryInfo(
            topology_file=layout.topology_path,
            trajectory_files=layout.trajectory_paths,
            n_segments=len(layout.trajectory_paths),
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
        import MDAnalysis as mda

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

    def find_topology(self, working_dir: Path) -> Path:
        """Find topology file in working directory.

        Delegates file discovery to the simulation engine.  The engine
        applies its own search order (e.g. PDB preference for GROMACS,
        ``solvated_system.pdb`` preference for OpenMM).

        This method is used by several plugins that pass an explicit
        ``working_dir`` unrelated to the current replicate.  The
        replicate index is inferred from the directory name when
        possible (``run_<N>``), falling back to ``1``.

        Parameters
        ----------
        working_dir : Path
            Directory to search for topology files.

        Returns
        -------
        Path
            Path to the topology file.

        Raises
        ------
        FileNotFoundError
            If no topology file is found.
        """
        layout = self._resolve_layout(working_dir, replicate=None)
        if layout.topology_path is None:
            raise FileNotFoundError(f"No topology file found in {working_dir}")
        return layout.topology_path

    def _find_trajectories(self, working_dir: Path) -> list[Path]:
        """Find trajectory files via the simulation engine.

        Parameters
        ----------
        working_dir : Path
            Working directory to search.

        Returns
        -------
        list[Path]
            Ordered trajectory file paths.

        Raises
        ------
        FileNotFoundError
            If no trajectory files are found.
        """
        layout = self._resolve_layout(working_dir, replicate=None)
        if not layout.trajectory_paths:
            raise FileNotFoundError(f"No production trajectory files found in {working_dir}")
        return layout.trajectory_paths


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
