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
import math
import re
from dataclasses import dataclass, field
from numbers import Real
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterator, Sequence

import numpy as np
from numpy.typing import NDArray
from pydantic import ValidationError

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from polyzymd.config.schema import SimulationConfig
    from polyzymd.engines.base import SimulationEngine, TrajectoryLayout

LOGGER = logging.getLogger(__name__)
_WARNED_GRO_TOPOLOGY_PATHS: set[Path] = set()
_WARNED_ELEMENT_ENRICHMENT_KEYS: set[str] = set()

_ELEMENT_SYMBOLS = frozenset(
    {
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Br",
        "I",
    }
)
_COMMON_ION_ELEMENTS = {
    "NA": "Na",
    "SOD": "Na",
    "CL": "Cl",
    "CLA": "Cl",
    "MG": "Mg",
    "CA": "Ca",
    "CAL": "Ca",
    "ZN": "Zn",
    "K": "K",
    "POT": "K",
    "FE": "Fe",
    "CU": "Cu",
    "MN": "Mn",
}
_BIOMOLECULAR_NAME_PREFIX_ELEMENTS = {"C", "N", "O", "S", "P", "F", "I"}
_STANDARD_BIOMOLECULAR_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "CYX",
    "GLN",
    "GLU",
    "GLY",
    "HID",
    "HIE",
    "HIP",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "A",
    "C",
    "G",
    "T",
    "U",
    "DA",
    "DC",
    "DG",
    "DT",
    "DU",
}


def _normalized_warning_path(path: Path) -> Path:
    """Return a stable key for topology warning de-duplication.

    Parameters
    ----------
    path : Path
        Topology path to normalize.

    Returns
    -------
    Path
        Expanded and resolved path suitable for process-wide warning gating.
    """
    expanded = Path(path).expanduser()
    try:
        return expanded.resolve(strict=False)
    except OSError:
        return expanded.absolute()


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


def _canonical_element_symbol(value: object) -> str | None:
    """Return a canonical element symbol for an unambiguous token.

    Parameters
    ----------
    value : object
        Candidate atom type or atom name token.

    Returns
    -------
    str or None
        Canonical element symbol, or ``None`` when the token is not a plain
        element symbol.
    """

    token = str(value).strip()
    if not token or not token.isalpha() or len(token) > 2:
        return None
    symbol = token[0].upper() + token[1:].lower()
    if symbol in _ELEMENT_SYMBOLS:
        return symbol
    return None


def _context_allows_type_element(symbol: str, name: object, resname: object | None) -> bool:
    """Return whether atom context supports an element inferred from type.

    Parameters
    ----------
    symbol : str
        Element symbol inferred from the atom type.
    name : object
        Atom name from the topology.
    resname : object or None
        Residue name used to reject ambiguous ion aliases.

    Returns
    -------
    bool
        ``True`` when the atom name and residue context support the type-derived
        element.
    """

    name_token = str(name).strip().upper()
    name_symbol = _canonical_element_symbol(name_token)
    residue = str(resname).strip().upper().rstrip("+-") if resname is not None else ""
    if (
        len(name_token) == 2
        and name_symbol == symbol
        and residue not in _STANDARD_BIOMOLECULAR_RESIDUES
        and residue not in _COMMON_ION_ELEMENTS
    ):
        return True

    name_symbol = _infer_element_from_atom_name(name, resname)
    return name_symbol == symbol


def _infer_elements_from_atom_types(universe: Any) -> tuple[list[str] | None, str]:
    """Infer topology elements from atom types when context agrees.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe or compatible test double.

    Returns
    -------
    tuple[list[str] | None, str]
        Inferred element symbols and a diagnostic message. The element list is
        ``None`` when any atom type is missing, ambiguous, or conflicts with
        atom name and residue context.
    """

    try:
        atom_types = list(universe.atoms.types)
    except (AttributeError, TypeError):
        return None, "atom types are unavailable"
    try:
        names = list(universe.atoms.names)
    except (AttributeError, TypeError):
        return None, "atom names are unavailable for type validation"
    try:
        resnames = list(universe.atoms.resnames)
    except (AttributeError, TypeError):
        resnames = [None] * len(names)

    if len(atom_types) != len(names) or len(atom_types) != len(resnames):
        return None, "atom type, name, and residue metadata lengths differ"

    inferred: list[str] = []
    for index, (atom_type, name, resname) in enumerate(
        zip(atom_types, names, resnames, strict=True)
    ):
        symbol = _canonical_element_symbol(atom_type)
        if symbol is None:
            return None, f"atom type {atom_type!r} is not element-like"
        if not _context_allows_type_element(symbol, name, resname):
            return None, (
                f"atom type {atom_type!r} at index {index} conflicts with "
                f"atom name {name!r} and residue {resname!r}"
            )
        inferred.append(symbol)
    return inferred, "all atom types are element-like and context-safe"


def _infer_element_from_atom_name(name: object, resname: object | None = None) -> str | None:
    """Infer an element from a conservative atom-name heuristic.

    Parameters
    ----------
    name : object
        Atom name from the topology.
    resname : object or None, optional
        Residue name used to distinguish ions such as calcium from protein
        alpha carbons named ``CA``.

    Returns
    -------
    str or None
        Inferred element symbol, or ``None`` when no conservative inference is
        possible.
    """

    raw_name = str(name).strip()
    token = re.sub(r"^\d+", "", raw_name).upper()
    if not token:
        return None
    base_token = re.sub(r"\d+$", "", token)
    has_numeric_suffix = base_token != token

    residue = str(resname).strip().upper().rstrip("+-") if resname is not None else ""
    residue_element = _COMMON_ION_ELEMENTS.get(residue)
    if residue_element is not None:
        name_element = _COMMON_ION_ELEMENTS.get(base_token)
        if name_element == residue_element:
            return name_element
        return None
    if base_token in {"CL", "BR"}:
        return _canonical_element_symbol(base_token)
    if (
        len(base_token) == 2
        and base_token == token
        and residue not in _STANDARD_BIOMOLECULAR_RESIDUES
        and _canonical_element_symbol(base_token)
    ):
        return None
    if (
        has_numeric_suffix
        and len(base_token) == 2
        and residue not in _STANDARD_BIOMOLECULAR_RESIDUES
        and _canonical_element_symbol(base_token)
    ):
        return None
    if token in _COMMON_ION_ELEMENTS and token != "CA" and not residue:
        return _COMMON_ION_ELEMENTS[token]
    if token.startswith("H"):
        return "H"
    first_letter = token[0]
    if first_letter in _BIOMOLECULAR_NAME_PREFIX_ELEMENTS:
        return first_letter
    return None


def _infer_elements_from_atom_names(universe: Any) -> tuple[list[str] | None, str]:
    """Infer topology elements from atom names when all atoms are recognized.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe or compatible test double.

    Returns
    -------
    tuple[list[str] | None, str]
        Inferred element symbols and a diagnostic message. The element list is
        ``None`` when any atom cannot be inferred conservatively.
    """

    try:
        names = list(universe.atoms.names)
    except (AttributeError, TypeError):
        return None, "atom names are unavailable"
    try:
        resnames = list(universe.atoms.resnames)
    except (AttributeError, TypeError):
        resnames = [None] * len(names)

    inferred: list[str] = []
    for index, (name, resname) in enumerate(zip(names, resnames, strict=True)):
        symbol = _infer_element_from_atom_name(name, resname)
        if symbol is None:
            return None, f"atom name {name!r} at index {index} is not safely inferable"
        inferred.append(symbol)
    return inferred, "all atom names matched conservative element rules"


def _universe_has_elements(universe: Any) -> bool:
    """Return whether a universe exposes complete element metadata.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe or compatible test double.

    Returns
    -------
    bool
        ``True`` when element metadata can be read for all atoms.
    """

    try:
        elements = list(universe.atoms.elements)
    except (AttributeError, TypeError):
        return False
    return len(elements) == len(universe.atoms)


def enrich_universe_elements(
    universe: Any, *, topology_key: str | Path | None = None
) -> dict[str, Any]:
    """Add missing MDAnalysis element metadata when it can be inferred safely.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe to enrich in place.
    topology_key : str or Path or None, optional
        Stable key for one-time logging, usually the topology path.

    Returns
    -------
    dict[str, Any]
        Lightweight enrichment metadata with ``applied``, ``source``, and
        diagnostic ``message`` or ``reason`` fields.
    """

    if _universe_has_elements(universe):
        metadata = {"applied": False, "source": "existing", "message": "elements already present"}
        universe._polyzymd_element_enrichment = metadata
        return metadata

    elements, type_reason = _infer_elements_from_atom_types(universe)
    source = "types"
    if elements is None:
        elements, name_reason = _infer_elements_from_atom_names(universe)
        source = "names"
        reason = f"type inference skipped: {type_reason}; name inference skipped: {name_reason}"
    else:
        reason = type_reason

    if elements is None:
        metadata = {"applied": False, "source": None, "reason": reason}
        universe._polyzymd_element_enrichment = metadata
        return metadata

    universe.add_TopologyAttr("elements", elements)
    metadata = {
        "applied": True,
        "source": source,
        "message": f"elements inferred from atom {source}",
    }
    universe._polyzymd_element_enrichment = metadata

    warning_key = str(topology_key or id(universe))
    if warning_key not in _WARNED_ELEMENT_ENRICHMENT_KEYS:
        _WARNED_ELEMENT_ENRICHMENT_KEYS.add(warning_key)
        LOGGER.info(
            "Inferred missing topology elements from atom %s for %s",
            source,
            topology_key or "MDAnalysis universe",
        )
    return metadata


def _trajectory_frame_index(trajectory: Any) -> int | None:
    """Return the current trajectory frame index when available.

    Parameters
    ----------
    trajectory : Any
        MDAnalysis trajectory reader or a compatible test double.

    Returns
    -------
    int | None
        Current frame index, or ``None`` when it cannot be determined.
    """

    for owner in (trajectory, getattr(trajectory, "ts", None)):
        frame = getattr(owner, "frame", None)
        if isinstance(frame, bool):
            continue
        try:
            return int(frame)
        except (TypeError, ValueError):
            continue
    return None


def _restore_trajectory_frame(trajectory: Any, frame_index: int | None) -> None:
    """Restore a trajectory reader to a previous frame when possible.

    Parameters
    ----------
    trajectory : Any
        MDAnalysis trajectory reader or a compatible test double.
    frame_index : int | None
        Frame index captured before a metadata probe.
    """

    if frame_index is None:
        return
    try:
        trajectory[frame_index]
    except (AttributeError, IndexError, TypeError, ValueError):
        return


def _finite_numeric_time(value: object) -> float | None:
    """Return a finite real-valued time or ``None``.

    Parameters
    ----------
    value : object
        Candidate time value from MDAnalysis metadata.

    Returns
    -------
    float | None
        Finite time value, or ``None`` when the value is missing, non-real, or
        non-finite.
    """

    if value is None or isinstance(value, bool) or not isinstance(value, Real):
        return None
    time_value = float(value)
    if not math.isfinite(time_value):
        return None
    return time_value


def _trajectory_raw_time(trajectory: Any) -> float | None:
    """Return raw timestamp metadata from a trajectory timestep.

    Parameters
    ----------
    trajectory : Any
        MDAnalysis trajectory reader or a compatible test double.

    Returns
    -------
    float | None
        Raw finite timestep time from ``ts.data['time']``, or ``None`` when it is
        unavailable.
    """

    ts = getattr(trajectory, "ts", None)
    data = getattr(ts, "data", None)
    if not isinstance(data, dict):
        return None
    return _finite_numeric_time(data.get("time"))


def _trajectory_time(trajectory: Any) -> float | None:
    """Return the best finite timestamp exposed by a trajectory reader.

    Parameters
    ----------
    trajectory : Any
        MDAnalysis trajectory reader or a compatible test double.

    Returns
    -------
    float | None
        Raw timestep time when available, otherwise ``trajectory.time`` when it
        is finite.
    """

    raw_time = _trajectory_raw_time(trajectory)
    if raw_time is not None:
        return raw_time
    return _finite_numeric_time(getattr(trajectory, "time", None))


class _TimestampPreservingTrajectory:
    """Proxy that exposes raw MDAnalysis timestep timestamps.

    Some MDAnalysis multi-DCD ``ChainReader`` instances normalize
    ``reader.time`` to a loaded-frame-relative origin even though each timestep
    keeps the absolute source timestamp in ``reader.ts.data['time']``. This proxy
    preserves the reader protocol while making ``time`` return that raw timestamp
    when available.
    """

    def __init__(self, reader: Any) -> None:
        """Store the wrapped trajectory reader.

        Parameters
        ----------
        reader : Any
            MDAnalysis trajectory reader to wrap.
        """

        object.__setattr__(self, "_reader", reader)

    def __len__(self) -> int:
        """Return the wrapped trajectory length."""

        return len(self._reader)

    def __iter__(self) -> Iterator[Any]:
        """Iterate over the wrapped trajectory reader."""

        return iter(self._reader)

    def __getitem__(self, item: Any) -> Any:
        """Delegate frame and slice access to the wrapped reader."""

        return self._reader[item]

    def __getattr__(self, name: str) -> Any:
        """Delegate unknown attributes to the wrapped reader."""

        return getattr(self._reader, name)

    def __setattr__(self, name: str, value: Any) -> None:
        """Delegate mutable reader attributes to the wrapped reader."""

        if name == "_reader":
            object.__setattr__(self, name, value)
            return
        setattr(self._reader, name, value)

    @property
    def time(self) -> float:
        """Return raw timestep time when available."""

        raw_time = _trajectory_raw_time(self._reader)
        if raw_time is not None:
            return raw_time
        return self._reader.time


def _wrap_timestamp_preserving_trajectory(trajectory: Any) -> Any:
    """Wrap trajectory readers that hide raw source timestamps.

    Parameters
    ----------
    trajectory : Any
        MDAnalysis trajectory reader.

    Returns
    -------
    Any
        Original reader when no correction is needed, otherwise a proxy that
        exposes raw timestamp metadata through ``time``.
    """

    previous_frame = _trajectory_frame_index(trajectory)
    try:
        trajectory[0]
        raw_time = _trajectory_raw_time(trajectory)
        reported_time = _finite_numeric_time(getattr(trajectory, "time", None))
    except (AttributeError, IndexError, TypeError, ValueError):
        return trajectory
    finally:
        _restore_trajectory_frame(trajectory, previous_frame)

    if raw_time is None or reported_time is None:
        return trajectory
    if math.isclose(raw_time, reported_time, rel_tol=1e-12, abs_tol=1e-12):
        return trajectory
    return _TimestampPreservingTrajectory(trajectory)


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
    topology_format : str or None, optional
        Engine-reported topology format, when available.
    trajectory_format : str or None, optional
        Engine-reported trajectory format, when available.
    warnings : list[str]
        Discovery warnings that should be preserved in downstream provenance.
    """

    topology_file: Path
    trajectory_files: list[Path] = field(default_factory=list)
    n_segments: int = 0
    working_directory: Path = field(default_factory=Path)
    replicate: int = 1
    topology_format: str | None = None
    trajectory_format: str | None = None
    warnings: list[str] = field(default_factory=list)

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
    call that needs it, so construction remains cheap. Engine resolution errors
    propagate unless an explicit ``engine_override`` is supplied.

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

    # ------------------------------------------------------------------
    # Engine delegation helpers
    # ------------------------------------------------------------------

    def _get_engine(self) -> "SimulationEngine":
        """Lazily create and cache the simulation engine.

        Engine resolution errors from ``create_engine()`` propagate to callers
        unless an explicit ``engine_override`` supplies a valid backend.

        Returns
        -------
        SimulationEngine
            Engine instance resolved from config.
        """
        if self._engine is None:
            from polyzymd.engines import create_engine

            self._engine = create_engine(self.config, override=self._engine_override)
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
            engine_dir = engine.resolve_engine_working_directory(working_dir)
        except (AttributeError, TypeError):
            engine_dir = working_dir
        try:
            layout = engine.resolve_trajectory_layout(engine_dir, replicate)
        except FileNotFoundError:
            raise
        except (TypeError, ValueError, ValidationError) as exc:
            # Invalid path-like inputs are treated as missing layouts so callers
            # can use the existing discovery fallback path
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

    def _gro_topology_warning(self, layout: "TrajectoryLayout") -> str | None:
        """Build the chain-ID warning for GRO topology layouts.

        GRO files do not reliably preserve chain identifiers, which can
        break chain-based selections (``chainid A/B/C``) used by many
        analysis plugins.

        Parameters
        ----------
        layout : TrajectoryLayout
            Resolved layout from the engine.

        Returns
        -------
        str or None
            Actionable warning text when the layout uses a GRO topology,
            otherwise ``None``.
        """
        if layout.topology_path is None or layout.topology_format.lower() != "gro":
            return None
        return (
            f"Using GRO topology {layout.topology_path} — GRO files may not preserve "
            "chain identifiers. Chain-based selections (chainid A/B/C) used by "
            "analysis plugins may be unreliable. Prefer a PDB topology when available."
        )

    def _warn_gro_topology(self, layout: "TrajectoryLayout") -> str | None:
        """Emit a one-time warning when the resolved topology is a GRO file.

        Parameters
        ----------
        layout : TrajectoryLayout
            Resolved layout from the engine.

        Returns
        -------
        str or None
            Actionable warning text when the layout uses a GRO topology,
            otherwise ``None``.
        """
        warning = self._gro_topology_warning(layout)
        if warning is not None and layout.topology_path is not None:
            topology_key = _normalized_warning_path(layout.topology_path)
            if topology_key in _WARNED_GRO_TOPOLOGY_PATHS:
                return warning
            _WARNED_GRO_TOPOLOGY_PATHS.add(topology_key)
            LOGGER.warning(warning)
        return warning

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
            raise FileNotFoundError(
                self._format_missing_data_message(
                    "Working directory not found",
                    working_dir=working_dir,
                    replicate=replicate,
                    available_replicates=available,
                    action=(
                        "Run or complete the simulation for this replicate, or verify "
                        "the config output/scratch paths before rerunning analysis."
                    ),
                )
            )

        # Delegate file discovery to the simulation engine
        layout = self._resolve_layout(working_dir, replicate=replicate)

        if layout.topology_path is None:
            raise FileNotFoundError(f"No topology file found in {working_dir}")
        if not layout.trajectory_paths:
            available = self._find_available_replicates()
            raise FileNotFoundError(
                self._format_missing_data_message(
                    "No production trajectory files found",
                    working_dir=working_dir,
                    replicate=replicate,
                    available_replicates=available,
                    action=(
                        "Run or complete the production simulation for this replicate, "
                        "then rerun analysis. Use --recompute only after trajectory files exist."
                    ),
                )
            )

        warnings = []
        gro_warning = self._gro_topology_warning(layout)
        if gro_warning is not None:
            warnings.append(gro_warning)

        return TrajectoryInfo(
            topology_file=layout.topology_path,
            trajectory_files=layout.trajectory_paths,
            n_segments=len(layout.trajectory_paths),
            working_directory=working_dir,
            replicate=replicate,
            topology_format=layout.topology_format,
            trajectory_format=layout.trajectory_format,
            warnings=warnings,
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
        enrich_universe_elements(u, topology_key=info.topology_file)
        u.trajectory = _wrap_timestamp_preserving_trajectory(u.trajectory)

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
        trajectory = u.trajectory

        previous_frame = _trajectory_frame_index(trajectory)
        try:
            times_ps = []
            for frame_index in range(len(trajectory)):
                trajectory[frame_index]
                time_ps = _trajectory_time(trajectory)
                if time_ps is None:
                    raise ValueError(
                        "Trajectory timestamps are unavailable for frame-time extraction"
                    )
                times_ps.append(time_ps)
        finally:
            _restore_trajectory_frame(trajectory, previous_frame)
        times = np.array(times_ps, dtype=np.float64)

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
        trajectory = u.trajectory

        # Get timestep from trajectory
        if len(trajectory) < 2:
            raise ValueError("Need at least 2 frames to determine timestep")

        previous_frame = _trajectory_frame_index(trajectory)
        try:
            trajectory[0]
            t0 = _trajectory_time(trajectory)
            trajectory[1]
            t1 = _trajectory_time(trajectory)
        finally:
            _restore_trajectory_frame(trajectory, previous_frame)
        if t0 is None or t1 is None:
            raise ValueError("Trajectory timestamps are unavailable for timestep detection")

        dt = t1 - t0  # in ps (MDAnalysis default)

        if unit == "ns":
            dt = dt / 1000.0
        elif unit != "ps":
            raise ValueError(f"Unknown time unit: {unit}")

        return float(dt)

    def get_first_frame_time(self, replicate: int, unit: str = "ps") -> float | None:
        """Return the first loaded frame timestamp when available.

        MDAnalysis reports trajectory times in picoseconds. This method probes
        cached Universe metadata without changing the caller-visible current
        frame when the reader exposes a restorable frame index.

        Parameters
        ----------
        replicate : int
            Replicate number.
        unit : str, optional
            Time unit for output. Options are ``"ps"`` and ``"ns"``, by default
            ``"ps"``.

        Returns
        -------
        float | None
            Finite first-frame timestamp in the requested unit, or ``None`` when
            the trajectory does not expose a usable timestamp.

        Raises
        ------
        ValueError
            Raised when ``unit`` is not ``"ps"`` or ``"ns"``.
        """

        if unit not in {"ps", "ns"}:
            raise ValueError(f"Unknown time unit: {unit}")

        u = self.load_universe(replicate)
        trajectory = u.trajectory
        previous_frame = _trajectory_frame_index(trajectory)
        try:
            trajectory[0]
            time_ps = _trajectory_time(trajectory)
        except (AttributeError, IndexError, TypeError, ValueError):
            return None
        finally:
            _restore_trajectory_frame(trajectory, previous_frame)

        if time_ps is None:
            return None
        if unit == "ns":
            return time_ps / 1000.0
        return time_ps

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
        discovered = self._discover_replicates_from_config()
        if discovered:
            return discovered

        scratch_dir = self._get_scratch_directory()
        if scratch_dir is None:
            return []
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

    def _discover_replicates_from_config(self) -> list[int]:
        """Discover replicate directories through SimulationConfig when available.

        Returns
        -------
        list[int]
            Sorted replicate numbers discovered by the config helper, or an
            empty list when unavailable.
        """
        discover = getattr(self.config, "discover_replicate_dirs", None)
        if not callable(discover):
            return []
        try:
            replicate_dirs = discover()
        except (AttributeError, TypeError, OSError, ValueError):
            return []

        replicates: list[int] = []
        for item in replicate_dirs:
            if isinstance(item, tuple) and item:
                try:
                    replicates.append(int(item[0]))
                    continue
                except (TypeError, ValueError):
                    item = item[-1]
            try:
                path = Path(item)
            except TypeError:
                continue
            match = re.search(r"run_?(\d+)$", path.name)
            if match:
                replicates.append(int(match.group(1)))
        return sorted(set(replicates))

    def _get_scratch_directory(self) -> Path | None:
        """Return configured scratch directory when it can be resolved.

        Returns
        -------
        Path or None
            Effective scratch directory, or ``None`` when config metadata is
            incomplete.
        """
        try:
            scratch_dir = self.config.output.effective_scratch_directory
        except AttributeError:
            return None
        if scratch_dir is None:
            return None
        try:
            return Path(scratch_dir)
        except TypeError:
            return None

    def _format_missing_data_message(
        self,
        headline: str,
        *,
        working_dir: Path,
        replicate: int,
        available_replicates: Sequence[int],
        action: str,
    ) -> str:
        """Build a user-facing message for missing replicate trajectory data.

        Parameters
        ----------
        headline : str
            Stable leading message used by existing tests.
        working_dir : Path
            Expected working directory for the replicate.
        replicate : int
            Requested replicate number.
        available_replicates : sequence of int
            Replicates discovered on disk.
        action : str
            Actionable hint for the user.

        Returns
        -------
        str
            Multi-line diagnostic message.
        """
        scratch_dir = self._get_scratch_directory()
        available = (
            ", ".join(str(rep) for rep in available_replicates)
            if available_replicates
            else "none found"
        )
        scratch_text = str(scratch_dir) if scratch_dir is not None else "not configured"
        return (
            f"{headline}: {working_dir}\n"
            f"Replicate: {replicate}\n"
            f"Expected working directory: {working_dir}\n"
            f"Scratch directory: {scratch_text}\n"
            f"Available replicates: {available}\n"
            f"Action: {action}"
        )

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
            replicate = self._infer_replicate(working_dir)
            raise FileNotFoundError(
                self._format_missing_data_message(
                    "No production trajectory files found",
                    working_dir=working_dir,
                    replicate=replicate,
                    available_replicates=self._find_available_replicates(),
                    action=(
                        "Run or complete the production simulation for this replicate, "
                        "then rerun analysis. Use --recompute only after trajectory files exist."
                    ),
                )
            )
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
