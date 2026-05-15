"""Universe loading and provenance helpers for the MDAnalysis extension layer."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Protocol

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from polyzymd.analyses.shared.loader import TrajectoryInfo
    from polyzymd.config.schema import SimulationConfig

LOGGER = logging.getLogger(__name__)


class _TrajectoryLoaderLike(Protocol):
    """Structural protocol for trajectory loaders used by ``UniverseProvider``."""

    def load_universe(self, replicate: int, cache: bool = True) -> Universe:
        """Load a universe for a replicate.

        Parameters
        ----------
        replicate : int
            Replicate index to load.
        cache : bool, optional
            Whether the loader may reuse its universe cache, by default True.

        Returns
        -------
        Universe
            Loaded MDAnalysis universe.
        """

    def get_trajectory_info(self, replicate: int) -> TrajectoryInfo:
        """Resolve trajectory files for a replicate without loading coordinates.

        Parameters
        ----------
        replicate : int
            Replicate index to inspect.

        Returns
        -------
        TrajectoryInfo
            Resolved trajectory metadata.
        """


LoaderFactory = Callable[..., _TrajectoryLoaderLike]

GRO_CHAIN_ID_WARNING_TEMPLATE = (
    "Using GRO topology {path} — GRO files may not preserve chain identifiers. "
    "Chain-based selections (chainid A/B/C) used by analysis plugins may be unreliable. "
    "Prefer a PDB topology when available."
)


@dataclass(frozen=True)
class FileIdentity:
    """Filesystem identity for an input topology or trajectory file."""

    path: Path
    format: str | None
    size_bytes: int
    mtime_ns: int

    @classmethod
    def from_path(cls, path: Path, file_format: str | None = None) -> FileIdentity:
        """Create file identity metadata from a filesystem path.

        Parameters
        ----------
        path : Path
            File path to identify.
        file_format : str or None, optional
            Format reported by the trajectory layout. When omitted, the file
            suffix is used without the leading dot.

        Returns
        -------
        FileIdentity
            Path, format, size, and modification-time metadata.
        """
        resolved_path = Path(path)
        stat = resolved_path.stat()
        inferred_format = file_format or resolved_path.suffix.removeprefix(".").lower() or None
        return cls(
            path=resolved_path,
            format=inferred_format,
            size_bytes=stat.st_size,
            mtime_ns=stat.st_mtime_ns,
        )

    def as_dict(self) -> dict[str, Any]:
        """Serialize the identity to JSON-compatible primitive values.

        Returns
        -------
        dict[str, Any]
            Dictionary representation with the path converted to a string.
        """
        return {
            "path": str(self.path),
            "format": self.format,
            "size_bytes": self.size_bytes,
            "mtime_ns": self.mtime_ns,
        }


@dataclass(frozen=True)
class UniverseProvenance:
    """Provenance for one replicate universe loaded from PolyzyMD outputs."""

    replicate: int
    working_directory: Path
    topology: FileIdentity
    trajectories: tuple[FileIdentity, ...]
    n_segments: int
    loader_class: str
    config_engine: str | None
    engine_override: str | None = None
    warnings: tuple[str, ...] = field(default_factory=tuple)

    def as_dict(self) -> dict[str, Any]:
        """Serialize provenance to JSON-compatible primitive values.

        Returns
        -------
        dict[str, Any]
            Dictionary representation suitable for manifests and tests.
        """
        return {
            "replicate": self.replicate,
            "working_directory": str(self.working_directory),
            "topology": self.topology.as_dict(),
            "trajectories": [trajectory.as_dict() for trajectory in self.trajectories],
            "n_segments": self.n_segments,
            "loader_class": self.loader_class,
            "config_engine": self.config_engine,
            "engine_override": self.engine_override,
            "warnings": list(self.warnings),
        }


@dataclass
class UniverseProvider:
    """Config-aware provider for MDAnalysis universes and input provenance."""

    config: SimulationConfig
    engine_override: str | None = None
    loader: _TrajectoryLoaderLike | None = None
    loader_factory: LoaderFactory | None = None
    _provenance_cache: dict[int, UniverseProvenance] = field(default_factory=dict, init=False)
    _warned_gro_topologies: set[Path] = field(default_factory=set, init=False)

    def __post_init__(self) -> None:
        """Validate loader injection settings after dataclass construction."""
        if self.loader is not None and self.loader_factory is not None:
            raise ValueError("Provide either loader or loader_factory, not both.")

    @classmethod
    def from_config(cls, config: SimulationConfig, **kwargs: Any) -> UniverseProvider:
        """Create a universe provider from a simulation configuration.

        Parameters
        ----------
        config : SimulationConfig
            PolyzyMD simulation configuration.
        **kwargs : Any
            Optional provider settings such as ``engine_override``, ``loader``,
            or ``loader_factory``.

        Returns
        -------
        UniverseProvider
            Provider that lazily instantiates the trajectory loader.
        """
        return cls(config=config, **kwargs)

    def load_universe(self, replicate: int, *, cache: bool = True) -> Universe:
        """Load an MDAnalysis universe for a replicate through the existing loader.

        Parameters
        ----------
        replicate : int
            Replicate index to load.
        cache : bool, optional
            Whether the underlying loader may reuse its universe cache, by
            default True.

        Returns
        -------
        Universe
            Loaded MDAnalysis universe from the underlying trajectory loader.
        """
        self.provenance_for(replicate, refresh=not cache)
        return self._get_loader().load_universe(replicate, cache=cache)

    def provenance_for(self, replicate: int, *, refresh: bool = False) -> UniverseProvenance:
        """Return provenance for a replicate, computing it when needed.

        Parameters
        ----------
        replicate : int
            Replicate index to inspect.
        refresh : bool, optional
            Recompute provenance even when cached, by default False.

        Returns
        -------
        UniverseProvenance
            Input file identity and loader metadata for the replicate.
        """
        if not refresh and replicate in self._provenance_cache:
            return self._provenance_cache[replicate]

        loader = self._get_loader()
        info = loader.get_trajectory_info(replicate)
        provenance = self._build_provenance(info=info, loader=loader)
        self._provenance_cache[replicate] = provenance
        return provenance

    def get_provenance(self, replicate: int) -> UniverseProvenance | None:
        """Return cached provenance without triggering trajectory discovery.

        Parameters
        ----------
        replicate : int
            Replicate index whose cached provenance should be returned.

        Returns
        -------
        UniverseProvenance or None
            Cached provenance when available, otherwise ``None``.
        """
        return self._provenance_cache.get(replicate)

    def _get_loader(self) -> _TrajectoryLoaderLike:
        """Return the lazily instantiated trajectory loader.

        Returns
        -------
        _TrajectoryLoaderLike
            Injected or default trajectory loader.
        """
        if self.loader is None:
            factory = self.loader_factory or self._default_loader_factory
            self.loader = factory(self.config, engine_override=self.engine_override)
        return self.loader

    @staticmethod
    def _default_loader_factory(
        config: SimulationConfig,
        *,
        engine_override: str | None = None,
    ) -> _TrajectoryLoaderLike:
        """Create the default shared trajectory loader lazily.

        Parameters
        ----------
        config : SimulationConfig
            PolyzyMD simulation configuration.
        engine_override : str or None, optional
            Engine override passed through to ``TrajectoryLoader``.

        Returns
        -------
        _TrajectoryLoaderLike
            Shared trajectory loader instance.
        """
        from polyzymd.analyses.shared.loader import TrajectoryLoader

        return TrajectoryLoader(config, engine_override=engine_override)

    def _build_provenance(
        self,
        *,
        info: TrajectoryInfo,
        loader: _TrajectoryLoaderLike,
    ) -> UniverseProvenance:
        """Build provenance metadata from shared loader trajectory info.

        Parameters
        ----------
        info : TrajectoryInfo
            Resolved trajectory metadata from the shared loader.
        loader : _TrajectoryLoaderLike
            Loader instance used to resolve the metadata.

        Returns
        -------
        UniverseProvenance
            Provenance with file identities and warnings.
        """
        topology_format = self._metadata_format(info, "topology_format", info.topology_file)
        trajectory_format = self._metadata_format(info, "trajectory_format", None)
        warnings = list(getattr(info, "warnings", []))
        gro_warning = self._gro_chain_id_warning(info, topology_format)
        if gro_warning is not None:
            if gro_warning not in warnings:
                warnings.append(gro_warning)
                self._log_gro_warning_once(info.topology_file, gro_warning)

        return UniverseProvenance(
            replicate=info.replicate,
            working_directory=Path(info.working_directory),
            topology=FileIdentity.from_path(info.topology_file, topology_format),
            trajectories=tuple(
                FileIdentity.from_path(path, trajectory_format) for path in info.trajectory_files
            ),
            n_segments=info.n_segments,
            loader_class=type(loader).__name__,
            config_engine=self._config_engine(),
            engine_override=self.engine_override,
            warnings=tuple(warnings),
        )

    def _config_engine(self) -> str | None:
        """Return the configured simulation engine name when it is concrete.

        Returns
        -------
        str or None
            String engine name from the config, otherwise ``None``.
        """
        engine = getattr(self.config, "engine", None)
        if isinstance(engine, str):
            return engine
        return None

    @staticmethod
    def _metadata_format(info: TrajectoryInfo, field_name: str, path: Path | None) -> str | None:
        """Return a format value from trajectory info with suffix fallback.

        Parameters
        ----------
        info : TrajectoryInfo
            Resolved trajectory metadata.
        field_name : str
            Name of the optional format field on ``TrajectoryInfo``.
        path : Path or None
            File path used to infer the format when metadata is absent.

        Returns
        -------
        str or None
            Lowercase format string or ``None`` when unavailable.
        """
        value = getattr(info, field_name, None)
        if isinstance(value, str) and value:
            return value.lower()
        if path is not None:
            suffix = Path(path).suffix.removeprefix(".").lower()
            if suffix:
                return suffix
        return None

    def _gro_chain_id_warning(
        self, info: TrajectoryInfo, topology_format: str | None
    ) -> str | None:
        """Return an actionable GRO chain-ID warning when applicable.

        Parameters
        ----------
        info : TrajectoryInfo
            Resolved trajectory metadata.
        topology_format : str or None
            Resolved topology format.

        Returns
        -------
        str or None
            Warning text for GRO topology inputs, otherwise ``None``.
        """
        topology_path = Path(info.topology_file)
        suffix_is_gro = topology_path.suffix.lower() == ".gro"
        format_is_gro = topology_format == "gro"
        if not suffix_is_gro and not format_is_gro:
            return None
        return GRO_CHAIN_ID_WARNING_TEMPLATE.format(path=topology_path)

    def _log_gro_warning_once(self, topology_path: Path, warning: str) -> None:
        """Log a GRO chain-ID warning once for each topology path.

        Parameters
        ----------
        topology_path : Path
            GRO topology path that triggered the warning.
        warning : str
            Warning text to log.
        """
        path = Path(topology_path)
        if path in self._warned_gro_topologies:
            return
        self._warned_gro_topologies.add(path)
        LOGGER.warning(warning)
