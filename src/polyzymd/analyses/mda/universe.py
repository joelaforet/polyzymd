"""Universe loading and provenance helpers for the MDAnalysis extension layer."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Protocol

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from polyzymd.analyses.shared.loader import TrajectoryInfo
    from polyzymd.config.schema import SimulationConfig


class _TrajectoryLoaderLike(Protocol):
    """Structural protocol for trajectory loaders used by ``UniverseProvider``."""

    def load_universe(self, replicate: int, cache: bool = True) -> Universe:
        """Load a universe for a replicate.

        Loaders that support periodic boundary policies also accept a
        ``pbc_policy`` keyword. The provider passes it only when a policy other
        than ``"as_is"`` is requested, so loaders without the keyword keep
        working.

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


def trajectory_variant(paths: "Sequence[Path]") -> str | None:
    """Name which GROMACS trajectory variant the engine picked, if any.

    The engine prefers the whole-molecule centered file when the run wrote one,
    and that changes what the coordinates mean, so the choice belongs in
    provenance rather than only in a filename.
    """
    names = [Path(path).name.lower() for path in paths]
    if not names:
        return None
    if any("centered" in name for name in names):
        return "centered"
    if any("nojump" in name for name in names):
        return "nojump"
    return "raw"


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
    """Provenance for one replicate universe loaded from PolyzyMD outputs.

    Attributes
    ----------
    pbc_policy : str
        Periodic boundary policy applied on load, ``"as_is"`` or
        ``"make_whole"``.
    topology_has_bonds : bool or None
        Whether the loaded topology carries bonds. ``None`` before a universe
        has been loaded.
    bond_source : str
        Where the bonds came from: ``"conect"``, ``"guessed"``, or ``"none"``.
    trajectory_variant : str or None
        Which trajectory the engine chose. GROMACS writes post-processed
        trajectories, so this is ``"centered"`` for ``prod_centered.xtc``,
        ``"nojump"`` for ``prod_nojump.xtc``, and ``"raw"`` otherwise. OpenMM
        segments are always ``"raw"``.
    """

    replicate: int
    working_directory: Path
    topology: FileIdentity
    trajectories: tuple[FileIdentity, ...]
    n_segments: int
    loader_class: str
    config_engine: str | None
    engine_override: str | None = None
    warnings: tuple[str, ...] = field(default_factory=tuple)
    excluded_segments: tuple[int, ...] = field(default_factory=tuple)
    segment_status: tuple[tuple[int, str], ...] = field(default_factory=tuple)
    pbc_policy: str = "as_is"
    topology_has_bonds: bool | None = None
    bond_source: str = "none"
    trajectory_variant: str | None = None

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
            "excluded_segments": list(self.excluded_segments),
            "segment_status": {str(index): status for index, status in self.segment_status},
            "pbc_policy": self.pbc_policy,
            "topology_has_bonds": self.topology_has_bonds,
            "bond_source": self.bond_source,
            "trajectory_variant": self.trajectory_variant,
        }


@dataclass
class UniverseProvider:
    """Config-aware provider for MDAnalysis universes and input provenance."""

    config: SimulationConfig
    engine_override: str | None = None
    require_complete: bool = True
    loader: _TrajectoryLoaderLike | None = None
    loader_factory: LoaderFactory | None = None
    pbc_policy: str = "as_is"
    _provenance_cache: dict[int, UniverseProvenance] = field(default_factory=dict, init=False)

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

    def load_universe(
        self,
        replicate: int,
        *,
        cache: bool = True,
        pbc_policy: str | None = None,
    ) -> Universe:
        """Load an MDAnalysis universe for a replicate through the existing loader.

        Parameters
        ----------
        replicate : int
            Replicate index to load.
        cache : bool, optional
            Whether the underlying loader may reuse its universe cache, by
            default True.
        pbc_policy : str or None, optional
            Periodic boundary policy for this call, overriding the provider
            setting. ``"as_is"`` (the default) leaves coordinates untouched;
            ``"make_whole"`` unwraps the protein and polymer selection and
            requires a topology with bonds.

        Returns
        -------
        Universe
            Loaded MDAnalysis universe from the underlying trajectory loader.
        """
        policy = self.pbc_policy if pbc_policy is None else str(pbc_policy)
        self.provenance_for(replicate, refresh=not cache)
        universe = self._get_loader().load_universe(
            replicate,
            cache=cache,
            **self._segment_kwargs(),
            **self._pbc_kwargs(policy),
        )
        self._record_universe_facts(replicate, universe, policy)
        return universe

    def _pbc_kwargs(self, policy: str) -> dict[str, Any]:
        """Forward the policy only when it differs from the loader default."""

        return {} if policy == "as_is" else {"pbc_policy": policy}

    def _record_universe_facts(self, replicate: int, universe: Any, policy: str) -> None:
        """Add what loading established to the cached provenance."""

        from polyzymd.analyses.shared.topology import topology_bond_source

        provenance = self._provenance_cache.get(replicate)
        if provenance is None:
            return
        has_bonds, bond_source = topology_bond_source(universe)
        self._provenance_cache[replicate] = replace(
            provenance,
            pbc_policy=policy,
            topology_has_bonds=has_bonds,
            bond_source=bond_source,
        )

    def provenance_for(self, replicate: int, *, refresh: bool = False) -> UniverseProvenance:
        """Return provenance for a replicate, computing it when needed.

        Parameters
        ----------
        replicate : int
            Replicate index to inspect.
        refresh : bool, optional
            Recompute provenance even when cached, by default False. Facts
            established by loading the universe, such as the bond source, are
            carried across a refresh because rediscovering the input files does
            not re-examine the topology.

        Returns
        -------
        UniverseProvenance
            Input file identity and loader metadata for the replicate.
        """
        previous = self._provenance_cache.get(replicate)
        if not refresh and previous is not None:
            return previous

        loader = self._get_loader()
        info = loader.get_trajectory_info(replicate, **self._segment_kwargs())
        provenance = self._build_provenance(info=info, loader=loader)
        if previous is not None and previous.topology_has_bonds is not None:
            # Discovery metadata can be refreshed without reloading the
            # universe, so keep what the last load established about it.
            provenance = replace(
                provenance,
                pbc_policy=previous.pbc_policy,
                topology_has_bonds=previous.topology_has_bonds,
                bond_source=previous.bond_source,
            )
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

    def _segment_kwargs(self) -> dict[str, Any]:
        """Return the segment-completeness keyword to forward to the loader.

        The keyword is forwarded only when it differs from the loader default,
        so loaders that predate the option keep working.

        Returns
        -------
        dict[str, Any]
            Either an empty mapping or ``{"require_complete": False}``.
        """

        return {} if self.require_complete else {"require_complete": False}

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
            excluded_segments=tuple(getattr(info, "excluded_segments", ()) or ()),
            segment_status=tuple(sorted((getattr(info, "segment_status", None) or {}).items())),
            pbc_policy=self.pbc_policy,
            trajectory_variant=trajectory_variant(info.trajectory_files),
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
