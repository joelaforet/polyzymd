"""Load the replicates of a set of simulation conditions as MDAnalysis universes.

:class:`Study` maps condition labels to simulation ``config.yaml`` files.
``study[label]`` gives a :class:`Condition`, whose replicates are the run
directories found on disk. Each :class:`Replicate` loads its production
trajectory through :class:`~polyzymd.analyses.mda.universe.UniverseProvider`
and gives the frame indices and times left after the equilibration window,
found by :func:`~polyzymd.analyses.shared.window.resolve_replicate_trajectory_window`.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from polyzymd.analyses.exceptions import ProtocolError

if TYPE_CHECKING:
    import numpy as np

    from polyzymd.analyses.shared.window import TrajectoryWindow


class Replicate:
    """One replicate of one condition, with its universe and production frames.

    Parameters
    ----------
    condition : Condition
        The condition this replicate belongs to.
    index : int
        Replicate number, as used in the run directory name.
    """

    def __init__(self, condition: Condition, index: int) -> None:
        self.condition = condition
        self.index = int(index)
        self._universe: Any = None
        self._window: TrajectoryWindow | None = None

    def __repr__(self) -> str:
        return f"Replicate(condition={self.condition.label!r}, index={self.index})"

    def universe(self) -> Any:
        """Load the production trajectory once and return the same universe after.

        Returns
        -------
        MDAnalysis.Universe
            Every production frame in segment order, including the
            equilibration window. Use :attr:`frames` to skip that window.
        """
        if self._universe is None:
            self._universe = self.condition._provider.load_universe(self.index)
        return self._universe

    def _production_window(self) -> TrajectoryWindow:
        """Measure the equilibration window on the loaded trajectory once."""
        if self._window is None:
            from polyzymd.analyses.shared.window import resolve_replicate_trajectory_window

            try:
                self._window = resolve_replicate_trajectory_window(
                    loader=self.condition._provider._get_loader(),
                    replicate=self.index,
                    equilibration=self.condition.equilibration,
                    n_frames_total=len(self.universe().trajectory),
                )
            except ValueError as exc:
                raise ProtocolError(
                    f"Cannot remove the equilibration window {self.condition.equilibration!r} "
                    f"from {self!r}: {exc}",
                    hint="Pass a shorter equilibration, or leave the replicate out with "
                    "replicates=[...].",
                ) from exc
        return self._window

    @property
    def frames(self) -> np.ndarray:
        """Indices of the production frames after the equilibration window.

        The first index is the first frame whose time is at or after the
        equilibration time and the last is the trajectory's last frame.
        """
        import numpy as np

        window = self._production_window()
        return np.arange(window.start, window.stop, window.step, dtype=np.int64)

    @property
    def times(self) -> np.ndarray:
        """Simulation time of each entry of :attr:`frames`, in ns.

        Each time is the first frame's timestamp, or zero when the trajectory
        has none, plus the frame index times the frame interval.
        """
        window = self._production_window()
        origin_ps = window.first_frame_time_ps or 0.0
        return (origin_ps + self.frames * window.timestep_ps) / 1000.0

    @property
    def identity(self) -> dict[str, Any]:
        """The config hash, equilibration window and input files of this replicate.

        Returns
        -------
        dict
            ``config_hash``, ``equilibration``, and the ``topology`` and
            ``trajectories`` file records (path, format, size and modification
            time) from :class:`~polyzymd.analyses.mda.universe.FileIdentity`.
        """
        provenance = self.condition._provider.provenance_for(self.index, refresh=True)
        return {
            "config_hash": self.condition.config_hash,
            "equilibration": self.condition.equilibration,
            "topology": provenance.topology.as_dict(),
            "trajectories": [item.as_dict() for item in provenance.trajectories],
        }


class Condition:
    """One simulation condition and the replicates chosen for it.

    Parameters
    ----------
    label : str
        Condition label.
    config_path : Path
        Absolute path of the simulation ``config.yaml``.
    equilibration : str
        Window removed from the start of every replicate, for example ``"10ns"``.
    replicates : sequence of int, optional
        Replicate numbers to use. Defaults to every run directory on disk.
    """

    def __init__(
        self,
        label: str,
        config_path: Path,
        equilibration: str,
        replicates: Sequence[int] | None = None,
    ) -> None:
        from polyzymd.analyses._framework.cache_identity import compute_config_hash
        from polyzymd.analyses.mda.universe import UniverseProvider
        from polyzymd.config.schema import SimulationConfig

        self.label = label
        self.config_path = config_path
        self.equilibration = equilibration
        try:
            self.config = SimulationConfig.from_yaml(config_path)
            found = sorted(int(index) for index, _ in self.config.discover_replicate_dirs())
        except (OSError, ValueError) as exc:
            raise ProtocolError(
                f"Condition {label!r}: cannot read {config_path}: {exc}",
                hint="Point the condition at a PolyzyMD simulation config.yaml.",
            ) from exc
        chosen = found if replicates is None else sorted({int(index) for index in replicates})
        missing = sorted(set(chosen) - set(found))
        if not chosen or missing:
            raise ProtocolError(
                f"Condition {label!r}: replicates {missing or chosen} have no run directory "
                f"under the scratch directory of {config_path}; found {found}.",
                hint="Pass replicates that exist on disk, or run the simulations first.",
            )
        self.config_hash = compute_config_hash(self.config)
        self._provider = UniverseProvider.from_config(self.config)
        self.replicates = [Replicate(self, index) for index in chosen]

    def __repr__(self) -> str:
        return f"Condition(label={self.label!r}, replicates={[r.index for r in self.replicates]})"


class Study:
    """A set of simulation conditions loaded as per-replicate universes.

    Build one with :meth:`from_configs`. Iterating yields the conditions in
    order, and the first one is the default control.
    """

    def __init__(self, conditions: Sequence[Condition]) -> None:
        self._conditions = {condition.label: condition for condition in conditions}

    @classmethod
    def from_configs(
        cls,
        configs: Mapping[str, str | Path] | Sequence[str | Path],
        *,
        equilibration: str,
        replicates: Sequence[int] | None = None,
    ) -> Study:
        """Build a study from simulation config paths.

        Parameters
        ----------
        configs : mapping of str to path, or sequence of paths
            Condition label to ``config.yaml``, control first. A sequence
            labels each condition by the name of the folder holding its config.
        equilibration : str
            Window removed from the start of every replicate's production
            trajectory, for example ``"100ns"``.
        replicates : sequence of int, optional
            Replicate numbers for every condition. Defaults to the run
            directories found on disk for each condition.

        Returns
        -------
        Study
            The study with every condition's replicates found.

        Raises
        ------
        ProtocolError
            If no config is given, a config is missing or unreadable, or a
            requested replicate has no run directory.
        """
        from polyzymd.analyses.protocols import _labels
        from polyzymd.analyses.shared.loader import parse_time_string

        if isinstance(configs, Mapping):
            labels = [str(label) for label in configs]
            paths = [Path(path).expanduser().resolve() for path in configs.values()]
        else:
            paths = [Path(path).expanduser().resolve() for path in configs]
            labels = _labels(paths, None)
        missing = [str(path) for path in paths if not path.is_file()]
        if not paths or missing:
            raise ProtocolError(
                f"Config file(s) not found: {', '.join(missing) or 'none given'}.",
                hint="Pass at least one simulation config.yaml; the first one is the control.",
            )
        try:
            parse_time_string(str(equilibration))
        except ValueError as exc:
            raise ProtocolError(
                f"Cannot read the equilibration window {equilibration!r}: {exc}",
                hint="Write it as a time such as '100ns', '500ps' or '0ns'.",
            ) from exc
        return cls(
            [
                Condition(label, path, str(equilibration), replicates)
                for label, path in zip(labels, paths, strict=True)
            ]
        )

    def __getitem__(self, label: str) -> Condition:
        if label not in self._conditions:
            raise ProtocolError(
                f"The study has no condition {label!r}.", hint=f"Use one of {self.labels}."
            )
        return self._conditions[label]

    def __iter__(self) -> Iterator[Condition]:
        return iter(self._conditions.values())

    def __len__(self) -> int:
        return len(self._conditions)

    def __repr__(self) -> str:
        return f"Study(labels={self.labels})"

    @property
    def labels(self) -> list[str]:
        """Condition labels in order, control first."""
        return list(self._conditions)

    @property
    def control(self) -> str:
        """Label of the first condition, the default control of a comparison."""
        return self.labels[0]

    def timeseries(self, function: Any, *args: Any, **kwargs: Any) -> Any:
        """Measure ``function`` on every production frame of every replicate.

        See :func:`polyzymd.analyses.timeseries.run_timeseries` for the
        arguments and what is stored.
        """
        from polyzymd.analyses.timeseries import run_timeseries

        return run_timeseries(self, function, *args, **kwargs)
