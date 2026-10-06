"""Load the replicates of a set of simulation conditions as MDAnalysis universes.

:class:`Study` maps condition labels to simulation ``config.yaml`` files.
``study[label]`` gives a :class:`Condition`, whose replicates are the run
directories found on disk. Each :class:`Replicate` loads its production
trajectory through :class:`~polyzymd.analyses.universe.UniverseProvider`
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

    from polyzymd.analyses.shared.loader import SegmentJoin
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
            equilibration window and any segment's last frame that the next
            segment records again. Use :attr:`frames` to skip both.
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

    def _segment_join(self) -> SegmentJoin | None:
        """Return the boundary repairs found when the universe was loaded, if any."""
        self.universe()
        get_join = getattr(self.condition._provider._get_loader(), "segment_join", None)
        join = get_join(self.index) if callable(get_join) else None
        return join if join is not None and join.repaired else None

    @property
    def frames(self) -> np.ndarray:
        """Indices of the production frames after the equilibration window.

        The first index is the first frame whose time is at or after the
        equilibration time and the last is the trajectory's last frame. When
        a segment's last frame is recorded again at the same time by the next
        segment (see :class:`~polyzymd.analyses.shared.loader.SegmentJoin`),
        that frame's index is left out. With the condition's ``stride`` above
        1, every ``stride``-th of these frames is kept, starting with the
        first. With the condition's ``until_ns``, frames after that time are
        left out, so conditions of unequal length can be compared over a
        common window.
        """
        frames = self._production_frames()
        if self.condition.until_ns is not None:
            frames = frames[self._times_of(frames) <= self.condition.until_ns + 1e-9]
        return frames

    @property
    def production_ns(self) -> float:
        """Time of the last production frame, in ns: how far the replicate was simulated."""
        window = self._production_window()
        last = window.n_frames_total - 1
        join = self._segment_join()
        if join is not None:
            return float(join.times_ps[last] / 1000.0)
        return float(((window.first_frame_time_ps or 0.0) + last * window.timestep_ps) / 1000.0)

    def _production_frames(self) -> np.ndarray:
        """The frames after the equilibration window, with the stride, before ``until_ns``."""
        import numpy as np

        from polyzymd.analyses.shared.window import FRAME_BOUNDARY_TOLERANCE

        window = self._production_window()
        join = self._segment_join()
        if join is None:
            frames = np.arange(window.start, window.stop, window.step, dtype=np.int64)
        else:
            cutoff_ps = window.equilibration_ps - FRAME_BOUNDARY_TOLERANCE * window.timestep_ps
            keep = join.times_ps >= cutoff_ps
            keep[list(join.dropped_frames)] = False
            frames = np.flatnonzero(keep)[:: window.step].astype(np.int64)
        return frames[:: self.condition.stride]

    @property
    def times(self) -> np.ndarray:
        """Simulation time of each entry of :attr:`frames`, in ns.

        Each time is the first frame's timestamp, or zero when the trajectory
        has none, plus the frame index times the frame interval. When a
        segment boundary was repaired, each time is instead the segment's
        first-frame time plus the frame's index within its segment times the
        frame interval, so frames after a dropped or missing frame keep their
        recorded times.
        """
        return self._times_of(self.frames)

    def _times_of(self, frames: np.ndarray) -> np.ndarray:
        """Simulation time of each of ``frames``, in ns; see :attr:`times`."""
        join = self._segment_join()
        if join is not None:
            return join.times_ps[frames] / 1000.0
        window = self._production_window()
        origin_ps = window.first_frame_time_ps or 0.0
        return (origin_ps + frames * window.timestep_ps) / 1000.0

    @property
    def identity(self) -> dict[str, Any]:
        """The config hash, equilibration window and input files of this replicate.

        Returns
        -------
        dict
            ``config_hash``, ``equilibration``, ``stride``, and the ``topology`` and
            ``trajectories`` file records: each file's path relative to the
            folder holding the condition's runs, its format, size and SHA-256
            (:func:`~polyzymd.analyses.shared.file_hashes.file_sha256`, from
            ``progress.json`` when the run recorded it).
        """
        from polyzymd.analyses.shared.file_hashes import file_sha256

        provenance = self.condition._provider.provenance_for(self.index, refresh=True)
        data_root = Path(self.condition.config.output.effective_scratch_directory).resolve()
        recorded = self.condition._provider.recorded_trajectory_hashes(self.index) or {}

        def located(item: Any) -> dict[str, Any]:
            # A record names its files by content and by their place among the
            # condition's runs, never by where this machine keeps them or when
            # they were copied, so a moved or downloaded study keeps its results.
            path = Path(item.path).resolve()
            return {
                "path": str(path.relative_to(data_root))
                if path.is_relative_to(data_root)
                else path.name,
                "format": item.format,
                "size_bytes": item.size_bytes,
                "sha256": file_sha256(path, recorded.get(path)),
            }

        identity = {
            "config_hash": self.condition.config_hash,
            "equilibration": self.condition.equilibration,
            "stride": self.condition.stride,
            "topology": located(provenance.topology),
            "trajectories": [located(item) for item in provenance.trajectories],
        }
        if self.condition.until_ns is not None:
            # Present only with a common window, so other stored records stay valid.
            identity["until_ns"] = self.condition.until_ns
        return identity


#: ``until`` value that ends every replicate at the shortest replicate's last time.
COMMON_UNTIL = "common"


def _parse_until(label: str, until: str | None) -> float | None:
    """Return the common window end ``until`` in ns, or ``None`` (also for ``common``, set later)."""
    if until is None or until == COMMON_UNTIL:
        return None
    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    try:
        value, unit = parse_time_string(str(until))
        return float(convert_time(value, unit, "ns"))
    except ValueError as exc:
        raise ProtocolError(
            f"Condition {label!r}: cannot read until {until!r}: {exc}",
            hint="Write it as a time such as '38ns'.",
        ) from exc


def resolve_common_until(conditions: Sequence[Condition]) -> None:
    """Replace ``until: common`` with the earliest last production time of any replicate.

    Every replicate of every condition then ends at the same time, so the
    replicates share their time points. The time is written back as each
    condition's ``until`` and ``until_ns``, so records name it.
    """
    if not any(condition.until == COMMON_UNTIL for condition in conditions):
        return
    ends = [r.production_ns for condition in conditions for r in condition.replicates]
    if not ends:
        raise ProtocolError(
            "until: common found no replicate to take the common end from.",
            hint="Check that the conditions' runs are where the configs or data.local.yaml say.",
        )
    end = f"{min(ends):.6g}ns"
    for condition in conditions:
        condition.until = end
        condition.until_ns = _parse_until(condition.label, end)


#: Relative difference in production length above which conditions are flagged.
PRODUCTION_LENGTH_TOLERANCE = 0.10


def production_length_warnings(study: Any, labels: Sequence[str] | None = None) -> list[str]:
    """Warn when the conditions were analysed over production of different lengths.

    A difference between conditions can then come from simulated time rather
    than from the conditions: a structure that drifts late appears only in the
    longer runs. Each condition's analysed span is the range of its
    replicates' last analysed times; conditions whose longest and shortest
    spans across the study differ by more than
    :data:`PRODUCTION_LENGTH_TOLERANCE` are named, with the remedy. This is a
    warning only; no value changes.
    """
    ends: dict[str, list[float]] = {}
    for label in list(labels) if labels is not None else list(study.labels):
        condition = study[label]
        spans = []
        for replicate in getattr(condition, "replicates", []):
            times = getattr(replicate, "times", None)
            if times is not None and len(times):
                spans.append(float(times[-1]))
        if spans:
            ends[label] = spans
    if len(ends) < 2:
        return []
    longest = max(max(v) for v in ends.values())
    shortest = min(min(v) for v in ends.values())
    if longest <= 0 or (longest - shortest) / longest <= PRODUCTION_LENGTH_TOLERANCE:
        return []

    def differ(values: list[float]) -> bool:
        return (max(values) - min(values)) / longest > PRODUCTION_LENGTH_TOLERANCE

    described = "; ".join(
        f"{label} {min(v):.4g}" + (f"-{max(v):.4g}" if max(v) != min(v) else "") + " ns"
        for label, v in ends.items()
    )
    remedy = (
        f"Compare over a common window with until {shortest:.4g}ns (--until, or until: in "
        "study.yaml)"
    )
    # Conditions differ when their shortest or their longest replicates do;
    # otherwise only replicates within each condition differ, the same way.
    if differ([min(v) for v in ends.values()]) or differ([max(v) for v in ends.values()]):
        return [
            f"the conditions were analysed up to different times ({described}); a difference "
            f"may come from simulated time rather than the condition. {remedy}"
        ]
    return [
        f"the replicates within each condition were analysed up to different times "
        f"({described}), the same range in every condition; a replicate's value may depend on "
        f"how long it ran. {remedy}"
    ]


def with_data_dir(config: Any, data_dir: Path | None) -> Any:
    """Return ``config`` with its run directories looked for under ``data_dir``.

    ``None`` returns ``config`` unchanged. Only the scratch directory, where
    the run directories are, changes; the config's own file is not edited.
    """
    if data_dir is None:
        return config
    output = config.output.model_copy(update={"scratch_directory": Path(data_dir)})
    return config.model_copy(update={"output": output})


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
    stride : int, optional
        Keep every ``stride``-th production frame of every replicate, 1 by
        default.
    until : str, optional
        End of a common analysis window, such as ``"38ns"``: production frames
        after it are left out, so conditions simulated for different lengths
        are compared over the same time.
    data_dir : Path, optional
        Where this machine keeps the condition's run directories, in place of
        the config's ``scratch_directory``, as ``data.local.yaml`` or
        ``--data`` give it. The config hash is that of the config as written,
        so moving the data does not change it.
    """

    def __init__(
        self,
        label: str,
        config_path: Path,
        equilibration: str,
        replicates: Sequence[int] | None = None,
        stride: int = 1,
        data_dir: Path | None = None,
        until: str | None = None,
    ) -> None:
        from polyzymd.analyses.identity import compute_config_hash
        from polyzymd.analyses.universe import UniverseProvider
        from polyzymd.config.schema import SimulationConfig

        self.label = label
        self.config_path = config_path
        self.equilibration = equilibration
        if isinstance(stride, bool) or not isinstance(stride, int) or stride < 1:
            raise ProtocolError(
                f"Condition {label!r}: stride must be a whole number of at least 1, got {stride!r}.",
                hint="Pass stride=5 to keep every fifth production frame.",
            )
        self.stride = stride
        self.data_dir = None if data_dir is None else Path(data_dir).expanduser().resolve()
        self.until = until
        self.until_ns = _parse_until(label, until)
        try:
            written = SimulationConfig.from_yaml(config_path)
        except (OSError, ValueError) as exc:
            raise ProtocolError(
                f"Condition {label!r}: cannot read {config_path}: {exc}",
                hint="Point the condition at a PolyzyMD simulation config.yaml.",
            ) from exc
        # The hash is of the config as written: where the data sits now is a
        # pointer, not part of what was simulated.
        self.config_hash = compute_config_hash(written)
        self.config = with_data_dir(written, self.data_dir)
        try:
            found = sorted(int(index) for index, _ in self.config.discover_replicate_dirs())
        except (OSError, ValueError) as exc:
            raise ProtocolError(
                f"Condition {label!r}: cannot read {config_path}: {exc}",
                hint="Point the condition at a PolyzyMD simulation config.yaml.",
            ) from exc
        chosen = found if replicates is None else sorted({int(index) for index in replicates})
        missing = sorted(set(chosen) - set(found))
        if not chosen or missing:
            where = self.config.output.effective_scratch_directory
            raise ProtocolError(
                f"Condition {label!r}: replicates {missing or chosen} have no run directory "
                f"under {where}; found {found}.",
                hint="Pass replicates that exist on disk, run the simulations first, or say "
                "where the runs are with data.local.yaml, polyzymd study locate DIR or --data.",
            )
        self._provider = UniverseProvider.from_config(self.config)
        self.replicates = [Replicate(self, index) for index in chosen]

    def __repr__(self) -> str:
        return f"Condition(label={self.label!r}, replicates={[r.index for r in self.replicates]})"


class Study:
    """A set of simulation conditions loaded as per-replicate universes.

    ``Study("study.yaml")`` (or the folder holding it) reads a study file:
    its conditions, equilibration window, stride, replicates and analysis
    settings (:mod:`polyzymd.analyses.study_file`). :meth:`from_configs`
    builds one from config paths instead. Iterating yields the conditions in
    order, and the first one is the default control.
    """

    def __init__(self, conditions: Sequence[Condition] | str | Path) -> None:
        self.protocol: Any = None
        self._built: dict[str, Condition] | None = None
        if isinstance(conditions, (str, Path)):
            from polyzymd.analyses.study_file import load_study_file

            # Conditions are built on first use, so settings() and results()
            # work on a study folder whose trajectories are not on this machine.
            self.protocol = load_study_file(conditions)
        else:
            self._built = {condition.label: condition for condition in conditions}

    @property
    def _conditions(self) -> dict[str, Condition]:
        if self._built is None:
            protocol = self.protocol
            self._built = {
                label: Condition(
                    label,
                    path,
                    protocol.equilibration,
                    protocol.replicates,
                    protocol.stride,
                    protocol.data.get(label),
                    protocol.until,
                )
                for label, path in protocol.conditions.items()
            }
            resolve_common_until(list(self._built.values()))
        return self._built

    @property
    def root(self) -> Path | None:
        """The study folder, when the study was read from a ``study.yaml``."""
        return None if self.protocol is None else self.protocol.root

    def settings(self, run: str) -> dict[str, Any]:
        """Return the settings ``study.yaml`` gives the analysis run ``run``.

        Raises
        ------
        ProtocolError
            If the study has no study file or no such run.
        """
        return dict(self._entry(run).settings)

    def results_dir(self, run: str) -> Path:
        """Return the folder of the run's results: ``<study>/results/<run>``."""
        self._entry(run)
        return self.protocol.results_dir(run)

    def results(self, run: str, *, folder: str | Path | None = None) -> Any:
        """Return the stored per-replicate values of ``run``, without loading any trajectory.

        See :func:`polyzymd.analyses.results.read_results`. This reads
        ``<study>/results/<run>``, where ``polyzymd analyze RUN --study``
        stores them, or ``folder`` when the run was given ``--output-dir``.
        """
        from polyzymd.analyses.results import read_results

        if folder is not None:
            return read_results(Path(folder).expanduser())
        default = self.results_dir(run)
        try:
            return read_results(default)
        except ProtocolError as exc:
            raise ProtocolError(
                str(exc),
                hint=f"Run polyzymd analyze {run} --study {self.protocol.path} first. Results "
                "written with --output-dir are read with Study.results(run, folder=<that folder>).",
            ) from exc

    def _entry(self, run: str) -> Any:
        if self.protocol is None:
            raise ProtocolError(
                "This study was built from config paths, so it has no analyses or results folder.",
                hint="Build it from a study file: pz.Study('study.yaml').",
            )
        if run not in self.protocol.analyses:
            raise ProtocolError(
                f"{self.protocol.path} lists no analysis run {run!r}.",
                hint=f"Use one of {', '.join(self.protocol.analyses) or 'none (add one under analyses:)'}.",
            )
        return self.protocol.analyses[run]

    @classmethod
    def from_configs(
        cls,
        configs: Mapping[str, str | Path] | Sequence[str | Path],
        *,
        equilibration: str,
        replicates: Sequence[int] | None = None,
        stride: int = 1,
        data: Mapping[str, str | Path] | None = None,
        until: str | None = None,
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
        stride : int, optional
            Keep every ``stride``-th production frame of every replicate,
            starting with the first after the equilibration window, 1 by
            default. Every measurement, reference and time then uses those
            frames only; a ``frame`` reference counts them from 1.
        data : mapping of str to path, optional
            Condition label to the directory holding its run directories on
            this machine, in place of its config's ``scratch_directory``. The
            key ``"*"`` applies to every condition the mapping does not name.
        until : str, optional
            End of a common analysis window for every condition, such as
            ``"38ns"``; see :class:`Condition`.

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
        conditions = [
            Condition(
                label,
                path,
                str(equilibration),
                replicates,
                stride,
                (data or {}).get(label, (data or {}).get("*")),
                until,
            )
            for label, path in zip(labels, paths, strict=True)
        ]
        resolve_common_until(conditions)
        return cls(conditions)

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
        if self._built is None:
            return list(self.protocol.conditions)
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

    def per_replicate(self, function: Any, *args: Any, **kwargs: Any) -> Any:
        """Compute one value, or one labelled array, per replicate with ``function``.

        See :func:`polyzymd.analyses.timeseries.run_per_replicate` for the
        arguments and what is stored.
        """
        from polyzymd.analyses.timeseries import run_per_replicate

        return run_per_replicate(self, function, *args, **kwargs)
