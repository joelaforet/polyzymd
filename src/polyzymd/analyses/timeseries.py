"""Run a per-frame function on every replicate of a study and compare conditions.

:func:`run_timeseries`, reached as ``study.timeseries(...)``, runs MDAnalysis
``AnalysisFromFunction`` on each replicate's production frames and stores each
replicate's series with a JSON record of the code, arguments and inputs that
produced it. :meth:`Timeseries.reduce` turns each series into one value per
replicate, and :class:`ReplicateValues` summarises and compares those values
with the replicate as the sampling unit, using
:mod:`polyzymd.analyses.shared.statistics` and
:mod:`polyzymd.analyses.shared.inferential_statistics`: Student t intervals,
Welch or Student t tests, Benjamini-Hochberg adjusted p values and Cohen's d,
with the pymbar statistical inefficiency of each series as a diagnostic.

References
----------
Grossfield, A., Patrone, P. N., Roe, D. R., Schultz, A. J., Siderius, D. W.,
    and Zuckerman, D. M. (2018). Best practices for quantifying sampling
    quality and uncertainty in molecular simulations. Living Journal of
    Computational Molecular Science 1:5067. doi:10.33011/livecoms.1.1.5067
Chodera, J. D., Swope, W. C., Pitera, J. W., Seok, C., and Dill, K. A. (2007).
    Use of the weighted histogram analysis method for the analysis of
    simulated and parallel tempering simulations. Journal of Chemical Theory
    and Computation 3:26-41. doi:10.1021/ct0502864
Welch, B. L. (1947). The generalization of Student's problem when several
    different population variances are involved. Biometrika 34:28-35.
    doi:10.1093/biomet/34.1-2.28
Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate:
    a practical and powerful approach to multiple testing. Journal of the Royal
    Statistical Society Series B 57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
"""

from __future__ import annotations

import hashlib
import inspect
import json
import math
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Sequence

from polyzymd.analyses.exceptions import ProtocolError

if TYPE_CHECKING:
    import numpy as np

    from polyzymd.analyses.protocols import ProtocolReport
    from polyzymd.analyses.study import Replicate, Study

RESULTS_DIR = "polyzymd_results"

#: A replicate is warned about when pymbar's detected start of the
#: equilibrated region falls later than this fraction of its production frames.
#: detect_equilibration picks the start that maximises the effective sample
#: size, and on a stationary series that maximum is flat, so the start moves
#: into the series by chance. On 200 stationary AR(1) series of 2000 frames
#: the start passed 5 percent in 5 to 10 percent of series with 100 or more
#: effective samples and passed 10 percent in 1 to 7 percent of them, while a
#: relaxation that decays over the first 10 percent was found. With about 10
#: effective samples the start passed 10 percent in 45 percent of series, so
#: there the warning says as much about the series length as about the window.
EQUILIBRATION_WARNING_FRACTION = 0.10


@dataclass(frozen=True)
class Select:
    """An MDAnalysis selection string, built into an ``AtomGroup`` per replicate."""

    selection: str


@dataclass(frozen=True)
class UniverseArgument:
    """Stands for the replicate's ``Universe`` in the arguments of a function."""


def select(selection: str) -> Select:
    """Stand for ``universe.select_atoms(selection)`` of each replicate.

    Parameters
    ----------
    selection : str
        MDAnalysis selection string. It is recorded with the result.

    Returns
    -------
    Select
        Placeholder that :func:`run_timeseries` replaces with an ``AtomGroup``.
    """
    return Select(str(selection))


def universe() -> UniverseArgument:
    """Stand for each replicate's ``Universe`` in the arguments of a function.

    Returns
    -------
    UniverseArgument
        Placeholder that :func:`run_timeseries` replaces with the ``Universe``.
    """
    return UniverseArgument()


def _function_record(function: Callable) -> dict[str, Any]:
    """Name a function and hash its source, or its bytecode when no source exists."""
    try:
        code, basis = inspect.getsource(function).encode(), "source"
    except (OSError, TypeError):
        warnings.warn(
            f"No source file for {function!r}, so its bytecode is hashed instead; a changed "
            "function can keep the same hash, and the stored record identifies it less reliably.",
            stacklevel=4,
        )
        body = getattr(function, "__code__", None)
        code = repr((body.co_code, body.co_consts) if body else function).encode()
        basis = "bytecode"
    return {
        "qualname": getattr(function, "__qualname__", repr(function)),
        "module": getattr(function, "__module__", None),
        "hash": hashlib.sha256(code).hexdigest(),
        "hash_of": basis,
    }


def _file_record(path: str | Path) -> dict[str, str]:
    """Record a file by its absolute path and the SHA-256 hash of its content."""
    path = Path(path).expanduser().resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _argument_record(value: Any) -> Any:
    """Describe one argument for the record; values JSON cannot hold are recorded by repr.

    A string or path naming an existing file is recorded with the SHA-256 hash
    of its content, so a changed file changes the record.
    """
    from polyzymd.analyses.reference import Reference

    if isinstance(value, Select):
        return {"select": value.selection}
    if isinstance(value, UniverseArgument):
        return {"universe": True}
    if isinstance(value, Reference):
        record = {key: getattr(value, key) for key in ("mode", "selection", "alignment", "frame")}
        return {"reference": {**record, "file": value.file and _file_record(value.file)}}
    if isinstance(value, (str, Path)):
        try:
            if Path(value).expanduser().is_file():
                return {"file": _file_record(value)}
        except (OSError, ValueError):
            pass
    try:
        return json.loads(json.dumps(value))
    except (TypeError, ValueError):
        return {"repr": repr(value)}


def _build(value: Any, universe_: Any) -> Any:
    """Replace a placeholder with the replicate's AtomGroup or Universe."""
    if isinstance(value, Select):
        atoms = universe_.select_atoms(value.selection)
        if len(atoms) == 0:
            raise ProtocolError(
                f"Selection {value.selection!r} matched no atoms.",
                hint="Check the selection string against the topology.",
            )
        return atoms
    return universe_ if isinstance(value, UniverseArgument) else value


def _build_arguments(
    arguments: dict[int | str, Any], replicate: Replicate
) -> tuple[dict[int | str, Any], dict[str, Any]]:
    """Build every placeholder for one replicate, with the frames that references chose."""
    from polyzymd.analyses.reference import Reference, build_reference

    u, built, chosen = replicate.universe(), {}, {}
    for where, value in arguments.items():
        if isinstance(value, Reference):
            built[where], chosen[str(where)] = build_reference(value, u, replicate.frames)
        else:
            built[where] = _build(value, u)
    return built, chosen


def _safe(text: str) -> str:
    """Turn a label into a folder name."""
    return re.sub(r"[^\w.+-]+", "_", text).strip("_") or "_"


def _versions() -> dict[str, str | None]:
    """Record the versions that do not decide whether a stored result is reused."""
    import MDAnalysis
    import numpy

    from polyzymd import __version__

    return {
        "polyzymd": __version__,
        "MDAnalysis": MDAnalysis.__version__,
        "numpy": numpy.__version__,
        "python": sys.version.split()[0],
    }


@dataclass
class ReplicateSeries:
    """One replicate's per-frame values with their frame indices and times in ns."""

    condition: str
    replicate: int
    values: np.ndarray
    frames: np.ndarray
    times: np.ndarray
    path: Path


def run_timeseries(
    study: Study,
    function: Callable,
    *args: Any,
    unit: str | None,
    name: str | None = None,
    recompute: bool = False,
    output_dir: str | Path | None = None,
    **kwargs: Any,
) -> Timeseries:
    """Measure ``function`` on every production frame of every replicate.

    Each replicate runs ``AnalysisFromFunction(function, *args, **kwargs)``
    with ``frames=replicate.frames``. Its values, frames and times go to
    ``series.npz`` and its record to ``record.json`` in
    ``<output_dir>/polyzymd_results/<name>/<condition>/replicate_<n>/``. The
    record holds the function's name, module and source hash, the arguments
    with their selection strings, the config hash, the input file records,
    the equilibration window, the frames, the times, the unit and the
    PolyzyMD, MDAnalysis, NumPy and Python versions. A stored series is read
    back instead of measured when every field of its record except the
    versions equals the new one. A
    :func:`~polyzymd.analyses.reference.reference` argument is built once per
    replicate before its frames are measured, and the production frame it
    chose, if any, is stored under ``chosen``, which is not compared either.

    Parameters
    ----------
    study : Study
        The study to measure.
    function : callable
        Called once per frame with ``*args`` and ``**kwargs``; returns a number.
    *args
        Arguments of ``function``. :func:`select` becomes the replicate's
        ``AtomGroup``, :func:`universe` its ``Universe`` and
        :func:`~polyzymd.analyses.reference.reference` the reference atoms.
    unit : str or None
        Unit of the returned number; ``None`` for a dimensionless quantity.
    name : str, optional
        Result name used for the folder and the report. Defaults to the
        function's ``__name__``.
    recompute : bool, optional
        Measure every replicate even when a matching stored series exists.
    output_dir : str or Path, optional
        Folder that holds ``polyzymd_results``. Defaults to the current directory.
    **kwargs
        Keyword arguments of ``function``, recorded like ``args``.

    Returns
    -------
    Timeseries
        Every replicate's series, by condition.

    Raises
    ------
    ProtocolError
        If ``function`` returns something other than one number per frame, or
        a selection matches no atoms.
    """
    import numpy as np
    from MDAnalysis.analysis.base import AnalysisFromFunction

    name = name or getattr(function, "__name__", "timeseries")
    root = Path(output_dir or Path.cwd()).expanduser().resolve() / RESULTS_DIR / _safe(name)
    base = {
        "name": name,
        "function": _function_record(function),
        "arguments": {
            "args": [_argument_record(value) for value in args],
            "kwargs": {key: _argument_record(value) for key, value in kwargs.items()},
        },
        "unit": unit,
    }
    series: dict[str, list[ReplicateSeries]] = {}
    for condition in study:
        series[condition.label] = []
        for replicate in condition.replicates:
            record = json.loads(json.dumps(_replicate_record(base, replicate)))
            folder = root / _safe(condition.label) / f"replicate_{replicate.index}"
            values = None if recompute else _stored_values(folder, record)
            if values is None:
                u = replicate.universe()
                built, chosen = _build_arguments({**dict(enumerate(args)), **kwargs}, replicate)
                analysis = AnalysisFromFunction(
                    function,
                    u.trajectory,
                    *(built[index] for index in range(len(args))),
                    **{key: built[key] for key in kwargs},
                ).run(frames=replicate.frames)
                values = np.asarray(analysis.results.timeseries, dtype=np.float64)
                if values.shape != (len(record["frames"]),):
                    raise ProtocolError(
                        f"{name}: {function!r} returned values of shape {values.shape[1:]} per "
                        "frame; this version takes one number per frame.",
                        hint="Return a single float from the function.",
                    )
                folder.mkdir(parents=True, exist_ok=True)
                np.savez(
                    folder / "series.npz",
                    values=values,
                    frames=replicate.frames,
                    times=replicate.times,
                )
                (folder / "record.json").write_text(
                    json.dumps({**record, "chosen": chosen, "versions": _versions()}, indent=1)
                )
            series[condition.label].append(
                ReplicateSeries(
                    condition.label,
                    replicate.index,
                    values,
                    replicate.frames,
                    replicate.times,
                    folder,
                )
            )
    return Timeseries(name, unit, study, series, root)


def _replicate_record(base: dict[str, Any], replicate: Replicate) -> dict[str, Any]:
    """Complete the shared record with one replicate's inputs, frames and times."""
    return {
        **base,
        "condition": replicate.condition.label,
        "replicate": replicate.index,
        **replicate.identity,
        "frames": replicate.frames.tolist(),
        "times_ns": replicate.times.tolist(),
    }


def _stored_values(folder: Path, record: dict[str, Any]) -> np.ndarray | None:
    """Return the stored values when the stored record matches, apart from versions and chosen."""
    import numpy as np

    try:
        stored = json.loads((folder / "record.json").read_text())
        stored.pop("versions", None)
        stored.pop("chosen", None)
        if stored != record:
            return None
        with np.load(folder / "series.npz") as data:
            return np.asarray(data["values"], dtype=np.float64)
    except (OSError, ValueError, KeyError):
        return None


class Timeseries:
    """Per-frame values of every replicate of a study, by condition.

    Attributes
    ----------
    name : str
        Result name.
    unit : str or None
        Unit of the values.
    series : dict of str to list of ReplicateSeries
        Each condition's replicate series, in replicate order.
    path : Path
        Folder the series are stored under.
    """

    def __init__(
        self,
        name: str,
        unit: str | None,
        study: Study,
        series: dict[str, list[ReplicateSeries]],
        path: Path,
    ) -> None:
        self.name, self.unit, self.study, self.series, self.path = name, unit, study, series, path

    def reduce(self, how: str | Callable = "mean", *, unit: Any = ...) -> ReplicateValues:
        """Turn each replicate's series into one value.

        Parameters
        ----------
        how : {"mean", "fraction", "std"} or callable
            ``"mean"`` averages the frames, ``"fraction"`` averages a series of
            0 and 1 and rejects any other value, ``"std"`` takes the sample
            standard deviation over frames (``ddof=1``). A callable receives
            the values and times in ns of one replicate and returns a number.
        unit : str or None, optional
            Unit of the reduced value. Defaults to the series unit, or
            ``None`` for ``"fraction"``.

        Returns
        -------
        ReplicateValues
            One value per replicate, with the statistical inefficiency and
            effective sample size of each series.
        """
        import numpy as np

        from polyzymd.analyses.shared.autocorrelation import (
            detect_equilibration,
            n_effective,
            statistical_inefficiency,
        )

        def fraction(values: np.ndarray, times: np.ndarray) -> float:
            if not np.isin(values, (0.0, 1.0)).all():
                raise ProtocolError(
                    f"{self.name}: 'fraction' needs a series of 0 and 1.",
                    hint="Use reduce('mean') for any other series.",
                )
            return float(np.mean(values))

        named = {
            "mean": lambda values, times: float(np.mean(values)),
            "std": lambda values, times: float(np.std(values, ddof=1)),
            "fraction": fraction,
        }
        if not callable(how) and how not in named:
            raise ProtocolError(
                f"Unknown reduction {how!r}.", hint="Use 'mean', 'fraction', 'std' or a function."
            )
        reducer = how if callable(how) else named[how]
        label = getattr(how, "__name__", "reduced") if callable(how) else how
        rows: dict[str, list[tuple]] = {}
        for condition, items in self.series.items():
            rows[condition] = []
            for item in items:
                g = statistical_inefficiency(item.values)
                start = detect_equilibration(item.values)
                rows[condition].append(
                    (
                        item.replicate,
                        float(reducer(item.values, item.times)),
                        g,
                        n_effective(len(item.values), g),
                        len(item.values),
                        start,
                        float(item.times[start]),
                    )
                )
        if unit is ...:
            unit = None if how == "fraction" else self.unit
        return ReplicateValues(self, f"{label}_{self.name}", unit, label == "fraction", rows)


class ReplicateValues:
    """One value per replicate of every condition, ready to summarise and compare."""

    def __init__(
        self,
        source: Timeseries,
        metric: str,
        unit: str | None,
        is_fraction: bool,
        rows: dict[str, list[tuple]],
    ) -> None:
        self.source, self.metric, self.unit = source, metric, unit
        self.is_fraction, self.rows = is_fraction, rows

    @property
    def values(self) -> dict[str, list[float]]:
        """Each condition's replicate values, in replicate order."""
        return {label: [row[1] for row in rows] for label, rows in self.rows.items()}

    def summary(self, conditions: Sequence[str] | None = None) -> ProtocolReport:
        """Give each condition's n, mean, standard error and 95 percent interval.

        The interval is the mean plus or minus the Student t coverage factor
        for n replicates times the standard error, from
        :func:`~polyzymd.analyses.shared.statistics.mean_sem_ci`. Every
        replicate value is listed with the statistical inefficiency and
        effective sample size of its series.

        Parameters
        ----------
        conditions : sequence of str, optional
            Conditions to include. Defaults to every condition of the study.

        Returns
        -------
        ProtocolReport
            One condition row each and no comparisons.
        """
        return self._report(self._chosen(conditions), [])

    def compare(
        self,
        control: str | None = None,
        conditions: Sequence[str] | None = None,
        test: str = "welch",
    ) -> ProtocolReport:
        """Compare every condition against the control.

        Each row gives ``mean(b) - mean(a)`` with its 95 percent interval from
        the same t test, the p value, the Benjamini-Hochberg adjusted p value
        across the rows of this call, and Cohen's d and Hedges' g oriented
        like the difference. A row is not testable when a condition has fewer
        than two replicates, or when both conditions have the same value in
        every replicate; it then takes no part in the correction.

        Parameters
        ----------
        control : str, optional
            Control condition. Defaults to the first condition of the study.
        conditions : sequence of str, optional
            Conditions to include, with the control. Defaults to all.
        test : {"welch", "student"}, optional
            Welch's t test (default) or Student's t test.

        Returns
        -------
        ProtocolReport
            The condition rows and one comparison row per non-control condition.
        """
        from polyzymd.analyses.protocols import PairwiseReport, _difference_ci
        from polyzymd.analyses.shared.inferential_statistics import (
            NO_SIGNIFICANT_CHANGE,
            benjamini_hochberg,
            cohens_d,
            independent_ttest,
        )

        control = control or self.source.study.control
        chosen = self._chosen(conditions)
        if control not in chosen:
            chosen = [control, *chosen]
        if test not in ("welch", "student"):
            raise ProtocolError(f"Unknown test {test!r}.", hint="Use 'welch' or 'student'.")
        values = self.values
        rows = []
        for label in chosen:
            if label == control:
                continue
            a, b = values[control], values[label]
            testable = min(len(a), len(b)) >= 2 and len(set(a)) + len(set(b)) > 2
            effect = cohens_d(b, a)
            rows.append(
                PairwiseReport(
                    a=control,
                    b=label,
                    delta=sum(b) / len(b) - sum(a) / len(a),
                    delta_ci95=_difference_ci(a, b, f"{test}_t"),
                    p=independent_ttest(b, a, method=test).p_value if testable else None,
                    test=f"{test}_t",
                    cohens_d=_finite(effect.cohens_d),
                    hedges_g=_finite(effect.hedges_g),
                    testable=testable,
                )
            )
        for row, corrected in zip(rows, benjamini_hochberg([row.p for row in rows]), strict=True):
            row.p_adjusted = corrected.adjusted_p_value
            row.significant = corrected.significant
            row.direction = (
                ("increased" if row.delta > 0 else "decreased")
                if row.significant
                else NO_SIGNIFICANT_CHANGE
            )
        return self._report(chosen, rows)

    def _chosen(self, conditions: Sequence[str] | None) -> list[str]:
        """Return the requested conditions in study order, rejecting unknown ones."""
        if conditions is None:
            return list(self.rows)
        unknown = [label for label in conditions if label not in self.rows]
        if unknown:
            raise ProtocolError(
                f"Unknown condition(s) {unknown}.", hint=f"Use labels from {list(self.rows)}."
            )
        return [label for label in self.rows if label in conditions]

    def _report(self, chosen: list[str], pairwise: list[Any]) -> ProtocolReport:
        """Build the report for the chosen conditions and comparison rows."""
        from polyzymd.analyses.protocols import (
            ProtocolProvenance,
            ProtocolReport,
            _condition,
            _verdict,
        )

        conditions, notes = [], []
        for label in chosen:
            rows = self.rows[label]
            item = _condition(label, [row[1] for row in rows], {}, self.metric)
            item.replicates = [row[0] for row in rows]
            item.statistical_inefficiency = [row[2] for row in rows]
            item.n_effective = [row[3] for row in rows]
            item.eq_detected_frame = [row[5] for row in rows]
            item.eq_detected_ns = [row[6] for row in rows]
            for row in rows:
                if row[5] > EQUILIBRATION_WARNING_FRACTION * row[4]:
                    notes.append(
                        f"condition {label} replicate {row[0]}: pymbar detect_equilibration puts "
                        f"the start of the equilibrated region at {row[6]:.4g} ns, production "
                        f"frame {row[5] + 1} of {row[4]}, after the equilibration window; the "
                        "window may be too short for it"
                    )
            if len(rows) > 1 and len({row[1] for row in rows}) == 1:
                item.ci95, item.ci_method = None, "not_estimable"
                notes.append(
                    f"condition {label} has the same {self.metric} in every replicate, so its "
                    "interval is not estimable"
                )
            elif len(rows) < 2:
                notes.append(f"condition {label} has one replicate, so it has no interval")
            if self.is_fraction and item.ci95 and (item.ci95[0] < 0 or item.ci95[1] > 1):
                notes.append(
                    f"the 95 percent interval of condition {label} extends past the fraction "
                    "bounds 0 and 1, where a t interval is not reliable"
                )
            conditions.append(item)
        for row in pairwise:
            if not row.testable:
                notes.append(
                    f"{row.a} vs {row.b} is not testable: a condition has fewer than two "
                    "replicates, or both have one value in every replicate"
                )
        study, versions = self.source.study, _versions()
        return ProtocolReport(
            analysis=self.source.name,
            protocol_version="2",
            metric=self.metric,
            unit=self.unit,
            equilibration=study[chosen[0]].equilibration,
            frames_per_replicate={label: [row[4] for row in self.rows[label]] for label in chosen},
            conditions=conditions,
            pairwise=pairwise,
            warnings=notes,
            provenance=ProtocolProvenance(
                polyzymd_version=versions["polyzymd"],
                mdanalysis_version=versions["MDAnalysis"],
                config_hashes={label: study[label].config_hash for label in chosen},
                output_paths={"results": str(self.source.path)},
            ),
            verdict=_verdict(self.metric, self.unit, conditions, pairwise),
        )


def _finite(value: float) -> float | None:
    """Return ``value``, or ``None`` when it is NaN or infinite."""
    return value if math.isfinite(value) else None
