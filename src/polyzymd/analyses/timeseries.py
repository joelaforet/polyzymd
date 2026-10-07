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
Benjamini, Y. (2010). Discovering the false discovery rate. Journal of the
    Royal Statistical Society Series B 72:405-416.
    doi:10.1111/j.1467-9868.2010.00746.x
"""

from __future__ import annotations

import functools
import hashlib
import inspect
import json
import logging
import math
import re
import sys
import warnings
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Mapping, Sequence

from polyzymd.analyses.exceptions import ProtocolError

LOGGER = logging.getLogger(__name__)

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

#: The detected start is judged only for a replicate with at least this many
#: effective samples. On the stationary AR(1) series above, the start passed
#: 10 percent in 1 to 7 percent of series with 100 or more effective samples
#: and in 45 percent of series with about 10. On the LipA 363 K rmsd
#: replicates, with 3 to 18 effective samples, pymbar put the start within
#: the last few percent of 28 of 30 runs, on a tail short enough that its
#: statistical inefficiency is near 1. Below 20 the start says more about the
#: length of the run than about the window, so it is reported without a warning.
EQUILIBRATION_MIN_N_EFFECTIVE = 20


@dataclass(frozen=True)
class Select:
    """An MDAnalysis selection string, built into an ``AtomGroup`` per replicate."""

    selection: str
    allow_empty: bool = False


@dataclass(frozen=True)
class UniverseArgument:
    """Stands for the replicate's ``Universe`` in the arguments of a function."""


def select(selection: str, *, allow_empty: bool = False) -> Select:
    """Stand for ``universe.select_atoms(selection)`` of each replicate.

    Parameters
    ----------
    selection : str
        MDAnalysis selection string. It is recorded with the result.
    allow_empty : bool, optional
        Give an empty ``AtomGroup`` for a replicate where ``selection``
        matches no atoms, instead of refusing it; for analyses that leave such
        replicates out themselves.

    Returns
    -------
    Select
        Placeholder that :func:`run_timeseries` replaces with an ``AtomGroup``.
    """
    return Select(str(selection), allow_empty)


def universe() -> UniverseArgument:
    """Stand for each replicate's ``Universe`` in the arguments of a function.

    Returns
    -------
    UniverseArgument
        Placeholder that :func:`run_timeseries` replaces with the ``Universe``.
    """
    return UniverseArgument()


def _function_record(function: Callable) -> dict[str, Any]:
    """Name a function and hash its source, or its bytecode when no source exists.

    A function loaded from a study's or project's own file
    (:func:`~polyzymd.analyses.user_functions.load_function`) is hashed with
    the Python files under that file's folder and the files of its ``data/``
    folder (:func:`code_files`), by path and content: the folder is on
    ``sys.path`` while the function runs, so a helper module or package it
    imports from there, or a data file it reads from ``data/``, is part of
    what produced the result. Editing any of those files changes the record.

    A function shipped with PolyzyMD is hashed with every Python file of
    ``polyzymd/analyses`` (:func:`_shipped_code_hash`), so a fix to any helper
    it calls, such as ``shared/centroid.py``, recomputes the results it made.
    """
    module_file = getattr(function, "__polyzymd_module_file__", None)
    module = getattr(function, "__module__", None) or ""
    try:
        if module_file:
            digest, basis = folder_hash(Path(module_file)), "module_folder"
        elif module.startswith("polyzymd."):
            digest, basis = _shipped_code_hash(), "polyzymd_modules"
        else:
            source = inspect.getsource(function).encode()
            digest, basis = hashlib.sha256(source).hexdigest(), "source"
    except (OSError, TypeError):
        warnings.warn(
            f"No source file for {function!r}, so its bytecode is hashed instead; a changed "
            "function can keep the same hash, and the stored record identifies it less reliably.",
            stacklevel=4,
        )
        body = getattr(function, "__code__", None)
        code = repr((body.co_code, body.co_consts) if body else function).encode()
        digest, basis = hashlib.sha256(code).hexdigest(), "bytecode"
    return {
        "qualname": getattr(function, "__qualname__", repr(function)),
        "module": getattr(function, "__module__", None),
        "hash": digest,
        "hash_of": basis,
    }


#: PolyzyMD modules whose text never changes a result: error messages and hints.
_NOT_RESULT_MODULES = {"polyzymd.analyses.exceptions"}


def _json_number(value: float) -> float | str:
    """Return ``value``, or its name when it is not finite, so a record compares equal to itself."""
    return value if math.isfinite(value) else str(value)


@functools.lru_cache(maxsize=None)
def _shipped_code_hash(package: Path = Path(__file__).parent) -> str:
    """Return the SHA-256 of every Python file of ``polyzymd/analyses``.

    Each file is hashed by its path inside the package and its content, in
    path order. :data:`_NOT_RESULT_MODULES` are left out. A change to any
    other analysis file changes the hash of every shipped analysis.
    """
    left_out = {name.rsplit(".", 1)[1] + ".py" for name in _NOT_RESULT_MODULES}
    digest = hashlib.sha256()
    for path in sorted(package.rglob("*.py")):
        name = path.relative_to(package).as_posix()
        if name not in left_out:
            digest.update(name.encode() + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()


#: Folders whose files are outputs or inputs of their own, never code a function uses.
_NOT_CODE = {
    "__pycache__",
    "results",
    "logs",
    "deposit",
    "conditions",
    "figures",
    "slurm",
    "slurm_logs",
    "runs",
    "site-packages",
}
#: Files beside code that say where data is on one machine, or that freeze writes.
_NOT_CODE_FILES = {
    "data.local.yaml",
    "manifest.json",
    "md_checklist.yaml",
    "system_summary.csv",
    "CITATION.cff",
    ".zenodo.json",
}


#: The folder, beside a function's file, whose files of any type count as its data.
DATA_FOLDER = "data"


def _candidate_files(folder: Path) -> list[tuple[Path, Path]]:
    """Return ``(path, relative path)`` of each file under ``folder`` outside skipped folders.

    Skipped: hidden files and folders (``.git``, ``.pixi``, ``.venv``),
    output and job folders (``results/``, ``logs/``, ``deposit/``,
    ``conditions/``, ``figures/``, ``slurm/``, ``runs/``, ``site-packages/``,
    ``__pycache__``), the folders of other studies and the files that
    ``freeze`` writes.
    """
    found = []
    for path in sorted(folder.rglob("*")):
        relative = path.relative_to(folder)
        if (
            not path.is_file()
            or path.name in _NOT_CODE_FILES
            or any(part.startswith(".") for part in relative.parts)
            or any(part in _NOT_CODE for part in relative.parts[:-1])
        ):
            continue
        if any(
            (folder / Path(*relative.parts[:i])).joinpath("study.yaml").is_file()
            for i in range(1, len(relative.parts))
        ):
            continue
        found.append((path, relative))
    return found


def code_files(folder: Path) -> list[Path]:
    """Return the files a function or plan in ``folder`` may use, in path order.

    Two kinds of files count, in subfolders too:

    - every Python file (``*.py``): the function's own file, and the helper
      modules and packages it imports from ``folder``;
    - every file under ``folder/data/``: the data files the function reads.

    When ``folder`` holds ``study.yaml`` or ``project.yaml`` (a function file
    beside it), only the Python files count: a ``data/`` folder there may
    hold trajectories. Keep data files in ``analyses/data/`` or
    ``stats/data/``. Other files, such as notes, figures or a copied
    trajectory, never count, so adding them recomputes nothing.
    :func:`ignored_files` names them.
    Skipped folders are listed in :func:`_candidate_files`.
    """
    root_of_study = (folder / "study.yaml").is_file() or (folder / "project.yaml").is_file()
    return [
        path
        for path, relative in _candidate_files(folder)
        if path.suffix == ".py" or (relative.parts[0] == DATA_FOLDER and not root_of_study)
    ]


def ignored_files(folder: Path) -> list[Path]:
    """Return the files under ``folder`` that :func:`code_files` leaves out of the hash.

    Compiled bytecode (``*.pyc``) is not listed. For a study's or project's
    own folder, where only Python files count by design, the list is empty.
    """
    if (folder / "study.yaml").is_file() or (folder / "project.yaml").is_file():
        return []
    counted = set(code_files(folder))
    return [
        path
        for path, _ in _candidate_files(folder)
        if path not in counted and path.suffix != ".pyc"
    ]


def folder_hash(module_file: Path) -> str:
    """Return the SHA-256 of the path and content of every file :func:`code_files` gives.

    The files are read in blocks, so a large data file beside the code is
    never held in memory whole.
    """
    folder = module_file.parent
    ignored = ignored_files(folder)
    if ignored:
        LOGGER.info(
            "The stored results of %s do not depend on %s: only Python files and the "
            "files in %s/ count. Move a file the function reads into %s/.",
            module_file.name,
            ", ".join(path.relative_to(folder).as_posix() for path in ignored),
            DATA_FOLDER,
            DATA_FOLDER,
        )
    digest = hashlib.sha256()
    for path in code_files(folder):
        digest.update(path.relative_to(folder).as_posix().encode() + b"\0")
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 22), b""):
                digest.update(block)
        digest.update(b"\0")
    return digest.hexdigest()


def _file_record(path: str | Path) -> dict[str, str]:
    """Record a file by its name and the SHA-256 hash of its content.

    Not by its location, so a study folder moved, copied or migrated to
    another machine reuses its stored results, while a changed file changes
    the record.
    """
    path = Path(path).expanduser().resolve()
    return {"name": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


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
        if len(atoms) == 0 and not value.allow_empty:
            raise ProtocolError(
                f"Selection {value.selection!r} matched no atoms.",
                hint="Check the selection string against the topology. If some conditions have "
                "no such atoms (a polymer selection in a no-polymer control), write "
                "'allow_empty: true' in the study.yaml entry to pass the function an empty "
                "AtomGroup there.",
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
    bounds: tuple[float | None, float | None] = (None, None),
    parts: Sequence[str] | None = None,
    **kwargs: Any,
) -> Timeseries | dict[str, Timeseries]:
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
    versions and the bounds equals the new one. The bounds change only how
    a distribution is drawn, so a reused series gets the bounds of this call
    written into its record. A
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
    bounds : tuple of (float or None, float or None), optional
        Lowest and highest value the quantity can take, ``None`` for no
        limit, such as ``(0.0, None)`` for a distance. Recorded with the
        result and used to correct distribution figures at the limits.
    parts : sequence of str, optional
        Names of several quantities ``function`` measures per frame, for
        example ``["area", "contacts"]``. ``function`` then returns, per
        frame, a dict with exactly these keys or a sequence in this order.
        They are stored as one ``(frames, parts)`` array; each part becomes
        its own series, named by its part.
    **kwargs
        Keyword arguments of ``function``, recorded like ``args``.

    Returns
    -------
    Timeseries or dict of str to Timeseries
        Every replicate's series, by condition, or with ``parts`` one such
        result per part.

    Raises
    ------
    ProtocolError
        If ``function`` returns something other than one number per frame, or
        a selection matches no atoms.
    """
    import numpy as np
    from MDAnalysis.analysis.base import AnalysisFromFunction

    from polyzymd.analyses.shared.loader import underlying_reader

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
        # The order of the parts names the stored columns, so it keys reuse.
        **({"parts": list(parts)} if parts is not None else {}),
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
                    function if parts is None else _per_part(function, list(parts), name),
                    # MDAnalysis passes anything that is not its own reader to
                    # the function as an argument, so hand it the unwrapped one.
                    underlying_reader(u.trajectory),
                    *(built[index] for index in range(len(args))),
                    **{key: built[key] for key in kwargs},
                ).run(frames=replicate.frames)
                n_frames = len(record["frames"])
                expected = (n_frames,) if parts is None else (n_frames, len(parts))
                several = (
                    "Return one float per frame, or name the quantities with parts: [...] in "
                    "the study.yaml entry (parts= in Python) and return a dict of them."
                )
                try:
                    values = np.asarray(analysis.results.timeseries, dtype=np.float64)
                except (TypeError, ValueError) as exc:
                    raise ProtocolError(
                        f"{name}: {function!r} returned something other than numbers per "
                        f"frame: {exc}",
                        hint=several,
                    ) from exc
                if values.shape != expected:
                    raise ProtocolError(
                        f"{name}: {function!r} returned values of shape {values.shape[1:]} per "
                        f"frame, not {expected[1:] or 'one number'}.",
                        hint=several,
                    )
                folder.mkdir(parents=True, exist_ok=True)
                np.savez(
                    folder / "series.npz",
                    values=values,
                    frames=replicate.frames,
                    times=replicate.times,
                )
                stored = {**record, "bounds": list(bounds), "chosen": chosen}
                (folder / "record.json").write_text(
                    json.dumps({**stored, "versions": _versions()}, indent=1)
                )
                if parts is not None:
                    # The part names, for polyzymd.analyses.results.
                    (folder / "parts.json").write_text(json.dumps(list(parts)))
            else:
                stored = json.loads((folder / "record.json").read_text())
                if stored.get("bounds") != list(bounds):
                    (folder / "record.json").write_text(
                        json.dumps({**stored, "bounds": list(bounds)}, indent=1)
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
    if parts is None:
        return Timeseries(name, unit, study, series, root, tuple(bounds))
    return {
        part: Timeseries(
            part,
            unit,
            study,
            {c: [replace(s, values=s.values[:, i]) for s in items] for c, items in series.items()},
            root,
            tuple(bounds),
        )
        for i, part in enumerate(parts)
    }


def _per_part(function: Callable, parts: list[str], name: str) -> Callable:
    """Wrap ``function`` so a dict it returns per frame becomes its values in ``parts`` order."""
    from collections.abc import Mapping

    def per_frame(*args: Any, **kwargs: Any) -> Any:
        out = function(*args, **kwargs)
        if not isinstance(out, Mapping):
            return out
        if set(out) != set(parts):
            raise ProtocolError(
                f"{name}: {function!r} returned the parts {sorted(out)}, not {parts}.",
                hint="Return a dict with exactly the names listed under parts.",
            )
        return [out[part] for part in parts]

    return per_frame


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


def _stored_values(
    folder: Path, record: dict[str, Any], file: str = "series.npz"
) -> np.ndarray | None:
    """Return the stored values when the stored record matches, apart from versions, chosen and bounds."""
    import numpy as np

    try:
        stored = json.loads((folder / "record.json").read_text())
        stored.pop("versions", None)
        stored.pop("chosen", None)
        stored.pop("bounds", None)
        if stored != record:
            return None
        with np.load(folder / file) as data:
            return np.asarray(data["values"], dtype=np.float64)
    except (OSError, ValueError, KeyError):
        return None


@dataclass
class Source:
    """The study, name, unit and folder of per-replicate values computed by a function."""

    name: str
    unit: str | None
    study: Study
    series: dict[str, list[Path]]
    path: Path


def run_per_replicate(
    study: Study,
    function: Callable,
    *args: Any,
    unit: str | None,
    labels: Sequence | Callable | None = None,
    missing: float | None = None,
    note_filled: bool = False,
    name: str | None = None,
    recompute: bool = False,
    output_dir: str | Path | None = None,
    bounds: tuple[float | None, float | None] = (None, None),
    parts: Sequence[str] | None = None,
    **kwargs: Any,
) -> ReplicateValues | dict[str, ReplicateValues]:
    """Compute one value, or one labelled array, per replicate with ``function``.

    ``function`` is called once per replicate with the built ``args``, the
    keyword ``frames`` holding the replicate's production frame indices and
    the built ``kwargs``. Arguments are built and recorded as in
    :func:`run_timeseries`. Each replicate's values go to ``values.npz`` and
    its record, which also holds the labels, to ``record.json`` in
    ``<output_dir>/polyzymd_results/<name>/<condition>/replicate_<n>/``, and
    a stored result is reused on the same terms.

    Parameters
    ----------
    study : Study
        The study to measure.
    function : callable
        Returns a number, or a one-dimensional array with one entry per label.
    *args
        Arguments of ``function``, as in :func:`run_timeseries`.
    unit : str or None
        Unit of the values.
    labels : sequence, callable or "returned", optional
        Name of each entry of the returned array, or a function of the
        replicate's ``Universe`` that returns them, such as residue IDs.
        With ``"returned"``, ``function`` returns ``(labels, values)``, for
        results whose entries are known only once measured, such as the
        residue pairs that formed a bond; pass ``missing`` for a label that
        only some replicates have. Replicates are lined up by label, never
        by position.
    missing : float, optional
        Value given to a label that one replicate lacks and another has.
        By default a missing label is an error.
    note_filled : bool, optional
        Name, in the report's warnings, every replicate given ``missing`` and
        its labels. For a value chosen to stand in for one that was not
        measured, as a study's ``missing:``; not for a shipped analysis whose
        ``missing`` is a definition (an unformed hydrogen bond is zero).
    name : str, optional
        Result name. Defaults to the function's ``__name__``.
    recompute : bool, optional
        Compute every replicate even when a matching stored result exists.
    output_dir : str or Path, optional
        Folder that holds ``polyzymd_results``. Defaults to the current directory.
    bounds : tuple of (float or None, float or None), optional
        Lowest and highest value the quantity can take, ``None`` for no
        limit, used to warn when a 95 percent interval extends past them.
        They change no value and are not part of the record.
    parts : sequence of str, optional
        Names of the rows of a two-dimensional array that ``function``
        returns with one column per label, for several quantities measured
        in one pass. Each row becomes its own result, named by its part.
    **kwargs
        Keyword arguments of ``function``, recorded like ``args``.

    Returns
    -------
    ReplicateValues or dict of str to ReplicateValues
        One value or one labelled array per replicate, or with ``parts`` one
        such result per part, sharing one record per replicate.

    Raises
    ------
    ProtocolError
        If the returned shape does not match the labels, labels repeat, or a
        label is missing from a replicate and ``missing`` is not given.
    """
    import numpy as np

    name = name or getattr(function, "__name__", "per_replicate")
    root = Path(output_dir or Path.cwd()).expanduser().resolve() / RESULTS_DIR / _safe(name)
    base = {
        "name": name,
        "function": _function_record(function),
        "arguments": {
            "args": [_argument_record(value) for value in args],
            "kwargs": {key: _argument_record(value) for key, value in kwargs.items()},
        },
        "unit": unit,
        # The order of the parts names the stored columns, so it keys reuse.
        **({"parts": list(parts)} if parts is not None else {}),
        # Readers fill a label a replicate lacks with it, as the report does.
        **({"missing": _json_number(missing)} if missing is not None else {}),
    }
    rows: dict[str, list[tuple]] = {}
    found: dict[str, list[Path]] = {}
    for condition in study:
        rows[condition.label], found[condition.label] = [], []
        for replicate in condition.replicates:
            returned = isinstance(labels, str) and labels == "returned"
            given = (
                None if returned else labels(replicate.universe()) if callable(labels) else labels
            )
            given = (
                None if given is None else [x.item() if hasattr(x, "item") else x for x in given]
            )
            record = json.loads(
                json.dumps(
                    {
                        **_replicate_record(base, replicate),
                        "labels": "returned" if returned else given,
                    }
                )
            )
            folder = root / _safe(condition.label) / f"replicate_{replicate.index}"
            values = None if recompute else _stored_values(folder, record, "values.npz")
            if values is not None and returned:
                try:
                    given = json.loads((folder / "labels.json").read_text())
                except (OSError, ValueError):
                    values = None
            if values is None:
                built, chosen = _build_arguments({**dict(enumerate(args)), **kwargs}, replicate)
                extra = {"times": replicate.times} if _takes(function, "times") else {}
                output = function(
                    *(built[index] for index in range(len(args))),
                    frames=replicate.frames,
                    **extra,
                    **{key: built[key] for key in kwargs},
                )
                if returned:
                    given, output = output
                    given = json.loads(
                        json.dumps([x.item() if hasattr(x, "item") else x for x in given])
                    )
                values = np.asarray(output, dtype=np.float64)
                expected = () if given is None else (len(given),)
                expected = expected if parts is None else (len(parts), *expected)
                if values.shape != expected:
                    raise ProtocolError(
                        f"{name}: {function!r} returned shape {values.shape} for "
                        f"{'no labels' if given is None else f'{len(given)} labels'}.",
                        hint="Return one number, or one value per label.",
                    )
                folder.mkdir(parents=True, exist_ok=True)
                np.savez(folder / "values.npz", values=values)
                if returned:
                    (folder / "labels.json").write_text(json.dumps(given))
                (folder / "record.json").write_text(
                    json.dumps({**record, "chosen": chosen, "versions": _versions()}, indent=1)
                )
                if parts is not None:
                    # The part names, for polyzymd.analyses.results; the record keys them.
                    (folder / "parts.json").write_text(json.dumps(list(parts)))
            if given is None:
                value = float(values) if parts is None else values.tolist()
            else:
                pairs = [dict(zip(given, row)) for row in np.atleast_2d(values).tolist()]
                value = pairs[0] if parts is None else pairs
            if given is not None and len(pairs[0]) != len(given):
                raise ProtocolError(
                    f"{name}: condition {condition.label} replicate {replicate.index} repeats "
                    "a label.",
                    hint="Give every entry its own label, for example segment and residue ID.",
                )
            rows[condition.label].append(
                (replicate.index, value, None, None, len(replicate.frames), None, None)
            )
            found[condition.label].append(folder)
    source = Source(name, unit, study, found, root)

    def result(table: dict[str, list[tuple]], metric: str) -> ReplicateValues:
        filled: list[str] = []
        order = None if labels is None else _label_order(table, missing, metric, filled)
        values = ReplicateValues(source, metric, unit, False, table, order)
        values.bounds = tuple(bounds)
        values.filled = filled if note_filled else []
        return values

    if parts is None:
        return result(rows, name)
    return {
        part: result(
            {c: [(r[0], r[1][i], *r[2:]) for r in items] for c, items in rows.items()}, part
        )
        for i, part in enumerate(parts)
    }


def _takes(function: Callable, name: str) -> bool:
    """Return whether ``function`` has a parameter called ``name``."""
    try:
        return name in inspect.signature(function).parameters
    except (TypeError, ValueError):
        return False


def _label_order(
    rows: dict[str, list[tuple]],
    missing: float | None,
    name: str,
    filled: list[str] | None = None,
) -> list:
    """Replace each replicate's label mapping with an array in the order labels first appear.

    A label a replicate lacks gets ``missing``; each replicate that got one is
    named in ``filled``, with its labels, for the report's warnings.
    """
    import numpy as np

    order = list(dict.fromkeys(key for items in rows.values() for row in items for key in row[1]))
    for condition, items in rows.items():
        for position, row in enumerate(items):
            absent = [key for key in order if key not in row[1]]
            if absent and missing is None:
                raise ProtocolError(
                    f"{name}: condition {condition} replicate {row[0]} has no value for "
                    f"labels {absent}, which other replicates have.",
                    hint="Measure the same labels in every replicate, or give the value a "
                    "replicate gets for a label it lacks: missing: .nan in the study.yaml entry "
                    "(missing=<value> in Python).",
                )
            if absent and filled is not None:
                filled.append(
                    f"condition {condition} replicate {row[0]}: {len(absent)} of {len(order)} "
                    f"labels had no value and were given missing={missing}: "
                    + ", ".join(str(key) for key in absent)
                )
            values = np.array([row[1].get(key, missing) for key in order], dtype=np.float64)
            items[position] = (row[0], values, *row[2:])
    return order


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
    bounds : tuple of (float or None, float or None)
        Lowest and highest value the quantity can take, ``None`` for no limit.
    """

    def __init__(
        self,
        name: str,
        unit: str | None,
        study: Study,
        series: dict[str, list[ReplicateSeries]],
        path: Path,
        bounds: tuple[float | None, float | None] = (None, None),
    ) -> None:
        self.name, self.unit, self.study, self.series, self.path = name, unit, study, series, path
        self.bounds = bounds

    def transform(
        self,
        function: Callable,
        *others: Timeseries,
        unit: Any = ...,
        name: str | None = None,
        bounds: Any = ...,
        **kwargs: Any,
    ) -> Timeseries:
        """Compute a new series from the stored values, without reading a trajectory.

        For every replicate, ``function`` receives the values of this series
        and of each series in ``others`` as NumPy arrays with one entry per
        frame, followed by ``kwargs``, and returns one value per frame, for
        example ``lambda d: d < 3.5``. Each replicate's result goes to
        ``series.npz`` and ``record.json`` in
        ``polyzymd_results/<name>/<condition>/replicate_<n>/``, next to this
        series. The record holds the function's name, module and source hash,
        ``kwargs``, the unit, and the path and SHA-256 hash of the record of
        every input series. It is written again on every call. A value that
        ``function`` reads from an enclosing scope is not recorded, so pass
        such values in ``kwargs``.

        Parameters
        ----------
        function : callable
            Called per replicate as ``function(values, *other_values, **kwargs)``.
        *others : Timeseries
            Further series of the same study with the same frames per replicate.
        unit : str or None, optional
            Unit of the new values. Defaults to the unit of this series.
        name : str, optional
            Result name. Defaults to the function's ``__name__`` and this name.
        bounds : tuple of (float or None, float or None), optional
            Lowest and highest value of the new quantity. Defaults to the
            bounds of this series.
        **kwargs
            Keyword arguments of ``function``, recorded with the result.

        Returns
        -------
        Timeseries
            The new per-frame values of every replicate.

        Raises
        ------
        ProtocolError
            If a series in ``others`` lacks a replicate or has other frames,
            or ``function`` returns other than one value per frame.
        """
        import numpy as np

        name = name or f"{getattr(function, '__name__', 'transform').strip('<>')}_{self.name}"
        unit = self.unit if unit is ... else unit
        bounds = self.bounds if bounds is ... else tuple(bounds)
        root = self.path.parent / _safe(name)
        base = {
            "name": name,
            "transform": _function_record(function),
            "kwargs": {key: _argument_record(value) for key, value in kwargs.items()},
            "unit": unit,
            "bounds": list(bounds),
        }
        series: dict[str, list[ReplicateSeries]] = {}
        for condition, items in self.series.items():
            series[condition] = []
            for index, item in enumerate(items):
                inputs = [item]
                for other in others:
                    match = other.series.get(condition, [])[index : index + 1]
                    if not match or not np.array_equal(match[0].frames, item.frames):
                        raise ProtocolError(
                            f"{name}: {other.name} has no series with the frames of {self.name} "
                            f"for condition {condition} replicate {item.replicate}.",
                            hint="Transform series measured on the same study and window.",
                        )
                    inputs.append(match[0])
                values = function(*(entry.values for entry in inputs), **kwargs)
                values = np.asarray(values, dtype=np.float64)
                if values.shape != item.values.shape:
                    raise ProtocolError(
                        f"{name}: {function!r} returned shape {values.shape} for "
                        f"{len(item.values)} frames.",
                        hint="Return one value per frame, for example lambda d: d < 3.5.",
                    )
                folder = root / _safe(condition) / f"replicate_{item.replicate}"
                folder.mkdir(parents=True, exist_ok=True)
                np.savez(folder / "series.npz", values=values, frames=item.frames, times=item.times)
                record = {
                    **base,
                    "condition": condition,
                    "replicate": item.replicate,
                    "inputs": [_file_record(entry.path / "record.json") for entry in inputs],
                    "versions": _versions(),
                }
                (folder / "record.json").write_text(json.dumps(record, indent=1))
                series[condition].append(
                    ReplicateSeries(
                        condition, item.replicate, values, item.frames, item.times, folder
                    )
                )
        return Timeseries(name, unit, self.study, series, root, bounds)

    def plot(
        self,
        output_dir: str | Path | None = None,
        name: str | None = None,
        plot_settings: Any = None,
    ) -> Path:
        """Draw every replicate's series against time, with each condition's mean.

        See :func:`polyzymd.analyses.figures.plot_timeseries`. The figure goes
        to ``<output_dir>/<name>.<format>``; ``output_dir`` defaults to the
        ``figures`` folder next to ``polyzymd_results`` and ``name`` to
        ``<name>_timeseries``. ``plot_settings`` is a
        :class:`~polyzymd.config.analysis_settings.PlotSettings`, by default its
        defaults. Returns the path of the figure file.
        """
        from polyzymd.analyses.figures import plot_timeseries

        folder = output_dir or self.path.parent.parent / "figures"
        return plot_timeseries(self, folder, name or f"{self.name}_timeseries", plot_settings)

    def plot_distribution(
        self,
        threshold: float | None = None,
        output_dir: str | Path | None = None,
        name: str | None = None,
        title: str | None = None,
        plot_settings: Any = None,
    ) -> Path:
        """Draw the pooled and per-replicate distribution of the values of each condition.

        See :func:`polyzymd.analyses.figures.plot_distribution`. ``threshold``,
        in the unit of the series, is drawn as a vertical line. ``name``
        defaults to ``<name>_distribution`` and ``title``, which also labels
        the x axis, to the series name; the other arguments are those of
        :meth:`plot`. Returns the path of the figure file.
        """
        from polyzymd.analyses.figures import plot_distribution

        folder = output_dir or self.path.parent.parent / "figures"
        name = name or f"{self.name}_distribution"
        return plot_distribution(self, folder, name, threshold, title, plot_settings)

    def reduce(
        self,
        how: str | Callable = "mean",
        *,
        unit: Any = ...,
        bounds: Any = ...,
        detect_equilibration: bool = True,
    ) -> ReplicateValues:
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
        bounds : tuple of (float or None, float or None), optional
            Lowest and highest value the reduced quantity can take, used to
            warn when a 95 percent interval extends past them. Defaults to
            the series bounds for ``"mean"``, ``(0, 1)`` for ``"fraction"``,
            ``(0, None)`` for ``"std"`` and no bounds for a callable.
        detect_equilibration : bool, optional
            Find the start of the equilibrated region of each series with
            pymbar, as a diagnostic that changes no value. ``False`` skips it.

        Returns
        -------
        ReplicateValues
            One value per replicate, with the statistical inefficiency and
            effective sample size of each series.
        """
        import numpy as np

        from polyzymd.analyses.shared.autocorrelation import detect_equilibration as detect_start
        from polyzymd.analyses.shared.autocorrelation import n_effective, statistical_inefficiency

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
                start = detect_start(item.values) if detect_equilibration else None
                rows[condition].append(
                    (
                        item.replicate,
                        float(reducer(item.values, item.times)),
                        g,
                        n_effective(len(item.values), g),
                        len(item.values),
                        start,
                        None if start is None else float(item.times[start]),
                    )
                )
        if unit is ...:
            unit = None if how == "fraction" else self.unit
        if bounds is ...:
            defaults = {"mean": self.bounds, "fraction": (0.0, 1.0), "std": (0.0, None)}
            bounds = (None, None) if callable(how) else defaults[how]
        values = ReplicateValues(self, f"{label}_{self.name}", unit, label == "fraction", rows)
        values.bounds = tuple(bounds)
        return values


class ReplicateValues:
    """One value, or one labelled array, per replicate of every condition.

    ``labels`` is ``None`` for one number per replicate. Otherwise each
    replicate's value is an array with one entry per label, in the order of
    ``labels``, and every summary and comparison is made per label.
    """

    def __init__(
        self,
        source: Timeseries | Source,
        metric: str,
        unit: str | None,
        is_fraction: bool,
        rows: dict[str, list[tuple]],
        labels: list | None = None,
    ) -> None:
        self.source, self.metric, self.unit = source, metric, unit
        self.is_fraction, self.rows, self.labels = is_fraction, rows, labels
        self.bounds: tuple[float | None, float | None] = (0.0, 1.0) if is_fraction else (None, None)
        #: One note per replicate whose missing labels were filled with ``missing``.
        self.filled: list[str] = []

    @property
    def values(self) -> dict[str, list]:
        """Each condition's replicate values, in replicate order, as numbers or arrays."""
        return {label: [row[1] for row in rows] for label, rows in self.rows.items()}

    def _column(self, position: int | None) -> dict[str, list[float]]:
        """Each condition's replicate values of the label at ``position``, or the values."""
        if position is None:
            return self.values
        return {
            label: [float(row[1][position]) for row in rows] for label, rows in self.rows.items()
        }

    def _entries(self, untested: Sequence = ()) -> list[tuple[int | None, str | None]]:
        """Return the position and text of each label to report, or one unlabelled entry."""
        if self.labels is None:
            return [(None, None)]
        skip = {str(key) for key in untested}
        return [(i, str(key)) for i, key in enumerate(self.labels) if str(key) not in skip]

    def over_labels(
        self,
        how: str | Callable = "mean",
        metric: str | None = None,
        labels: Sequence | None = None,
    ) -> ReplicateValues:
        """Turn each replicate's labelled array, or the entries of some labels, into one number.

        Parameters
        ----------
        how : "mean" or callable
            ``"mean"`` averages the entries of each replicate. A callable
            receives one replicate's array in the order of :attr:`labels`.
        metric : str, optional
            Name of the new values. Defaults to ``<how>_<source name>``.
            ``"mean"`` keeps the bounds of these values, a callable has none.
        labels : sequence, optional
            Labels whose entries ``how`` receives, in this order. Defaults to
            every label.

        Returns
        -------
        ReplicateValues
            One number per replicate.
        """
        import numpy as np

        function = np.mean if how == "mean" else how
        if self.labels is None or not callable(function):
            raise ProtocolError(
                f"{self.metric}: over_labels needs labelled values and 'mean' or a function.",
                hint="Use over_labels('mean') on a result computed with labels=.",
            )
        name = how if isinstance(how, str) else getattr(how, "__name__", "reduced")
        known = {str(key): i for i, key in enumerate(self.labels)}
        chosen = list(self.labels) if labels is None else list(labels)
        unknown = [key for key in chosen if str(key) not in known]
        if unknown or not chosen:
            raise ProtocolError(
                f"{self.metric}: labels {unknown or chosen} are not labels of this result.",
                hint="Pass labels that the result has, such as residue IDs of its selection.",
            )
        where = [known[str(key)] for key in chosen]
        rows = {
            label: [(row[0], float(function(row[1][where])), *row[2:]) for row in items]
            for label, items in self.rows.items()
        }
        metric = metric or f"{name}_{self.metric}"
        values = ReplicateValues(self.source, metric, self.unit, False, rows)
        values.bounds = self.bounds if how == "mean" else (None, None)
        return values

    def plot(
        self,
        output_dir: str | Path | None = None,
        name: str | None = None,
        title: str | None = None,
        plot_settings: Any = None,
        highlight: Sequence = (),
        xlabel: str = "label",
    ) -> Path:
        """Draw each condition's mean with its 95 percent interval and every replicate value.

        See :func:`polyzymd.analyses.figures.plot_condition_values`, or for
        labelled values :func:`polyzymd.analyses.figures.plot_profile`, which
        marks the ``highlight`` labels and names the x axis ``xlabel``. The
        figure goes to ``<output_dir>/<name>.<format>``; ``output_dir``
        defaults to the ``figures`` folder next to ``polyzymd_results``,
        ``name`` to ``<source name>_<metric>`` and ``title`` to the source
        name. Returns the path of the figure file.
        """
        from polyzymd.analyses.figures import plot_condition_values, plot_profile

        folder = output_dir or self.source.path.parent.parent / "figures"
        name = name or f"{self.source.name}_{self.metric}"
        if self.labels is not None:
            return plot_profile(self, folder, name, title, plot_settings, highlight, xlabel)
        return plot_condition_values(self, folder, name, title, plot_settings)

    def summary(self, conditions: Sequence[str] | None = None) -> ProtocolReport:
        """Give each condition's n, mean, standard error and 95 percent interval.

        The interval is the mean plus or minus the Student t coverage factor
        for n replicates times the standard error, from
        :func:`~polyzymd.analyses.shared.statistics.mean_sem_ci`. Every
        replicate value is listed, with the statistical inefficiency and
        effective sample size of its series when it came from one. Labelled
        values give one row per condition and label, with the label in
        ``entry``.

        Parameters
        ----------
        conditions : sequence of str, optional
            Conditions to include. Defaults to every condition of the study.

        Returns
        -------
        ProtocolReport
            The condition rows and no comparisons.
        """
        return self._report(self._chosen(conditions), [])

    def compare(
        self,
        control: str | Mapping[str, Any] | None = None,
        conditions: Sequence[str] | None = None,
        test: str = "welch",
        untested: Sequence = (),
        within: str | Sequence[str] | None = None,
    ) -> ProtocolReport:
        """Compare every condition against the control, per label for labelled values.

        Without ``within`` every condition is compared with one control.
        With ``within``, one factor name or a list, the conditions that share
        their values of those factors form one stratum, and each condition is
        compared with the control of its own stratum
        (:func:`~polyzymd.analyses.study_statistics.comparison_pairs`); each
        row names its control in ``a`` and its stratum in ``stratum``. The
        factors are those of the study's ``study.yaml``. ``within`` and
        ``control`` default to its ``comparison:`` block, and ``within``
        given without ``control`` takes the block's control too;
        ``within=[]`` turns the strata off. The control of a stratum is never
        replaced by another condition: a control without replicate values,
        or not among the conditions, is refused.

        Each row gives ``mean(b) - mean(a)`` with its 95 percent interval from
        the same t test, the p value, the Benjamini-Hochberg adjusted p value,
        and Cohen's d and Hedges' g oriented like the difference. A row is not
        testable when a condition has fewer than two replicates, or when both
        conditions have the same value in every replicate; it then takes no
        part in the correction.

        The correction family is every tested row of this call. For one
        number per replicate that is the conditions compared with their
        control, over every stratum with ``within``, because together they
        answer one question, whether the condition changes the metric in any
        stratum, and a family is the set of tests behind one conclusion
        (Bender and Lange 2001, doi:10.1016/S0895-4356(00)00314-0).
        For labelled values it is every label of every condition
        compared, because scanning a profile for the labels that changed is a
        selective search over that set, and a false discovery rate holds for
        the discoveries of a search only when the family is the set searched
        (Benjamini 2010, doi:10.1111/j.1467-9868.2010.00746.x). Labels in
        ``untested`` are summarised but take no part in the tests or the
        correction.

        Parameters
        ----------
        control : str or mapping, optional
            Control condition. Defaults to the first condition of the study,
            or with ``within`` to the control of the ``comparison:`` block.
            With ``within``, the control of each stratum is the condition
            whose other factors equal this condition's, or, for a mapping
            such as ``{"polymer": "none"}``, whose factors have those values.
        conditions : sequence of str, optional
            Conditions to include, with their controls. Defaults to all.
        test : {"welch", "student"}, optional
            Welch's t test (default) or Student's t test.
        untested : sequence, optional
            Labels left out of the tests, for example one whose value the
            others fix.
        within : str or sequence of str, optional
            Factors whose values form each stratum, such as
            ``"temperature_K"``.

        Returns
        -------
        ProtocolReport
            The condition rows and one comparison row per non-control
            condition and tested label.

        Raises
        ------
        ProtocolError
            If a stratum has no control or more than one, a control or a
            compared condition has no replicate values, a control is not
            among the conditions, or a condition does not declare a
            ``within`` factor.
        """
        chosen = self._chosen(conditions)
        return self._compare(chosen, self._pairs(chosen, control, within), test, untested)

    def _pairs(
        self,
        chosen: list[str],
        control: str | Mapping[str, Any] | None = None,
        within: str | Sequence[str] | None = None,
    ) -> list[tuple[str, str, dict[str, Any] | None]]:
        """Return ``(control, condition, stratum)`` for each condition of ``chosen``; see :meth:`compare`.

        Controls are looked up among every condition of the study, also those
        with no replicate values, so a missing control is never replaced.
        """
        from polyzymd.analyses.study_statistics import comparison_pairs

        study = self.source.study
        default = getattr(study, "comparison", None) or {}
        if within is None:
            within = default.get("within") or []
        within = [within] if isinstance(within, str) else list(within)
        if within:
            control = control or default.get("control")
        factors = getattr(study, "factors", None) or {}
        # The study file's factors also name the conditions left out with --label.
        labels = list(dict.fromkeys([*self.rows, *factors]))
        return comparison_pairs(chosen, labels, factors, within, control or study.control)

    def _compare(
        self,
        chosen: list[str],
        pairs: list[tuple[str, str, dict[str, Any] | None]],
        test: str = "welch",
        untested: Sequence = (),
    ) -> ProtocolReport:
        """Report the conditions of ``chosen`` and their controls, and test ``pairs``."""
        from polyzymd.analyses.protocols import PairwiseReport, _difference_ci
        from polyzymd.analyses.shared.inferential_statistics import (
            NO_SIGNIFICANT_CHANGE,
            benjamini_hochberg,
            cohens_d,
            independent_ttest,
        )

        problems = {}
        for control, label, stratum in pairs:
            where = ", ".join(f"{name} {value}" for name, value in (stratum or {}).items())
            whose = f"the stratum {where} has no control values: its" if where else "the"
            if control not in self.rows:
                problems[control] = f"{whose} control {control} is not among the conditions"
            elif not self.rows[control]:
                problems[control] = f"{whose} control {control} has no replicate values"
            if not self.rows[label]:
                problems[label] = f"condition {label} has no replicate values"
        if problems:
            raise ProtocolError(
                f"Cannot compare: {'; '.join(problems.values())}.",
                hint="Give conditions= the conditions with replicate values, and include the "
                "control of each (with polyzymd analyze --label, add a --label for it). A "
                "condition is compared only with its own control.",
            )
        controls = dict.fromkeys(pair[0] for pair in pairs)
        chosen = [label for label in controls if label not in chosen] + chosen
        if test not in ("welch", "student"):
            raise ProtocolError(f"Unknown test {test!r}.", hint="Use 'welch' or 'student'.")
        rows = []
        for position, entry in self._entries(untested):
            values = self._column(position)
            for control, label, stratum in pairs:
                a, b = values[control], values[label]
                testable = min(len(a), len(b)) >= 2 and len(set(a)) + len(set(b)) > 2
                effect = cohens_d(b, a)
                rows.append(
                    PairwiseReport(
                        a=control,
                        b=label,
                        stratum=stratum,
                        entry=entry,
                        delta=sum(b) / len(b) - sum(a) / len(a),
                        delta_ci95=_difference_ci(a, b, f"{test}_t"),
                        p=independent_ttest(b, a, method=test).p_value if testable else None,
                        test=f"{test}_t",
                        cohens_d=_finite(effect.cohens_d),
                        hedges_g=_finite(effect.hedges_g),
                        testable=testable,
                    )
                )
        # One Benjamini-Hochberg family per call: the conditions compared with
        # their control, over every stratum, for every label compared. A
        # family is the set of tests behind one conclusion (Bender and Lange 2001), FDR controlled in
        # separate families stays controlled overall (Benjamini and Yekutieli
        # 2001), and a scan of labels for changes is corrected over the labels
        # scanned (Benjamini 2010).
        family_size = sum(row.p is not None and math.isfinite(row.p) for row in rows)
        for row, corrected in zip(rows, benjamini_hochberg([row.p for row in rows]), strict=True):
            row.family_size = family_size if row.p is not None and math.isfinite(row.p) else None
            row.p_adjusted = corrected.adjusted_p_value
            row.significant = corrected.significant
            row.direction = (
                ("increased" if row.delta > 0 else "decreased")
                if row.significant
                else NO_SIGNIFICANT_CHANGE
            )
        report = self._report(chosen, rows)
        left_out = [str(key) for key in untested if str(key) in map(str, self.labels or [])]
        if left_out:
            report.warnings.append(
                f"labels {', '.join(left_out)} are summarised but left out of the tests and the "
                "correction, as untested= asks"
            )
        return report

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
        from polyzymd.analyses.study_freeze import group_warnings

        conditions, notes = [], list(self.filled)
        for position, entry in self._entries():
            at = "" if entry is None else f" at {entry}"
            values = self._column(position)
            for label in chosen:
                rows = self.rows[label]
                item = _condition(label, values[label])
                item.entry = entry
                item.replicates = [row[0] for row in rows]
                if rows and rows[0][2] is not None:
                    item.statistical_inefficiency = [row[2] for row in rows]
                    item.n_effective = [row[3] for row in rows]
                checked = [row for row in rows if row[5] is not None]
                item.eq_detected_frame = [row[5] for row in checked]
                item.eq_detected_ns = [row[6] for row in checked]
                few = [str(row[0]) for row in checked if row[3] < EQUILIBRATION_MIN_N_EFFECTIVE]
                if few:
                    notes.append(
                        f"condition {label}: replicates {', '.join(few)} have fewer than "
                        f"{EQUILIBRATION_MIN_N_EFFECTIVE} effective samples, so the start of an "
                        "equilibrated region cannot be detected reliably; values and statistics "
                        "are unaffected"
                    )
                for row in checked:
                    if (
                        row[3] >= EQUILIBRATION_MIN_N_EFFECTIVE
                        and row[5] > EQUILIBRATION_WARNING_FRACTION * row[4]
                    ):
                        notes.append(
                            f"condition {label} replicate {row[0]}: pymbar detect_equilibration "
                            f"puts the start of the equilibrated region at {row[6]:.4g} ns, "
                            f"production frame {row[5] + 1} of {row[4]}, after the equilibration "
                            "window; the window may be too short for it"
                        )
                if len(rows) > 1 and len(set(values[label])) == 1:
                    item.ci95, item.ci_method = None, "not_estimable"
                    notes.append(
                        f"condition {label} has the same {self.metric}{at} in every replicate, "
                        "so its interval is not estimable"
                    )
                elif len(rows) < 2 and position in (None, 0):
                    notes.append(f"condition {label} has one replicate, so it has no interval")
                low, high = self.bounds
                if item.ci95 and (
                    (low is not None and item.ci95[0] < low)
                    or (high is not None and item.ci95[1] > high)
                ):
                    # Grossfield et al. (2018): a bounded quantity is not Gaussian,
                    # so a t interval that crosses the bound is not reliable.
                    notes.append(
                        f"the 95 percent interval of condition {label}{at} extends past the "
                        f"bounds {'-inf' if low is None else format(low, 'g')} to "
                        f"{'inf' if high is None else format(high, 'g')} of {self.metric}, "
                        "where a t interval is not reliable"
                    )
                conditions.append(item)
        untestable: dict[tuple[str, str], list] = {}
        for row in pairwise:
            untestable.setdefault((row.a, row.b), [])
            if not row.testable:
                untestable[(row.a, row.b)].append(row)
        for (a, b), rows in untestable.items():
            total = sum(1 for row in pairwise if (row.a, row.b) == (a, b))
            if rows and len(rows) == total and total > 1:
                # One note for the pair; its rows stay in the report's pairwise table.
                notes.append(
                    f"{a} vs {b} is not testable at any of its {total} labels: a condition has "
                    "fewer than two replicates, or both have one value in every replicate"
                )
                continue
            for row in rows:
                at = "" if row.entry is None else f" at {row.entry}"
                notes.append(
                    f"{row.a} vs {row.b}{at} is not testable: a condition has fewer than two "
                    "replicates, or both have one value in every replicate"
                )
        study, versions = self.source.study, _versions()
        for label in chosen:
            # What loading found about each replicate's files, such as segments
            # missing on disk, belongs in the report, not only in the log.
            condition = study[label]
            provider = getattr(condition, "_provider", None)
            for replicate in getattr(condition, "replicates", []):
                if provider is None:
                    break
                for text in provider.provenance_for(replicate.index).warnings:
                    note = f"condition {label} replicate {replicate.index}: {text}"
                    if note not in notes:
                        notes.append(note)
        if len(chosen) > 1:
            from polyzymd.analyses.study import production_length_warnings

            notes.extend(production_length_warnings(study, chosen))
        # One line for a warning that several conditions share, as freeze writes it.
        notes = group_warnings(notes, [f"condition {label}" for label in chosen])
        if self.labels is None:
            verdict = _verdict(self.metric, self.unit, conditions, pairwise)
        else:
            verdict = _labelled_verdict(self.metric, self.unit, chosen, conditions, pairwise)
        return ProtocolReport(
            analysis=self.source.name,
            protocol_version="2",
            metric=self.metric,
            unit=self.unit,
            equilibration=study[chosen[0]].equilibration,
            stride=getattr(study[chosen[0]], "stride", 1),
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
            verdict=verdict,
        )


def _labelled_verdict(
    metric: str, unit: str | None, chosen: list[str], conditions: list, pairwise: list
) -> list[str]:
    """Write one sentence per condition, or per comparison, naming every significant label."""
    from polyzymd.analyses.protocols import VERDICT_NOT_TESTABLE, _num

    unit_text = f" {unit}" if unit else ""
    if not pairwise:
        sentences = []
        for label in chosen:
            means = [item.mean for item in conditions if item.label == label]
            finite = [mean for mean in means if math.isfinite(mean)]
            if not finite:
                sentences.append(f"{label} {metric} over {len(means)} labels, no finite label mean")
                continue
            sentences.append(
                f"{label} {metric} over {len(means)} labels, label means from "
                f"{_num(min(finite))} to {_num(max(finite))}{unit_text}"
                + (
                    f" ({len(means) - len(finite)} labels without a mean)"
                    if len(finite) < len(means)
                    else ""
                )
            )
        return sentences
    replicates = {item.label: item.n_replicates for item in conditions}
    sentences = []
    for label in dict.fromkeys(row.b for row in pairwise):
        rows = [row for row in pairwise if row.b == label]
        tested = [row for row in rows if row.p_adjusted is not None]
        if not tested:
            n_text = f"n {replicates.get(rows[0].a, 0)} vs {replicates.get(label, 0)}"
            sentences.append(
                f"{VERDICT_NOT_TESTABLE}: {metric} for {rows[0].a} vs {label} needs at least two "
                f"replicates per condition and values that vary ({n_text})"
            )
            continue
        larger = [row.entry for row in tested if row.significant and row.delta > 0]
        smaller = [row.entry for row in tested if row.significant and row.delta < 0]
        family = tested[0].family_size if tested else 0
        sentences.append(
            f"{label} vs {rows[0].a}: {len(larger) + len(smaller)} of {len(tested)} tested labels "
            f"differ in {metric} after BH over a family of {family}; larger at "
            f"{', '.join(larger) or 'none'}; smaller at {', '.join(smaller) or 'none'}"
        )
    return sentences


def _finite(value: float) -> float | None:
    """Return ``value``, or ``None`` when it is NaN or infinite."""
    return value if math.isfinite(value) else None
