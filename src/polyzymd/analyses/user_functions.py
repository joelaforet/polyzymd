"""Run a study's own analysis functions, named by file and function in ``study.yaml``.

An ``analyses:`` entry with ``function: analyses/lid.py:lid_distance`` runs
that function through :meth:`~polyzymd.analyses.study.Study.per_replicate`
or :meth:`~polyzymd.analyses.study.Study.timeseries`, with the same stored
records, statistics, report and figures as a shipped analysis. Its stored
results are keyed on the hash of the whole file, not only the function's
source, so editing a helper the function calls recomputes them.
"""

from __future__ import annotations

import shutil
import sys
import tempfile
import types
from pathlib import Path
from typing import Any, Callable

from polyzymd.analyses.exceptions import ProtocolError

#: Attribute that tells the stored record to hash the function's whole file.
MODULE_FILE_ATTRIBUTE = "__polyzymd_module_file__"
#: Folders that function files were loaded from in this process.
_FOLDERS: set[Path] = set()


def load_function(file: Path, qualname: str) -> Callable:
    """Import ``qualname`` from the Python file ``file`` and return it.

    The file is imported by :func:`load_module`. The function is marked so
    its stored results are keyed on the hash of the whole file.

    Raises
    ------
    ProtocolError
        If the file cannot be imported or has no such callable.
    """
    file = Path(file).resolve()
    function: Any = load_module(file)
    for part in qualname.split("."):
        function = getattr(function, part, None)
    if not callable(function):
        raise ProtocolError(
            f"{file} has no function {qualname!r}.",
            hint="Name a function defined in the file, as 'file.py:function_name'.",
        )
    try:
        setattr(function, MODULE_FILE_ATTRIBUTE, str(file))
    except (AttributeError, TypeError):
        pass
    return function


def _owner_name(file: Path) -> str:
    """Return the name of the study or project folder that holds ``file``, as an identifier.

    The nearest folder above ``file`` with a ``study.yaml`` or ``project.yaml``,
    or the file's own folder when there is none.
    """
    import re

    owner = file.parent
    for folder in file.parents:
        if (folder / "study.yaml").is_file() or (folder / "project.yaml").is_file():
            owner = folder
            break
    name = re.sub(r"\W", "_", owner.name) or "study"
    return f"_{name}" if name[0].isdigit() else name


def load_module(file: Path) -> types.ModuleType:
    """Import the Python file ``file`` from its current text and return the module.

    The file's folder is put first on ``sys.path`` while it imports, so it
    can import helper modules beside it; helpers imported earlier from that
    folder or another study's are imported again, so edits since then take
    effect and a ``helper.py`` of another study is never used.

    Raises
    ------
    ProtocolError
        If the file cannot be imported.
    """
    file = Path(file).resolve()
    # The record keeps the module name, so it must not depend on where the
    # study folder sits; it names the study or project folder the file
    # belongs to, so two studies' f.py are two modules. The module and its
    # parent packages are put in sys.modules, where dataclasses and pickle
    # look up the classes the file defines.
    package = f"polyzymd_study.{_owner_name(file)}"
    module_name = f"{package}.{file.stem}"
    module = types.ModuleType(module_name)
    module.__file__ = str(file)
    sys.modules.setdefault("polyzymd_study", types.ModuleType("polyzymd_study"))
    sys.modules.setdefault(package, types.ModuleType(package))
    sys.modules[module_name] = module
    # A helper module imported earlier in this process from this folder may
    # have been edited since, and one from another study's folder may have
    # the same name; drop both so the import reads this folder's current
    # file, whose content is what the stored results' hash covers. Only
    # modules imported through such a folder are dropped (helper.py,
    # util/k.py): an installed package that merely lies under it, in a
    # .pixi environment or a source checkout, is left alone.
    _FOLDERS.add(file.parent)
    for name, loaded in list(sys.modules.items()):
        if any(_imported_from(loaded, name, folder) for folder in _FOLDERS):
            del sys.modules[name]
    sys.path.insert(0, str(file.parent))
    # Helpers it imports are compiled into a fresh cache: a .pyc beside them is
    # trusted by modification second and size, so an edit of the same size
    # within one second would otherwise run the old helper.
    saved_prefix, sys.pycache_prefix = sys.pycache_prefix, tempfile.mkdtemp(prefix="pzpyc")
    try:
        # Compile the file's current text rather than importing it, so a
        # cached .pyc that predates an edit is never run: the code that runs
        # is the code whose hash keys the stored results.
        exec(compile(file.read_bytes(), str(file), "exec"), module.__dict__)  # noqa: S102
    except Exception as exc:  # noqa: BLE001 - any error in user code is reported the same way
        raise ProtocolError(
            f"Importing {file} failed: {type(exc).__name__}: {exc}",
            hint="Run 'python FILE' to see the full error, and fix it in the file.",
        ) from exc
    finally:
        sys.path.remove(str(file.parent))
        shutil.rmtree(sys.pycache_prefix, ignore_errors=True)
        sys.pycache_prefix = saved_prefix
    return module


def _imported_from(module: Any, name: str, folder: Path) -> bool:
    """Return whether ``module`` was imported through ``folder`` on ``sys.path``.

    True when its top-level name is a file or package directly in ``folder``
    and its file lies there: ``helper`` from ``folder/helper.py``, or
    ``util.k`` from ``folder/util/k.py``.
    """
    origin = getattr(module, "__file__", None)
    if not origin:
        return False
    try:
        relative = Path(origin).resolve().relative_to(folder)
    except ValueError:
        return False
    top = name.split(".")[0]
    first = relative.parts[0]
    return first == top or first == f"{top}.py"


def run_user_analysis(
    study: Any,
    run: str,
    user: Any,
    *,
    settings: dict[str, Any] | None = None,
    output_dir: Path,
    recompute: bool = False,
    plots: bool = True,
    part: str | None = None,
) -> Any:
    """Run a study's user function over every replicate and return its ProtocolReport.

    ``user`` is the :class:`~polyzymd.analyses.study_file.UserFunction` of
    the entry ``run``; ``settings`` (from ``--set``) override its
    ``settings``. With one condition the report summarises it, and with
    several it compares each with the first. With ``parts``, every part is
    stored and plotted, and the report covers ``part`` (``--run``), the first
    part by default, naming the others in ``all_runs``.
    """
    from polyzymd.analyses.timeseries import select, universe

    function = load_function(user.file, user.qualname)
    # With allow_empty, a selection matching no atoms reaches the function as
    # an empty AtomGroup, so a no-polymer control is measured, not dropped.
    kwargs: dict[str, Any] = {
        key: select(value, allow_empty=user.allow_empty) for key, value in user.selections.items()
    }
    if user.universe:
        kwargs[user.universe] = universe()
    kwargs.update(user.settings)
    kwargs.update(settings or {})
    parts = user.parts
    if part is not None and part not in (parts or []):
        raise ProtocolError(
            f"{run} has no part {part!r}.",
            hint=f"Use --run with one of {', '.join(parts)}."
            if parts
            else "Leave --run out; this function measures one quantity.",
        )
    extra = {} if parts is None else {"parts": parts}
    if user.kind == "timeseries":
        measured = study.timeseries(
            function,
            unit=user.unit,
            name=run,
            recompute=recompute,
            output_dir=output_dir,
            **extra,
            **kwargs,
        )
        every = measured if parts is not None else {None: measured}
        if plots:
            for series in every.values():
                series.plot()
        results = {key: series.reduce(user.reduce) for key, series in every.items()}
    else:
        measured = study.per_replicate(
            function,
            unit=user.unit,
            labels=user.labels,
            missing=user.missing,
            note_filled=True,
            name=run,
            recompute=recompute,
            output_dir=output_dir,
            **extra,
            **kwargs,
        )
        results = measured if parts is not None else {None: measured}
    if plots:
        for values in results.values():
            values.plot()
    chosen = part or (parts[0] if parts else None)
    values = results[chosen]
    report = values.compare() if len(study) > 1 else values.summary()
    if parts is None:
        return report
    return report.model_copy(update={"analysis": run, "run": chosen, "all_runs": list(parts)})
