"""Run a study's own analysis functions, named by file and function in ``study.yaml``.

An ``analyses:`` entry with ``function: analyses/lid.py:lid_distance`` runs
that function through :meth:`~polyzymd.analyses.study.Study.per_replicate`
or :meth:`~polyzymd.analyses.study.Study.timeseries`, with the same stored
records, statistics, report and figures as a shipped analysis. Its stored
results are keyed on the hash of the whole file, not only the function's
source, so editing a helper the function calls recomputes them.
"""

from __future__ import annotations

import sys
import types
from pathlib import Path
from typing import Any, Callable

from polyzymd.analyses.exceptions import ProtocolError

#: Attribute that tells the stored record to hash the function's whole file.
MODULE_FILE_ATTRIBUTE = "__polyzymd_module_file__"


def load_function(file: Path, qualname: str) -> Callable:
    """Import ``qualname`` from the Python file ``file`` and return it.

    The file's folder is put first on ``sys.path`` while it imports, so it
    can import helper modules beside it. The function is marked so its
    stored results are keyed on the hash of the whole file.

    Raises
    ------
    ProtocolError
        If the file cannot be imported or has no such callable.
    """
    file = Path(file).resolve()
    # The record keeps the module name, so it must not depend on where the
    # study folder sits; the module is not added to sys.modules.
    module_name = f"polyzymd_study.{file.stem}"
    module = types.ModuleType(module_name)
    module.__file__ = str(file)
    sys.path.insert(0, str(file.parent))
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
    function: Any = module
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


def run_user_analysis(
    study: Any,
    run: str,
    user: Any,
    *,
    settings: dict[str, Any] | None = None,
    output_dir: Path,
    recompute: bool = False,
    plots: bool = True,
) -> Any:
    """Run a study's user function over every replicate and return its ProtocolReport.

    ``user`` is the :class:`~polyzymd.analyses.study_file.UserFunction` of
    the entry ``run``; ``settings`` (from ``--set``) override its
    ``settings``. With one condition the report summarises it, and with
    several it compares each with the first.
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
    if user.kind == "timeseries":
        series = study.timeseries(
            function,
            unit=user.unit,
            name=run,
            recompute=recompute,
            output_dir=output_dir,
            **kwargs,
        )
        if plots:
            series.plot()
        values = series.reduce(user.reduce)
    else:
        values = study.per_replicate(
            function,
            unit=user.unit,
            labels=user.labels,
            name=run,
            recompute=recompute,
            output_dir=output_dir,
            **kwargs,
        )
    if plots:
        values.plot()
    return values.compare() if len(study) > 1 else values.summary()
