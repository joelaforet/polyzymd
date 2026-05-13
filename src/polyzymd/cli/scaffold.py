"""Scaffold generator for new analysis plugins.

Creates the minimal file set needed for a new analysis plugin using Jinja
templates stored as package resources. The default scaffold is a single-file
measurement plugin. Advanced scaffolds keep the runner-backed package layout:
lifecycle wiring lives in ``__init__.py`` while trajectory logic lives in
``_runner.py``. Pydantic style scaffolds also place result models in
``_results.py``.

Usage::

    polyzymd new-analysis my_analysis
    polyzymd new-analysis my_analysis --advanced
    polyzymd new-analysis my_analysis --style pydantic
    polyzymd new-analysis my_analysis --dry-run
    polyzymd new-analysis my_analysis --force
"""

from __future__ import annotations

import keyword
import re
from pathlib import Path

from polyzymd.cli._scaffold.models import ADVANCED_STYLES, VALID_STYLES, ScaffoldSpec
from polyzymd.cli._scaffold.renderer import render_scaffold

# ---------------------------------------------------------------------------
# Name helpers
# ---------------------------------------------------------------------------

_NAME_RE = re.compile(r"^[a-z][a-z0-9_]*$")

_RESERVED_NAMES = frozenset(
    {
        "base",
        "cli",
        "config",
        "discovery",
        "orchestrator",
        "runner",
        "shared",
        "stats",
    }
)


def validate_name(name: str, *, check_existing: bool = True) -> str | None:
    """Return an error message if *name* is invalid, otherwise ``None``.

    Parameters
    ----------
    name : str
        Proposed plugin name in snake_case.
    check_existing : bool, optional
        If True, also reject names that collide with already-registered
        analysis plugins or their aliases, by default True.

    Returns
    -------
    str or None
        Validation error text, or ``None`` when the name is valid.
    """
    if not _NAME_RE.match(name):
        return (
            f"'{name}' is not a valid plugin name. Use lowercase snake_case (e.g. 'my_analysis')."
        )
    if keyword.iskeyword(name):
        return f"'{name}' is a Python keyword."
    if name in _RESERVED_NAMES:
        return f"'{name}' is reserved for framework infrastructure."
    if check_existing:
        try:
            from polyzymd.analyses.discovery import list_all_names

            if name in list_all_names():
                return f"'{name}' already exists as a registered analysis plugin."
        except ModuleNotFoundError as exc:
            # Only suppress if the discovery module itself is unavailable
            if exc.name not in ("polyzymd.analyses", "polyzymd.analyses.discovery"):
                raise
    return None


def validate_class_name(class_name: str) -> str | None:
    """Return an error message if *class_name* is invalid, otherwise ``None``.

    Parameters
    ----------
    class_name : str
        Proposed PascalCase class prefix, for example ``SolventShell``.

    Returns
    -------
    str or None
        Validation error text, or ``None`` when the class prefix is valid.
    """
    if not class_name.isidentifier():
        return (
            f"'{class_name}' is not a valid Python identifier. "
            "Use PascalCase (e.g. 'SolventShell')."
        )
    if keyword.iskeyword(class_name):
        return f"'{class_name}' is a Python keyword and cannot be used as a class name."
    if not class_name[0].isupper():
        return f"'{class_name}' should start with an uppercase letter (PascalCase convention)."
    return None


def to_pascal_case(snake: str) -> str:
    """Convert a snake_case name to PascalCase.

    Parameters
    ----------
    snake : str
        Snake-case plugin name.

    Returns
    -------
    str
        PascalCase class prefix.
    """
    return "".join(part.capitalize() for part in snake.split("_"))


def _build_spec(
    name: str, class_name: str | None, style: str, *, advanced: bool = False
) -> ScaffoldSpec:
    """Validate user inputs and build a scaffold render specification.

    Parameters
    ----------
    name : str
        Plugin name in snake_case.
    class_name : str or None
        Optional PascalCase class prefix.
    style : str
        Scaffold style name.
    advanced : bool, optional
        If True, use the runner-backed package scaffold. When no explicit
        advanced style is requested, this selects the legacy ``"dict"`` style.

    Returns
    -------
    ScaffoldSpec
        Validated scaffold rendering specification.

    Raises
    ------
    ValueError
        If any input is invalid.
    """
    effective_style = "dict" if advanced and style == "measurement" else style
    if effective_style not in VALID_STYLES:
        raise ValueError(f"Invalid style '{style}'. Choose from: {', '.join(VALID_STYLES)}")
    if advanced and effective_style not in ADVANCED_STYLES:
        raise ValueError(
            "Advanced scaffolds require style 'dict' or 'pydantic', " f"got '{effective_style}'."
        )

    name_error = validate_name(name, check_existing=True)
    if name_error:
        raise ValueError(name_error)

    cls = class_name or to_pascal_case(name)
    cls_error = validate_class_name(cls)
    if cls_error:
        raise ValueError(cls_error)

    return ScaffoldSpec(name=name, class_name=cls, style=effective_style)


def _check_layout_conflicts(spec: ScaffoldSpec, project_root: Path) -> None:
    """Reject single-file/package layout collisions for a scaffold.

    Parameters
    ----------
    spec : ScaffoldSpec
        Validated scaffold rendering specification.
    project_root : Path
        Repository root directory containing ``src/`` and ``tests/``.

    Raises
    ------
    FileExistsError
        If the opposite plugin layout already exists.
    """
    analyses_root = project_root / "src" / "polyzymd" / "analyses"
    module_path = analyses_root / f"{spec.name}.py"
    package_path = analyses_root / spec.name

    if spec.uses_measurement_api and package_path.exists():
        raise FileExistsError(
            f"Cannot create {module_path}: package layout {package_path} already exists."
        )
    if not spec.uses_measurement_api and module_path.exists():
        raise FileExistsError(
            f"Cannot create {package_path}: single-file layout {module_path} already exists."
        )


def generate_scaffold(
    name: str,
    project_root: Path,
    *,
    class_name: str | None = None,
    style: str = "measurement",
    advanced: bool = False,
    force: bool = False,
    dry_run: bool = False,
) -> list[Path]:
    """Create scaffold files for a new analysis plugin.

    Parameters
    ----------
    name : str
        Plugin name in snake_case, for example ``"solvent_shell"``.
    project_root : Path
        Repository root directory containing ``src/`` and ``tests/``.
    class_name : str or None, optional
        PascalCase class prefix. Auto-derived from *name* when omitted, by
        default None.
    style : str, optional
        ``"measurement"`` for a single-file scalar measurement plugin,
        ``"dict"`` for advanced plain-dict results, or ``"pydantic"`` for
        advanced typed result models, by default ``"measurement"``.
    advanced : bool, optional
        Use the advanced runner-backed package scaffold. If ``style`` is left
        as ``"measurement"``, this selects the legacy ``"dict"`` scaffold, by
        default False.
    force : bool, optional
        Overwrite existing files, by default False.
    dry_run : bool, optional
        Return paths without writing files, by default False.

    Returns
    -------
    list[Path]
        Paths of created, or would-be-created, files.

    Raises
    ------
    FileExistsError
        If a target path exists and ``force`` is False.
    ValueError
        If the plugin name, class prefix, or style is invalid.
    """
    spec = _build_spec(name=name, class_name=class_name, style=style, advanced=advanced)
    _check_layout_conflicts(spec=spec, project_root=project_root)
    files = render_scaffold(spec=spec, project_root=project_root)

    created: list[Path] = []
    for path, content in files.items():
        if path.exists() and not force:
            raise FileExistsError(f"{path} already exists. Use --force to overwrite.")
        if not dry_run:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
        created.append(path)

    return created
