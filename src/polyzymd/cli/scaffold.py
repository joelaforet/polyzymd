"""Scaffold generator for new analysis plugins.

Creates the minimal file set needed for a new analysis plugin using Jinja
templates stored as package resources. The default scaffold is a single-file
MDAnalysis-native plugin. Advanced scaffolds generate packages that keep
MDAnalysis ``AnalysisBase`` code in a private, lazy-imported module.

Usage::

    polyzymd new-analysis my_analysis
    polyzymd new-analysis my_analysis --advanced
    polyzymd new-analysis my_analysis --style dict
    polyzymd new-analysis my_analysis --dry-run
    polyzymd new-analysis my_analysis --force
"""

from __future__ import annotations

import keyword
import re
from pathlib import Path

from polyzymd.cli._scaffold.models import ADVANCED_STYLES, DEFAULT_STYLE, VALID_STYLES, ScaffoldSpec
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
        "mda",
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
        Scaffold style name. ``"simple"`` requests the default single-file
        MDAnalysis-native scaffold; ``"dict"`` requests an advanced package
        scaffold using canonical artifact payloads.
    advanced : bool, optional
        If True, request an advanced package scaffold.

    Returns
    -------
    ScaffoldSpec
        Validated scaffold rendering specification.

    Raises
    ------
    ValueError
        If any input is invalid.
    """
    effective_style = "dict" if advanced and style == DEFAULT_STYLE else style
    if effective_style not in VALID_STYLES:
        raise ValueError(f"Invalid style '{style}'. Choose from: {', '.join(VALID_STYLES)}")

    name_error = validate_name(name, check_existing=False)
    if name_error:
        raise ValueError(name_error)

    cls = class_name or to_pascal_case(name)
    cls_error = validate_class_name(cls)
    if cls_error:
        raise ValueError(cls_error)

    return ScaffoldSpec(name=name, class_name=cls, style=effective_style)


def _format_paths(paths: list[Path]) -> str:
    """Return a concise, deterministic path list for errors.

    Parameters
    ----------
    paths : list[Path]
        Paths to format.

    Returns
    -------
    str
        Comma-separated path list.
    """
    return ", ".join(str(path) for path in sorted(paths))


def _source_target_paths(files: dict[Path, str], project_root: Path) -> set[Path]:
    """Return generated source paths that can define a plugin module.

    Parameters
    ----------
    files : dict[Path, str]
        Rendered scaffold files keyed by output path.
    project_root : Path
        Repository root for the scaffold operation.

    Returns
    -------
    set[Path]
        Resolved generated paths under ``src/polyzymd/analyses``.
    """
    analyses_root = (project_root / "src" / "polyzymd" / "analyses").resolve()
    paths: set[Path] = set()
    for path in files:
        resolved = path.resolve()
        try:
            resolved.relative_to(analyses_root)
        except ValueError:
            continue
        paths.add(resolved)
    return paths


def _path_has_scaffold_signature(path: Path) -> bool:
    """Return whether an existing file looks like scaffold output.

    Parameters
    ----------
    path : Path
        Existing source path to inspect.

    Returns
    -------
    bool
        True when the file contains scaffold-specific placeholder text.
    """
    if not path.exists() or not path.is_file():
        return False

    text = path.read_text(encoding="utf-8")
    # Keep legacy signature text so --force can still identify old scaffold output
    signatures = (
        "Replace this placeholder with domain-specific measurement logic",
        "Replace this placeholder with domain-specific MDAnalysis logic",
        "Generated MDAnalysis-native analysis plugin",
        "Generated advanced MDAnalysis-native analysis package",
        "Uses plain dicts for result containers",
        "Uses typed Pydantic result models for validation",
    )
    return any(signature in text for signature in signatures)


def _check_registered_name_conflict(
    spec: ScaffoldSpec,
    files: dict[Path, str],
    project_root: Path,
    *,
    force: bool,
) -> None:
    """Reject registered-name collisions unless force targets scaffold files.

    Parameters
    ----------
    spec : ScaffoldSpec
        Validated scaffold rendering specification.
    files : dict[Path, str]
        Rendered scaffold files keyed by output path.
    project_root : Path
        Repository root for the scaffold operation.
    force : bool
        Whether existing scaffold files may be overwritten.

    Raises
    ------
    ValueError
        If the requested name collides with a built-in, external, or alias-only
        analysis registration.
    """
    try:
        from polyzymd.analyses.discovery import get_analysis, list_all_names

        if spec.name not in list_all_names():
            return

        if not force:
            raise ValueError(f"'{spec.name}' already exists as a registered analysis plugin.")

        analysis_cls = get_analysis(spec.name)
    except ModuleNotFoundError as exc:
        if exc.name not in ("polyzymd.analyses", "polyzymd.analyses.discovery"):
            raise
        return

    if getattr(analysis_cls, "name", None) != spec.name:
        raise ValueError(
            f"'{spec.name}' already exists as an analysis alias and cannot be overwritten."
        )

    module_file = Path(
        __import__(analysis_cls.__module__, fromlist=["__file__"]).__file__
    ).resolve()
    source_paths = _source_target_paths(files=files, project_root=project_root)
    if module_file not in source_paths or not _path_has_scaffold_signature(module_file):
        raise ValueError(
            f"'{spec.name}' already exists as a registered analysis plugin and is not a "
            "matching scaffold target. Choose a different name."
        )


def _check_target_conflicts(files: dict[Path, str], *, force: bool) -> None:
    """Preflight all scaffold output paths before writing.

    Parameters
    ----------
    files : dict[Path, str]
        Rendered scaffold files keyed by output path.
    force : bool
        Whether existing files may be overwritten.

    Raises
    ------
    FileExistsError
        If any parent path is not a directory, any target path has the wrong
        type, or any target path exists and ``force`` is False.
    """
    parent_conflicts = sorted(
        {
            parent
            for path in files
            for parent in path.parents
            if parent.exists() and not parent.is_dir()
        }
    )
    if parent_conflicts:
        raise FileExistsError(
            "Scaffold parent path is not a directory: "
            f"{_format_paths(parent_conflicts)}. Resolve the path collision before scaffolding."
        )

    target_type_conflicts = [path for path in files if path.exists() and not path.is_file()]
    if target_type_conflicts:
        raise FileExistsError(
            "Scaffold target path exists but is not a file: "
            f"{_format_paths(target_type_conflicts)}. Resolve the path collision before scaffolding."
        )

    existing_files = [path for path in files if path.exists()]
    if existing_files and not force:
        raise FileExistsError(
            f"Scaffold target path already exists: {_format_paths(existing_files)}. "
            "Use --force to overwrite."
        )


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

    if spec.uses_single_file_layout and package_path.exists():
        raise FileExistsError(
            f"Cannot create {module_path}: package layout {package_path} already exists."
        )
    if spec.uses_package_layout and module_path.exists():
        raise FileExistsError(
            f"Cannot create {package_path}: single-file layout {module_path} already exists."
        )


def generate_scaffold(
    name: str,
    project_root: Path,
    *,
    class_name: str | None = None,
    style: str = DEFAULT_STYLE,
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
        ``"simple"`` for the default single-file MDAnalysis-native plugin, or
        ``"dict"`` for advanced canonical artifacts, by default ``"simple"``.
    advanced : bool, optional
        Request an advanced package scaffold, by default False.
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
    _check_registered_name_conflict(spec=spec, files=files, project_root=project_root, force=force)
    _check_target_conflicts(files=files, force=force)

    created: list[Path] = []
    for path, content in files.items():
        if not dry_run:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
        created.append(path)

    return created
