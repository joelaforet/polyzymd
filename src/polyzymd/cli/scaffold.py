"""Scaffold generator for new analysis plugins.

One analysis is one module holding a settings model and a ``compute()`` that
returns observables, plus one test module. Both come from Jinja templates
stored as package resources.

Usage::

    polyzymd new-analysis my_analysis
    polyzymd new-analysis my_analysis --dry-run
    polyzymd new-analysis my_analysis --force
"""

from __future__ import annotations

import keyword
import re
from pathlib import Path

from polyzymd.cli._scaffold.models import ScaffoldSpec
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
        analysis plugins, by default True.

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


def _build_spec(name: str, class_name: str | None, *, in_study: bool = False) -> ScaffoldSpec:
    """Validate user inputs and build a scaffold render specification.

    Parameters
    ----------
    name : str
        Plugin name in snake_case.
    class_name : str or None
        Optional PascalCase class prefix.
    in_study : bool, optional
        Whether the plugin goes in a study's analyses folder, by default False.

    Returns
    -------
    ScaffoldSpec
        Validated scaffold rendering specification.

    Raises
    ------
    ValueError
        If the name or the class prefix is invalid.
    """
    name_error = validate_name(name, check_existing=False)
    if name_error:
        raise ValueError(name_error)

    cls = class_name or to_pascal_case(name)
    cls_error = validate_class_name(cls)
    if cls_error:
        raise ValueError(cls_error)

    return ScaffoldSpec(name=name, class_name=cls, module=name if in_study else None)


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


def _check_registered_name_conflict(spec: ScaffoldSpec) -> None:
    """Reject registered-name collisions uniformly.

    Parameters
    ----------
    spec : ScaffoldSpec
        Validated scaffold rendering specification.

    Raises
    ------
    ValueError
        If the requested name collides with a built-in or external analysis
        registration.
    """
    try:
        from polyzymd.analyses.discovery import list_all_names

        if spec.name in list_all_names():
            raise ValueError(f"'{spec.name}' already exists as a registered analysis plugin.")
    except ModuleNotFoundError as exc:
        if exc.name not in ("polyzymd.analyses", "polyzymd.analyses.discovery"):
            raise


def _check_target_conflicts(files: dict[Path, str], project_root: Path, *, force: bool) -> None:
    """Preflight all scaffold output paths before writing.

    Parameters
    ----------
    files : dict[Path, str]
        Rendered scaffold files keyed by output path.
    project_root : Path
        Repository root for the scaffold operation.
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
    existing_source_files = [
        path
        for path in existing_files
        if path.resolve() in _source_target_paths(files, project_root)
    ]
    if existing_source_files:
        raise FileExistsError(
            f"Scaffold source target already exists: {_format_paths(existing_source_files)}. "
            "Choose a different analysis name."
        )

    if existing_files and not force:
        raise FileExistsError(
            f"Scaffold target path already exists: {_format_paths(existing_files)}. "
            "Use --force to overwrite."
        )


def _check_layout_conflicts(spec: ScaffoldSpec, project_root: Path) -> None:
    """Reject a new module whose name is already taken by a plugin package.

    Parameters
    ----------
    spec : ScaffoldSpec
        Validated scaffold rendering specification.
    project_root : Path
        Repository root directory containing ``src/`` and ``tests/``.

    Raises
    ------
    FileExistsError
        If a package of that name already exists.
    """
    analyses_root = project_root / "src" / "polyzymd" / "analyses"
    package_path = analyses_root / spec.name
    if package_path.exists():
        raise FileExistsError(
            f"Cannot create {analyses_root / f'{spec.name}.py'}: "
            f"package layout {package_path} already exists."
        )


def generate_scaffold(
    name: str,
    project_root: Path | None = None,
    *,
    class_name: str | None = None,
    force: bool = False,
    dry_run: bool = False,
    study_root: Path | None = None,
) -> list[Path]:
    """Create scaffold files for a new analysis plugin.

    Parameters
    ----------
    name : str
        Plugin name in snake_case, for example ``"solvent_shell"``.
    project_root : Path or None, optional
        Repository root directory containing ``src/`` and ``tests/``, for a
        built-in analysis. Not used when ``study_root`` is given.
    class_name : str or None, optional
        PascalCase class prefix. Auto-derived from *name* when omitted, by
        default None.
    force : bool, optional
        Overwrite existing files, by default False.
    dry_run : bool, optional
        Return paths without writing files, by default False.
    study_root : Path or None, optional
        Root of a study. When given, the plugin and its test are written to the
        study's analyses folder instead of the PolyzyMD source tree.

    Returns
    -------
    list[Path]
        Paths of created, or would-be-created, files.

    Raises
    ------
    FileExistsError
        If a target path exists and ``force`` is False.
    ValueError
        If the plugin name or the class prefix is invalid, or neither root is
        given.
    """
    spec = _build_spec(name=name, class_name=class_name, in_study=study_root is not None)
    if study_root is not None:
        from polyzymd.config.study import STUDY_FILE, StudyConfig

        analyses_dir = (
            Path(study_root) / StudyConfig.from_yaml(Path(study_root) / STUDY_FILE).analyses
        )
        _check_study_name_conflict(spec=spec, analyses_dir=analyses_dir)
        files = render_scaffold(spec=spec, project_root=Path(study_root), analyses_dir=analyses_dir)
        _check_target_conflicts(files=files, project_root=Path(study_root), force=force)
    else:
        if project_root is None:
            raise ValueError("generate_scaffold() needs project_root or study_root")
        _check_layout_conflicts(spec=spec, project_root=project_root)
        files = render_scaffold(spec=spec, project_root=project_root)
        _check_registered_name_conflict(spec=spec)
        _check_target_conflicts(files=files, project_root=project_root, force=force)

    created: list[Path] = []
    for path, content in files.items():
        if not dry_run:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
        created.append(path)

    return created


def _check_study_name_conflict(spec: ScaffoldSpec, analyses_dir: Path) -> None:
    """Reject a name a built-in analysis or another study module already uses.

    Raises
    ------
    ValueError
        If the name is a built-in analysis or an importable module, since the
        generated test imports the plugin by that name.
    FileExistsError
        If a module or package of that name already exists in the study's
        analyses folder.
    """
    import importlib.util

    from polyzymd.analyses.discovery import _cached_registry

    if spec.name in _cached_registry():
        raise ValueError(
            f"'{spec.name}' is a built-in PolyzyMD analysis; choose another name for the "
            "study's analysis."
        )
    if importlib.util.find_spec(spec.name) is not None:
        raise ValueError(
            f"'{spec.name}' is already an importable Python module, so the generated test "
            "could not import the analysis by that name. Choose another name."
        )
    for existing in (analyses_dir / f"{spec.name}.py", analyses_dir / spec.name):
        if existing.exists():
            raise FileExistsError(f"{existing} already exists. Choose a different analysis name.")
