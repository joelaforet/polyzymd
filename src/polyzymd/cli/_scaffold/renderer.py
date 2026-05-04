"""Jinja renderer for analysis scaffold package resources."""

from __future__ import annotations

from pathlib import Path

from jinja2 import Environment

from polyzymd.cli._scaffold.models import ScaffoldSpec
from polyzymd.utils.templates import create_package_environment


def create_environment() -> Environment:
    """Create the Jinja environment used for scaffold templates.

    Returns
    -------
    Environment
        Configured Jinja environment using package-resource templates.
    """
    return create_package_environment("polyzymd.cli._scaffold", "templates")


def render_template(template_name: str, spec: ScaffoldSpec) -> str:
    """Render one scaffold template with a specification.

    Parameters
    ----------
    template_name : str
        Template filename within the package-resource template directory.
    spec : ScaffoldSpec
        Scaffold rendering specification.

    Returns
    -------
    str
        Rendered file content.
    """
    env = create_environment()
    template = env.get_template(template_name)
    return template.render(spec=spec)


def render_scaffold(spec: ScaffoldSpec, project_root: Path) -> dict[Path, str]:
    """Render all files for one analysis scaffold.

    Parameters
    ----------
    spec : ScaffoldSpec
        Scaffold rendering specification.
    project_root : Path
        Repository root that will contain generated ``src`` and ``tests`` files.

    Returns
    -------
    dict[Path, str]
        Mapping of output paths to rendered content.
    """
    analyses_dir = project_root / "src" / "polyzymd" / "analyses" / spec.name
    tests_dir = project_root / "tests" / "analyses" / "plugins"

    files = {
        analyses_dir / "__init__.py": render_template("plugin_init.py.jinja", spec),
        analyses_dir / "_runner.py": render_template("runner.py.jinja", spec),
        tests_dir / f"test_{spec.name}.py": render_template("test_plugin.py.jinja", spec),
    }
    if spec.uses_pydantic_results:
        files[analyses_dir / "_results.py"] = render_template("results.py.jinja", spec)
    return files
