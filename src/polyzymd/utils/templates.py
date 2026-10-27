"""Shared Jinja helpers for package-resource templates."""

from __future__ import annotations

import json
import shlex
from collections.abc import Mapping
from typing import Any

from jinja2 import Environment, PackageLoader, StrictUndefined


def yaml_quote(value: object) -> str:
    """Return a YAML-safe double-quoted string scalar.

    Parameters
    ----------
    value : object
        String value to quote for YAML output.

    Returns
    -------
    str
        JSON-compatible quoted scalar, which is also valid YAML.
    """
    return json.dumps(str(value), ensure_ascii=False)


def shell_quote(value: object) -> str:
    """Return a POSIX shell-safe single argument.

    Parameters
    ----------
    value : object
        Value to quote for shell command interpolation.

    Returns
    -------
    str
        Shell-escaped value that is parsed as one argument.
    """
    return shlex.quote(str(value))


def create_package_environment(package_name: str, template_dir: str = "templates") -> Environment:
    """Create the shared package-resource Jinja environment.

    Parameters
    ----------
    package_name : str
        Importable package containing the template directory.
    template_dir : str, optional
        Directory within ``package_name`` that contains templates, by default
        ``"templates"``.

    Returns
    -------
    Environment
        Configured Jinja environment for package resources.
    """
    env = Environment(
        loader=PackageLoader(package_name, template_dir),
        undefined=StrictUndefined,
        autoescape=False,
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    env.filters["yaml_quote"] = yaml_quote
    env.filters["shell_quote"] = shell_quote
    return env


def render_package_template(
    package_name: str,
    template_name: str,
    context: Mapping[str, Any] | None = None,
    *,
    template_dir: str = "templates",
) -> str:
    """Render a package-resource Jinja template.

    Parameters
    ----------
    package_name : str
        Importable package containing the template directory.
    template_name : str
        Template filename within ``template_dir``.
    context : Mapping[str, Any] or None, optional
        Values exposed to the template, by default ``None``.
    template_dir : str, optional
        Directory within ``package_name`` that contains templates, by default
        ``"templates"``.

    Returns
    -------
    str
        Rendered template content.
    """
    env = create_package_environment(package_name, template_dir)
    template = env.get_template(template_name)
    return template.render(**dict(context or {}))
