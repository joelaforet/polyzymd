"""Report a simulation config that cannot be loaded as ``error:`` and ``fix:`` lines.

Every command that reads a ``config.yaml`` reports a bad one the same way:
:func:`load_config_or_exit` for a command that stops, and
:func:`report_config_error` for one that reports each condition and goes on.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

import click
import yaml
from pydantic import ValidationError

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig

#: What loading a config raises for a file the user must correct.
CONFIG_LOAD_ERRORS = (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError)


def config_load_lines(error: BaseException, path: str | Path) -> list[tuple[str, str]]:
    """Explain why the config at ``path`` cannot be used, as ``error`` and ``fix`` lines.

    Parameters
    ----------
    error : BaseException
        One of :data:`CONFIG_LOAD_ERRORS`, raised while loading or checking the config.
    path : str or Path
        The config file, named in the fix line.

    Returns
    -------
    list of tuple of str
        ``("error", text)`` and ``("fix", text)`` pairs, in printing order. The
        last fix line always says to run ``polyzymd validate`` again.
    """
    from polyzymd.config.loader import ConfigFileError, config_error_lines

    if isinstance(error, ValidationError):
        return config_error_lines(error, path)
    if isinstance(error, ConfigFileError):
        message, fix = error.message, error.fix
    else:
        message, _, fix = str(error).partition("\nfix: ")
    lines = [("error", " ".join(message.split()))]
    if fix:
        lines.append(("fix", " ".join(fix.split())))
    if "polyzymd validate" not in fix:
        lines.append(("fix", f"correct {path} and run polyzymd validate again."))
    return lines


def report_config_error(
    error: BaseException, path: str | Path, *, prefix: str = "", err: bool = True
) -> None:
    """Print the ``error:`` and ``fix:`` lines of a config that cannot be used.

    Parameters
    ----------
    error : BaseException
        The error raised while loading or checking the config.
    path : str or Path
        The config file.
    prefix : str, optional
        Text put before each error, such as ``"condition B: cannot read X: "``.
    err : bool, optional
        Print to stderr (default) or stdout.
    """
    for kind, text in config_load_lines(error, path):
        if kind == "error":
            click.echo(click.style(f"error: {prefix}{text}", fg="red"), err=err)
        else:
            click.echo(f"fix: {text}", err=err)


def load_config_or_exit(path: str | Path, exit_code: int = 1) -> SimulationConfig:
    """Load the config at ``path``, or report why it cannot be loaded and exit.

    Parameters
    ----------
    path : str or Path
        The config file.
    exit_code : int, optional
        The exit code of the command for a config it cannot load.

    Returns
    -------
    SimulationConfig
        The loaded config.
    """
    from polyzymd.config.schema import SimulationConfig

    try:
        return SimulationConfig.from_yaml(path)
    except CONFIG_LOAD_ERRORS as error:
        report_config_error(error, path)
        sys.exit(exit_code)
