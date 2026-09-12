"""Shared utilities for compare subcommands in ``polyzymd.cli``.

Provides decorators and helper functions that reduce boilerplate across
the ``polyzymd compare`` subcommands (run, run-all, plot-all).
"""

from __future__ import annotations

import functools
from pathlib import Path

import click
import yaml
from pydantic import ValidationError

DEFAULT_COMPARE_FORMATS = ("table", "markdown", "json")


def common_compare_options(func=None, *, formats=DEFAULT_COMPARE_FORMATS):
    """Decorator adding shared Click options for compare subcommands.

    Adds: -f/--file, --eq-time, --format, -o/--output, -q/--quiet, --debug

    Note: --recompute is not included here because subcommands decide whether
    a recompute flag is meaningful for their workflow.

    Parameters
    ----------
    func : callable, optional
        The command function, when the decorator is used bare.
    formats : sequence of str, optional
        Choices offered by ``--format``. A subcommand that can render an extra
        format passes its own tuple, for example ``compare run`` which also
        accepts ``agent``.

    Returns
    -------
    callable
        The decorated command, or a decorator when ``func`` is ``None``.
    """
    if func is None:
        return functools.partial(common_compare_options, formats=formats)

    @click.option(
        "--debug",
        is_flag=True,
        help="Enable DEBUG logging for troubleshooting.",
    )
    @click.option(
        "-q",
        "--quiet",
        is_flag=True,
        help="Suppress INFO messages, show warnings/errors only.",
    )
    @click.option(
        "-o",
        "--output",
        "output_path",
        type=click.Path(path_type=Path),
        default=None,
        help="Save output to file. JSON comparison results are cached under comparison/.",
    )
    @click.option(
        "--format",
        "output_format",
        type=click.Choice(list(formats)),
        default="table",
        help=f"Output format: table (default), {', '.join(formats[1:])}.",
    )
    @click.option(
        "--eq-time",
        default=None,
        help="Override equilibration time (e.g., '10ns', '5000ps').",
    )
    @click.option(
        "-f",
        "--file",
        "config_file",
        type=click.Path(path_type=Path),
        default="comparison.yaml",
        help="Path to comparison.yaml config file.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def load_comparison_config(config_file: Path) -> ComparisonConfig:
    """Load and return a ComparisonConfig, raising ClickException on error.

    Checks file existence with a friendly error message, then loads
    via ComparisonConfig.from_yaml(). Raises click.ClickException on errors.

    Parameters
    ----------
    config_file : Path
        Path to the comparison.yaml config file.

    Returns
    -------
    ComparisonConfig
        The loaded configuration object.
    """
    from polyzymd.config.comparison import ComparisonConfig

    config_file = Path(config_file).resolve()
    if not config_file.exists():
        raise click.ClickException(
            f"Config file not found: {config_file}\n"
            "Run 'polyzymd compare init -n <name>' to create a comparison project."
        )

    click.echo(f"Loading config: {config_file}")
    try:
        config = ComparisonConfig.from_yaml(config_file)
    except FileNotFoundError as e:
        raise click.ClickException(
            f"{e}\nRun 'polyzymd compare init -n <name>' to create a comparison project."
        )
    except (yaml.YAMLError, ValidationError, ValueError) as e:
        raise click.ClickException(f"Error loading config: {e}") from e

    return config


def validate_and_report(config) -> bool:
    """Validate config and print errors if any.

    Calls config.validate_config() and raises a ClickException
    if there are errors.

    Parameters
    ----------
    config : ComparisonConfig
        The configuration to validate.

    Returns
    -------
    bool
        True if validation passed (no errors).
    """
    errors = config.validate_config()
    if errors:
        msg = "Configuration errors:\n" + "\n".join(f"  - {e}" for e in errors)
        raise click.ClickException(msg)
    return True
