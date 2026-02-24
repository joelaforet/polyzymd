"""Shared utilities for compare CLI subcommands.

Provides decorators and helper functions that reduce boilerplate across
the ``polyzymd compare`` subcommands (rmsf, triad, contacts, exposure,
run, binding-free-energy).
"""

from __future__ import annotations

import functools
import sys
from pathlib import Path

import click


def common_compare_options(func):
    """Decorator adding shared Click options for compare subcommands.

    Adds: -f/--file, --eq-time, --format, -o/--output, -q/--quiet, --debug

    Note: --recompute is NOT included because exposure uses two separate
    recompute flags (--recompute-sasa, --recompute-exposure).
    """

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
        help="Save output to file. Also saves JSON result to results/ directory.",
    )
    @click.option(
        "--format",
        "output_format",
        type=click.Choice(["table", "markdown", "json"]),
        default="table",
        help="Output format: table (default), markdown, or json.",
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
    """Load and return a ComparisonConfig, exiting on error.

    Checks file existence with a friendly error message, then loads
    via ComparisonConfig.from_yaml(). Exits with sys.exit(1) on any error.

    Parameters
    ----------
    config_file : Path
        Path to the comparison.yaml config file.

    Returns
    -------
    ComparisonConfig
        The loaded configuration object.
    """
    from polyzymd.compare.config import ComparisonConfig

    config_file = Path(config_file).resolve()
    if not config_file.exists():
        click.echo(f"Error: Config file not found: {config_file}", err=True)
        click.echo(
            "Run 'polyzymd compare init -n <name>' to create a comparison project.",
            err=True,
        )
        sys.exit(1)

    click.echo(f"Loading config: {config_file}")
    try:
        config = ComparisonConfig.from_yaml(config_file)
    except FileNotFoundError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo(
            "Run 'polyzymd compare init -n <name>' to create a comparison project.",
            err=True,
        )
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error loading config: {e}", err=True)
        sys.exit(1)

    return config


def save_and_display_result(
    result,
    formatter,
    output_format: str,
    config_file: Path,
    config_name: str,
    result_prefix: str,
    output_path: Path | None = None,
) -> None:
    """Save JSON result, format/display output, and optionally save formatted text.

    Performs the three-step result output pattern shared by compare subcommands:
    1. Save structured JSON result to ``results/`` directory
    2. Format and echo the result to stdout
    3. If ``output_path`` is provided, save formatted text to that file

    Parameters
    ----------
    result : object
        Comparison result object with a ``.save(path)`` method.
    formatter : callable
        Formatter function, called as ``formatter(result, format=output_format)``.
    output_format : str
        Output format string (e.g., ``"table"``, ``"markdown"``, ``"json"``).
    config_file : Path
        Path to the comparison.yaml config file (used to locate ``results/`` dir).
    config_name : str
        Comparison project name (used in the JSON filename).
    result_prefix : str
        Prefix for the JSON filename (e.g., ``"rmsf"``, ``"triad"``).
    output_path : Path or None, optional
        If provided, save formatted text output to this file.
    """
    # Save JSON result
    results_dir = config_file.parent / "results"
    results_dir.mkdir(exist_ok=True)
    json_path = results_dir / f"{result_prefix}_comparison_{config_name.replace(' ', '_')}.json"
    result.save(json_path)
    click.echo(f"Saved result: {json_path}")
    click.echo()

    # Format and display output
    formatted = formatter(result, format=output_format)
    click.echo(formatted)

    # Save formatted output if requested
    if output_path:
        output_path = Path(output_path)
        output_path.write_text(formatted)
        click.echo(f"Saved output: {output_path}")


def validate_and_report(config) -> bool:
    """Validate config and print errors if any.

    Calls config.validate_config() and prints any errors to stderr.
    Exits with sys.exit(1) if there are errors.

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
        click.echo("Configuration errors:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)
    return True
