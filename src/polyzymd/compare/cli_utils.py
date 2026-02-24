"""Shared utilities for compare CLI subcommands.

Provides decorators and helper functions that reduce boilerplate across
the ``polyzymd compare`` subcommands (rmsf, triad, contacts, exposure,
run, binding-free-energy).
"""

from __future__ import annotations

import functools
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
