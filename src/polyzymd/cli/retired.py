"""The hidden ``polyzymd init`` command, which names its replacement.

It accepts any arguments, prints where the workflow went to stderr and exits
with status 2, the status ``polyzymd analyze`` uses for a refused analysis.
"""

from __future__ import annotations

import sys

import click

EXIT_RETIRED = 2

_ANY_ARGUMENTS = {"ignore_unknown_options": True, "allow_extra_args": True, "help_option_names": []}


@click.command("init", hidden=True, context_settings=_ANY_ARGUMENTS)
@click.argument("arguments", nargs=-1, type=click.UNPROCESSED)
def init(arguments: tuple[str, ...]) -> None:
    """Retired: point to project init and study add-condition --new, and exit 2."""
    del arguments
    click.echo(
        "error: polyzymd init is retired: make a project with polyzymd project init PATH "
        "--study LABEL, then a condition with polyzymd study add-condition LABEL --new",
        err=True,
    )
    sys.exit(EXIT_RETIRED)
