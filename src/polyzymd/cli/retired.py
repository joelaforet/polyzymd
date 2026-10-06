"""Hidden ``polyzymd compare``, ``new-analysis`` and ``init`` commands that name their replacements.

Each command accepts any arguments, print where the workflow went to stderr
and exit with status 2, the status ``polyzymd analyze`` uses for a refused
analysis.
"""

from __future__ import annotations

import sys

import click

EXIT_RETIRED = 2

_ANY_ARGUMENTS = {"ignore_unknown_options": True, "allow_extra_args": True, "help_option_names": []}


def compare_message() -> str:
    """Return the lines ``polyzymd compare`` prints: the replacement commands, docs and skill."""
    from polyzymd.analyses.protocols import ANALYZE_AGENT_SKILL, ANALYZE_PROTOCOL_URL

    return "\n".join(
        [
            "error: polyzymd compare is retired, and comparison.yaml is no longer read.",
            "run an analysis: polyzymd analyze NAME -c config.yaml [-c other/config.yaml ...] "
            "--eq 10ns (the first -c is the control)",
            "run it on SLURM: polyzymd analyze NAME -c config.yaml ... --submit --preset <cluster>",
            f"docs: {ANALYZE_PROTOCOL_URL}",
            f"for an agent: point it at {ANALYZE_AGENT_SKILL} or that page",
        ]
    )


def new_analysis_message() -> str:
    """Return the lines ``polyzymd new-analysis`` prints: the study API, docs and skill."""
    from polyzymd.analyses.protocols import ANALYSIS_API_URL, ANALYZE_AGENT_SKILL

    return "\n".join(
        [
            "error: polyzymd new-analysis is retired: an analysis is now a Python function "
            "of an MDAnalysis Universe, and there is no plugin to scaffold.",
            "write the function and run it with Study.per_replicate (one value per replicate) "
            "or Study.timeseries (one value per frame)",
            f"docs: {ANALYSIS_API_URL} (docs/source/how_to/study_api.md)",
            f"for an agent: point it at {ANALYZE_AGENT_SKILL} or that page",
        ]
    )


@click.command("compare", hidden=True, context_settings=_ANY_ARGUMENTS)
@click.argument("arguments", nargs=-1, type=click.UNPROCESSED)
def compare(arguments: tuple[str, ...]) -> None:
    """Retired: print the polyzymd analyze commands that replace compare, and exit 2."""
    del arguments
    click.echo(compare_message(), err=True)
    sys.exit(EXIT_RETIRED)


@click.command("new-analysis", hidden=True, context_settings=_ANY_ARGUMENTS)
@click.argument("arguments", nargs=-1, type=click.UNPROCESSED)
def new_analysis(arguments: tuple[str, ...]) -> None:
    """Retired: point to the study API that replaces analysis plugins, and exit 2."""
    del arguments
    click.echo(new_analysis_message(), err=True)
    sys.exit(EXIT_RETIRED)


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
