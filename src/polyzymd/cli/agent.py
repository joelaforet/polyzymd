"""``polyzymd agent``: the guidance for AI coding agents in a project."""

from __future__ import annotations

import sys
from pathlib import Path

import click

EXIT_AGENT_ERROR = 2


@click.group("agent")
def agent_group() -> None:
    """Manage the guidance for AI coding agents in a project or study.

    See https://polyzymd.readthedocs.io/en/latest/how_to/use_with_ai_agents.html.
    """


@agent_group.command("sync")
@click.argument("path", type=click.Path(path_type=Path), default=Path("."))
@click.option("--force", is_flag=True, help="Also overwrite files that were edited.")
@click.option("--dry-run", is_flag=True, help="Say what would change; write nothing.")
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["text", "agent"]),
    default="text",
    show_default=True,
    help="agent prints only one line per file.",
)
def sync_command(path: Path, force: bool, dry_run: bool, output_format: str) -> None:
    """Rewrite the agent files of a project or study to this PolyzyMD version.

    PATH is a folder inside the project or study (default: here); the
    project holding it is used, or else the study. The files are AGENTS.md,
    CLAUDE.md, .claude/skills/<name>/SKILL.md and
    .agents/skills/<name>/SKILL.md. Each line says whether a file was
    created, updated, unchanged or skipped. A file you edited is skipped
    unless --force is given. Nothing is committed.

    \b
    Examples:
        polyzymd agent sync
        polyzymd agent sync Paper_1 --dry-run
    """
    from polyzymd.agent_files import find_root, sync
    from polyzymd.analyses.exceptions import ProtocolError

    try:
        root = find_root(path)
        actions = sync(root, force=force, dry_run=dry_run)
    except ProtocolError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}", err=True)
        if exc.hint:
            click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
        sys.exit(EXIT_AGENT_ERROR)
    for action in actions:
        click.echo(action.line())
    if output_format == "agent":
        return
    changed = [action.path for action in actions if action.action in ("created", "updated")]
    if dry_run:
        click.echo(f"dry run: nothing written in {root}")
    elif changed:
        click.echo(f"next: review and commit them in {root}: git add {' '.join(changed)}")
