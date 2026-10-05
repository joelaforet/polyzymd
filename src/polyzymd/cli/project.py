"""``polyzymd project``: the studies of one paper, one per protein."""

from __future__ import annotations

import sys
from pathlib import Path

import click

EXIT_PROJECT_ERROR = 2


@click.group("project")
def project_group() -> None:
    """Work with a project folder: one paper's studies, one study per protein.

    See https://polyzymd.readthedocs.io/en/latest/explanation/projects.html.
    """


@project_group.command("check")
@click.argument("path", type=click.Path(path_type=Path), default=Path("."))
@click.pass_context
def check_command(ctx: click.Context, path: Path) -> None:
    """Check a project.yaml and each of its studies without loading any trajectory.

    PATH is the project.yaml or the folder holding it (default: here). Prints
    which studies run each project analysis, then the polyzymd study check of
    every study. Exits 2 when the project file or any study cannot be read.
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.project import Project
    from polyzymd.cli.study import check_command as study_check

    try:
        project = Project(path)
        for run in project.protocol.analyses:
            click.echo(f"analysis {run}: studies {', '.join(project.runs_in(run))}")
    except ProtocolError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}", err=True)
        if exc.hint:
            click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
        sys.exit(EXIT_PROJECT_ERROR)
    failed = []
    for label, folder in project.protocol.studies.items():
        click.echo(f"== study {label}")
        try:
            ctx.invoke(study_check, path=folder)
        except SystemExit as exit_:
            if exit_.code:
                failed.append(label)
    if failed:
        click.echo(f"error: studies {', '.join(failed)} did not check", err=True)
        sys.exit(EXIT_PROJECT_ERROR)
