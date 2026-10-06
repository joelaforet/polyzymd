"""``polyzymd stats``: run the statistical plan a project or study names."""

from __future__ import annotations

import sys
from pathlib import Path

import click

EXIT_STATS_ERROR = 2


@click.command("stats")
@click.argument("path", type=click.Path(path_type=Path), default=Path("."))
def stats_command(path: Path) -> None:
    """Run the stats: plan of a project.yaml or study.yaml on its stored results.

    PATH is a project or study folder, or its file (default: here). The plan's
    function receives the project (pz.Project) or study (pz.Study) and returns
    a dict of tables and values; they are written to results/stats/<function>/
    with a record of the plan's code and of the reports it read. No trajectory
    is loaded.

    \b
    Examples:
        polyzymd stats Paper_1
        polyzymd stats Paper_1/lipa363
    """
    from polyzymd.analyses.exceptions import AnalysisError
    from polyzymd.analyses.project_file import PROJECT_FILE
    from polyzymd.analyses.statistics_plan import run_stats_plan
    from polyzymd.cli.study import _study_logging

    path = Path(path).expanduser().resolve()
    _study_logging(path, "stats")
    is_project = path.name == PROJECT_FILE or (path / PROJECT_FILE).is_file()
    try:
        if is_project:
            from polyzymd.analyses.project import Project

            target = Project(path)
        else:
            from polyzymd.analyses.study import Study
            from polyzymd.analyses.study_file import find_study_file

            target = Study(find_study_file(path))
        plan = target.protocol.stats
        if plan is None:
            click.echo(f"error: {target.protocol.path} has no stats: plan.", err=True)
            click.echo("fix: Add 'stats: {plan: stats/plan.py:plan}' to it.", err=True)
            sys.exit(EXIT_STATS_ERROR)
        folder = run_stats_plan(target, plan)
    except AnalysisError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}", err=True)
        if getattr(exc, "hint", None):
            click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
        sys.exit(EXIT_STATS_ERROR)
    written = sorted(p.name for p in folder.iterdir())
    click.echo(f"stats {plan.qualname}: wrote {', '.join(written)} to {folder}")
