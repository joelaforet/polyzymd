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
    from polyzymd.cli.study import _study_logging
    from polyzymd.cli.study import check_command as study_check

    _study_logging(path, "project-check")

    try:
        project = Project(path)
        for run in project.protocol.analyses:
            click.echo(f"analysis {run}: studies {', '.join(project.runs_in(run))}")
    except ProtocolError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}", err=True)
        if exc.hint:
            click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
        sys.exit(EXIT_PROJECT_ERROR)
    if project.protocol.stats is not None:
        from polyzymd.analyses.statistics_plan import stats_status

        plan = project.protocol.stats
        click.echo(f"stats {plan.qualname}: {stats_status(project, plan)}")
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


def _fail(exc: Exception) -> None:
    """Print the error and its hint to stderr and exit with status 2."""
    click.echo(f"error: {' '.join(str(exc).split())}", err=True)
    if getattr(exc, "hint", None):
        click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
    sys.exit(EXIT_PROJECT_ERROR)


@project_group.command("init")
@click.argument("path", type=click.Path(path_type=Path))
@click.option(
    "--study",
    "studies",
    multiple=True,
    required=True,
    help="Label of a study of the project, one per protein; also its folder name. Repeatable.",
)
@click.option("--holder", default=None, help="Copyright holder for the licence files.")
@click.option("--no-git", is_flag=True, help="Do not make the project a git repository.")
def init_command(path: Path, studies: tuple[str, ...], holder: str | None, no_git: bool) -> None:
    """Create a project folder at PATH: project.yaml and one empty study per protein.

    Writes project.yaml listing the studies, analyses/, stats/, figures/,
    licences and a README, and one study folder per --study with a
    study.yaml to fill in. To bring existing studies in, follow the tutorial
    "Move existing studies into a project".

    \b
    Examples:
        polyzymd project init Paper_1 --study lipa363 --study calb343 --study rml333
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.project_scaffold import create_project

    try:
        created = create_project(path, list(studies), holder=holder, git=not no_git)
    except ProtocolError as exc:
        _fail(exc)
    click.echo(f"created project {created.root} with studies {', '.join(created.studies)}")
    click.echo(
        "next: fill description:, structures:, regions: and conditions: in each study.yaml "
        "(polyzymd study add-condition adds a condition), the analyses and metadata: in "
        "project.yaml, then polyzymd project check"
    )


@project_group.command("freeze")
@click.argument("path", type=click.Path(path_type=Path), default=Path("."))
@click.option(
    "--tag", default=None, help="Git tag of the frozen project. Default: project-v1, -v2, ..."
)
def freeze_command(path: Path, tag: str | None) -> None:
    """Freeze the project at PATH and every study in it, for one publication.

    Freezes each study (its manifest, checklist, system summary, engine
    inputs and final frames), then writes the project's manifest.json,
    CITATION.cff and .zenodo.json from project.yaml's metadata, commits and
    tags the project, and lays out deposit/ for one upload, with
    deposit/UPLOAD.md. Every gap is a warning; PolyzyMD uploads nothing.
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.project_freeze import freeze_project
    from polyzymd.cli.study import _study_logging

    _study_logging(path, "project-freeze")

    try:
        result = freeze_project(path, tag=tag)
    except ProtocolError as exc:
        _fail(exc)
    conditions = result.manifest["conditions"]
    replicates = sum(len(c["replicates"]) for c in conditions.values())
    click.echo(
        f"froze {result.root}"
        + (f" as {result.tag} ({result.commit[:12]})" if result.tag else " without a git tag")
    )
    click.echo(
        f"manifest: {len(result.manifest['studies'])} studies, {len(conditions)} conditions, "
        f"{replicates} replicates hashed"
    )
    click.echo(f"deposit: {result.deposit}; files to upload in {result.upload}")
    for warning in result.warnings:
        click.echo(f"warning: {warning}")
    click.echo(f"next: follow {result.guide}; PolyzyMD uploads nothing")
