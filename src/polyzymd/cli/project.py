"""``polyzymd project``: the studies of one paper."""

from __future__ import annotations

import sys
from pathlib import Path

import click

EXIT_PROJECT_ERROR = 2


@click.group("project")
def project_group() -> None:
    """Work with a project folder: the studies of one paper.

    See https://polyzymd.readthedocs.io/en/latest/explanation/projects.html.
    """


@project_group.command("check")
@click.argument("path", type=click.Path(path_type=Path), default=Path("."))
@click.option(
    "--production",
    is_flag=True,
    help="Also give each condition's production length, as study check --production does; "
    "this reads every run's trajectory.",
)
@click.option(
    "--publish",
    is_flag=True,
    help="Also list every publishing gap of the project and its studies, the freeze step and "
    "the citation, as study check --publish does.",
)
@click.pass_context
def check_command(
    ctx: click.Context, path: Path, production: bool = False, publish: bool = False
) -> None:
    """Check a project.yaml and each of its studies.

    PATH is the project.yaml or the folder holding it (default: here). Prints
    which studies run each project analysis, then the polyzymd study check of
    every study. No trajectory is loaded unless --production is given. The
    last line counts the publishing gaps of the project and its studies;
    --publish lists them. Exits 2 when the project file or any study cannot
    be read.
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.project import Project
    from polyzymd.analyses.study_git import describe, git_state
    from polyzymd.citation import citation_line
    from polyzymd.cli.study import (
        PROJECT_CHECK,
        PROJECT_GAPS,
        _study_logging,
        metadata_gaps_of,
        publishing_summary,
    )
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
    failed = []
    ctx.meta[PROJECT_CHECK] = True
    for label, folder in project.protocol.studies.items():
        click.echo(f"== study {label}")
        try:
            ctx.invoke(study_check, path=folder, production=production, publish=publish)
        except SystemExit as exit_:
            if exit_.code:
                failed.append(label)
    click.echo("== project")
    click.echo(describe(git_state(project.root)))
    gaps = metadata_gaps_of(project.protocol.metadata, "project.yaml")
    metadata_read = gaps is not None
    if publish:
        for gap in gaps or []:
            click.echo(f"publishing gap: {gap}")
        click.echo(f"cite: {citation_line()}")
    else:
        count = ctx.meta.get(PROJECT_GAPS, 0) + len(gaps or [])
        click.echo(publishing_summary(count, "project check"))
    if failed:
        click.echo(f"error: studies {', '.join(failed)} did not check", err=True)
    if failed or not metadata_read:
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
    help="Label of a study of the project; also its folder name. Repeatable.",
)
@click.option("--holder", default=None, help="Copyright holder for the licence files.")
@click.option("--no-git", is_flag=True, help="Do not make the project a git repository.")
def init_command(path: Path, studies: tuple[str, ...], holder: str | None, no_git: bool) -> None:
    """Create a project folder at PATH: project.yaml and one empty study per --study.

    Writes project.yaml listing the studies, analyses/, stats/, figures/,
    licences and a README, and one study folder per --study with a
    study.yaml to fill in. Add a study later with polyzymd project
    add-study. To bring existing studies in, follow the how-to "Move
    existing studies into a project".

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


@project_group.command("add-study")
@click.argument("label")
@click.option(
    "--project",
    "project_path",
    type=click.Path(path_type=Path),
    default=Path("."),
    show_default=True,
    help="project.yaml, or the folder holding it.",
)
def add_study_command(label: str, project_path: Path) -> None:
    """Add the study LABEL to a project: a study folder named LABEL, listed in project.yaml.

    The folder gets a study.yaml to fill in, as project init writes. LABEL
    must be a folder name: lower case, digits and _. Nothing is committed.

    \b
    Examples:
        polyzymd project add-study calb343 --project Paper_1
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.project_scaffold import add_study

    try:
        folder = add_study(project_path, label)
    except ProtocolError as exc:
        _fail(exc)
    click.echo(f"study {label}: {folder}, listed in project.yaml")
    click.echo(
        f"next: fill in {folder / 'study.yaml'}, add conditions with "
        f"polyzymd study add-condition LABEL --new --study {folder}, commit, "
        "and run polyzymd project check"
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
        f"{replicates} replicates' files hashed for the manifest"
    )
    click.echo(f"deposit: {result.deposit}; files to upload in {result.upload}")
    for warning in result.warnings:
        click.echo(f"warning: {warning}")
    click.echo(f"next: follow {result.guide}; PolyzyMD uploads nothing")
    if result.git_failed:
        click.echo(
            "error: git could not commit and tag the project; fix that and freeze again", err=True
        )
        sys.exit(EXIT_PROJECT_ERROR)
