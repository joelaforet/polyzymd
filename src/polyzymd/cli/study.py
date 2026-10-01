"""``polyzymd study``: commands on a study folder."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import click

EXIT_STUDY_ERROR = 2


@click.group("study")
def study_group() -> None:
    """Work with a study folder: one MD study's conditions, protocol and results.

    See https://polyzymd.readthedocs.io/en/latest/explanation/study_folders.html.
    """


@study_group.command("check")
@click.argument("path", type=click.Path(path_type=Path), default=Path("."))
def check_command(path: Path) -> None:
    """Check a study.yaml without loading any trajectory.

    PATH is the study.yaml or the folder holding it (default: here). Prints
    one line per condition, saying where its runs were found, and one per
    analysis run, saying whether it has stored results. Exits 2 when the file
    or a condition's config cannot be read; missing runs are reported but are
    not errors, so a study folder without its trajectories still checks.
    """
    import polyzymd
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.results import REPORT_FILE
    from polyzymd.analyses.study import with_data_dir
    from polyzymd.analyses.study_file import load_study_file
    from polyzymd.citation import citation_line
    from polyzymd.config.schema import SimulationConfig

    try:
        protocol = load_study_file(path)
    except ProtocolError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}", err=True)
        if exc.hint:
            click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
        sys.exit(EXIT_STUDY_ERROR)

    click.echo(
        f"study {protocol.path}  equilibration {protocol.equilibration}  stride {protocol.stride}"
        + (f"  replicates {protocol.replicates}" if protocol.replicates else "")
    )
    if protocol.polyzymd and protocol.polyzymd != polyzymd.__version__:
        click.echo(
            f"warning: written for PolyzyMD {protocol.polyzymd}; this is {polyzymd.__version__}"
        )
    failed = False
    for index, (label, config_path) in enumerate(protocol.conditions.items()):
        role = "control" if index == 0 else "condition"
        try:
            config = SimulationConfig.from_yaml(config_path)
        except (OSError, ValueError) as exc:
            click.echo(
                f"error: {role} {label}: cannot read {config_path}: {' '.join(str(exc).split())}"
            )
            failed = True
            continue
        source = "data.local.yaml" if label in protocol.data else "config"
        config = with_data_dir(config, protocol.data.get(label))
        found = sorted(int(i) for i, _ in config.discover_replicate_dirs())
        where = f"{config.output.effective_scratch_directory} (from {source})"
        if not found:
            click.echo(
                f"{role} {label}: no runs found under {where}; stored results can still be read, "
                "and polyzymd study locate DIR finds downloaded runs"
            )
            continue
        missing = sorted(set(protocol.replicates or []) - set(found))
        click.echo(
            f"{role} {label}: runs {found} under {where}"
            + (f"; missing replicates {missing}" if missing else "")
        )
    from polyzymd.analyses.user_functions import load_function

    for run, entry in protocol.analyses.items():
        folder = protocol.results_dir(run)
        stored = any(folder.glob("polyzymd_results/*/*/replicate_*/record.json"))
        if entry.function is not None:
            user = entry.function
            try:
                load_function(user.file, user.qualname)
            except ProtocolError as exc:
                click.echo(f"error: analysis {run}: {' '.join(str(exc).split())}")
                failed = True
                continue
            relative = (
                user.file.relative_to(protocol.root)
                if user.file.is_relative_to(protocol.root)
                else user.file
            )
            what = f"{run} ({relative}:{user.qualname}, {user.kind})"
            settings = (
                ", ".join(
                    [f"{k}={v!r}" for k, v in user.selections.items()]
                    + [f"{k}={v}" for k, v in user.settings.items()]
                )
                or "no arguments"
            )
        else:
            what = entry.analysis if entry.analysis == run else f"{entry.analysis} as {run}"
            settings = ", ".join(f"{k}={v}" for k, v in entry.settings.items()) or "defaults"
        click.echo(
            f"analysis {what}: {settings}; "
            + (
                f"stored results in {folder}"
                + (" with its report" if (folder / REPORT_FILE).is_file() else "")
                if stored
                else "no stored results"
            )
        )
    from polyzymd.analyses.study_git import describe, git_state

    click.echo(describe(git_state(protocol.root)))
    click.echo(f"cite: {citation_line()}")
    if failed:
        sys.exit(EXIT_STUDY_ERROR)


def find_run_parents(config: Any, root: Path, max_depth: int = 6) -> dict[Path, list[int]]:
    """Return each directory under ``root`` that holds run directories of ``config``, with their replicates.

    A run directory is named by the config's ``naming_template`` with a
    replicate number, as :meth:`SimulationConfig.discover_replicate_dirs`
    matches it. Directories deeper than ``max_depth`` below ``root`` are not
    searched.
    """
    import os
    import re

    pattern = config.format_run_directory_name(replicate="*")
    regex = re.compile("^" + re.escape(pattern).replace(r"\*", r"(?P<replicate>\d+)") + "$")
    parents: dict[Path, list[int]] = {}
    root = root.resolve()
    base_depth = len(root.parts)
    for current, folders, _ in os.walk(root):
        here = Path(current)
        if len(here.parts) - base_depth >= max_depth:
            folders[:] = []
            continue
        for folder in list(folders):
            match = regex.match(folder)
            if match:
                parents.setdefault(here, []).append(int(match.group("replicate")))
                folders.remove(folder)  # do not search inside a run directory
    return {parent: sorted(found) for parent, found in parents.items()}


@study_group.command("locate")
@click.argument("directory", type=click.Path(path_type=Path, exists=True, file_okay=False))
@click.option(
    "--study",
    "study_path",
    type=click.Path(path_type=Path),
    default=Path("."),
    show_default=True,
    help="study.yaml, or the folder holding it.",
)
def locate_command(directory: Path, study_path: Path) -> None:
    """Find each condition's runs under DIRECTORY and record where they are in data.local.yaml.

    Use it after downloading or moving trajectories. For every condition,
    the directory under DIRECTORY holding the most of its run directories
    (named by its config's naming_template) is written to data.local.yaml
    beside study.yaml, which is never committed or published; entries for
    conditions not found are kept as they were. Moving data never changes
    the study or its stored results' config hashes.
    """
    import yaml

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study_file import DATA_FILE, load_study_file
    from polyzymd.config.schema import SimulationConfig

    try:
        protocol = load_study_file(study_path)
    except ProtocolError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}", err=True)
        if exc.hint:
            click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
        sys.exit(EXIT_STUDY_ERROR)
    located = {label: str(folder) for label, folder in protocol.data.items()}
    missing = []
    for label, config_path in protocol.conditions.items():
        try:
            config = SimulationConfig.from_yaml(config_path)
        except (OSError, ValueError) as exc:
            click.echo(f"error: {label}: cannot read {config_path}: {' '.join(str(exc).split())}")
            missing.append(label)
            continue
        parents = find_run_parents(config, directory)
        if not parents:
            click.echo(
                f"{label}: no run directories named {config.format_run_directory_name('*')} "
                f"under {directory.resolve()}"
            )
            missing.append(label)
            continue
        best = max(parents, key=lambda parent: (len(parents[parent]), -len(parent.parts)))
        located[label] = str(best)
        others = f" ({len(parents) - 1} other folders also hold some)" if len(parents) > 1 else ""
        click.echo(f"{label}: runs {parents[best]} under {best}{others}")
    target = protocol.root / DATA_FILE
    target.write_text(
        "# Where this machine keeps each condition's runs. Written by polyzymd study locate;\n"
        "# never commit or publish it.\n"
        + yaml.safe_dump(
            {k: located[k] for k in protocol.conditions if k in located}, sort_keys=False
        )
    )
    click.echo(f"wrote {target}")
    if missing:
        sys.exit(EXIT_STUDY_ERROR)


@study_group.command("init")
@click.argument("directory", type=click.Path(path_type=Path))
@click.option(
    "--condition",
    "conditions",
    multiple=True,
    metavar="LABEL=CONFIG",
    help="Copy an existing config.yaml, with the input files it names, into "
    "conditions/<label>/. Repeatable; the first is the control.",
)
@click.option(
    "--new-condition",
    "new_conditions",
    multiple=True,
    metavar="LABEL",
    help="Create conditions/<label>/ with polyzymd init, to fill in. Repeatable.",
)
@click.option("--equilibration", default=None, help="The study's equilibration window, e.g. 100ns.")
@click.option(
    "--holder",
    default=None,
    help="Copyright holder in LICENSE-data and LICENSE-code. Default: git user.name.",
)
@click.option("--no-git", is_flag=True, help="Do not make the folder a git repository.")
def init_command(
    directory: Path,
    conditions: tuple[str, ...],
    new_conditions: tuple[str, ...],
    equilibration: str | None,
    holder: str | None,
    no_git: bool,
) -> None:
    """Create a study folder at DIRECTORY, ready to version, analyse and publish.

    Writes study.yaml, conditions/, structures/, analyses/, figures/,
    results/, environment/, a README, LICENSE-data (CC-BY-4.0) and
    LICENSE-code (MIT), which you can replace, data.example.yaml and a
    .gitignore, then makes the folder a git repository and commits it.

    \b
    Examples:
        polyzymd study init lipase_363K --condition "No polymer=runs/noPoly/config.yaml" \\
            --condition "SBMA 50%=runs/SBMA50/config.yaml" --equilibration 100ns
        polyzymd study init new_study --new-condition "No polymer" --new-condition "SBMA 50%"
    """
    import subprocess

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study_scaffold import create_study

    parsed: dict[str, Path] = {}
    for entry in conditions:
        label, separator, config = entry.partition("=")
        if not separator or not label.strip() or not config.strip():
            click.echo(f"error: cannot read --condition {entry!r}", err=True)
            click.echo('fix: Write it as --condition "LABEL=path/to/config.yaml".', err=True)
            sys.exit(EXIT_STUDY_ERROR)
        parsed[label.strip()] = Path(config.strip())
    if holder is None:
        try:
            holder = (
                subprocess.run(
                    ["git", "config", "user.name"], capture_output=True, text=True, timeout=10
                ).stdout.strip()
                or None
            )
        except (OSError, subprocess.SubprocessError):
            holder = None
    try:
        created = create_study(
            directory,
            conditions=parsed,
            new_conditions=list(new_conditions),
            equilibration=equilibration,
            holder=holder,
            git=not no_git,
        )
    except ProtocolError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}", err=True)
        if exc.hint:
            click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
        sys.exit(EXIT_STUDY_ERROR)
    click.echo(f"created {created.root}")
    for label, config in created.conditions.items():
        copied = created.copied.get(label, [])
        click.echo(
            f"condition {label}: {config.relative_to(created.root)}"
            + (f", with {len(copied)} input files copied to structures/" if copied else "")
        )
        for path in created.left_absolute.get(label, []):
            click.echo(f"warning: {label}: {path} does not exist, so it was left as it was")
    if not no_git:
        click.echo(
            f"git: committed {created.commit[:12]}"
            if created.commit
            else "warning: git could not commit; run git init and commit yourself"
        )
    if equilibration is None:
        click.echo("note: set equilibration in study.yaml before analysing")
    click.echo(f"next: polyzymd study check {created.root}")
