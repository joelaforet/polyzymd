"""``polyzymd study``: commands on a study folder."""

from __future__ import annotations

import sys
from pathlib import Path

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
        found = sorted(int(i) for i, _ in config.discover_replicate_dirs())
        where = config.output.effective_scratch_directory
        if not found:
            click.echo(
                f"{role} {label}: no runs found under {where} (stored results can still be read)"
            )
            continue
        missing = sorted(set(protocol.replicates or []) - set(found))
        click.echo(
            f"{role} {label}: runs {found} under {where}"
            + (f"; missing replicates {missing}" if missing else "")
        )
    for run, entry in protocol.analyses.items():
        folder = protocol.results_dir(run)
        stored = any(folder.glob("polyzymd_results/*/*/replicate_*/record.json"))
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
    click.echo(f"cite: {citation_line()}")
    if failed:
        sys.exit(EXIT_STUDY_ERROR)
