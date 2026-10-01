"""``polyzymd hash-trajectories``: record trajectory hashes for runs that predate them."""

from __future__ import annotations

import sys
from pathlib import Path

import click

EXIT_CONFLICT = 2


@click.command("hash-trajectories")
@click.option(
    "-c",
    "--config",
    "configs",
    multiple=True,
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    help="Simulation config.yaml; every run it finds is hashed. Repeatable.",
)
@click.option(
    "--study",
    "study_path",
    type=click.Path(path_type=Path),
    default=None,
    help="study.yaml or its folder: hash the runs of every condition (and data.local.yaml).",
)
@click.option(
    "--replicates", "replicate_spec", default=None, help="Replicates, e.g. 1-5. Default: all found."
)
@click.option(
    "--verify", is_flag=True, help="Also rehash segments that have a recorded hash and compare."
)
@click.option(
    "--dry-run", "dry_run", is_flag=True, help="Report what would be hashed; write nothing."
)
@click.option(
    "--force",
    is_flag=True,
    help="Also hash runs recorded as running, once their jobs have stopped.",
)
def hash_trajectories_command(
    configs: tuple[Path, ...],
    study_path: Path | None,
    replicate_spec: str | None,
    verify: bool,
    dry_run: bool,
    force: bool,
) -> None:
    """Record each production segment's trajectory SHA-256 in progress.json.

    For runs that finished before PolyzyMD recorded segment hashes. Running
    it again changes nothing: segments with a recorded hash are skipped
    without reading their files, and a recorded hash is never overwritten;
    a disagreement is printed as a conflict and exits 2. Reading the
    trajectories takes about a second per gigabyte, so on a cluster run it
    in a batch job.

    \b
    Examples:
        polyzymd hash-trajectories -c LipA_SBMA/config.yaml
        polyzymd hash-trajectories --study lipase_363K --dry-run
    """
    from polyzymd.analyses.study import with_data_dir
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.simulation.progress import record_trajectory_hashes
    from polyzymd.utils.replicates import parse_replicate_range

    targets: list[tuple[str, SimulationConfig]] = []
    if study_path is not None:
        from polyzymd.analyses.study_file import load_study_file

        protocol = load_study_file(study_path)
        for label, path in protocol.conditions.items():
            targets.append(
                (label, with_data_dir(SimulationConfig.from_yaml(path), protocol.data.get(label)))
            )
    for path in configs:
        targets.append((path.resolve().parent.name, SimulationConfig.from_yaml(path)))
    if not targets:
        click.echo("error: give -c config.yaml or --study.", err=True)
        sys.exit(EXIT_CONFLICT)
    wanted = set(parse_replicate_range(replicate_spec)) if replicate_spec else None
    conflicts = 0
    for label, config in targets:
        for index, working_dir in config.discover_replicate_dirs():
            if wanted is not None and int(index) not in wanted:
                continue
            report = record_trajectory_hashes(
                working_dir, verify=verify, dry_run=dry_run, force=force
            )
            name = f"{label} replicate {index}"
            if report["skipped"]:
                click.echo(f"{name}: skipped, {report['skipped']}")
                continue
            parts = [
                f"{'would hash' if dry_run else 'hashed'} {len(report['hashed'])}",
                f"already recorded {len(report['recorded'])}",
            ]
            if verify:
                parts.append(f"verified {len(report['verified'])}")
            if report["missing"]:
                parts.append(f"no file for segments {report['missing']}")
            click.echo(f"{name}: " + ", ".join(parts))
            for message in report["conflicts"]:
                click.echo(f"conflict: {name} {message}")
                conflicts += 1
    if conflicts:
        sys.exit(EXIT_CONFLICT)
