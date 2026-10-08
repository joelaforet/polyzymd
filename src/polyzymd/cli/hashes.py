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
    "--verify", is_flag=True, help="Also rehash files that have a recorded hash and compare."
)
@click.option(
    "--dry-run", "dry_run", is_flag=True, help="Report what would be hashed; write nothing."
)
@click.option(
    "--rehash-changed",
    "rehash_changed",
    is_flag=True,
    help="Hash again, and replace, entries this command recorded for files whose size has "
    "changed since, such as a GROMACS run that was extended.",
)
def hash_trajectories_command(
    configs: tuple[Path, ...],
    study_path: Path | None,
    replicate_spec: str | None,
    verify: bool,
    dry_run: bool,
    rehash_changed: bool,
) -> None:
    """Record the SHA-256 of each run's finished trajectory files in trajectory_hashes.json.

    For runs whose runner did not record trajectory hashes (older OpenMM
    runs, downsampled copies, GROMACS runs), with any simulation engine: each
    config's engine (OpenMM, GROMACS) says which files are its finished
    trajectories. It writes only trajectory_hashes.json in each run's engine
    working directory, never progress.json. Running it again changes
    nothing: files with a recorded hash are skipped without being read, and
    a recorded hash is never overwritten; a disagreement is printed as a
    conflict and exits 2. Reading the trajectories takes about a second per
    gigabyte, so on a cluster run it in a batch job.

    \b
    Examples:
        polyzymd hash-trajectories -c LipA_SBMA/config.yaml
        polyzymd hash-trajectories --study lipase_363K --dry-run
    """
    from polyzymd.analyses.study import with_data_dir
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.engines import create_engine
    from polyzymd.utils.replicates import parse_replicate_range

    targets: list[tuple[str, SimulationConfig]] = []
    if study_path is not None:
        from polyzymd.analyses.exceptions import ProtocolError
        from polyzymd.analyses.study_file import load_study_file
        from polyzymd.cli.study import EXIT_STUDY_ERROR

        try:
            protocol = load_study_file(study_path)
        except ProtocolError as exc:
            click.echo(f"error: {' '.join(str(exc).split())}", err=True)
            if exc.hint:
                click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
            sys.exit(EXIT_STUDY_ERROR)
        for label, path in protocol.conditions.items():
            targets.append(
                (label, with_data_dir(SimulationConfig.from_yaml(path), protocol.data.get(label)))
            )
    for path in configs:
        targets.append((path.resolve().parent.name, SimulationConfig.from_yaml(path)))
    if not targets:
        raise click.UsageError("give -c config.yaml or --study.")
    wanted = set(parse_replicate_range(replicate_spec)) if replicate_spec else None
    conflicts = 0
    for label, config in targets:
        # Hashing reads files only, so no engine binary is needed.
        engine = create_engine(config, defer_binary=True)
        runs = [
            (index, root)
            for index, root in config.discover_replicate_dirs()
            if wanted is None or int(index) in wanted
        ]
        if not runs:
            scratch = config.output.effective_scratch_directory
            click.echo(f"{label}: no replicates found under {scratch}")
        for index, root in runs:
            report = engine.record_trajectory_hashes(
                engine.resolve_engine_working_directory(root),
                verify=verify,
                dry_run=dry_run,
                rehash_changed=rehash_changed,
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
            if report["rehashed"]:
                parts.append(
                    f"{'would rehash' if dry_run else 'rehashed'} {len(report['rehashed'])}"
                )
            click.echo(f"{name} ({engine.name}): " + ", ".join(parts))
            for message in report["conflicts"]:
                click.echo(f"conflict: {name} {message}")
                conflicts += 1
    if conflicts:
        sys.exit(EXIT_CONFLICT)
