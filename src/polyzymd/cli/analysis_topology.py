"""The ``polyzymd analysis-topology`` command.

Writes ``system.prmtop`` for runs built before PolyzyMD wrote it, so their
analyses read a topology with every bond instead of a PDB whose CONECT
records are unreadable above 99,999 atoms.
"""

from __future__ import annotations

import sys
from pathlib import Path

import click


@click.command("analysis-topology")
@click.argument("run_dirs", nargs=-1, required=True, type=click.Path(path_type=Path))
@click.option("--overwrite", is_flag=True, help="Rewrite system.prmtop where it already exists.")
def analysis_topology_command(run_dirs: tuple[Path, ...], overwrite: bool) -> None:
    """Write system.prmtop for existing runs from their PDB and system.xml.

    Each RUN_DIR is a replicate working directory holding solvated_system.pdb
    and system.xml. New builds write the file themselves; this command is for
    runs that predate it. Exits 1 if any directory could not be converted.

    \b
    Examples:
        polyzymd analysis-topology /scratch/campaign/condition_a/run_1
        polyzymd analysis-topology /scratch/campaign/*/run_*
    """
    from polyzymd.simulation.analysis_topology import rebuild_analysis_topology

    failures = 0
    for run_dir in run_dirs:
        try:
            written = rebuild_analysis_topology(Path(run_dir), overwrite=overwrite)
        except (FileNotFoundError, ValueError, OSError) as exc:
            click.echo(f"{run_dir}: {exc}", err=True)
            failures += 1
            continue
        if written is None:
            click.echo(f"{run_dir}: ParmEd could not convert the system; see the log", err=True)
            failures += 1
        else:
            click.echo(f"{run_dir}: wrote {written.name}")
    if failures:
        sys.exit(1)
