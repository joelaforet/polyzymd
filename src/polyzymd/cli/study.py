"""CLI commands for the study folder.

A study is one directory that holds every condition, every comparison and every
analysis of a piece of work, so it can be zipped with a paper or kept in git.
See :mod:`polyzymd.config.study`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import click
import yaml

from polyzymd.config.study import ANALYSES_DIR, STUDY_FILE

#: Folders ``polyzymd study init`` creates, each with a ``.gitkeep``.
STUDY_FOLDERS = ("conditions", "comparisons", ANALYSES_DIR, "structures", "workflows")

_README = """\
# {name}

This is a PolyzyMD study. Everything needed to reproduce it lives in this
folder; trajectories live wherever each condition's `output.scratch_directory`
points and are archived separately.

```text
study.yaml      name of the study and where its analyses live
conditions/     one folder per simulated condition, each made by `polyzymd init`
comparisons/    one folder per comparison, each made by `polyzymd compare init`
analyses/       analyses written for this study, shared by every comparison
structures/     inputs shared across conditions (enzyme PDB, docked ligand)
workflows/      scripts that make the paper's tables and figures from results
```

Add a condition and a comparison from this folder:

```bash
polyzymd init -n conditions/<condition>
polyzymd compare init -n <comparison> -o comparisons
```

List each condition's `config.yaml` in the comparison's `comparison.yaml` by a
path relative to that file, for example `../../conditions/<condition>/config.yaml`.

Write an analysis with `polyzymd new-analysis <name>`, run from anywhere in the
study. It creates `analyses/<name>.py` and `analyses/test_<name>.py`. Name the
analysis under `plugins:` in a `comparison.yaml` to run it.
"""

_GITIGNORE = """\
# Trajectories, checkpoints and caches are archived separately, not versioned.
*.dcd
*.xtc
*.trr
*.nc
*.chk
*.cpt
__pycache__/
.pytest_cache/
"""


@click.group()
def study() -> None:
    """Create and manage a study folder."""


@study.command("init")
@click.option("-n", "--name", required=True, help="Name of the study directory to create.")
@click.option("--description", default=None, help="One sentence on what the study asks.")
def init(name: str, description: str | None) -> None:
    """Create a study folder.

    \b
    Example:
        polyzymd study init -n lipase_thermal_stability
    """
    root = Path(name)
    if root.exists():
        click.echo(f"Error: {root} already exists.", err=True)
        sys.exit(1)

    root.mkdir(parents=True)
    study_config = {"name": root.name, "description": description, "analyses": ANALYSES_DIR}
    (root / STUDY_FILE).write_text(yaml.safe_dump(study_config, sort_keys=False))
    (root / "README.md").write_text(_README.format(name=root.name))
    (root / ".gitignore").write_text(_GITIGNORE)
    for folder in STUDY_FOLDERS:
        (root / folder).mkdir()
        (root / folder / ".gitkeep").touch()

    click.echo(f"Created study {root}/")
    for entry in (STUDY_FILE, "README.md", *(f"{folder}/" for folder in STUDY_FOLDERS)):
        click.echo(f"  {entry}")
    click.echo()
    click.echo("Next steps:")
    click.echo(f"  cd {root}")
    click.echo("  polyzymd init -n conditions/<condition>")
    click.echo("  polyzymd compare init -n <comparison> -o comparisons")


@study.command("results")
@click.option(
    "-o",
    "--output",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory for the CSV files. Default: results/ in the study root.",
)
@click.option("--analysis", "analyses", multiple=True, help="Only this analysis; repeatable.")
def results(output: Path | None, analyses: tuple[str, ...]) -> None:
    """Write every comparison's numbers in the study as three CSV tables.

    Run anywhere inside a study. Writes conditions.csv (one row per condition
    and observable), comparisons.csv (one row per pairwise test) and
    profiles.csv (one row per residue or bin of a profile observable).
    """
    from polyzymd.analyses.results import load_results
    from polyzymd.config.study import find_study_root

    root = find_study_root(Path.cwd())
    if root is None:
        click.echo("Error: not inside a study (no study.yaml found above this folder).", err=True)
        sys.exit(1)
    loaded = load_results(root, analyses=analyses or None)
    for warning in loaded.warnings:
        click.secho(f"Warning: {warning}", fg="yellow", err=True)
    written = loaded.to_csv(output or root / "results")
    click.echo(
        f"{len(loaded.conditions)} condition rows, {len(loaded.comparisons)} comparison rows, "
        f"{len(loaded.profiles)} profile rows"
    )
    for path in written:
        click.echo(f"  {path}")


@study.command("export")
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Zip file to write. Default: <study>.zip next to the study folder.",
)
def export(output: Path | None) -> None:
    """Package the study for publication, without its trajectories.

    Run anywhere inside a study. Condition folders no comparison.yaml lists are
    left out. The zip holds bundle_manifest.json, which records every file's
    SHA-256 and every trajectory the results were computed from, to archive
    separately.
    """
    from polyzymd.config.study import find_study_root
    from polyzymd.study_bundle import export_study

    root = find_study_root(Path.cwd())
    if root is None:
        click.echo("Error: not inside a study (no study.yaml found above this folder).", err=True)
        sys.exit(1)
    target = output or root.parent / f"{root.name}.zip"
    plan = export_study(root, target)
    click.echo(
        f"Wrote {target}: {len(plan.files)} files, {plan.size_bytes / 2**20:.1f} MiB, "
        f"{len(plan.conditions)} condition(s)"
    )
    click.echo(f"  {len(plan.trajectories)} trajectory file(s) to archive separately")
    for path in plan.unreferenced:
        click.secho(f"  Left out {path.relative_to(root)}: no comparison lists it", fg="yellow")
    for entry in plan.outside:
        click.secho(f"  Not packaged, outside the study: {entry}", fg="yellow")
    for line in plan.partial:
        click.secho(f"  Partial result packaged: {line}", fg="yellow")


@study.command("verify")
@click.argument("path", type=click.Path(file_okay=False, exists=True, path_type=Path), default=".")
def verify(path: Path) -> None:
    """Check an unpacked study against its bundle_manifest.json.

    Every packaged file is checked by SHA-256. Every trajectory the manifest
    lists is checked by size and content fingerprint if it has been downloaded
    to where its condition's config expects it. Exits non-zero when a file is
    missing or changed, or a trajectory changed.
    """
    from polyzymd.study_bundle import MANIFEST_NAME, verify_study

    if not (path / MANIFEST_NAME).is_file():
        click.echo(f"Error: {path} holds no {MANIFEST_NAME}.", err=True)
        sys.exit(1)
    report = verify_study(path)
    click.echo(f"Checked {report.checked_files} files")
    for name in report.missing_files:
        click.secho(f"  Missing: {name}", fg="red")
    for name in report.changed_files:
        click.secho(f"  Changed: {name}", fg="red")
    click.echo(
        f"Trajectories: {len(report.trajectories_ok)} match, "
        f"{len(report.trajectories_absent)} not downloaded, "
        f"{len(report.trajectories_changed)} changed"
    )
    for name in report.trajectories_changed:
        click.secho(f"  Changed: {name}", fg="red")
    if not report.ok:
        sys.exit(1)
