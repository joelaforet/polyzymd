"""``polyzymd study``: commands on a study folder."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import click

EXIT_STUDY_ERROR = 2
#: Context ``meta`` key set by project check, which prints the git and metadata lines once.
PROJECT_CHECK = "polyzymd.project_check"


def _report_status(path: Path) -> str:
    """Return how ``study check`` describes a run's ``report.json``: absent, complete or partial."""
    import json

    if not path.is_file():
        return ""
    try:
        report = json.loads(path.read_text())
    except (OSError, ValueError):
        return " with an unreadable report"
    problems = report.get("problems") or []
    if report.get("status", "complete") == "complete":
        return " with its report"
    return f" with a partial report ({len(problems)} problems): " + "; ".join(problems)


@click.group("study")
def study_group() -> None:
    """Work with a study folder: one MD study's conditions, protocol and results.

    See https://polyzymd.readthedocs.io/en/latest/explanation/study_folders.html.
    """


@study_group.command("check")
@click.argument("path", type=click.Path(path_type=Path), default=Path("."))
@click.option(
    "--production",
    is_flag=True,
    help="Also give each condition's production length, read from every run's trajectory "
    "headers and segments: seconds for a few runs, minutes for long restarted chains.",
)
def check_command(path: Path, production: bool = False) -> None:
    """Check a study.yaml without loading any trajectory.

    PATH is the study.yaml or the folder holding it (default: here). Prints
    one line per condition, saying where its runs were found, and one per
    analysis run, saying whether it has stored results. With --production,
    each condition line also gives how long its replicates were simulated.
    Exits 2 when the file or a condition's config cannot be read; missing
    runs are reported but are not errors, so a study folder without its
    trajectories still checks.
    """
    _study_logging(path, "study-check")
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
    if protocol.description:
        click.echo(f"system: {protocol.description}")
    if protocol.project is not None:
        click.echo(f"project {protocol.project.path} as study {protocol.project_label}")
    for name, where in protocol.structures.items():
        click.echo(f"structure {name}: {where}")
    for name, selection in protocol.regions.items():
        click.echo(f"region {name}: {selection}")
    failed = False
    from polyzymd.analyses.study_file import outside_configs

    outside = outside_configs(protocol)
    labels = list(protocol.conditions)
    controls = labels[:1]
    if protocol.comparison is not None:
        from polyzymd.analyses.study_statistics import comparison_pairs

        within, control = protocol.comparison["within"], protocol.comparison["control"]
        pairs = comparison_pairs(labels, labels, protocol.factors, within, control)
        controls = [label for label in labels if label not in {pair[1] for pair in pairs}]
    for label, config_path in protocol.conditions.items():
        role = "control" if label in controls else "condition"
        if label in outside:
            click.echo(
                f"warning: {role} {label}: the config is outside the study folder, so "
                "polyzymd study freeze will refuse it; copy it in with polyzymd study add-condition"
            )
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
                f"{role} {label}: no replicates found under {where}; stored results can still be read, "
                "and polyzymd study locate DIR finds downloaded runs"
            )
            continue
        missing = sorted(set(protocol.replicates or []) - set(found))
        factors = protocol.factors.get(label)
        click.echo(
            f"{role} {label}: replicates {found} under {where}"
            + (f"; factors {', '.join(f'{k}={v}' for k, v in factors.items())}" if factors else "")
            + (f"; missing replicates {missing}" if missing else "")
            + (_production_summary(label, config_path, protocol) if production else "")
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
            from polyzymd.analyses.study_file import missing_files

            for value in missing_files(user.settings):
                click.echo(
                    f"warning: analysis {run}: setting {value} names no file; give a path "
                    "relative to the file that lists the analysis, or a structure <name>"
                )
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
        equilibration, until = protocol.window(run)
        own = any(x is not None for x in (entry.equilibration, entry.until, entry.stride))
        stride = protocol.stride_of(run)
        window = (
            f"window eq {equilibration}"
            + (f" until {until}" if until else "")
            + (f" stride {stride}" if stride != 1 else "")
            + (" (its own)" if own else "")
        )
        click.echo(
            f"analysis {what}: {settings}; {window}; "
            + (
                f"stored results in {folder}" + _report_status(folder / REPORT_FILE)
                if stored
                else "no stored results"
            )
        )
    from polyzymd.analyses.study_git import describe, git_state

    # A study of a project is committed with it, so project.yaml and the shared
    # analyses/ count among its inputs. project check prints the git line and the
    # project's metadata line once, for all its studies.
    project = protocol.project
    in_project = click.get_current_context().meta.get(PROJECT_CHECK, False)
    if not in_project:
        click.echo(describe(git_state(project.root if project is not None else protocol.root)))
    if not (in_project and protocol.metadata == project.metadata):
        # A study's own metadata, differing from its project's, is its study.yaml's.
        own = project is None or protocol.metadata != project.metadata
        owner = "study.yaml" if own else "project.yaml"
        if not _echo_metadata(protocol.metadata, owner):
            failed = True
    guide = protocol.root / "deposit" / "UPLOAD.md"
    frozen = protocol.root / "manifest.json"
    if guide.is_file():
        click.echo(f"publish: follow {guide}")
    elif frozen.is_file():
        click.echo(
            "reproduce: this study was frozen (manifest.json); point it at downloaded "
            "trajectories with polyzymd study locate DIR --verify, then rerun polyzymd analyze "
            "--study, or redraw figures from results/ without trajectories"
        )
    elif protocol.analyses:
        command = "project freeze" if protocol.project is not None else "study freeze"
        click.echo(f"publish: when the analyses are final, run polyzymd {command}")
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
@click.option(
    "--verify",
    is_flag=True,
    help="Also check every located file's SHA-256 against manifest.json, not only its size.",
)
def locate_command(directory: Path, study_path: Path, verify: bool) -> None:
    """Find each condition's runs under DIRECTORY and record where they are in data.local.yaml.

    Use it after downloading or moving trajectories. For every condition,
    the directory under DIRECTORY holding the most of its run directories
    (named by its config's naming_template), preferring one whose files have the sizes
    manifest.json records, is written to data.local.yaml
    beside study.yaml, which is never committed or published; entries for
    conditions not found are kept as they were. Runs named alike go to the
    folder named for the condition (no_polymer/), and one folder is never
    written for two conditions. Moving data never changes the study or its
    stored results' config hashes.
    """
    _study_logging(study_path, "study-locate")
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
    from polyzymd.analyses.study_scaffold import condition_folder

    located: dict[str, Path] = {}
    missing = []
    configs = {}
    for label, config_path in protocol.conditions.items():
        try:
            configs[label] = SimulationConfig.from_yaml(config_path)
        except (OSError, ValueError) as exc:
            click.echo(f"error: {label}: cannot read {config_path}: {' '.join(str(exc).split())}")
            missing.append(label)
    # Run directories are found by name only, so two conditions whose configs
    # name their runs alike cannot be told apart, except by manifest.json.
    named: dict[str, list[str]] = {}
    for label, config in configs.items():
        named.setdefault(config.format_run_directory_name("*"), []).append(label)
    for label, config in configs.items():
        parents = find_run_parents(config, directory)
        if not parents:
            click.echo(
                f"{label}: no run directories named {config.format_run_directory_name('*')} "
                f"under {directory.resolve()}"
            )
            missing.append(label)
            continue
        candidates = list(parents)
        alike = [other for other in named[config.format_run_directory_name("*")] if other != label]
        by_name = [p for p in parents if condition_folder(label) in p.parts]
        if alike and len(by_name) == 1:
            # One folder per condition, named as in the study (no_polymer/,
            # sds/), as a deposit lays them out, tells them apart before sizes do.
            candidates = by_name
        matching = [p for p in candidates if _matches_manifest(protocol.root, label, p, verify)]
        if alike and len(by_name) != 1 and not matching:
            click.echo(
                f"error: {label}: its runs are named like those of {', '.join(alike)} "
                f"({config.format_run_directory_name('*')}), so their folders cannot be told "
                f"apart: {', '.join(str(p) for p in parents)}"
            )
            click.echo(
                f"fix: Write the folder of {label}'s runs into {protocol.root / DATA_FILE} by "
                f"hand ('{label}: /path/to/folder'), or run study locate on a directory that "
                "holds only its runs."
            )
            missing.append(label)
            continue
        best = max(
            matching or candidates, key=lambda parent: (len(parents[parent]), -len(parent.parts))
        )
        located[label] = best
        others = f" ({len(parents) - 1} other folders also hold some)" if len(parents) > 1 else ""
        click.echo(f"{label}: replicates {parents[best]} under {best}{others}")
        for line in _check_against_manifest(protocol.root, label, best, verify):
            click.echo(line)
            if line.startswith("error"):
                missing.append(label)
    # One folder may hold the runs of several conditions when their runs are
    # named apart. Conditions whose runs are named alike and land in one
    # folder would read the same runs, so neither is written.
    for folder, pattern in {
        (where, configs[label].format_run_directory_name("*")) for label, where in located.items()
    }:
        shared = [
            label
            for label, where in located.items()
            if where == folder and configs[label].format_run_directory_name("*") == pattern
        ]
        if len(shared) > 1:
            click.echo(
                f"error: {' and '.join(shared)} were both located in {folder}; neither is written"
            )
            click.echo(
                f"fix: Write each condition's folder into {protocol.root / DATA_FILE} by hand "
                "('<label>: /path/to/folder'), or run study locate on a directory that holds "
                "each condition's runs in a folder named for it."
            )
            missing.extend(shared)
            for label in shared:
                del located[label]
    if located:
        # Entries this run did not locate are kept as written, hand-written ones included.
        target = protocol.root / DATA_FILE
        written = (yaml.safe_load(target.read_text()) if target.is_file() else None) or {}
        written.update({label: str(folder) for label, folder in located.items()})
        target.write_text(
            "# Where this machine keeps each condition's runs. Written by polyzymd study locate;\n"
            "# never commit or publish it.\n" + yaml.safe_dump(written, sort_keys=False)
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
    help="Create conditions/<label>/ with a template config.yaml to fill in. Repeatable.",
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
        polyzymd study init lipase_363K --condition "No polymer=sims/noPoly/config.yaml" \\
            --condition "SBMA 50%=sims/SBMA50/config.yaml" --equilibration 100ns
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


def _matches_manifest(root: Path, label: str, folder: Path, verify: bool = False) -> bool:
    """Return whether every run file manifest.json lists for ``label`` is in ``folder``.

    Files are compared by size, and with ``verify`` also by SHA-256, which
    tells apart conditions whose runs have the same names and sizes.
    """
    import json

    from polyzymd.analyses.shared.file_hashes import file_sha256
    from polyzymd.analyses.study_freeze import MANIFEST

    try:
        replicates = json.loads((root / MANIFEST).read_text())["conditions"][label]["replicates"]
    except (OSError, ValueError, KeyError):
        return False
    files = [item for replicate in replicates.values() for item in replicate.get("files", [])]
    return bool(files) and all(
        (folder / item["path"]).is_file()
        and (folder / item["path"]).stat().st_size == item["size"]
        and (not verify or file_sha256(folder / item["path"], use_cache=False) == item["sha256"])
        for item in files
    )


def _echo_metadata(metadata: dict, owner: str) -> bool:
    """Print the metadata line of ``owner``; return False when the metadata cannot be read."""
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study_metadata import check_metadata

    try:
        _, gaps = check_metadata(metadata)
    except ProtocolError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}")
        return False
    click.echo(
        f"metadata ({owner}): {len(gaps)} gaps for publishing: {'; '.join(gaps)}"
        if gaps
        else f"metadata ({owner}): complete"
    )
    return True


def _production_summary(label: str, config_path: Path, protocol: Any) -> str:
    """Return ``; production X ns`` for a condition, the range when its replicates differ.

    Choosing an equilibration window or a common window (``until``) needs
    how far each condition was simulated. This reads trajectory headers, not
    frames; when the runs cannot be read, nothing is added.
    """
    import warnings

    from polyzymd.analyses.study import Condition

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            condition = Condition(
                label, config_path, "0ns", protocol.replicates, 1, protocol.data.get(label)
            )
            lengths = [r.production_ns for r in condition.replicates]
    except Exception:  # noqa: BLE001 - the summary is a convenience; check reports real errors
        return ""
    low, high = min(lengths), max(lengths)
    return (
        f"; production {low:.4g} ns"
        if high - low < 1e-6
        else f"; production {low:.4g}-{high:.4g} ns"
    )


def _study_logging(path: Path, command: str) -> None:
    """Keep the console to warnings and log everything to ``<study>/logs/``; see analysis_logging."""
    import logging

    from polyzymd.cli.logging_utils import analysis_logging

    if any(isinstance(h, logging.FileHandler) for h in logging.getLogger().handlers):
        return
    root = Path(path)
    root = root if root.is_dir() else root.parent
    context = click.get_current_context(silent=True)
    verbose = bool(context and context.find_root().params.get("verbose"))
    click.echo(f"log: {analysis_logging(root / 'logs', command, verbose=verbose)}", err=True)


def _check_against_manifest(root: Path, label: str, folder: Path, verify: bool) -> list[str]:
    """Compare located run files with the sizes (and with ``verify``, SHA-256) in manifest.json."""
    import json

    from polyzymd.analyses.shared.file_hashes import file_sha256
    from polyzymd.analyses.study_freeze import MANIFEST

    try:
        manifest = json.loads((root / MANIFEST).read_text())
    except (OSError, ValueError):
        return []
    replicates = manifest.get("conditions", {}).get(label, {}).get("replicates", {})
    if not replicates:
        return []
    checked, problems = 0, []
    for index, replicate in replicates.items():
        for item in replicate.get("files", []):
            path = folder / item["path"]
            if not path.is_file():
                problems.append(f"replicate {index}: {item['path']} is missing")
                continue
            if path.stat().st_size != item["size"]:
                problems.append(f"replicate {index}: {item['path']} has another size")
                continue
            if verify:
                if file_sha256(path, use_cache=False) != item["sha256"]:
                    problems.append(f"replicate {index}: {item['path']} has another SHA-256")
                    continue
            checked += 1
    if problems:
        return [f"error: {label}: {problem}" for problem in problems]
    return [
        f"{label}: {checked} files match manifest.json" + (" (SHA-256)" if verify else " (size)")
    ]


@study_group.command("freeze")
@click.argument("path", type=click.Path(path_type=Path), default=Path("."))
@click.option(
    "--tag", default=None, help="Git tag of the frozen study. Default: study-v1, study-v2, ..."
)
def freeze_command(path: Path, tag: str | None) -> None:
    """Freeze the study at PATH for publication.

    Checks the metadata, the git state and whether each analysis's stored
    results still match the study; hashes the trajectories and writes their
    engine inputs and final frames to deposit/; writes manifest.json,
    md_checklist.yaml, system_summary.csv, CITATION.cff and .zenodo.json;
    commits those and results/, tags the commit, and lays out deposit/ for
    upload, with deposit/UPLOAD.md saying how to publish it on Zenodo; PolyzyMD
    uploads and publishes nothing. Freeze refuses while the study's input files
    have uncommitted changes, so the deposit matches a commit; every other gap,
    such as a missing DOI, is a warning. Refreeze after filling a gap.
    """
    _study_logging(path, "study-freeze")
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study_freeze import freeze

    try:
        result = freeze(path, tag=tag)
    except ProtocolError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}", err=True)
        if exc.hint:
            click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
        sys.exit(EXIT_STUDY_ERROR)
    conditions = result.manifest["conditions"]
    replicates = sum(len(c["replicates"]) for c in conditions.values())
    click.echo(
        f"froze {result.root}"
        + (f" as {result.tag} ({result.commit[:12]})" if result.tag else " without a git tag")
    )
    click.echo(
        f"manifest: {len(result.manifest['files'])} study files, {len(conditions)} conditions, "
        f"{replicates} replicates' files hashed for the manifest"
    )
    click.echo(f"deposit: {result.deposit}; files to upload in {result.upload}")
    for warning in result.warnings:
        click.echo(f"warning: {warning}")
    click.echo(
        f"next: follow {result.guide}, which says how to reserve the DOI, upload and publish "
        "on Zenodo; PolyzyMD uploads nothing"
    )
    if result.git_failed:
        click.echo(
            "error: git could not commit and tag the study; fix that and freeze again", err=True
        )
        sys.exit(EXIT_STUDY_ERROR)


@study_group.command("add-condition")
@click.argument("label")
@click.option(
    "--config",
    "config",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Copy this config.yaml, with the input files it names, into conditions/<label>/.",
)
@click.option(
    "--from",
    "source",
    default=None,
    metavar="OTHER_LABEL",
    help="Copy the config of this condition of the study, with the input files it names.",
)
@click.option(
    "--new", is_flag=True, help="Create conditions/<label>/ with a template config.yaml to fill in."
)
@click.option(
    "--study",
    "study_path",
    type=click.Path(path_type=Path),
    default=Path("."),
    show_default=True,
    help="study.yaml, or the folder holding it.",
)
def add_condition_command(
    label: str, config: Path | None, source: str | None, new: bool, study_path: Path
) -> None:
    """Add the condition LABEL to an existing study and list it in study.yaml.

    Give one of --config, --from and --new. The new config.yaml writes its
    runs into the git-ignored runs/ folder of the project (of the study, for
    a study in no project) unless you set its scratch_directory. With
    --config, when the config's scratch_directory holds its runs, that
    folder is written to data.local.yaml.

    \b
    Examples:
        polyzymd study add-condition "No polymer" --new
        polyzymd study add-condition "SBMA 100%" --from "No polymer"
        polyzymd study add-condition "SBMA 100%" --config sims/SBMA100/config.yaml
    """
    import os

    import yaml

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study_file import DATA_FILE, find_study_file
    from polyzymd.analyses.study_scaffold import RUNS_WARNING, add_condition, runs_folder

    try:
        path = add_condition(study_path, label, config=config, new=new, source=source)
    except ProtocolError as exc:
        click.echo(f"error: {' '.join(str(exc).split())}", err=True)
        if exc.hint:
            click.echo(f"fix: {' '.join(exc.hint.split())}", err=True)
        sys.exit(EXIT_STUDY_ERROR)
    click.echo(f"condition {label}: {path}, listed in study.yaml")
    root = find_study_file(study_path).parent
    data = root / DATA_FILE
    where = (yaml.safe_load(data.read_text()) or {}).get(label) if data.is_file() else None
    if config is not None and where is not None:
        click.echo(f"data {label}: {where} (from the config's scratch_directory)")
    runs = os.path.relpath(runs_folder(root, label))
    click.echo(f"warning: {RUNS_WARNING.format(runs=runs)}")
    click.echo(
        "next: "
        + ("fill in the config, check it with polyzymd validate, " if new else "")
        + "commit, and run polyzymd study check"
    )
