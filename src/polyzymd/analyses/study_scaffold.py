"""Create a study folder: the layout ``polyzymd study init`` writes.

See the "Study folders" explanation page for what each part is for. Each
condition lives in ``conditions/<name>/``: either a new ``polyzymd init``
project, or a copy of an existing ``config.yaml`` with the input files it
names copied into ``conditions/<name>/structures/`` and its paths made
relative, so the folder can be moved and published as a whole. The
config's scratch and projects directories stay as they were: they say
where the data lives, which ``data.local.yaml`` can override.
"""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError

#: Folders of a new study, with what each holds.
FOLDERS = {
    "conditions": "One folder per condition: its config.yaml and structures/",
    "structures": "Reference structures for analyses, such as a crystal structure",
    "analyses": "Your measurement functions, listed in study.yaml as file.py:function",
    "figures": "Notebooks and scripts that turn stored results into the paper's figures",
    "results": "Stored per-replicate values, reports and figures, written by polyzymd analyze",
    "environment": "The software environment that produced the results",
}

GITIGNORE = """\
# Where this machine keeps the trajectories: never commit or publish it.
data.local.yaml
# Python and notebook caches.
__pycache__/
*.pyc
.ipynb_checkpoints/
# What polyzymd study freeze lays out for upload.
deposit/
# Full logs of polyzymd commands; the console shows only warnings.
logs/
# Job scripts and logs of polyzymd analyze --submit, which name this
# machine's paths, and SLURM logs of the simulations.
results/*/slurm/
conditions/*/slurm_logs/
# Software environments, which are rebuilt, never committed.
.pixi/
.venv/
"""

DATA_EXAMPLE = """\
# data.local.yaml says where THIS machine keeps each condition's run
# directories. It is gitignored and never published, because a path is only
# a pointer to where the data lives now: moving data must not change the study.
#
# Without data.local.yaml, each condition's config.yaml says where its runs are
# (its scratch_directory), which suits whoever ran the simulations. After
# downloading trajectories, for example from Zenodo, write it with
#
#     polyzymd study locate DOWNLOAD_DIR
#
# or by hand, as condition label -> directory holding its run directories:
#
# No polymer: /path/to/downloaded/no_polymer
# SBMA 50%: /pl/active/my_lab/polyzymd_sims/LipA_363K
"""

MIT = """\
MIT License

Copyright (c) {year} {holder}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

CC_BY = """\
Copyright (c) {year} {holder}

The data, results and figures of this study (everything except the code in
analyses/ and figures/, which LICENSE-code covers) are licensed under the
Creative Commons Attribution 4.0 International License (CC BY 4.0):
https://creativecommons.org/licenses/by/4.0/legalcode

SPDX-License-Identifier: CC-BY-4.0

Replace this file and LICENSE-code to publish under other licences, and set
metadata.license in study.yaml to match.
"""


@dataclass
class CreatedStudy:
    """What :func:`create_study` wrote."""

    root: Path
    conditions: dict[str, Path]
    copied: dict[str, list[str]] = field(default_factory=dict)
    left_absolute: dict[str, list[str]] = field(default_factory=dict)
    commit: str | None = None


def condition_folder(label: str) -> str:
    """Return the folder name of a condition label: lower case, words joined by ``_``."""
    return re.sub(r"[^\w.+-]+", "_", label.strip().lower()).strip("_") or "condition"


def copy_condition(config: Path, folder: Path) -> tuple[list[str], list[str]]:
    """Copy ``config`` to ``folder/config.yaml`` with the input files it names.

    Every path key of the config (``config.loader.PATH_KEYS``) that names an
    existing file or directory is copied into ``folder/structures/`` and
    written as a path relative to the new config. The directories that say
    where one machine keeps its runs and job files are taken out, as in a
    deposited config: ``projects_directory`` becomes ``.``,
    ``scratch_directory`` becomes ``data`` (where the runs are goes into
    ``data.local.yaml`` instead, :func:`record_data_location`), and a polymer
    ``cache_directory`` is left out so the default is used. The config hash
    leaves these out, so stored results still match.

    Returns
    -------
    tuple of list of str
        The copied paths and the paths left absolute because nothing exists
        there.
    """
    import yaml

    from polyzymd.config.loader import PATH_KEYS, _expand_paths
    from polyzymd.config.schema import SimulationConfig

    config = Path(config).expanduser().resolve()
    try:
        raw = yaml.safe_load(config.read_text()) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ProtocolError(f"Cannot read {config}: {exc}", hint="Give a config.yaml.") from exc
    data = _expand_paths(raw, config.parent)
    structures = folder / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    missing: list[str] = []
    seen: dict[Path, str] = {}

    def place(source: Path) -> str:
        if source in seen:
            return seen[source]
        name, stem, n = source.name, source.stem, 1
        while (structures / name).exists():
            n += 1
            name = f"{stem}_{n}{source.suffix}"
        target = structures / name
        if source.is_dir():
            shutil.copytree(source, target)
        else:
            shutil.copy2(source, target)
        seen[source] = f"structures/{name}"
        copied.append(str(source))
        return seen[source]

    def walk(key: str | None, value: Any) -> Any:
        if isinstance(value, dict):
            return {k: walk(k, v) for k, v in value.items()}
        if isinstance(value, list):
            return [walk(key, item) for item in value]
        if key in PATH_KEYS and isinstance(value, str) and value.lower().strip() != "default":
            source = Path(value)
            if key != "cache_directory" and source.exists():
                return place(source)
            if key != "cache_directory":
                missing.append(value)
        return value

    rewritten = walk(None, data)
    polymers = rewritten.get("polymers")
    if isinstance(polymers, dict):
        polymers.pop("cache_directory", None)
    from polyzymd.analyses.study_freeze import without_machine_paths

    target = folder / "config.yaml"
    target.write_text(
        without_machine_paths(
            f"# Copied by polyzymd from {config.name}\n"
            + yaml.safe_dump(rewritten, sort_keys=False)
        )
    )
    try:
        SimulationConfig.from_yaml(target)
    except (OSError, ValueError) as exc:
        raise ProtocolError(
            f"The copy of {config} in {target} does not load: {exc}",
            hint="Check that the original config loads with polyzymd validate.",
        ) from exc
    return copied, missing


def record_data_location(study_root: Path, label: str, config: Path) -> None:
    """Write where the runs of ``config`` are into the study's ``data.local.yaml``, for ``label``.

    The directory is the config's own scratch directory, read before the
    copy takes it out; a relative one is relative to the config's folder.
    Entries for other conditions are kept.
    """
    import yaml

    from polyzymd.analyses.study_file import DATA_FILE
    from polyzymd.config.schema import SimulationConfig

    try:
        where = SimulationConfig.from_yaml(Path(config)).output.effective_scratch_directory
    except (OSError, ValueError):
        return
    where = Path(where).expanduser()
    if not where.is_absolute():
        where = Path(config).resolve().parent / where
    file = Path(study_root) / DATA_FILE
    current = yaml.safe_load(file.read_text()) if file.is_file() else None
    current = dict(current or {})
    current[label] = str(where.resolve())
    file.write_text(
        "# Where this machine keeps each condition's runs; never committed or published.\n"
        + yaml.safe_dump(current, sort_keys=False)
    )


def _study_yaml(conditions: dict[str, str], equilibration: str | None) -> str:
    import yaml

    listed = (
        yaml.safe_dump({"conditions": conditions}, sort_keys=False)
        if conditions
        else ("conditions: {}   # add each condition: label: conditions/<name>/config.yaml\n")
    )
    window = equilibration or "0ns"
    note = "" if equilibration else "   # set the window to discard from every replicate"
    return f"""\
# The analysis protocol of this study: see
# https://polyzymd.readthedocs.io/en/latest/how_to/study_yaml.html
equilibration: {window}{note}
{listed}analyses: {{}}
# analyses:
#   rg: {{}}
#   contacts:
#     method: occlusion
#   lid_opening:
#     function: analyses/lid.py:lid_distance
#     kind: timeseries
#     unit: A
#     selections:
#       lid: "protein and resid 140-150 and name CA"
#       core: "protein and resid 4-120 and name CA"
"""


def _readme(name: str) -> str:
    from polyzymd.citation import citation_line

    rows = "\n".join(f"| `{folder}/` | {what} |" for folder, what in FOLDERS.items())
    return f"""\
# {name}

A molecular dynamics study made with PolyzyMD. `study.yaml` holds the
analysis protocol: the conditions, the equilibration window and every
analysis setting. Each condition's simulation is in `conditions/<name>/config.yaml`.

| Path | Holds |
|---|---|
| `study.yaml` | The analysis protocol |
{rows}
| `data.example.yaml` | How to say where this machine keeps the trajectories (`data.local.yaml`) |
| `LICENSE-data`, `LICENSE-code` | The licences of the data and of the code |

## Reproduce

Install the PolyzyMD version the results were made with (each report and
`manifest.json` record it; see `environment/`), then:

1. **Figures, from the stored results only:** run the notebooks and scripts in
   `figures/`, which read `pz.Study("study.yaml").results(run)`.
2. **Analyses, from the trajectories:** download them, run
   `polyzymd study locate DOWNLOAD_DIR`, check with `polyzymd study check`,
   and run `polyzymd analyze --study study.yaml`.
3. **Simulations:** each `conditions/<name>/config.yaml` builds and runs its
   condition with `polyzymd`; the replicate number is the random seed.
   Results agree within the statistical noise of MD, not bit for bit.

## Publish

Fill in `metadata:` in `study.yaml` and run `polyzymd study freeze`, which
writes `manifest.json`, `CITATION.cff`, `.zenodo.json` and
`md_checklist.yaml`, commits and tags the study, and lays out `deposit/` for
upload.

## How to cite

Cite this study's paper, and PolyzyMD, which produced the analyses:

> {citation_line()}
"""


def _environment_readme() -> str:
    import polyzymd

    return f"""\
# Environment

This study was started with PolyzyMD {polyzymd.__version__}. To reproduce
it, install the version its results were made with, for example from
https://github.com/joelaforet/polyzymd at that tag, with its pixi environment:

    pixi install -e analysis

Every report records the PolyzyMD version that made it, and
`polyzymd study freeze` records PolyzyMD and its main dependencies in
`manifest.json` and warns when stored results were made with another version.
To pin the whole environment, copy the `pixi.toml` and `pixi.lock` you ran
with into this folder and commit them.
"""


def create_study(
    root: str | Path,
    *,
    conditions: dict[str, Path] | None = None,
    new_conditions: list[str] | None = None,
    equilibration: str | None = None,
    holder: str | None = None,
    git: bool = True,
) -> CreatedStudy:
    """Write a new study folder at ``root``.

    ``conditions`` maps labels to existing configs, copied with their input
    files (:func:`copy_condition`); ``new_conditions`` are labels for which a
    ``polyzymd init`` project is created. With ``git``, the folder is made a
    git repository and everything in it is committed; nothing is committed
    afterwards.

    Raises
    ------
    ProtocolError
        If ``root`` already holds a ``study.yaml``, a label is given twice, or
        a config cannot be copied.
    """
    from polyzymd.analyses.study_file import STUDY_FILE
    from polyzymd.analyses.study_git import init_repository

    root = Path(root).expanduser().resolve()
    if (root / STUDY_FILE).exists():
        raise ProtocolError(
            f"{root} already holds a {STUDY_FILE}.",
            hint="Choose a new folder, or edit the existing study.yaml.",
        )
    labels = [*(conditions or {}), *(new_conditions or [])]
    folders = [condition_folder(label) for label in labels]
    if len(set(labels)) != len(labels) or len(set(folders)) != len(folders):
        raise ProtocolError(
            "Each condition needs its own label and folder name.",
            hint=f"The labels {labels} give the folders {folders}; rename the clashing ones.",
        )
    root.mkdir(parents=True, exist_ok=True)
    for folder in FOLDERS:
        (root / folder).mkdir(exist_ok=True)
    created = CreatedStudy(root, {})
    for label, config in (conditions or {}).items():
        folder = root / "conditions" / condition_folder(label)
        created.copied[label], created.left_absolute[label] = copy_condition(Path(config), folder)
        record_data_location(root, label, Path(config))
        created.conditions[label] = folder / "config.yaml"
    for label in new_conditions or []:
        folder = root / "conditions" / condition_folder(label)
        _polyzymd_init(folder)
        created.conditions[label] = folder / "config.yaml"

    year, holder = date.today().year, holder or "the study's authors"
    listed = {label: str(path.relative_to(root)) for label, path in created.conditions.items()}
    (root / STUDY_FILE).write_text(_study_yaml(listed, equilibration))
    (root / "README.md").write_text(_readme(root.name))
    (root / "data.example.yaml").write_text(DATA_EXAMPLE)
    (root / ".gitignore").write_text(GITIGNORE)
    (root / "LICENSE-code").write_text(MIT.format(year=year, holder=holder))
    (root / "LICENSE-data").write_text(CC_BY.format(year=year, holder=holder))
    (root / "environment" / "README.md").write_text(_environment_readme())
    for folder, what in FOLDERS.items():
        readme = root / folder / "README.md"
        if folder != "environment" and not readme.exists():
            readme.write_text(f"# {folder}/\n\n{what}.\n")
    if git:
        created.commit = init_repository(root, "Create study with polyzymd study init")
    return created


def _polyzymd_init(folder: Path) -> None:
    """Create a ``polyzymd init`` project in ``folder``."""
    from polyzymd.cli.main import init

    folder.parent.mkdir(parents=True, exist_ok=True)
    try:
        init.callback(name=str(folder))
    except SystemExit as exc:
        raise ProtocolError(
            f"polyzymd init could not create {folder}.", hint="Check that the folder is new."
        ) from exc


def add_condition(
    study: str | Path, label: str, *, config: Path | None = None, new: bool = False
) -> Path:
    """Add a condition to an existing study folder and list it in ``study.yaml``.

    With ``config``, the config and its input files are copied as
    :func:`copy_condition` does; with ``new``, ``conditions/<name>/`` is a new
    ``polyzymd init`` project to fill in. The condition is added as one line
    under ``conditions:``, so the rest of ``study.yaml``, comments included, is
    kept as written. Returns the new condition's config path.

    Raises
    ------
    ProtocolError
        If neither or both of ``config`` and ``new`` are given, the label or
        its folder is taken, or the config cannot be copied.
    """
    from polyzymd.analyses.study_file import find_study_file

    if (config is None) == (not new):
        raise ProtocolError(
            "Give either a config to copy or new, not both or neither.",
            hint="polyzymd study add-condition LABEL --config path/to/config.yaml, or --new.",
        )
    file = find_study_file(study)
    # Only the existing labels are read, so a new study with no condition yet,
    # or one still being filled in, takes its first conditions.
    listed = (_read_yaml(file).get("conditions") or {}).keys()
    root = file.parent
    folder = root / "conditions" / condition_folder(label)
    if label in listed or folder.exists():
        raise ProtocolError(
            f"{file} already has a condition {label!r} or a folder {folder}.",
            hint="Choose another label.",
        )
    if config is not None:
        copy_condition(Path(config), folder)
        record_data_location(root, label, Path(config))
    else:
        _polyzymd_init(folder)
    _list_condition(file, label, str((folder / "config.yaml").relative_to(root)))
    if label not in (_read_yaml(file).get("conditions") or {}):
        raise ProtocolError(
            f"{file}: the condition {label!r} was not listed under conditions:.",
            hint="Add it under conditions: by hand.",
        )
    return folder / "config.yaml"


def _read_yaml(file: Path) -> dict:
    """Read a YAML mapping, or raise ProtocolError saying the file is not one."""
    import yaml

    try:
        raw = yaml.safe_load(file.read_text()) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ProtocolError(
            f"Cannot read {file}: {exc}", hint="Check that it is valid YAML."
        ) from exc
    if not isinstance(raw, dict):
        raise ProtocolError(f"{file} must be a YAML mapping.", hint="See the documented example.")
    return raw


def _list_condition(file: Path, label: str, path: str) -> None:
    """Insert ``label: path`` as the last entry of ``conditions:`` in ``file``, keeping the rest."""
    import json
    import re

    lines = file.read_text().splitlines()
    entry = (
        f"  {json.dumps(label) if re.search(r'[:#{}\\[\\],&*!|>%@`]', label) else label}: {path}"
    )
    for index, line in enumerate(lines):
        if re.match(r"^conditions:\s*(\{\s*\})?\s*(#.*)?$", line):
            if "{" in line:
                lines[index] = "conditions:"
                lines.insert(index + 1, entry)
                break
            end = index + 1
            while end < len(lines) and (
                lines[end].startswith((" ", "\t")) or not lines[end].strip()
            ):
                end += 1
            while end > index + 1 and not lines[end - 1].strip():
                end -= 1
            lines.insert(end, entry)
            break
    else:
        raise ProtocolError(
            f"{file} has no top-level 'conditions:' block to add to.",
            hint="Add the condition under conditions: by hand.",
        )
    file.write_text("\n".join(lines) + "\n")
