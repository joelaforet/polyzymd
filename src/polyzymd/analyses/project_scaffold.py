"""Create a project folder, one study per protein, and move existing studies into it.

:func:`create_project` writes what ``polyzymd project init`` makes: a
``project.yaml``, shared ``analyses/``, ``stats/`` and ``figures/``, and one
study folder per protein, either new or moved from an existing
``study.yaml`` (:func:`migrate_study`). See the "Projects and studies"
explanation page.
"""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError

PROJECT_FOLDERS = {
    "analyses": "Measurement functions every study uses, listed in project.yaml as file.py:function",
    "stats": "The statistical plan, named by stats: in project.yaml and run by polyzymd stats",
    "figures": "Notebooks and scripts that read results with pz.Project('.').results(run)",
}

PROJECT_GITIGNORE = """\
# Where this machine keeps the trajectories: never commit or publish it.
data.local.yaml
# Python and notebook caches.
__pycache__/
*.pyc
.ipynb_checkpoints/
# What polyzymd project freeze and study freeze lay out for upload.
deposit/
# Full logs of polyzymd commands.
logs/
# SLURM logs of polyzymd analyze --submit and of the simulations.
*/results/*/slurm/*/logs/
*/conditions/*/slurm_logs/
"""


@dataclass
class MigratedStudy:
    """What :func:`migrate_study` did with one existing study."""

    label: str
    root: Path
    structures: dict[str, str] = field(default_factory=dict)
    copied_results: bool = False
    left_absolute: list[str] = field(default_factory=list)


@dataclass
class CreatedProject:
    """What :func:`create_project` wrote."""

    root: Path
    studies: dict[str, Path]
    migrated: dict[str, MigratedStudy] = field(default_factory=dict)
    shared: list[str] = field(default_factory=list)
    commit: str | None = None


def _project_yaml(studies: dict[str, str], analyses: dict[str, Any], metadata: dict) -> str:
    import yaml

    import polyzymd

    head = (
        "# The analyses and publishing metadata of this project: one paper or thesis\n"
        "# chapter. Each study below is one protein (or other system) with its own\n"
        "# conditions, structures and regions; the analyses listed here run in every\n"
        "# study. See https://polyzymd.readthedocs.io/en/latest/how_to/project.html\n"
        f"polyzymd: {polyzymd.__version__}\n\n"
        + yaml.safe_dump({"studies": studies}, sort_keys=False)
    )
    if analyses:
        body = "\n" + yaml.safe_dump({"analyses": analyses}, sort_keys=False)
    else:
        body = """
analyses: {}                  # run in every study, with that study's regions and structures
# analyses:
#   native_contacts:
#     selection: protein and not element H
#     reference_file: structure reference      # each study's structures: reference
#   native_contacts_full:                      # the same analysis over its own window
#     analysis: native_contacts
#     selection: protein and not element H
#     reference_file: structure reference
#     equilibration: 0ns
#   lid_opening:                               # only the studies that have a lid
#     function: analyses/lid.py:lid_distance
#     kind: timeseries
#     unit: A
#     studies: [calb343, rml333]
#     selections:
#       lid: region lid
#       core: region core
"""
    stats = """
# stats:                      # your own statistics, run with polyzymd stats
#   plan: stats/plan.py:plan  # receives the project; use project.replicate_table(run)
"""
    if metadata:
        meta = "\n" + yaml.safe_dump({"metadata": metadata}, sort_keys=False)
    else:
        meta = """
metadata:                     # TODO before polyzymd project freeze
  title: TODO
  authors:
    - {name: TODO, orcid: TODO, affiliation: TODO}
  license: {code: MIT, data: CC-BY-4.0}
  keywords: []
"""
    return head + body + stats + meta


def _study_yaml(conditions: dict[str, str], equilibration: str | None) -> str:
    import yaml

    listed = (
        yaml.safe_dump({"conditions": conditions}, sort_keys=False)
        if conditions
        else "conditions: {}                # control first; add with polyzymd study add-condition\n"
    )
    return f"""\
# One protein (or other system) of the project: its conditions, structures
# and regions. Analyses come from ../project.yaml; add this protein's own
# analyses under analyses:.
description: TODO             # e.g. B. subtilis lipase A (1ISP) at 363 K
equilibration: {equilibration or "0ns"}            # the burn-in to discard from every replicate

structures: {{}}                # name: file under structures/, used as 'structure <name>'
# structures:
#   reference: structures/1ISP_clean.pdb

regions: {{}}                   # name: selection, used as 'region <name>'
# regions:
#   core: resid 5-8 15-27 32-37
#   catalytic_triad: resid 76 132 155

{listed}# conditions:
#   No Polymer: conditions/no_polymer
#   SBMA-EGMA 50:50:
#     config: conditions/sbma_50
#     factors: {{sbma_fraction: 0.5}}           # optional: what varies, for trend tests and plots

analyses: {{}}                  # optional: analyses only this protein runs
"""


def _structure_name(path: Path, taken: set[str]) -> str:
    name = re.sub(r"[^\w-]+", "_", path.stem) or "structure"
    while name in taken:
        name += "_2"
    return name


def migrate_study(source: str | Path, root: Path, label: str) -> tuple[MigratedStudy, dict]:
    """Move the study of ``source`` (a ``study.yaml`` or its folder) into a new study at ``root``.

    The conditions' configs and input files are copied into ``conditions/``
    (:func:`~polyzymd.analyses.study_scaffold.copy_condition`), where each
    condition's runs are found today goes into ``data.local.yaml``, the
    study's ``analyses/`` code and ``results/`` are copied, and every setting
    that names an existing file is copied into ``structures/`` and written as
    ``structure <name>`` (``reference`` when the study names only one file).
    Settings naming other absolute paths are kept and listed in
    ``left_absolute``. The source is only read.

    Returns
    -------
    tuple
        What was done, and the study's raw ``metadata:`` (the project's to keep).
    """
    import yaml

    from polyzymd.analyses.study import with_data_dir
    from polyzymd.analyses.study_file import DATA_FILE, find_study_file, load_study_file
    from polyzymd.analyses.study_scaffold import create_study
    from polyzymd.config.schema import SimulationConfig

    file = find_study_file(source)
    old = load_study_file(file)
    raw = yaml.safe_load(file.read_text()) or {}
    created = create_study(
        root, conditions=dict(old.conditions), equilibration=old.equilibration, git=False
    )
    for name in ("LICENSE-code", "LICENSE-data", "data.example.yaml"):
        (root / name).unlink(missing_ok=True)
    data = {}
    for condition, config in old.conditions.items():
        loaded = with_data_dir(SimulationConfig.from_yaml(config), old.data.get(condition))
        data[condition] = str(loaded.output.effective_scratch_directory)
    (root / DATA_FILE).write_text(
        "# Where this machine keeps each condition's runs (written by polyzymd project init).\n"
        + yaml.safe_dump(data, sort_keys=False)
    )
    if (file.parent / "analyses").is_dir():
        shutil.copytree(file.parent / "analyses", root / "analyses", dirs_exist_ok=True)
    migrated = MigratedStudy(label, root)
    if (file.parent / "results").is_dir():
        shutil.copytree(
            file.parent / "results",
            root / "results",
            dirs_exist_ok=True,
            ignore=shutil.ignore_patterns("logs", "__pycache__"),
        )
        migrated.copied_results = True

    files: dict[Path, str] = {}

    def collect(value: Any) -> None:
        if isinstance(value, dict):
            for item in value.values():
                collect(item)
        elif isinstance(value, list):
            for item in value:
                collect(item)
        elif isinstance(value, str) and Path(value).is_absolute() and Path(value).is_file():
            files.setdefault(Path(value), "")

    entries = dict(raw.get("analyses") or {})
    collect(entries)
    taken: set[str] = set()
    for path in files:
        name = "reference" if len(files) == 1 else _structure_name(path, taken)
        taken.add(name)
        target = root / "structures" / path.name
        shutil.copy2(path, target)
        files[path] = name
        migrated.structures[name] = f"structures/{path.name}"

    def rewrite(value: Any) -> Any:
        if isinstance(value, dict):
            return {key: rewrite(item) for key, item in value.items()}
        if isinstance(value, list):
            return [rewrite(item) for item in value]
        if isinstance(value, str) and Path(value).is_absolute():
            if Path(value) in files:
                return f"structure {files[Path(value)]}"
            migrated.left_absolute.append(value)
        if isinstance(value, str) and value.endswith(".py") is False and ".py:" in value:
            file_part, _, function = value.rpartition(":")
            return f"analyses/{Path(file_part).name}:{function}"
        return value

    analyses = {run: rewrite(entry) for run, entry in entries.items()}
    conditions = {
        label_: str(path.relative_to(root)) for label_, path in created.conditions.items()
    }
    study: dict[str, Any] = {"description": raw.get("description") or "TODO"}
    for key in ("equilibration", "stride", "until", "replicates"):
        if key in raw:
            study[key] = raw[key]
    study["structures"] = migrated.structures
    study["regions"] = {}
    study["conditions"] = conditions
    study["analyses"] = analyses
    (root / "study.yaml").write_text(
        f"# Moved into a project by polyzymd project init from {file}\n"
        + yaml.safe_dump(study, sort_keys=False)
    )
    return migrated, dict(raw.get("metadata") or {})


def create_project(
    root: str | Path,
    studies: dict[str, Path | None],
    *,
    holder: str | None = None,
    git: bool = True,
) -> CreatedProject:
    """Write a new project folder at ``root`` with one study per entry of ``studies``.

    ``studies`` maps each study's label to ``None`` for a new, empty study, or
    to an existing ``study.yaml`` (or its folder) to move in
    (:func:`migrate_study`). Analyses that every moved study defines the same
    way are moved into ``project.yaml``. With ``git``, the project is made a
    git repository and everything is committed.

    Raises
    ------
    ProtocolError
        If ``root`` already holds a ``project.yaml``, no study is given, or a
        label cannot be a folder name.
    """
    import yaml

    from polyzymd.analyses.project_file import PROJECT_FILE
    from polyzymd.analyses.study_git import init_repository
    from polyzymd.analyses.study_scaffold import CC_BY, MIT, condition_folder, create_study

    root = Path(root).expanduser().resolve()
    if (root / PROJECT_FILE).exists():
        raise ProtocolError(f"{root} already holds a {PROJECT_FILE}.", hint="Choose a new folder.")
    if not studies:
        raise ProtocolError(
            "A project needs at least one study.",
            hint="Give --study LABEL, or --study LABEL=path/to/old/study.yaml to move one in.",
        )
    for label in studies:
        if condition_folder(label) != label:
            raise ProtocolError(
                f"The study label {label!r} is not a folder name.",
                hint=f"Use lower case, digits and _ only, such as {condition_folder(label)!r}.",
            )
    root.mkdir(parents=True, exist_ok=True)
    for folder, what in PROJECT_FOLDERS.items():
        (root / folder).mkdir(exist_ok=True)
        readme = root / folder / "README.md"
        if not readme.exists():
            readme.write_text(f"# {folder}/\n\n{what}.\n")
    created = CreatedProject(root, {})
    metadata: dict = {}
    for label, source in studies.items():
        folder = root / label
        if source is None:
            create_study(folder, git=False)
            for name in ("LICENSE-code", "LICENSE-data"):
                (folder / name).unlink(missing_ok=True)
            (folder / "study.yaml").write_text(_study_yaml({}, None))
        else:
            migrated, meta = migrate_study(source, folder, label)
            created.migrated[label] = migrated
            metadata = metadata or meta
        created.studies[label] = folder

    shared: dict[str, Any] = {}
    moved = [root / label / "study.yaml" for label in created.migrated]
    if len(moved) > 1:
        texts = [yaml.safe_load(path.read_text()) for path in moved]
        common = set.intersection(*(set(t.get("analyses") or {}) for t in texts))
        for run in sorted(common):
            if all(t["analyses"][run] == texts[0]["analyses"][run] for t in texts):
                shared[run] = texts[0]["analyses"][run]
        if shared:
            for path, text in zip(moved, texts, strict=True):
                text["analyses"] = {k: v for k, v in text["analyses"].items() if k not in shared}
                head = path.read_text().splitlines()[0]
                path.write_text(head + "\n" + yaml.safe_dump(text, sort_keys=False))
            for run, entry in shared.items():
                function = entry.get("function") if isinstance(entry, dict) else None
                if function:
                    name = Path(function.rpartition(":")[0]).name
                    first = root / next(iter(created.migrated)) / "analyses" / name
                    if first.is_file():
                        shutil.copy2(first, root / "analyses" / name)
    created.shared = sorted(shared)

    (root / PROJECT_FILE).write_text(
        _project_yaml({label: label for label in studies}, shared, metadata)
    )
    year, holder = date.today().year, holder or "the project's authors"
    (root / ".gitignore").write_text(PROJECT_GITIGNORE)
    (root / "LICENSE-code").write_text(MIT.format(year=year, holder=holder))
    (root / "LICENSE-data").write_text(CC_BY.format(year=year, holder=holder))
    (root / "README.md").write_text(
        f"# {root.name}\n\nA PolyzyMD project: one paper, one study per protein "
        f"({', '.join(studies)}).\n`project.yaml` lists the studies and the analyses each "
        "runs; each study's `study.yaml` holds that protein's conditions, structures and "
        "regions.\n\n```bash\npolyzymd project check .\npolyzymd analyze --project .\n"
        "polyzymd stats .\npolyzymd project freeze .\n```\n"
    )
    if git:
        created.commit = init_repository(root, "Create project with polyzymd project init")
    return created
