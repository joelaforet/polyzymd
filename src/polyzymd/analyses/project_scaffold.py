"""Create a project folder: the layout ``polyzymd project init`` writes.

A project holds one paper: ``project.yaml``, the shared ``analyses/``,
``stats/`` and ``figures/`` folders, and one study folder for each study
label, each with a ``study.yaml`` to fill in. See the "Projects and studies" explanation
page and the tutorial on moving existing studies into a project.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from pathlib import Path

from polyzymd.analyses.exceptions import ProtocolError

#: Folders of a new project, with what each holds.
PROJECT_FOLDERS = {
    "analyses": "Measurement functions every study uses, listed in project.yaml as file.py:function",
    "stats": "Statistics scripts that read pz.Project('.').replicate_table(run); freeze publishes them",
    "figures": "Notebooks and scripts that read results with pz.Project('.').results(run)",
}

PROJECT_GITIGNORE = """\
# Where this machine keeps the trajectories: never commit or publish it.
data.local.yaml
# The runs of the conditions' configs, unless they set scratch_directory.
runs/
# Python and notebook caches.
__pycache__/
*.pyc
.ipynb_checkpoints/
# What polyzymd project freeze and study freeze lay out for upload.
deposit/
# Full logs of polyzymd commands.
logs/
# Job scripts and logs of polyzymd analyze --submit, which name this
# machine's paths, and SLURM logs of the simulations.
**/results/*/slurm/
**/conditions/*/slurm_logs/
# Software environments, which are rebuilt, never committed.
.pixi/
.venv/
"""

PROJECT_YAML = """\
# The analyses and publishing metadata of this project: one paper or thesis
# chapter. Each study below is a set of conditions compared with each other,
# with its own conditions, structures and regions; the analyses listed here
# run in every study. See https://polyzymd.readthedocs.io/en/latest/how_to/project.html

{studies}
# Run in every study, with that study's regions and structures.
analyses: {{}}
# analyses:
#   native_contacts:
#     selection: protein and not element H
#     reference_file: structure reference      # each study's structures: reference
#   native_contacts_full:                      # the same analysis over its own window
#     analysis: native_contacts
#     selection: protein and not element H
#     reference_file: structure reference
#     equilibration: 0ns
#   lid_opening:                               # only the studies whose protein has a lid
#     function: analyses/lid.py:lid_distance
#     kind: timeseries
#     unit: A
#     studies: [calb343, rml333]
#     selections:
#       lid: region lid
#       core: region core

metadata:                     # TODO before polyzymd project freeze
  title: TODO
  authors:
    - {{family-names: TODO, given-names: TODO, orcid: TODO, affiliation: TODO}}
  license: {{code: MIT, data: CC-BY-4.0}}
  keywords: []
"""

STUDY_YAML = """\
# One study of the project: conditions compared with each other. They share
# one residue numbering, these structures and regions, one equilibration
# window and the control (the first condition). Analyses come from
# ../project.yaml; add this study's own analyses under analyses:.
description: TODO             # e.g. B. subtilis lipase A (1ISP) at 363 K
equilibration: 0ns            # TODO: the burn-in to discard from every replicate

# name: file under structures/, used as 'structure <name>'.
structures: {}
# structures:
#   reference: structures/1ISP_clean.pdb

# name: selection, used as 'region <name>'.
regions: {}
# regions:
#   core: resid 5-8 15-27 32-37
#   catalytic_triad: resid 76 132 155

# Control first; add each with polyzymd study add-condition.
conditions: {}
# conditions:
#   No Polymer: conditions/no_polymer
#   SBMA-EGMA 50:50:
#     config: conditions/sbma_50
#     factors: {sbma_fraction: 0.5}           # optional: what varies, for trend tests and plots

# Optional: analyses only this study runs.
analyses: {}
"""


@dataclass
class CreatedProject:
    """What :func:`create_project` wrote.

    Attributes
    ----------
    root : Path
        The project folder.
    studies : dict of str to Path
        Each study's folder, by label.
    commit : str or None
        The first commit, or ``None`` without git.
    """

    root: Path
    studies: dict[str, Path]
    commit: str | None = None


def create_project(
    root: str | Path,
    studies: list[str],
    *,
    holder: str | None = None,
    git: bool = True,
) -> CreatedProject:
    """Write a new project folder at ``root`` with one empty study per label.

    Writes ``project.yaml`` listing the studies, ``analyses/``, ``stats/``
    and ``figures/`` with a README each, ``.gitignore``, ``LICENSE-code``
    (MIT), ``LICENSE-data`` (CC-BY-4.0) and ``README.md``, and for each label
    a study folder made by
    :func:`~polyzymd.analyses.study_scaffold.create_study` with a
    ``study.yaml`` to fill in (description, equilibration, structures,
    regions, conditions).

    Parameters
    ----------
    root : str or Path
        The new project folder.
    studies : list of str
        Study labels; each is also its folder name.
    holder : str, optional
        Copyright holder written into the licence files.
    git : bool, optional
        Make the project a git repository and commit everything.

    Returns
    -------
    CreatedProject
        The folders written and the first commit.

    Raises
    ------
    ProtocolError
        If ``root`` exists and is not empty, no study is given, a
        label is given twice, or a label is not a folder name (lower case,
        digits and ``_``).
    """
    from polyzymd.analyses.project_file import PROJECT_FILE
    from polyzymd.analyses.study_git import init_repository
    from polyzymd.analyses.study_scaffold import CC_BY, MIT

    root = Path(root).expanduser().resolve()
    if root.exists() and any(root.iterdir()):
        raise ProtocolError(
            f"{root} is not empty.",
            hint="Choose a new folder: project init writes and commits everything in it.",
        )
    if not studies or len(set(studies)) != len(studies):
        raise ProtocolError(
            "A project needs at least one study, each with its own label.",
            hint="Give --study LABEL once per study, such as --study lipa363 --study rml333.",
        )
    for label in studies:
        _check_study_label(label)
    root.mkdir(parents=True, exist_ok=True)
    for folder, what in PROJECT_FOLDERS.items():
        (root / folder).mkdir(exist_ok=True)
        (root / folder / "README.md").write_text(f"# {folder}/\n\n{what}.\n")
    created = CreatedProject(root, {})
    for label in studies:
        created.studies[label] = _write_study(root / label)
    listed = "studies:                      # label: folder holding its study.yaml\n" + "".join(
        f"  {label}: {label}\n" for label in studies
    )
    (root / PROJECT_FILE).write_text(PROJECT_YAML.format(studies=listed))
    year, holder = date.today().year, holder or "the project's authors"
    (root / ".gitignore").write_text(PROJECT_GITIGNORE)
    (root / "LICENSE-code").write_text(MIT.format(year=year, holder=holder))
    (root / "LICENSE-data").write_text(CC_BY.format(year=year, holder=holder))
    (root / "README.md").write_text(
        f"# {root.name}\n\nA PolyzyMD project: the studies of one paper "
        f"({', '.join(studies)}).\n`project.yaml` lists the studies and the analyses each "
        "runs; each study's `study.yaml` holds its conditions, structures and "
        "regions.\n\n```bash\npolyzymd project check .\npolyzymd analyze --project .\n"
        "polyzymd project freeze .\n```\n"
    )
    if git:
        created.commit = init_repository(root, "Create project with polyzymd project init")
    return created


def _check_study_label(label: str) -> None:
    """Refuse a study label that is not a folder name (lower case, digits and ``_``) or is reserved."""
    from polyzymd.analyses.study_scaffold import check_label, condition_folder

    if condition_folder(label) != label:
        raise ProtocolError(
            f"The study label {label!r} is not a folder name.",
            hint=f"Use lower case, digits and _ only, such as {condition_folder(label)!r}.",
        )
    check_label(label, "study")


def _write_study(folder: Path) -> Path:
    """Write an empty study of a project into ``folder``, without licences, and return it."""
    from polyzymd.analyses.study_scaffold import create_study

    create_study(folder, git=False)
    # The project holds the licences, once.
    for name in ("LICENSE-code", "LICENSE-data"):
        (folder / name).unlink(missing_ok=True)
    (folder / "study.yaml").write_text(STUDY_YAML)
    return folder


def add_study(project: str | Path, label: str) -> Path:
    """Write an empty study ``label`` into the project at ``project`` and list it in ``project.yaml``.

    The study folder is ``<project>/<label>``, as :func:`create_project`
    writes it. The study is added as one line under ``studies:``, so the
    rest of ``project.yaml``, comments included, is kept. ``runs/`` is
    added to the project's ``.gitignore`` if it lacks it. Nothing is
    committed. Returns the study folder.

    Raises
    ------
    ProtocolError
        If there is no ``project.yaml``, the label is not a folder name or is
        reserved, or the project already lists the label or holds its folder.
    """
    import yaml

    from polyzymd.analyses.project_file import find_project_file
    from polyzymd.analyses.study_scaffold import _list_entry, ignore_runs

    file = find_project_file(project)
    _check_study_label(label)
    listed = (yaml.safe_load(file.read_text()) or {}).get("studies") or {}
    folder = file.parent / label
    if label in listed or folder.exists():
        raise ProtocolError(
            f"{file} already has a study {label!r} or a folder {folder}.",
            hint="Choose another label.",
        )
    _write_study(folder)
    ignore_runs(file.parent)
    _list_entry(file, "studies", label, label)
    return folder
