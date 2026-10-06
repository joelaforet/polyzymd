"""Create a project folder: the layout ``polyzymd project init`` writes.

A project holds one paper: ``project.yaml``, the shared ``analyses/``,
``stats/`` and ``figures/`` folders, and one study folder per protein, each
with a ``study.yaml`` to fill in. See the "Projects and studies" explanation
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
# chapter. Each study below is one protein (or other system) with its own
# conditions, structures and regions; the analyses listed here run in every
# study. See https://polyzymd.readthedocs.io/en/latest/how_to/project.html

{studies}
analyses: {{}}                  # run in every study, with that study's regions and structures
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

metadata:                     # TODO before polyzymd project freeze
  title: TODO
  authors:
    - {{family-names: TODO, given-names: TODO, orcid: TODO, affiliation: TODO}}
  license: {{code: MIT, data: CC-BY-4.0}}
  keywords: []
"""

STUDY_YAML = """\
# One protein (or other system) of the project: its conditions, structures
# and regions. Analyses come from ../project.yaml; add this protein's own
# analyses under analyses:.
description: TODO             # e.g. B. subtilis lipase A (1ISP) at 363 K
equilibration: 0ns            # TODO: the burn-in to discard from every replicate

structures: {}                # name: file under structures/, used as 'structure <name>'
# structures:
#   reference: structures/1ISP_clean.pdb

regions: {}                   # name: selection, used as 'region <name>'
# regions:
#   core: resid 5-8 15-27 32-37
#   catalytic_triad: resid 76 132 155

conditions: {}                # control first; add with polyzymd study add-condition
# conditions:
#   No Polymer: conditions/no_polymer
#   SBMA-EGMA 50:50:
#     config: conditions/sbma_50
#     factors: {sbma_fraction: 0.5}           # optional: what varies, for trend tests and plots

analyses: {}                  # optional: analyses only this protein runs
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
        Study labels, one per protein; each is also its folder name.
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
    from polyzymd.analyses.study_scaffold import CC_BY, MIT, condition_folder, create_study

    root = Path(root).expanduser().resolve()
    if root.exists() and any(root.iterdir()):
        raise ProtocolError(
            f"{root} is not empty.",
            hint="Choose a new folder: project init writes and commits everything in it.",
        )
    if not studies or len(set(studies)) != len(studies):
        raise ProtocolError(
            "A project needs at least one study, each with its own label.",
            hint="Give --study LABEL once per protein, such as --study lipa363 --study rml333.",
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
        (root / folder / "README.md").write_text(f"# {folder}/\n\n{what}.\n")
    created = CreatedProject(root, {})
    for label in studies:
        folder = root / label
        create_study(folder, git=False)
        # The project holds the licences, once.
        for name in ("LICENSE-code", "LICENSE-data"):
            (folder / name).unlink(missing_ok=True)
        (folder / "study.yaml").write_text(STUDY_YAML)
        created.studies[label] = folder
    listed = "studies:                      # label: folder holding its study.yaml\n" + "".join(
        f"  {label}: {label}\n" for label in studies
    )
    (root / PROJECT_FILE).write_text(PROJECT_YAML.format(studies=listed))
    year, holder = date.today().year, holder or "the project's authors"
    (root / ".gitignore").write_text(PROJECT_GITIGNORE)
    (root / "LICENSE-code").write_text(MIT.format(year=year, holder=holder))
    (root / "LICENSE-data").write_text(CC_BY.format(year=year, holder=holder))
    (root / "README.md").write_text(
        f"# {root.name}\n\nA PolyzyMD project: one paper, one study per protein "
        f"({', '.join(studies)}).\n`project.yaml` lists the studies and the analyses each "
        "runs; each study's `study.yaml` holds that protein's conditions, structures and "
        "regions.\n\n```bash\npolyzymd project check .\npolyzymd analyze --project .\n"
        "polyzymd project freeze .\n```\n"
    )
    if git:
        created.commit = init_repository(root, "Create project with polyzymd project init")
    return created
