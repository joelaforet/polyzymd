# Study folders: publishing a reproducible MD study

```{admonition} Status: design, partly implemented
:class: warning
This page is the agreed design of PolyzyMD study folders. It is implemented in
slices: A (the protocol file), B (the folder), C (publishing) and D (upload to
Zenodo), plus a small slice that records trajectory hashes when segments
finish. Each section names its slice. Slices A, B and C are implemented:
`study.yaml`, `polyzymd analyze --study`, `polyzymd study check` and
`Study.results` ({doc}`../how_to/study_yaml`); `polyzymd study init`,
`data.local.yaml`, `polyzymd study locate`, `--data` and git provenance
({doc}`../how_to/study_folder`); and `polyzymd study freeze`
({doc}`../how_to/study_freeze`), which also prepares the upload to Zenodo
and writes the steps (slice D: PolyzyMD uploads and publishes nothing). The
segment hashes recorded when simulations finish and identity by content do
not exist yet.
```

A **study folder** holds one MD study: every condition's simulation config,
the analysis protocol, the analysis and figure code, and the stored results.
You can version it with git, zip it, and publish it with the paper, with the
trajectories deposited on Zenodo. Someone who installs the PolyzyMD version
the study names can then reproduce every reported figure from the folder
alone, and every analysis from the folder plus the trajectories. The design
follows the FAIR principles (Wilkinson et al. 2016; Barker et al. 2022) and the
TRUE principles for molecular simulation (Thompson et al. 2020).

## Three levels of reproduction

| Level | The reproducer has | Command | Slice |
|---|---|---|---|
| 1. Figures | The study folder only | Run the scripts in `figures/`, which read stored results with `pz.Study("study.yaml").results(name)` and need no trajectories | A |
| 2. Analyses | The folder plus the trajectories | `polyzymd study locate DIR`, then `polyzymd analyze --study study.yaml` | A, B |
| 3. Simulations | The folder plus compute | Each `conditions/<label>/config.yaml`, with the replicate number as the random seed | existing |

Level 3 reproduces results within the statistical noise of MD, not bit for
bit: floating-point arithmetic, parallel reduction order and hardware differ
between machines (Thompson et al. 2020).

## Layout

```
my_study/
├── study.yaml                 # the analysis protocol (committed)
├── data.example.yaml          # explains data.local.yaml (committed)
├── data.local.yaml            # where this machine keeps the trajectories (gitignored)
├── README.md                  # generated: how to reproduce at each level, how to cite
├── LICENSE-data, LICENSE-code # CC-BY-4.0 and MIT by default; replace with your own
├── conditions/<label>/
│   ├── config.yaml            # the simulation, as written by polyzymd init
│   └── structures/            # that condition's input structures
├── structures/                # analysis references (crystal, catalytically competent frame)
├── analyses/                  # measurement functions: Universe -> value
├── figures/                   # notebooks and scripts: stored results -> paper figures
├── results/                   # stored per-replicate values, reports, default figures
├── environment/               # pixi.toml and pixi.lock pinning PolyzyMD and its stack
└── .gitignore
```

`polyzymd study init DIR` (slice B) writes this layout, runs `git init` and
makes the first commit. Conditions point at their configs in place; nothing
is copied. Structures are duplicated across conditions that share them,
which costs little next to the trajectories.

## `study.yaml`, the analysis protocol (slice A)

`config.yaml` describes one simulation. `study.yaml` holds everything else
needed to regenerate the paper's figures from the trajectories: the
conditions, the equilibration window, and the settings of every analysis.
Analysis settings change far more often than simulation settings, which is
why they live in their own file.

```yaml
polyzymd: 1.3.0                 # the version that produced the results
equilibration: 100ns            # one window for every condition and replicate
stride: 1                       # optional
replicates: [1, 2, 3, 4, 5]     # optional; default: every run found
conditions:                     # control first; paths relative to this file
  No polymer: conditions/no_polymer/config.yaml
  SBMA 50%: conditions/sbma50/config.yaml
analyses:
  contacts:                     # a shipped analysis, named by its key
    method: occlusion
  contacts_4A:                  # the same analysis with other settings
    analysis: contacts
    method: distance
    cutoff: 4.0
  lid_opening:                  # your own function
    function: analyses/lid.py:lid_distance
    kind: per_replicate         # or timeseries
    unit: A
    selections: {lid: "resid 140-150 and name CA", core: "resid 4-120 and name CA"}
    settings: {cutoff: 8.0}
metadata: {}                    # see "Publishing metadata"
```

| Rule | Why |
|---|---|
| One equilibration window for the whole study | Replicates and conditions are compared on equal footing |
| Command-line options override `study.yaml`, which overrides defaults | A one-off change needs no edit |
| Unknown keys are errors with a "did you mean" hint | A misspelt setting must not fall back to its default silently |
| Results go to `results/` beside `study.yaml` unless `--output-dir` says otherwise | The protocol and its results travel together |
| Each `selections:` entry becomes an `AtomGroup` keyword argument built for each replicate; `settings:` are passed as they are | User functions keep the `compute(universe_or_groups, **settings)` form people already write |
| A user function's stored results are keyed on the hash of its whole module file | Editing a helper the function calls must recompute its results |

| Command | Does |
|---|---|
| `polyzymd analyze NAME --study study.yaml` | Runs one listed analysis, shipped or user-defined |
| `polyzymd analyze --study study.yaml` | Runs every listed analysis; with `--submit`, one SLURM array per analysis |
| `polyzymd study check` | Validates paths, configs, setting names and that every function imports, without loading trajectories; prints where each condition's runs were found |

## Where the trajectories are: `data.local.yaml` (slice B)

A trajectory's path is a pointer to where the file lives now, not part of the
study. Moving data from cluster storage to a Zenodo download must not change
the study or invalidate its results, so machine paths stay out of
`study.yaml`:

```yaml
# data.local.yaml: where this machine keeps each condition's runs (gitignored)
No polymer: /pl/active/shirts_archive/LaforetJoe/polyzymd_sims/LipA_363K_REDO
SBMA 50%: /pl/active/shirts_archive/LaforetJoe/polyzymd_sims/LipA_363K_REDO
```

| Situation | What PolyzyMD reads |
|---|---|
| No `data.local.yaml` (someone actively doing the science) | Each config's own `scratch_directory`, as today |
| `data.local.yaml` present | Its directory for each condition it names |
| `--data DIR` | `DIR` for this command only |
| Downloaded trajectories | `polyzymd study locate DIR` finds each condition's runs under `DIR`, checks them against the freeze manifest, and writes `data.local.yaml` |

This follows tools that separate what a project is from where its files sit
on one machine: DVC keeps content hashes in committed files and machine
locations in an uncommitted `config.local`; Git LFS commits a pointer holding
a content hash; the twelve-factor convention keeps deployment settings in a
gitignored `.env` beside a committed `.env.example`.

## Identity rests on content

Stored results are reused only when what produced them is unchanged, and
that test must not depend on where the files sit, so a moved or downloaded
study keeps its results and its numbers can be checked against the published
ones.

| Input | Identified by | Status |
|---|---|---|
| A condition's config | Its simulation content: input structures by the SHA-256 of their content; projects and scratch directories left out | Implemented |
| A function | Its source, or its whole module file for user functions | Implemented |
| Trajectory and topology files | Their path relative to the folder holding the runs, size and modification time; the freeze manifest adds their SHA-256 | Implemented; records still include modification time |
| Trajectory files by content in records | SHA-256 and size, cached by path, size and modification time, so copied or downloaded data is not measured again | Planned, with the runner slice |

The simulation runner will record each segment's SHA-256 when the segment
finishes (its own small slice). That is bookkeeping only: trajectories and
results do not change.

## Git and provenance (slice B)

- Reuse is decided by content, never by commit, so committing a fix to a
  figure script recomputes nothing.
- Every stored record and report records the study's commit and whether the
  working tree had uncommitted changes, and which files.
- Uncommitted changes give a warning, never a refusal.
- PolyzyMD never commits for you after `study init`.

## Publishing: `polyzymd study freeze` (slice C)

`freeze` prepares the folder for deposit:

1. Checks that the working tree is clean and that every listed analysis has
   stored results matching its current code and settings; warns otherwise.
2. Collects, without user input:
   - the serialized engine inputs of every replicate (OpenMM `system.xml`, or
     GROMACS `.top`, `.mdp` and `.tpr`), which freeze the parameters the
     toolkit actually assigned, including generated polymer parameters
     (Thompson et al. 2020);
   - each replicate's final-frame coordinates (Reliability and reproducibility checklist 2023);
   - force-field names and versions, the charge method and its toolkit
     version, cutoffs and the long-range method;
   - platform, precision, production length and frames analysed for each
     replicate, and any post-processing such as a stride or stripped solvent
     (Tiemann et al. 2024);
   - a system-setup table: box, atom count, waters, salt and composition
     (Reliability and reproducibility checklist 2023).
3. Writes `manifest.json`: file hashes and sizes, package versions, the
   provenance above, trajectory DOIs and the git commit.
4. Writes `md_checklist.yaml`, the Communications Biology reliability and
   reproducibility checklist (2023), filled in from the manifest.
   It is informational and can accompany a journal submission.
5. Writes `CITATION.cff` and `.zenodo.json` from the one `metadata:` block, so
   they cannot disagree (Zenodo ignores `CITATION.cff` when `.zenodo.json`
   exists).
6. Tags the commit and lays out `deposit/`, with the manifest, README and
   `CITATION.cff` at the top level, unzipped so that they stay indexed and
   previewable (Tiemann et al. 2024).

Slice D prepares the upload rather than doing it: publishing on Zenodo is
permanent and mints a DOI, so it stays the author's step. `freeze` writes
`deposit/upload/`, exactly the files to add to a Zenodo record within its
limits (100 files and 50 GB by default), with the study, engine inputs and
final frames as one zip each; `deposit/trajectories.csv`, the trajectory
files grouped into batches that each fit a record; and `deposit/UPLOAD.md`,
the steps for this study from reserving its DOI (`metadata.doi`) to
publishing. Zenodo's REST API remains available to anyone who wants to
script the upload.

### Publishing metadata

```yaml
metadata:
  title: "..."
  description: "..."
  purpose: "..."                  # the main purpose of the simulations
  keywords: [molecular dynamics, enzyme, polymer]
  system_type: [protein, polymer]
  authors:
    - {family-names: ..., given-names: ..., orcid: "https://orcid.org/...", affiliation: "..."}
  license: {data: CC-BY-4.0, code: MIT}   # SPDX identifiers
  funding: [{funder: "...", award: "..."}]
  related:
    paper: {doi: "10.XXXX/...", status: in-preparation, title: "..."}
    trajectories: [{doi: "10.5281/zenodo.NNNN", conditions: [...]}]
  zenodo: {communities: [], access_right: open}
```

| FAIR or TRUE requirement | Met by |
|---|---|
| Persistent identifier (FAIR F1, F3) | The study's Zenodo version and concept DOIs, written into the manifest and both citation files |
| Rich metadata, purpose, licence (F2, R1, R1.1) | `metadata:`, with SPDX licence IDs and LICENSE files |
| Detailed provenance (R1.2) | `manifest.json`, the git tag, the stored records |
| Community standards (R1.3) | `md_checklist.yaml`; MDverse sharing guidelines (Tiemann et al. 2024) |
| Qualified references (I3) | `.zenodo.json` relations: `isSupplementTo` the paper, `requires` PolyzyMD, `references` the trajectories |
| All engine inputs, versioned, scripted (TRUE) | `conditions/`, serialized engine inputs, `analyses/`, `figures/`, `environment/`, git |
| Statistical rather than exact reproducibility (TRUE) | Seeds (replicate numbers), platform and precision in the manifest |

## Citing PolyzyMD

Every published study cites PolyzyMD, so its authors get credit and agents
reading the study know which framework produced it:

| Where | What |
|---|---|
| `CITATION.cff` | The study's paper as `preferred-citation`; PolyzyMD (software) and the PolyzyMD paper under `references`, copied from the installed package's own `CITATION.cff` |
| `.zenodo.json` | A `requires` relation to PolyzyMD's DOI |
| Generated `README.md` | A "How to cite" section |
| Reports and figures | The PolyzyMD version line and figure watermark |
| `polyzymd study check` | One line naming PolyzyMD and how to cite it |

While a DOI is still a placeholder, such as an unpublished paper, `freeze`
warns; refreezing once it is known writes the citations again under the next
tag. Updating PolyzyMD's own `CITATION.cff` when its paper is
published updates every study frozen afterwards.

## Sources

- Wilkinson, M. D. et al. (2016). The FAIR Guiding Principles for scientific
  data management and stewardship. Scientific Data 3:160018.
  doi:10.1038/sdata.2016.18
- Barker, M. et al. (2022). Introducing the FAIR Principles for research
  software. Scientific Data 9:622. doi:10.1038/s41597-022-01710-x
- Thompson, M. W. et al. (2020). Towards molecular simulations that are
  transparent, reproducible, usable by others, and extensible (TRUE).
  Molecular Physics 118:e1742938. doi:10.1080/00268976.2020.1742938
- Reliability and reproducibility checklist for molecular dynamics
  simulations (2023). Communications Biology 6:268.
  doi:10.1038/s42003-023-04653-0 (an unsigned editorial, so cited by its
  title)
- Tiemann, J. K. S. et al. (2024). MDverse, shedding light on the dark matter
  of molecular dynamics simulations. eLife 12:RP90061.
  doi:10.7554/eLife.90061
- Amaro, R. et al. (2025). The need to implement FAIR principles in
  biomolecular simulations. Nature Methods 22:641-645.
  doi:10.1038/s41592-025-02635-0 (its minimum metadata, including the
  purpose of the simulations, shaped the `metadata:` block)
- Citation File Format 1.2.0. doi:10.5281/zenodo.1003149
