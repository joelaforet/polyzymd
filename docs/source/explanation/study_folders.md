# Study folders: publishing a reproducible MD study

A {term}`study` folder holds one MD study: a set of conditions that you
compare with each other, in one analysis frame. For what the frame is and
which conditions share it, see {doc}`projects`. The folder contains these
items:

- the simulation config of each condition;
- the analysis protocol, `study.yaml`;
- your analysis and figure code;
- the stored results.

You can version the folder with git, zip it, and publish it with the paper.
The trajectories go to Zenodo. A reader installs the PolyzyMD version that
the manifest records. The reader can then make every figure from the folder
alone, and every analysis from the folder and the trajectories.

The design follows the FAIR principles (Wilkinson et al. 2016; Barker et al.
2022) and the TRUE principles for molecular simulation (Thompson et al.
2020).

To make a study folder, see {doc}`../how_to/study_folder`.

## Three levels of reproduction

| Level | The reader has | The reader runs |
|---|---|---|
| 1. Figures | The study folder only | The scripts in `figures/`. They read stored results with `pz.Study("study.yaml").results(name)` and need no trajectories |
| 2. Analyses | The folder and the trajectories | `polyzymd study locate DIR`, then `polyzymd analyze --study study.yaml` |
| 3. Simulations | The folder and compute time | Each `conditions/<label>/config.yaml`. The replicate number seeds the starting structure and the dynamics |

Level 3 gives the same results within the statistical noise of MD, not bit
for bit. Floating-point arithmetic, the order of parallel sums and the
hardware differ between machines (Thompson et al. 2020).

## Layout

```
my_study/
├── study.yaml                 # the analysis protocol (committed)
├── data.example.yaml          # explains data.local.yaml (committed)
├── data.local.yaml            # where this machine keeps the trajectories (git ignores it)
├── README.md                  # generated: how to reproduce at each level, how to cite
├── LICENSE-data, LICENSE-code # CC-BY-4.0 and MIT by default; replace them with your own
├── conditions/<label>/
│   ├── config.yaml            # the simulation config of the condition
│   └── structures/            # the input structures of that condition
├── structures/                # analysis references (crystal, catalytically competent frame)
├── analyses/                  # measurement functions: Universe -> value
├── figures/                   # notebooks and scripts: stored results -> paper figures
├── results/                   # stored per-replicate values, reports, default figures
├── environment/               # how to install the PolyzyMD version; add pixi.toml and pixi.lock
└── .gitignore
```

`polyzymd study init DIR` writes this layout, runs `git init` and makes the
first commit. For each condition, `study init` and `study add-condition` copy
the config into `conditions/<label>/`, with the input files that it names.
Conditions that share a structure each get a copy. These copies are small
compared with the trajectories.

## `study.yaml`, the analysis protocol

`config.yaml` describes one simulation. `study.yaml` holds the rest of what
you need to make the figures of the paper from the trajectories:

- the conditions;
- the equilibration window;
- the settings of each analysis.

Analysis settings change much more often than simulation settings. So they
have their own file.

```yaml
equilibration: 100ns            # one window for every condition and replicate
stride: 1                       # optional
replicates: [1, 2, 3, 4, 5]     # optional; default: every replicate found
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
    kind: timeseries            # or per_replicate
    unit: A
    selections: {lid: "resid 140-150 and name CA", core: "resid 4-120 and name CA"}
    settings: {cutoff: 8.0}
metadata: {}                    # see "Publishing metadata"
```

The rules of the file, and the reason for each:

| Rule | Reason |
|---|---|
| Each analysis uses the window of the study, unless it sets its own | Replicates and conditions are compared on equal terms. A time-resolved analysis can still start where a steady-state analysis cannot |
| Command-line options override `study.yaml`, and `study.yaml` overrides the defaults | A one-time change needs no edit |
| An unknown key is an error with a "did you mean" hint | A misspelled setting must not fall back to its default |
| Results go to `results/` beside `study.yaml`, unless you give `--output-dir` | The protocol and its results stay together |
| Each `selections:` entry becomes an `AtomGroup` keyword argument for each replicate. `settings:` are passed unchanged | Your function keeps the form `compute(groups, **settings)` |
| The stored results of your function depend on every file in the folder of the function | A change to a helper that the function calls must recompute its results |

| Command | What it does |
|---|---|
| `polyzymd analyze NAME --study study.yaml` | Runs one listed analysis, shipped or your own |
| `polyzymd analyze --study study.yaml` | Runs every listed analysis. With `--submit`, it submits one SLURM array for each analysis |
| `polyzymd study check` | Checks paths, configs and setting names, and imports every function. It loads no trajectory. It prints where it found the replicates of each condition |

## Where the trajectories are: `data.local.yaml`

The path of a trajectory points to where the file is now. It is not part of
the study. When you move data from cluster storage, or download it from
Zenodo, the study and its results must not change. So machine paths stay out
of `study.yaml`:

```yaml
# data.local.yaml: where this machine keeps the replicates of each condition (git ignores it)
No polymer: /data/me/polyzymd_sims/LipA_363K
SBMA 50%: /data/me/polyzymd_sims/LipA_363K
```

| Situation | What PolyzyMD reads |
|---|---|
| No entry in `data.local.yaml` for a condition | The `scratch_directory` of the config |
| An entry in `data.local.yaml` | The folder of that entry |
| `--data DIR` | `DIR`, for this command only |
| Downloaded trajectories | `polyzymd study locate DIR` finds the replicates of each condition under `DIR`. It prefers files whose sizes match the freeze manifest (`--verify` also checks the SHA-256). It writes `data.local.yaml` |

Other tools separate what a project is from where its files are on one
machine in the same way:

- DVC keeps content hashes in committed files, and machine locations in an
  uncommitted `config.local`.
- Git LFS commits a pointer that holds a content hash.
- The twelve-factor convention keeps deployment settings in a `.env` file
  that git ignores, beside a committed `.env.example`.

## Identity rests on content

PolyzyMD reuses a stored result only when its inputs are unchanged. This test
does not depend on where the files are. So a moved or downloaded study keeps
its results, and you can check its numbers against the published numbers.

| Input | Identified by |
|---|---|
| The config of a condition | Its simulation content. Input structures count by the SHA-256 of their content. The projects and scratch directories do not count |
| A shipped function | Its module file, and the PolyzyMD modules that the file imports |
| Your function | Every file in the folder of its file |
| Trajectory and topology files | The path relative to the folder of the replicate folders, the size and the SHA-256 |
| A file that an analysis reads, such as a reference structure | Its name and its SHA-256, not its location |

PolyzyMD gets the SHA-256 of a trajectory from the first of these sources:

1. `progress.json`. The simulation runner records the hash there when each
   production segment finishes.
2. `trajectory_hashes.json`. `polyzymd hash-trajectories` records the hash
   there for simulations whose runner did not.
3. A cache in `~/.cache/polyzymd/hashes`, or `$POLYZYMD_CACHE_DIR/hashes`.
   PolyzyMD computes the hash once and keeps it while the size and the
   modification time of the file do not change.

So PolyzyMD hashes copied or downloaded data once, and reuses its stored
results.

## Git and provenance

- PolyzyMD decides reuse by content, never by commit. A commit that fixes a
  figure script recomputes nothing.
- Each stored record and report records the commit of the study, and the
  uncommitted files, if there are any.
- During analysis, uncommitted changes give a warning, never an error.
- `freeze` publishes, so it refuses uncommitted inputs.
- PolyzyMD commits for you only in `study init` and in `freeze`.

## Publishing: `polyzymd study freeze`

`freeze` prepares the folder for deposit. It does these steps:

1. It refuses uncommitted inputs, and names them. It warns about each listed
   analysis whose stored results do not match its current code and settings.
2. It collects these items without user input:
   - the engine inputs of each replicate (OpenMM `system.xml`, or GROMACS
     `.top`, `.mdp` and `.tpr`). They record the parameters that the toolkit
     assigned, generated polymer parameters included (Thompson et al. 2020);
   - the final-frame coordinates of each replicate (Reliability and
     reproducibility checklist 2023);
   - the force-field names and versions, the charge method and its toolkit
     version, the cutoffs and the long-range method;
   - the platform, the precision and the production length of each
     replicate, and any post-processing such as a stride or stripped solvent
     (Tiemann et al. 2024);
   - a system-setup table, `system_summary.csv`: the box, the atom count, the
     waters, the salt and the composition (Reliability and reproducibility
     checklist 2023).
3. It writes `manifest.json`: the file hashes and sizes, the package versions,
   the provenance above, the trajectory DOIs and the parent of the tagged
   commit.
4. It writes `md_checklist.yaml`, the reliability and reproducibility
   checklist of Communications Biology (2023), filled in from the manifest.
   The checklist is for information. You can send it with a journal
   submission.
5. It writes `CITATION.cff` and `.zenodo.json` from the one `metadata:`
   block, so the two files always agree. Zenodo ignores `CITATION.cff` when
   `.zenodo.json` exists.
6. It commits these files and `results/`, and tags the commit.
7. It lays out `deposit/`. The manifest, the README and `CITATION.cff` are at
   the top level and not zipped, so Zenodo indexes them and shows a preview
   (Tiemann et al. 2024).

`freeze` prepares the upload but does not upload. A Zenodo publication is
permanent and gets a DOI, so the author does that step. `freeze` writes these
files:

- `deposit/upload/`: the files to add to one Zenodo record, within its limits
  (100 files and 50 GB by default). The study, the engine inputs and the
  final frames are one zip file each.
- `deposit/trajectories.csv`: the trajectory files, in batches that each fit
  one record.
- `deposit/UPLOAD.md`: the steps for this study, from the reservation of its
  DOI (`metadata.doi`) to publication.

You can also script the upload with the Zenodo REST API.

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
| Persistent identifier (FAIR F1, F3) | The version DOI and concept DOI of the study on Zenodo, written into the manifest and both citation files |
| Rich metadata, purpose, license (F2, R1, R1.1) | `metadata:`, with SPDX license IDs and LICENSE files |
| Detailed provenance (R1.2) | `manifest.json`, the git tag, the stored records |
| Community standards (R1.3) | `md_checklist.yaml`; the MDverse sharing guidelines (Tiemann et al. 2024) |
| Qualified references (I3) | The `.zenodo.json` relations: `isSupplementTo` the paper, `requires` PolyzyMD, `references` the trajectories |
| All engine inputs, versioned, scripted (TRUE) | `conditions/`, the engine inputs, `analyses/`, `figures/`, `environment/`, git |
| Statistical rather than exact reproducibility (TRUE) | The replicate numbers, which seed the starting structures, and the platform and precision in the manifest |

## Citing PolyzyMD

Each published study cites PolyzyMD. Its authors get credit, and an agent
that reads the study knows which software produced it.

| Where | What |
|---|---|
| `CITATION.cff` | The paper of the study as `preferred-citation`. PolyzyMD (software) and the PolyzyMD paper under `references`, copied from the `CITATION.cff` of the installed package |
| `.zenodo.json` | A `requires` relation to the DOI of PolyzyMD |
| Generated `README.md` | A "How to cite" section |
| Reports and figures | The PolyzyMD version line and the figure watermark |
| `polyzymd study check` | One line that names PolyzyMD and how to cite it |

While a DOI is a placeholder, for example for an unpublished paper, `freeze`
warns. When you know the DOI, freeze again. The citations are then written
again under the next tag. When the `CITATION.cff` of PolyzyMD changes, for
example when its paper is published, every study frozen after that cites the
new entry.

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
