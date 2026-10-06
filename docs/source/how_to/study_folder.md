# Create a study folder

A {term}`study` folder holds one protein (or other system) under its
conditions. It holds the config of each condition, the analysis protocol,
your analysis and figure code, and the stored results. You can version it,
share it and publish it. To write `study.yaml` and run the analyses, see
{doc}`study_yaml`. For the reasons behind the layout, see
{doc}`../explanation/study_folders`.

:::{admonition} Environment Setup
:class: tip

The commands on this page assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## Create the folder

Give `polyzymd study init` the config of each condition that you set up or
ran. Give the control first.

```bash
polyzymd study init lipase_363K \
  --condition "No polymer=runs/noPoly/config.yaml" \
  --condition "SBMA 50%=runs/SBMA50/config.yaml" \
  --equilibration 100ns
```

```
created /home/me/lipase_363K
condition No polymer: conditions/no_polymer/config.yaml, with 1 input files copied to structures/
condition SBMA 50%: conditions/sbma_50/config.yaml, with 1 input files copied to structures/
git: committed 53cd9f37d690
next: polyzymd study check /home/me/lipase_363K
```

For each `--condition`, `study init` does these steps:

1. It copies the config to `conditions/<name>/config.yaml`.
2. It copies each input file that the config names (the protein PDB, a
   substrate SDF, polymer building blocks) to `conditions/<name>/structures/`.
   The paths in the copied config are relative.
3. It removes the machine paths from the copied config. `projects_directory`
   becomes `.`, and a `scratch_directory` becomes `data`. These paths are not
   part of the {term}`config hash`, so stored results still match.
4. It writes the folder that holds the simulations of the condition to
   `data.local.yaml`.

The original configs do not change.

For a condition that you have not simulated yet, use
`--new-condition "No polymer"`. It creates `conditions/no_polymer/` with
`polyzymd init`, for you to fill in and simulate.

### Add a condition later

```bash
polyzymd study add-condition "SBMA 100%" --config runs/SBMA100/config.yaml --study lipase_363K
polyzymd study add-condition "SBMA 100%" --new --study lipase_363K   # an empty simulation folder
```

`add-condition --config` copies the config and its inputs as `study init`
does, and writes the location of the simulations to `data.local.yaml`. Both
forms add one line under `conditions:` in `study.yaml`. The rest of the
file, with its comments, does not change.

### What the folder holds

| Path | Holds |
|---|---|
| `study.yaml` | The analysis protocol: conditions, equilibration window, analyses |
| `conditions/<name>/` | The `config.yaml` and `structures/` of each condition |
| `structures/` | Reference structures for analyses, such as a crystal structure |
| `analyses/` | Your measurement functions, listed in `study.yaml` |
| `figures/` | Notebooks and scripts that make figures from stored results |
| `results/` | What `polyzymd analyze --study` stores |
| `environment/` | The PolyzyMD version and how to install it |
| `README.md` | How to reproduce the study and how to cite it |
| `LICENSE-data`, `LICENSE-code` | CC-BY-4.0 and MIT by default |
| `data.example.yaml` | How to write `data.local.yaml` |
| `data.local.yaml` | Where this machine keeps the simulations. Git ignores it |
| `.gitignore` | Leaves out `data.local.yaml`, caches and SLURM logs |

`study init` also makes the folder a git repository with one commit. Use
`--no-git` to skip git. To use a different license, replace a LICENSE file.
`--holder NAME` sets the copyright holder. The default is your
`git config user.name`.

## Say where the trajectories are

The path of a trajectory says where the file is now. It is not part of the
study. So the path goes in `data.local.yaml`, which git ignores and freeze
never publishes.

| Situation | What to do |
|---|---|
| You made the study with `study init --condition` or `add-condition --config` | Nothing. The command wrote `data.local.yaml` |
| A condition has no entry in `data.local.yaml` | Nothing. PolyzyMD uses the `scratch_directory` of the config |
| You moved the data, for example to PetaLibrary | Edit `data.local.yaml`, or run `polyzymd study locate NEW_DIR` |
| You downloaded the trajectories, for example from Zenodo | Run `polyzymd study locate DOWNLOAD_DIR` |
| You want one command to read a different copy | Run `polyzymd analyze RUN --study S --data DIR` |

```bash
polyzymd study locate ~/Downloads/zenodo_1234567 --study lipase_363K
```

```
No polymer: runs [1, 2, 3, 4, 5] under /home/me/Downloads/zenodo_1234567/no_polymer
SBMA 50%: runs [1, 2, 3, 4, 5] under /home/me/Downloads/zenodo_1234567/sbma_50
wrote /home/me/lipase_363K/data.local.yaml
```

`data.local.yaml` maps each condition label to the folder that holds its
{term}`replicate folders <replicate folder>`:

```yaml
No polymer: /pl/active/my_lab/polyzymd_sims/LipA_363K
SBMA 50%: /pl/active/my_lab/polyzymd_sims/LipA_363K
```

`polyzymd study check` prints where it found the simulations of each
condition. It also says whether the location came from `data.local.yaml` or
from the config.

When you move data, the config hash of a condition does not change. The hash
identifies input structures by their content, and it leaves out the
projects and scratch directories. Stored records name each trajectory file
relative to the folder that holds the replicate folders.

Stored records identify each trajectory and topology file by its SHA-256 and
size. So moved, copied or downloaded data reuses every stored result. At a
new location, the first analysis reads each file once to hash it. This takes
about one second per gigabyte. It is not necessary when the hashes of the
replicate are recorded.

### Record the hashes of older simulations

Some simulations have no recorded hashes:

- OpenMM simulations that finished before PolyzyMD recorded hashes;
- downsampled copies;
- GROMACS simulations.

Each machine that analyzes them hashes them again. To prevent this, record
the hashes once, on the machine that holds the data. The command writes only
`trajectory_hashes.json` beside the files of each replicate.

```bash
polyzymd hash-trajectories --study lipase_363K --dry-run   # what it would hash
polyzymd hash-trajectories --study lipase_363K             # record them
```

```
No polymer replicate 1 (openmm): hashed 3, already recorded 0
SBMA 50% replicate 1 (gromacs): hashed 2, already recorded 0
```

- A second run prints `hashed 0` and changes nothing.
- The command never overwrites a recorded hash. If a hash does not agree, it
  prints `conflict:` and exits with code 2.
- `--verify` hashes the files again and compares.
- `--rehash-changed` records the hash again for a file that grew after it
  was hashed.
- `-c config.yaml` works for simulations outside a study.

On a cluster, run the command in a batch job, because it reads every
trajectory once. See {doc}`../reference/cli_reference`.

## Commit as you go

Each report of `polyzymd analyze --study` records, under `provenance.study`
in `results/<run>/report.json`:

- the SHA-256 of `study.yaml`;
- the git commit of the folder;
- the uncommitted files, if there are any.

An uncommitted input, such as an edited `study.yaml` or analysis function,
gives a warning. Changes in `results/` and `data.local.yaml` do not, because
they are outputs and pointers.

```
warning: the study has uncommitted changes (study.yaml); the report records them, and committing them makes it reproducible.
```

PolyzyMD reuses stored results by content, never by commit. So a commit never
causes a recomputation. Before you freeze, PolyzyMD commits only once, in
`study init`. After that, you decide what to commit. `polyzymd study check` prints the commit and
each uncommitted input.

Commit every input before you publish. `polyzymd study freeze` stops while an
input is not committed; see {doc}`study_freeze`.
