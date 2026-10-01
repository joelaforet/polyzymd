# Create a study folder

Use this to start a study you will version, share and publish: one folder
with every condition's config, the analysis protocol, your analysis and
figure code, and the stored results. For writing `study.yaml` and running
the analyses, see {doc}`study_yaml`; for why the folder is laid out this way,
see {doc}`../explanation/study_folders`.

:::{admonition} Environment Setup
:class: tip

The commands on this page assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## Create it

From simulations you have already set up or run:

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

Each config is copied to `conditions/<name>/config.yaml`, and every input
file it names (the protein PDB, a substrate SDF, polymer building blocks) to
`conditions/<name>/structures/`, with the path in the config made relative.
The config's scratch and projects directories are kept: they say where the
trajectories are. The originals are not changed.

For a new study, `--new-condition "No polymer"` creates
`conditions/no_polymer/` with `polyzymd init` instead, for you to fill in
and simulate.

To add a condition later:

```bash
polyzymd study add-condition "SBMA 100%" --config runs/SBMA100/config.yaml --study lipase_363K
polyzymd study add-condition "SBMA 100%" --new --study lipase_363K   # a polyzymd init project
```

It copies the config and its inputs as `study init` does, and adds one line
under `conditions:` in `study.yaml`, leaving the rest of the file, comments
included, as it was.

| Path | Holds |
|---|---|
| `study.yaml` | The analysis protocol: conditions, equilibration window, analyses |
| `conditions/<name>/` | Each condition's `config.yaml` and `structures/` |
| `structures/` | Reference structures for analyses, such as a crystal structure |
| `analyses/` | Your measurement functions, listed in `study.yaml` |
| `figures/` | Notebooks and scripts that turn stored results into figures |
| `results/` | What `polyzymd analyze --study` stores |
| `environment/` | The PolyzyMD version and how to install it |
| `README.md` | How to reproduce the study and how to cite it |
| `LICENSE-data`, `LICENSE-code` | CC-BY-4.0 and MIT by default |
| `data.example.yaml` | How to say where this machine keeps the trajectories |
| `.gitignore` | Leaves out `data.local.yaml`, caches and SLURM logs |

`study init` also makes the folder a git repository with one commit. Replace
either LICENSE file to use another licence. `--holder NAME` sets the
copyright holder, by default your `git config user.name`.

## Say where the trajectories are

A trajectory's path is a pointer to where the file is now, not part of the
study, so it lives in `data.local.yaml`, which is gitignored and never
published:

| Situation | Do |
|---|---|
| You ran the simulations, and the configs point at the data | Nothing: each config's `scratch_directory` is used |
| You moved the data, for example to PetaLibrary | Write `data.local.yaml`, or run `polyzymd study locate NEW_DIR` |
| You downloaded the trajectories, for example from Zenodo | `polyzymd study locate DOWNLOAD_DIR` |
| One command against another copy | `polyzymd analyze RUN --study S --data DIR` |

```bash
polyzymd study locate ~/Downloads/zenodo_1234567 --study lipase_363K
```

```
No polymer: runs [1, 2, 3, 4, 5] under /home/me/Downloads/zenodo_1234567/no_polymer
SBMA 50%: runs [1, 2, 3, 4, 5] under /home/me/Downloads/zenodo_1234567/sbma_50
wrote /home/me/lipase_363K/data.local.yaml
```

`data.local.yaml` maps each condition label to the directory holding its run
directories:

```yaml
No polymer: /pl/active/my_lab/polyzymd_sims/LipA_363K
SBMA 50%: /pl/active/my_lab/polyzymd_sims/LipA_363K
```

`polyzymd study check` then says, for each condition, where its runs were
found and whether that came from `data.local.yaml` or from its config.
Moving data does not change a condition's config hash, which is that of the
config as written.

Moving data never changes a condition's config hash: input structures are
identified by their content, and the projects and scratch directories are not
part of it. Stored records name their trajectory files relative to the
folder holding the runs.

Stored records identify each trajectory and topology file by its SHA-256 and
size, so moved, copied or downloaded data reuses every stored result. The
first analysis at a new location reads each file once to hash it, about a
second per gigabyte, unless the run recorded the hash in `progress.json`.

### Record the hashes of older runs

Runs that finished before PolyzyMD recorded segment hashes have none in
`progress.json`, so each machine that analyses them hashes them again. Record
them once, where the runs are:

```bash
polyzymd hash-trajectories --study lipase_363K --dry-run   # what it would hash
polyzymd hash-trajectories --study lipase_363K             # record them
```

```
No polymer replicate 1: hashed 3, already recorded 0
SBMA 50% replicate 1: hashed 12, already recorded 0
```

Running it again prints `hashed 0` and changes nothing; a hash already in
`progress.json` is never overwritten, and a disagreement is reported as
`conflict:` with exit code 2. `--verify` rehashes and compares. On a cluster,
run it in a batch job: it reads every trajectory once. `-c config.yaml` works
for runs outside a study. See {doc}`../reference/cli_reference`.

## Commit as you go

Every report from `polyzymd analyze --study` records the study file's
SHA-256 and the folder's git commit, with any uncommitted files, under
`provenance.study` in `results/<run>/report.json`. Uncommitted inputs, such
as an edited `study.yaml` or analysis function, give a warning; changes in
`results/` and `data.local.yaml` do not, because those are outputs and
pointers:

```
warning: the study has uncommitted changes (study.yaml); the report records them, and committing them makes it reproducible.
```

Stored results are reused by content, never by commit, so committing never
recomputes anything. PolyzyMD commits only once, in `study init`; what to
commit after that is up to you. `polyzymd study check` prints the commit and
any uncommitted inputs.
