# Move existing studies into a project

This tutorial starts from a blank `polyzymd project init` and moves studies
you already have, each with its own `study.yaml`, into it. You end with one
project folder for the paper, whose analyses are written once and run in every
study, and whose stored results are kept.

The example is a paper about three lipases, each simulated under six polymer
conditions, with an old study folder per enzyme:

```
old/
├── lipa363/   study.yaml  analyses/  results/
├── calb343/   study.yaml  analyses/  results/
└── rml333/    study.yaml  analyses/  results/
```

In this paper, each study is one enzyme at one temperature: `lipa363` is
lipase A at 363 K. That is one layout, not a rule. A study holds every
condition that shares one analysis frame (residue numbering, reference
structures, regions, equilibration window and control). So one study can
vary the temperature and the polymer composition together, with factors such
as `temperature_K` and `sbma_fraction`. See {doc}`../explanation/projects`.

Each old `study.yaml` lists its conditions by absolute config path and names
its reference structure by absolute path. Nothing in `old/` is changed.

:::{admonition} Environment Setup
:class: tip

The commands on this page assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

**You need:** the old study folders, and their simulation configs where the old
`study.yaml` files point. The trajectories need not be on this machine.

## 1. Create the project

One study for each old study folder. A study label is also its folder name,
so use lower case, digits and `_`:

```bash
polyzymd project init Paper_1 --study lipa363 --study calb343 --study rml333
```

```
created project Paper_1 with studies lipa363, calb343, rml333
```

`Paper_1/` now holds `project.yaml`, `analyses/`, `stats/`, `figures/`, the
licenses, and one folder per study, each with a `study.yaml` full of `TODO`s.

## 2. Add the conditions of each study

For each old study, add its conditions in the order the old `study.yaml` lists
them: the first is the control. `--config` copies the config, with the input
files it names, into the new study:

```bash
polyzymd study add-condition "No Polymer" --study Paper_1/lipa363 \
    --config old_runs/lipa_no_polymer/config.yaml
polyzymd study add-condition "SBMA-EGMA 0:100" --study Paper_1/lipa363 \
    --config old_runs/lipa_sbma_egma_0_100/config.yaml
# ... and so on for every condition of every study
```

With many conditions, let Python read the old file and run the commands:

```python
import subprocess
import yaml

for study in ("lipa363", "calb343", "rml333"):
    old = yaml.safe_load(open(f"old/{study}/study.yaml"))
    for label, config in old["conditions"].items():
        subprocess.run(
            ["polyzymd", "study", "add-condition", label,
             "--study", f"Paper_1/{study}", "--config", config],
            check=True,
        )
```

**Check:** `Paper_1/lipa363/conditions/` has one folder per condition, and
`conditions:` in `Paper_1/lipa363/study.yaml` lists them in order.

To record what varies between conditions, for trend tests and plots, give
each polymer condition its factor (the control has none):

```yaml
conditions:
  No Polymer: conditions/no_polymer/config.yaml
  SBMA-EGMA 0:100: {config: conditions/sbma-egma_0_100/config.yaml, factors: {sbma_fraction: 0.0}}
  SBMA-EGMA 25:75: {config: conditions/sbma-egma_25_75/config.yaml, factors: {sbma_fraction: 0.25}}
```

## 3. Describe the study and its window

In each new `study.yaml`, replace the `TODO`s, copying the window from the old
file:

```yaml
description: Bacillus subtilis lipase A (1ISP) at 363 K
equilibration: 200ns
```

Copy `stride:`, `until:` and `replicates:` too if the old file has them.

## 4. Name the structures

Copy each reference structure the old analyses name into the study's
`structures/` folder, and name it:

```bash
cp old_runs/structures/1ISP_clean.pdb Paper_1/lipa363/structures/
```

```yaml
structures:
  reference: structures/1ISP_clean.pdb
```

Use the same name for the equivalent structure in every study (here
`reference` for each crystal structure), so one analysis can say
`structure reference` for all of them. A protein with a second structure,
such as RML's closed conformation, names it too: `closed: structures/3TGL....pdb`.

## 5. Name the regions

Every residue set the old analyses list (a core, the catalytic triad, the
active site) becomes a named region, in the residue numbering of that study:

```yaml
regions:
  core: resid 5-8 15-27 32-37 47-65 70-75 77-88 91-93 95-101 105-107 111 123-129 137-140 143 146-150 156-160 162-175 177
  catalytic_triad: resid 76 132 155
  oxyanion_hole: resid 11 77
  active_site_5A: resid 10 11 74-80 101-105 128-137 153-159
```

Give the same region the same name in every study: `core` is LipA's core in
`lipa363` and RML's in `rml333`.

## 6. Move the analyses

Each old analysis entry goes to one of three places:

- **every study runs it:** `analyses:` in `Paper_1/project.yaml`, written
  once, with `structure <name>` and `region <name>` where the old entry had a
  path or a residue list;
- **only some studies run it:** the same, with `studies: [...]` listing them;
- **only one study runs it** (RML's lid gate): `analyses:` in that study's
  `study.yaml`.

An old entry

```yaml
  rmsf:
    selection: protein and name CA
    alignment_selection: protein and name CA and resid 5 6 7 8 15 16 ...
    reference_mode: external
    reference_file: old_runs/structures/1ISP_clean.pdb
```

becomes, in `project.yaml`,

```yaml
analyses:
  rmsf:
    selection: protein and name CA
    alignment_selection: protein and name CA and region core
    reference_mode: external
    reference_file: structure reference
```

Copy the Python files the entries name: shared functions into
`Paper_1/analyses/`, the functions of one study into `Paper_1/<study>/analyses/`. A
`function:` path is relative to the file that lists it, so
`function: analyses/dssp.py:dssp_core` in `project.yaml` is
`Paper_1/analyses/dssp.py`.

Two things the old results depend on must go into the entries:

- **Options given on the command line** when the results were made, such as
  `--stride 10` or `--eq 0ns`, become `stride: 10` or `equilibration: 0ns` in
  the entry. Otherwise the stored results are stale against the new file.
- **Settings that pointed into the old study folder**, such as a sidecar
  directory, become paths in the new one.

## 7. Bring the results

Copy each old study's results into its new study:

```bash
cp -r old/lipa363/results Paper_1/lipa363/
```

## 8. Point at the trajectories

The copied configs still say where their runs are. If that is right on this
machine, there is nothing to do. Otherwise, for each study:

```bash
polyzymd study locate /path/to/the/runs --study Paper_1/lipa363
```

which writes `data.local.yaml` (never committed or published).

## 9. Check

```bash
polyzymd project check Paper_1
```

**Success looks like:** one `analysis <run>: studies ...` line per project
analysis, naming the studies that run it, then each study's check with its
structures, regions and conditions, and no `error:` line. A name a study does
not define is an error naming the study and the analysis: define it in that
study, or limit the analysis with `studies:`.

Then rerun the analyses:

```bash
polyzymd analyze --project Paper_1
```

A stored result whose record matches is read back, not measured again.
Results stored by a PolyzyMD version that identified input files differently
are measured once more.

## 10. Finish

- Fill `metadata:` in `project.yaml` (title, authors with `family-names`,
  `given-names` and ORCID, licenses).
- Move any standalone statistics script into `stats/`, reading results with
  `pz.Project(".").replicate_table(run)` ({doc}`../how_to/project`).
- Commit: `git -C Paper_1 add -A && git -C Paper_1 commit -m "Move the studies into a project"`.

Read the results of every study in one table:

```python
import polyzymd as pz

pz.Project("Paper_1").results("rmsf").table
```

## Next steps

- {doc}`../how_to/project`: analyse, run statistics and publish a project.
- {doc}`../explanation/projects`: why the conditions of a study share one analysis frame.
