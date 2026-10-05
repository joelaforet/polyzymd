# Analyse several proteins as one project

Use this when a paper studies more than one protein, such as three lipases,
each simulated under the same polymer conditions, and every protein must be
analysed the same way. For the design, see {doc}`../explanation/projects`.

:::{admonition} Environment Setup
:class: tip

The commands on this page assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## One study per protein

**A study is tied to its protein.** Each protein gets its own study folder,
with its own `study.yaml`, because everything below depends on the protein:

- its residue numbering, and so the selections of its core, active site or lid;
- its reference structure;
- how long it takes to equilibrate;
- its control: each study's conditions are compared with that study's first
  condition, never with another protein's.

The project holds what the paper asks of every protein: the analyses, with
their settings written once.

```
Paper_1/
├── project.yaml
├── analyses/           # functions every study uses
├── lipa363/study.yaml  # lipase A at 363 K
├── calb343/study.yaml
└── rml333/study.yaml
```

## Start a project

```bash
polyzymd project init Paper_1 --study lipa363 --study calb343 --study rml333
```

writes `project.yaml`, `analyses/`, `stats/`, `figures/`, licences and one
study folder per protein, each with a `study.yaml` to fill in, and makes the
project a git repository. Add each protein's conditions with
`polyzymd study add-condition Paper_1/lipa363 ...`.

### Move existing studies in

Studies made before projects existed, one `study.yaml` per protein, move in
with `LABEL=path`:

```bash
polyzymd project init Paper_1 \
    --study lipa363=Paper_1_REDO/lipa363 \
    --study calb343=Paper_1_REDO/calb343 \
    --study rml333=Paper_1_REDO/rml333
```

Each old study is only read. Its conditions' configs and input files are
copied into `conditions/`, where its runs are today goes into the study's
`data.local.yaml`, and its `analyses/` code and `results/` are copied. A
setting that names a file, such as a `reference_file`, is copied into the
study's `structures/` and written as `structure reference` (or the file's
name, when the study names several). Analyses that every study defines alike
move into `project.yaml`. The command lists settings it left as absolute
paths. Stored results stay valid: the next `polyzymd analyze --project`
reuses them, unless they were written by a PolyzyMD version that identified
input files differently, in which case they are measured again once.

Then name each protein's regions and replace residue lists that differ only
by numbering with `region <name>`, so the analysis can move into
`project.yaml`.

## Write each protein's `study.yaml`

Name the protein's structures and residue regions, so the project's analyses
can use them:

```yaml
description: Bacillus subtilis lipase A (1ISP) at 363 K
equilibration: 200ns
structures:
  reference: structures/1ISP_clean.pdb
regions:                       # this protein's numbering
  core: resid 5-8 15-27 32-37
  catalytic_triad: resid 76 132 155
conditions:                    # control first
  No Polymer: conditions/no_polymer
  SBMA-EGMA 50:50: {config: conditions/sbma_50, factors: {sbma_fraction: 0.5}}
analyses: {}                   # analyses only this protein runs
```

A condition is the folder holding its `config.yaml` (or the file itself), or
`{config: ..., factors: {...}}` when you want to record what varies between
conditions; factors become columns of the results table.

## Write `project.yaml`

```yaml
polyzymd: 1.3.0
studies:                       # label: folder
  lipa363: lipa363
  calb343: calb343
  rml333: rml333
analyses:
  rmsf:
    selection: protein and name CA
    alignment_selection: protein and name CA and region core
    reference_file: structure reference
  lid_opening:                 # only the proteins that have a lid
    function: analyses/lid.py:lid_distance
    kind: timeseries
    studies: [calb343, rml333]
    selections: {lid: region lid, core: region core}
metadata:
  title: ...
```

- `region <name>` in a selection becomes that study's region, so `rmsf`
  aligns on each protein's own core. `structure <name>` as a value becomes the
  path of that study's structure.
- A name a study does not define stops the command with a message naming the
  study. Limit an analysis to the proteins that have the region with
  `studies: [...]`, or put a protein's own analysis in its `study.yaml`.
- Records and reports keep the resolved selections, so each result says which
  residues it measured.

## Check and analyse

```bash
polyzymd project check Paper_1                    # which studies run each analysis, then each study's check
polyzymd analyze rmsf --project Paper_1           # every study that runs rmsf, one after another
polyzymd analyze rmsf --study Paper_1/lipa363     # one protein
```

Each study's results go to its own `results/<run>`, so a study inside a
project can also be analysed and read on its own.

## Read every protein's results

```python
import polyzymd as pz

paper = pz.Project("Paper_1")
rmsf = paper.results("rmsf")
rmsf.table          # every study's rows, with a study column and a column per factor
rmsf.reports        # each study's report: comparisons stay within a protein
```

A study that runs the analysis but has no stored results is named in the
error, so no protein is left out of a figure silently.

## Statistics

**One row per replicate.** The replicate is the sampling unit of every test
(Grossfield et al. 2018). `replicate_table` gives exactly that, from stored
results, with a column per factor:

```python
paper.replicate_table("rmsf")          # every study; or pz.Study(...).replicate_table("rmsf")
```

Per-frame results are reduced over each replicate's frames as the analysis
reduces them (the entry's `reduce`, otherwise the mean).

**Trend over a factor.** When a study's conditions declare a numeric factor,
such as `sbma_fraction`, every report of that study adds the slope of the
replicate values against it, over the conditions that declare it (a control
without the factor is left out), with a 95 percent interval and a t test of
zero slope; several numeric factors are one Benjamini-Hochberg family:

```
trend sbma_fraction  slope 2  ci95 1.977 to 2.023  p 1.7e-15  p_adj 1.9e-14  r2 0.99  n 15  conditions 5
verdict: core_rmsf falls with sbma_fraction (slope ...)
```

It is a straight-line fit: look at the per-condition means before reading it
as a dose response. Labelled results (one value per residue) get none.

**Your own statistical plan.** A plan that goes further, such as Paper 1's
gatekeeping (Welch against the control, then a trend, then bootstrap
intervals), is a function in the project (or study) folder:

```yaml
stats:
  plan: stats/plan.py:plan
```

```python
# stats/plan.py
def plan(project):
    table = project.replicate_table("core_rmsf")
    ...
    return {"tier1": tier1_table, "tier2": tier2_table, "alpha": 0.05}
```

```bash
polyzymd stats Paper_1
```

It receives the `Project` (or `Study`) and returns a dict of tables and
values, written to `results/stats/<function>/` (`<name>.csv` per table,
`values.json` for the rest) with `record.json`: the SHA-256 of the plan's
file and of every report it could read. `polyzymd project check` then says
whether it is up to date, or stale because its code or the analyses changed.
The plan lives in the folder, so it is committed and published with the
paper.

## Publish

```bash
polyzymd project freeze Paper_1
```

freezes every study (its manifest, checklist, system summary, engine inputs
and final frames), then writes the project's `manifest.json`, which lists
each study's manifest by SHA-256 and every condition as `<study> /
<condition>`, and one `CITATION.cff` and `.zenodo.json` from
`project.yaml`'s `metadata:`. It commits and tags the project
(`project-v1`, ...) and lays out one `deposit/` with `deposit/UPLOAD.md`, as
`polyzymd study freeze` does for a study ({doc}`study_freeze`): one dataset,
one DOI, for the paper.
