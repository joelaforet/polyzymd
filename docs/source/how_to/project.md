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
project a git repository. Add each protein's conditions, control first:

```bash
polyzymd study add-condition "No Polymer" --config path/to/no_polymer/config.yaml --study Paper_1/lipa363
polyzymd study add-condition "SBMA 50%" --config path/to/sbma50/config.yaml --study Paper_1/lipa363
```

Then give each polymer condition its factor by editing its line in
`study.yaml` (there is no command for it):

```yaml
conditions:
  No Polymer: conditions/no_polymer/config.yaml
  SBMA 50%: {config: conditions/sbma_50/config.yaml, factors: {sbma_fraction: 0.5}}
```

To bring studies you already have into a project, follow the tutorial
{doc}`../tutorials/move_studies_into_project`.

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

Give the control a factor value only when it is a point on the same axis.
For polymer loading (grams of polymer, or polymer chains per protein), the
control is loading 0: write `factors: {polymer_loading: 0}` and the trend
includes it. For the polymer's composition, such as `sbma_fraction`, a
control without polymer has no composition: leave its factor out, and the
trend is fitted over the polymer conditions only.

## Write `project.yaml`

```yaml
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
metadata:                      # title, authors, license, keywords, related paper and data
  title: ...
```

- `region <name>` in a selection becomes that study's region, so `rmsf`
  aligns on each protein's own core. `structure <name>` as a value becomes the
  path of that study's structure.
- A `function:` path is relative to the file that lists it: in `project.yaml`,
  `analyses/lid.py` is `Paper_1/analyses/lid.py`.
- A name a study does not define stops the command with a message naming the
  study. Limit an analysis to the proteins that have the region with
  `studies: [...]`, or put a protein's own analysis in its `study.yaml`.
- Records and reports keep the resolved selections, so each result says which
  residues it measured.

## Check and analyse

```bash
polyzymd project check Paper_1                    # which studies run each analysis, then each study's check
polyzymd analyze rmsf --project Paper_1           # every study that runs rmsf, one after another
polyzymd analyze --project Paper_1                # every analysis of every study
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
paper.replicate_table("rmsf")   # one row per replicate: the values every test uses

# Every study's comparisons with its control, and its trend tests, side by side:
for study, report in rmsf.reports.items():
    for pair in report.pairwise:
        print(study, pair.b, "vs", pair.a, pair.delta, pair.p_adjusted, pair.significant)
    for trend in report.trends:
        print(study, trend.factor, trend.slope, trend.p_adjusted, trend.reason)
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
such as `sbma_fraction`, every report of that study adds a straight line
fitted through the condition means against the factor, over the conditions
that declare it (a control without the factor is left out), with a 95
percent interval and a t test of zero slope; several numeric factors are one
Benjamini-Hochberg family:

```
trend sbma_fraction  slope -0.4  ci95 -0.61 to -0.19  p 0.009  p_adj 0.009  r2 0.92  condition_means 5  replicates 15
verdict: core_rmsf falls with sbma_fraction (slope -0.4 A per unit sbma_fraction, ..., fitted on 5 condition means of 15 replicates)
```

Each condition is one point, the mean of its replicates: the factor varies
only between conditions, so with five conditions the test has three degrees
of freedom, and only a clear, steady change across conditions is called a
trend. It needs at least three factor levels (with two it would restate a
pairwise comparison) and is reported as not testable, with the reason, when
a replicate value is not finite. Labelled results (one value per residue) get
none. A different model, such as a dose response that is not a straight
line, belongs in your statistical plan.

**Your own statistical plan.** A plan that goes further is a function in the project (or study) folder:

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
each study's manifest by SHA-256, every condition as `<study> /
<condition>`, and the size and SHA-256 of every project file outside the
studies (`project.yaml`, the shared `analyses/` and `stats/` code, figures,
the stats plan's output), and one `CITATION.cff` and `.zenodo.json` from
`project.yaml`'s `metadata:` (the keys are those of a study's metadata,
listed in {ref}`study-metadata`). It commits and tags the project
(`project-v1`, ...) and lays out one `deposit/` with `deposit/UPLOAD.md`, as
`polyzymd study freeze` does for a study ({doc}`study_freeze`): one dataset,
one DOI, for the paper.

Before writing anything, freeze warns about every run whose stored results no
longer match the project (a changed config, window, stride, function or
anything in its folder, any key of the analysis entry, setting, selection,
file content, factor, condition, or replicate set),
every partial report, and a stats plan that is stale or was never run. Each
warning names what to rerun; freezing anyway is your call.
