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
