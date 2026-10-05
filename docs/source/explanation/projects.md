# Projects and studies: one paper, one study per protein

```{admonition} Status: design, being implemented
:class: warning
This page is the agreed design for grouping studies into a project. It is
built in three slices: P1 (the project file, shared analyses, regions and
structures, `Project.results`), P2 (statistics: `replicate_table`, trend
tests and the `stats:` hook) and P3 (`project init`, `project freeze` and
moving an existing study into a project). Each section names its slice.
P1, P2 and P3 are implemented ({doc}`../how_to/project`).
```

## The hierarchy

| Level | Is | Holds | File |
|---|---|---|---|
| Project | One paper or thesis chapter | The analyses every study runs, the statistical plan, publishing metadata, figures | `project.yaml` |
| **Study** | **One protein (or other system) under one set of conditions**, such as lipase A at 363 K | Its equilibration window, its structures, its named residue regions, its conditions | `study.yaml` |
| Condition | One simulated variant of that protein, such as a polymer composition | A simulation `config.yaml`, optionally `factors:` | `conditions/<name>/config.yaml` |
| Replicate | One run of a condition | Trajectories; its number is its random seed | the run directory |

**A study is tied to its protein.** Everything that depends on the protein
lives in its study: the residue numbering, the reference structure, which
residues form the core, the active site or a lid, and how long the protein
takes to equilibrate. Comparisons are made within a study, against the
study's own control. Comparing one protein's polymer condition with another
protein's control is never a default.

The project holds what the paper asks of every protein: the analyses, with
their settings written once, and the statistics. A study can be used on its
own; a project is needed only to analyse several proteins the same way.

## A project folder (P1)

```
Paper_1/
├── project.yaml
├── analyses/          # functions shared by every study
├── stats/             # the statistical plan (P2)
├── figures/           # read results through pz.Project
├── lipa363/           # one study per protein
│   ├── study.yaml
│   ├── structures/
│   ├── conditions/
│   └── results/       # this protein's stored results
├── calb343/
└── rml333/
```

`project.yaml` lists the studies and the analyses they all run:

```yaml
polyzymd: 1.3.0
studies:                  # label: folder holding its study.yaml
  lipa363: lipa363
  rml333: rml333
analyses:
  native_contacts:
    selection: protein and not element H
    reference_file: structure reference   # each study's structures: reference
  lid_opening:
    function: analyses/lid.py:lid_distance
    kind: timeseries
    studies: [rml333]                      # only the proteins that have a lid
    selections: {lid: region lid, core: region core}
metadata: {title: ..., authors: [...], license: {...}}
```

Each `study.yaml` holds what is specific to its protein:

```yaml
description: Bacillus subtilis lipase A (1ISP) at 363 K
equilibration: 200ns
structures:
  reference: structures/1ISP_clean.pdb
regions:
  core: resid 5-8 15-27 32-37
  catalytic_triad: resid 76 132 155
conditions:               # control first
  No Polymer: conditions/no_polymer
  SBMA-EGMA 50:50: {config: conditions/sbma_50, factors: {sbma_fraction: 0.5}}
analyses: {}              # analyses only this protein runs
```

Rules:

- **Names in selections.** `region <name>` in any selection string becomes
  that study's region, in parentheses, and a value `structure <name>` becomes
  the path of that study's structure. A name the study does not define is an
  error naming the study, never a silent skip. Records and reports keep the
  resolved selection, so a result says which residues it measured.
- **Which studies run an analysis.** A project analysis runs in every study
  unless it lists `studies:`. An analysis that only one protein has goes in
  that study's own `analyses:`. A study may not redefine a project analysis
  under the same name.
- **Results stay with the protein.** `polyzymd analyze RUN --project Paper_1`
  runs RUN in every study that runs it, and stores each study's results in
  that study's `results/RUN`, exactly as `polyzymd analyze RUN --study
  Paper_1/lipa363` would. A study inside a project can still be analysed on
  its own.
- **One table for the paper.** `pz.Project("Paper_1").results(RUN).table` is
  every study's table with a `study` column, and a column for each factor.
  Each study's report stays separate.

## Statistics (P2)

- `study.replicate_table(RUN)` and `project.replicate_table(RUN)` give one
  row per replicate: study, condition, factors, replicate and its value after
  the analysis window. This is the sampling unit for every test.
- When a study's conditions declare a numeric factor, the report adds a trend
  test: the slope of the replicate values against that factor, over the
  conditions that declare it (a control without the factor is left out).
- `stats: {plan: stats/plan.py:plan}` in `project.yaml` or `study.yaml` names
  a function that receives those tables and returns its results. They are
  stored with the reports, keyed on the function's source like any analysis,
  and deposited with the project, so the statistical plan ships with the
  paper.

## Publishing (P3)

- `polyzymd project init NAME --study A --study B` writes the folder above.
- `polyzymd project freeze` checks every study as `study freeze` does and
  writes one manifest, `CITATION.cff` and deposit for the project, so the
  paper has one DOI. Publishing metadata lives in `project.yaml`, once.
- `polyzymd project init --from` moves existing studies into a project:
  their configs into `conditions/`, machine paths into `data.local.yaml`,
  and their results kept. Stored results recompute once on the next
  `analyze` when the PolyzyMD version that wrote them identified files
  differently.
