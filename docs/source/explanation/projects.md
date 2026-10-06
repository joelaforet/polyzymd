# Projects and studies: one paper, one study per protein

PolyzyMD organizes simulations in four levels:

| Level | What it is | What it holds | File |
|---|---|---|---|
| {term}`Project <project>` | One paper | The analyses that every study runs, the statistics and figure code, the publishing metadata | `project.yaml` |
| **{term}`Study <study>`** | **One protein (or other system) under its conditions** | Its equilibration window, its structures, its named residue regions, its conditions | `study.yaml` |
| {term}`Condition <condition>` | One simulated variant of that protein, such as a polymer composition, a co-solvent or a temperature | A simulation `config.yaml`, and optionally its `factors:` | `conditions/<name>/config.yaml` |
| {term}`Replicate <replicate>` | One independent simulation of a condition | Its trajectories. Its number seeds its starting structure and its dynamics | The {term}`replicate folder` |

## A study is tied to its protein

Each protein differs in everything that an analysis must know about it:

- **Residue numbering.** The catalytic serine of lipase A is residue 76 in the
  trajectory. The numbering of RML is offset from its crystal structure. A
  selection written for one protein selects the wrong atoms in another.
- **Reference structures.** Native contacts and RMSF use the crystal
  structure of the protein. Some proteins need two, such as an open lid and a
  closed lid.
- **Regions.** The fold of the protein sets which residues form the core, the
  active site or a lid.
- **Equilibration.** Each protein relaxes at its own rate, at its own
  temperature.
- **The control.** A polymer condition has meaning only against the same
  protein without polymer.

So a study holds exactly one protein. PolyzyMD makes every comparison within
a study, against the control of that study. By default, PolyzyMD never
compares a polymer condition of one protein with the control of a different
protein.

## What the project adds

The question stays the same across the proteins of a paper. A project writes
each analysis once, in `project.yaml`, and runs it in every study. So the
three enzymes of a paper get the same measurement by design. You do not keep
three copies of one file in step.

The project also holds what belongs to the paper as a whole:

- the statistical plan;
- the figures;
- the publishing metadata, published as one dataset with one DOI.

A study can still stand on its own. You need a project only to analyze
several proteins the same way. Every study command works the same on a study
inside a project and on a study alone.

## How one analysis fits every protein

One analysis must select different residues and read different reference
files in each protein. So each study names its protein-specific parts, and
the analysis uses the names:

```yaml
# project.yaml
analyses:
  rmsf:
    alignment_selection: protein and name CA and region core
    reference_file: structure reference
```

```yaml
# lipa363/study.yaml
structures: {reference: structures/1ISP_clean.pdb}
regions: {core: resid 5-8 15-27 32-37}
```

In each study, `region core` becomes the core selection of that protein.
`structure reference` becomes its crystal structure. The analysis runs with
the resolved selection, not with the name. The stored record holds the
resolved selection too, so each result says which residues it measured.

If a study does not define a name, the analysis stops with an error that
names the study. PolyzyMD never skips a protein without a message, because a
figure that is missing one enzyme looks complete.

Some analyses apply only to some proteins, such as a lid opening when one
enzyme has no lid. The project then lists the proteins with
`studies: [...]`. `polyzymd project check` lists which studies run each
analysis, so you can see each protein that is left out. Put an analysis that
only one protein has in the `study.yaml` of that protein.

Results stay with the protein. Each study stores its results in its own
`results/` folder. When you read them through the project, the rows of every
study go into one table with a `study` column. Cross-protein figures need
this table. The report of each study keeps its own comparisons.

(project-statistics)=
## Statistics follow the same rule

The sampling unit is the replicate, the independent simulation. It is not the
frame (Grossfield et al. 2018). Every test in PolyzyMD uses one value per
replicate. `replicate_table` gives that table, for a project or for a study.
It reduces the per-frame values of each replicate the way the analysis does:
with the `reduce` of the analysis entry, or else with the mean.

Conditions can declare factors, such as the SBMA fraction of a copolymer.
Factors describe what varies between conditions. They do not group anything.
Factors need not match across studies: each study tests the factors that its
conditions declare.

### Why a trend uses one point per condition

A numeric factor lets the report of a study test for a trend, in addition to
the comparison of each condition with the control. PolyzyMD fits the line
through the condition means, with one point per condition.

The factor varies only between conditions. More replicates make the mean of
a condition more precise, but they do not add points to the line. If the
replicates were independent points, the test would claim a confidence that
the conditions cannot give (Hurlbert 1984; Lazic 2010).

With one point per condition, the test has few degrees of freedom: n - 2 for
n condition means. So the test finds only a clear, steady change across the
conditions. The test needs at least three factor levels. With two levels it
would repeat a pairwise comparison.

### The statistical plan of the paper

The plan of a paper often goes further than the report. Examples are a
hierarchy of tests with multiplicity control across them, or comparisons
across proteins that you choose on purpose. Write that plan as a script in
the `stats/` folder of the project. The script reads the stored results with
`replicate_table`. So the plan is committed and published with the paper.

## Reproducing and publishing

A project is the unit of publication. When you freeze a project, PolyzyMD
does these steps:

1. It freezes every study.
2. It records the hash of the manifest of each study in one project manifest.
3. It writes one citation and one deposit.

A reader then reproduces the paper from one dataset.

Stored results identify their inputs by content:

- trajectories, topologies and argument files by name, size and SHA-256;
- configs by the {term}`config hash`, which leaves out where the data is.

So when you move a study into a project, or move a project to a different
machine, every stored result with unchanged inputs is reused.

## See also

- {doc}`../how_to/project`: analyze, test and publish a project.
- {doc}`../how_to/move_studies_into_project`: bring existing studies into a project.
- {doc}`study_folders`: what a study folder holds and how it is published.
