# Projects and studies: one paper, one analysis frame per study

PolyzyMD organizes simulations in four levels:

| Level | What it is | What it holds | File |
|---|---|---|---|
| {term}`Project <project>` | The studies of one paper | The analyses that the studies share, the `stats/` scripts, the figure code, the publishing metadata | `project.yaml` |
| **{term}`Study <study>`** | **A set of conditions that you compare with each other** | Its analysis frame: its equilibration window, its structures, its named residue regions, its control. Its conditions | `study.yaml` |
| {term}`Condition <condition>` | One point in the space of independent variables: one value of each variable, such as the protein variant, the polymer composition, the co-solvent and the temperature | A simulation `config.yaml`, and optionally its `factors:`, which name its coordinates | `conditions/<name>/config.yaml` |
| {term}`Replicate <replicate>` | One independent simulation of a condition | Its trajectories. Its number seeds its starting structure and its dynamics | The {term}`replicate folder` |

## A study shares one analysis frame

A study is a set of conditions that you compare with each other: a part of
the space of independent variables. The study holds fixed each variable that
it does not vary. All its conditions share one analysis frame: the same atom
and residue numbering, so one selection means the same atoms in every
condition; the same reference structures; the same named regions; one
equilibration window; and a control to compare against. PolyzyMD compares
conditions within one study, never across studies. Conditions that cannot
share one frame go in separate studies. The most common reason is a different
protein.

Each part of the frame is something that an analysis must know:

- **Residue numbering.** The catalytic serine of lipase A is residue 76 in the
  trajectory. The numbering of RML is offset from its crystal structure. A
  selection written for one protein selects the wrong atoms in another.
- **Reference structures.** Native contacts and RMSF use one crystal
  structure for every condition of the study. Some proteins need two, such as
  an open lid and a closed lid.
- **Regions.** The fold of the protein sets which residues form the core, the
  active site or a lid. A region names the same residues in every condition.
- **Equilibration.** PolyzyMD removes the same start of production from
  every condition of the study.
- **The control.** PolyzyMD compares each condition with the first condition
  of the study, or, with a `comparison:` block, with the control of its own
  stratum. A polymer condition has meaning only against the same system
  without polymer.

So two conditions go in one study when they share all of these parts. A
study can vary several variables together. For example, one enzyme at 300 K,
330 K and 360 K, each with and without SBMA, is one study with six
conditions:

```yaml
conditions:                    # control first
  No polymer 300 K: {config: conditions/none_300, factors: {temperature_K: 300}}
  No polymer 330 K: {config: conditions/none_330, factors: {temperature_K: 330}}
  No polymer 360 K: {config: conditions/none_360, factors: {temperature_K: 360}}
  SBMA 300 K: {config: conditions/sbma_300, factors: {temperature_K: 300, polymer: SBMA}}
  SBMA 330 K: {config: conditions/sbma_330, factors: {temperature_K: 330, polymer: SBMA}}
  SBMA 360 K: {config: conditions/sbma_360, factors: {temperature_K: 360, polymer: SBMA}}
```

By default, PolyzyMD compares each condition with the first condition, the
control. In this study, `SBMA 360 K` is then compared with `No polymer 300 K`,
so that difference holds the effect of the polymer and the effect of the
temperature.

To compare each polymer condition with the no-polymer condition at the same
temperature, add a `comparison:` block:

```yaml
comparison:
  within: temperature_K        # one factor name, or a list
```

The conditions with the same `temperature_K` form one stratum. Each condition
is compared with the control of its stratum: the condition whose other
factors equal those of the first condition. Here the first condition has no
factor other than `temperature_K`, so each control is the condition with no
other factor: `No polymer 300 K`, `No polymer 330 K` and `No polymer 360 K`.
`SBMA 360 K` is compared with `No polymer 360 K`, so the difference holds the
effect of the polymer only. `control: {polymer: none}` names the control's
factor values instead, when every condition declares `polymer`. PolyzyMD
refuses a stratum with no control or with two.

The Benjamini-Hochberg correction covers every comparison of the report,
over all strata. The trend tests are unchanged. To fit both factors in one
model, write a script in `stats/`.

A point mutant keeps the residue numbering of its wild type. It can share the
frame of the wild type when one reference structure fits both. A different
protein usually cannot: its numbering, structures and regions differ. So it
goes in its own study. By default, PolyzyMD never compares a condition of one study with the
control of a different study.

The examples on these pages come from a paper with one enzyme at one
temperature in each study: `lipa363` is lipase A at 363 K. That layout is one
choice, not a rule.

## What the project adds

The question stays the same across the studies of a paper. A project writes
each analysis once, in `project.yaml`, and runs it in every study. So the
three enzymes of a paper get the same measurement by design. You do not keep
three copies of one file in step.

The project also holds what belongs to the paper as a whole:

- the statistical plan;
- the figures;
- the publishing metadata, published as one dataset with one DOI.

A study can still stand on its own. You need a project only to analyze
several studies the same way. Every study command works the same on a study
inside a project and on a study alone.

## How one analysis fits every study

One analysis must select different residues and read different reference
files in each study. So each study names the parts of its frame, and the
analysis uses the names:

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

In each study, `region core` becomes the core selection of that study.
`structure reference` becomes its crystal structure. The analysis runs with
the resolved selection, not with the name. The stored record holds the
resolved selection too, so each result says which residues it measured.

If a study does not define a name, the analysis stops with an error that
names the study. PolyzyMD never skips a study without a message, because a
figure that is missing one study looks complete.

Some analyses apply only to some studies, such as a lid opening when one
enzyme has no lid. The project then lists the studies with
`studies: [...]`. `polyzymd project check` lists which studies run each
analysis, so you can see each study that is left out. Put an analysis that
only one study runs in the `study.yaml` of that study.

Results stay with the study. Each study stores its results in its own
`results/` folder. When you read them through the project, the rows of every
study go into one table with a `study` column. Figures across studies need
this table. The report of each study keeps its own comparisons.

(project-statistics)=
## Statistics follow the same rule

The sampling unit is the replicate, the independent simulation. It is not the
frame (Grossfield et al. 2018). Every test in PolyzyMD uses one value per
replicate. `replicate_table` gives that table, for a project or for a study.
It reduces the per-frame values of each replicate the way the analysis does:
with the `reduce` of the analysis entry, or else with the mean.

Conditions can declare factors, such as the SBMA fraction of a copolymer or
the temperature. The factors of a condition are its coordinates in the space
of independent variables. They describe what varies between conditions. They
do not group anything. Factors need not match across studies: each study
tests the factors that its conditions declare.

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

PolyzyMD fits each numeric factor on its own, through every condition that
declares it. In a study that varies the temperature and the polymer together,
the temperature line goes through the condition means of every polymer. To
separate the two effects, fit a model with both factors in a `stats/` script.

### The statistical plan of the paper

The plan of a paper often goes further than the report. Examples are a
hierarchy of tests with multiplicity control across them, a model with two
factors together, or comparisons across studies that you choose on purpose.
Write that plan as a script in the `stats/` folder of the project. The script
reads the stored results with `replicate_table`. So the plan is committed and published with the paper.

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
