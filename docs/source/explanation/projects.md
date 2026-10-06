# Projects and studies: one paper, one study per protein

PolyzyMD organizes experiments in four levels:

| Level | Is | Holds | File |
|---|---|---|---|
| Project | One paper | The analyses every study runs, statistics and figures code, publishing metadata | `project.yaml` |
| **Study** | **One protein (or other system) under its conditions**, such as protein X | Its equilibration window, its structures, its named residue regions, its conditions | `study.yaml` |
| Condition | One simulated variant of that protein, such as a polymer composition or a different temperature | A simulation `config.yaml`, optionally its `factors:` | `conditions/<name>/config.yaml` |
| Replicate | One run of a condition | Trajectories; its number seeds its starting structure | the run directory |

## A study is tied to its protein

Everything an analysis needs to know about a protein differs from one protein
to the next:

- **Residue numbering.** The catalytic serine of lipase A is residue 76 in the
  trajectory; RML's numbering is offset from its crystal structure. A
  selection written for one protein selects the wrong atoms in another.
- **Reference structures.** Native contacts and RMSF are measured against a
  protein's own crystal structure, and some proteins need two (an open and a
  closed lid).
- **Regions.** Which residues form the core, the active site or a lid is a
  property of the fold.
- **Equilibration.** Each protein relaxes at its own pace, at its own
  temperature.
- **The control.** A polymer condition is meaningful only against the same
  protein without polymer.

So a study holds exactly one protein, and every comparison is made within a
study, against that study's control. Comparing one protein's polymer
condition with another protein's control is never something PolyzyMD does by
default.

## What the project adds

What stays the same across proteins is the question. A project writes each
analysis once, in `project.yaml`, and runs it in every study, so the paper's
three enzymes are measured the same way by construction rather than by
keeping three copies of a file in step. The project also holds what belongs
to the paper as a whole: the statistical plan, the figures, and the
publishing metadata, which is written once and published as one dataset with
one DOI.

A study can still stand on its own. A project is only needed to analyse
several proteins the same way, and every study command works on a study
inside a project as it does on a study alone.

## How one analysis fits every protein

An analysis written once has to select different residues and read different
reference files in each protein. Studies therefore name their protein-specific
parts, and the analysis refers to the names:

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

In each study, `region core` becomes that protein's core selection and
`structure reference` its crystal structure. The resolved selection, not the
name, is what the analysis runs with and what its stored record holds, so a
result always says which residues it measured.

A name a study does not define stops the analysis with an error naming the
study. PolyzyMD never skips a protein silently, because a figure missing one
enzyme looks complete. When an analysis applies only to some proteins (a lid
opening, where one enzyme has no lid), the project says so explicitly with
`studies: [...]`, and `polyzymd project check` lists which studies run each
analysis, so a protein left out is visibly left out. An analysis that only one protein has belongs in that study's
own `study.yaml`.

Results stay with the protein: each study stores its results in its own
`results/` folder. Reading them through the project puts every study's rows
in one table with a `study` column, which is what cross-protein figures need,
while each study's report keeps its own comparisons.

## Statistics follow the same rule

The sampling unit is the replicate, the independent simulation, not the
frame (Grossfield et al. 2018). Every test in PolyzyMD works on one value per
replicate, and the project and its studies give exactly that table.

Conditions may declare factors, such as the SBMA fraction of a copolymer.
Factors describe what varies between conditions; they do not group anything.
A numeric factor lets a study's report test for a trend, in addition to
comparing each condition with the control. The trend is fitted through the
condition means, one point per condition, because the factor varies only
between conditions: more replicates of a condition make its mean more
precise, but they are not more points on the line. Treating them as
independent points would claim a confidence the conditions cannot give
(Hurlbert 1984; Lazic 2010). Factors need not
match across studies: each study tests the factors its conditions declare.

A paper's statistical plan often goes further, for example a hierarchy of
tests with multiplicity control across them, or comparisons across proteins
made deliberately. That plan is a script in the project's `stats/` folder,
reading the stored results with `replicate_table`, so it is committed and
published with the paper instead of living beside it.

## Reproducing and publishing

A project is the unit of publication. Freezing it freezes every study,
records each study's manifest by its hash in one project manifest, and
writes one citation and one deposit, so a reader reproduces the paper from
one dataset.

Stored results identify their inputs by content: trajectories, topologies
and argument files by name, size and SHA-256, and configs by a hash that
leaves out where data lives. A study moved into a project, or a project
moved to another machine, therefore reuses every stored result whose inputs
are unchanged.

## See also

- {doc}`../how_to/project`: analyse, test and publish a project.
- {doc}`../how_to/move_studies_into_project`: bring existing studies into a project.
- {doc}`study_folders`: what a study folder holds and how it is published.
