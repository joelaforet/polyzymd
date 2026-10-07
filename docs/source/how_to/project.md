# Analyze several studies as one project

A {term}`project` holds the studies of one paper. It runs the same analyses
with the same settings in every {term}`study`. Use a project when a paper has
two or more studies, such as three lipases under the same polymer conditions,
one study for each lipase. For the reasons behind
this design, see {doc}`../explanation/projects`.

:::{admonition} Environment Setup
:class: tip

The commands on this page assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## Checklist

Do these steps in this order. The sections below give the details.

1. Make the project: `polyzymd project init Paper_1 --study lipa363 --study calb343`.
   Add a study later with `polyzymd project add-study rml333 --project Paper_1`.
2. Add the conditions of each study, control first, with `polyzymd study add-condition`.
3. Edit each `study.yaml`: add the structures, the regions and the factors.
4. Check the project: `polyzymd project check Paper_1`. To see how long each
   condition ran, also run `polyzymd study check Paper_1/lipa363 --production`.
5. Commit the project: `git -C Paper_1 add -A && git -C Paper_1 commit -m "Inputs"`.
6. Run the analyses: `polyzymd analyze --project Paper_1`.
7. Write your own statistics scripts in `Paper_1/stats/`, and run them.
8. Fill in `metadata:` in `project.yaml`.
9. Commit again.
10. Freeze the project: `polyzymd project freeze Paper_1`.

Commit before you freeze. `polyzymd project freeze` stops when an input file
is not committed, and it names the files.

## Make the project

A study is a set of conditions that you compare with each other. All its
conditions share one analysis frame:

- the residue numbering, and so the selections of the core, the active site
  or the lid;
- the reference structure;
- the equilibration window;
- the control. PolyzyMD compares the conditions of a study with the first
  condition of that study. It never compares them with the control of a
  different study.

Conditions that cannot share one frame go in separate studies, each with its
own folder and `study.yaml`. A different protein is the most common case. A
study can vary several variables together, such as the temperature and the
polymer composition.

The project holds what the paper asks of every study: the analyses and their
settings, written once.

In this example paper, each study is one enzyme at one temperature:

```
Paper_1/
├── project.yaml
├── analyses/           # functions that every study uses
├── stats/              # your statistics scripts
├── figures/            # your figure scripts and notebooks
├── lipa363/study.yaml  # lipase A at 363 K
├── calb343/study.yaml  # CALB at 343 K
└── rml333/study.yaml   # RML at 333 K
```

Make the project with one `--study` for each study:

```bash
polyzymd project init Paper_1 --study lipa363 --study calb343 --study rml333
```

The command writes these files and folders:

- `project.yaml`, which lists the studies;
- `analyses/`, `stats/` and `figures/`, each with a `README.md`;
- `LICENSE-code` (MIT), `LICENSE-data` (CC-BY-4.0), `README.md` and
  `.gitignore`;
- one study folder for each `--study`, with a `study.yaml` to fill in.

It then makes the project a git repository and commits these files. Use
`--no-git` to skip git, and `--holder NAME` to set the copyright holder of
the licenses.

To add a study later, give its label:

```bash
polyzymd project add-study rml333 --project Paper_1
```

`add-study` writes the study folder `Paper_1/rml333/` with a `study.yaml` to
fill in, and adds one line under `studies:` in `project.yaml`. It commits
nothing. The label is also the folder name, so use lower case letters,
digits and `_`.

To bring studies that you already have into a project, see
{doc}`move_studies_into_project`.

## Add the conditions

Add the conditions of each study. Add the control first. For a condition
that you have not simulated yet, write a template config, fill it in, and
copy it for the next condition:

```bash
polyzymd study add-condition "No polymer" --new --study Paper_1/lipa363
polyzymd study add-condition "SBMA 50%" --from "No polymer" --study Paper_1/lipa363
```

The tutorial {doc}`../tutorials/own_system` shows how to fill in the
template. For a condition that you already simulated, copy its config:

```bash
polyzymd study add-condition "No polymer" --config path/to/no_polymer/config.yaml --study Paper_1/lipa363
```

For each condition, `add-condition` does these steps:

1. It writes or copies the config to `conditions/<name>/config.yaml`.
2. It copies each input file that the config names into
   `conditions/<name>/structures/`.
3. It points the config at `Paper_1/runs/<study>/<name>/` for new runs.
4. With `--config`, if the config's `scratch_directory` holds its runs, it
   writes that folder into the study's `data.local.yaml`. Git ignores this
   file.
5. It adds one line under `conditions:` in `study.yaml`.

:::{warning}
The runs go into the project folder, in `runs/`, unless you set
`scratch_directory` in the config. Trajectories can use a lot of disk space.
On a cluster, set `scratch_directory` to scratch storage. Git ignores
`runs/`, and freeze never publishes it.
:::

## Write the `study.yaml` of each study

Name the structures and the residue regions of the study. The analyses of
the project use these names.

```yaml
description: Bacillus subtilis lipase A (1ISP) at 363 K
equilibration: 200ns
structures:
  reference: structures/1ISP_clean.pdb
regions:                       # this study's numbering
  core: resid 5-8 15-27 32-37
  catalytic_triad: resid 76 132 155
conditions:                    # control first
  No polymer: conditions/no_polymer
  SBMA 50%: {config: conditions/sbma_50, factors: {sbma_fraction: 0.5}}
analyses: {}                   # analyses that only this study runs
```

A condition is one of these:

- the folder that holds its `config.yaml`;
- the `config.yaml` file itself;
- `{config: ..., factors: {...}}`, to record what varies between conditions.

Each factor becomes a column of the results table. No command sets factors.
Edit the line of the condition in `study.yaml`.

### Factors without polymers

A factor can be any setting that varies between conditions. This study of
one protein in water varies the temperature:

```yaml
conditions:                    # control first
  300 K: {config: conditions/t300, factors: {temperature_K: 300}}
  330 K: {config: conditions/t330, factors: {temperature_K: 330}}
  360 K: {config: conditions/t360, factors: {temperature_K: 360}}
```

This study varies the concentration of a co-solvent:

```yaml
conditions:
  Water: {config: conditions/water, factors: {urea_M: 0}}
  Urea 2 M: {config: conditions/urea_2m, factors: {urea_M: 2}}
  Urea 4 M: {config: conditions/urea_4m, factors: {urea_M: 4}}
```

### Factors of two variables

A study can vary two variables together. This study varies the temperature
and the polymer. Each condition names both of its coordinates:

```yaml
conditions:                    # control first
  No polymer 300 K: {config: conditions/none_300, factors: {temperature_K: 300}}
  No polymer 330 K: {config: conditions/none_330, factors: {temperature_K: 330}}
  No polymer 360 K: {config: conditions/none_360, factors: {temperature_K: 360}}
  SBMA 300 K: {config: conditions/sbma_300, factors: {temperature_K: 300, polymer: SBMA}}
  SBMA 330 K: {config: conditions/sbma_330, factors: {temperature_K: 330, polymer: SBMA}}
  SBMA 360 K: {config: conditions/sbma_360, factors: {temperature_K: 360, polymer: SBMA}}
```

The report compares each condition with the control, `No polymer 300 K`. The
trend test of `temperature_K` goes through all six condition means. `polymer`
is text, so it gets no trend test; it is a column of the results table. To
fit both factors in one model, write a script in `stats/`.

### Give the control a factor value only when it is on the same axis

- For polymer loading (grams of polymer, or chains per protein), the
  control has loading 0. Write `factors: {polymer_loading: 0}`. The trend
  then includes the control.
- For polymer composition, such as `sbma_fraction`, a control without polymer
  has no composition. Do not give it the factor. The trend then uses the
  polymer conditions only.

## Write `project.yaml`

```yaml
studies:                       # label: folder directly inside the project
  lipa363: lipa363
  calb343: calb343
  rml333: rml333
analyses:
  rmsf:
    selection: protein and name CA
    alignment_selection: protein and name CA and region core
    reference_file: structure reference
  lid_opening:                 # only the studies whose protein has a lid
    function: analyses/lid.py:lid_distance
    kind: timeseries
    studies: [calb343, rml333]
    selections: {lid: region lid, core: region core}
metadata:                      # title, authors, license, keywords, related paper and data
  title: ...
```

The function of `lid_opening` has one keyword argument for each key of
`selections:`. At each frame, PolyzyMD passes each selection as an
MDAnalysis `AtomGroup` of the replicate:

```python
# Paper_1/analyses/lid.py
import numpy as np


def lid_distance(lid, core):
    """Distance between the centers of the lid and the core, in Å."""
    return float(np.linalg.norm(lid.center_of_geometry() - core.center_of_geometry()))
```

These rules apply to the names in `project.yaml`:

- In a selection, `region <name>` becomes the region of that name in each
  study. So `rmsf` aligns on the core of each study.
- As a value, `structure <name>` becomes the path of that structure in each
  study.
- A `function:` path is relative to the file that lists it. In
  `project.yaml`, `analyses/lid.py` is `Paper_1/analyses/lid.py`.
- If a study does not define a name, the command stops. The message names the
  study.
- To run an analysis in some studies only, list them in `studies: [...]`.
  Put an analysis of one study only in the `study.yaml` of that study.
- Records and reports keep the resolved selections. Each result says which
  residues it measured.

## Check and analyze

```bash
polyzymd project check Paper_1                    # which studies run each analysis, then each study's check
polyzymd study check Paper_1/lipa363 --production # the production length of each condition of one study
polyzymd analyze rmsf --project Paper_1           # every study that runs rmsf, one after another
polyzymd analyze --project Paper_1                # every analysis of every study
polyzymd analyze rmsf --study Paper_1/lipa363     # one study
```

`project check` reads no trajectory. `study check --production` reads the
trajectory headers, and gives the simulated time of each condition. Use it to
choose the equilibration window.

`analyze --project` writes one log to `Paper_1/logs/`. The results of each
study go to `results/<run>/` in the folder of that study. You can analyze and
read a study of a project on its own.

## Read the results of every study

```python
import polyzymd as pz

paper = pz.Project("Paper_1")
rmsf = paper.results("rmsf")
rmsf.table          # every study's stored values, with a study column and a column per factor
rmsf.reports        # each study's report; comparisons stay within one study
```

`results(run).table` has one row for each stored value. For a `timeseries`
run, such as `lid_opening`, that is one row for each frame of each replicate.
To get one value for each replicate, use `replicate_table`:

```python
table = paper.replicate_table("lid_opening")
table.groupby(["study", "condition"])["value"].mean()   # the mean of each condition
```

`replicate_table` reduces the per-frame values of each replicate the way the
analysis does: with the `reduce` of the analysis entry, or else with the mean.
For one study, use `pz.Study("Paper_1/lipa363").replicate_table("lid_opening")`.

`replicate_table` has the columns `study`, `condition`, `replicate`, `name`,
`part`, `label`, `value` and `unit`, then one column for each factor. `name`
is the stored quantity. A run that stores two quantities, such as two
hydrogen-bond summaries, has one row for each quantity. In that case, group
by `name` too.

Print the comparisons with the control and the trend tests of each study:

```python
for study, report in rmsf.reports.items():
    for pair in report.pairwise:
        print(study, pair.b, "vs", pair.a, pair.delta, pair.p_adjusted, pair.significant)
    for trend in report.trends:
        print(study, trend.factor, trend.slope, trend.p_adjusted, trend.reason)
```

If a study runs the analysis but has no stored results, the call stops and
names the study. A figure therefore never leaves out a study without a
message.

## Trend tests

When the conditions of a study declare a numeric factor, each report of that
study adds a trend test. PolyzyMD fits a straight line through the condition
means against the factor. It then tests whether the slope is zero. For the
reasons behind this test, see {ref}`project-statistics`.

This study has a control without polymer and four SBMA fractions (0.25, 0.5,
0.75 and 1.0), with five replicates each. The control has no
`sbma_fraction`, so the line uses four condition means:

```
trend sbma_fraction  slope -0.4  ci95 -0.71 to -0.09  p 0.03  p_adj 0.03  r2 0.94  condition_means 4  replicates 20
verdict: core_rmsf falls with sbma_fraction (slope -0.4 A per unit sbma_fraction, 95% CI -0.71 to -0.09, p_adj 0.03, fitted on 4 condition means of 20 replicates)
```

- The test needs at least three factor levels.
- With four condition means, the test has two degrees of freedom.
- Several numeric factors of one study form one {term}`Benjamini-Hochberg`
  family.
- If a replicate value is not finite, the trend is "not testable", with the
  reason.
- A factor with a level such as `PEG` or `true` gets no trend test.
- A factor with a level that is a number written as text is "not testable".
  YAML reads `1e-3` as text; write `1.0e-3`.
- Labeled results, such as one value per residue, get no trend test.

## Your own statistics

Write each test that the report does not make as a script in `stats/`. Read
one row per replicate with `replicate_table`:

```python
# Paper_1/stats/tiers.py
import polyzymd as pz

table = pz.Project("Paper_1").replicate_table("core_rmsf")
...
tier1.to_csv("stats/tier1.csv", index=False)
```

A different model, such as a dose response that is not a straight line,
belongs in such a script. The script and its output files are project files.
`project freeze` hashes them and publishes them with the paper.

## Publish

1. Fill in `metadata:` in `project.yaml`. The keys are the same as the keys of
   the metadata of a study, listed in {ref}`study-metadata`.
2. Commit every input file.
3. Freeze the project:

   ```bash
   polyzymd project freeze Paper_1
   ```

`polyzymd project freeze` does these steps:

1. It freezes each study. Each study gets its manifest, its MD checklist, its
   system summary, and the engine inputs and final frames of its replicates.
   The manifest of each study lists the files that git tracks and the commit
   of the project.
2. It writes the `manifest.json` of the project. The manifest lists:
   - the SHA-256 of the manifest of each study;
   - each condition, as `<study> / <condition>`;
   - the size and SHA-256 of each project file outside the studies:
     `project.yaml`, the shared `analyses/` and `stats/` code, and the
     figures and the files that they wrote.
3. It writes `CITATION.cff` and `.zenodo.json` from `metadata:`.
4. It commits and tags the project as `project-v1` (then `project-v2`, and so
   on). `--tag NAME` sets another tag.
5. It writes `deposit/` and `deposit/UPLOAD.md`. The paper gets one dataset
   and one DOI. The upload steps are the same as for a study; see
   {doc}`study_freeze`.

The freeze stops before it writes anything in these cases:

- An input file is not committed. Commit it, then freeze again.
- Git has no user name and email. Set them with `git config`.

The freeze also warns about each partial report, and about each run whose
stored results no longer match the project. Each warning names the run to
analyze again. A change to any of these items makes stored results stale:

- a config;
- the equilibration window or the stride;
- a function, or a file in its folder;
- a key of the analysis entry, or a setting;
- a selection;
- the content of a file that an analysis reads;
- a factor;
- a condition;
- the set of replicates.

You decide whether to freeze with warnings.

`polyzymd study freeze` refuses a study that is part of a project. Freeze the
whole project instead.
