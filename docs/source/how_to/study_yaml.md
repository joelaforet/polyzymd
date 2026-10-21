# Run a study from `study.yaml`

A {term}`study` is a set of conditions that you compare with each other.
All its conditions share one analysis frame: the same residue numbering,
reference structures, named regions, equilibration window and control.
`study.yaml` holds the analysis protocol of the study: the conditions, the
equilibration window and the settings of each analysis. The simulation of
each condition stays in its own `config.yaml`.

Use `study.yaml` when you analyze a study again and again, or when someone
else must reproduce it. To analyze several studies the same way, such as one
study for each protein, put them in a {doc}`project <project>`. For the
design of study folders and how they are published, see
{doc}`../explanation/study_folders`.

:::{admonition} Environment Setup
:class: tip

The commands on this page assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## Write `study.yaml`

Put `study.yaml` at the top of the study folder. Paths in it are relative to
the file, so you can move the folder as a whole.

```yaml
equilibration: 100ns            # the window for every analysis that sets none
stride: 1                       # optional
replicates: 1-5                 # optional; default: every replicate found
conditions:                     # control first
  No polymer: conditions/no_polymer/config.yaml
  SBMA 50%: conditions/sbma_50/config.yaml
analyses:
  rg: {}                        # a shipped analysis with its defaults
  contacts:
    method: occlusion
  contacts_4A:                  # the same analysis with other settings
    analysis: contacts
    method: distance
    cutoff: 4.0
  rg_full:                      # the same analysis over its own window
    analysis: rg
    equilibration: 0ns          # from the start of production
```

| Key | Meaning |
|---|---|
| `equilibration` | Required. The time that each analysis removes from the start of the production trajectory of each replicate |
| `conditions` | Required. Condition label to `config.yaml`. The first condition is the control. Each label is also a folder name, so two labels that differ only in case or punctuation (`SBMA 50`, `SBMA 50%`) are an error. Analysis reads a config anywhere, but `polyzymd study freeze` refuses one outside the study folder (or its project folder); `polyzymd study add-condition --config` copies one in |
| `stride`, `replicates` | Optional. The same as `--stride` and `--replicates`. `analyze` deletes the stored results of a replicate that `replicates` does not list and the run did not use |
| `until` | Optional. The end of a common analysis window, such as `38ns`, the same as `--until`. See {ref}`study-until` |
| `analyses` | Run name to settings. Each entry is one {term}`run`. The settings are the keys that `--set` takes. `analysis:` names the shipped analysis when the run name is not its name |
| `comparison` | Optional. `within:` names the factor, or a list of factors, whose values form each stratum, such as `temperature_K`. Each condition is then compared with the control of its stratum, not with the first condition. `control:` gives the control's factor values, such as `{polymer: none}`; without it, the control of a stratum is the condition whose other factors equal those of the first condition. See {ref}`study-comparison` |
| `metadata` | Publishing metadata for `polyzymd study freeze`; see {doc}`study_freeze` |

(study-comparison)=
### Compare each condition with the control of its stratum

A study that varies the temperature and the polymer gives each condition
both factors. With `comparison:`, each polymer condition is compared with the
no-polymer condition at its own temperature:

```yaml
conditions:
  none 300 K: {config: conditions/none_300_k, factors: {temperature_K: 300, polymer: none}}
  SBMA 300 K: {config: conditions/sbma_300_k, factors: {temperature_K: 300, polymer: SBMA}}
  none 360 K: {config: conditions/none_360_k, factors: {temperature_K: 360, polymer: none}}
  SBMA 360 K: {config: conditions/sbma_360_k, factors: {temperature_K: 360, polymer: SBMA}}
comparison:
  within: temperature_K
  control: {polymer: none}
```

- Every condition must give each `within` factor.
- Each stratum must have exactly one control. PolyzyMD refuses the file
  otherwise and names the stratum.
- Without `control:`, PolyzyMD takes the control's factor values from the
  first condition when it reads the file, and stores them in each report.
  A different first condition can change the control, and then makes the
  stored reports stale for `polyzymd study freeze`.
- When the control of a stratum has no values, for example no polymer atoms
  for a polymer analysis, the conditions of that stratum are summarised and
  not compared, with a warning. The other strata are compared.
- `--label` with conditions of several strata must include the control of
  each; PolyzyMD names the controls to add.
- `polyzymd study check` prints the control of each stratum as `control`.
- Each comparison line of the report names its control and its stratum, for
  example `none 360 K vs SBMA 360 K  temperature_K 360  delta ...`. In the
  JSON report, `a` is the control and `stratum` holds the `within` values.
- The Benjamini-Hochberg correction covers every comparison of the report,
  over all strata.
- Without `comparison:`, every condition is compared with the first one.
- A changed `comparison:` block makes the stored report stale for
  `polyzymd study freeze`.

PolyzyMD refuses a key that it does not know, and gives the nearest known
spelling. So a typo never falls back to a default:

```
error: study.yaml: analyses.contacts has an unknown key 'methd'.
fix: Did you mean 'method'? The keys it takes are analysis, method, ...
```

## Check the study

```bash
polyzymd study check my_study --production
```

```
study my_study/study.yaml  equilibration 100ns  stride 1
control No polymer: replicates [1, 2, 3, 4, 5] under /data/me/LipA_363K (from data.local.yaml); production 1000 ns
condition SBMA 50%: replicates [1, 2, 3, 4, 5] under /data/me/LipA_363K (from data.local.yaml); production 1000 ns
analysis rg: defaults; no stored results
analysis contacts as contacts_4A: method=distance, cutoff=4.0; no stored results
git: commit 53cd9f37d690; inputs committed
metadata (study.yaml): 3 gaps for publishing: metadata.title is missing; ...
publish: when the analyses are final, run polyzymd study freeze
cite: Laforet, Joseph R., Jr. PolyzyMD: ... (version 1.3.0). https://github.com/joelaforet/polyzymd
```

- Without `--production`, `study check` reads no trajectory.
- With `--production`, it reads the trajectory headers and segments of each
  replicate, but no frames. This takes seconds for a few replicates and
  minutes for long chains of segments. When it cannot read a replicate's
  trajectory, the condition line has no production part and no warning.
- The production length of each condition gives the longest equilibration
  window that you can use. It also shows whether the conditions were
  simulated for the same time.
- `study check` reports missing replicates, but they are not errors. So a
  study folder without its trajectories still checks.
- If a file or a config cannot be read, `study check` exits with code 2.

## Run the analyses

```bash
polyzymd analyze contacts_4A --study my_study     # one run
polyzymd analyze --study my_study                 # every run, in order
```

| Option with `--study` | Effect |
|---|---|
| `--eq`, `--stride`, `--until`, `--replicates`, `--output-dir` | Override `study.yaml` for this command |
| `--set KEY=VALUE` | Overrides one setting of the run |
| `--label LABEL` | Analyzes only that condition. Repeatable. The command prints the report but does not save it as the `report.json` of the run |
| `--run NAME` | Selects the {term}`result` to report, when the analysis reports several |
| `--submit`, `--dry-run` | One SLURM array task per condition and replicate, then a report job, as in {doc}`hpc_execution` |
| `-c` | Refused, because `study.yaml` names the conditions |

If you name a shipped analysis that `study.yaml` does not list, it runs with
its defaults, and the command prints a note.

The results of each run go to `my_study/results/<run>/`:

| Path | Holds |
|---|---|
| `polyzymd_results/` | The stored per-replicate values and their records |
| `report.json` | The full report |
| `figures/` | The figures |
| `slurm/` | The scripts and logs of `--submit` |

Two runs of one analysis keep separate stored results. So a change from
`contacts` to `contacts_4A` recomputes nothing.

The console shows only the report and the warnings. The full log, with
library messages, goes to `my_study/logs/polyzymd-analyze-<time>.log`. The
first line of the output gives its path. To show the full log on the
console, use `polyzymd -v analyze ...`.

(study-until)=
### Compare conditions over the same time

The conditions can have different simulated lengths, or different lengths of
data on this machine. A difference between them can then come from the
simulated time, not from the condition. For example, a structure that drifts
late shows only in the longer simulations. The report then warns:

```
warning: the conditions were analysed up to different times (No polymer 456.4 ns; SBMA 50% 38.4 ns); a difference may come from simulated time rather than the condition. Compare over a common window with until 38.4ns (--until, or until: in study.yaml)
```

- `until: 38ns` in `study.yaml`, or `--until 38ns`, removes the production
  after 38 ns from every condition.
- `until: common`, or `--until common`, stops every replicate at the last
  frame of the shortest replicate. All replicates then have the same time
  points, also when segment restarts left them a frame or two apart.
  PolyzyMD finds that time once, over all conditions, and records it in each
  result.

### Give one analysis its own window

An analysis entry can set its own `equilibration:`, `until:` and `stride:`.
Use a stride for an analysis that costs too much to run on every frame.

Use a separate window, for example, for a system that still changes at the
end of production. A time-resolved analysis can then start at 0 ns, while a
steady-state analysis keeps the window of the study.

- Two entries of one analysis with different windows store their results
  side by side.
- Each record names its window.
- `polyzymd study check` prints the window of each analysis.

The window comes from the first of these that is set:

1. the command line (`--eq`, `--until`, `--stride`);
2. the analysis entry;
3. the study.

When a window is part of the protocol, give it in the entry, not on the
command line. Otherwise the stored results do not match `study.yaml`, and
freeze reports them as stale.

The report also names each replicate whose `progress.json` lists production
segments that are not on disk. This happens with a copy that kept only some
segments.

## Run your own function

First check whether a shipped analysis measures your quantity.
`polyzymd analyze --list` prints each shipped analysis, what it measures and
its settings with their defaults. For example, `contacts` with
`method: occlusion` is the buried-surface definition of a contact. A shipped
analysis is verified and documented, and its settings go straight into an
entry of `study.yaml`.

To run your own function, write it in a Python file of the study and list it
in `study.yaml`. It then runs like a shipped analysis, with stored records,
statistics, a report and figures:

```python
# my_study/analyses/lid.py
import numpy as np


def lid_distance(lid, core):
    """Distance between the centers of the lid and the core, in Å."""
    return float(np.linalg.norm(lid.center_of_geometry() - core.center_of_geometry()))
```

```yaml
analyses:
  lid_opening:
    function: analyses/lid.py:lid_distance
    kind: timeseries
    unit: A
    selections:
      lid: "protein and resid 140-150 and name CA"
      core: "protein and resid 4-120 and name CA"
```

Each key of `selections:` is a keyword argument of the function. Here,
`lid` and `core` become the parameters `lid` and `core` of `lid_distance`.

| Key | Meaning |
|---|---|
| `function` | `file.py:function_name`, relative to `study.yaml` |
| `kind` | `timeseries` or `per_replicate`; see the next table |
| `unit` | The unit of the values |
| `selections` | Keyword argument to selection. PolyzyMD passes each as an `AtomGroup` of the replicate |
| `universe` | The keyword argument that receives the `Universe` of the replicate |
| `settings` | Keyword arguments that PolyzyMD passes unchanged. `--set` overrides them |
| `labels: returned` | For a `per_replicate` function that returns `(labels, values)`, such as one value per residue |
| `parts: [area, contacts]` | For a function that measures several quantities in one pass; see below |
| `missing: .nan` | With `labels: returned`, the value for a label that other replicates returned and this replicate did not; see below |
| `allow_empty: true` | Passes a selection that matches no atoms as an empty `AtomGroup`; see below |

| `kind` | PolyzyMD calls the function | The function returns |
|---|---|---|
| `timeseries` | On each production frame | One number. PolyzyMD then averages the series of each replicate (`reduce: mean`) |
| `per_replicate` | Once per replicate, with `frames=` the production frame indices, and `times=` their times in ns if the function has a `times` parameter | One value |

**`parts:`** A `timeseries` function returns a dict with these keys at each
frame, or a sequence in this order. A `per_replicate` function returns one
row per part. PolyzyMD stores and plots each part as its own result. Each
part has its own `part` value in `Study.results().table`. The report covers
the part that `--run` names, or the first part.

**`missing:`** An example is a frame index past the end of a shorter
replicate. The report warnings name each replicate that gets the `missing`
value, with its labels. Without `missing:`, such a replicate stops the
report, with a message.

**`allow_empty: true`** An example is a polymer selection in a control
without polymer. The function then decides the value, for example `0.0` when
`len(polymer) == 0`. Without `allow_empty`, such a replicate stops the run,
with a message. The shipped `contacts` and `hydrogen_bonds` analyses need no
such setting: they report 0 for a replicate without polymer.

### Label per-frame values by time

For a per-frame label, such as one value per residue and frame, name the
frame by its time, not by its index. Restarts can leave duplicate frames,
which PolyzyMD removes. After that, frame index *k* is not the same time in
every replicate.

1. Give the function a `times` parameter.
2. Label each value as `f"{resid}|{t:.3f}"`.
3. Add `until: common`, so that every replicate covers the same times.

For a few quantities per frame, use `kind: timeseries` with `parts:`. It
keeps the time axis for you.

### When stored results are recomputed

The stored results of a function depend on the function and on these files
in the folder of the function, in subfolders too. A change to any of them
recomputes the results at the next analysis:

- every Python file: the file of the function, and a helper module or
  package that the file imports from that folder;
- every file in a `data/` folder beside the function file, such as
  `analyses/data/reference_distances.csv`. Put each data file that a
  function reads there.

No other file counts. Notes, figures, PDFs or a copied trajectory in
`analyses/` change no result, and PolyzyMD reads none of them. The analysis
log names the files that it leaves out, so a data file outside `data/` is
easy to find. Hidden folders (`.git`, `.pixi`, `.venv`), `results/`,
`logs/`, `deposit/`, `conditions/`, `figures/`, job folders and the folders
of other studies never count.

A function file beside `study.yaml` or `project.yaml` depends only on the
Python files of that folder, because a `data/` folder there can hold
trajectories. Keep the data files of such a function in `analyses/data/`, or
pass them as arguments. PolyzyMD compiles the files from their current text
each time, never from a cached `.pyc`.

`polyzymd study check` imports each listed function. So it reports a broken
file before it reads any trajectory.

## Read results back for figures

`Study("study.yaml")` reads the stored results without loading any
trajectory. So a figure script or notebook in `figures/` works on a copy of
the study folder that has no trajectories:

```python
import polyzymd as pz

study = pz.Study("my_study/study.yaml")
results = study.results("lid_opening")

results.report                 # the ProtocolReport, as polyzymd analyze printed it
results.table                  # one row per stored value
study.replicate_table("lid_opening").groupby("condition")["value"].mean()
```

`study.results` reads `results/<run>` in the study folder. To read results
that `polyzymd analyze --output-dir` wrote elsewhere, use
`study.results("lid_opening", folder="that/folder")`. `analyze` warns when
`--output-dir` puts results outside the study folder, because `study check`
and `study freeze` see only `results/`.

`results.table` has the columns `name`, `condition`, `replicate`, `part`,
`label`, `frame`, `time_ns`, `value` and `unit`:

- A per-frame series fills `frame` and `time_ns`. The table then has one row
  per frame.
- A per-residue result fills `label`.
- A result with several parts, such as the contact fraction and its
  components, fills `part`.

`study.replicate_table(run)` has one row per replicate, quantity (`name`),
part and label. Every test of the report uses these values.

`study.settings("contacts_4A")` returns the settings of a run. When the
trajectories are present, the conditions load as usual
(`study["SBMA 50%"].replicates`).

### The Study object at a glance

| You want | Write |
|---|---|
| A study from its file | `study = pz.Study("my_study/study.yaml")` (`polyzymd.Study`) |
| Condition labels, control first | `study.labels`, `study.control` |
| The conditions | `study.conditions`, `for condition in study`, `study["SBMA 50%"]` |
| The replicates of a condition | `study["SBMA 50%"].replicates`. Each has `.index`, `.frames`, `.times`, `.universe()` |
| The stored results of a run | `study.results("rg")` (`.table`, `.report`, `.warnings`) |
| One row per replicate | `study.replicate_table("rg")`: one row per replicate, quantity (`name`), part and label |
| The settings of a run | `study.settings("rg")` |
| A file of the study, from a relative path in its settings | `study.path(study.settings("rmsf")["reference_file"])`. `study.root` is the folder |
| The Python code of the study, as the analyses import it | `study.module("interface")` for `analyses/interface.py` (also looked up in the project folder) |
| Several studies | `pz.Project("Paper_1")`, with the same `results` and `replicate_table`, and a `study` column |
