# Study API

The Python interface that runs a function on every replicate of a study. For
a worked example, see {doc}`../how_to/study_api`. For the reasons behind the
design, see {doc}`../explanation/analysis_api`. The module documentation is
in {doc}`../api/analyses`.

`import polyzymd as pz` exports `Study`, `Project`, `select`, `universe`,
`reference`, `plot_values` and `plot_distributions`.

## Load a study

| Call | Returns |
|---|---|
| `pz.Study.from_configs(configs, *, equilibration, replicates=None, stride=1, data=None, until=None)` | A `Study` from a mapping of condition label to `config.yaml`, control first. A list of paths labels each condition by the folder of its config |
| `pz.Study("study.yaml")` | A `Study` from a study file, or the folder that holds it |
| `pz.Project("Paper_1")` | A `Project`: every study of `project.yaml`. Give the project folder or its `project.yaml` |

| `from_configs` argument | Meaning |
|---|---|
| `equilibration` | The time removed from the start of the production trajectory of each replicate, such as `"100ns"` |
| `replicates` | The replicate numbers of every condition. Default: the replicate folders found |
| `stride` | Keep every `stride`-th production frame, from the first frame after the window |
| `data` | Condition label to the folder that holds its replicate folders, in place of the `scratch_directory` of the config. The key `"*"` applies to every condition that the mapping does not name |
| `until` | The end of a common analysis window, such as `"38ns"` or `"common"` |

## Study, condition and replicate

| Attribute or method | Gives |
|---|---|
| `study.labels`, `study.control` | The condition labels, control first, and the control label |
| `study.conditions`, `study[label]`, `for condition in study` | The conditions |
| `study[label].replicates` | The replicates of a condition |
| `replicate.index` | The replicate number |
| `replicate.universe()` | The MDAnalysis `Universe`, with the production segments in order |
| `replicate.frames` | The trajectory frame indices after the window and the stride, counted from 0 |
| `replicate.times` | The times of those frames, in ns |
| `replicate.production_ns` | The simulated production time of the replicate, in ns |
| `study.results(run, *, folder=None)` | The stored results of a run of a study file, read without trajectories (`.table`, `.report`, `.warnings`) |
| `study.replicate_table(run)` | One row per replicate, quantity (`name`), part and label, with a column for each factor |
| `study.settings(run)` | The settings of a run of `study.yaml` |
| `study.path(relative)`, `study.root` | A file of the study, and the study folder |
| `study.module(name)` | A Python module of the study, as the analyses import it |

A frame that you name yourself, such as the `frame` of `pz.reference` or
`--set reference_frame=N`, is a production frame counted from 1 after the
window. So `1` is the first production frame. The record stores the
trajectory frame that PolyzyMD used.

## Project

| Attribute or method | Gives |
|---|---|
| `project.root` | The project folder |
| `project.labels` | The study labels, in the order of `project.yaml` |
| `project[label]`, `for study in project`, `len(project)` | The studies, each a `Study` |
| `project.runs_in(run)` | The labels of the studies that run the analysis `run` |
| `project.results(run)` | The stored results of `run` in every study that runs it, read without trajectories: `.table` with a `study` column first, `.reports` and `.folders` by study label |
| `project.replicate_table(run)` | One row per replicate of `run` in every study that runs it, with a `study` column |

`project.results` and `project.replicate_table` stop with an error that names
each study that runs `run` but has no stored results.

## Arguments that stand for each replicate

| Placeholder | Becomes, in each replicate |
|---|---|
| `pz.select(selection, *, allow_empty=False)` | `universe.select_atoms(selection)`. An empty selection is an error unless `allow_empty=True` |
| `pz.universe()` | The `Universe` |
| `pz.reference(mode, selection, *, frame=None, file=None, alignment=None)` | The reference atoms, in a separate one-frame universe |

| `reference` mode | The reference is |
|---|---|
| `external` | The `selection` atoms of `file`. The record holds the SHA-256 of the file |
| `frame` | Production frame `frame`, counted from 1 |
| `average` | The mean positions. PolyzyMD superposes the `alignment` atoms of every frame on the first frame, averages, superposes again on that average, and takes the mean |
| `centroid` | The production frame whose `alignment` atoms have the smallest RMSD to their iterative average structure |

`alignment` defaults to `selection`. Giving `file` with a mode other than
`external` is an error.

## Measure

| Method | Calls the function | The function returns | Returns |
|---|---|---|---|
| `study.timeseries(function, *args, unit, name=None, recompute=False, output_dir=None, bounds=(None, None), parts=None, **kwargs)` | Once per production frame, through MDAnalysis `AnalysisFromFunction` | One number, or with `parts` a dict or sequence of one number per part | `Timeseries`, or with `parts` a dict of `Timeseries` |
| `study.per_replicate(function, *args, unit, labels=None, missing=None, note_filled=False, name=None, recompute=False, output_dir=None, bounds=(None, None), parts=None, **kwargs)` | Once per replicate, with the keyword argument `frames` | One number, a one-dimensional array with one entry per label, or with `parts` one row per part | `ReplicateValues`, or with `parts` a dict of `ReplicateValues` |

| Argument | Meaning |
|---|---|
| `unit` | The unit of the value. `None` for a dimensionless quantity |
| `name` | The result name, used for the folder and the report. Default: the name of the function |
| `recompute` | Measure every replicate, even when a stored result matches |
| `output_dir` | The folder that holds `polyzymd_results/`. Default: the current folder |
| `bounds` | The lowest and highest possible value, `None` for no limit. Used for distribution figures and for a warning when an interval passes a limit. Not part of the record |
| `parts` | The names of several quantities measured in one pass. Each part becomes its own result |
| `labels` | `per_replicate` only. The name of each entry: a list, a function of the `Universe`, or `"returned"` when the function returns `(labels, values)` |
| `missing` | `per_replicate` only. The value for a label that one replicate lacks. Default: a missing label is an error |
| `note_filled` | `per_replicate` only. Name each replicate given `missing` in the report warnings |
| `**kwargs` | Keyword arguments of the function, recorded with the result |

## Timeseries

| Method | Returns |
|---|---|
| `reduce(how="mean", *, unit=..., bounds=..., detect_equilibration=True)` | `ReplicateValues`, one value per replicate |
| `transform(function, *others, unit=..., name=None, bounds=..., **kwargs)` | A new `Timeseries` computed from the stored values, frame by frame. It reads no trajectory |
| `plot(output_dir=None, name=None, plot_settings=None)` | The path of a figure of every replicate's series against time, with the mean and 95 % interval of each condition |
| `plot_distribution(threshold=None, output_dir=None, name=None, title=None, plot_settings=None)` | The path of a figure of the distribution of frame values of each condition, pooled and per replicate |

| `how` | One value per replicate | Default `unit` | Default `bounds` |
|---|---|---|---|
| `"mean"` | The mean over frames | The series unit | The series bounds |
| `"fraction"` | The mean of a series of 0 and 1. Other values are refused | `None` | `(0, 1)` |
| `"std"` | The sample standard deviation over frames (`ddof=1`) | The series unit | `(0, None)` |
| a function `f(values, times)` | Its return value. `times` are in ns | The series unit | None |

`detect_equilibration=True` runs `pymbar.timeseries.detect_equilibration` on
each series as a diagnostic. It changes no value.

`transform` calls `function(values, *other_values, **kwargs)` for each
replicate. Each series in `others` must have the same frames. Pass a value
such as a threshold in `kwargs`: a value from an enclosing scope is not
recorded.

## ReplicateValues

| Method | Returns |
|---|---|
| `summary(conditions=None)` | A `ProtocolReport` with n, the mean, the standard error, the 95 % interval and every replicate value of each condition |
| `compare(control=None, conditions=None, test="welch", untested=(), within=None)` | A `ProtocolReport` with each condition against the control, or against the control of its stratum |
| `over_labels(how="mean", metric=None, labels=None)` | One number per replicate from a labelled array, over every label or over `labels` |
| `plot(output_dir=None, name=None, title=None, plot_settings=None, highlight=(), xlabel="label")` | The path of a figure of the mean and interval of each condition, with every replicate value. For a labelled result, a profile |
| `values` | The stored values |

`compare` details:

- `control` defaults to the first condition. `test` is `"welch"` or
  `"student"`.
- `within`, a factor name or a list, compares each condition with the
  control of its stratum: the conditions with the same values of those
  factors. The control of a stratum is the condition whose other factors
  equal those of `control`, or, for `control={"polymer": "none"}`, whose
  factors have those values. The factors come from `study.yaml`. `within`
  and `control` default to its `comparison:` block, and `within=` given
  without `control=` takes the block's control too; `within=[]` turns it
  off. Without a block, `control` defaults to the first condition of the
  study, also when that condition has no values.
- A stratum with no control or two is refused. So is a control or a compared
  condition with no replicate values, and a control that is not among the
  conditions. A condition is never compared with another stratum's control.
- Each row gives `mean(b) - mean(a)` with its 95 % interval, the p value, the
  {term}`Benjamini-Hochberg` adjusted p value, and Cohen's d and Hedges' g.
- The correction family is every tested row of the call, over every stratum
  with `within`. For a labelled result, that is every label of every compared
  condition.
- A row is not testable when a condition has fewer than two replicates, or
  when both conditions have the same value in every replicate. It then takes
  no part in the correction.
- Labels in `untested` are summarized but not tested.

`ProtocolReport.to_agent_text()` prints the report in the format of
{ref}`polyzymd analyze <cli-analyze>`. `model_dump_json()` gives the JSON. For
every field, see {doc}`analysis_protocol_report`.

## Figures of several results

| Function | Draws |
|---|---|
| `pz.plot_values(results, labels=None, output_dir=None, name="values", title=None, plot_settings=None)` | One group of bars per `ReplicateValues`, one bar per condition, with intervals and replicate points. All results must share a unit |
| `pz.plot_distributions(series, thresholds=None, titles=None, output_dir=None, name="distributions", quantity="value", plot_settings=None)` | One panel per `Timeseries`, with a shared axis and a threshold per panel |
| `polyzymd.analyses.figures.plot_differences(values, report, output_dir, name, title=None, plot_settings=None, xlabel="label")` | One panel per condition: the difference from the control at each label, its 95 % interval from the `compare()` report, and a point on each significant label |
| `polyzymd.analyses.figures.plot_decomposition(parts, output_dir, name, title=None, plot_settings=None, xlabel="label")` | Several labelled results of one unit, one panel per condition |

Every figure goes to `figures/` beside `polyzymd_results/`, or to
`output_dir`. A figure with error bars or a band has a footnote that names
the interval and the replicates it is computed across.

## Stored results

Each replicate of a result is stored in
`<output_dir>/polyzymd_results/<name>/<condition>/replicate_<n>/`:

| File | Holds |
|---|---|
| `series.npz` | `timeseries`: the values, frames and times |
| `values.npz` | `per_replicate`: the values |
| `labels.json` | `per_replicate` with `labels="returned"`: the labels of the replicate |
| `record.json` | What produced the values. See below |

`record.json` holds:

- the name and module of the function, and the hash of its code;
- every argument, with the selection strings;
- the config hash;
- the relative path, size and SHA-256 of every topology and trajectory file;
- the equilibration window, the frames and the times;
- the unit and the labels;
- the PolyzyMD, MDAnalysis, NumPy and Python versions.

A new call reuses a stored result when every field except the versions and
the bounds is equal. The hash of the code depends on the function:

| Function | Hash of |
|---|---|
| A function of a study or project file (`function:` in `study.yaml`) | Every file in the folder of its file |
| A shipped function of `polyzymd.analyses` | Its module file and the PolyzyMD modules that the file imports |
| Any other function | Its source code |
| A function without a source file, such as a lambda or a notebook cell | Its bytecode, with a warning |
