# Run a study from `study.yaml`

Use this when a study has several conditions and analyses that you rerun
often, or that someone else must reproduce. `study.yaml` holds the analysis
protocol: the conditions, the equilibration window and every analysis
setting. Each condition's simulation stays in its own `config.yaml`. For the
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

Put it at the top of the study folder. Paths are relative to the file, so
the folder can be moved as a whole.

```yaml
polyzymd: 1.3.0                 # the version that produced the results
equilibration: 100ns            # one window for every condition and replicate
stride: 1                       # optional
replicates: 1-5                 # optional; default: every run found
conditions:                     # control first
  No polymer: conditions/no_polymer/config.yaml
  SBMA 50%: conditions/sbma50/config.yaml
analyses:
  rg: {}                        # a shipped analysis with its defaults
  contacts:
    method: occlusion
  contacts_4A:                  # the same analysis with other settings
    analysis: contacts
    method: distance
    cutoff: 4.0
```

| Key | Meaning |
|---|---|
| `equilibration` | Required. Window removed from the start of every replicate's production trajectory |
| `conditions` | Required. Condition label to `config.yaml`, control first |
| `stride`, `replicates` | Optional, as `--stride` and `--replicates` |
| `until` | Optional end of a common analysis window, such as `38ns`, as `--until`: production after it is left out for every condition |
| `analyses` | Run name to settings. The settings are those `--set` takes; `analysis:` names the shipped analysis when the run name is not one |
| `polyzymd` | The PolyzyMD version the study was run with; a different version gives a warning |
| `metadata` | Publishing metadata, read by `polyzymd study freeze`; see {doc}`study_freeze` |

A key PolyzyMD does not know is refused with the nearest known spelling, so a
typo never falls back to a default:

```
error: study.yaml: analyses.contacts has an unknown key 'methd'.
fix: Did you mean 'method'? The keys it takes are analysis, method, ...
```

## Check it

```bash
polyzymd study check my_study
```

```
study my_study/study.yaml  equilibration 100ns  stride 1
control No polymer: runs [1, 2, 3, 4, 5] under /pl/active/.../LipA_363K_REDO (from config); production 1000 ns
condition SBMA 50%: runs [1, 2, 3, 4, 5] under /pl/active/.../LipA_363K_REDO (from config); production 1000 ns
analysis rg: defaults; no stored results
analysis contacts as contacts_4A: method=distance, cutoff=4.0; no stored results
git: commit 53cd9f37d690; inputs committed
metadata: 3 gaps for publishing; polyzymd study freeze lists them
publish: when the analyses are final, run polyzymd study freeze
cite: Laforet, Joseph R., Jr. PolyzyMD: ... (version 1.3.0). https://github.com/joelaforet/polyzymd
```

It reads trajectory headers, not frames, so it takes seconds. Each
condition's production length tells you how long an equilibration window can
be, and whether the conditions were simulated for the same time. Missing runs
are reported but are not errors, so a study folder without its trajectories
still checks. An unreadable file or config exits 2.

## Run the analyses

```bash
polyzymd analyze contacts_4A --study my_study     # one run
polyzymd analyze --study my_study                 # every run, in order
```

| Option with `--study` | Effect |
|---|---|
| `--eq`, `--stride`, `--until`, `--replicates`, `--output-dir` | Override the file for this command |
| `--set KEY=VALUE` | Overrides one setting of the run |
| `--label LABEL` | Runs only that condition; repeatable |
| `--submit`, `--dry-run` | One SLURM array task per condition and replicate, then a report job, as in {doc}`hpc_execution` |
| `-c` | Refused: the file names the conditions |

A shipped analysis the file does not list runs with its defaults and a note.

Each run's results go to `my_study/results/<run>/`:

| Path | Holds |
|---|---|
| `polyzymd_results/` | The stored per-replicate values and their records |
| `report.json` | The report, in full |
| `figures/` | The figures |
| `slurm/` | The scripts and logs of `--submit` |

Two runs of one analysis keep separate stored results, so switching between
`contacts` and `contacts_4A` recomputes nothing.

The console shows only the report and warnings. The full log, with library
messages, goes to `my_study/logs/polyzymd-analyze-<time>.log`, whose path is
the first line printed; `polyzymd -v analyze ...` shows it on the console.

### Compare conditions over the same time

When the conditions were simulated, or are on this machine, for different
lengths, a difference can come from simulated time rather than the
condition: a structure that drifts late appears only in the longer runs. The
report then warns:

```
warning: the conditions were analysed up to different times (No polymer 456.4 ns; SBMA 50% 38.4 ns); a difference may come from simulated time rather than the condition. Compare over a common window with until 38.4ns (--until, or until: in study.yaml)
```

`until: 38ns` in `study.yaml`, or `--until 38ns`, leaves out production after
that time for every condition. A replicate whose `progress.json` records
production segments missing from disk, as in a copy that kept only some
segments, is named in the report too.

## Run your own function

List a function from a Python file of the study, and it runs like a shipped
analysis, with stored records, statistics, a report and figures:

```python
# my_study/analyses/lid.py
import numpy as np


def lid_distance(lid, core):
    """Distance between the centres of the lid and the core, in Å."""
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

| Key | Meaning |
|---|---|
| `function` | `file.py:function_name`, relative to `study.yaml` |
| `kind` | `timeseries`: called on every production frame, returning one number; each replicate's series is then averaged (`reduce: mean`). `per_replicate`: called once per replicate with `frames=` the production frame indices, returning one value |
| `unit` | Unit of the values |
| `selections` | Keyword argument to selection; each is passed as the replicate's `AtomGroup` |
| `universe` | Keyword argument that receives the replicate's `Universe` |
| `settings` | Keyword arguments passed as they are; `--set` overrides them |
| `labels: returned` | For a `per_replicate` function that returns `(labels, values)`, such as one value per residue |
| `allow_empty: true` | Leave out, with a warning in the report, every replicate where a selection matches no atoms, such as a polymer selection in a no-polymer control. Without it such a replicate stops the run, with a message saying so |

The stored results are keyed on the whole file, not only the function: edit
any helper in it and the next run recomputes. The file is compiled from its
current text every time, never from a cached `.pyc`.

`polyzymd study check` imports every listed function, so a broken file is
reported before any trajectory is read.

## Read results back for figures

`Study("study.yaml")` reads the stored results without loading any
trajectory, so a figure script or notebook in `figures/` works on a copy of
the study folder that has no trajectories:

```python
import polyzymd as pz

study = pz.Study("my_study/study.yaml")
results = study.results("lid_opening")

results.report                 # the ProtocolReport, as polyzymd analyze printed it
results.table                  # one row per stored value
results.table.groupby(["condition", "replicate"])["value"].mean()
```

`results.table` has the columns `name`, `condition`, `replicate`, `part`,
`label`, `frame`, `time_ns`, `value` and `unit`. A per-frame series fills
`frame` and `time_ns`; a per-residue result fills `label`; a result with
several parts, such as the contact fraction and its components, fills `part`.

`study.settings("contacts_4A")` returns a run's settings, and the conditions
load as usual (`study["SBMA 50%"].replicates`) when the trajectories are
present.
