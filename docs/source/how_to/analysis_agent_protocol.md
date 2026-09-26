# Get a validated number with one command

Use this when you want the value of a metric for one simulation, or want to
know whether two simulations differ, and you do not want to write a comparison
YAML first. One command gives a number that carries its unit, its uncertainty,
its sample size and its provenance.

## Summarize one condition

```bash
pixi run -e analysis polyzymd analyze rg -c enzyme_water/config.yaml --eq 10ns
```

```
# polyzymd analyze rg  metric mean_rg  unit A  eq 10ns  conditions 1  replicates 3  protocol rg/2
enzyme_water  n 3  mean 18.42  sem 0.05  ci95 18.2 to 18.64  values 18.4, 18.5, 18.36  replicates 1,2,3  g 473.7  n_eff 19
verdict: enzyme_water mean_rg 18.42 A (95% CI 18.2 to 18.64, n 3)
```

The replicates come from the directories the config points at. Name them
explicitly with `--replicates 1-3` when only some of them should be used.

## Compare two or more conditions

The first `-c` is the control; every other condition is compared against it.

```bash
pixi run -e analysis polyzymd analyze rg \
  -c noPoly/config.yaml \
  -c SBMA50/config.yaml \
  --label "no polymer" --label "50% SBMA" \
  --eq 10ns
```

```
# polyzymd analyze rg  metric mean_rg  unit A  eq 10ns  conditions 2  replicates 3,3  protocol rg/2
no polymer  n 3  mean 18.42  sem 0.05  ci95 18.2 to 18.64  values 18.4, 18.5, 18.36  replicates 1,2,3  g 473.7  n_eff 19
50% SBMA  n 3  mean 18.73  sem 0.06  ci95 18.47 to 18.99  values 18.71, 18.8, 18.68  replicates 1,2,3  g 402.1  n_eff 22
no polymer vs 50% SBMA  delta +0.31  ci95 0.02 to 0.6  p 0.041  p_adj 0.041  test welch_t  correction BH  d 1.9  significant
verdict: 50% SBMA larger mean_rg than no polymer (delta +0.31 A, 95% CI 0.02 to 0.6, p_adj 0.041, n 3 vs 3)
```

Without `--label`, each condition is named after the directory holding its
config.

## Change plugin settings

Pass `--set key=value` once per setting. Values are read as YAML, so numbers
and booleans arrive with the right type. `--set` takes only top-level settings;
put nested ones in a comparison.yaml and pass it with `-f`.

```bash
pixi run -e analysis polyzymd analyze rmsf \
  -c A/config.yaml -c B/config.yaml \
  --set selection='name CA' --set reference_mode=average
```

## Pick a selection when the analysis measures several

Some analyses report one metric for several selections: sasa measures four
contexts and distances measures each atom pair. The report covers one of them at a time, names it in `run`, and
lists the rest in `all_runs`. Pick another with `--run`:

```bash
pixi run -e analysis polyzymd analyze distances -c A/config.yaml -c B/config.yaml \
  --run "Substrate-Ser77"
```

## Get the full record

`--format json` prints the whole `ProtocolReport`, including every metric the
plugin reported, the frames each replicate contributed, the package versions
and a SHA-256 of each simulation config. Every field is listed in
{doc}`../reference/analysis_protocol_report`.

```bash
pixi run -e analysis polyzymd analyze rg -c A/config.yaml -c B/config.yaml \
  --format json -o rg_report.json
```

In Python the same call is:

```python
from polyzymd.analyses.protocols import analyze

report = analyze("rg", ["A/config.yaml", "B/config.yaml"], equilibration="10ns")
print(report.to_agent_text())
print(report.conditions[0].ci95, report.unit)
```

## Reuse an existing comparison project

When a `comparison.yaml` already exists, analyze it directly instead of
repeating its conditions on the command line:

```bash
pixi run -e analysis polyzymd analyze rmsf -f comparison.yaml
```

The same summary is available from the full comparison workflow, which also
writes plots:

```bash
pixi run -e analysis polyzymd compare run rmsf -f comparison.yaml --format agent
```

## Read the result

- `n` is the number of replicates. It is the sample size for every test; frames
  are never treated as independent samples.
- `ci95` on a condition is the Student t interval on its mean. On a comparison
  it is the interval on the difference and carries no multiplicity correction,
  so it can exclude zero while `p_adj` does not clear alpha.
- `not testable` means a condition has fewer than two replicates, so no test
  exists. It does not mean the conditions are the same.
- `no test recorded` means the plugin stored a raw p value but no
  multiplicity-corrected one, so the line describes a difference without
  deciding it. Distances is the analysis that does this today.
- A plugin that stores only means and standard errors, such as contacts, gets
  its condition intervals rebuilt from the standard error and the replicate
  count; `ci_method` then reads `student_t_from_sem` and no interval is given
  on a difference.
- Every `warning:` line is part of the answer. A warning that a condition has
  two replicates changes how wide the interval really is.

## Where the outputs land

The run writes `analysis/`, `comparison/` and `figures/` under the current
directory, or under `--output-dir` when you give one. Cached replicate results
are reused; pass `--recompute` to rebuild them.

## Troubleshooting

| Message | Fix |
|---|---|
| `error: Unknown analysis 'rgyr'` | Use a name from `polyzymd compare run --list`. |
| `error: ... no replicate directories under the scratch directory` | Pass `--replicates 1-3` to state which replicates to use. |
| `error: Config file(s) not found` | Point `-c` at a simulation `config.yaml`, not a `comparison.yaml`. |
| `polyzymd: command not found` | Run through `pixi run -e analysis`. |

The command exits 0 on success and 2 on any of the errors above, printing the
message and the fix on one line each.

## See also

- {doc}`../explanation/analysis_entry_points` for which entry point to use.
- {doc}`analysis_compare_conditions` for the full comparison workflow.
- {doc}`../reference/cli_reference` for every flag.
