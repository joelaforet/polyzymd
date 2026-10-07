# Get one number with its uncertainty

Run one `polyzymd analyze` command to get the value of a metric for one
condition, or to find out whether two conditions differ. The command reads the
`config.yaml` of each condition. Its result carries the unit, the 95 %
confidence interval, the number of replicates and the provenance.

## Summarize one condition

```bash
pixi run -e analysis polyzymd analyze rg -c enzyme_water/config.yaml --eq 10ns
```

```
# polyzymd analyze rg  metric mean_rg  unit A  eq 10ns  conditions 1  replicates 3  protocol rg/2
enzyme_water  n 3  mean 18.42  sem 0.05  ci95 18.2 to 18.64  values 18.4, 18.5, 18.36  replicates 1,2,3  g 473.7  n_eff 19
verdict: enzyme_water mean_rg 18.42 A (95% CI 18.2 to 18.64, n 3)
```

PolyzyMD finds the replicates in the folders that the config names. To use
only some replicates, add `--replicates 1-3`.

## Compare two or more conditions

The first `-c` is the control. PolyzyMD compares each other condition with
it.

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

Without `--label`, each condition takes the name of the folder that holds its
config.

## Change analysis settings

Give `--set key=value` once for each setting. PolyzyMD reads each value as
YAML, so numbers and booleans keep their type. `--set` takes only top-level
settings. Give a nested setting as one YAML mapping, for example
`--set groups='{protein: chainid A, polymer: chainid C}'`.

```bash
pixi run -e analysis polyzymd analyze rmsf \
  -c A/config.yaml -c B/config.yaml \
  --set selection='name CA' --set reference_mode=average
```

## Select one result

Some analyses report several results:

- `sasa` reports the total and the per-residue SASA of each context.
- `distances` reports a mean distance and a fraction below the threshold for
  each atom pair.
- `rmsf` reports core, region and mean values, and the per-residue profiles.

The report covers one result at a time. It names that result in `run` and
lists the others in `all_runs`. To report a different result, use `--run`:

```bash
pixi run -e analysis polyzymd analyze distances -c A/config.yaml -c B/config.yaml \
  --set pairs=pairs.yaml --run "Substrate-Ser76"
```

## Get the full record

`--format json` prints the whole `ProtocolReport`. It includes these items:

- every metric that the analysis reported;
- the frames that each replicate contributed;
- the package versions;
- the SHA-256 of each simulation config.

For every field, see {doc}`../reference/analysis_protocol_report`.

```bash
pixi run -e analysis polyzymd analyze rg -c A/config.yaml -c B/config.yaml \
  --format json -o rg_report.json
```

In Python the same call is:

```python
from polyzymd.analyses import analyze

report = analyze("rg", ["A/config.yaml", "B/config.yaml"], equilibration="10ns")
print(report.to_agent_text())
print(report.conditions[0].ci95, report.unit)
```

The catalytic triad is not an analysis of this command. Measure it with the
routine in {doc}`analysis_triad_quickstart`. For the triad distances alone, run
`polyzymd analyze distances --set pairs=<pairs.yaml>`.

An agent learns this protocol from
`.claude/skills/polyzymd-analyze/SKILL.md` in the repository, or from this
page.

## Read the result

For the format of each line and the verdict words, see
{ref}`polyzymd analyze <cli-analyze>`. Read these points first:

- `n` is the number of replicates. It is the sample size of every test.
  PolyzyMD never treats frames as independent samples.
- `ci95` of a condition is the Student t interval of its mean. `ci95` of a
  comparison is the interval of the difference. It has no correction for
  multiple tests. So it can exclude zero while `p_adj` is above alpha.
- `not testable` means that a condition has fewer than two replicates, so no
  test is possible. It does not mean that the conditions are the same.
- `no test recorded` means that a stored comparison holds a raw p value but no
  corrected p value. The line gives the difference but no decision.
- A stored comparison that holds only means and standard errors gets condition
  intervals rebuilt from the standard error and the number of replicates.
  `ci_method` is then `student_t_from_sem`, and the difference has no interval.
- Each `warning:` line is part of the answer. For example, a warning that a
  condition has two replicates changes how much the interval tells you.

## Where the outputs land

The command writes two folders in the current folder, or in `--output-dir`:

- `polyzymd_results/` holds the values of each replicate and their record.
- `figures/<name>/` holds the figures.

PolyzyMD reuses a stored replicate result if these items are unchanged: the
function, the settings, the config, the input files, the equilibration window
and the frames. To measure it again, add `--recompute`.

## Troubleshooting

| Message | Fix |
|---|---|
| `error: No analysis named 'rgyr'.` | Use one of the names that the `fix:` line lists. |
| `error: Condition ...: replicates [...] have no run directory under ...` | Give replicates that exist with `--replicates`, or run the simulations first. If the replicate folders are in another place, name it with `data.local.yaml`, `polyzymd study locate DIR` or `--data`. |
| `error: Config file(s) not found` | Point `-c` at a simulation `config.yaml`. For a study, use `--study` with its `study.yaml`. |
| `error: No analysis named 'catalytic_triad'.` | The catalytic triad is a routine on the study API: follow {doc}`analysis_triad_quickstart`, or run `polyzymd analyze distances --set pairs=<pairs.yaml>` for the distances. |
| `polyzymd: command not found` | Run through `pixi run -e analysis`. |

The command exits with 0 on success. On each error above, it exits with 2 and
prints the message and the fix, one line each.

## See also

- {doc}`../explanation/analysis_entry_points` for which entry point to use.
- {doc}`analysis_compare_conditions` for comparing conditions step by step.
- {doc}`../reference/cli_reference` for every flag.
