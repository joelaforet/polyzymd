---
name: polyzymd-analyze
description: Get a validated number out of PolyzyMD trajectories with ONE command, `polyzymd analyze <name> -c config.yaml [-c other.yaml]`. Use whenever the question is "what is the Rg / RMSF / SASA / contact count of this system" or "does condition A differ from condition B". Do not write your own MDAnalysis loop for an analysis that already exists.
---

# polyzymd-analyze, one command for one validated number

## 1. The command

```bash
pixi run -e analysis polyzymd analyze rg -c A/config.yaml -c B/config.yaml --eq 10ns
```

One `-c` gives a per-condition summary. Two or more give pairwise comparisons
with the first config as the control. `--format agent` is the default and
prints at most 25 lines. Use `--format json` when you need the full report, and
`-f comparison.yaml` when a comparison project already exists.

Environment: only the `analysis` and `sim-cuda-12-4` pixi envs have the CLI. The
bare `polyzymd` on PATH points at a system Python without click. Always go
through `pixi run -e analysis`.

List the nine analysis names with
`pixi run -e analysis polyzymd compare run --list`.

## 2. Reading the output

```
# polyzymd analyze rg  metric mean_rg  unit A  eq 10ns  conditions 2  replicates 3,3  protocol rg/1
A  n 3  mean 18.42  sem 0.05  ci95 18.2 to 18.64  values 18.4, 18.5, 18.36
B  n 3  mean 18.73  sem 0.06  ci95 18.47 to 18.99  values 18.71, 18.8, 18.68
A vs B  delta +0.31  ci95 0.02 to 0.6  p 0.041  p_adj 0.041  test student_t  correction BH  d 1.9  significant
verdict: B larger mean_rg than A (delta +0.31 A, 95% CI 0.02 to 0.6, p_adj 0.041, n 3 vs 3)
```

- `n` is the replicate count, which is the sample size for every test. Frames
  are never the sampling unit.
- `ci95` on a condition is the Student t interval on its mean. On a comparison
  it is the interval on the difference, uncorrected for multiplicity.
- `delta` is `mean(b) - mean(a)`, and `d` has the same sign.
- `significant` uses `p_adj` against the configured alpha, 0.05 by default.
- `run` in the header names the selection when an analysis measures several
  (rg does protein and polymer); `--run LABEL` picks another.

Verdict vocabulary, fixed so you can branch on it:

| Word | Means |
|---|---|
| `larger` / `smaller` | the second condition differs from the control after correction |
| `no significant difference` | the test ran and did not clear alpha; read the CI before calling it "the same" |
| `no test recorded` | the plugin stored no corrected p value; the line describes, it does not decide |
| `not testable` | a condition has fewer than two replicates, so no test exists |

Report the verdict sentence verbatim, with the unit and the replicate counts.

## 3. Rules

- Do not write your own MDAnalysis loop for an analysis that exists. Run the
  protocol and cite its provenance: the `analysis`, the `protocol_version`, the
  equilibration window, and the `config_hashes` from `--format json`.
- Do not read `_mda.py` to find out what an error bar means. The report says.
- Warnings are part of the answer. A `warning:` line about two replicates
  changes how the number should be read.
