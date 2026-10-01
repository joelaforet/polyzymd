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
prints one line per condition and comparison, every one of them. Use
`--format json` when you need every field of the report.
`polyzymd analyze` does not read `comparison.yaml`: `-f comparison.yaml` exits 2
and prints the equivalent `-c ... --label ... --replicates ... --eq ...` command
built from the file, which is the command to run. `polyzymd compare` and
`polyzymd new-analysis` are retired: with any arguments they exit 2 and print
the replacement. An unknown analysis name exits 2 with the list of names.

Environment: only the `analysis` and `sim-cuda-12-4` pixi envs have the CLI. The
bare `polyzymd` on PATH points at a system Python without click. Always go
through `pixi run -e analysis`.

`pixi run -e analysis polyzymd analyze --help` lists the analysis names: rg,
rmsd, rmsf, rmsd_per_residue, sasa, secondary_structure, contacts, native_contacts,
hydrogen_bonds and distances, the keys of
`polyzymd.analyses.protocols.FUNCTION_ANALYSES`. Every one reads `-c config.yaml`.
For contacts, `--run mean_lifetime` reports how long contacts last; for
hydrogen_bonds, `--run protein_polymer_mean_lifetime`, `protein_polymer_residues`
and `protein_polymer_pairs` report how long the bonds last and how often each
residue and residue pair is bonded.

`rmsd` is one number per frame; `rmsd_per_residue`, `rmsf` and `offset` are one
number per residue over all frames, measured from the reference, from each
atom's own mean, and as the mean's distance from the reference, with
`rmsd_per_residue² = rmsf² + offset²`. Report `rmsf` for flexibility and
`offset` for distance from a reference state; read the "Fluctuation, offset
and deviation" section of `docs/source/explanation/analysis_rmsf_best_practices.md`
before choosing.

The catalytic triad is not an analysis name: `polyzymd analyze catalytic_triad`
exits 2. It is a routine on the study API that counts each triad hydrogen bond
with `functions.hbond_count` and combines them with `Timeseries.transform`:
follow `docs/source/how_to/analysis_triad_quickstart.md`
(https://polyzymd.readthedocs.io/en/latest/how_to/analysis_triad_quickstart.html),
and run `polyzymd analyze distances -c config.yaml --set pairs=pairs.yaml` for
the triad distances. For any other question, write a function of an MDAnalysis
`Universe` and run it with `Study.timeseries` or `Study.per_replicate`; see
`docs/source/reference/analysis_functions.md` and
`docs/source/how_to/analysis_agent_protocol.md`.

Many replicates: add `--submit --preset <cluster>` (`alpine-cpu`,
`blanca-shirts`, `blanca-chbe-rdi`, `bridges2-rm`; `--partition`, `--account`,
`--qos`, `--time`, `--mem`, `--cpus` override) to the same command. It submits
one SLURM array task per condition and replicate and a report job that starts
after them all (`afterany`), from `<output-dir>/slurm/<name>_<time>/`, and the
report lands in `report.txt` there (`report.json` with `--format json`, or
`-o PATH`). `--dry-run` writes the
scripts without submitting. Load the cluster's SLURM module first, such as
`module load slurm/blanca`. See `docs/source/how_to/hpc_execution.md`.

A study folder with a `study.yaml` needs no `-c` list. Run
`pixi run -e analysis polyzymd study check STUDY` first: it reads no
trajectory and prints each condition's runs and each analysis run, with
whether it has stored results. Then `polyzymd analyze RUN --study STUDY`
runs one entry of its `analyses:`, and `polyzymd analyze --study STUDY` runs
them all; command-line options override the file, and results go to
`STUDY/results/RUN/`. To read results back without trajectories, use
`pz.Study("STUDY/study.yaml").results(RUN).table`. See
`docs/source/how_to/study_yaml.md`. When `study check` finds no runs for a
condition, the trajectories are elsewhere: run `polyzymd study locate DIR`
(writes the gitignored `data.local.yaml`) or pass `--data DIR`; never edit the
configs' paths to point at moved data. Read `study check`'s production lengths
before choosing `--eq`; if a report warns that conditions were analysed up to
different times, rerun with `--until <shortest>` before comparing them. The
console shows only reports and warnings; the full log is the `log:` path, so
read it only when a run fails unexpectedly. To publish, fill `metadata:` and run
`polyzymd study freeze STUDY`; its `warning:` lines list what is missing or stale
(`docs/source/how_to/study_freeze.md`). Then hand the author `deposit/UPLOAD.md`:
uploading and publishing on Zenodo are theirs, never an agent's.

## 2. Reading the output

```
# polyzymd analyze rg  metric mean_rg  unit A  eq 10ns  conditions 2  replicates 3,3  protocol rg/2
A  n 3  mean 18.42  sem 0.05  ci95 18.2 to 18.64  values 18.4, 18.5, 18.36
B  n 3  mean 18.73  sem 0.06  ci95 18.47 to 18.99  values 18.71, 18.8, 18.68
A vs B  delta +0.31  ci95 0.02 to 0.6  p 0.041  p_adj 0.041  test welch_t  correction BH  d 1.9  significant
verdict: B larger mean_rg than A (delta +0.31 A, 95% CI 0.02 to 0.6, p_adj 0.041, n 3 vs 3)
```

- `n` is the replicate count, which is the sample size for every test. Frames
  are never the sampling unit.
- `ci95` on a condition is the Student t interval on its mean. On a comparison
  it is the interval on the difference, uncorrected for multiplicity.
- `delta` is `mean(b) - mean(a)`, and `d` has the same sign.
- `significant` means `p_adj` is at most 0.05 (Benjamini-Hochberg over the
  comparisons of the report).
- `run` in the header names the selection when an analysis measures several
  (rg does protein and polymer); `--run LABEL` picks another.

Verdict vocabulary, fixed so you can branch on it:

| Word | Means |
|---|---|
| `larger` / `smaller` | the second condition differs from the control after correction |
| `no significant difference` | the test ran and did not clear alpha; read the CI before calling it "the same" |
| `no test recorded` | the row has no corrected p value; the line describes, it does not decide |
| `not testable` | a condition has fewer than two replicates, so no test exists |

Report the verdict sentence verbatim, with the unit and the replicate counts.

## 3. Rules

- Do not write your own MDAnalysis loop for an analysis that exists. Run the
  protocol and cite its provenance: the `analysis`, the `protocol_version`, the
  equilibration window, and the `config_hashes` from `--format json`.
- Do not read the source to find out what an error bar means. The report says,
  and `docs/source/reference/analysis_protocol_report.md` defines every field.
- Warnings are part of the answer. A `warning:` line about two replicates
  changes how the number should be read.
