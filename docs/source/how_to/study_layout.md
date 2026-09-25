# Organize a study

A study is one folder that holds every condition, every comparison and every
analysis of a piece of work. You open that one folder to work on the study, and
when it is finished you zip it or push it to GitHub and publish it with the
paper. The trajectories are archived separately, for example on Zenodo, and
anyone with the folder can regenerate them within the precision of molecular
dynamics.

## Create the study

```bash
polyzymd study init -n lipase_thermal_stability
cd lipase_thermal_stability
```

```text
lipase_thermal_stability/
├── study.yaml      name of the study and where its analyses live
├── README.md       this layout, for the next person or agent
├── conditions/     one folder per simulated condition
├── comparisons/    one folder per comparison
├── analyses/       analyses written for this study
├── structures/     inputs shared across conditions
└── workflows/      scripts that make the paper's tables and figures
```

`study.yaml` marks the root. Every PolyzyMD command run anywhere below it finds
the study by walking up to that file.

## Add conditions and comparisons

Each condition is its own `polyzymd init` project, and each comparison is its
own `polyzymd compare init` project:

```bash
polyzymd init -n conditions/CALB_noPoly_343K
polyzymd init -n conditions/CALB_SBMA_EGMA_50_50_343K
polyzymd compare init -n CALB_343K -o comparisons
```

A paper usually holds several comparisons, for example one per enzyme and
temperature. Keeping each in its own folder keeps its `analysis/`,
`comparison/` and `figures/` results apart from the others.

List conditions in `comparisons/CALB_343K/comparison.yaml` by paths relative to
that file, so the study still works after it is copied to another machine:

```yaml
conditions:
  - label: "No Polymer"
    config: "../../conditions/CALB_noPoly_343K/config.yaml"
    replicates: [1, 2, 3]
  - label: "SBMA-EGMA 50:50"
    config: "../../conditions/CALB_SBMA_EGMA_50_50_343K/config.yaml"
    replicates: [1, 2, 3]
```

## Write an analysis for the study

Run `polyzymd new-analysis` anywhere inside the study:

```bash
polyzymd new-analysis lid_opening
```

It writes `analyses/lid_opening.py` and `analyses/test_lid_opening.py`. The
plugin is a settings model and a `compute()` that returns observables; the
{doc}`contributor guide <../contributor_guide/analysis_plugins/index>` explains
the contract. Run its tests from the study root with `pytest analyses/`.

To run it, name it under `plugins:` in a comparison, exactly like a built-in
analysis:

```yaml
plugins:
  lid_opening:
    selection: "resid 139-147"
```

```bash
cd comparisons/CALB_343K
polyzymd compare run lid_opening
```

Every comparison in the study sees every analysis in `analyses/`, so an
analysis is written once, not copied into each condition. A study analysis may
not reuse the name of a built-in one, so `rmsf` in any published comparison
means the same code.

Files in `analyses/` whose names start with `_` or `test_`, and a `tests/`
folder, are not loaded as analyses. Put shared helpers in a file such as
`analyses/_geometry.py` and import them with `from ._geometry import angle`.

## Collect every number

`polyzymd study results`, run anywhere inside the study, writes three tables to
`results/`:

- `conditions.csv`: one row per comparison, analysis, condition and scalar
  observable, with the mean, SEM, 95 percent interval over replicates, the
  replicate values and their count.
- `comparisons.csv`: one row per pairwise test, with the difference, the test,
  the adjusted p-value and the effect size.
- `profiles.csv`: one row per residue or bin of a profile observable.

The same tables are available in Python, which is how a figure script in
`workflows/` should read them:

```python
from polyzymd.analyses import load_results

results = load_results(".")                      # the study root
temperatures = results.to_dataframe("conditions")
```

Only the comparison results are read, so this also works on a published study
without its trajectories.

## Look at a study before it finishes

Nothing waits for every replicate. If a comparison lists replicates 1 to 5 and
only 1 to 3 have trajectories, `polyzymd compare run` computes the condition
from those three. A replicate still running contributes the segments that have
finished. Every result records what it used: the replicates listed and the
ones used, why any were left out, and for each replicate the frames read, the
simulated time they cover and that time as a fraction of the planned
production length.

A result from partial data says so wherever it is read:

- the text report and the agent report start with `PARTIAL:` lines, for
  example `PARTIAL: SBMA: replicates 1-3 of 1-5; replicate 2 at 64% of 100 ns (still running)`;
- the figure footnote ends with the same statement;
- rows of `polyzymd study results` have `complete` false, and condition rows
  name the replicates listed and the smallest production fraction reached;
- `polyzymd study export` names every partial result it packages.

To look at a replicate while it is still being written, for example to see
whether the protein is unfolding, also read the segment in progress:

```bash
polyzymd analyze rmsd -c conditions/CALB_noPoly_343K/config.yaml --replicates 2 \
    --eq 10ns --include-running
```

Run it again later and the replicate is recomputed with the new frames, since
the growing trajectory no longer matches the cached result.

## Work outside the framework

The framework's statistics apply to analyses written against the contract. For
anything else, write a script in `workflows/` and read the data yourself.
`load_replicate` returns the universe and the production window an analysis
would receive, with restart segments joined and checked:

```python
from polyzymd.analyses import iter_frames, load_replicate

universe, frames = load_replicate(
    "conditions/CALB_noPoly_343K/config.yaml", 1, equilibration="10ns"
)
protein = universe.select_atoms("protein")
rg = [protein.radius_of_gyration() for _ in iter_frames(universe, frames)]
```

`frames.run_kwargs()` gives the same window for an MDAnalysis
`AnalysisBase.run()` call. Scripts in `workflows/` may also read the JSON
results under each comparison and style figures however the paper needs.

## Share or move the study

Every replicate result records the files it was computed from by their path
relative to the replicate's working directory, their size and a content
fingerprint, and the config hash leaves out where the study lives. So a study
copied to another machine, unzipped or cloned from git keeps its cached
results: a file whose modification time changed still matches when its content
does, and an edited trajectory does not.

A published study usually arrives without its trajectories. Its comparisons
still aggregate, compare and plot from the replicate results in `analysis/`.
The inputs those results were computed from cannot be checked, so each
aggregate records a warning saying so rather than refusing. Recomputing a
replicate needs the trajectories again, downloaded to the location each
condition's `output.scratch_directory` names.

Run `polyzymd compare validate` in each comparison before publishing. It warns
about any condition named by an absolute path, which stops the study working
once it moves.

To package the study, run `polyzymd study export` anywhere inside it. It writes
`<study>.zip` next to the study folder. Only condition folders that some
`comparison.yaml` lists are packaged, so reruns and tests that sit beside them
stay out, and the command names each one it left out. Trajectories,
checkpoints, caches and SLURM logs are never packaged. The zip holds
`bundle_manifest.json`, which records the SHA-256 of every packaged file and
every trajectory the published results were computed from, with its size and
content fingerprint; that list is what to deposit on Zenodo.

Whoever receives the zip unpacks it and runs `polyzymd study verify <study>`.
Every packaged file is checked, and every listed trajectory that has been
downloaded to where its condition's config expects it is checked by content.

## Studies that predate the study folder

A comparison outside any study still works. It finds analyses in an
`analyses/` folder next to its `comparison.yaml`, and it resolves condition
paths relative to that file as before. Absolute condition paths still load,
but `compare validate` warns about them.
