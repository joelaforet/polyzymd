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

## Studies that predate the study folder

A comparison outside any study still works. It finds analyses in an
`analyses/` folder next to its `comparison.yaml`, and it resolves condition
paths relative to that file as before. Absolute condition paths still load,
but a study that uses them will not run on another machine.
