# Analysis module rules

## The real tree

Verified against `src/polyzymd/analyses/` on 2026-09-13. Every file named here
exists. If you add or delete a module, update this list in the same commit.

```
src/polyzymd/analyses/
├── contract.py          # Observable, the four kinds, aggregation, testing,
│                        # and contract_analysis(), which builds the Analysis
├── base.py              # The Analysis lifecycle and the framework contexts
├── contract_plots.py    # One figure per observable kind
├── discovery.py         # pkgutil auto-discovery of plugins
├── orchestrator.py      # Engine: compute, aggregate, compare, plot
├── protocols.py         # analyze() and ProtocolReport, the agent surface
├── stats.py             # interpret_direction, format_pct
├── exceptions.py        # Typed analysis errors
├── _framework/          # aggregate_validation, cache_identity,
│                        # comparison_models, contexts, io, lifecycle
├── mda/                 # artifacts, base, frame_selection, job, lifecycle,
│                        # pair_distance, store, universe
├── shared/              # aa_classification, alignment, autocorrelation,
│                        # centroid, diagnostics, inferential_statistics,
│                        # loader, paths, plotting, selections, statistics,
│                        # topology, window
├── catalytic_triad.py
├── contacts/__init__.py
├── distances.py
├── hydrogen_bonds.py
├── rg.py
├── rmsd.py
├── rmsf/__init__.py
├── sasa.py
└── secondary_structure.py
```

Every plugin is one module. There is no plugin package holding `_mda.py`,
`_plotters.py`, `_models.py`, `_formatters.py`, `_comparison.py` or
`_events.py`. Result models, persistence, statistics and figures belong to the
framework, not to a plugin.

## Loading trajectories

`polyzymd.analyses.shared.loader.TrajectoryLoader` is the canonical universe
loader. It resolves topology and trajectory files for a replicate, checks
segment lineage, and builds the MDAnalysis universe.

`polyzymd.analyses.mda.universe.UniverseProvider` wraps `TrajectoryLoader`. It
takes a `SimulationConfig`, instantiates the loader lazily, and adds input
provenance (`UniverseProvenance`) to each load. The replicate lifecycle uses
`UniverseProvider`. Nothing else builds a `Universe` directly, and no plugin
constructs file paths by hand.

## Public import surface

A plugin imports from `polyzymd.analyses.contract` (`Observable`,
`iter_frames`, `contract_analysis`) and from `polyzymd.analyses.shared` for a
helper that already exists. `polyzymd.analyses.base` holds `Analysis` and the
four framework contexts, which a plugin receives rather than constructs. Do not
import `polyzymd.analyses._framework` from a plugin.

## Adding a plugin

1. Run `polyzymd new-analysis <name>`. It writes one module and one test file.
2. Define `Settings` as a pydantic v2 `BaseModel`.
3. Define a plugin class with `name`, `Settings`, `references` and
   `compute(universe, frames, settings)` returning a sequence of `Observable`.
4. End the module with `NameAnalysis = contract_analysis(Name)`. Discovery
   finds that class; there is no registry to edit.
5. Never aggregate across replicates, run a test, write a file or import
   matplotlib inside `compute()`. The framework does all of it from `kind`.

## The four observable kinds

| kind | values are | replicate value | reported as |
|---|---|---|---|
| `mean_of_timeseries` | one number per frame | mean | mean, SEM, 95 percent CI across replicates |
| `fluctuation` | one number per frame | sample standard deviation | the same, on the fluctuation |
| `fraction` | 0 or 1 per frame | occupancy | the same, on the fraction |
| `profile` | one number per index | the whole vector | per-index mean and SEM |

A `profile` gets a comparable scalar through `reduce="mean_over_index"` or
`"sum_over_index"`. An observable that is a function of others the plugin
reports is declared `tested=False`, which keeps it out of the pairwise tests and
out of the Benjamini-Hochberg family.

## One lifecycle

`Analysis` is concrete. `contract_analysis()` builds one subclass per plugin,
differing only in `name`, `Settings` and `plugin`. The lifecycle writes the
replicate artifact with a framework-written identity block, reuses a replicate
whose identity still matches, aggregates by kind, compares under one
Benjamini-Hochberg family, plots per kind, and formats the result. There are no
`extract_metrics()`, `build_mda_collector()` or `compare()` hooks to override,
and no second comparison path.

## Results and plotting

Replicate, condition and comparison artifacts persist through `ArtifactStore`.
The per-frame series of every observable goes to an NPZ sidecar; a plugin that
also produces a raw table returns it as an extra sidecar from `compute()`. Do
not invent a plugin-specific cache filename scheme.

`plot()` reads cached artifacts and sidecars. It must not reload a trajectory or
rerun an analysis.

## Statistical contract

The replicate is the sampling unit for every cross-condition test and every
condition-level uncertainty. Equilibration is one global value applied
uniformly, and no diagnostic is allowed to select data. Every metric carries a
unit and a stated uncertainty. Invoke the `livecoms-check` skill in
`.claude/skills/` before committing analysis code.
