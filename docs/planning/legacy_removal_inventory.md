# Legacy removal inventory

Every module under `src/polyzymd/analyses/` as it stands on `analyses_ported2`
at 677c86e4, with its line count, the modules and test files that import it,
and a verdict. The module is 26,722 lines across 58 files at this point. The
target after the removal is under 20,000.

This file is a plan, not part of the Sphinx documentation. It records what the
wave 4 branch `analyses/remove-legacy-framework` intends to delete and why, so
the plan can be checked before any code moves.

## How the importers were found

`ast` parses every `.py` file under `src/` and `tests/`, resolves relative
imports against the importing package, and records an edge whenever a module or
a name inside it is imported. Documentation and `.opencode` references were
found by literal search for the dotted module path. A second pass lists every
public top-level name whose only references outside its own module come from
tests, which is the grep equivalent of a `vulture` run. `vulture` itself is not
in any pixi environment here.

One caveat on that pass. The nine plugin classes and their `Settings` models
look test-only because `polyzymd.analyses.discovery` imports plugin modules by
name at runtime rather than through a static import. They are live.

## Verdicts

`keep` means the module stays roughly as it is. `trim` means the module stays
but loses named dead symbols. `collapse` means its content moves into another
module and the file disappears. `delete` means the whole file goes.

### Top level, 8,472 lines

| Module | Lines | Imported by | Verdict |
| --- | --- | --- | --- |
| `__init__.py` | 101 | `discovery`, 6 test files | keep, prune re-exports of deleted names |
| `base.py` | 630 | 14 source modules including `config/comparison.py` and `workflow/analysis_slurm.py`, 21 test files | collapse, step 3 folds the generated contract subclass in and drops every hook only `contract_runner` overrode |
| `catalytic_triad.py` | 163 | discovery at runtime, 3 test files | keep, contract plugin |
| `contacts/__init__.py` | 483 | discovery at runtime, 4 test files | keep, contract plugin |
| `contract.py` | 814 | 12 plugin and framework modules, 11 test files | keep, this is the surface the refactor is built on |
| `contract_plots.py` | 612 | `contract_runner`, `config/comparison.py` | keep |
| `contract_runner.py` | 622 | all nine plugins | collapse into `base.py` in step 3, widen the cache identity in step 4 |
| `discovery.py` | 285 | `__init__`, `orchestrator`, `protocols`, `cli/compare.py`, `cli/scaffold.py`, `config/comparison.py` | keep |
| `distances.py` | 163 | discovery at runtime, 3 test files | keep, contract plugin |
| `exceptions.py` | 124 | 23 source modules, 25 test files | keep |
| `hydrogen_bonds.py` | 484 | discovery at runtime, 2 test files | keep, contract plugin |
| `orchestrator.py` | 679 | `__init__`, `protocols`, `cli/compare.py`, `workflow/analysis_slurm.py` | keep |
| `protocols.py` | 1135 | `__init__`, `cli/analyze.py`, `cli/compare.py` | collapse to the one contract result shape in step 5, which removes `_read`, the legacy summary reader and the metric-key plumbing |
| `rg.py` | 288 | discovery at runtime, 6 test files | keep, contract plugin |
| `rmsd.py` | 376 | discovery at runtime, 3 test files | keep, contract plugin |
| `rmsf/__init__.py` | 273 | discovery at runtime, 3 test files | keep, contract plugin |
| `sasa.py` | 338 | discovery at runtime, 2 test files | keep, contract plugin |
| `secondary_structure.py` | 210 | discovery at runtime, 2 test files | keep, contract plugin |
| `stats.py` | 1448 | `_framework/compare`, `mda/comparison`, `protocols`, `shared/multi_run_formatting` | trim to `interpret_direction`, which `protocols` calls, and `format_pct`, which it uses for the same sentence. The other 1,300 lines serve only modules step 2 deletes |

### `_framework/`, 4,108 lines

| Module | Lines | Imported by | Verdict |
| --- | --- | --- | --- |
| `__init__.py` | 1 | package marker | keep |
| `aggregate_validation.py` | 477 | `_framework/io`, `_framework/lifecycle`, `base`, `mda/store` | trim, the contract lifecycle needs the fingerprint and the not-outdated check, not the per-field legacy validators |
| `cache_identity.py` | 478 | `_framework/aggregate_validation`, `_framework/lifecycle`, `contract_runner` | trim, `compute_cache_identity`, `extract_settings_fingerprint_from_path` and `validate_settings_fingerprint` have no caller outside tests. Step 4 rewrites the rest |
| `compare.py` | 152 | `base` only | delete, step 2. It is the default `compare()` for plugins that report `MetricValue` dictionaries, and no plugin does any more |
| `comparison_models.py` | 321 | `_framework/compare`, `base`, `contract_plots`, `mda/comparison` | delete most, step 2. Keep `BasePlotSettings`, which `ContractPlotSettings` subclasses, and `SlurmResourceHint`, which the orchestrator reads off a plugin class. `MetricValue`, `PairwiseResult`, `ANOVAResult`, `ComparisonResult`, `ConditionSummary`, `BaseConditionSummary` and `BaseComparisonResult` have no reader left once `compare.py` and `mda/comparison.py` go, and `ProtocolReport` builds its own models |
| `contexts.py` | 160 | 6 source modules | keep, minus the `MDA*Context` references step 3 removes |
| `contract.py` | 111 | `base` only | delete, step 3. It is the removed-hook police, it exists to raise when a subclass defines `run_replicate`, a method deleted two releases ago |
| `io.py` | 477 | `base` only | trim, the artifact readers stay, the per-plugin deserialize fallbacks go with the legacy result classes |
| `lifecycle.py` | 1919 | `orchestrator`, `protocols` | keep and trim. This is the runner. `AnalysisLifecycleAdapter` has no importer at all |
| `results_base.py` | 12 | none in `src`, one test | delete, it re-exports `polyzymd.utils.version.get_polyzymd_version` for backwards compatibility with nothing |

### `mda/`, 5,157 lines

| Module | Lines | Imported by | Verdict |
| --- | --- | --- | --- |
| `__init__.py` | 102 | 6 source modules, 22 test files | collapse, the re-export list loses every name step 2 and step 3 delete |
| `aggregation.py` | 659 | `mda/__init__`, `mda/comparison` | delete, step 3. It aggregates `AggregatedMetric` scalars off legacy `MDAAnalysisJob` collectors. `contract.aggregate_observables` replaced it |
| `artifacts.py` | 533 | 11 source modules, 10 test files | keep, the artifact envelopes are the on-disk format |
| `base.py` | 53 | 6 mda modules | keep |
| `comparison.py` | 773 | `_framework/compare`, `mda/__init__` | delete, step 2. It is the metric-dictionary comparison engine, superseded by `contract.compare_observables` |
| `frame_selection.py` | 581 | `contract`, 5 mda modules, 14 test files | keep |
| `job.py` | 539 | `_framework/lifecycle`, `contract_runner`, 3 mda modules | keep, `MDAAnalysisJob`, `MDABackendPolicy` and `MDAUniversePolicy` stay |
| `lifecycle.py` | 388 | `_framework/lifecycle`, `base`, `contract_runner`, `mda/__init__` | collapse, step 3 folds the parts the contract lifecycle needs into it and drops `build_mda_replicate_job_context`, which nothing imports |
| `pair_distance.py` | 181 | `catalytic_triad`, `distances`, `mda/__init__` | keep |
| `plugin.py` | 322 | `contract_runner`, `mda/__init__`, `mda/lifecycle` | collapse, step 3. `frame_selection_payload` and `MDACollectorContext` move to `mda/lifecycle.py`. `StrictJSONMDAResultCollector` and `strict_json_payload` exist to serialize an arbitrary plugin return value, and a contract plugin returns `Observable` objects |
| `store.py` | 522 | 7 source modules, 6 test files | keep |
| `universe.py` | 504 | `contract_runner`, `mda/__init__` | keep, `UniverseProvider` and `FileIdentity` carry the input identity |

### `shared/`, 8,229 lines

| Module | Lines | Imported by | Verdict |
| --- | --- | --- | --- |
| `__init__.py` | 168 | 31 references | keep |
| `aa_classification.py` | 279 | `sasa` | trim, `AAClass`, `get_aa_class`, `get_residues_for_class` and `get_selection_for_class` have no caller |
| `alignment.py` | 433 | `contract_runner`, `rmsd`, `rmsf`, `shared/centroid` | keep |
| `autocorrelation.py` | 723 | `contract`, `contract_runner`, `mda/aggregation`, `shared/__init__` | keep, the `mda/aggregation` edge goes with that file |
| `centroid.py` | 270 | `rmsd`, `shared/alignment` | trim, `get_reference_mode_description` has no caller |
| `diagnostics.py` | 384 | `mda/pair_distance`, `shared/centroid`, `shared/window` | trim, three of its public helpers have no caller |
| `inferential_statistics.py` | 793 | `contract`, `protocols`, `shared/multi_run_comparison`, `stats` | keep, trim the typing protocols nothing implements |
| `loader.py` | 1891 | 7 source modules | keep |
| `multi_run_comparison.py` | 239 | nothing in `src` | delete, step 2 |
| `multi_run_formatting.py` | 249 | `stats` only | delete, step 2, it formats the text report of the comparison engine that goes with it |
| `paths.py` | 59 | `_framework/io`, `_framework/lifecycle`, `cli/compare.py`, `workflow/analysis_slurm.py` | keep, trim `format_replicate_cache_token` |
| `plotting.py` | 1535 | `contract_plots`, `shared/__init__` | keep, trim `ArtifactPlotData` and the replicate-scatter helpers `contract_plots` does not call |
| `selections.py` | 226 | `mda/pair_distance` | keep, trim `translate_selection` and `ParsedSelection` |
| `statistics.py` | 468 | 10 source modules | keep |
| `topology.py` | 159 | `contacts`, `mda/universe`, `rg`, `shared/loader` | keep |
| `window.py` | 353 | `base`, `mda/frame_selection`, `shared/__init__` | keep |

## Outside the analyses package

The scaffold templates under `src/polyzymd/cli/_scaffold/templates/` are part of
the same removal. `simple_mda_plugin.py.jinja` (337), `advanced_plugin_init.py.jinja`
(241), `advanced_mda.py.jinja` (68), `test_simple_mda_plugin.py.jinja` (277) and
`test_advanced_plugin.py.jinja` (247) teach the deleted hooks. Step 6 deletes
them and makes the contract template the only style, which also removes the
`--advanced` flag and the `--style` choices from `polyzymd new-analysis`.

## What the deletions add up to

| Step | Modules | Lines removed, approximate |
| --- | --- | --- |
| 2, comparison engines | `_framework/compare`, `mda/comparison`, `shared/multi_run_comparison`, `shared/multi_run_formatting`, most of `_framework/comparison_models`, most of `stats` | 2,900 |
| 3, one lifecycle | `_framework/contract`, `_framework/results_base`, `mda/aggregation`, `mda/plugin`, `contract_runner` folded into `base`, hooks off `Analysis` | 1,900 |
| 4, cache identity | net small, the identity walk replaces the per-plugin version constants | 0 |
| 5, one result shape | `protocols` readers and their real-artifact fixtures | 500 |
| 6, one scaffold style | five templates and their tests | 1,200 source, of which 1,170 are templates |

That lands the analyses module near 21,000 lines before the trims, and the
named dead symbols account for the rest.

## Risks

Three modules carry the on-disk format and must not move in the same commit as
a behaviour change. `mda/artifacts.py` defines the artifact envelopes every
campaign tree already holds. `mda/store.py` reads and writes them.
`_framework/cache_identity.py` decides whether a stored replicate is still
valid, and step 4 changes it on purpose, so the parity run after step 4 is the
check that it still reproduces frozen numbers rather than silently recomputing
everything.

## What the removal actually did

Recorded on 2026-09-13, after the branch finished. The module went from 26,722
lines in 57 files to 20,569 lines in 47 files, a drop of 6,153. The 20,000
target was missed by 569 lines; the reasons are below.

| Step | Commit subject | Lines after | Change |
| --- | --- | --- | --- |
| 1 | inventory the legacy framework before removal | 26,722 | 0 |
| 2 | delete the scalar comparison engine | 23,619 | -3,103 |
| 3 | collapse the two analysis lifecycles | 22,077 | -1,542 |
| 4 | hash the code a cached replicate depends on | 22,179 | +102 |
| 5 | read one comparison shape in the report | 21,981 | -198 |
| 6 | delete public names with no caller | 21,479 | -502 |
| 7 | delete the job seam and the review's remaining findings | 20,569 | -910 |

The scaffold change the inventory lists as its own step landed inside step 2,
because the simple and advanced templates instantiate the result models step 2
deletes and no commit in between would have had a passing test suite.

### Where the remaining lines are

| Area | Lines | Why they stay |
| --- | --- | --- |
| `shared/` | 7,126 | `loader.py` (1,891) resolves and validates trajectory inputs, `plotting.py` (1,331) draws every figure the framework generates, `inferential_statistics.py` (586) and `autocorrelation.py` (714) are the cited statistics. All live. |
| top level | 6,741 | The nine plugins are 2,770 of it. `base.py` (1,024) is the lifecycle and the cache identity, `protocols.py` (960) the agent report, `contract.py` (876) the contract itself, `contract_plots.py` (612) and `orchestrator.py` (669) the figures and the scheduler. |
| `_framework/` | 3,221 | `lifecycle.py` (1,743) is the runner that walks conditions and replicates, resolves settings, writes caches and reports failures. `io.py` (477) and `aggregate_validation.py` (477) read and check stored aggregates. |
| `mda/` | 2,727 | The artifact envelopes, the store, frame selection, the universe provider and the replicate lifecycle. |

### Judged too risky to delete

The `AggregatedResultClass` and `ReplicateResultClass` hooks in
`_framework/io.py`. They are `None` for every plugin, but the generic reload
path still reads a plain JSON aggregate through them, which the lifecycle tests
rely on. Tightening the reader to accept only a `ConditionArtifact` is a
behaviour change to the stored-result contract, not a deletion, and it belongs
with a migration note.

`shared/inferential_statistics.py`'s remaining result models. `TTestResult`,
`EffectSize` and `BHResult` look unreferenced from outside the module, and they
are its own return annotations.

### Reconsidered after review

`mda/job.py` was kept in the first pass on the grounds that
`mda_backend_policy` in `comparison.yaml` is applied there. The review pointed
out that the key was already inert, because the contract lifecycle never passes
a policy to anything, so keeping 539 lines of job machinery preserved nothing
but the appearance of a feature. `MDAAnalysisJob`, `MDAFunctionAdapter` and
`MDABackendPolicy` are deleted, `MDAUniversePolicy` and `MDAJobResult` moved to
`mda/lifecycle.py`, and the settings key now parses, warns that it is ignored,
and will be rejected next release.
