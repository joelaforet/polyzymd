# API map for analysis plugin contributors

This reference page maps common contributor needs to the stable import surface
and the nearest detailed documentation. Use it as a lookup bridge when writing or
reviewing analysis plugins.

```{warning}
`polyzymd.analyses._framework` is private implementation infrastructure, not a
contributor import surface. Leading-underscore helper modules inside built-in
plugin packages are examples of package organization only; they are not public
APIs for contributor plugins.
```

## Public API map

| Need | Import from | Docs link |
|------|-------------|-----------|
| Subclass the plugin base class, receive lifecycle contexts, describe scalar metrics, or return comparison models | `polyzymd.analyses.base` | {doc}`../../api/analyses_base` |
| Build MDAnalysis-native replicate jobs and function adapters | `polyzymd.analyses.mda` | {doc}`../../api/analyses_mda` |
| Select frames for trajectory work | `polyzymd.analyses.mda.FrameSelection` | {doc}`../../api/analyses_mda` |
| Convert completed MDAnalysis jobs into canonical replicate artifacts | `polyzymd.analyses.mda.MDACollectorContext`, `polyzymd.analyses.mda.MDAArtifactCollector`, `polyzymd.analyses.mda.ReplicateArtifact` | {doc}`../../api/analyses_mda` |
| Write and read artifact JSON, manifests, and validated sidecars | `polyzymd.analyses.mda.ArtifactStore`, `polyzymd.analyses.mda.ArtifactSidecarRef` | {doc}`../../api/analyses_mda` |
| Use default replicate-to-condition aggregation or condition-artifact comparison | `polyzymd.analyses.mda.aggregate_replicate_artifacts`, `polyzymd.analyses.mda.compare_condition_artifacts` | {doc}`../../api/analyses_mda` |
| Discover plugin classes by name in tests, CLI helpers, or developer diagnostics | `polyzymd.analyses.get_analysis`, `polyzymd.analyses.list_analyses`; detailed helpers in `polyzymd.analyses.discovery` | {doc}`../../api/analyses` |
| Use reusable scalar comparison and CLI formatting helpers | `polyzymd.analyses.stats` | {doc}`../../api/analyses` |
| Locate trajectories and resolve analysis windows | `polyzymd.analyses.shared.loader`, `polyzymd.analyses.shared.window` | {doc}`../../api/analyses_shared` |
| Use documented autocorrelation, summary-statistic, inferential-test, or convergence utilities | `polyzymd.analyses.shared.autocorrelation`, `polyzymd.analyses.shared.statistics`, `polyzymd.analyses.shared.inferential_statistics`, `polyzymd.analyses.shared.convergence` | {doc}`../../api/analyses_shared` |
| Apply shared plotting themes, figure paths, axis helpers, legends, grouped bars, or matrix annotations | `polyzymd.analyses.shared.plotting` | {doc}`../../api/analyses_shared` |
| Reuse documented selection parsing, selector implementations, residue-classification, or centroid helpers | `polyzymd.analyses.shared.selections`; selector submodules `polyzymd.analyses.shared.selectors.base`, `polyzymd.analyses.shared.selectors.protein`, `polyzymd.analyses.shared.selectors.polymer`, `polyzymd.analyses.shared.selectors.solvent`; `polyzymd.analyses.shared.aa_classification`; `polyzymd.analyses.shared.centroid` | {doc}`../../api/analyses_shared` |
| Use multi-run comparison or formatting helpers for analyses with several named runs per condition | `polyzymd.analyses.shared.multi_run_comparison`, `polyzymd.analyses.shared.multi_run_formatting` | {doc}`../../api/analyses_shared` |
| Scaffold a new analysis plugin and matching tests | Command: `pixi run -e build polyzymd new-analysis <name>` | {ref}`CLI new-analysis reference <polyzymd-new-analysis>` |
| Look up built-in plugin settings in `comparison.yaml` | `plugins.<plugin_name>` settings keys in YAML | {doc}`../../reference/analysis_plugin_settings` |

## Stable import boundaries

- Prefer `polyzymd.analyses.base` for plugin lifecycle classes and comparison
  result vocabulary.
- Prefer `polyzymd.analyses.mda` for trajectory-native jobs, collectors,
  artifacts, sidecars, stores, aggregation, and comparison.
- Use only utilities documented on {doc}`../../api/analyses_shared` from
  `polyzymd.analyses.shared`.
- Treat built-in plugin package roots, such as `polyzymd.analyses.contacts`, as
  public locations for their documented analysis classes and settings. Treat
  their leading-underscore helper modules as implementation examples, not
  contributor dependencies.

## Related contributor pages

- {doc}`architecture` explains why these public surfaces exist.
- {doc}`first_scaffold` shows the scaffold command in a guided first pass.
- {doc}`simple_scalar_plugin` shows a small MDAnalysis-native plugin using the
  public base and MDAnalysis facades.
- {doc}`sidecars` shows artifact-owned sidecar storage patterns.
- {doc}`checklist` summarizes pre-review API-boundary checks.
