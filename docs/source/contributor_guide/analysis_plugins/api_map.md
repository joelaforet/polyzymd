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

Use the public facades below as your default imports. They keep contributor
plugins insulated from internal framework changes while still covering the
common plugin lifecycle tasks.

### Create a plugin class

Import lifecycle classes and comparison vocabulary from
`polyzymd.analyses.base`.

```python
from polyzymd.analyses.base import Analysis, MetricValue
```

Use this surface when you need to:

- subclass `Analysis`;
- receive lifecycle contexts such as `AggregateContext`, `ComparisonContext`,
  or `PlotContext`;
- describe scalar metrics with `MetricValue`;
- return framework comparison models such as `ComparisonResult`.

See {doc}`../../api/analyses_base`.

### Build MDAnalysis-native replicate jobs

Import trajectory-native job helpers from `polyzymd.analyses.mda`.

```python
from polyzymd.analyses.mda import FrameSelection, MDAAnalysisJob
```

Use this surface when you need to:

- build replicate jobs around MDAnalysis `AnalysisBase`-compatible work;
- adapt simple trajectory functions into framework jobs;
- select frames with `FrameSelection`;
- receive MDAnalysis job context objects from the framework.

See {doc}`../../api/analyses_mda`.

### Emit artifacts and sidecars

Import artifact, collector, and sidecar helpers from `polyzymd.analyses.mda`.

```python
from polyzymd.analyses.mda import ArtifactStore, MDAArtifactCollector
from polyzymd.analyses.mda import MDACollectorContext, ReplicateArtifact
```

Use this surface when you need to:

- convert completed MDAnalysis jobs into canonical replicate artifacts;
- write and read artifact JSON and manifests;
- attach validated sidecars with `ArtifactSidecarRef`;
- use default replicate-to-condition aggregation;
- compare condition artifacts without inventing plugin-specific cache files.

See {doc}`../../api/analyses_mda`.

### Discover and test plugins

Import plugin discovery helpers from `polyzymd.analyses`.

```python
from polyzymd.analyses import get_analysis, list_analyses
```

Use this surface in tests, CLI helpers, and developer diagnostics when you need
to confirm that a plugin is auto-discovered or retrieve a plugin class by name.

See {doc}`../../api/analyses`.

### Compare scalar metrics

Import reusable scalar comparison and formatting helpers from
`polyzymd.analyses.stats`.

```python
from polyzymd.analyses.stats import default_scalar_comparison
from polyzymd.analyses.stats import format_scalar_comparison
```

Use this surface when a plugin reports one or more scalar `MetricValue` entries
and can use the framework's default t-tests, ANOVA, ranking, and CLI summary
formatting.

See {doc}`../../api/analyses`.

### Reuse shared analysis helpers

Use only the shared helpers documented on {doc}`../../api/analyses_shared`.
These modules are public contributor surfaces when imported from
`polyzymd.analyses.shared.*`.

For trajectory setup and analysis windows:

```python
from polyzymd.analyses.shared.loader import TrajectoryLoader
from polyzymd.analyses.shared.window import resolve_replicate_trajectory_window
```

For selections, selectors, and molecular group helpers:

```python
from polyzymd.analyses.shared import aa_classification, centroid, selections
from polyzymd.analyses.shared.selectors import base, polymer, protein, solvent
```

For statistics and convergence checks:

```python
from polyzymd.analyses.shared import autocorrelation, convergence
from polyzymd.analyses.shared import inferential_statistics, statistics
```

For plotting and multi-run summaries:

```python
from polyzymd.analyses.shared import multi_run_comparison
from polyzymd.analyses.shared import multi_run_formatting, plotting
```

### Start from scaffolding and configuration reference

Use the scaffold command when you want a working plugin skeleton and matching
tests before filling in analysis-specific logic.

```bash
pixi run -e build polyzymd new-analysis solvent_shell
```

Replace `solvent_shell` with the desired plugin name.

- See {ref}`CLI new-analysis reference <polyzymd-new-analysis>` for command
  options.
- See {doc}`../../reference/analysis_plugin_settings` for built-in
  `comparison.yaml` plugin settings under `plugins.<plugin_name>`.

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
