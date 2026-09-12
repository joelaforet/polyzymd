# RMSF plugin reference

For a step-by-step guide to running RMSF analysis, see
{doc}`../how_to/analysis_rmsf_quickstart`.

## Settings

RMSF settings live under `plugins.rmsf` in `comparison.yaml`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `selection` | `str` | `"protein and name CA"` | MDAnalysis selection whose residues carry the profile |
| `alignment_selection` | `str` | `"protein and name CA"` | MDAnalysis selection superposed before the fluctuation is measured |
| `centroid_selection` | `str` | `"protein"` | MDAnalysis selection used to pick the representative frame in centroid mode |
| `reference_mode` | `str` | `"centroid"` | Alignment reference: `centroid`, `average`, `frame` or `external` |
| `reference_frame` | `int \| null` | `null` | One-indexed frame, required when `reference_mode: frame` |
| `reference_file` | `str \| null` | `null` | Structure file, required when `reference_mode: external` and it must exist |

### Minimal plugin block

```yaml
plugins:
  rmsf:
    selection: "protein and name CA"
    reference_mode: "centroid"
```

### External reference example

```yaml
plugins:
  rmsf:
    selection: "protein and name CA"
    reference_mode: "external"
    reference_file: "/path/to/crystal_structure.pdb"
```

## Observables

Every frame of the production window is used. The unit is angstrom throughout.

| Name | Kind | Reported per replicate | When |
|---|---|---|---|
| `rmsf` | `profile` | Per-residue fluctuation about the mean structure of the aligned window, indexed by residue ID | always |
| `rmsf_mean` | `fluctuation` | Mean of `rmsf` over residues | always |
| `rmsd_about_reference_per_residue` | `profile` | Per-residue deviation from the external reference structure, indexed by residue ID | `reference_mode: external` |
| `rmsd_about_reference_per_residue_mean` | `mean_of_timeseries` | Mean of that profile over residues | `reference_mode: external` |

The framework reduces each observable by kind. A profile keeps its per-index
vector, and the condition artifact reports a per-residue mean and SEM across
replicates. A `fluctuation` gives a mean, a SEM and a Student t 95 percent
interval across replicates, and it is the quantity the cross-condition test
uses. Profiles are not tested pairwise.

In external mode the trajectory is superposed on a structure that the
simulation did not produce, so the deviation from that structure is not a
fluctuation about the trajectory mean. The two are reported as separate
observables rather than one number whose meaning depends on a setting, and only
`rmsf_mean` is labelled a fluctuation.

The external reference file is part of the replicate identity, so replacing it
in place recomputes rather than reusing the cached profile.

`reference_mode: centroid` and `reference_mode: average` build their reference
from a contiguous slice of the trajectory, so the framework may not hand them a
non-uniform list of frames; doing so raises `PluginContractError`. Use
`reference_mode: frame` or `external` with such a window.

## Output files

```text
<comparison_workspace>/
├── analysis/
│   └── <condition>/
│       └── rmsf/
│           ├── run_1/
│           │   ├── result.json
│           │   └── observables.npz
│           ├── run_2/ ...
│           └── aggregated/
│               └── result.json
└── comparison/
    └── rmsf/
        └── result.json
```

| Level | Artifact | Path |
|-------|----------|------|
| Per replicate | `ReplicateArtifact` | `analysis/<condition>/rmsf/run_<replicate>/result.json` |
| Per condition | `ConditionArtifact` | `analysis/<condition>/rmsf/aggregated/result.json` |
| Cross condition | `ComparisonArtifact` | `comparison/rmsf/result.json` |
| Per-residue arrays | NPZ sidecar | `analysis/<condition>/rmsf/run_<replicate>/observables.npz` |

`payload["observables"]` is a list of observable records: `ObservableEstimate`
in a replicate artifact, `ObservableAggregate` in a condition artifact. The
comparison artifact adds `payload["comparisons"]`, one `ObservableComparison`
per observable and non-control condition. The field names are those of
{doc}`../api/analyses`.

```python
from pathlib import Path

from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.mda import ArtifactStore

store = ArtifactStore(Path("analysis/PEGylated/rmsf/aggregated"))
for payload in store.read_condition_result().payload["observables"]:
    aggregate = ObservableAggregate.model_validate(payload)
    print(aggregate.name, aggregate.kind, aggregate.mean, aggregate.unit)
```

## Figures

RMSF has no plotter of its own. Figures come from the framework, keyed on
observable kind, and are not yet implemented, so `compare run --plot` produces
no RMSF figures. A `plot_settings.rmsf` block in an existing `comparison.yaml`
still loads and raises a `DeprecationWarning` saying it does nothing.

The optional secondary-structure annotation bar under the profile plot was
dropped with the old plotter. It duplicated what the `secondary_structure`
analysis reports and needed an external reference file to exist.

## Common CLI options

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Path to comparison configuration |
| `--eq-time` | `0ns` | Equilibration time to skip |
| `--recompute` | off | Ignore cached results and recompute |
| `--format` | `table` | Output format (`table` or `json`) |
| `-o, --output` | (none) | Save formatted output to file |
| `-q, --quiet` | off | Suppress INFO messages |
| `--debug` | off | Enable DEBUG logging |

## Troubleshooting

### `rmsf: selection ... matched no atoms`

The MDAnalysis selection matches nothing in the topology. Check residue
numbering and atom names, and start from `selection: "protein and name CA"`.

### `reference_file does not exist`

`reference_mode: external` is set and the path is wrong. Use an absolute path
or one relative to the working directory.

### `rmsf: external reference ... gives N atoms over M residues`

The selection resolves to a different atom or residue set in the trajectory and
in the reference. Use a reference holding the same selected residues in the
same order.

### Very high values, above 10 angstrom

Usually an alignment mismatch, an overly broad selection, or genuine
instability. Check `alignment_selection` against `selection`, and cross-check
with `reference_mode: "average"`.

For interpretation guidance, see
{doc}`../explanation/analysis_rmsf_best_practices`.
