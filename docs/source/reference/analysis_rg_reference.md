# Rg Plugin Reference

For a step-by-step guide to running Rg analysis, see
{doc}`../how_to/analysis_rg_quickstart`.

## Configuration Reference

All fields for `RgRunSettings`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Human-readable run label (must be unique) |
| `selection` | `str` | *required* | MDAnalysis selection for Rg calculation |
| `calculation_mode` | `str` | `"selection"` | `"selection"` for whole-group Rg, `"fragments"` for per-fragment reduction |
| `fragment_weighting` | `str` | `"equal"` | `"equal"` (arithmetic mean) or `"mass"` (mass-weighted mean). Only valid when `calculation_mode="fragments"` |
| `save_fragment_distribution` | `bool` | `true` | Save per-fragment Rg values in NPZ sidecar for distribution analysis |
| `histogram_bins` | `int` | `50` | Number of bins for fragment/reduced distribution histograms (minimum 2) |

Top-level `RgSettings` contains a single field:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | `list[RgRunSettings]` | *required* | One or more named Rg runs (at least one required) |

```{note}
Run labels must be unique within a single `comparison.yaml`. Duplicate labels
raise a validation error.
```

```{warning}
Unlike RMSD (which defaults `selection` to `"protein and name CA"`), Rg has
**no default selection**. You must always specify `selection` explicitly for
each run.
```

## Fragment Mode Reference

```{versionadded} 1.3.0
Fragment-aware Rg calculation was added in PolyzyMD 1.3.0.
```

Standard selection mode computes one Rg value per frame for the full atom
group matched by `selection`. Fragment mode first computes Rg for each
disconnected topological fragment, then reduces those per-fragment values to
one per-frame value.

Use fragment mode when your selection includes many independent molecules and
you care about the average fragment size rather than the size of the entire
multi-molecule cloud.

### Fragment mode configuration

```yaml
plugins:
  rg:
    runs:
      - label: protein_rg
        selection: protein

      - label: polymer_blob_rg
        selection: "resname SBM or resname EGM or resname EGP"
        calculation_mode: fragments
        fragment_weighting: equal
```

| Setting | Meaning |
|---------|---------|
| `calculation_mode: "selection"` | Whole-group Rg |
| `calculation_mode: "fragments"` | Per-fragment Rg with reduction |
| `fragment_weighting: "equal"` | Arithmetic mean over fragments |
| `fragment_weighting: "mass"` | Mass-weighted mean over fragments |

### How fragment mode reduction works

1. Identify disconnected fragments within the selected atom group
2. Compute fragment-level Rg for each frame
3. Reduce fragment values to one per-frame value using equal or mass weighting
4. Use the reduced timeseries for summary statistics and comparisons
5. Optionally save pooled fragment values in NPZ sidecars for distribution plots

## Why Rg Has No Alignment or Reference Fields

Rg is based on mass-weighted distances from the center of mass, so it is
translation and rotation invariant.

Because of that, Rg runs do not use RMSD-style fields such as
`alignment_selection`, `reference_mode`, `reference_file`, or
`reference_frame`.

## Output Files

Results are saved as canonical v1.3 artifacts. JSON files use framework-owned
artifact envelopes, and per-frame or distribution arrays are stored as NPZ
sidecars.

```text
<comparison_workspace>/
├── analysis/
│   └── <condition>/
│       └── rg/
│           ├── run_1/
│           │   ├── result.json
│           │   └── sidecars/
│           │       ├── rg_protein_rg_timeseries.npz
│           │       └── rg_polymer_blob_rg_timeseries.npz
│           ├── run_2/
│           │   └── ...
│           ├── run_3/
│           │   └── ...
│           └── aggregated/
│               ├── result.json
│               └── sidecars/
│                   └── rg_polymer_blob_rg_distribution.npz
└── comparison/
    └── rg/
        └── result.json
```

The canonical paths are:

| Level | Artifact | Path |
|-------|----------|------|
| Per replicate | `ReplicateArtifact` | `analysis/<condition>/rg/run_<replicate>/result.json` |
| Per condition | `ConditionArtifact` | `analysis/<condition>/rg/aggregated/result.json` |
| Cross condition | Comparison result | `comparison/rg/result.json` |
| Large arrays | NPZ sidecars | `analysis/<condition>/rg/**/sidecars/*.npz` |

Each replicate artifact contains JSON summaries for all configured runs. NPZ
sidecars store per-frame Rg timeseries and optional fragment distributions.

### Artifact envelope fields

| Field | Description |
|-------|-------------|
| `payload` | Rg run summaries, scalar metrics, fragment statistics, and relative sidecar paths |
| `metadata` | Run settings, calculation modes, equilibration labels, and units |
| `provenance` | Input topology/trajectory identity and workflow details |
| `sidecars` | Validated references to `sidecars/*.npz` arrays with hashes and sizes |

Use `ArtifactStore` for programmatic access:

```python
from pathlib import Path

from polyzymd.analyses.mda import ArtifactStore

replicate = ArtifactStore(Path("analysis/PEGylated/rg/run_1")).read_replicate_result()
condition = ArtifactStore(Path("analysis/PEGylated/rg/aggregated")).read_condition_result()
print(replicate.payload["runs"][0]["mean_rg"])
print(condition.payload["runs"][0]["metrics"]["mean_rg"])
```

### NPZ sidecar arrays

Each `rg_<label>_timeseries.npz` may include:

| Array | Mode | Description |
|-------|------|-------------|
| `rg_values` | always | Per-frame reduced Rg timeseries (Å) |
| `time_ns` | always | Time axis in ns |
| `frames` | always | 0-indexed frame indices |
| `fragment_rg_values` | fragments only | Pooled fragment-level Rg values across all frames |
| `fragment_counts_per_frame` | fragments only | Number of fragments detected per frame |
| `fragment_masses` | fragments + mass weighting | Fragment masses used for weighted reduction |

### JSON result structures

Per-replicate result (`ReplicateArtifact`), representative structure:

```python
{
    "schema_version": "1",
    "artifact_type": "replicate",
    "analysis_name": "rg",
    "condition_label": "PEGylated",
    "replicate": 1,
    "payload": {
        "runs": [
            {
                "run_label": "protein_rg",
                "selection": "protein",
                "calculation_mode": "selection",
                "mean_rg": 18.234,
                "sem_rg": 0.098,
                "timeseries_sidecar": "sidecars/rg_protein_rg_timeseries.npz"
            },
            {
                "run_label": "polymer_blob_rg",
                "selection": "resname SBM or resname EGM or resname EGP",
                "calculation_mode": "fragments",
                "fragment_weighting": "equal",
                "mean_rg": 8.412,
                "sem_rg": 0.054,
                "mean_fragments_per_frame": 50.0,
                "timeseries_sidecar": "sidecars/rg_polymer_blob_rg_timeseries.npz"
            }
        ]
    },
    "metadata": {"equilibration": "10ns", "time_unit": "ns"},
    "provenance": {"trajectory_files": ["..."], "n_frames_used": 9000},
    "sidecars": [
        {"path": "sidecars/rg_protein_rg_timeseries.npz", "metadata": {"kind": "timeseries"}},
        {"path": "sidecars/rg_polymer_blob_rg_timeseries.npz", "metadata": {"kind": "timeseries"}}
    ]
}
```

Aggregated result (`ConditionArtifact`), representative structure:

```python
{
    "schema_version": "1",
    "artifact_type": "condition",
    "analysis_name": "rg",
    "condition_label": "PEGylated",
    "replicates": [1, 2, 3],
    "payload": {
        "runs": [
            {
                "run_label": "protein_rg",
                "calculation_mode": "selection",
                "metrics": {
                    "mean_rg": {"values": [18.234, 18.291, 18.244], "mean": 18.256, "sem": 0.044}
                },
                "distribution_sidecar": "sidecars/rg_protein_rg_distribution.npz"
            },
            {
                "run_label": "polymer_blob_rg",
                "calculation_mode": "fragments",
                "metrics": {
                    "mean_rg": {"values": [8.412, 8.445, 8.439], "mean": 8.432, "sem": 0.021}
                },
                "fragment_summary": {"overall_mean_fragments_per_frame": 50.0},
                "distribution_sidecar": "sidecars/rg_polymer_blob_rg_distribution.npz"
            }
        ]
    },
    "metadata": {"equilibration": "10ns"},
    "provenance": {"source_replicates": [1, 2, 3]},
    "sidecars": [
        {"path": "sidecars/rg_polymer_blob_rg_distribution.npz", "metadata": {"kind": "distribution"}}
    ]
}
```

## Plot Types

The Rg plugin generates figures via `polyzymd compare plot-all`.

| Plot output | Description |
|-------------|-------------|
| `rg_timeseries_<run>.png` | Mean Rg vs time with SEM shading, one figure per run |
| `rg_comparison_<run>.png` | Grouped bar chart of mean Rg across conditions, one figure per run |
| `rg_distribution_<run>.png` | Distribution view. Selection mode: reduced distribution panel only. Fragment mode: reduced + pooled fragment distributions |

Distribution plots are generated for runs that include histogram data in
aggregated results.

Plot settings in `comparison.yaml`:

```yaml
plot_settings:
  rg:
    show_per_replicate: false    # Overlay individual replicate traces
    figsize: [10, 6]             # Default figure size (bar charts)
    timeseries_figsize: [12, 5]  # Timeseries figure size
```

## Common CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Path to comparison configuration |
| `--eq-time` | `0ns` | Equilibration time to skip |
| `--recompute` | off | Ignore cached results and recompute |
| `--format` | `table` | Output format (`table` or `json`) |
| `-o, --output` | (none) | Save formatted output to file |
| `-q, --quiet` | off | Suppress INFO messages |
| `--debug` | off | Enable DEBUG logging |

Replicates are configured per condition in `comparison.yaml`:

```yaml
conditions:
  - label: "no_polymer"
    config: "configs/no_polymer.yaml"
    replicates: [1, 3, 5]
```

## Troubleshooting

### "Selection matched no atoms"

When a run selection matches zero atoms, that run is skipped with a warning.
Analysis continues for runs and conditions with valid data.

If unexpected:
- Confirm residue numbering and atom naming in your topology
- Check selection syntax directly against your system
- Re-run with `--debug` for detailed diagnostics

### "At least one Rg run must be defined"

Cause: `plugins.rg.runs` is missing or empty.

Fix: add at least one run with `label` and `selection`.

### "Rg run labels must be unique"

Cause: duplicate run labels.

Fix: assign a unique `label` to each run.

### "Equilibration removed all frames"

Cause: `--eq-time` exceeds trajectory duration.

Fix: lower equilibration time or verify simulation completion.

### Very large Rg fluctuations (> 5 Å std)

Cause: often unfolding, large flexibility, or an overly broad selection.

Fix:
- Validate the selection (whole protein vs backbone vs core)
- Check trajectory integrity
- Inspect timeseries for transitions or discontinuities

### "Low statistical reliability" warning

Cause: correlation time is long relative to trajectory length.

This is informational. Consider longer simulations or more replicates.

### Missing replicate data

Message: `Skipping replicate N: trajectory data not found`

Cause: trajectory files are missing or replicate is incomplete.

Fix: verify simulation output paths and replicate status.

### Control condition missing for a run in fragment-mode workflows

Message: control has no data for a run and comparison falls back to all-vs-all
pairwise testing for that run.

Cause: control selection matched no atoms for that run.

Behavior: expected in mixed run sets (for example, polymer-only runs with a
no-polymer control).

## Rg vs Other Metrics

### Rg vs RMSD

| Feature | Rg | RMSD |
|---------|----|------|
| **Measures** | Structural compactness (mass-weighted size) | Deviation from a reference structure |
| **Output** | One value per frame (timeseries) | One value per frame (timeseries) |
| **Reference** | None required | Required (`centroid`, `average`, `frame`, or `external`) |
| **Alignment** | Not required | Required |
| **Configuration** | `label` + `selection` | `label` + `selection` + alignment + reference |
| **Direction labels** | `compaction` / `expansion` / `unchanged` | `stabilizing` / `destabilizing` / `unchanged` |
| **Best for** | Compaction, swelling, folding state shifts | Drift and reference-relative stability |

### Rg vs RMSF

| Feature | Rg | RMSF |
|---------|----|------|
| **Measures** | Global compactness over time | Per-residue positional fluctuation |
| **Primary output** | Timeseries | Residue profile |
| **Best question** | Is the structure compacting or expanding? | Which regions are flexible or rigid? |
