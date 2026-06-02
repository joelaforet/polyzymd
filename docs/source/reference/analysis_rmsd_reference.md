# RMSD Plugin Reference

For a step-by-step guide to running RMSD analysis, see
{doc}`../how_to/analysis_rmsd_quickstart`.

## Configuration Reference

All fields for `RMSDRunSettings`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Human-readable run label (must be unique) |
| `selection` | `str` | `"protein and name CA"` | MDAnalysis selection for RMSD calculation |
| `alignment_selection` | `str` | `"protein and name CA"` | MDAnalysis selection for trajectory alignment |
| `reference_mode` | `str` | `"centroid"` | Reference mode: `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | `int` | `0` | 0-indexed frame for `reference_mode: frame` |
| `reference_file` | `str \| None` | `null` | Path to external PDB for `reference_mode: external` |
| `centroid_selection` | `str \| None` | `null` | Selection for centroid finding; defaults to `alignment_selection` |
| `convergence_window_size_ns` | `float` | `15.0` | Sliding window size for convergence detection (ns) |
| `convergence_step_size_ns` | `float` | `5.0` | Step between successive windows (ns) |
| `convergence_slope_threshold` | `float` | `0.0005` | Max absolute slope to qualify as "flat" (Å/ns) |
| `convergence_sustained_for_ns` | `float` | `15.0` | Required sustained duration below threshold (ns) |

Top-level `RMSDSettings` contains a single field:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | `list[RMSDRunSettings]` | *required* | One or more named RMSD runs (at least one required) |

```{note}
Run labels must be unique within a single `comparison.yaml`. Duplicate labels
raise a validation error.
```

## Output Files

Results are saved as canonical v1.3 artifacts. JSON files are stable artifact
envelopes, while per-frame RMSD arrays are stored in NPZ sidecars.

```text
<comparison_workspace>/
├── analysis/
│   └── <condition>/
│       └── rmsd/
│           ├── run_1/
│           │   ├── result.json
│           │   └── sidecars/
│           │       ├── rmsd_Protein_Backbone_timeseries.npz
│           │       └── rmsd_Active_Site_timeseries.npz
│           ├── run_2/
│           │   └── ...
│           ├── run_3/
│           │   └── ...
│           └── aggregated/
│               ├── result.json
│               └── sidecars/
│                   └── rmsd_Protein_Backbone_timeseries.npz
└── comparison/
    └── rmsd/
        └── result.json
```

The canonical paths are:

| Level | Artifact | Path |
|-------|----------|------|
| Per replicate | `ReplicateArtifact` | `analysis/<condition>/rmsd/run_<replicate>/result.json` |
| Per condition | `ConditionArtifact` | `analysis/<condition>/rmsd/aggregated/result.json` |
| Cross condition | Comparison result | `comparison/rmsd/result.json` |
| Large arrays | NPZ sidecars | `analysis/<condition>/rmsd/**/sidecars/*.npz` |

Each replicate artifact contains JSON-compatible summaries for all configured
runs. Raw per-frame RMSD timeseries are sidecars referenced from `payload` and
listed in `sidecars` with recorded size and hash metadata.

### Artifact envelope fields

| Field | Description |
|-------|-------------|
| `payload` | RMSD run summaries, scalar metrics, convergence diagnostics, and sidecar paths |
| `metadata` | Settings such as selections, reference modes, equilibration labels, and units |
| `provenance` | Input topology/trajectory identity and workflow details |
| `sidecars` | Validated references to `sidecars/*.npz` arrays for timeseries and aggregate profiles |

Use `ArtifactStore` for programmatic access:

```python
from pathlib import Path

from polyzymd.analyses.mda import ArtifactStore

replicate = ArtifactStore(Path("analysis/PEGylated/rmsd/run_1")).read_replicate_result()
condition = ArtifactStore(Path("analysis/PEGylated/rmsd/aggregated")).read_condition_result()
print(replicate.payload["runs"][0]["mean_rmsd"])
print(condition.payload["runs"][0]["metrics"]["mean_rmsd"])
```

### JSON result structure

Per-replicate result (`ReplicateArtifact`), representative structure:

```python
{
    "schema_version": "1",
    "artifact_type": "replicate",
    "analysis_name": "rmsd",
    "condition_label": "PEGylated",
    "replicate": 1,
    "payload": {
        "runs": [
            {
                "run_label": "Protein Backbone",
                "selection": "protein and name CA",
                "alignment_selection": "protein and name CA",
                "reference_mode": "centroid",
                "mean_rmsd": 1.823,
                "std_rmsd": 0.312,
                "median_rmsd": 1.791,
                "sem_rmsd": 0.078,
                "converged": true,
                "convergence_time_ns": 12.5,
                "timeseries_sidecar": "sidecars/rmsd_Protein_Backbone_timeseries.npz"
            }
        ]
    },
    "metadata": {"equilibration": "10ns", "time_unit": "ns"},
    "provenance": {"trajectory_files": ["..."], "n_frames_used": 9000},
    "sidecars": [
        {
            "path": "sidecars/rmsd_Protein_Backbone_timeseries.npz",
            "metadata": {"kind": "timeseries", "run_label": "Protein Backbone"}
        }
    ]
}
```

Aggregated result (`ConditionArtifact`), representative structure:

```python
{
    "schema_version": "1",
    "artifact_type": "condition",
    "analysis_name": "rmsd",
    "condition_label": "PEGylated",
    "replicates": [1, 2, 3],
    "payload": {
        "runs": [
            {
                "run_label": "Protein Backbone",
                "selection": "protein and name CA",
                "metrics": {
                    "mean_rmsd": {"values": [1.823, 1.891, 1.854], "mean": 1.856, "sem": 0.034}
                },
                "convergence": {
                    "n_converged_replicates": 3,
                    "convergence_fraction": 1.0,
                    "mean_convergence_time_ns": 13.2
                }
            }
        ]
    },
    "metadata": {"equilibration": "10ns"},
    "provenance": {"source_replicates": [1, 2, 3]}
}
```

## Plot Types

The RMSD plugin generates figures through `polyzymd compare plot-all`:

| Plot output | Description |
|-------------|-------------|
| `rmsd_timeseries_<run>.png` | Mean RMSD vs time with SEM shading, one per run |
| `rmsd_comparison_<run>.png` | Grouped bar chart of mean RMSD across conditions, one per run |
| `rmsd_convergence_<condition>_<run>.png` | Dual-axis plot: RMSD timeseries with sliding-window slope and convergence marker (requires `show_convergence_plots: true`) |

**Timeseries plot features:**
- Mean RMSD curve per condition with SEM shading
- Legend placed outside the plot area (`bbox_to_anchor=(1.02, 0.5)`)
- Optional per-replicate traces via `show_per_replicate: true`

RMSD plot behavior can be customized in `comparison.yaml`:

```yaml
plot_settings:
  rmsd:
    show_per_replicate: false    # Overlay individual replicate traces
    figsize: [10, 6]             # Default figure size (bar charts)
    timeseries_figsize: [12, 5]  # Timeseries figure size (wider)
    show_convergence_plots: false  # Generate per-replicate convergence diagnostics
    convergence_figsize: [12, 5]   # Convergence panel figure size
```

## Convergence Detection

Convergence detection is always on — every RMSD run automatically applies a
sliding-window slope heuristic to determine whether the RMSD timeseries has
plateaued. This is a purely additive diagnostic: it does not affect ranking,
statistical tests, or any other comparison output. Convergence results appear
as additional fields in per-replicate and aggregated JSON files, and optional
convergence plots can be enabled via `show_convergence_plots: true`.

For a conceptual explanation of the algorithm, its parameters, and its
limitations, see {doc}`../explanation/convergence_detection`.

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

## Troubleshooting

### "Selection matched no atoms"

**Cause:** MDAnalysis selection doesn't match any atoms in your topology.

**Fix:**
- Check residue numbering in your PDB vs. MDAnalysis (0-indexed vs 1-indexed)
- Verify atom names match your topology
- Use `polyzymd --debug compare run rmsd -f comparison.yaml ...` for detailed
  diagnostics

### "At least one RMSD run must be defined"

**Cause:** The `runs` list in `plugins.rmsd` is empty or missing.

**Fix:** Add at least one run entry with a `label` field:

```yaml
plugins:
  rmsd:
    runs:
      - label: "Protein Backbone"
```

### "reference_file does not exist"

**Cause:** Using `reference_mode: external` but the PDB path is invalid.

**Fix:** Provide an absolute path or a path relative to the working directory:

```yaml
reference_mode: "external"
reference_file: "/absolute/path/to/crystal.pdb"
```

### "atom count mismatch between trajectory and external PDB"

**Cause:** The `selection` string matches different numbers of atoms in the
trajectory vs. the external reference PDB.

**Fix:**
- Ensure both systems use the same atom naming convention
- Check that the external PDB contains the same residues as your simulation
- Use a more specific selection if topologies differ

### Very high RMSD values (> 10 Å)

**Cause:** Usually indicates alignment issues, wrong selection, or unfolding.

**Fix:**
- Check that `alignment_selection` matches atoms in your system
- Try `reference_mode: "average"` to compare
- Verify trajectory files are complete
- Check for protein unfolding or large conformational changes

### "Low statistical reliability" warning

**Cause:** Long correlation time relative to trajectory length.

**This is informational, not an error.** Results are still valid but
uncertainties may be underestimated.

**Mitigation:**
- Use multiple replicates (aggregated SEM is more reliable)
- Run longer simulations
- Results are still useful for qualitative comparisons

### Missing replicate data

**Message:** `Skipping replicate N: trajectory data not found`

**Cause:** The requested replicate hasn't completed or path is incorrect.

**Fix:** This is informational — analysis continues with available replicates.
Check simulation status if unexpected.

## RMSD vs RMSF Comparison

| Feature | RMSD | RMSF |
|---------|------|------|
| **Measures** | Global deviation from reference | Per-residue fluctuation |
| **Output** | One value per frame (timeseries) | One value per residue (profile) |
| **Reference** | Fixed structure (centroid/average/external) | Time-averaged position |
| **Detects** | Conformational drift, unfolding | Flexible loops, rigid core |
| **Multi-run** | Yes (`runs` list with different selections) | Single selection |
| **Best for** | Equilibration assessment, stability comparison | Flexibility mapping |

```{tip}
Use RMSD first to assess overall stability and choose equilibration time,
then use RMSF to identify which regions drive flexibility differences.
```
