# RMSF Plugin Reference

For a step-by-step guide to running RMSF analysis, see
{doc}`../how_to/analysis_rmsf_quickstart`.

## Configuration Reference

RMSF settings live under `plugins.rmsf` in `comparison.yaml`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `selection` | `str` | `"protein and name CA"` | MDAnalysis selection used for RMSF calculation |
| `reference_mode` | `str` | `"centroid"` | Alignment reference mode: `centroid`, `average`, `frame`, `external` |
| `reference_frame` | `int \| null` | `null` | 1-indexed frame when `reference_mode: frame` |
| `reference_file` | `str \| null` | `null` | Path to external PDB when `reference_mode: external` |
| `alignment_selection` | `str` | `"protein and name CA"` | Selection used for trajectory alignment |
| `centroid_selection` | `str` | `"protein"` | Selection used to find centroid reference frame |

```{note}
Validation rules:

- `reference_mode: frame` requires `reference_frame`
- `reference_mode: external` requires `reference_file`
- `reference_file` must point to an existing PDB file
```

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

```{important}
When `reference_mode: external` or another explicit reference mode is used,
RMSF values are deviations about that specified reference structure. They are
not necessarily fluctuations about the trajectory mean. Use this mode when the
scientific question is reference deviation, such as motion away from a crystal
or catalytically competent conformation.
```

## Output Files

RMSF writes canonical v1.3 artifact outputs under each condition's analysis
directory. Artifact JSON files use stable envelope fields instead of exposing
plugin-private result model classes as public output schemas.

```text
<comparison_workspace>/
├── analysis/
│   └── <condition>/
│       └── rmsf/
│           ├── run_1/
│           │   ├── result.json
│           │   └── sidecars/
│           │       └── rmsf_profile.npz
│           ├── run_2/
│           │   └── ...
│           ├── run_3/
│           │   └── ...
│           └── aggregated/
│               ├── result.json
│               └── sidecars/
│                   └── rmsf_profile.npz
└── comparison/
    └── rmsf/
        └── result.json
```

The canonical paths are:

| Level | Artifact | Path |
|-------|----------|------|
| Per replicate | `ReplicateArtifact` | `analysis/<condition>/rmsf/run_<replicate>/result.json` |
| Per condition | `ConditionArtifact` | `analysis/<condition>/rmsf/aggregated/result.json` |
| Cross condition | Comparison result | `comparison/rmsf/result.json` |
| Large arrays | NPZ sidecars | `analysis/<condition>/rmsf/**/sidecars/*.npz` |

### Artifact envelope fields

All RMSF artifact JSON files use common envelope fields:

| Field | Description |
|-------|-------------|
| `payload` | JSON-compatible RMSF metrics, summaries, and relative sidecar references |
| `metadata` | Plugin settings, equilibration labels, units, and other cache metadata |
| `provenance` | Input identity and workflow details needed to audit how the artifact was produced |
| `sidecars` | Validated references to large arrays or profiles in `sidecars/*.npz`, with hashes and sizes |

Use the public artifact store API to inspect artifacts programmatically:

```python
from pathlib import Path

from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact, ReplicateArtifact

rep_store = ArtifactStore(Path("analysis/PEGylated/rmsf/run_1"))
replicate: ReplicateArtifact = rep_store.read_replicate_result()
print(replicate.payload["metrics"])

condition_store = ArtifactStore(Path("analysis/PEGylated/rmsf/aggregated"))
condition: ConditionArtifact = condition_store.read_condition_result()
print(condition.payload["metrics"])
```

### Per-replicate JSON (`ReplicateArtifact`)

Representative structure:

```python
{
    "schema_version": "1",
    "artifact_type": "replicate",
    "analysis_name": "rmsf",
    "condition_label": "PEGylated",
    "replicate": 1,
    "payload": {
        "metrics": {"mean_rmsf": 0.621},
        "residue_ids": [1, 2, 3],
        "residue_names": ["MET", "ALA", "SER"],
        "profile_sidecar": "sidecars/rmsf_profile.npz",
        "summary": {
            "std_rmsf": 0.215,
            "min_rmsf": 0.248,
            "max_rmsf": 3.160,
            "n_frames_used": 9000
        }
    },
    "metadata": {
        "selection": "protein and name CA",
        "reference_mode": "centroid",
        "alignment_selection": "protein and name CA",
        "equilibration": "10ns"
    },
    "provenance": {"trajectory_files": [".../prod_1.xtc"]},
    "sidecars": [
        {"path": "sidecars/rmsf_profile.npz", "metadata": {"kind": "rmsf_profile"}}
    ]
}
```

The sidecar stores arrays such as `rmsf_values`, residue indices, and any
profile data that should not be duplicated into JSON.

### Aggregated JSON (`ConditionArtifact`)

Representative structure:

```python
{
    "schema_version": "1",
    "artifact_type": "condition",
    "analysis_name": "rmsf",
    "condition_label": "PEGylated",
    "replicates": [1, 2, 3],
    "payload": {
        "metrics": {
            "mean_rmsf": {"values": [0.64, 0.59, 0.63], "mean": 0.62, "sem": 0.02}
        },
        "profile_sidecar": "sidecars/rmsf_profile.npz",
        "summary": {"overall_min_rmsf": 0.30, "overall_max_rmsf": 4.21}
    },
    "metadata": {"equilibration": "10ns"},
    "provenance": {"source_replicates": [1, 2, 3]},
    "sidecars": [
        {"path": "sidecars/rmsf_profile.npz", "metadata": {"kind": "mean_sem_profile"}}
    ]
}
```

### Comparison JSON (`result.json`)

```python
{
    "metric": "rmsf",
    "conditions": [
        {
            "label": "No Polymer",
            "n_replicates": 3,
            "mean_rmsf": 0.715,
            "sem_rmsf": 0.020,
            "replicate_values": [0.755, 0.693, 0.696]
        },
        {
            "label": "With Polymer",
            "n_replicates": 3,
            "mean_rmsf": 0.551,
            "sem_rmsf": 0.034,
            "replicate_values": [0.590, 0.520, 0.542]
        }
    ],
    "pairwise_comparisons": [
        {
            "condition_a": "No Polymer",
            "condition_b": "With Polymer",
            "percent_change": -22.9,
            "p_value": 0.0211,
            "cohens_d": 4.06,
            "significant": true,
            "direction": "stabilizing"
        }
    ],
    "ranking": ["With Polymer", "No Polymer"]
}
```

## Plot Types

RMSF plots are generated by `polyzymd compare plot-all -f comparison.yaml`.

| Plot output | Description |
|-------------|-------------|
| `rmsf_profile.png` | Per-residue RMSF profile by condition; optional SEM shading |
| `rmsf_comparison.png` | Horizontal bar chart of condition-level mean RMSF with SEM |

The profile plot can include a reference secondary-structure annotation row
when `reference_mode: external` uses a readable `reference_file`. The
annotation is computed during RMSF aggregation and stored in the condition
artifact, so plotting remains artifact-only and does not reload the reference
PDB. Recompute RMSF if older artifacts were created before this annotation was
available.

### RMSF plot settings

```yaml
plot_settings:
  rmsf:
    show_error: true                # Show SEM band/bars
    show_reference_secondary_structure: true  # Show external-reference structure strip
    highlight_residues: [77, 133]   # Vertical guide lines in profile plot
    figsize_profile: [14, 4]        # Profile figure size
    figsize_comparison: [8, 6]      # Comparison figure size
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `show_error` | `bool` | `true` | Show SEM shading/bars |
| `show_reference_secondary_structure` | `bool` | `true` | Show the cached external-reference secondary-structure strip when present |
| `highlight_residues` | `list[int]` | `[]` | Residue numbers to mark on profile plot |
| `figsize_profile` | `tuple[float, float]` | `[14, 4]` | Profile figure size |
| `figsize_comparison` | `tuple[float, float]` | `[8, 6]` | Comparison figure size |

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

**Cause:** The MDAnalysis selection does not match atoms in the topology.

**Fix:**

- Check residue numbering and atom names in your input structure
- Start with `selection: "protein and name CA"`
- Re-run with `--debug` for detailed selection diagnostics

### "reference_file does not exist"

**Cause:** `reference_mode: external` is set, but the path is invalid.

**Fix:** Use an absolute path or a path relative to your working directory.

### "External PDB atom count does not match trajectory selection"

**Cause:** The `selection` string resolves to different atom counts in
trajectory and external reference.

**Fix:**

- Ensure both structures use compatible atom naming
- Use a stricter selection that matches in both systems
- Confirm the external PDB contains the same residue set

### Very high RMSF values (> 10 Å)

**Cause:** Usually alignment mismatch, overly broad selection, or genuine
structural instability.

**Fix:**

- Verify `alignment_selection` and `selection`
- Try `reference_mode: "average"` as a cross-check
- Confirm trajectory files are complete

### "Low statistical reliability" warning

**Cause:** Correlation time is large relative to available production data.

**Fix:**

- Use more replicates
- Extend simulation length
- Treat the result as qualitative if uncertainty is large

For interpretation guidance, see
{doc}`../explanation/analysis_rmsf_best_practices`.

### Missing replicate data

**Message:** `Skipping replicate N: trajectory data not found`

**Cause:** Replicate output is missing or path configuration is incorrect.

**Fix:** Analysis continues with available replicates. Verify simulation
completion and file paths if missing replicates are unexpected.
