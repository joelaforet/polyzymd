# Catalytic Triad Plugin Reference

For a step-by-step workflow, see
{doc}`../how_to/analysis_triad_quickstart`.

## Configuration Reference

Top-level plugin key in `comparison.yaml`:

```yaml
plugins:
  catalytic_triad:
    name: "LipA_catalytic_triad"
    description: "Ser-His-Asp catalytic triad"
    threshold: 3.5
    pairs:
      - label: "Asp133-His156"
        selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
        selection_b: "protein and resid 156 and name ND1"
      - label: "His156-Ser77"
        selection_a: "protein and resid 156 and name NE2"
        selection_b: "protein and resid 77 and name OG"
```

### `CatalyticTriadSettings` fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `name` | `str` | `"catalytic_triad"` | Triad identifier saved in output files |
| `description` | `str \| null` | `null` | Optional human-readable triad description |
| `threshold` | `float` | `3.5` | Contact cutoff in angstroms for each pair |
| `pairs` | `list[TriadPair]` | *required* | Pair definitions used to compute distances and contacts |

### `TriadPair` fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Pair label used in tables and JSON |
| `selection_a` | `str` | *required* | First MDAnalysis selection |
| `selection_b` | `str` | *required* | Second MDAnalysis selection |

### Selection syntax

| Syntax | Description | Example |
|--------|-------------|---------|
| Standard | MDAnalysis atom selection | `protein and resid 77 and name OG` |
| `midpoint()` | Geometric midpoint of selected atoms | `midpoint(protein and resid 133 and name OD1 OD2)` |
| `com()` | Center of mass of selected atoms | `com(protein and resid 133 and name OD1 OD2)` |

```{important}
Use chain-aware selections. Residue IDs restart per chain in PolyzyMD systems,
so bare `resid X` selections can match non-protein atoms. Prefer
`protein and resid X` for catalytic residues.
```

## Output Files

Triad outputs are written as canonical v1.3 artifacts under each condition plus
a comparison result. JSON files use framework-owned artifact envelopes instead
of exposing plugin-private result model classes as public schemas.

```text
comparison_workspace/
├── analysis/
│   └── <condition>/
│       └── catalytic_triad/
│           ├── run_1/
│           │   ├── result.json
│           │   └── sidecars/
│           │       └── triad_distances.npz
│           ├── run_2/
│           │   └── ...
│           ├── run_3/
│           │   └── ...
│           └── aggregated/
│               ├── result.json
│               └── sidecars/
│                   └── triad_distance_profiles.npz
└── comparison/
    └── catalytic_triad/
        └── result.json
```

The canonical paths are:

| Level | Artifact | Path |
|-------|----------|------|
| Per replicate | `ReplicateArtifact` | `analysis/<condition>/catalytic_triad/run_<replicate>/result.json` |
| Per condition | `ConditionArtifact` | `analysis/<condition>/catalytic_triad/aggregated/result.json` |
| Cross condition | Comparison result | `comparison/catalytic_triad/result.json` |
| Large arrays | NPZ sidecars | `analysis/<condition>/catalytic_triad/**/sidecars/*.npz` |

### Artifact envelope fields

| Field | Description |
|-------|-------------|
| `payload` | Triad contact metrics, pair summaries, and relative sidecar paths |
| `metadata` | Triad settings, cutoff threshold, equilibration labels, and units |
| `provenance` | Input topology/trajectory identity and workflow details |
| `sidecars` | Validated references to `sidecars/*.npz` arrays with hashes and sizes |

Use the public artifact API to load saved triad artifacts:

```python
from pathlib import Path

from polyzymd.analyses.mda import ArtifactStore

replicate = ArtifactStore(Path("analysis/PEGylated/catalytic_triad/run_1"))
condition = ArtifactStore(Path("analysis/PEGylated/catalytic_triad/aggregated"))
print(replicate.read_replicate_result().payload["metrics"])
print(condition.read_condition_result().payload["metrics"])
```

### Per-replicate JSON structure (`ReplicateArtifact`)

Representative structure:

```python
{
    "schema_version": "1",
    "artifact_type": "replicate",
    "analysis_name": "catalytic_triad",
    "condition_label": "PEGylated",
    "replicate": 1,
    "payload": {
        "metrics": {"simultaneous_contact_fraction": 0.741},
        "pair_results": [
            {"pair_label": "Asp133-His156", "mean_distance": 2.91, "fraction_below_threshold": 0.932},
            {"pair_label": "His156-Ser77", "mean_distance": 4.07, "fraction_below_threshold": 0.551}
        ],
        "distance_sidecar": "sidecars/triad_distances.npz",
        "summary": {"n_frames_used": 2000, "n_frames_total": 3000}
    },
    "metadata": {"triad_name": "LipA_catalytic_triad", "threshold": 3.5, "equilibration": "100ns"},
    "provenance": {"trajectory_files": ["..."]},
    "sidecars": [
        {"path": "sidecars/triad_distances.npz", "metadata": {"kind": "distance_timeseries"}}
    ]
}
```

### Aggregated JSON structure (`ConditionArtifact`)

Representative structure:

```python
{
    "schema_version": "1",
    "artifact_type": "condition",
    "analysis_name": "catalytic_triad",
    "condition_label": "PEGylated",
    "replicates": [1, 2, 3],
    "payload": {
        "metrics": {
            "simultaneous_contact_fraction": {
                "values": [0.741, 0.512, 0.245],
                "mean": 0.499,
                "sem": 0.273
            }
        },
        "pair_results": [
            {"pair_label": "Asp133-His156", "overall_mean": 3.09, "overall_sem": 0.21},
            {"pair_label": "His156-Ser77", "overall_mean": 4.03, "overall_sem": 1.07}
        ],
        "distance_profile_sidecar": "sidecars/triad_distance_profiles.npz"
    },
    "metadata": {"triad_name": "LipA_catalytic_triad", "threshold": 3.5, "equilibration": "100ns"},
    "provenance": {"source_replicates": [1, 2, 3]},
    "sidecars": [
        {"path": "sidecars/triad_distance_profiles.npz", "metadata": {"kind": "distance_profiles"}}
    ]
}
```

## Plot Types

Catalytic triad plots are generated by `polyzymd compare plot-all`.

| Plot output | Description |
|-------------|-------------|
| `triad_kde_panel.png` | Multi-row KDE of per-pair distance distributions across conditions |
| `triad_threshold_bars.png` | Grouped bar chart of the fraction below threshold for each pair |

Optional historical plotting helpers in the code can generate 2D KDE figures,
but the current plugin `plot()` lifecycle emits the KDE panel and threshold bar
figures listed above.

Use top-level `plot_settings` in `comparison.yaml` to tune sizes, styles,
and output resolution for triad plots.

## Common CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Path to comparison config |
| `--eq-time` | from YAML defaults | Equilibration time to skip |
| `--recompute` | off | Ignore cached results and recompute |
| `--format` | `table` | Output format: `table` or `json` |
| `-q, --quiet` | off | Suppress INFO logging |
| `--debug` | off | Enable DEBUG logging |

Use either plugin name:

```bash
polyzymd compare run triad -f comparison.yaml --eq-time 100ns
polyzymd compare run catalytic_triad -f comparison.yaml --eq-time 100ns
```

## Troubleshooting

### "No catalytic_triad section found"

**Cause:** `plugins.catalytic_triad` is missing from `comparison.yaml`

**Fix:** Add a complete `catalytic_triad` block with `name`, `threshold`, and
at least one pair in `pairs`

### "Selection ... returned 0 atoms"

**Cause:** Residue ID or atom name does not match the loaded topology

**Fix:**
- Verify residue numbering in the topology used for analysis
- Verify atom names (`OD1`/`OD2`, `ND1`/`NE2`, `OG`)
- Use chain-aware selections like `protein and resid ...`

### Very high distances (> 10 Å)

**Cause:** Usually incorrect selections or mismatched residue numbering

**Fix:** Validate selections in VMD or PyMOL and confirm the exact residues

### Very low simultaneous contact (< 10%)

**Cause:** Could be a genuine disrupted triad or a too-strict threshold

**Fix:**
- Inspect pair-level fractions to identify the limiting pair
- Re-check selections
- If scientifically justified, test a slightly larger threshold (for example 4.0 Å)

### "Low statistical reliability" warning

**Cause:** Long autocorrelation time relative to trajectory length

**Fix:** Informational warning only. Prefer more replicates and/or longer
simulations for tighter uncertainty

### "Skipping replicate N: trajectory data not found"

**Cause:** Replicate output is missing or path mapping is incorrect

**Fix:** Analysis continues with available data. Verify simulation completion and
project paths if this is unexpected
