# Distances Plugin Reference

For a step-by-step task guide, see
{doc}`../how_to/analysis_distances_quickstart`.

## Configuration Reference

All fields for `plugins.distances`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `threshold` | `float` | `3.5` | Global default threshold in Angstroms |
| `pairs` | `list[DistancePair]` | *required* | One or more named distance pairs |
| `use_pbc` | `bool` | `true` | Apply periodic boundary conditions using minimum-image distance |
| `align_trajectory` | `bool` | `true` | Align trajectory before distance calculation |
| `alignment_selection` | `str` | `"protein and name CA"` | MDAnalysis selection used for alignment |
| `alignment_mode` | `str` | `"centroid"` | Alignment reference mode: `centroid` or `frame` |
| `alignment_frame` | `int \| None` | `null` | Reference frame index when `alignment_mode: frame` |

Each entry in `pairs`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Human-readable pair label |
| `selection_a` | `str` | *required* | Selection for group A |
| `selection_b` | `str` | *required* | Selection for group B |
| `threshold` | `float` | global `threshold` | Per-pair threshold override |
| `below_label` | `str` | `"Below {threshold}Å"` | Label for distance `<= threshold` state |
| `above_label` | `str` | `"Above {threshold}Å"` | Label for distance `> threshold` state |

### Selection syntax extensions

Distances supports MDAnalysis selections plus helper wrappers:

| Syntax | Description | Typical use |
|--------|-------------|-------------|
| `midpoint(selection)` | Geometric midpoint of selected atoms | Carboxylate oxygens (Asp/Glu) |
| `com(selection)` | Center of mass of selected atoms | Whole residues, ligands, or rings |
| `pdbindex N` | Atom by PDB serial index (1-indexed) | Copying atom IDs from PDB/PyMOL |

Example pair definitions:

```yaml
plugins:
  distances:
    pairs:
      - label: "His156-Asp133"
        selection_a: "protein and resid 156 and name ND1"
        selection_b: "midpoint(protein and resid 133 and name OD1 OD2)"
      - label: "Ligand-COM to Ser77"
        selection_a: "com(resname LIG)"
        selection_b: "protein and resid 77 and name OG"
      - label: "Restraint atom check"
        selection_a: "pdbindex 2740"
        selection_b: "pdbindex 3011"
```

```{important}
Residue indices restart by chain in PolyzyMD systems. For protein residues,
prefer `protein and resid ...` to avoid accidental multi-chain matches.
```

### PBC and alignment behavior

- `use_pbc: true` computes minimum-image distances for wrapped trajectories
- `align_trajectory: true` removes global rotation/translation before analysis
- Orthorhombic boxes are fully supported for PBC correction
- Triclinic boxes trigger a warning and fall back to Euclidean distance

Alignment reference options:

| Mode | Description | Best fit |
|------|-------------|----------|
| `centroid` | Align to most populated conformation | General default |
| `frame` | Align to a specific frame index | Reproducible fixed-reference comparisons |

### Cache invalidation by settings

Distances cache keys include equilibration and geometry settings. Changing
PBC/alignment settings produces new cache filenames automatically, for example:

```text
distances_Ser77-His156_eq10ns_pbc_align-centroid.json
distances_Ser77-His156_eq10ns_nopbc_noalign.json
```

## Output Files

Results are saved under your project analysis directory:

```text
<projects_directory>/
└── analysis/
    └── distances/
        ├── run_1/
        │   ├── distances_Ser77-His156_eq10ns_pbc_align-centroid.json
        │   └── distances_His156-Asp133_eq10ns_pbc_align-centroid.json
        ├── run_2/
        │   └── ...
        ├── run_3/
        │   └── ...
        └── aggregated/
            └── distances_reps1-3_eq10ns.json
```

Per-replicate result structure (representative):

```python
{
    "config_hash": "abc123...",
    "replicate": 1,
    "equilibration_time": 10.0,
    "equilibration_unit": "ns",
    "n_frames_total": 10000,
    "n_frames_used": 9000,
    "pair_results": [
        {
            "pair_label": "Ser77-His156",
            "selection1": "protein and resid 77 and name OG",
            "selection2": "protein and resid 156 and name NE2",
            "mean_distance": 3.42,
            "std_distance": 0.87,
            "sem_distance": 0.15,
            "median_distance": 3.31,
            "min_distance": 2.61,
            "max_distance": 5.87,
            "kde_peak": 3.18,
            "threshold": 3.5,
            "fraction_below_threshold": 0.624,
            "correlation_time": 245.3,
            "n_independent_frames": 34,
            "histogram_edges": [...],
            "histogram_counts": [...],
            "kde_x": [...],
            "kde_y": [...]
        }
    ]
}
```

Aggregated result structure (representative):

```python
{
    "replicates": [1, 2, 3],
    "n_replicates": 3,
    "pair_summaries": [
        {
            "pair_label": "Ser77-His156",
            "overall_mean": 3.39,
            "overall_sem": 0.11,
            "overall_median": 3.28,
            "per_replicate_means": [3.42, 3.31, 3.45],
            "per_replicate_sems": [0.15, 0.13, 0.16],
            "threshold": 3.5,
            "overall_fraction_below_threshold": 0.61
        }
    ]
}
```

## Plot Types

Generate figures with:

```bash
polyzymd compare plot-all -f comparison.yaml
```

Distances plot outputs:

| Plot output | Description |
|-------------|-------------|
| `distance_kde_<pair>.png` | Distribution overlays across conditions for one pair |
| `distance_threshold_bars.png` | Grouped bars of fraction below threshold |
| `distance_state_<pair>.png` | Per-pair below/above threshold state summary |

`plot_settings.distances` options:

| Field | Default | Description |
|-------|---------|-------------|
| `show_threshold` | `true` | Draw threshold line on distributions |
| `use_kde` | `true` | Use KDE overlays (else histogram emphasis) |
| `generate_state_bars` | `true` | Generate per-pair state bar figures |

## Common CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Comparison config path |
| `--eq-time` | `0ns` | Equilibration time to discard |
| `--recompute` | off | Ignore cached results and recompute |
| `--format` | `table` | Output format (`table` or `json`) |
| `-o, --output` | (none) | Write formatted output to a file |
| `-q, --quiet` | off | Suppress INFO logs |
| `--debug` | off | Enable DEBUG logging |

Typical run commands:

```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
polyzymd compare run distances -f comparison.yaml --eq-time 10ns --recompute --format json
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```

## Troubleshooting

### "Selection matched no atoms"

**Cause:** Selection string does not match topology atoms.

**Fix:**
- Add chain-aware qualifiers, such as `protein and resid ...`
- Verify atom names and residue IDs in your topology
- Re-run with `--debug` for expanded selection diagnostics

### Very wide or multimodal distance distribution

**Cause:** Selection may include flexible groups or unintended atoms.

**Fix:**
- Confirm each selection resolves to the intended atom/group
- Use `midpoint(...)` or `com(...)` where chemically appropriate
- Inspect atom selections in a molecular viewer

### Apparent long-distance outliers near box boundaries

**Cause:** PBC handling disabled or unsupported box geometry.

**Fix:**
- Ensure `use_pbc: true`
- Check logs for triclinic fallback warnings
- Compare with aligned trajectories to reduce rigid-body artifacts

### "Low statistical reliability" warning

**Cause:** Correlation time is large relative to trajectory length.

**Fix:**
- Add replicates and compare aggregated results
- Extend production trajectory length
- Treat uncertainty estimates as conservative qualitative guidance

### Missing replicate data

**Message:** `Skipping replicate N: trajectory data not found`

**Cause:** Replicate output is missing or incomplete.

**Fix:**
- Confirm the requested replicate finished simulation
- Verify scratch/project path mapping in config
- Re-run analysis after data is available
