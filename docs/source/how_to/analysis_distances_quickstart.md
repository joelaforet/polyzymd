# Distance Analysis: Quick Start

Compute inter-atomic distances with proper statistical handling in under 5 minutes.

```{note}
**Want to understand the statistics?** This guide focuses on getting results quickly.
For proper uncertainty quantification (autocorrelation correction, SEM vs. SD),
see the [Statistical Best Practices Guide](../explanation/analysis_statistics_best_practices.md).
```

## TL;DR

```bash
# Configure distance pairs in comparison.yaml, then run:
polyzymd compare run distances -f comparison.yaml --eq-time 10ns

# Run all enabled analyses in the same workflow
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute and machine-readable output
polyzymd compare run distances -f comparison.yaml --eq-time 10ns --recompute --format json
```

## Prerequisites

Before running distance analysis, you need:

1. **Completed production simulation(s)** - at least one replicate
2. **Comparison config** - `comparison.yaml` with conditions and plugin settings
3. **Trajectory files** - in the scratch directory specified in config

Verify your setup:

```bash
# Check that trajectories exist
ls $(polyzymd info -c config.yaml --scratch-dir)/production_*/
```

## What Distance Analysis Provides

The distance analysis module computes:

| Feature | Description |
|---------|-------------|
| **Mean distance** | Average distance over trajectory (equilibrated portion) |
| **SEM** | Autocorrelation-corrected standard error of the mean |
| **Mode (KDE peak)** | Most probable distance from kernel density estimation |
| **Contact fraction** | \% of frames below a distance threshold |
| **Distribution** | Full histogram and KDE for visualization |

```{tip}
**When to use distances vs. contacts vs. triad:**
- **Distances**: Specific atom pairs with continuous distance values
- **Contacts**: All residue-residue contacts at an interface (binary count)
- **Triad**: Pre-defined catalytic geometry with simultaneous contact analysis
```

## Basic Usage

`````{tab-set}

````{tab-item} YAML (Recommended)
For reproducible analysis, define distance pairs in `comparison.yaml`:

```yaml
# comparison.yaml
name: "distance_quickstart"
control: "no_polymer"

conditions:
  - label: "no_polymer"
    config: "configs/no_polymer.yaml"
    replicates: [1, 2, 3]
  - label: "with_polymer"
    config: "configs/with_polymer.yaml"
    replicates: [1, 2, 3]

plugins:
  distances:
    enabled: true
    pairs:
      - label: "Ser77-His156"
        selection_a: "protein and resid 77 and name OG"
        selection_b: "protein and resid 156 and name NE2"
      - label: "His156-Asp133"
        selection_a: "protein and resid 156 and name ND1"
        selection_b: "midpoint(protein and resid 133 and name OD1 OD2)"
```

Then run:

```bash
# Run distance analysis only
polyzymd compare run distances -f comparison.yaml --eq-time 10ns

# Run all enabled plugins in comparison.yaml
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute
polyzymd compare run distances -f comparison.yaml --eq-time 10ns --recompute
```

**Benefits:**
- Version-controlled, reproducible
- Self-documenting experiment setup
- Easy to re-run with different parameters
````

````{tab-item} CLI
### Single analysis run

```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
```

**Expected behavior:**

```
Loading comparison config from: comparison.yaml
Running plugin: distances
  Equilibration: 10ns
  Conditions: no_polymer, with_polymer
Distance comparison complete
```

### Multiple pairs

Define multiple pairs in `comparison.yaml`, then run the same command:

```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
```

### All enabled analyses

Run distances plus any other enabled plugins:

```bash
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```
````

````{tab-item} Python
Use the analysis plugin workflow for reproducible results.

```python
from pathlib import Path
import subprocess

# Distance analysis is configured in comparison.yaml and executed
# through the compare workflow
comparison_file = Path("comparison.yaml")
subprocess.run(
    [
        "polyzymd",
        "compare",
        "run",
        "distances",
        "-f",
        str(comparison_file),
        "--eq-time",
        "10ns",
    ],
    check=True,
)
```

`DistanceCalculator` is still part of `polyzymd.analyses.distances`, but it is an
internal implementation detail used by the distances and catalytic triad plugins.
For user workflows, prefer `comparison.yaml` + `polyzymd compare run`.
````

`````

## Contact Threshold Analysis

Set `threshold` in `comparison.yaml` to compute the fraction of frames where the
distance is below a cutoff (useful for hydrogen bond analysis, active site proximity, etc.):

`````{tab-set}

````{tab-item} YAML (Recommended)
```yaml
# comparison.yaml
plugins:
  distances:
    enabled: true
    threshold: 3.5  # Angstroms (H-bond cutoff)
    pairs:
      - label: "Ser77-His156"
        selection_a: "protein and resid 77 and name OG"
        selection_b: "protein and resid 156 and name NE2"
```

```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
```
````

````{tab-item} CLI
```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
```

**Output:**

```
Distance Analysis Complete
  resid77_OG-resid133_NE2:
    Mean: 3.42 ± 0.15 Å
    Min:  2.61 Å
    Max:  5.87 Å
    Contact fraction (<3.5Å): 62.4%
```
````

````{tab-item} Python
```python
import subprocess

# Threshold is read from plugins.distances.threshold in comparison.yaml
subprocess.run(
    [
        "polyzymd",
        "compare",
        "run",
        "distances",
        "-f",
        "comparison.yaml",
        "--eq-time",
        "10ns",
        "--format",
        "table",
    ],
    check=True,
)
```
````

`````

## Special Selection Syntax

PolyzyMD extends MDAnalysis selections with special position modes and keywords:

:::{warning}
**Chain-Aware Selections Required**

Residue numbers restart at 1 for each chain in PolyzyMD systems. A selection like
`resid 141-148` will match residues from **all chains** (protein, polymer, and water).

For protein residues, always use `protein and resid X`:

```yaml
# INCORRECT - selects from all chains, causing wrong distances
selection_a: "com(resid 141-148)"

# CORRECT - restricts to protein chain only
selection_a: "com(protein and resid 141-148)"
```

PolyzyMD will emit a runtime warning if your selection spans multiple chains,
but it's best to write correct selections from the start.
:::

### Position Modes

| Syntax | Description | Use Case |
|--------|-------------|----------|
| `midpoint(selection)` | Geometric midpoint of selected atoms | Carboxylate groups (Asp, Glu) |
| `com(selection)` | Center of mass of selected atoms | Entire residues, aromatic rings |

### PolyzyMD Keywords

| Keyword | Description | Example |
|---------|-------------|---------|
| `pdbindex N` | Atom by PDB serial number (1-indexed) | `pdbindex 2740 and name CA` |

The `pdbindex` keyword lets you reference atoms by their PDB ATOM serial number
(the number displayed in PyMOL as "id"). This is especially useful when copying
atom selections from restraint definitions in `config.yaml`.

```{tip}
**Consistency with restraints:** You can use the same `pdbindex` selections in
both your restraint configuration (config.yaml) and analysis commands. This
makes it easy to verify that restrained distances match observed distances.
```

### Examples

```yaml
# Midpoint of Asp carboxylate oxygens (protein residue)
selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"

# Center of mass of entire ligand (non-protein, no chain restriction needed)
selection_b: "com(resname LIG)"

# Standard single atom (protein residue)
selection_a: "protein and resid 77 and name OG"

# Atom by PDB serial number (unique, no chain restriction needed)
selection_a: "pdbindex 2740"
```

```{tip}
Use `midpoint()` for carboxylate groups (Asp, Glu) where either oxygen can
accept a hydrogen bond. This gives a single representative point instead of
choosing arbitrarily between OD1/OD2 or OE1/OE2.
```

### CLI execution

Configure quoted selection strings in YAML, then run the plugin:

```yaml
plugins:
  distances:
    enabled: true
    pairs:
      - label: "His156-Asp133"
        selection_a: "protein and resid 156 and name ND1"
        selection_b: "midpoint(protein and resid 133 and name OD1 OD2)"
```

```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
```

## Output Files

Results are saved in your project's analysis directory:

```
<projects_directory>/
└── analysis/
    └── distances/
        ├── run_1/
        │   └── distances_resid77_OG-resid133_NE2_eq10ns.json
        ├── run_2/
        │   └── distances_resid77_OG-resid133_NE2_eq10ns.json
        ├── run_3/
        │   └── distances_resid77_OG-resid133_NE2_eq10ns.json
        └── aggregated/
            └── distances_reps1-3_eq10ns.json
```

### JSON Result Structure

```python
{
    "pair_results": [
        {
            "pair_label": "resid77_OG-resid133_NE2",
            "selection1": "resid 77 and name OG",
            "selection2": "resid 133 and name NE2",
            "mean_distance": 3.42,
            "std_distance": 0.87,
            "sem_distance": 0.15,  # Autocorrelation-corrected
            "median_distance": 3.31,
            "min_distance": 2.61,
            "max_distance": 5.87,
            "kde_peak": 3.18,  # Mode from KDE
            "threshold": 3.5,
            "fraction_below_threshold": 0.624,
            "correlation_time": 245.3,  # ps
            "n_independent_frames": 34,
            "histogram_edges": [...],
            "histogram_counts": [...],
            "kde_x": [...],
            "kde_y": [...]
        }
    ],
    "n_frames_total": 10000,
    "n_frames_used": 9000,
    "equilibration_time": 10.0,
    "equilibration_unit": "ns",
    # ... additional metadata
}
```

## Visualization

`````{tab-set}

````{tab-item} CLI
Generate plugin plots after running comparisons:

```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
polyzymd compare plot-all -f comparison.yaml
```

Plots are saved to `<projects_directory>/plots/distances/`.
````

````{tab-item} Python
Plotting is handled by the distances plugin `plot()` method, which is invoked
through the compare plotting workflow:

```python
import subprocess

subprocess.run(
    ["polyzymd", "compare", "plot-all", "-f", "comparison.yaml"],
    check=True,
)
```
````

`````

### Available plot types

The distances plugin generates figures through `polyzymd compare plot-all`:

| Plot output | Description |
|-------------|-------------|
| `distance_kde_*.png` | KDE distribution overlay per configured distance pair |
| `distance_threshold_bars.png` | Grouped bars for fraction below threshold across conditions |
| `distance_state_*.png` | Per-pair above/below-threshold state summaries |

## Common Options

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Path to comparison configuration |
| `--eq-time` | `0ns` | Equilibration time to skip |
| `--recompute` | off | Ignore cached results and recompute |
| `--format` | `table` | Output format (`table` or `json`) |

### Replicate Specification

Replicates are configured per condition in `comparison.yaml`:

```yaml
conditions:
  - label: "no_polymer"
    config: "configs/no_polymer.yaml"
    replicates: [1, 3, 5]
```

## PBC-Aware Distances and Trajectory Alignment

```{versionadded} 0.3.0
Distance calculations now include PBC-aware distances and trajectory alignment
by default.
```

### Periodic Boundary Conditions (PBC)

By default, distances are computed using the **minimum image convention**, which
correctly handles molecules near periodic boundaries. This prevents artificially
large distances (60-70Å) when atoms are actually close but on opposite sides of
the simulation box.

```{note}
**When does PBC matter?** PBC correction is critical when:
- Molecules diffuse across box boundaries
- Long polymers span the periodic box
- Active sites are near the box edge

For well-centered proteins in large boxes, PBC usually has minimal effect, but
it's always safer to keep it enabled (the default).
```

**Supported box types:**
- ✅ Orthorhombic boxes (cubic, rectangular): Fully supported
- ⚠️ Triclinic boxes: Warning issued, falls back to Euclidean distance

`````{tab-set}

````{tab-item} YAML
```yaml
# comparison.yaml
plugins:
  distances:
    use_pbc: true  # Default, can be omitted
    pairs:
      - label: "Ser77-His156"
        selection_a: "protein and resid 77 and name OG"
        selection_b: "protein and resid 156 and name NE2"
```
````

`````

### Trajectory Alignment

By default, trajectories are **aligned to a reference structure** before computing
distances. This removes rotational drift and center-of-mass motion that can add
noise to distance measurements.

**Why alignment matters:** MD simulations allow the entire system to rotate and
translate. Without alignment, even a rigid protein will show larger fluctuations
in inter-atomic distances due to this global motion.

**Reference modes:**

| Mode | Description | Best for |
|------|-------------|----------|
| `centroid` (default) | Align to most populated conformational cluster (K-Means) | General use |
| `average` | Align to mathematical average structure | Pure thermal fluctuation analysis |
| `frame` | Align to a specific frame number | Comparing to known functional conformation |

```{note}
When alignment is performed, an INFO-level log message notifies you. This
ensures you're aware that trajectory coordinates have been modified in-memory.
```

`````{tab-set}

````{tab-item} YAML
```yaml
# comparison.yaml
plugins:
  distances:
    align_trajectory: true  # Default
    alignment_mode: centroid  # Default
    alignment_selection: "protein and name CA"  # Default
    pairs:
      - label: "Ser77-His156"
        selection_a: "protein and resid 77 and name OG"
        selection_b: "protein and resid 156 and name NE2"
```
````

`````

### Cache Invalidation

The result filename includes PBC and alignment settings, so changing these
parameters automatically invalidates the cache:

```
distances_resid77_OG-resid133_NE2_eq10ns_pbc_align-centroid.json
distances_resid77_OG-resid133_NE2_eq10ns_nopbc_noalign.json
```

This means you can safely experiment with different settings without manually
clearing cached results.

## Troubleshooting

### "Selection matched no atoms"

**Cause:** MDAnalysis selection doesn't match any atoms in your topology.

**Fix:**
- Check residue numbering in your PDB vs. MDAnalysis (0-indexed vs 1-indexed)
- Verify atom names match your topology: `protein and resid 77` to see available atoms
- Use `polyzymd --debug compare run distances -f comparison.yaml ...` for detailed diagnostics

### Very wide distance distribution

**Cause:** The selected atoms may be flexible or the selection is too broad.

**Fix:**
- Ensure selections resolve to single atoms (or use `midpoint()`/`com()`)
- Check that `selection1` and `selection2` are correctly specified
- Visualize the selections in a molecular viewer

### "Low statistical reliability" warning

**Cause:** Long correlation time relative to trajectory length.

**This is informational, not an error.** Results are still valid but uncertainties
may be underestimated.

**Mitigation:**
- Use multiple replicates (aggregated SEM is more reliable)
- Run longer simulations
- Results are still useful for qualitative comparisons

### Missing replicate data

**Message:** `Skipping replicate N: trajectory data not found`

**Cause:** The requested replicate hasn't completed or path is incorrect.

**Fix:** This is informational—analysis continues with available replicates.
Check simulation status if unexpected.

## Comparison with Catalytic Triad Analysis

Distance analysis and [Catalytic Triad Analysis](analysis_triad_quickstart.md)
both measure atom-pair distances, but serve different purposes:

| Feature | Distances | Catalytic Triad |
|---------|-----------|-----------------|
| **Focus** | Any atom pairs | Pre-defined catalytic geometry |
| **Configuration** | `comparison.yaml` (`plugins.distances`) | `comparison.yaml` with conditions |
| **Multi-condition** | Via `compare run distances` | Built-in condition comparison |
| **Simultaneous contacts** | Not computed | Key metric (all pairs < threshold) |
| **Use case** | Ad-hoc distance measurements | Structured enzyme comparisons |

```{tip}
Use **distances** for exploratory analysis of specific interactions.
Use **catalytic triad** when comparing enzyme integrity across conditions.
```

## Comparing Distances Across Conditions

To statistically compare distances across multiple simulation conditions (e.g.,
different polymer compositions), use the `compare run distances` command:

```bash
# Add distances section to comparison.yaml, then:
polyzymd compare run distances -f comparison.yaml
```

This provides:
- **Dual-metric ranking**: By mean distance (primary) and fraction below threshold (secondary)
- **Statistical tests**: t-tests, Cohen's d effect sizes, ANOVA
- **Per-pair summaries**: Distance statistics for each defined pair

See [Comparing Distances Across Conditions](analysis_compare_conditions.md)
for full documentation.

## Next Steps

- **Compare distances across conditions**: [Comparing Conditions Guide](analysis_compare_conditions.md)
- **Catalytic triad analysis**: [Triad Quick Start](analysis_triad_quickstart.md)
- **Understand statistics**: [Statistical Best Practices](../explanation/analysis_statistics_best_practices.md)
- **Contact analysis**: [Contacts Quick Start](analysis_contacts_quickstart.md)
