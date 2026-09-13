# Catalytic triad plugin reference

For a step-by-step workflow, see
{doc}`../how_to/analysis_triad_quickstart`.

## Settings

```yaml
plugins:
  catalytic_triad:
    name: "LipA Catalytic Triad"
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

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `pairs` | `list[PairSelection]` | *required* | Pairs to monitor, at least one |
| `threshold` | `float` | `3.5` | Contact cutoff in angstrom, applied to every pair |
| `name` | `str` | `"catalytic_triad"` | Name of the active site |
| `description` | `str \| null` | `null` | What the active site is |

Each entry in `pairs`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Name the pair is reported under |
| `selection_a` | `str` | *required* | First endpoint selection |
| `selection_b` | `str` | *required* | Second endpoint selection |

An endpoint is one point, so a multi-atom selection must be wrapped in
`midpoint(...)` or `com(...)`. The syntax is the same as for
{doc}`analysis_distances_reference`.

```{important}
Use chain-aware selections. Residue IDs restart per chain in PolyzyMD systems,
so a bare `resid X` can match non-protein atoms. Prefer `protein and resid X`.
```

## Observables

A triad of *n* pairs produces `2n + 1` observables.

| Observable | Kind | Unit | Meaning |
|------------|------|------|---------|
| `<label>` | `mean_of_timeseries` | `A` | Distance for that pair, per frame |
| `<label> within <threshold> A` | `fraction` | `fraction` | Frames in which that pair is strictly within the cutoff |
| `simultaneous_contact_fraction` | `fraction` | `fraction` | Frames in which every pair is within the cutoff at once |

The charge-relay system only works when every link is short at the same time, so
the simultaneous fraction is the conjunction of the per-pair indicators and not
the average of the per-pair fractions. It is stored as a fraction in `[0, 1]`; a
report that wants a percentage multiplies at display time. A pair counts as
within the cutoff when its distance is strictly less than it, and every fraction
records that as `threshold_operator: strict_less_than` in its metadata, beside
the name of the active site and the selections it was measured from.

The framework reduces each replicate to one number per observable, reports the
mean over replicates with its SEM and a Student t interval, and tests conditions
against the control on replicate-level values. The per-frame series of every
observable is kept in the replicate's `observables.npz` sidecar.

The pair distances and the simultaneous fraction are tested; the per-pair
contact fractions are not. A per-pair fraction is a monotone function of the
same series as that pair's mean distance, so it would add no information to the
Benjamini-Hochberg family that every test in the run shares, while the
simultaneous fraction says something no single pair does. The untested
observables still carry their mean, SEM and interval.

## Periodic boundaries

Distances are measured on the coordinates as the trajectory stores them, with
the minimum image convention applied against the box of each frame, read from
`Timestep.dimensions`. The `pbc_policy` chosen at load time, the trajectory
variant and the bond source are recorded in provenance; see
{doc}`analysis_plugin_settings`.

```{versionchanged} 1.3.0
The triad no longer aligns the trajectory before measuring. Distances are
invariant under rigid-body motion, so alignment changed nothing that was
correct, while the in-memory alignment rotated coordinates without rotating the
box vectors, which corrupted minimum-image distances for pairs separated by more
than half a box length. For the reasoning, see
{doc}`../explanation/analysis_concepts`.
```

```{versionchanged} 1.3.0
The plugin is written against the observable contract. Contact fractions are
stored as fractions rather than percentages, the payload holds observables
rather than the former `metrics` and `pair_results` records, the KDE panel is
replaced by the shared contract figures, and a `plot_settings.catalytic_triad`
block is accepted for one release and ignored.
```

## Output files

```text
comparison_workspace/
├── analysis/
│   └── <condition>/
│       └── catalytic_triad/
│           ├── run_1/
│           │   ├── result.json
│           │   └── observables.npz
│           ├── run_2/
│           ├── run_3/
│           └── aggregated/result.json
└── comparison/
    └── catalytic_triad/result.json
```

| Level | Artifact | Path |
|-------|----------|------|
| Per replicate | `ReplicateArtifact` | `analysis/<condition>/catalytic_triad/run_<replicate>/result.json` |
| Per condition | `ConditionArtifact` | `analysis/<condition>/catalytic_triad/aggregated/result.json` |
| Cross condition | `ComparisonArtifact` | `comparison/catalytic_triad/result.json` |
| Per-frame series | NPZ sidecar | `analysis/<condition>/catalytic_triad/run_<replicate>/observables.npz` |

Read saved artifacts through the public store:

```python
from pathlib import Path

from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.mda import ArtifactStore

condition = ArtifactStore(Path("analysis/PEGylated/catalytic_triad/aggregated"))
for payload in condition.read_condition_result().payload["observables"]:
    aggregate = ObservableAggregate.model_validate(payload)
    print(aggregate.name, aggregate.mean, aggregate.sem, aggregate.unit)
```

## Commands

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Path to comparison config |
| `--eq-time` | from YAML defaults | Equilibration time to skip |
| `--recompute` | off | Ignore cached replicates and measure again |
| `--format` | `table` | Output format, `table` or `json` |

```bash
polyzymd compare run catalytic_triad -f comparison.yaml --eq-time 200ns
```

## Troubleshooting

### "selection ... matched no atoms"

The residue ID or atom name does not match the loaded topology. Check the
numbering, check the atom names (`OD1`/`OD2`, `ND1`/`NE2`, `OG`), and prefer
chain-aware selections such as `protein and resid 132`.

### "selection ... matched N atoms, and a pair endpoint is one point"

Wrap the selection in `midpoint(...)` or `com(...)`.

### Distances above 10 angstrom

Usually a wrong selection or mismatched residue numbering. Confirm the residues
in a viewer before reading the result as a disrupted triad.

### A simultaneous contact fraction near zero

Read the per-pair fractions first: one limiting pair distinguishes a genuinely
disrupted relay from a cutoff that is too strict for this system.
