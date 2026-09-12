# Distances plugin reference

For a step-by-step task guide, see
{doc}`../how_to/analysis_distances_quickstart`.

## Settings

All fields for `plugins.distances`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `pairs` | `list[DistancePair]` | *required* | One or more named pairs, at least one |
| `threshold` | `float \| null` | `3.5` | Threshold in angstrom for pairs that do not set their own. `null` reports distances only |
| `use_pbc` | `bool` | `true` | Take minimum-image distances against the box of the measured frame |

Each entry in `pairs`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Name the pair is reported under |
| `selection_a` | `str` | *required* | First endpoint selection |
| `selection_b` | `str` | *required* | Second endpoint selection |
| `threshold` | `float \| null` | global `threshold` | Threshold for this pair |
| `below_label` | `str \| null` | `"below <threshold> A"` | Name of the below-threshold state |

An endpoint is one point. A selection that matches several atoms is rejected
unless it is wrapped in `midpoint(...)` or `com(...)`, and a selection that
matches no atoms raises `SelectionError` with the topology diagnostics.

### Selection syntax

| Syntax | Meaning | Typical use |
|--------|---------|-------------|
| `midpoint(selection)` | Centre of geometry of the selected atoms | Carboxylate oxygens of Asp or Glu |
| `com(selection)` | Centre of mass of the selected atoms | A whole residue, ligand or domain |
| `pdbindex N` | Atom by PDB serial number, 1-indexed | Copying atom IDs from PyMOL |

```yaml
plugins:
  distances:
    threshold: 3.5
    pairs:
      - label: "His156-Asp133"
        selection_a: "protein and resid 156 and name ND1"
        selection_b: "midpoint(protein and resid 133 and name OD1 OD2)"
      - label: "Ser77(OG)-Substrate"
        selection_a: "protein and resid 76 and name OG"
        selection_b: "resname RBY and name C13x"
        threshold: 10.0
        below_label: "Within 10 Angstrom"
```

```{important}
Residue indices restart by chain in PolyzyMD systems. For protein residues,
prefer `protein and resid ...` to avoid an accidental multi-chain match.
```

## Observables

One pair produces one or two observables.

| Observable | Kind | Unit | Meaning |
|------------|------|------|---------|
| `<label>` | `mean_of_timeseries` | `A` | Distance between the two endpoints, per frame |
| `<label> <below_label>` | `fraction` | `fraction` | Frames in which the distance is strictly below the threshold |

The framework reduces each replicate to one number per observable, then reports
the mean over replicates with its SEM and a Student t interval, and tests
conditions against the control on replicate-level values. Nothing averages one
pair into another. The per-frame series of every observable is kept in the
replicate's `observables.npz` sidecar.

## Periodic boundaries

With `use_pbc: true` the distance uses the minimum image convention against the
box stored in the frame being measured, read from `Timestep.dimensions`. A frame
whose box is missing or degenerate is measured without periodicity, and the run
logs one warning saying so.

```{versionchanged} 1.3.0
Distances are measured without alignment. `align_trajectory` and the
`alignment_*` fields are still accepted so existing `comparison.yaml` files keep
loading, but they change nothing and raise a `DeprecationWarning`. A distance is
invariant under rigid-body motion, so alignment could never improve it, while
aligning in memory rotated the coordinates without rotating the box vectors and
so corrupted minimum-image distances for pairs separated by more than half a box
length. For the reasoning, see {doc}`../explanation/analysis_concepts`.
```

```{versionchanged} 1.3.0
The plugin is written against the observable contract. The per-pair `above_label`
field and the `plot_settings.distances` block are accepted for one release and
ignored, the KDE distribution figures are replaced by the shared contract
figures, and the per-replicate and aggregated files hold observables rather than
the former `pair_results` records.
```

## Output files

```text
<projects_directory>/
└── analysis/
    └── distances/
        ├── run_1/
        │   ├── replicate.json
        │   └── observables.npz
        ├── run_2/
        ├── run_3/
        └── aggregated/
            └── condition.json
```

Each `replicate.json` holds one `ObservableEstimate` per observable, with its
kind, unit, frame count, statistical inefficiency and effective sample size. The
condition artifact holds one `ObservableAggregate` per observable, with the
replicate values, the mean, the SEM, the interval and its coverage.

## Commands

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Comparison config path |
| `--eq-time` | `0ns` | Equilibration time to discard |
| `--recompute` | off | Ignore cached replicates and measure again |
| `--format` | `table` | Output format, `table` or `json` |
| `-o, --output` | (none) | Write formatted output to a file |

```bash
polyzymd compare run distances -f comparison.yaml --eq-time 200ns
polyzymd compare run distances -f comparison.yaml --eq-time 200ns --format json
```

## Troubleshooting

### "selection ... matched no atoms"

The selection does not match the topology. Add chain-aware qualifiers such as
`protein and resid ...`, check atom names and residue IDs against the topology,
and read the diagnostics block in the error, which lists what the topology does
contain.

### "selection ... matched N atoms, and a pair endpoint is one point"

Wrap the selection in `midpoint(...)` or `com(...)` to say which point of the
group you mean.

### Long-distance outliers near the box boundary

Check that `use_pbc` is true, and check whether a measured group can straddle a
periodic boundary. The centre of mass of a split molecule is meaningless, and
only a whole-molecule trajectory or a `make_whole` load fixes that.
