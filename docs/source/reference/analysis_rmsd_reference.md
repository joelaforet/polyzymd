# RMSD plugin reference

For a step-by-step guide to running RMSD analysis, see
{doc}`../how_to/analysis_rmsd_quickstart`.

The plugin measures the deviation of a selection from a reference structure,
frame by frame, in angstrom. It reports one observable of kind
`mean_of_timeseries` per configured run and leaves aggregation, uncertainty,
cross-condition tests, storage and formatting to the analysis framework.

## Configuration reference

Top-level `RMSDSettings`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | `list[RMSDRunSettings]` | required, at least one | Named RMSD runs to compute |

All fields of `RMSDRunSettings`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Run label, unique within the file |
| `selection` | `str` | `"protein and name CA"` | MDAnalysis selection whose deviation is measured |
| `alignment_selection` | `str` | `"protein and name CA"` | MDAnalysis selection that superposition minimises over; see below when it differs from `selection` |
| `reference_mode` | `str` | `"centroid"` | `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | `int` | `0` | 0-indexed frame used when `reference_mode: frame` |
| `reference_file` | `str \| None` | `null` | External PDB used when `reference_mode: external` |
| `centroid_selection` | `str \| None` | `null` | Selection used to find the representative frame; defaults to `alignment_selection` |

```{note}
Run labels must be unique within one `comparison.yaml`. A duplicate label
raises a validation error, because two runs would otherwise write one
observable name.
```

### Settings removed in 1.3.0

`convergence_window_size_ns`, `convergence_step_size_ns`,
`convergence_slope_threshold` and `convergence_sustained_for_ns` no longer do
anything. They still parse for one release and raise a `DeprecationWarning`, so
an existing `comparison.yaml` keeps working, and they will be rejected in the
release after that. The sliding-window convergence flag they configured was
removed because its default slope threshold sat below the scatter of successive
window means, so the flag tracked noise.

## Observable names

Each run reports one observable named `rmsd_<label slug>_ref_<reference_mode>`,
where the slug is the label lowercased with every run of non-alphanumeric
characters replaced by an underscore. A run labelled `Protein Backbone` with
`reference_mode: centroid` reports `rmsd_protein_backbone_ref_centroid`.

The reference mode is part of the name on purpose. A centroid or average
reference differs between replicates and between conditions, so such a run
measures spread within a replicate, not deviation from a shared structure. Only
`external` mode measures every condition against one structure. Keeping the
mode in the name stops a comparison table from presenting the two as one
quantity.

| Mode | Reference structure |
|------|---------------------|
| `centroid` | The frame of the production window closest to the aligned mean |
| `average` | The mean position of each selected atom over the aligned window |
| `frame` | The frame named by `reference_frame` |
| `external` | The structure in `reference_file` |

## What `alignment_selection` does

The measurement superimposes each frame on the reference itself, so the plugin
runs no separate alignment pass over the trajectory. `alignment_selection`
names the atoms that superposition minimises over.

When `alignment_selection` equals `selection`, the reported value is the
minimised RMSD of those atoms, which is the usual global stability number.

When the two differ, each frame is superimposed on `alignment_selection` and
the deviation is then reported for `selection`. That is how a loop, a lid or a
bound ligand is measured against a rigid core. Superimpose on the core, then
ask how far the other group moved. The two settings answer different questions,
and the second is not a refinement of the first. Superimposing a group on itself always
hides its own displacement, so a run whose `selection` is a small flexible
group and whose `alignment_selection` is the same group reports almost nothing.

A superposition group of fewer than three atoms leaves the rotation
undetermined, and the plugin raises `SelectionError` rather than reporting the
NaN that MDAnalysis returns.

`average` mode is the one mode that superimposes the trajectory in place,
because a mean structure has no meaning until every frame shares a frame of
reference. The value it then reports is unaffected by that pass, for the same
reason as above.

## Output files

Results are canonical v1.3 artifacts. The JSON files carry the reduced
observables; the full per-frame series is an NPZ sidecar.

```text
<comparison_workspace>/
├── analysis/
│   └── <condition>/
│       └── rmsd/
│           ├── run_1/
│           │   ├── result.json
│           │   └── observables.npz
│           ├── run_2/
│           ├── run_3/
│           └── aggregated/
│               └── result.json
└── comparison/
    └── rmsd/
        └── result.json
```

| Level | Artifact | Path |
|-------|----------|------|
| Per replicate | `ReplicateArtifact` | `analysis/<condition>/rmsd/run_<replicate>/result.json` |
| Per-frame series | NPZ sidecar | `analysis/<condition>/rmsd/run_<replicate>/observables.npz` |
| Per condition | `ConditionArtifact` | `analysis/<condition>/rmsd/aggregated/result.json` |
| Cross condition | `ComparisonArtifact` | `comparison/rmsd/result.json` |

The NPZ sidecar holds one array per observable, keyed by observable name, with
one value per analysed frame.

### Replicate payload

`payload["observables"]` is a list of `ObservableEstimate` records:

| Field | Description |
|-------|-------------|
| `name` | Observable name, reference mode included |
| `kind` | `"mean_of_timeseries"` |
| `unit` | `"A"` |
| `value` | Mean RMSD of this replicate |
| `n_frames` | Frames in the production window |
| `statistical_inefficiency` | Correlation diagnostic g, or `null` for a short or constant series |
| `n_eff` | Effective sample size within the replicate, reported but never used to shrink an error bar |
| `higher_is_better` | `false` |

`provenance["identity"]` carries the polyzymd version, the plugin source hash,
the settings fingerprint, the config hash, the equilibration setting and the
identity of every input file. A replicate is recomputed when any of these
changes and reused when none does.

### Condition payload

`payload["observables"]` is a list of `ObservableAggregate` records, one per
observable, computed across replicates:

| Field | Description |
|-------|-------------|
| `replicate_values` | One mean RMSD per replicate, the sample every statistic is computed from |
| `mean`, `sem` | Mean and standard error across replicates, `ddof = 1` |
| `ci95_low`, `ci95_high`, `ci_method`, `coverage` | Student t interval on the mean, `null` for a single replicate |
| `n_replicates` | Size of the replicate sample |
| `n_eff_min` | Smallest effective sample size among the replicates |

### Comparison payload

`payload["comparisons"]` holds one `ObservableComparison` per observable and
non-control condition, with `delta`, `percent_change`, `p_value`, `p_adjusted`,
`correction`, `cohens_d`, `significant`, `testable` and `note`. Tests run on
replicate-level values, and every test in the run forms one
Benjamini-Hochberg family unless the config asks for Tukey's test with three or
more conditions.

Read an artifact with `ArtifactStore`:

```python
from pathlib import Path

from polyzymd.analyses.mda import ArtifactStore

replicate = ArtifactStore(Path("analysis/PEGylated/rmsd/run_1")).read_replicate_result()
condition = ArtifactStore(Path("analysis/PEGylated/rmsd/aggregated")).read_condition_result()
print(replicate.payload["observables"][0]["value"])
print(condition.payload["observables"][0]["mean"])
```

## Figures

Figures come from the framework, keyed on the observable kind, not from the
plugin, so there are no rmsd-specific plot settings any more. A
`mean_of_timeseries` observable gives a comparison bar chart with the replicate
points overlaid and a per-replicate time series panel, both footnoted with the
interval, the replicate count and the production window. Until the generic
figures land, `polyzymd compare run-all --plot` writes no figure for rmsd; the
per-frame series in the NPZ sidecar is the input for a plot of your own. A
`plot_settings.rmsd` block in an existing comparison file still loads and warns
that it does nothing.

## Common CLI options

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Path to comparison configuration |
| `--eq-time` | none, falls back to `defaults.equilibration_time` in the file | Override the equilibration time to skip |
| `--recompute` | off | Ignore cached results and recompute |
| `--format` | `table` | Output format (`table`, `markdown`, `json` or `agent`) |
| `-o, --output` | (none) | Save formatted output to file |
| `-q, --quiet` | off | Suppress INFO messages |
| `--debug` | off | Enable DEBUG logging |

## Troubleshooting

### "selection ... matched no atoms"

The MDAnalysis selection matches nothing in the topology. Check residue
numbering in the PDB against MDAnalysis, verify atom names, and rerun with
`polyzymd --debug compare run rmsd -f comparison.yaml` for the full context.
An empty selection raises `SelectionError`; it never yields an RMSD of zero.

### "reference_file ... does not exist on this machine"

External mode was requested without a readable PDB. Give an absolute path, or a
path relative to the working directory. A path that does not exist where the
settings are parsed is a warning, because a comparison file often names a
cluster path and is validated on a laptop; it becomes a `ReplicateError` on the
machine that reads the trajectory.

### "selection ... matches N atoms in the trajectory but M in ..."

The selection string picks different atom counts in the trajectory and in the
external reference. Make the atom naming consistent, check that the external
PDB covers the same residues, or narrow the selection.

### "rmsd needs a contiguous production window"

The frame selection listed explicit frame indices. The centroid and average
references are defined on a window, so the plugin needs a start, stop and step.
Use an equilibration time rather than an explicit frame list.

### Very high RMSD values, above 10 Å

Usually the alignment selection, the reference mode, or the system itself.
Check that `alignment_selection` matches atoms in the system, compare against
`reference_mode: "average"`, verify the trajectory files are complete, and
inspect the structures for unfolding or a large conformational change.

## RMSD compared with RMSF

| Feature | RMSD | RMSF |
|---------|------|------|
| Measures | Global deviation from a reference | Per-residue fluctuation |
| Output | One value per frame | One value per residue |
| Reference | Fixed structure, centroid, average, frame or external | Time-averaged position |
| Detects | Conformational drift, unfolding | Flexible loops, rigid core |
| Multi-run | Yes, a `runs` list with different selections | Single selection |
| Best for | Stability comparison, equilibration assessment | Flexibility mapping |

```{tip}
Use RMSD first to judge overall stability and choose an equilibration time,
then use RMSF to find which regions drive a flexibility difference.
```
