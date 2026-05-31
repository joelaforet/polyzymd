# Plugin settings in comparison.yaml

This reference lists the plugin-specific YAML settings accepted under the
top-level `plugins:` key in `comparison.yaml`. Use it to look up field names,
types, defaults, and meanings for each analysis plugin.

```yaml
plugins:
  rmsf:
    selection: "protein and name CA"
  contacts:
    cutoff: 4.5
```

FDR thresholds for comparison workflows are configured through top-level
`defaults.fdr_alpha` in `comparison.yaml` where supported, not through
plugin-local settings unless a plugin explicitly lists its own `fdr_alpha` field
below.

Related pages: {doc}`comparison_yaml`, {doc}`analysis_comparison_reference`, and
the per-plugin reference pages in this section.

## Plugin index

Structure and stability
: [RMSF](#rmsf-settings), [RMSD](#rmsd-settings), [Radius of gyration](#radius-of-gyration-settings), [Secondary structure](#secondary-structure-settings)

Interactions and contacts
: [Contacts](#contacts-settings), [Hydrogen bonds](#hydrogen-bonds-settings), [Distances](#distances-settings), [Catalytic triad](#catalytic-triad-settings)

Solvent exposure
: [SASA](#sasa-settings)

## RMSF settings

Path: `plugins.rmsf`

### Fields

#### `selection`

Type
: `str`

Default
: `"protein and name CA"`

MDAnalysis selection used for RMSF calculation.

#### `reference_mode`

Type
: `str`

Default
: `"centroid"`

Reference mode: `centroid`, `average`, `frame`, or `external`.

#### `reference_frame`

Type
: `int | null`

Default
: `null`

Frame number used when `reference_mode: frame` (1-indexed).

#### `reference_file`

Type
: `str | null`

Default
: `null`

External PDB path used when `reference_mode: external`.

#### `alignment_selection`

Type
: `str`

Default
: `"protein and name CA"`

Selection used for trajectory alignment.

#### `centroid_selection`

Type
: `str`

Default
: `"protein"`

Selection used to find centroid frame.

## Catalytic triad settings

Path: `plugins.catalytic_triad`

### Fields

#### `name`

Type
: `str`

Default
: `"catalytic_triad"`

Name of the triad/active-site definition.

#### `pairs`

Type
: `list[TriadPairSettings]`

Default
: required

Distance pairs to monitor.

#### `threshold`

Type
: `float`

Default
: `3.5`

Contact threshold in Å.

#### `description`

Type
: `str | null`

Default
: `null`

Optional human-readable description.

### Triad pair fields

Path: `plugins.catalytic_triad.pairs[]`

#### `label`

Type
: `str`

Default
: required

Human-readable pair name.

#### `selection_a`

Type
: `str`

Default
: required

First atom/point selection.

#### `selection_b`

Type
: `str`

Default
: required

Second atom/point selection.

## Distances settings

Path: `plugins.distances`

### Fields

#### `threshold`

Type
: `float | null`

Default
: `3.5`

Global distance threshold (Å) for contact-style state analysis.

#### `pairs`

Type
: `list[DistancePairSettings]`

Default
: `[]` (must be non-empty)

Distance pairs to monitor.

#### `use_pbc`

Type
: `bool`

Default
: `true`

Use minimum-image PBC-aware distances.

#### `align_trajectory`

Type
: `bool`

Default
: `true`

Align trajectory before distance calculations.

#### `alignment_selection`

Type
: `str`

Default
: `"protein and name CA"`

Selection for alignment.

#### `alignment_mode`

Type
: `str`

Default
: `"centroid"`

Alignment reference mode: `centroid`, `average`, or `frame`.

#### `alignment_frame`

Type
: `int | null`

Default
: `null`

Frame index (1-indexed) when `alignment_mode: frame`.

### Distance pair fields

Path: `plugins.distances.pairs[]`

#### `label`

Type
: `str`

Default
: required

Human-readable pair name.

#### `selection_a`

Type
: `str`

Default
: required

First atom/point selection.

#### `selection_b`

Type
: `str`

Default
: required

Second atom/point selection.

#### `threshold`

Type
: `float | null`

Default
: `null`

Per-pair threshold override; falls back to global `threshold`.

#### `below_label`

Type
: `str | null`

Default
: `null`

Display label for below-threshold state.

#### `above_label`

Type
: `str | null`

Default
: `null`

Display label for above-threshold state.

## Contacts settings

Path: `plugins.contacts`

### Fields

#### `polymer_selection`

Type
: `str`

Default
: `"chainid C"`

MDAnalysis selection for polymer atoms.

#### `protein_selection`

Type
: `str`

Default
: `"chainid A"`

MDAnalysis selection for protein atoms.

#### `cutoff`

Type
: `float`

Default
: `4.5`

Contact cutoff distance in Å.

#### `polymer_types`

Type
: `list[str] | null`

Default
: `null`

Optional polymer residue-name filter.

#### `grouping`

Type
: `str`

Default
: `"aa_class"`

Grouping mode: `aa_class`, `secondary_structure`, or `none`.

#### `compute_residence_times`

Type
: `bool`

Default
: `true`

Compute aggregate residence-time summaries and plots; per-replicate contact
events remain stored when disabled.

#### `protein_groups`

Type
: `dict[str, list[int]] | null`

Default
: `null`

Custom residue groups, for example `{name: [resid, ...]}`.

#### `protein_partitions`

Type
: `dict[str, list[str]] | null`

Default
: `null`

Partition definitions over custom `protein_groups`.

#### `fdr_alpha`

Type
: `float`

Default
: `0.05`

Benjamini-Hochberg false-discovery-rate alpha.

#### `min_effect_size`

Type
: `float`

Default
: `0.5`

Minimum Cohen's d highlighted in output.

#### `top_residues`

Type
: `int`

Default
: `10`

Number of top-contact residues shown in summaries.

## Secondary structure settings

Path: `plugins.secondary_structure`

### Fields

#### `chain_id`

Type
: `str`

Default
: `"A"`

Protein chain letter to analyze (PolyzyMD convention: chain A).

## SASA settings

Path: `plugins.sasa`

### Fields

#### `runs`

Type
: `list[SASARunSettings]`

Default
: `[]` (must be non-empty)

SASA runs to compute.

#### `probe_radius_nm`

Type
: `float`

Default
: `0.14`

MDTraj Shrake-Rupley probe radius (nm).

#### `n_sphere_points`

Type
: `int`

Default
: `960`

MDTraj Shrake-Rupley sphere point count.

#### `chunk_size`

Type
: `int`

Default
: `100`

Frames per chunk for SASA computation.

### SASA run fields

Path: `plugins.sasa.runs[]`

#### `label`

Type
: `str`

Default
: required

Human-readable run label.

#### `target_selection`

Type
: `str`

Default
: required

Selection whose SASA is reported.

#### `context_selection`

Type
: `str | null`

Default
: `null`

Environment/context selection for SASA computation. `null` defaults to
`target_selection`.

#### `stride`

Type
: `int`

Default
: `1`

Frame stride; `1` means every frame.

## Hydrogen bonds settings

Path: `plugins.hydrogen_bonds`

### Fields

#### `groups`

Type
: `dict[str, str]`

Default
: `{protein: "chainid A", polymer: "chainid C"}`

Named MDAnalysis selections used by summaries.

#### `summaries`

Type
: `list[HydrogenBondSummarySettings]`

Default
: `[{name: protein_polymer, between: [protein, polymer]}]`

Summary definitions to compute.

#### `distance_cutoff`

Type
: `float`

Default
: `3.0`

Donor-acceptor cutoff (Å).

#### `angle_cutoff`

Type
: `float`

Default
: `150.0`

D-H...A angle cutoff (degrees).

#### `update_selections`

Type
: `bool`

Default
: `true`

Re-evaluate selections each frame.

#### `top_n_pairs`

Type
: `int`

Default
: `15`

Number of top residue pairs reported.

#### `allow_empty_groups`

Type
: `bool`

Default
: `true`

Warn/skip empty groups instead of raising.

#### `allow_overlapping_composition`

Type
: `bool`

Default
: `false`

Allow overlapping composition partitions; otherwise raise.

#### `composition`

Type
: `HydrogenBondCompositionSettings | null`

Default
: `null`

Optional partitioning for composition analysis.

#### `timestep_ps`

Type
: `float | null`

Default
: `null`

Optional timestep override (ps) for time-axis plots.

Time-axis plots assume uniformly saved frames. PolyzyMD maps frame index to time
as `frame_index * timestep_ps`; variable-timestep concatenated trajectories are
not supported.

Hydrogen detection uses MDAnalysis `HydrogenBondAnalysis` with hydrogens
selected as `(<group union>) and element H`; explicit hydrogens and reliable
element metadata are required for meaningful counts.

### Hydrogen-bond summary fields

Path: `plugins.hydrogen_bonds.summaries[]`

Exactly one of `between` or `within` must be set for each summary.

#### `name`

Type
: `str`

Default
: required

Unique summary name.

#### `between`

Type
: `tuple[str, str] | null`

Default
: `null`

Cross-group summary mode.

#### `within`

Type
: `str | null`

Default
: `null`

Intra-group summary mode.

### Hydrogen-bond composition fields

Path: `plugins.hydrogen_bonds.composition`

#### `partitions`

Type
: `dict[str, str]`

Default
: `{}`

Named composition partitions as MDAnalysis selections.

## Radius of gyration settings

Path: `plugins.rg`

### Fields

#### `runs`

Type
: `list[RgRunSettings]`

Default
: `[]` (must be non-empty)

Named Rg runs to compute.

### Radius-of-gyration run fields

Path: `plugins.rg.runs[]`

#### `label`

Type
: `str`

Default
: required

Human-readable run label.

#### `selection`

Type
: `str`

Default
: required

MDAnalysis selection for Rg calculation.

#### `calculation_mode`

Type
: `"selection" | "fragments"`

Default
: `"selection"`

Whole-selection vs fragment-reduced Rg mode.

#### `fragment_weighting`

Type
: `"equal" | "mass"`

Default
: `"equal"`

Fragment reduction weighting (fragment mode).

#### `save_fragment_distribution`

Type
: `bool`

Default
: `true`

Save per-fragment distribution sidecar outputs.

#### `histogram_bins`

Type
: `int`

Default
: `50`

Histogram bins for fragment distribution summaries.

## RMSD settings

Path: `plugins.rmsd`

### Fields

#### `runs`

Type
: `list[RMSDRunSettings]`

Default
: `[]` (must be non-empty)

Named RMSD runs to compute.

### RMSD run fields

Path: `plugins.rmsd.runs[]`

#### `label`

Type
: `str`

Default
: required

Human-readable run label.

#### `selection`

Type
: `str`

Default
: `"protein and name CA"`

Selection used for RMSD calculation.

#### `alignment_selection`

Type
: `str`

Default
: `"protein and name CA"`

Selection used for alignment.

#### `reference_mode`

Type
: `str`

Default
: `"centroid"`

Reference mode: `centroid`, `average`, `frame`, or `external`.

#### `reference_frame`

Type
: `int`

Default
: `0`

0-indexed frame index for `reference_mode: frame`.

#### `reference_file`

Type
: `str | null`

Default
: `null`

External PDB path for `reference_mode: external`.

#### `centroid_selection`

Type
: `str | null`

Default
: `null`

Optional centroid-mode selection; falls back to `alignment_selection`.

#### `convergence_window_size_ns`

Type
: `float`

Default
: `15.0`

Window size for convergence detection.

#### `convergence_step_size_ns`

Type
: `float`

Default
: `5.0`

Step size for sliding convergence windows.

#### `convergence_slope_threshold`

Type
: `float`

Default
: `0.0005`

Maximum slope threshold in Å/ns for a window to be considered converged.

#### `convergence_sustained_for_ns`

Type
: `float`

Default
: `15.0`

Sustained duration required before declaring convergence.
