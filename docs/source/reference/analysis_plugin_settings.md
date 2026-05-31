# Analysis Plugin Settings Reference

This page is a complete YAML key reference for plugin settings under
`comparison.yaml`:

```yaml
plugins:
  <plugin_name>:
    ...settings...
```

Use this as a lookup table for field names, types, defaults, and meanings.

FDR thresholds for comparison workflows are configured through top-level
`defaults.fdr_alpha` in `comparison.yaml` where supported, not through
plugin-local settings unless a plugin explicitly lists its own `fdr_alpha` field
below.

## `rmsf`

| Key | Type | Default | Description |
|---|---|---|---|
| `selection` | `str` | `"protein and name CA"` | MDAnalysis selection used for RMSF calculation |
| `reference_mode` | `str` | `"centroid"` | Reference mode: `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | `int \| null` | `null` | Frame number used when `reference_mode: frame` (1-indexed) |
| `reference_file` | `str \| null` | `null` | External PDB path used when `reference_mode: external` |
| `alignment_selection` | `str` | `"protein and name CA"` | Selection used for trajectory alignment |
| `centroid_selection` | `str` | `"protein"` | Selection used to find centroid frame |

## `catalytic_triad`

| Key | Type | Default | Description |
|---|---|---|---|
| `name` | `str` | `"catalytic_triad"` | Name of the triad/active-site definition |
| `pairs` | `list[TriadPairSettings]` | required | Distance pairs to monitor |
| `threshold` | `float` | `3.5` | Contact threshold in Å |
| `description` | `str \| null` | `null` | Optional human-readable description |

`TriadPairSettings` entries in `pairs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable pair name |
| `selection_a` | `str` | required | First atom/point selection |
| `selection_b` | `str` | required | Second atom/point selection |

## `distances`

| Key | Type | Default | Description |
|---|---|---|---|
| `threshold` | `float \| null` | `3.5` | Global distance threshold (Å) for contact-style state analysis |
| `pairs` | `list[DistancePairSettings]` | `[]` (must be non-empty) | Distance pairs to monitor |
| `use_pbc` | `bool` | `true` | Use minimum-image PBC-aware distances |
| `align_trajectory` | `bool` | `true` | Align trajectory before distance calculations |
| `alignment_selection` | `str` | `"protein and name CA"` | Selection for alignment |
| `alignment_mode` | `str` | `"centroid"` | Alignment reference mode: `centroid`, `average`, or `frame` |
| `alignment_frame` | `int \| null` | `null` | Frame index (1-indexed) when `alignment_mode: frame` |

`DistancePairSettings` entries in `pairs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable pair name |
| `selection_a` | `str` | required | First atom/point selection |
| `selection_b` | `str` | required | Second atom/point selection |
| `threshold` | `float \| null` | `null` | Per-pair threshold override (falls back to global `threshold`) |
| `below_label` | `str \| null` | `null` | Display label for below-threshold state |
| `above_label` | `str \| null` | `null` | Display label for above-threshold state |

## `contacts`

| Key | Type | Default | Description |
|---|---|---|---|
| `polymer_selection` | `str` | `"chainid C"` | MDAnalysis selection for polymer atoms |
| `protein_selection` | `str` | `"chainid A"` | MDAnalysis selection for protein atoms |
| `cutoff` | `float` | `4.5` | Contact cutoff distance in Å |
| `polymer_types` | `list[str] \| null` | `null` | Optional polymer residue-name filter |
| `grouping` | `str` | `"aa_class"` | Grouping mode: `aa_class`, `secondary_structure`, or `none` |
| `compute_residence_times` | `bool` | `true` | Compute aggregate residence-time summaries and plots; per-replicate contact events remain stored when disabled |
| `protein_groups` | `dict[str, list[int]] \| null` | `null` | Custom residue groups, e.g. `{name: [resid, ...]}` |
| `protein_partitions` | `dict[str, list[str]] \| null` | `null` | Partition definitions over custom `protein_groups` |
| `fdr_alpha` | `float` | `0.05` | Benjamini-Hochberg false-discovery-rate alpha |
| `min_effect_size` | `float` | `0.5` | Minimum Cohen's d highlighted in output |
| `top_residues` | `int` | `10` | Number of top-contact residues shown in summaries |

## `secondary_structure`

| Key | Type | Default | Description |
|---|---|---|---|
| `chain_id` | `str` | `"A"` | Protein chain letter to analyze (PolyzyMD convention: chain A) |

## `sasa`

| Key | Type | Default | Description |
|---|---|---|---|
| `runs` | `list[SASARunSettings]` | `[]` (must be non-empty) | SASA runs to compute |
| `probe_radius_nm` | `float` | `0.14` | MDTraj Shrake-Rupley probe radius (nm) |
| `n_sphere_points` | `int` | `960` | MDTraj Shrake-Rupley sphere point count |
| `chunk_size` | `int` | `100` | Frames per chunk for SASA computation |

`SASARunSettings` entries in `runs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable run label |
| `target_selection` | `str` | required | Selection whose SASA is reported |
| `context_selection` | `str \| null` | `null` | Environment/context selection for SASA computation (`null` defaults to `target_selection`) |
| `stride` | `int` | `1` | Frame stride (1 means every frame) |

## `hydrogen_bonds`

| Key | Type | Default | Description |
|---|---|---|---|
| `groups` | `dict[str, str]` | `{protein: "chainid A", polymer: "chainid C"}` | Named MDAnalysis selections used by summaries |
| `summaries` | `list[HydrogenBondSummarySettings]` | `[{name: protein_polymer, between: [protein, polymer]}]` | Summary definitions to compute |
| `distance_cutoff` | `float` | `3.0` | Donor-acceptor cutoff (Å) |
| `angle_cutoff` | `float` | `150.0` | D-H...A angle cutoff (degrees) |
| `update_selections` | `bool` | `true` | Re-evaluate selections each frame |
| `top_n_pairs` | `int` | `15` | Number of top residue pairs reported |
| `allow_empty_groups` | `bool` | `true` | Warn/skip empty groups instead of raising |
| `allow_overlapping_composition` | `bool` | `false` | Allow overlapping composition partitions (otherwise raise) |
| `composition` | `HydrogenBondCompositionSettings \| null` | `null` | Optional partitioning for composition analysis |
| `timestep_ps` | `float \| null` | `null` | Optional timestep override (ps) for time-axis plots |

Time-axis plots assume uniformly saved frames. PolyzyMD maps frame index to time
as `frame_index * timestep_ps`; variable-timestep concatenated trajectories are
not supported.

`HydrogenBondSummarySettings` entries in `summaries`:

| Key | Type | Default | Description |
|---|---|---|---|
| `name` | `str` | required | Unique summary name |
| `between` | `tuple[str, str] \| null` | `null` | Cross-group summary mode |
| `within` | `str \| null` | `null` | Intra-group summary mode |

Exactly one of `between` or `within` must be set for each summary.

Hydrogen detection uses MDAnalysis `HydrogenBondAnalysis` with hydrogens
selected as `(<group union>) and element H`; explicit hydrogens and reliable
element metadata are required for meaningful counts.

`HydrogenBondCompositionSettings`:

| Key | Type | Default | Description |
|---|---|---|---|
| `partitions` | `dict[str, str]` | `{}` | Named composition partitions as MDAnalysis selections |

## `rg`

| Key | Type | Default | Description |
|---|---|---|---|
| `runs` | `list[RgRunSettings]` | `[]` (must be non-empty) | Named Rg runs to compute |

`RgRunSettings` entries in `runs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable run label |
| `selection` | `str` | required | MDAnalysis selection for Rg calculation |
| `calculation_mode` | `"selection" \| "fragments"` | `"selection"` | Whole-selection vs fragment-reduced Rg mode |
| `fragment_weighting` | `"equal" \| "mass"` | `"equal"` | Fragment reduction weighting (fragment mode) |
| `save_fragment_distribution` | `bool` | `true` | Save per-fragment distribution sidecar outputs |
| `histogram_bins` | `int` | `50` | Histogram bins for fragment distribution summaries |

## `rmsd`

| Key | Type | Default | Description |
|---|---|---|---|
| `runs` | `list[RMSDRunSettings]` | `[]` (must be non-empty) | Named RMSD runs to compute |

`RMSDRunSettings` entries in `runs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable run label |
| `selection` | `str` | `"protein and name CA"` | Selection used for RMSD calculation |
| `alignment_selection` | `str` | `"protein and name CA"` | Selection used for alignment |
| `reference_mode` | `str` | `"centroid"` | Reference mode: `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | `int` | `0` | 0-indexed frame index for `reference_mode: frame` |
| `reference_file` | `str \| null` | `null` | External PDB path for `reference_mode: external` |
| `centroid_selection` | `str \| null` | `null` | Optional centroid-mode selection (falls back to `alignment_selection`) |
