# `comparison.yaml` Schema Reference

The `comparison.yaml` file defines a cross-condition analysis project. It
specifies which simulation conditions to compare, which analysis plugins to
run, and how to visualize results. Create one with `polyzymd compare init -n
<name>` and place it at the root of your comparison project directory.

Source of truth: {func}`polyzymd.config.comparison.ComparisonConfig` in
`src/polyzymd/config/comparison.py`.

For CLI commands that consume this file, see
{doc}`analysis_comparison_reference`. For directory layout and data
expectations, see {doc}`data_requirements`.

---

## Minimal Working Example

```yaml
name: "polymer_stability_study"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]
  - label: "100% SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    selection: "protein and name CA"
```

---

## Top-Level Fields

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `name` | string | **yes** | — | Human-readable project name |
| `description` | string | no | `null` | Description of what is being compared |
| `control` | string | no | `null` | Label of the control condition. Must match a `label` in `conditions`. Used for relative comparisons (e.g., Δ from control). |
| `conditions` | list | **yes** | — | List of condition entries (min 1 required) |
| `defaults` | mapping | no | see below | Default analysis parameters |
| `plugins` | mapping | no | `{}` | Analysis plugin settings — **what** to compute |
| `plot_settings` | mapping | no | see below | Plot customization — **how** to visualize |

**Legacy key handling:**
- `analysis_settings:` is accepted as a backward-compatible alias for `plugins:` (emits deprecation warning).
- Unknown top-level keys are logged as warnings and ignored.

---

## `conditions[*]`

Each entry describes one simulation condition to include in the comparison.

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `label` | string | **yes** | — | Display name (must be unique across all conditions) |
| `config` | path | **yes** | — | Path to the simulation's `config.yaml`. Relative paths resolved from `comparison.yaml` location. |
| `replicates` | list of int | **yes** | — | Replicate numbers to include. A single `int` is auto-wrapped to a list. |

---

## `defaults`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `equilibration_time` | string | `"10ns"` | Time to discard as equilibration (e.g., `"10ns"`, `"5000ps"`) |
| `fdr_alpha` | float (0, 1] | `0.05` | Significance threshold for pairwise comparisons and ANOVA. Used as the Benjamini-Hochberg FDR threshold when `posthoc_method` is `"ttest_bh"`, and as the family-wise alpha threshold when `posthoc_method` is `"tukey_hsd"`. |
| `posthoc_method` | `"ttest_bh"` or `"tukey_hsd"` | `"ttest_bh"` | Post-hoc pairwise comparison method. See {doc}`posthoc_testing` for details. |
| `ttest_method` | `"student"` or `"welch"` | `"student"` | Two-sample t-test variance assumption. Only used when `posthoc_method` is `"ttest_bh"`. |

---

## `plugins`

Presence of a key **enables** that analysis. The value is a mapping of that
plugin's settings. An empty mapping (`rmsf: {}`) enables the plugin with all
defaults.

### `plugins.rmsf`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `selection` | string | `"protein and name CA"` | MDAnalysis selection string for RMSF computation |
| `reference_mode` | string | `"centroid"` | Reference structure: `"centroid"`, `"average"`, or `"frame"` |
| `reference_frame` | int | `null` | Required when `reference_mode` is `"frame"` |
| `reference_file` | path | `null` | Crystal/input PDB for secondary structure annotation bar on profiles |

### `plugins.secondary_structure`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `chain_id` | string | `"A"` | Chain letter for the protein to analyze via DSSP |

### `plugins.sasa`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | list | **(required)** | List of SASA run definitions (see sub-fields) |
| `probe_radius_nm` | float | `0.14` | SASA probe radius in nanometers |
| `n_sphere_points` | int | `960` | Number of sphere points for Shrake-Rupley SASA |
| `chunk_size` | int | `100` | Frames per chunk for memory management |

Each entry in `runs`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | string | **(required)** | Name for this SASA computation |
| `target_selection` | string | **(required)** | MDAnalysis selection for the target surface |
| `context_selection` | string | same as `target_selection` | Atoms to include in SASA context (affects shadowing) |
| `stride` | int | `1` | Frame stride |

### `plugins.catalytic_triad`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `name` | string | `"catalytic_triad"` | Display name for the triad analysis |
| `threshold` | float | `3.5` | Distance threshold in Angstroms (H-bond cutoff) |
| `pairs` | list | **(required)** | List of atom pair definitions |

Each entry in `pairs`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | string | **(required)** | Display label (e.g., `"Asp-His"`) |
| `selection_a` | string | **(required)** | MDAnalysis selection for atom/group A. Supports `midpoint(...)` syntax. |
| `selection_b` | string | **(required)** | MDAnalysis selection for atom/group B |

### `plugins.distances`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `threshold` | float | `3.5` | Global default threshold in Angstroms |
| `pairs` | list | **(required)** | List of distance pair definitions |

Each entry in `pairs`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | string | **(required)** | Display label (e.g., `"Ser77-Substrate"`) |
| `selection_a` | string | **(required)** | MDAnalysis selection for group A. Supports `com(...)` syntax. |
| `selection_b` | string | **(required)** | MDAnalysis selection for group B |
| `threshold` | float | global `threshold` | Per-pair threshold override |
| `below_label` | string | `"Below {threshold}Å"` | Display text for d ≤ threshold |
| `above_label` | string | `"Above {threshold}Å"` | Display text for d > threshold |

### `plugins.contacts`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `polymer_selection` | string | `"chainID C"` | MDAnalysis selection for polymer atoms |
| `protein_selection` | string | `"protein"` | MDAnalysis selection for protein atoms |
| `cutoff` | float | `4.5` | Contact distance cutoff in Angstroms |
| `grouping` | string | `"aa_class"` | Residue grouping: `"aa_class"`, `"secondary_structure"`, or `"none"` |
| `compute_residence_times` | bool | `true` | Whether to compute contact residence times |
| `compute_binding_preference` | bool | `false` | **Experimental.** Enable enrichment by residue group |
| `surface_exposure_threshold` | float | `0.2` | Relative SASA cutoff defining "surface exposed" (for binding preference) |
| `enzyme_pdb_for_sasa` | path | `null` | Path to enzyme PDB for standalone SASA computation (relative to `comparison.yaml`) |
| `include_default_aa_groups` | bool | `true` | Include built-in amino acid groups (aromatic, polar, nonpolar, charged) |
| `protein_groups` | mapping | `null` | Custom residue groups: `{group_name: [resid, ...]}` |
| `protein_partitions` | mapping | `null` | Mutually exclusive partitions for coverage plots: `{partition_name: [group_name, ...]}` |
| `fdr_alpha` | float | `0.05` | Per-plugin FDR threshold |
| `min_effect_size` | float | `0.5` | Minimum Cohen's d for practical significance |
| `top_residues` | int | `10` | Max residues shown per condition in formatted output |

### `plugins.rmsd`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | list | **(required)** | List of RMSD run definitions |

Each entry in `runs`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | string | **(required)** | Name for this RMSD computation (e.g., `"backbone"`) |
| `selection` | string | **(required)** | MDAnalysis selection |
| `align_selection` | string | same as `selection` | MDAnalysis selection for alignment |
| `stride` | int | `1` | Frame stride |

### `plugins.rg`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | list | **(required)** | List of Rg run definitions |

Each entry in `runs`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | string | **(required)** | Name for this Rg computation |
| `selection` | string | **(required)** | MDAnalysis selection |
| `stride` | int | `1` | Frame stride |

### `plugins.hydrogen_bonds`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `groups` | mapping | `{"protein": "chainid A", "polymer": "chainid C"}` | Named atom groups: `{name: "MDAnalysis selection"}` |
| `summaries` | list or mapping | one default summary (`protein_polymer` between `protein` and `polymer`) | Named H-bond summaries (see below) |
| `distance_cutoff` | float | `3.0` | H-bond distance cutoff in Angstroms |
| `angle_cutoff` | float | `150` | H-bond angle cutoff in degrees |
| `update_selections` | bool | `true` | Update atom selections every frame |
| `top_n_pairs` | int | `15` | Number of top residue pairs to report |
| `allow_empty_groups` | bool | `true` | Warn and skip summaries when a group selection matches no atoms (`false` = raise error) |
| `allow_overlapping_composition` | bool | `false` | Whether overlapping composition partitions are allowed |
| `composition` | mapping | `null` | Composition analysis settings |
| `timestep_ps` | float | `null` | Override trajectory timestep in picoseconds for time-axis plots |

Each summary entry in `summaries` has:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | yes | Unique summary name |
| `between` | `[group_a, group_b]` | exactly one of `between` / `within` | Inter-group H-bonds |
| `within` | `group_name` | exactly one of `between` / `within` | Intra-group H-bonds |

For mapping-form input, keys are treated as `name` values.

`composition` sub-fields:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `partitions` | mapping | — | Named partitions: `{name: "MDAnalysis selection"}` |

### `plugins.exposure`

```{admonition} Experimental
:class: warning

Exposure dynamics is an experimental analysis. Results should be interpreted
with caution and are subject to change.
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `exposure_threshold` | float | `0.20` | Fraction SASA defining "exposed" |
| `transient_lower` | float | `0.20` | Lower bound for transient classification |
| `transient_upper` | float | `0.80` | Upper bound for transient classification |
| `min_event_length` | int | `1` | Minimum consecutive frames for an event |
| `protein_chain` | string | `"A"` | Chain ID for protein |
| `protein_selection` | string | `"protein"` | MDAnalysis selection for protein |
| `polymer_selection` | string | `"chainID C"` | MDAnalysis selection for polymer |
| `polymer_resnames` | list of string | `null` | Residue names for enrichment analysis |
| `probe_radius_nm` | float | `0.14` | SASA probe radius (nm) |
| `n_sphere_points` | int | `960` | Number of sphere points |

### `plugins.binding_free_energy`

```{admonition} Experimental
:class: warning

Binding free energy decomposition is experimental and under active development.
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `units` | string | `"kT"` | Energy units: `"kT"`, `"kcal/mol"`, or `"kJ/mol"` |
| `compute_binding_preference` | bool | `true` | Recompute binding preference from contacts if no cache is available |
| `surface_exposure_threshold` | float | `0.2` | Minimum relative SASA for surface-exposed |
| `enzyme_pdb_for_sasa` | path | `null` | Enzyme PDB for SASA computation |
| `include_default_aa_groups` | bool | `true` | Include built-in amino acid class groups |
| `protein_groups` | mapping | `null` | Custom residue groups: `{name: [resid, ...]}` |
| `protein_partitions` | mapping | `null` | Mutually exclusive protein-group partitions |
| `polymer_type_selections` | mapping | `null` | Custom MDAnalysis selections per polymer type |
| `polymer_chain` | string | `"C"` | Chain ID used for polymer auto-detection |
| `fdr_alpha` | float | `0.05` | FDR threshold |

### `plugins.polymer_affinity`

```{admonition} Experimental
:class: warning

Polymer affinity scoring is experimental and under active development.
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `compute_binding_preference` | bool | `true` | Recompute binding preference from contacts if no cache is available |
| `surface_exposure_threshold` | float | `0.2` | Minimum relative SASA |
| `enzyme_pdb_for_sasa` | path | `null` | Enzyme PDB for SASA computation |
| `include_default_aa_groups` | bool | `true` | Use built-in AA groups |
| `protein_groups` | mapping | `null` | Custom residue groups |
| `protein_partitions` | mapping | `null` | Mutually exclusive partitions |
| `polymer_type_selections` | mapping | `null` | Custom MDAnalysis selections per polymer type |
| `polymer_chain` | string | `"C"` | Chain ID used for polymer auto-detection |
| `fdr_alpha` | float | `0.05` | FDR threshold |

### `plugins.polymer_bridging`

```{admonition} Experimental
:class: warning

Polymer bridging detection is experimental and under active development.
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `protein_selection` | string | `"protein"` | MDAnalysis selection for protein |
| `polymer_selection` | string | `"chainID C"` | MDAnalysis selection for polymer |
| `cutoff` | float | `4.5` | Contact distance cutoff in Angstroms for oligomer-protein contact detection |
| `min_ca_distance_angstrom` | float | `0.0` | Minimum frame-wise CA-CA distance to count as multisite (`0.0` disables geometric filtering) |

---

## `plot_settings`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `output_dir` | path | `"figures/"` | Directory for generated plots (relative to `comparison.yaml`) |
| `format` | string | `"png"` | Image format: `"png"`, `"pdf"`, or `"svg"` |
| `dpi` | int | `300` | Resolution for raster formats. Range: 50–600. |
| `style` | string | `"publication"` | Style preset: `"publication"`, `"presentation"`, or `"minimal"` |
| `color_palette` | string | `"tab10"` | Seaborn/matplotlib color palette name |
| `theme` | mapping | from style preset | Visual theme overrides (see below) |

### `plot_settings.theme`

All fields are optional — defaults are drawn from the selected `style` preset.

**Font sizes:**

| Field | publication | presentation | Description |
|-------|------------|-------------|-------------|
| `title_fontsize` | 13 | 18 | Axes title font size |
| `suptitle_fontsize` | 14 | 20 | Figure suptitle font size |
| `label_fontsize` | 11 | 15 | Axis label font size |
| `tick_fontsize` | 9 | 12 | Tick label font size |
| `legend_fontsize` | 9 | 12 | Legend entry font size |
| `annotation_fontsize` | 9 | 12 | Heatmap annotation font size |
| `small_fontsize` | 8 | 10 | Secondary annotation font size |
| `tiny_fontsize` | 7 | 9 | Fine-grained annotation font size |

**Bar chart:**

| Field | Default | Description |
|-------|---------|-------------|
| `bar_alpha` | `0.85` | Bar fill opacity |
| `bar_edgecolor` | `"black"` | Bar edge color |
| `bar_linewidth` | `0.5` | Bar edge line width |
| `bar_capsize` | `4` | Error bar cap size in points |

**Replicate dots:**

| Field | Default | Description |
|-------|---------|-------------|
| `dot_size` | `18` (minimal: `0`) | Scatter marker size |
| `dot_alpha` | `0.7` (minimal: `0`) | Dot opacity |
| `dot_color` | `"black"` | Dot color |

**Lines:**

| Field | Default | Description |
|-------|---------|-------------|
| `line_alpha` | `0.8` | Line plot opacity |
| `fill_alpha` | `0.25` (presentation: `0.3`, minimal: `0.15`) | fill_between band opacity |
| `reference_line_color` | `"black"` | Reference line color |
| `reference_line_style` | `"--"` | Reference line style |
| `reference_line_width` | `1.5` (presentation: `2.0`, minimal: `1.0`) | Reference line width |
| `highlight_line_alpha` | `0.5` | Vertical highlight line opacity |

**Axes chrome:**

| Field | Default | Description |
|-------|---------|-------------|
| `hide_top_spine` | `true` | Hide top axis spine |
| `hide_right_spine` | `true` | Hide right axis spine |

**Title & legend:**

| Field | Default | Description |
|-------|---------|-------------|
| `title_fontweight` | `"bold"` | Title font weight |
| `legend_loc` | `"center left"` | Matplotlib legend location |
| `legend_bbox` | `[1.02, 0.5]` | bbox_to_anchor for legend placement |
| `show_watermark` | `true` | Render "Made by PolyzyMD" watermark |

### Per-Analysis Plot Settings

Per-analysis plot customization keys go under `plot_settings:` at the same
level as `style`, `dpi`, etc.

**`plot_settings.rmsf`:**

| Field | Default | Description |
|-------|---------|-------------|
| `show_error` | `true` | Show SEM fill_between bands |
| `highlight_residues` | `[]` | Residue IDs for vertical reference lines |
| `figsize_profile` | `[14, 4]` | Per-residue profile figure size |
| `figsize_comparison` | `[8, 6]` | Bar comparison figure size |

**`plot_settings.catalytic_triad`:**

| Field | Default | Description |
|-------|---------|-------------|
| `generate_kde_panel` | `true` | Multi-row KDE panel |
| `generate_bars` | `true` | Threshold bar chart |
| `generate_2d_kde` | `false` | 2D joint KDE |
| `kde_xlim` | `[0, 7]` | X-axis range for KDE (Angstroms) |

**`plot_settings.distances`:**

| Field | Default | Description |
|-------|---------|-------------|
| `show_threshold` | `true` | Threshold line on distributions |
| `use_kde` | `true` | KDE vs histogram |
| `generate_state_bars` | `true` | Above/below threshold bars |

**`plot_settings.contacts`:**

| Field | Default | Description |
|-------|---------|-------------|
| `generate_enrichment_heatmap` | `true` | Binding preference heatmap |
| `generate_enrichment_bars` | `true` | Enrichment bar chart |
| `generate_system_coverage_heatmap` | `true` | System coverage heatmap |
| `generate_system_coverage_bars` | `true` | System coverage bar chart |
| `generate_contact_fraction_profile` | `true` | Per-residue contact fraction profile |
| `generate_residence_time_profile` | `true` | Per-residue residence time profile |

**`plot_settings.binding_free_energy`:**

| Field | Default | Description |
|-------|---------|-------------|
| `generate_heatmap` | `true` | ΔG_sel heatmap |
| `generate_bars` | `true` | ΔG_sel bar chart |
| `colormap` | `"RdBu_r"` | Diverging colormap for heatmap |

**`plot_settings.polymer_affinity`:**

| Field | Default | Description |
|-------|---------|-------------|
| `generate_stacked_bars` | `true` | Total score by condition |
| `generate_group_bars` | `true` | Per-group contributions |

**`plot_settings.secondary_structure`:**

| Field | Default | Description |
|-------|---------|-------------|
| `generate_timeline` | `true` | Residue × time SS heatmap |
| `generate_content_bars` | `true` | Helix/strand/coil fraction bars |
| `generate_individual_bars` | `true` | One bar chart per SS type |
| `generate_diff_heatmap` | `true` | Δ(helix persistence) vs control |
| `diff_colormap` | `"RdBu_r"` | Diverging colormap for diff heatmap |

---

```{tip}
**Common tips:**

- Run `polyzymd compare validate` to check your `comparison.yaml` for errors
  before launching a full analysis run.
- Relative paths in `config:` are resolved from the directory containing
  `comparison.yaml`, not from your working directory.
- An empty plugin mapping (e.g., `rmsf: {}`) enables the analysis with all
  default settings — you only need to specify fields you want to override.
- Set `control:` to match one of your condition labels to get Δ-from-control
  columns in comparison tables and plots.
```
