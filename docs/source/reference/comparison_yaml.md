# `comparison.yaml` Schema Reference

The `comparison.yaml` file defines a cross-condition analysis project. It
specifies which simulation conditions to compare, which analysis plugins to
run, and how to visualize results. Create one with `polyzymd compare init -n
<name>` and place it at the root of your comparison project directory.

Source of truth: {func}`polyzymd.config.comparison.ComparisonConfig` in
`src/polyzymd/config/comparison.py`.

```{important}
Plugin settings path fields are resolved relative to the directory containing
`comparison.yaml`.

For example, in:

`plugins.rmsf.reference_file`, condition `config` paths, and other
plugin-declared path fields, a relative path like `structures/enzyme.pdb` is
interpreted as:

`<comparison_yaml_parent>/structures/enzyme.pdb`
```

For CLI commands that consume this file, see
{doc}`analysis_comparison_reference`. For directory layout and data
expectations, see {doc}`data_requirements`.

Typical local workflow:

```bash
pixi run -e build polyzymd compare validate -f comparison.yaml
pixi run -e build polyzymd compare run rmsf -f comparison.yaml
pixi run -e build polyzymd compare plot-all -f comparison.yaml
```

Typical SLURM workflow:

```bash
pixi run -e build polyzymd compare submit sasa -f comparison.yaml --dry-run
pixi run -e build polyzymd compare submit sasa -f comparison.yaml --partition <part>
pixi run -e build polyzymd compare status sasa -f comparison.yaml
pixi run -e build polyzymd compare finalize sasa -f comparison.yaml
pixi run -e build polyzymd compare plot-all -f comparison.yaml
```

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
| `mda_backend_policy` | mapping | no | `{}` | Optional MDAnalysis internal backend policy for job-backed analyses |
| `plot_settings` | mapping | no | see below | Plot customization — **how** to visualize |

Unknown top-level keys raise a `ValueError` listing the invalid keys and valid
alternatives. Use `plugins:` for analysis plugin settings; unsupported keys such
as `analysis_settings:` are rejected.

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

`equilibration_time` is interpreted as an absolute MDAnalysis trajectory
timestamp when the loaded trajectory exposes finite frame times. This handles
continuation runs where the first loaded segment may begin after 0 ps. If frame
timestamps are unavailable, PolyzyMD treats the first loaded frame as time zero.

## `mda_backend_policy`

The default policy is empty and forwards no backend-related keyword arguments to
MDAnalysis. This avoids nested oversubscription: PolyzyMD schedules work across
conditions/replicates, while each replicate remains serial unless you explicitly
opt into an MDAnalysis backend.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `backend` | string | `null` | Backend name forwarded to `AnalysisBase.run()`, such as `"multiprocessing"` or `"dask"` |
| `n_workers` | positive int | `null` | Worker count forwarded only when `backend` is set |
| `n_parts` | positive int | `null` | Optional partition count forwarded only when `backend` is set |

Example opt-in for local MDAnalysis internal parallelism:

```yaml
mda_backend_policy:
  backend: "multiprocessing"
  n_workers: 2
  n_parts: 2
```

Function-adapter jobs generated by the simple scaffold reject non-default
backend policies; use an `AnalysisBase`-compatible job for MDAnalysis internal
parallelism.

---

## `plugins`

Presence of a key **enables** that analysis. The value is a mapping of that
plugin's settings. An empty mapping (`rmsf: {}`) enables the plugin with all
defaults.

### `plugins.rmsf`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `selection` | string | `"protein and name CA"` | MDAnalysis selection string for RMSF computation |
| `reference_mode` | string | `"centroid"` | Reference structure: `"centroid"`, `"average"`, `"frame"`, or `"external"` |
| `reference_frame` | int | `null` | Required when `reference_mode` is `"frame"` |
| `reference_file` | path | `null` | Path to external PDB reference structure. Required when `reference_mode` is `"external"`. Also used for secondary structure annotation on profile plots. |
| `alignment_selection` | string | `"protein and name CA"` | MDAnalysis selection used for trajectory alignment before RMSF calculation |
| `centroid_selection` | string | `"protein"` | MDAnalysis selection used to compute the centroid reference structure when `reference_mode` is `"centroid"` |

### `plugins.secondary_structure`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `chain_id` | string | `"A"` | Chain letter for the protein to analyze via DSSP |
| `selection` | string | `null` | Explicit MDAnalysis protein-residue selection. Overrides `chain_id` when provided. |

The default secondary-structure selection is `protein and chainid A`, preserving
the PolyzyMD/PDB convention that chain A is the protein. GROMACS `.gro`
topologies may not preserve chain IDs, so use `selection` for those files, for
example `selection: "protein"`, `selection: "protein and resid 1:269"`, or
`selection: "protein and resindex 0:268"`. DSSP needs complete residues; do not
use CA-only selections such as `protein and name CA`.

### `plugins.sasa`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | list | **(required)** | List of SASA run definitions (see sub-fields) |
| `probe_radius_nm` | float | `0.14` | MDTraj Shrake-Rupley probe radius in nanometers |
| `n_sphere_points` | int | `960` | Number of sphere points for MDTraj Shrake-Rupley SASA |
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
| `description` | string | `null` | Optional description of the triad (e.g., `"Ser-His-Asp catalytic triad"`) |
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
| `use_pbc` | bool | `true` | Apply periodic boundary conditions to distance calculations |
| `align_trajectory` | bool | `true` | Align trajectory before computing distances |
| `alignment_selection` | string | `"protein and name CA"` | MDAnalysis selection used for trajectory alignment |
| `alignment_mode` | string | `"centroid"` | Alignment reference mode: `"centroid"` or `"frame"` |
| `alignment_frame` | int | `null` | Frame index to use as reference when `alignment_mode` is `"frame"` |

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
| `polymer_selection` | string | `"chainid C"` | MDAnalysis selection for polymer atoms |
| `protein_selection` | string | `"chainid A"` | MDAnalysis selection for protein atoms |
| `cutoff` | float | `4.5` | Contact distance cutoff in Angstroms |
| `grouping` | string | `"aa_class"` | Residue grouping: `"aa_class"`, `"secondary_structure"`, or `"none"` |
| `compute_residence_times` | bool | `true` | Whether to compute aggregate residence-time summaries and plots. When `false`, per-replicate contact events are still stored and the canonical artifact identity changes. |
| `protein_groups` | mapping | `null` | Custom residue groups: `{group_name: [resid, ...]}` |
| `protein_partitions` | mapping | `null` | Mutually exclusive partitions for contact-fraction and residence-time plots: `{partition_name: [group_name, ...]}` |
| `polymer_types` | list of string | `null` | Explicit polymer type labels. If `null`, types are auto-detected from topology. |
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
| `selection` | string | **(required)** | MDAnalysis selection for RMSD atoms |
| `alignment_selection` | string | same as `selection` | MDAnalysis selection for alignment |
| `reference_mode` | string | `"centroid"` | Reference structure mode: `"centroid"` or `"frame"` |
| `reference_frame` | int | `0` | Frame index to use as reference when `reference_mode` is `"frame"` |
| `reference_file` | path | `null` | Path to external PDB reference structure |
| `centroid_selection` | string | `null` | MDAnalysis selection for centroid computation. If `null`, uses `alignment_selection`. |
| `convergence_window_size_ns` | float | `15.0` | Rolling window size in nanoseconds for convergence detection |
| `convergence_step_size_ns` | float | `5.0` | Step size in nanoseconds between convergence windows |
| `convergence_slope_threshold` | float | `0.0005` | Maximum slope (Å/ns) for a window to be considered converged |
| `convergence_sustained_for_ns` | float | `15.0` | Duration in nanoseconds that convergence must be sustained |

### `plugins.rg`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | list | **(required)** | List of Rg run definitions |

Each entry in `runs`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | string | **(required)** | Name for this Rg computation |
| `selection` | string | **(required)** | MDAnalysis selection for Rg atoms |
| `calculation_mode` | string | `"selection"` | Computation mode: `"selection"` (single Rg for the whole selection) or `"fragments"` (per-fragment Rg) |
| `fragment_weighting` | string | `"equal"` | How to weight fragments when `calculation_mode` is `"fragments"`: `"equal"` or `"mass"` |
| `save_fragment_distribution` | bool | `true` | Save per-frame fragment Rg distributions |
| `histogram_bins` | int | `50` | Number of bins for Rg distribution histograms |

### `plugins.hydrogen_bonds`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `groups` | mapping | `{"protein": "chainid A", "polymer": "chainid C"}` | Named atom groups: `{name: "MDAnalysis selection"}` |
| `summaries` | list or mapping | one default summary (`protein_polymer` between `protein` and `polymer`) | Named H-bond summaries (see below) |
| `distance_cutoff` | float | `3.0` | H-bond distance cutoff in Angstroms |
| `angle_cutoff` | float | `150` | H-bond angle cutoff in degrees |
| `update_selections` | bool | `true` | Update atom selections every frame |
| `top_n_pairs` | int | `15` | Number of top residue pairs to report |
| `allow_empty_groups` | bool | `true` | Allow empty group selections: `true` = warn and skip summaries when a group matches no atoms; `false` = raise error |
| `allow_overlapping_composition` | bool | `false` | Whether overlapping composition partitions are allowed |
| `composition` | mapping | `null` | Composition analysis settings |
| `timestep_ps` | float | `null` | Override trajectory timestep in picoseconds for time-axis plots |

Time-axis plots assume uniformly saved frames. PolyzyMD converts frame index to
time as `frame_index * timestep_ps`; variable-timestep concatenated
trajectories are not supported.

Each summary entry in `summaries` has:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | yes | Unique summary name |
| `between` | `[group_a, group_b]` | exactly one of `between` / `within` | Inter-group H-bonds |
| `within` | `group_name` | exactly one of `between` / `within` | Intra-group H-bonds |

For mapping-form input, keys are treated as `name` values.

Hydrogen detection uses MDAnalysis `HydrogenBondAnalysis` with hydrogens
selected as `(<group union>) and element H`; topologies need explicit hydrogens
and usable element metadata.

`composition` sub-fields:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `partitions` | mapping | — | Named partitions: `{name: "MDAnalysis selection"}` |

---

## `plot_settings`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `output_dir` | path | `"figures/"` | Directory for generated plots (relative to `comparison.yaml`) |
| `format` | string | `"png"` | Image format: `"png"`, `"pdf"`, or `"svg"` |
| `dpi` | int | `300` | Resolution for raster formats. Range: 50–600. |
| `style` | string | `"compact"` | PolyzyMD theme preset: `"compact"`, `"large_elements"`, or `"low_ink"` |
| `color_palette` | string | `"tab10"` | Seaborn/matplotlib color palette name |
| `semantic_colors` | mapping | disabled | Optional condition-label color and display-order rules for condition-series plots |
| `theme` | mapping | from style preset | Visual theme overrides (see below) |

`style` selects a PolyzyMD built-in theme preset for standard analysis plots. It
is not a matplotlib or seaborn stylesheet, and it does not control `format`,
`dpi`, per-analysis figure sizes, or color palettes.

`theme` values are merged on top of the selected preset, so you can choose a
base style and override only the fields that need project-specific changes.

### `plot_settings.semantic_colors`

Semantic colors let a comparison project encode condition meaning directly in
figures. The settings are optional and disabled by default; when disabled,
plots keep using `color_palette` and each plotter's existing category colors.

Semantic ordering is **plot-only**. It changes the display order of conditions
in figures, but it does not mutate comparison statistics, rankings, cached
artifacts, or JSON result files.

Top-level fields:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `false` | Opt in to semantic condition colors and plot ordering |
| `order` | list of string | `[]` | Explicit plot display order by condition label. Labels not present keep their relative order after condition-level `order` sorting. |
| `manual_colors` | mapping | `{}` | Direct color overrides by exact condition label. Highest precedence color rule. |
| `conditions` | mapping | `{}` | Per-condition semantic metadata keyed by exact condition label |
| `families` | mapping | `{}` | Family-level colormap rules keyed by family name |
| `control_color` | color | `"black"` | Color used for the configured `control` condition or a condition with `role: control` |
| `missing_color` | color | `"lightgray"` | Fallback color for conditions with incomplete semantic metadata |
| `default_color` | color or `null` | `null` | Fallback for labels missing from `conditions`. If `null`, the regular palette color is used. |

`conditions.<label>` fields:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `color` | color or `null` | `null` | Direct color for this condition, after `manual_colors` and before control/family rules |
| `family` | string or `null` | `null` | Semantic family name used to look up `families.<family>` |
| `value` | scalar or `null` | `null` | Numeric or ordinal value mapped through the family color rule |
| `order` | int or `null` | `null` | Plot-only display order used after explicit `semantic_colors.order` |
| `role` | string or `null` | `null` | Optional semantic role. Use `control` to apply `control_color`. |

`families.<family>` fields:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `colormap` | string | `"viridis"` | Matplotlib colormap name for values in this family |
| `scale` | `"linear"` or `"ordinal"` | `"linear"` | Map numeric values continuously (`linear`) or ordered categories/steps discretely (`ordinal`) |
| `value_order` | list | `[]` | Explicit value order for `ordinal` mapping. If omitted, observed values are used in label order. |
| `vmin` | float or `null` | `null` | Lower bound for `linear` normalization. If omitted, observed values set the bound. |
| `vmax` | float or `null` | `null` | Upper bound for `linear` normalization. If omitted, observed values set the bound. |
| `colormap_range` | two floats | `[0.0, 1.0]` | Fractional colormap interval to sample, useful for avoiding colors that are too pale or too dark |
| `reverse` | bool | `false` | Reverse the value-to-colormap direction |
| `value_colors` | mapping | `{}` | Explicit color overrides by value. These override the family colormap for matching values. |

Color precedence for each condition label is:

1. `semantic_colors.manual_colors.<label>`
2. `semantic_colors.conditions.<label>.color`
3. `semantic_colors.control_color` when the label is the top-level `control` or
   the condition has `role: control`
4. `families.<family>.value_colors.<value>`
5. `families.<family>` colormap mapping
6. `semantic_colors.missing_color` for incomplete condition metadata
7. `semantic_colors.default_color` or the regular `color_palette` for labels
   missing from `semantic_colors.conditions`

Semantic colors apply to plots where colors represent comparison conditions.
Non-condition categories, such as secondary-structure states or residue classes,
may still use categorical palettes or plot-specific colormaps.

### `plot_settings.theme`

All fields are optional. Defaults are drawn from the selected `style` preset,
then any values under `theme:` override individual fields.

#### Theme presets

| Preset | Use when | Notes |
|--------|----------|-------|
| `compact` | You want the default compact print-style output. | Uses moderate fonts, replicate dots, bar edges, and reference lines. |
| `large_elements` | You need slides, posters, or high-visibility figures. | Increases font sizes, replicate dot size, bar line width, error-bar caps, reference-line width, and fill opacity. |
| `low_ink` | You want simpler, lower-ink plots. | Hides replicate dots, removes bar edges, and reduces reference-line width and fill opacity. |

#### Tweakable `PlotTheme` fields

| Field | `compact` | `large_elements` | `low_ink` | Description |
|-------|-----------|------------------|-----------|-------------|
| `title_fontsize` | `13` | `18` | `13` | Axes title font size |
| `suptitle_fontsize` | `14` | `20` | `14` | Figure suptitle font size |
| `label_fontsize` | `11` | `15` | `11` | Axis label font size |
| `tick_fontsize` | `9` | `12` | `9` | Tick label font size |
| `legend_fontsize` | `9` | `12` | `9` | Legend entry font size |
| `annotation_fontsize` | `9` | `12` | `9` | Heatmap annotation font size |
| `small_fontsize` | `8` | `10` | `8` | Secondary annotation font size |
| `tiny_fontsize` | `7` | `9` | `7` | Fine-grained annotation font size |
| `bar_alpha` | `0.85` | `0.85` | `0.85` | Bar fill opacity |
| `bar_edgecolor` | `"black"` | `"black"` | `"none"` | Bar edge color |
| `bar_linewidth` | `0.5` | `0.8` | `0.0` | Bar edge line width |
| `bar_capsize` | `4` | `5` | `3` | Error bar cap size in points |
| `dot_size` | `18` | `30` | `0` | Scatter marker size for replicate dots |
| `dot_alpha` | `0.7` | `0.7` | `0.0` | Replicate dot opacity |
| `dot_color` | `"black"` | `"black"` | `"black"` | Replicate dot color |
| `line_alpha` | `0.8` | `0.8` | `0.8` | Line plot opacity |
| `fill_alpha` | `0.25` | `0.3` | `0.15` | `fill_between` band opacity |
| `reference_line_color` | `"black"` | `"black"` | `"black"` | Reference line color |
| `reference_line_style` | `"--"` | `"--"` | `"--"` | Reference line style |
| `reference_line_width` | `1.5` | `2.0` | `1.0` | Reference line width |
| `highlight_line_alpha` | `0.5` | `0.5` | `0.5` | Vertical highlight line opacity |
| `hide_top_spine` | `true` | `true` | `true` | Hide top axis spine |
| `hide_right_spine` | `true` | `true` | `true` | Hide right axis spine |
| `title_fontweight` | `"bold"` | `"bold"` | `"bold"` | Title font weight |
| `legend_loc` | `"center left"` | `"center left"` | `"center left"` | Matplotlib legend location |
| `legend_bbox` | `[1.02, 0.5]` | `[1.02, 0.5]` | `[1.02, 0.5]` | `bbox_to_anchor` for legend placement |
| `show_watermark` | `true` | `true` | `true` | Render the "Made by PolyzyMD" watermark |

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
| `generate_contact_fraction_profile` | `true` | Per-residue contact fraction profile |
| `generate_residence_time_profile` | `true` | Per-residue residence time profile |
| `generate_cf_by_aa_class_bars` | `true` | Contact fraction by amino acid class bar chart |
| `generate_cf_by_partition_bars` | `true` | Contact fraction by user partition bar charts |
| `generate_rt_by_aa_class_bars` | `true` | Residence time by amino acid class bar chart |
| `generate_rt_by_partition_bars` | `true` | Residence time by user partition bar charts |

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
