# `comparison.yaml` Schema Reference

The `comparison.yaml` file lists the simulation conditions of a comparison,
with their replicates, the control and the default equilibration window.
Create one with `polyzymd compare init -n <name>`. No analysis is configured
in it: every shipped analysis runs with `polyzymd analyze NAME -c
<config.yaml> ...` on the simulation configs, and `polyzymd analyze` does not
read `comparison.yaml`. The `plugins` and `plot_settings` sections are
retired.

`comparison.yaml` is still read by:

- `polyzymd compare validate`, which checks the conditions and their config
  paths;
- `polyzymd compare run NAME` for a shipped analysis name and `polyzymd
  analyze NAME -f comparison.yaml`, which run nothing and print the
  equivalent `polyzymd analyze NAME -c <config> --label <label> ...
  --replicates ... --eq ...` command built from the file's conditions;
- the plugin framework commands (`polyzymd compare run`, `run-all`,
  `plot-all`, `submit`, `status`, `finalize`) for analysis plugins you
  register yourself. No shipped analysis is a plugin any more, and the
  framework is being removed.

Source of truth: {func}`polyzymd.config.comparison.ComparisonConfig` in
`src/polyzymd/config/comparison.py`.

```{important}
Condition `config` paths are resolved relative to the directory containing
`comparison.yaml`. For example, `../sbma_100/config.yaml` is interpreted as
`<comparison_yaml_parent>/../sbma_100/config.yaml`.
```

For the commands that run the analyses, see
{doc}`../how_to/analysis_agent_protocol` and {doc}`cli_reference`. For
directory layout and data expectations, see {doc}`data_requirements`.

Typical workflow:

```bash
pixi run -e analysis polyzymd compare validate -f comparison.yaml
pixi run -e analysis polyzymd analyze hydrogen_bonds \
  -c ../no_polymer/config.yaml -c ../sbma_100/config.yaml \
  --label "No Polymer" --label "100% SBMA" --eq 10ns
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
| `plugins` | mapping | no | `{}` | Retired; see {ref}`comparison-yaml-retired`. Settings of analysis plugins you register yourself |
| `mda_backend_policy` | mapping | no | `{}` | Optional MDAnalysis internal backend policy for the jobs of registered plugins |
| `plot_settings` | mapping | no | see below | Retired; see {ref}`comparison-yaml-retired`. Plot customization of registered plugins |

Unknown top-level keys raise a `ValueError` listing the invalid keys and valid
alternatives; unsupported keys such as `analysis_settings:` are rejected.

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

(comparison-yaml-retired)=
## Retired `plugins` and `plot_settings` blocks

A `plugins.<name>` or `plot_settings.<name>` block for a shipped analysis still
loads, so an old `comparison.yaml` keeps working for `polyzymd compare
validate`, and each block is ignored with one `UserWarning`. The retired names
and the function each warning names are:

| Block | `polyzymd analyze` command | Python |
|---|---|---|
| `rg` | `polyzymd analyze rg` | `study.timeseries` with `functions.radius_of_gyration` |
| `rmsd` | `polyzymd analyze rmsd` | `study.timeseries` with `functions.rmsd` |
| `rmsf` | `polyzymd analyze rmsf` | `study.per_replicate` with `functions.rmsf` |
| `sasa` | `polyzymd analyze sasa` | `study.timeseries` with `functions.sasa` |
| `secondary_structure` | `polyzymd analyze secondary_structure` | `study.per_replicate` with `functions.dssp_occupancy` |
| `distances` | `polyzymd analyze distances --set pairs=<pairs.yaml>` | `study.timeseries` with `functions.pair_distance` |
| `contacts` | `polyzymd analyze contacts` | `study.per_replicate` with `functions.residue_occlusion`, `residue_contacts` for `--set method=distance` and `contact_lifetimes` |
| `hydrogen_bonds` | `polyzymd analyze hydrogen_bonds` | `study.per_replicate` with `functions.hydrogen_bonds`, `hbond_lifetimes`, `residue_hbond_occupancy` and `residue_pair_hbond_occupancy` |
| `catalytic_triad` | `polyzymd analyze distances --set pairs=<pairs.yaml>` for the distances | the triad routine, which counts each triad hydrogen bond with `functions.hbond_count` and combines them with `Timeseries.transform` |

Every warning names the replacement command and the Python function, and
ends with the address of {doc}`../how_to/analysis_agent_protocol`
(<https://polyzymd.readthedocs.io/en/latest/how_to/analysis_agent_protocol.html>)
and the agent skill `.claude/skills/polyzymd-analyze/SKILL.md`, which an agent
can be pointed at to learn the protocol. The `catalytic_triad` warning points
to the triad routine instead, {doc}`../how_to/analysis_triad_quickstart`. The
settings of a retired block are not carried over: pass them to `polyzymd
analyze` with `--set`, as the how-to page of each analysis shows.

Any other key in `plugins` must name an analysis plugin you registered
yourself; an unknown key raises a `ValueError` that names the registered
plugins and the same address and skill.

## `PlotSettings` fields

`polyzymd analyze` does not read `plot_settings`: its figures use the default
`PlotSettings`. In Python, every figure of the study API takes a
`plot_settings=` argument, a {class}`polyzymd.config.comparison.PlotSettings`
with the fields below, for example
`values.plot(plot_settings=PlotSettings(format="pdf", style="large_elements"))`.
The plugin framework commands also read these fields from a `plot_settings:`
block for registered plugins.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `output_dir` | path | `"figures/"` | Directory for generated plots of registered plugins (relative to `comparison.yaml`); the study figures take the folder as their own `output_dir` argument |
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

### `semantic_colors`

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

### `theme`

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

### Per-analysis plot settings

No shipped analysis has per-analysis plot settings. A `plot_settings.<name>`
block for a retired name is ignored with the warning of
{ref}`comparison-yaml-retired`. Mark residues on the RMSF profiles with
`polyzymd analyze rmsf --set highlight_residues='[...]'`.

---

```{tip}
**Common tips:**

- Run `polyzymd compare validate` to check the conditions and their config
  paths.
- Relative paths in `config:` are resolved from the directory containing
  `comparison.yaml`, not from your working directory.
- In `polyzymd analyze`, the first `-c` is the control.
```
