# Customizing Plots for Publication

PolyzyMD generates plots automatically when you run analyses. This guide shows
you how to control output format, resolution, figure size, and style using
`plot_settings` in `comparison.yaml`.

:::{admonition} Environment Setup
:class: tip

All commands below assume you have activated the PolyzyMD pixi environment:

```bash
pixi shell -e build
```

Alternatively, prefix each command with `pixi run -e build`.
:::

## The `plot_settings` block

Add this block to `comparison.yaml`:

```yaml
plot_settings:
  format: "png"              # or "pdf", "svg"
  dpi: 300
  style: "publication"       # or "presentation", "minimal"
```

These fields are defined in PolyzyMD's `PlotSettings` model:

- `format`
  - Output file format
  - Allowed values: `png`, `pdf`, `svg`
- `dpi`
  - Plot resolution
  - Allowed range: `50` to `600`
  - Primarily affects raster output (`png`)
- `style`
  - Global style preset
  - Allowed values: `publication`, `presentation`, `minimal`

You can also set these optional global fields:

- `output_dir` (default: `figures/`)
- `color_palette` (default: `tab10`)
- `theme` (fine-grained visual overrides)

Example with all global fields:

```yaml
plot_settings:
  output_dir: "figures/"
  format: "png"
  dpi: 300
  style: "publication"
  color_palette: "tab10"
```

## Use semantic colors for condition series

Use `plot_settings.semantic_colors` when the condition colors should carry
experimental meaning, such as polymer family, composition, or dose. Semantic
colors are opt-in and apply to condition-series plots. Plot elements that are
not conditions, such as residue classes or secondary-structure categories, may
keep their categorical colors.

```yaml
conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]
  - label: "100% SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3]
  - label: "100% EGMA"
    config: "../egma_100/config.yaml"
    replicates: [1, 2, 3]
  - label: "1% EGPMA"
    config: "../egpma_1/config.yaml"
    replicates: [1, 2, 3]
  - label: "2% EGPMA"
    config: "../egpma_2/config.yaml"
    replicates: [1, 2, 3]
  - label: "5% EGPMA"
    config: "../egpma_5/config.yaml"
    replicates: [1, 2, 3]
  - label: "10% EGPMA"
    config: "../egpma_10/config.yaml"
    replicates: [1, 2, 3]

control: "No Polymer"

plot_settings:
  style: "publication"
  semantic_colors:
    enabled: true

    # Plot-only order. This does not change statistics, rankings, or artifacts.
    order:
      - "No Polymer"
      - "100% SBMA"
      - "100% EGMA"
      - "1% EGPMA"
      - "2% EGPMA"
      - "5% EGPMA"
      - "10% EGPMA"

    control_color: "#222222"
    missing_color: "#bdbdbd"

    conditions:
      "No Polymer":
        role: control
        order: 0
      "100% SBMA":
        family: sbma
        value: 100
        order: 10
      "100% EGMA":
        family: egma
        value: 100
        order: 20
      "1% EGPMA":
        family: egpma
        value: 1
        order: 30
      "2% EGPMA":
        family: egpma
        value: 2
        order: 40
      "5% EGPMA":
        family: egpma
        value: 5
        order: 50
      "10% EGPMA":
        family: egpma
        value: 10
        order: 60

    families:
      sbma:
        colormap: "Blues"
        scale: linear
        vmin: 0
        vmax: 100
        colormap_range: [0.55, 0.85]
      egma:
        colormap: "Greens"
        scale: linear
        vmin: 0
        vmax: 100
        colormap_range: [0.55, 0.85]
      egpma:
        colormap: "Purples"
        scale: ordinal
        value_order: [1, 2, 5, 10]
        colormap_range: [0.35, 0.9]
```

Labels containing `%` should be quoted in YAML. Every key under
`semantic_colors.conditions`, every entry in `semantic_colors.order`, and the
top-level `control` value must match a `conditions[*].label` exactly.

For small percentage series such as `1`, `2`, `5`, and `10`, `scale: ordinal`
often gives more readable figures than a linear scale because each configured
level receives a visually distinct shade. For numeric gradients with many
levels, `scale: linear` is usually appropriate.

For multi-component chemistry, avoid naive additive RGB mixing: colors that add
component hues together can become hard to interpret and hard to reproduce.
Prefer a simpler mapping where hue identifies the family and shade identifies
the ordered composition or dose within that family.

After changing only plot settings, regenerate figures from cached comparison
results:

```bash
polyzymd compare plot-all -f comparison.yaml
```

## Changing output format

Switch file format by changing `plot_settings.format`.

### PNG

```yaml
plot_settings:
  format: "png"
  dpi: 300
  style: "publication"
```

Use PNG for quick drafts, sharing in chat/slides, and embedding in docs.

### PDF

```yaml
plot_settings:
  format: "pdf"
  dpi: 300
  style: "publication"
```

PDF is vector output, so it scales cleanly.

### SVG

```yaml
plot_settings:
  format: "svg"
  dpi: 300
  style: "publication"
```

SVG is vector output and is convenient for web use and post-editing.

## Changing DPI

Set `plot_settings.dpi` to control raster resolution.

```yaml
plot_settings:
  format: "png"
  dpi: 150
  style: "publication"
```

Common ranges:

- `72-150`: screen and draft output
- `300`: print/publication output
- `600`: high-resolution print output

## Per-plugin plot settings

In addition to global settings, plugins can define their own plot options under
`plot_settings.<plugin_name>`.

### RMSD example

```yaml
plot_settings:
  format: "pdf"
  style: "publication"
  rmsd:
    show_per_replicate: true
    figsize: [10, 6]
    timeseries_figsize: [12, 5]
    show_convergence_plots: true
    convergence_figsize: [12, 5]
```

What changes:

- `show_per_replicate: true` overlays each replicate trace
- `timeseries_figsize` changes the width/height of RMSD time plots
- `show_convergence_plots: true` adds convergence diagnostics

### RMSF example

```yaml
plot_settings:
  style: "publication"
  rmsf:
    show_error: true
    highlight_residues: [77, 133, 156]
    figsize_profile: [14, 4]
    figsize_comparison: [8, 6]
```

What changes:

- `show_error` turns error shading/bars on or off
- `highlight_residues` adds vertical markers at selected residue IDs
- `figsize_*` controls profile and comparison figure sizes

### Contacts example

```yaml
plot_settings:
  style: "publication"
  contacts:
    generate_contact_fraction_profile: true
    generate_residence_time_profile: true
    generate_cf_by_aa_class_bars: true
    generate_cf_by_partition_bars: true
    figsize_contact_fraction_profile: [16, 5]
    show_contact_fraction_profile_error: true
    highlight_residues: [77, 133, 156]
```

What changes:

- `generate_*` flags turn specific plot families on/off
- `figsize_contact_fraction_profile` sets profile dimensions
- `show_contact_fraction_profile_error` toggles profile error bands

## Re-generating plots after changing settings

After editing `comparison.yaml`, regenerate figures with:

```bash
polyzymd compare plot-all -f comparison.yaml
```

This command re-draws figures from existing comparison results. It does not
recompute per-replicate or aggregated analysis data.

If a plugin has no cached comparison result yet, run that comparison first
(for example, `polyzymd compare run rmsf -f comparison.yaml`).

## See Also

- [How to Compare Simulation Conditions](analysis_compare_conditions.md)
- [Comparison and Plotting Reference](../reference/analysis_comparison_reference.md)
- [RMSD Analysis Reference](../reference/analysis_rmsd_reference.md)
- [RMSF Analysis Reference](../reference/analysis_rmsf_reference.md)
