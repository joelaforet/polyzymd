# Customizing Plots for Publication

PolyzyMD generates plots automatically when you run analyses. This guide shows
you how to control output format, resolution, figure size, and style using
`plot_settings` in `comparison.yaml`.

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
    generate_enrichment_heatmap: true
    generate_system_coverage_heatmap: false
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
