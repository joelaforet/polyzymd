# Customizing Plots for Publication

`polyzymd analyze` draws each analysis's figures with the default plot
settings, 300 dpi PNG files in the `compact` style. This guide shows you how
to redraw them with another output format, resolution, PolyzyMD theme preset
or condition colors, by passing a `PlotSettings` to the `plot` methods of the
study API.

:::{admonition} Environment Setup
:class: tip

All plotting commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## Write the settings and redraw a figure

Write the settings to a YAML file, for example `plot_settings.yaml`:

```yaml
format: "pdf"              # or "png", "svg"
dpi: 300
style: "compact"           # or "large_elements", "low_ink"
```

Load it into a `PlotSettings` and pass it as `plot_settings` when you draw a
result of the study API. Run in the directory where `polyzymd analyze
hydrogen_bonds -c ... -c ...` ran, and with the condition labels given there
with `--label`, the measurement below reads back the per-replicate values
stored then, so no trajectory is loaded again:

```python
from pathlib import Path

import yaml

import polyzymd as pz
from polyzymd.analyses.functions import HBOND_PARTS, hydrogen_bonds
from polyzymd.config.comparison import PlotSettings

settings = PlotSettings(**yaml.safe_load(Path("plot_settings.yaml").read_text()))
study = pz.Study.from_configs(
    {"No Polymer": "noPoly/config.yaml", "100% SBMA": "SBMA100/config.yaml"},
    equilibration="10ns",
)
rows = study.per_replicate(
    hydrogen_bonds,
    pz.select("chainid A"),
    pz.select("chainid C"),
    unit=None,
    name="hydrogen_bonds_protein_polymer",
    parts=list(HBOND_PARTS),
)
rows["mean_hbonds"].plot("figures/hydrogen_bonds", "mean_hbonds", plot_settings=settings)
```

`plot` writes `figures/hydrogen_bonds/mean_hbonds.pdf`. `Timeseries.plot`,
`Timeseries.plot_distribution`, `pz.plot_values` and `pz.plot_distributions`
take `plot_settings` the same way; {doc}`../reference/analysis_functions`
lists the figures each analysis draws.

These fields are defined in PolyzyMD's `PlotSettings` model:

- `format`
  - Output file format
  - Allowed values: `png`, `pdf`, `svg`
- `dpi`
  - Plot resolution
  - Allowed range: `50` to `600`
  - Primarily affects raster output (`png`)
- `style`
  - PolyzyMD built-in theme preset for standard analysis plots
  - Allowed values: `compact`, `large_elements`, `low_ink`

`style` is a PolyzyMD theme preset selector, not a matplotlib or seaborn
stylesheet name. It does not choose the output `format`, `dpi`, figure sizes, or
condition color palettes. Use `format`, `dpi`, and `color_palette` or
`semantic_colors` for those choices.

The built-in presets are:

- `compact` — the default compact print-style preset.
- `large_elements` — larger fonts, markers, lines, error-bar caps, and fill
  opacity for slides or high-visibility output.
- `low_ink` — removes or reduces replicate dots, bar edges, and heavy reference
  lines for simpler figures.

Use `theme` to override individual visual values on top of the selected preset.
For example, this starts from `large_elements` and then changes only the listed
fields:

```yaml
format: "png"
dpi: 300
style: "large_elements"
theme:
  dot_size: 40
  title_fontsize: 20
  show_watermark: false
```

Standard PolyzyMD plots are designed for consistent screening and reporting.
Final manuscript figures may still need custom plotting from the stored
per-replicate values, {doc}`custom_artifact_plotting`, when a journal, panel
layout, or statistical annotation requires bespoke styling.

You can also set these optional global fields:

- `output_dir` (default: `figures/`), which the `plot` methods of the study
  API do not read: they take the folder as their first argument
- `color_palette` (default: `tab10`)
- `theme` (fine-grained visual overrides)

Example with all global fields:

```yaml
output_dir: "figures/"
format: "png"
dpi: 300
style: "compact"
color_palette: "tab10"
```

## Use semantic colors for condition series

Use `semantic_colors` when the condition colors should carry
experimental meaning, such as polymer family, composition, or dose. Semantic
colors are opt-in and apply to condition-series plots. Plot elements that are
not conditions, such as residue classes or secondary-structure categories, may
keep their categorical colors.

The conditions of this example are "No Polymer", the control and first
condition of the study, "100% SBMA", "100% EGMA", and "1% EGPMA" to "10% EGPMA":

```yaml
style: "compact"
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
`semantic_colors.conditions` and every entry in `semantic_colors.order` must
match a condition label of the study exactly.

For small percentage series such as `1`, `2`, `5`, and `10`, `scale: ordinal`
often gives more readable figures than a linear scale because each configured
level receives a visually distinct shade. For numeric gradients with many
levels, `scale: linear` is usually appropriate.

For multi-component chemistry, avoid naive additive RGB mixing: colors that add
component hues together can become hard to interpret and hard to reproduce.
Prefer a simpler mapping where hue identifies the family and shade identifies
the ordered composition or dose within that family.

After changing only plot settings, redraw the figures as in
[Write the settings and redraw a figure](#write-the-settings-and-redraw-a-figure);
the stored values are read back, not measured again.

## Changing output format

Switch file format by changing `format`.

### PNG

```yaml
format: "png"
dpi: 300
style: "compact"
```

Use PNG for quick drafts, sharing in chat/slides, and embedding in docs.

### PDF

```yaml
format: "pdf"
dpi: 300
style: "compact"
```

PDF is vector output, so it scales cleanly.

### SVG

```yaml
format: "svg"
dpi: 300
style: "compact"
```

SVG is vector output and is convenient for web use and post-editing.

## Changing DPI

Set `dpi` to control raster resolution.

```yaml
format: "png"
dpi: 150
style: "compact"
```

Common ranges:

- `72-150`: screen and draft output
- `300`: print/publication output
- `600`: high-resolution print output

## Figure options of `polyzymd analyze`

`polyzymd analyze` draws with the default settings. `--no-plots` turns its
figures off, and a few analyses take figure options with `--set`, such as
`highlight_residues` for `rmsf`. For any other setting on this page, redraw
the figure in Python with a `PlotSettings`.

## See Also

- [How to Compare Simulation Conditions](analysis_compare_conditions.md)
- [Create Custom Plots from Study Results](custom_artifact_plotting.md)
- [Shipped analysis functions](../reference/analysis_functions.md)
