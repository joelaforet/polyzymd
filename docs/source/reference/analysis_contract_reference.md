# Observable contract reference

A plugin written against the observable contract is a settings model plus a
`compute(universe, frames, settings)` function that returns `Observable`
objects. Everything after that is framework work, including the figures. This
page describes what the framework generates from an observable and what a
plugin may change.

## Observable kinds

| kind | `values` hold | replicate value | figures |
|---|---|---|---|
| `mean_of_timeseries` | one number per frame | mean over frames | comparison bars, time series |
| `fluctuation` | one number per frame | sample standard deviation over frames | comparison bars, time series |
| `fraction` | 0 or 1 per frame | occupancy over frames | comparison bars on a `[0, 1]` axis |
| `profile` | one number per index | the vector itself | one line per condition with a band, or grouped bars |

## Generated figures

`polyzymd.analyses.contract_plots` renders the figures. The runner calls it
from the `plot(ctx)` method of the generated analysis class, so a plugin gets
them without writing any plotting code.

**Comparison bars.** One figure per scalar observable, one bar per condition,
the bar height being the mean over replicates. The error bar is the interval
named by `error_bar`, and the per-replicate values are drawn as jittered points
on top of it. A `fraction` observable is drawn on an axis pinned to `[0, 1]`.

**Time series.** One figure per `mean_of_timeseries` and `fluctuation`
observable, reading the per-frame values from the `observables.npz` sidecar
beside each replicate artifact. Each replicate contributes a faint trace in its
condition's colour, and the mean across replicates is drawn as a solid line.
The x axis is the frame index inside the production window.

**Profiles.** One figure per `profile` observable. When the index is integral
and holds no more than `max_categories_for_bars` entries, which is how residue
labels and pair labels arrive, the profile is drawn as grouped bars, one group
per category. Otherwise it is drawn as one line per condition over the index,
with a shaded band at the chosen interval and a faint trace per replicate.

Every figure that draws an error bar or a band carries a footnote naming the
interval, the number of replicates behind it and the production window, as
Grossfield et al. (2018) ask. The time-series panel draws no interval, so it
carries no such footnote.

### File names

Figures are written to the figures directory the comparison resolves, under the
analysis name, using the format in `plot_settings.format`.

| figure | file |
|---|---|
| comparison bars, and any profile | `<plugin>_<observable>_comparison.png` |
| time series | `<plugin>_<observable>_timeseries.png` |

## `ContractPlotSettings`

Every contract plugin exposes this model as its `PlotSettingsModel`, so a block
under the plugin's name in the `plot_settings` section of the comparison YAML
file is parsed into it.

| field | type | default | meaning |
|---|---|---|---|
| `error_bar` | `"ci95"` or `"sem"` | `"ci95"` | Interval drawn on bars and bands. `"ci95"` is the Student t interval across replicates, 4.303 times the standard error at n = 3. |
| `figsize` | `(float, float)` | `(10.0, 6.0)` | Width and height of every generated figure, in inches. |
| `show_replicates` | `bool` | `true` | Draw the per-replicate points on bars and the per-replicate traces on lines. |
| `max_categories_for_bars` | `int` | `30` | Longest categorical profile still drawn as grouped bars. |

```yaml
plot_settings:
  format: png
  rg2:
    error_bar: ci95
    show_replicates: false
    max_categories_for_bars: 20
```

A plugin that wants more fields declares a subclass and attaches it as
`PlotSettings`, which the runner uses in place of the default:

```python
class RgPlotSettings(ContractPlotSettings):
    show_histogram: bool = False


class RgContract:
    name: ClassVar[str] = "rg2"
    Settings: ClassVar[type[BaseModel]] = RgSettings
    PlotSettings: ClassVar[type[BasePlotSettings]] = RgPlotSettings
```

## References

Grossfield, A., Patrone, P. N., Roe, D. R., Schultz, A. J., Siderius, D. W. &
Zuckerman, D. M. (2018). Best practices for quantifying the uncertainty in
molecular simulations. *Living Journal of Computational Molecular Science*,
1(1), 5067. doi:10.33011/livecoms.1.1.5067
