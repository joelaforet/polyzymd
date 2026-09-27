# Shipped analysis functions

These functions ship in `polyzymd.analyses.functions`. Each one measures one
frame and runs through `Study.timeseries` like a function you write yourself;
see {doc}`../explanation/analysis_api`.

| Function | Arguments | Returns | Measurement | Command |
|---|---|---|---|---|
| `radius_of_gyration` | `atoms`, an `AtomGroup` | float, Å | MDAnalysis `AtomGroup.radius_of_gyration()`: mass-weighted, coordinates as loaded, molecules split across periodic boundaries not unwrapped | `polyzymd analyze rg` |
| `rmsd` | `atoms` and `reference`, `AtomGroup`s with the same atoms | float, Å | MDAnalysis `rms.rmsd(center=True, superposition=True)`: unweighted, `atoms` superposed on `reference` | `polyzymd analyze rmsd` |
| `pair_distance` | `atoms_a`, `atoms_b`, `mode_a` and `mode_b` (`single`, `midpoint`, `centroid` or `com`), `pbc` | float, Å | MDAnalysis `lib.distances.calc_bonds` between the two points, with the minimum image of the frame's box when `pbc` is true | `polyzymd analyze distances`, `polyzymd analyze catalytic_triad` |
| `all_below` | distance series, `thresholds` | 1 or 0 per frame | 1 for each frame in which every distance is strictly below its threshold; used with `Timeseries.transform` | `polyzymd analyze catalytic_triad` |

`polyzymd analyze rg` measures `radius_of_gyration` of `--set selection=...`
(default `protein`) on every production frame, reduces each replicate to its
mean, and compares every condition with the first by Welch's t test. See
{doc}`../how_to/analysis_rg_quickstart`.

`polyzymd analyze rmsd` measures `rmsd` of `--set selection=...` (default
`protein and name CA`) against
`pz.reference(reference_mode, selection, frame=reference_frame, file=reference_file, alignment=alignment_selection)`,
with `centroid` as the default mode. `reference_frame` counts production frames
from 1 after the equilibration window, and `alignment_selection` is the set of
atoms superposed to build the `average` and `centroid` references. Each frame
is superposed on the `selection` atoms themselves. See
{doc}`../how_to/analysis_rmsd_quickstart`.

`polyzymd analyze distances --set pairs=pairs.yaml` reads a YAML or JSON list of
pairs, each with `label`, `selection_a`, `selection_b` and optionally
`threshold` and `below_label`. It reports `<label>`, the mean distance, and
`<label> <below_label>`, the fraction of frames below the pair's threshold (the
analysis `threshold`, 3.5 Å by default, when the pair sets none).
`polyzymd analyze catalytic_triad` adds `simultaneous`, the fraction of frames in
which every pair is below its threshold, computed from the stored distances,
and reports it first. Pick a result with `--run`. See
{doc}`../how_to/analysis_distances_quickstart` and
{doc}`../how_to/analysis_triad_quickstart`.

## Figures

`polyzymd analyze` writes these figures to `<output-dir>/figures/<analysis>/`
and records the folder under `output_paths.figures` in the JSON report;
`--no-plots`, or `plots=False` in `analyze`, skips them.

| Analysis | Figures |
|---|---|
| `rg` | `rg_timeseries` (every replicate against time), `rg_comparison` (condition means with every replicate value), `rg_distribution` (per-frame distribution) |
| `rmsd` | `rmsd_timeseries`, `rmsd_comparison` |
| `distances` | `distance_kde_<pair>` for each pair with its threshold, and `distance_fraction_<result>` for each fraction below threshold |
| `catalytic_triad` | `triad_kde_<pair>` for each pair with its threshold, and `triad_fraction_<result>` for `simultaneous` and each pair's fraction |
