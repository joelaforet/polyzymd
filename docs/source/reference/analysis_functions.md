# Shipped analysis functions

These functions ship in `polyzymd.analyses.functions`. Each one measures one
frame and runs through `Study.timeseries` like a function you write yourself;
see {doc}`../explanation/analysis_api`.

| Function | Arguments | Returns | Measurement | Command |
|---|---|---|---|---|
| `radius_of_gyration` | `atoms`, an `AtomGroup` | float, Å | MDAnalysis `AtomGroup.radius_of_gyration()`: mass-weighted, coordinates as loaded, molecules split across periodic boundaries not unwrapped | `polyzymd analyze rg` |
| `rmsd` | `atoms` and `reference`, `AtomGroup`s with the same atoms | float, Å | MDAnalysis `rms.rmsd(center=True, superposition=True)`: unweighted, `atoms` superposed on `reference` | `polyzymd analyze rmsd` |

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
