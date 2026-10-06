# Analysis System Concepts

PolyzyMD analyzes finished trajectories in four steps, for its shipped
analyses and for your own functions. This page explains the steps, the terms
of the output, and the reasons behind the choices.

## The analysis pipeline

Every analysis in PolyzyMD runs on the study API, in the same four steps:

```text
load replicates  →  measure  →  reduce per replicate  →  summarise or compare
```

| Step | Scope | What it does |
|------|-------|--------------|
| **load replicates** | Every replicate of every condition | `Study.from_configs` finds the replicates of each simulation `config.yaml`, loads each one as an MDAnalysis `Universe` with its production segments in order, and removes the equilibration window |
| **measure** | One replicate | `study.timeseries` runs a per-frame function on every production frame; `study.per_replicate` runs a function once on the replicate's production frames. Each replicate's result is written to `polyzymd_results/<name>/<condition>/replicate_<n>/` with a `record.json` of how it was made |
| **reduce per replicate** | One replicate | A time series becomes one number per replicate, for example its mean or the fraction of frames that pass a test |
| **summarise or compare** | All conditions | The replicate values give each condition's mean, standard error and 95 percent Student t interval, and each condition is compared with the control |

`polyzymd analyze NAME -c A/config.yaml -c B/config.yaml` runs these steps for
a shipped analysis and prints a report in which every number states its unit,
its uncertainty and its sample size. The first `-c` config is the control. The
same steps are open to your own functions in Python; see
{doc}`../how_to/study_api`.

Figures are drawn from the stored values into `figures/<name>/`. They do not
reload trajectories or rerun the measurement.

## Conditions, replicates and the control

A **condition** is one simulation setup, such as "No polymer", "SBMA 100%"
or "Urea 2 M". Each condition has its own `config.yaml`, which defines the
system: the protein, the polymer or co-solvent, the force field and the
temperature.

A **replicate** is one independent simulation of a condition. Replicates
have the numbers 1, 2, 3 and so on, and use the same `config.yaml`. The
replicate number seeds the starting structure. The starting velocities and
the thermostat noise are random in each simulation. These choices separate
the trajectories, but they do not make them statistically independent by
themselves. Independence also depends on equilibration, stationarity and
whether the simulated time is long enough for the process that you measure.

The replicate is the sampling unit: every interval and test uses one value
per replicate. With 2 conditions of 3 replicates each, a comparison rests on
3 values against 3 values, whatever the number of frames.

`polyzymd analyze` takes one `-c config.yaml` per condition. `--label` gives
each condition the name of the reports and figures. The default label is the
name of the folder of the config. `--replicates 1-3` selects the replicates;
the default is every replicate found. `--eq 10ns` sets the equilibration
window, removed from the start of the production trajectory of every
replicate. The first condition is the control. PolyzyMD compares every other
condition with it.

In Python the same inputs go to `Study.from_configs`:

```python
import polyzymd as pz

study = pz.Study.from_configs(
    {"No polymer": "no_polymer/config.yaml", "SBMA 100%": "sbma_100/config.yaml"},
    equilibration="10ns",
)
```

In folder names, PolyzyMD replaces characters that a file system cannot hold.
For example, `SBMA 100%` can become `SBMA_100_` in a path, and stays
`SBMA 100%` in reports and figures.

## The shipped analyses

Each shipped analysis is a function in `polyzymd.analyses.functions`, and
`polyzymd analyze` runs it with settings given as `--set key=value`:

| Name | What it measures |
|------|------------------|
| `rg` | Radius of gyration |
| `rmsd` | RMSD from a reference structure, one value per frame |
| `rmsf` | Per-residue fluctuations |
| `rmsd_per_residue` | Per-residue RMS deviation from a reference over the frames, as `gmx rmsf -od` gives |
| `sasa` | Solvent-accessible surface area |
| `secondary_structure` | DSSP secondary-structure fractions |
| `contacts` | Polymer-protein contacts and how long they last |
| `native_contacts` | Fraction of native contacts |
| `hydrogen_bonds` | Hydrogen-bond counts, lifetimes and per-residue occupancy |
| `distances` | Pair distances and the fraction of frames below a threshold |

The catalytic triad is not a separate analysis but a routine that combines
these functions; see {doc}`../how_to/analysis_triad_quickstart`. For what each
function measures, see {doc}`../reference/analysis_functions`.

## Periodic boundaries, whole molecules, and alignment

Molecular dynamics runs in a periodic box, and the coordinates written to a
trajectory are usually wrapped back into the primary cell. A molecule that
drifts across a face of the box is then stored in pieces, with some atoms at one
edge and the rest at the opposite edge. Nothing in the file says whether that
happened, so an observable computed from raw coordinates can be right on one
frame and badly wrong on the next.

Which observables care depends on what they measure. A radius of gyration is a
spread about a center of mass, so a molecule stored in two pieces reports a
radius of roughly half a box length rather than its real size. Alignment and
RMSD have the same problem, because fitting a structure that is split in two
fits a shape the molecule never adopted. A pair distance is different. The
minimum image convention gives the distance to the nearest periodic copy of the
second point, which is the physically meaningful separation whether or not the
molecules were wrapped, so distances need the box rather than whole molecules.
PolyzyMD therefore makes wholeness an explicit choice at load time through
`pbc_policy`, records the choice in provenance, and refuses to unwrap a topology
that has no bonds, because unwrapping walks the bond graph.

Alignment does not belong anywhere in a distance calculation, for two reasons.
The first is that a distance is invariant under rotation and translation, so
removing rigid-body motion cannot change a correct answer, and paying for a full
in-memory copy of the trajectory to do it buys nothing. The second is that it
was actively harmful here. MDAnalysis applies the minimum image convention by
building box vectors from the stored `(a, b, c, alpha, beta, gamma)` and
assuming those vectors describe the lattice of the coordinates it is handed.
Aligning in memory rotates every coordinate but copies the box unchanged, so the
two no longer agree. For a pair separated by less than half a box length in
every Cartesian component the wrapping is a no-op and the answer survives, which
is why catalytic-triad distances looked fine. For a long pair, such as a domain
center of mass against a polymer center of mass in a 90 Angstrom box, a
component can be folded against the wrong lattice vector and the reported
distance is wrong. Distances are therefore measured on the coordinates as the
trajectory stores them.

### References

- Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., and Beckstein, O. (2011).
  MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
  *Journal of Computational Chemistry*, 32(10), 2319-2327.
  doi:10.1002/jcc.21787
- MDAnalysis `MDAnalysis.lib.distances` documentation, on the `box` argument and
  the minimum image convention:
  <https://docs.mdanalysis.org/stable/documentation_pages/lib/distances.html>
- MDAnalysis `MDAnalysis.transformations.wrap` documentation, on `unwrap` and
  its bond requirement:
  <https://docs.mdanalysis.org/stable/documentation_pages/transformations/wrap.html>

## Why hydrogen bonds exclude carbon

MDAnalysis finds hydrogen bonds from geometry alone: it pairs each hydrogen with
a nearby heavy atom, then keeps the triplets whose donor-acceptor distance and
D-H...A angle pass the cutoffs. Which atoms are allowed to be donors and
acceptors is therefore a scientific choice, not a detail. If that choice is the
whole selection, every aliphatic and aromatic carbon that carries a hydrogen
becomes a donor and every atom becomes an acceptor, so C-H...O and N-H...C
contacts are counted alongside real hydrogen bonds.

The IUPAC definition requires the donor to be more electronegative than
hydrogen and the acceptor to carry a lone pair or a pi cloud, which carbon
generally does not (Arunan et al. 2011). `functions.hbond_atoms` therefore
takes as donors only the hydrogens covalently bonded to nitrogen, oxygen or
sulfur, and as acceptors every oxygen and every nitrogen or sulfur bonded to at
most two atoms, which keeps a lone pair. The practical reason matters as much
as the formal one: the share of short C-H...O geometries depends on polymer
chemistry, so counting them biases one condition relative to another instead of
shifting every condition by the same amount. How the rule was checked is in
{doc}`analysis_hydrogen_bonds_verification`.

Arunan, E., et al. (2011). Definition of the hydrogen bond (IUPAC
Recommendations 2011). Pure and Applied Chemistry, 83(8), 1637-1641.
doi:10.1351/PAC-REC-10-01-02

## Statistical comparison

When you have two or more conditions, every condition is compared with the
control, on one value per replicate:

- **Welch's t test** by default, or Student's t test with `test="student"` in
  Python, with the difference of the means and its 95 percent interval from
  the same test.
- **Benjamini-Hochberg correction** over one family per call: every tested row
  of that comparison. For a per-residue result the family is every residue of
  every condition compared, because scanning a profile for the residues that
  changed is a search over that set.
- **Effect sizes**, Cohen's d and Hedges' g, oriented like the difference, so
  you can see not just whether a difference is significant but how large it is.

A row is not testable when a condition has fewer than two replicates, or when
both conditions have the same value in every replicate; it then takes no part
in the correction. The report states a verdict per comparison in plain words.
For the report's fields, see {doc}`../reference/analysis_protocol_report`, and
for the choices behind the tests, {doc}`analysis_statistics_best_practices`.

## Output structure

`polyzymd analyze NAME -c ... --output-dir out` writes:

```text
out/
├── polyzymd_results/
│   └── <name>/
│       └── <condition>/
│           └── replicate_<n>/
│               ├── series.npz   # per-frame values, frames and times, or values.npz
│               └── record.json  # how the values were made
└── figures/
    └── <name>/
        └── *.<format>           # png by default
```

A per-frame measurement writes `series.npz`, a per-replicate one
`values.npz`. The `record.json` beside it holds the function's name, module and
source hash, the arguments with their selection strings, the config hash, the
topology and trajectory file records, the equilibration window, the frames, the
times, the unit and the PolyzyMD, MDAnalysis, NumPy and Python versions.
`--no-plots` skips the figures, and `--format json` prints the whole report.

## Why a cached result is checked against its inputs

A replicate that gained three segments since its result was computed will keep
reporting the old window for as long as the cache is reused, so a figure
regenerated a week later still describes last week's data. The failure is quiet
in the same way an unfinished segment is quiet: the cached number is a real
average over the frames that were read, and the provenance describes the run
rather than the subset, so nothing in the output says the two have drifted
apart.

Reuse is therefore conditional rather than automatic. A stored result is read
back only when every field of its `record.json`, except the versions and the
plot bounds, equals the record the new call would write: the same function
source, arguments, config, equilibration window, stride and frames, and a
topology and set of trajectory files with the recorded relative paths, sizes
and SHA-256 hashes. A file that grew, or changed in any other way, therefore
does not match. PolyzyMD reads a trajectory hash from `progress.json` or
`trajectory_hashes.json` when the run recorded it, or else computes it once
and caches it. Any other result is measured again, and `--recompute` measures
every replicate again whatever its record says.

## Why an unfinished segment is not read

A production segment that is still being written has a trajectory file that
ends at the last flush. If the loader reads it, the analysis window is shorter
than the run directory suggests, and the result records nothing about the
difference. The mean is a real average over fewer nanoseconds, the provenance
describes the whole run, and no later check can separate the two cases.

The OpenMM engine writes a status for each segment in `progress.json`.
The loader reads that status and skips segments marked `running` or
`failed`, and it lists the skipped indices in the provenance. Segments marked
`interrupted` are kept, because the continuation chain restarts from an
interrupted segment's saved state and its frames belong to the same time line.
The frame count of an interrupted segment is not recorded, so the last segment
of a chain may be partial. Dropping it instead would leave a gap between the
segments that were kept, and the lineage check would refuse the trajectory.

## See also

- {doc}`../tutorials/first_analysis` — Hands-on tutorial for running your
  first analysis
- {doc}`../how_to/analysis_compare_conditions` — Practical guide to comparing
  several conditions
- {doc}`../how_to/study_api` — Running your own functions on a study
- {doc}`analysis_api` — How the study API works
- {doc}`../reference/analysis_functions` — What each shipped function measures
