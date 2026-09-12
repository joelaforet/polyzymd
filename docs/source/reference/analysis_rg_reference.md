# Rg plugin reference

For a step-by-step guide to running Rg analysis, see
{doc}`../how_to/analysis_rg_quickstart`.

The plugin is written against the observable contract, so it measures each
replicate and the framework owns aggregation, uncertainty, testing and
persistence. See {doc}`analysis_comparison_reference` for what the framework
does with each observable kind.

## Settings

`RgSettings` has one field, `runs`, a list of at least one `RgRunSettings`.
Labels must stay distinct after slugging, because the label names the
observable.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Run label, slugged into the observable name |
| `selection` | `str` | *required* | MDAnalysis selection string, no default |
| `calculation_mode` | `str` | `"selection"` | `"selection"` measures the whole group, `"fragments"` measures each bonded fragment |
| `fragment_weighting` | `str` | `"equal"` | How the per-frame mean over fragments weights each fragment, `"equal"` or `"mass"`. Only valid in fragment mode |
| `save_fragment_distribution` | `bool` | `true` | Report the distribution of fragment values. Fragment mode only |
| `histogram_bins` | `int` | `50` | Number of bins in that distribution, at least 2 |
| `histogram_range` | `[float, float]` or unset | `[0.0, 50.0]` | Lowest and highest fragment Rg the distribution covers, in Å |
| `allow_single_fragment_fallback` | `bool` | `false` | Measure the whole selection as one fragment when the topology has no bonds, instead of raising `TopologyBondsMissingError` |

```yaml
plugins:
  rg:
    runs:
      - label: "Protein"
        selection: "protein"
        calculation_mode: "selection"
      - label: "Polymer Oligomers"
        selection: "resname SBM EGM"
        calculation_mode: "fragments"
        fragment_weighting: "equal"
        save_fragment_distribution: true
        histogram_bins: 50
```

## Observables

A run labelled `Polymer Oligomers` slugs to `polymer_oligomers`. Every value is
the mass weighted radius of gyration MDAnalysis computes for an `AtomGroup`.

| Observable | Kind | Unit | Modes | Values |
|---|---|---|---|---|
| `rg_<label>` | `mean_of_timeseries` | `A` | both | Selection mode: the Rg of the whole group per frame. Fragment mode: the weighted mean over fragments per frame |
| `rg_<label>_fragments` | `profile` | `A` | fragments | Mean Rg of each fragment over the window, indexed by fragment number |
| `rg_<label>_distribution` | `profile` | `1/A` | fragments, unless `save_fragment_distribution` is off | Probability density of the fragment values, indexed by bin centre |

Each observable carries three metadata entries: `pbc_policy`, always
`as_loaded`; `topology_has_bonds`; and `bond_source`, one of `conect`,
`guessed` or `none`. A run that used the single-fragment fallback also carries
`fragment_fallback` with the reason. The framework copies them into the
replicate artifact.

The distribution uses fixed bin edges from `histogram_range`, not edges derived
from the data, because every replicate of a condition must report the same
profile index and a replicate cannot see its neighbours. A fragment value
outside the range raises rather than being dropped from the density; widen
`histogram_range` when that happens.

## Fragment mode requires topology bonds

Fragments are connected components of the bond graph, so fragment mode needs
the selected atoms to be bonded. The check is made against the selection, not
the topology as a whole, because a topology can carry bonds for the protein and
none for the polymer; in that case `AtomGroup.fragments` raises nothing and
returns one singleton fragment per atom. When the selected atoms have no usable
bonds the plugin raises
`polyzymd.analyses.exceptions.TopologyBondsMissingError`, naming the topology
file, its atom count and the two fixes: load a topology that carries bonds,
such as the OpenMM system XML read through ParmEd, or guess bonds for the
selection with `MDAnalysis.Universe(..., guess_bonds=True)`.

MDAnalysis skips CONECT records when a PDB holds atom serials above 99999,
which OpenMM writes in hexadecimal, so a solvated system above that size loads
without bonds even though the PDB contains CONECT lines.

Set `allow_single_fragment_fallback: true` to measure the whole selection as
one fragment instead, and say so when reporting the number.

## Coordinates and periodic boundaries

Rg is computed on the coordinates as they are loaded. No unwrap, centering or
make-whole transformation is applied, so a molecule split across a periodic
boundary reports an Rg of roughly half a box length rather than its real size.
Check whether your trajectory stores whole molecules. The GROMACS engine
prefers a whole-molecule trajectory when the run wrote one; OpenMM writes
wrapped coordinates.

## Why Rg has no alignment or reference fields

Rg is based on mass weighted distances from the centre of mass, so it is
translation and rotation invariant. Rg runs therefore have no
`alignment_selection`, `reference_mode`, `reference_file` or `reference_frame`.

## Output files

```text
<comparison_workspace>/
├── analysis/
│   └── <condition>/
│       └── rg/
│           ├── run_1/
│           │   ├── result.json
│           │   └── observables.npz
│           └── aggregated/
│               └── result.json
└── comparison/
    └── rg/
        └── result.json
```

| Level | Artifact | Path |
|-------|----------|------|
| Per replicate | `ReplicateArtifact` | `analysis/<condition>/rg/run_<replicate>/result.json` |
| Per replicate series | NPZ sidecar | `analysis/<condition>/rg/run_<replicate>/observables.npz` |
| Per condition | `ConditionArtifact` | `analysis/<condition>/rg/aggregated/result.json` |
| Cross condition | `ComparisonArtifact` | `comparison/rg/result.json` |

The replicate payload holds one `ObservableEstimate` per observable, the
condition payload one `ObservableAggregate`, and the comparison payload the
aggregates plus one `ObservableComparison` per observable and non-control
condition. The NPZ sidecar holds the full per-frame or per-index array of every
observable, keyed by observable name.

```python
from pathlib import Path

from polyzymd.analyses.mda import ArtifactStore

condition = ArtifactStore(Path("analysis/PEGylated/rg/aggregated")).read_condition_result()
for observable in condition.payload["observables"]:
    print(observable["name"], observable["mean"], observable["sem"], observable["unit"])
```

A replicate is recomputed when the plugin source, the settings, the simulation
config, the equilibration or the input files change, and reused otherwise.

## Plots

Figures come from the generic plotters keyed on observable kind, not from the
plugin. Until those land, `polyzymd compare run rg --plot` writes no figures;
the numbers and the artifacts are unaffected. There is no `plot_settings.rg`
block any more.

## Errors

| Message | Cause | Fix |
|---|---|---|
| `selection ... matched no atoms` | A run selection matches nothing | Check the selection against the topology. The run is no longer skipped with a warning |
| `in fragment mode needs topology bonds` | Fragment mode on an unbonded selection | Load a bonded topology, guess bonds, or opt into the fallback |
| `outside histogram_range` | A fragment Rg falls outside the distribution range | Widen `histogram_range` |
| `run labels must be unique after slugging` | Two labels slug to the same observable name | Rename one run |
| `asked for mass weighting but the topology gives fragment masses` | Zero or non-finite masses | Load a topology with masses |

## Rg against other metrics

| Feature | Rg | RMSD | RMSF |
|---|---|---|---|
| Measures | Compactness, mass weighted size | Deviation from a reference | Per-residue fluctuation |
| Reference | None | Required | Required for alignment |
| Main observable | `mean_of_timeseries` | `mean_of_timeseries` | `profile` over residues |
| Best question | Is the structure compacting or expanding? | Is it drifting from the reference? | Which regions are flexible? |
