# Run RMSF analysis

Measure three quantities for each residue of each replicate:

- the RMSF: how much the residue fluctuates about its mean position;
- the offset: how far its mean position is from a reference structure;
- `rmsd_per_residue`: its root mean square deviation from the reference over
  the frames.

Then compare the conditions, with one value per replicate.

For what RMSF can and cannot show, see
{doc}`../explanation/analysis_rmsf_best_practices`. For the choice of
reference, see {doc}`../explanation/analysis_reference_selection`.

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## What is measured

PolyzyMD superposes each production frame on a reference structure by the
`alignment_selection` atoms. Then it measures these quantities for each atom
of `selection`:

| Quantity | Meaning | GROMACS equivalent |
|---|---|---|
| `rmsf` | The fluctuation about the mean position of the atom | `gmx rmsf -o` |
| `offset` | The distance of the mean position from the reference | none |
| `rmsd_per_residue` | The root mean square over frames of the distance from the reference | `gmx rmsf -od` |

For each residue, `rmsd_per_residue² = rmsf² + offset²`. PolyzyMD gives each
atom of a residue equal weight. GROMACS weights atoms by mass. The two give the
same result for Cα-only selections. The values agree with `gmx rmsf` to the
4 decimals in nm that GROMACS writes. See
{doc}`../explanation/analysis_rmsf_verification`.

`rmsd_per_residue` is one number per residue. `rmsd` is one number per frame.
For the definitions, and for which quantity answers which question, see
{ref}`Fluctuation, offset and deviation <rmsf-fluctuation-offset-deviation>`.

## From the command line

```bash
polyzymd analyze rmsf -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. One pass over each replicate gives all three
quantities. By default, `polyzymd analyze rmsf` reports the core RMSF.
`polyzymd analyze rmsd_per_residue` makes the same measurement, and by default
it reports the core `rmsd_per_residue`. To report a different result, use
`--run`:

| `--run` | One value per replicate |
|---|---|
| `core_rmsf`, `core_offset`, `core_rmsd_per_residue` | Root of the mean square over the core residues (below) |
| `<region>_rmsf`, `<region>_offset`, `<region>_rmsd_per_residue` | The same over one named region |
| `mean_rmsf`, `mean_offset`, `mean_rmsd_per_residue` | Plain mean over every selected residue |
| `rmsf`, `offset`, `rmsd_per_residue` | The per-residue profile, compared residue by residue |

PolyzyMD compares the results in two ways:

- A one-value result: PolyzyMD compares each condition with the control by
  Welch's t test. It corrects the p values over the compared conditions with
  the {term}`Benjamini-Hochberg` method.
- A profile: PolyzyMD compares each residue. The correction family is every
  residue of every compared condition.

For each condition, the text report gives the number of residues that are
significantly lower and higher than in the control, and lists them. The JSON
report keeps every per-residue row.

Change what is measured and what it is measured against with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `selection` | `protein and name CA` | Atoms measured |
| `alignment_selection` | `protein and name CA` | Atoms superposed on the reference |
| `reference_mode` | `external` if `reference_file` is set, else `centroid` | `centroid`, `average`, `frame` or `external` |
| `reference_frame` | `1` | Production frame used by `frame` mode, counted from 1 after the equilibration window |
| `reference_file` | none | Structure file used by `external` mode |
| `core` | all selected residues | MDAnalysis selection of the residues combined into the `core_*` values |
| `regions` | none | Mapping of region names to selections, each combined into `<region>_*` values |
| `highlight_residues` | none | Residue IDs marked on the profile figures |

For example, to measure against a crystal structure, with the lid as its own
region and the termini left out of the core:

```bash
polyzymd analyze rmsd_per_residue -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set reference_file=structures/1ISP.pdb \
  --set "core=resid 10:170" \
  --set "regions={lid: resid 70-90, termini: resid 1-9 or resid 171-181}"
```

Rules for `core` and `regions`:

- PolyzyMD intersects `core` and each region with `selection`.
- Each must select at least one residue.
- Each must select the same residues in every replicate.
- `core` and `mean` cannot be region names.

The JSON report records the residues of the core and of each region, the
reference mode used and every setting under `provenance.settings`.

```{note}
PolyzyMD superposes the frames by `alignment_selection`, not by the core. To
measure motion within the core, fit on the core atoms. For the core above, add
`--set "alignment_selection=name CA and resid 10:170"`. If you do not, a mobile
terminus in the fit adds apparent motion to the core.
```

```{note}
In `external` mode, the structure file must contain the atoms of `selection`
and `alignment_selection`. PolyzyMD stops with an error if the atom counts of
the trajectory and the file differ. PolyzyMD records the SHA-256 of the file
with the result. If you edit the file, the next command measures the
replicates again.
```

```{note}
`rmsf` is always the fluctuation about the mean position of each atom. The
deviation of a residue from a reference structure is `rmsd_per_residue`.
`reference_frame` counts production frames from 1, after the equilibration
window.
```

Useful options:

- `--format json` prints the full report.
- `--replicates 1-3` uses only some replicates.
- `--recompute` ignores stored results.

For the line format and the verdict words, see
{ref}`polyzymd analyze <cli-analyze>`.

## Figures

`polyzymd analyze rmsf` and `rmsd_per_residue` write these figures to
`<output-dir>/figures/<analysis>/`. Add `--no-plots` to skip them.

| Figure | What it shows |
|---|---|
| `rmsf_profile`, `offset_profile`, `rmsd_per_residue_profile` | The value of each replicate at each residue as a thin line, and the mean of each condition with its 95 % interval. The `highlight_residues` are marked |
| `rms_decomposition` | One panel per condition, with its mean `rmsd_per_residue`, RMSF and offset at every residue |
| `rmsf_difference`, `offset_difference`, `rmsd_per_residue_difference` | With several conditions, one panel per condition. Each panel shows the difference from the control at each residue, its 95 % interval, and a point on each significant residue |
| `rmsf_comparison` | The three core values of every condition, with every replicate value |

## From Python

`rms_decomposition` returns six rows per replicate. `RMS_PARTS` and `MS_PARTS`
hold the names of the rows:

- `rmsd_per_residue`, `rmsf` and `offset`, in Å;
- `ms_deviation`, `msf` and `ms_offset`, in Å². Each is the mean over the atoms
  of the residue of the squared deviation, RMSF or offset.

`parts=` gives each row its name, so `rows["rmsf"]` is one result.

```python
import polyzymd as pz
from polyzymd.analyses.functions import MS_PARTS, RMS_PARTS, rms_decomposition

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
ca = "protein and name CA"
rows = study.per_replicate(
    rms_decomposition,
    pz.select(ca),
    pz.select(ca),
    pz.reference("external", ca, file="structures/1ISP.pdb", alignment=ca),
    unit="Å",
    labels=lambda u: u.select_atoms(ca).residues.resids,
    parts=RMS_PARTS + MS_PARTS,
    bounds=(0.0, None),
)
profile = rows["rmsf"]                    # one value per residue per replicate
print(profile.compare(control="No polymer").to_agent_text())

mean_rmsf = profile.over_labels("mean")   # plain mean over the residues
lid = [r for r in profile.labels if 70 <= r <= 90]
lid_rmsf = rows["msf"].over_labels(lambda v: float(v.mean() ** 0.5), "lid_rmsf", lid)
```

A core or region value is the square root of the mean of a mean-square row
over its residues. This keeps `rmsd_per_residue² = rmsf² + offset²` true for
the whole set. `over_labels(how, metric, labels)` turns the profile of each
replicate into one number. It uses every label, or only the labels you give.
`polyzymd.analyses.functions.rmsf` and
`polyzymd.analyses.functions.rmsd_per_residue` each return one profile.

## Before interpreting the numbers

A lower RMSF shows that a region moves less. It does not show that the protein
is more stable. A replicate that partly unfolds contains two states. Its RMSF
then depends on when the unfolding happened, not only on equilibrium
fluctuation. For which claim each quantity supports, and how to choose the
core, see {doc}`../explanation/analysis_rmsf_best_practices`.

## Next steps

- **Understand RMSF interpretation**: {doc}`../explanation/analysis_rmsf_best_practices`
- **Choose a reference**: {doc}`../explanation/analysis_reference_selection`
- **RMSD analysis**: {doc}`analysis_rmsd_quickstart`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **Distance analysis**: {doc}`analysis_distances_quickstart`
