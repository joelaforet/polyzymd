# RMSF analysis: quick start

Measure, for every residue, how much it fluctuates about its mean position
(RMSF), how far its mean position sits from a reference structure (the
offset), and its root mean square deviation from that reference, on every
replicate. Compare conditions with the replicate as the sampling unit.

```{versionadded} 1.3.0
RMSF analysis was added in PolyzyMD 1.3.0.
```

```{note}
**Want to understand the statistics?** This guide focuses on getting results
quickly. For what RMSF can and cannot support, how to combine residues into one
number, and which reference to choose, see
{doc}`../explanation/analysis_rmsf_best_practices` and
{doc}`../explanation/analysis_reference_selection`. For what each shipped
function measures, see {doc}`../reference/analysis_functions`.
```

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

```{tip}
**RMSD vs RMSF vs Distances, when to use which:**
- **RMSD**: Global structural deviation over time, "is the protein drifting?"
- **RMSF**: Per-residue fluctuation around average, "which residues are flexible?"
- **Distances**: Specific atom-pair distances, "is this H-bond intact?"
```

## What is measured

For each replicate, every production frame is superposed on a reference
structure by the `alignment_selection` atoms. Then, for each atom *i* of
`selection`:

| Quantity | Definition | GROMACS equivalent |
|---|---|---|
| `rmsf` | $\sqrt{\langle \lvert x_i(t) - \langle x_i \rangle \rvert^2 \rangle}$, the fluctuation about the atom's mean position | `gmx rmsf -o` |
| `offset` | $\lvert \langle x_i \rangle - x_i^{\mathrm{ref}} \rvert$, how far the mean position sits from the reference | none |
| `rms_deviation` | $\sqrt{\langle \lvert x_i(t) - x_i^{\mathrm{ref}} \rvert^2 \rangle}$, the deviation from the reference | `gmx rmsf -od` |

The values agree with `gmx rmsf` to the 4 decimals in nm that GROMACS writes;
{doc}`../explanation/analysis_rmsf_verification` gives the comparison and the
script that reruns it. The three are exact parts of one another: the squared deviation equals the
squared RMSF plus the squared offset for every atom. Each residue's value is
the mean over its atoms in `selection`. The reference decides what the frames
are superposed on and what the offset and deviation are measured from; `rmsf`
is always the fluctuation about the mean.

## From the command line

```bash
polyzymd analyze rmsf -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. One pass over each replicate gives all three
quantities. `polyzymd analyze rmsf` reports the core RMSF by default, and
`polyzymd analyze rms_deviation`, the same measurement, reports the core
deviation by default. Pick another result with `--run`:

| `--run` | One value per replicate |
|---|---|
| `core_rmsf`, `core_offset`, `core_rms_deviation` | Root of the mean square over the core residues (below) |
| `<region>_rmsf`, `<region>_offset`, `<region>_rms_deviation` | The same over one named region |
| `mean_rmsf`, `mean_offset`, `mean_rms_deviation` | Plain mean over every selected residue |
| `rmsf`, `offset`, `rms_deviation` | The per-residue profile, compared residue by residue |

A one-value result is summarised per condition, and every other condition is
compared with the control by Welch's t test with the Benjamini-Hochberg
correction over the compared conditions. A profile is compared at every
residue, and its correction family is every residue of every compared
condition. The text report gives, for each condition, how many residues are
significantly lower and higher than in the control, and lists them. Every
per-residue row is kept in the JSON report.

Change what is measured and what it is measured against with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `selection` | `protein and name CA` | Atoms measured |
| `alignment_selection` | `protein and name CA` | Atoms superposed on the reference |
| `reference_mode` | `centroid` for `rmsf`; for `rms_deviation`, `external` when `reference_file` is set and `centroid` otherwise | `centroid`, `average`, `frame` or `external` |
| `reference_frame` | `1` | Production frame used by `frame` mode, counted from 1 after the equilibration window |
| `reference_file` | none | Structure file used by `external` mode |
| `core` | all selected residues | MDAnalysis selection of the residues combined into the `core_*` values |
| `regions` | none | Mapping of region names to selections, each combined into `<region>_*` values |
| `highlight_residues` | none | Residue IDs marked on the profile figures |

For example, to measure against a crystal structure, with the lid as its own
region and the termini left out of the core:

```bash
polyzymd analyze rms_deviation -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set reference_file=structures/1ISP.pdb \
  --set "core=resid 10:170" \
  --set "regions={lid: resid 70-90, termini: resid 1-9 or resid 171-181}"
```

`core` and each region are intersected with `selection`, must pick at least
one residue, and must pick the same residues in every replicate. `core` and
`mean` cannot be region names. The residues of the core and of each region,
the resolved reference mode and every setting are recorded under
`provenance.settings` in the JSON report.

```{note}
The frames are superposed by `alignment_selection`, whatever the core. To
measure motion within the core, fit on the same atoms, for example
`--set "alignment_selection=name CA and resid 10:170"` together with the core
above. Otherwise a mobile terminus in the fit adds apparent motion to the core.
```

```{note}
When using `external` reference mode, the structure file must contain the atoms
of `selection` and `alignment_selection`. PolyzyMD checks that the atom counts
match between the trajectory and the file and raises an error on a mismatch.
The file's SHA-256 hash is recorded with the result, so editing the file
measures the replicates again.
```

```{note}
`rmsf` is always the fluctuation about each atom's mean position; each
residue's deviation from a reference structure is `rms_deviation`.
`reference_frame` counts production frames from 1, after the equilibration
window.
```

Add `--format json` for the full report, `--replicates 1-3` to use only some
replicates, and `--recompute` to ignore stored results. The line format and the
verdict words are described under {ref}`polyzymd analyze <cli-analyze>`.

## Figures

`polyzymd analyze rmsf` and `rms_deviation` write these figures to
`<output-dir>/figures/<analysis>/`; `--no-plots` skips them.

| Figure | What it shows |
|---|---|
| `rmsf_profile`, `offset_profile`, `rms_deviation_profile` | Each replicate's value at every residue as a thin line, and each condition's mean with its 95 percent interval; `highlight_residues` are marked |
| `rms_decomposition` | One panel per condition, with its mean deviation, RMSF and offset at every residue |
| `rmsf_difference`, `offset_difference`, `rms_deviation_difference` | With several conditions, one panel per condition: its difference from the control at every residue, the 95 percent interval of the difference, and a point on each significant residue |
| `rmsf_comparison` | The three core values of every condition, with every replicate value |

## From Python

`rms_decomposition` returns six rows per replicate. `RMS_PARTS` and `MS_PARTS`
are the names of those rows: `rms_deviation`, `rmsf` and `offset` in Å, then
their per-residue mean squares `ms_deviation`, `msf` and `ms_offset` in Å².
`parts=` gives each row its name, so `rows["rmsf"]` is one result. What the
three quantities are, and how they add up, is explained in
[Fluctuation, offset and deviation](../explanation/analysis_rmsf_best_practices.md#fluctuation-offset-and-deviation).

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
    unit="A",
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

`rms_decomposition` returns six rows per replicate: the per-residue means
`rms_deviation`, `rmsf` and `offset`, then their per-residue mean squares
`ms_deviation`, `msf` and `ms_offset`. A core or region value is the root of
the mean of a mean-square row over its residues, which keeps the identity
deviation² = RMSF² + offset² for the whole set. `over_labels(how, metric,
labels)` turns each replicate's profile into one number, over every label or
only the ones you pass. `polyzymd.analyses.functions.rmsf` and
`polyzymd.analyses.functions.rms_deviation` return a single profile each.

## Before interpreting the numbers

A lower RMSF says a region moves less. It does not by itself say the protein
is more stable. A replicate that partly unfolds mixes states, which raises its
RMSF according to when the unfolding happened rather than through an
equilibrium fluctuation. {doc}`../explanation/analysis_rmsf_best_practices`
covers which claim each quantity supports and how to choose the core.

## Next steps

- **Understand RMSF interpretation**: {doc}`../explanation/analysis_rmsf_best_practices`
- **Choose a reference**: {doc}`../explanation/analysis_reference_selection`
- **RMSD analysis**: {doc}`analysis_rmsd_quickstart`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **Distance analysis**: {doc}`analysis_distances_quickstart`
