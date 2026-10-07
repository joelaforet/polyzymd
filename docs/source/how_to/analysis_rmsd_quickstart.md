# Run RMSD analysis

Measure the RMSD of a selection from a reference structure on each production
frame of each replicate. Then compare the conditions, with one value per
replicate.

For how to read RMSD curves and choose a reference, see
{doc}`../explanation/analysis_rmsd_best_practices`. `rmsd` gives one number per
frame. For the per-residue quantities `rmsd_per_residue`, `rmsf` and `offset`,
see {ref}`Fluctuation, offset and deviation <rmsf-fluctuation-offset-deviation>`.

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## From the command line

In a {term}`study`, give the study folder. The study names the conditions,
the control and the equilibration window. The results go to
`<study>/results/rmsd/`:

```bash
polyzymd analyze rmsd --study my_study
```

For a quick look without a study, give the `config.yaml` of each condition
instead. The results then go to the current folder:

```bash
polyzymd analyze rmsd -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. The command does these steps:

1. It superposes the `protein and name CA` atoms of each production frame on a
   reference structure.
2. It measures the RMSD of each frame.
3. It takes the mean over frames of each replicate.
4. It compares each condition with the control by Welch's t test. It corrects
   the p values with the {term}`Benjamini-Hochberg` method.

Change what is measured and what it is measured against with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `selection` | `protein and name CA` | The atoms to measure. PolyzyMD superposes each frame on these atoms |
| `alignment_selection` | `protein and name CA` | Atoms superposed to build the `average` and `centroid` references |
| `reference_mode` | `external` if `reference_file` is set, else `centroid` | `centroid`, `average`, `frame` or `external` |
| `reference_frame` | `1` | Production frame used by `frame` mode, counted from 1 after the equilibration window |
| `reference_file` | none | Structure file used by `external` mode |

For example, to measure deviation from a crystal structure:

```bash
polyzymd analyze rmsd -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set reference_mode=external --set reference_file=structures/1ISP.pdb
```

```{note}
In `external` mode, the structure file must contain the atoms of `selection`.
PolyzyMD stops with an error if the atom counts of the trajectory and the file
differ. PolyzyMD records the SHA-256 of the file with the result. If you edit
the file, the next command measures the replicates again.
```

**Which reference to use:**

| Mode | Question answered |
|------|-------------------|
| `centroid` (default) | How much does the structure deviate from its most representative production frame? |
| `average` | How much does the structure deviate from its time-averaged conformation? |
| `frame` | How much does the structure deviate from a specific production frame? |
| `external` | How much does the structure deviate from a known functional geometry? |

The `centroid` reference is the production frame closest to the mean
structure. MDAnalysis `align.iterative_average` computes the mean structure on
the `alignment_selection` atoms. PolyzyMD builds the `average` and `centroid`
references for each replicate from the production frames of that replicate.

```{tip}
For an enzyme, run RMSD two times. Use `centroid` mode to measure overall
stability. Use `external` mode with the crystal structure to measure the
distance from the known active geometry.
```

Useful options:

- `--format json` prints the full report.
- `--replicates 1-3` uses only some replicates.
- `--recompute` ignores stored results.

For the line format and the verdict words, see
{ref}`polyzymd analyze <cli-analyze>`.

```{note}
`reference_frame` counts production frames from 1, after the equilibration
window, not trajectory frames from 0.
```

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import rmsd

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
ca = "protein and name CA"
values = study.timeseries(
    rmsd, pz.select(ca), pz.reference("centroid", ca, alignment=ca), unit="A"
)
print(values.reduce("mean").compare(control="No polymer").to_agent_text())
```

`pz.reference(mode, selection, frame=None, file=None, alignment=None)` gives
the reference atoms. PolyzyMD builds the reference once per replicate in a
separate universe, so the trajectory is not changed. The record of each
replicate holds the mode, the selections and the reference frame. For
`external`, it also holds the SHA-256 of the file.

To measure your own quantity, see {doc}`study_api`.

## Before interpreting the numbers

Look at the time series of each replicate before you choose the equilibration
window. If the RMSD still rises after the window, the replicate has not
relaxed. For the equilibration diagnostic in the report, see
{doc}`../explanation/convergence_detection`.

## Next steps

- **Understand RMSD interpretation**: {doc}`../explanation/analysis_rmsd_best_practices`
- **RMSF analysis**: {doc}`analysis_rmsf_quickstart`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **Distance analysis**: {doc}`analysis_distances_quickstart`
- **Contact analysis**: {doc}`analysis_contacts_quickstart`
