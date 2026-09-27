# RMSD analysis: quick start

Measure the RMSD of a selection from a reference structure on every production
frame of every replicate, and compare conditions with the replicate as the
sampling unit.

```{versionadded} 1.3.0
RMSD analysis was added in PolyzyMD 1.3.0.
```

```{note}
**Want to understand the statistics?** This guide focuses on getting results
quickly. For interpretation of RMSD curves and the choice of reference, see
{doc}`../explanation/analysis_rmsd_best_practices`. For what each shipped
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

## From the command line

```bash
polyzymd analyze rmsd -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. For each replicate, every production frame of
the `protein and name CA` atoms is superposed on a reference structure and its
RMSD is measured. The replicate means are then summarised per condition, and
every other condition is compared with the control by Welch's t test with the
Benjamini-Hochberg correction.

Change what is measured and what it is measured against with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `selection` | `protein and name CA` | Atoms whose RMSD is measured; each frame is superposed on these atoms |
| `alignment_selection` | `protein and name CA` | Atoms superposed to build the `average` and `centroid` references |
| `reference_mode` | `centroid` | `centroid`, `average`, `frame` or `external` |
| `reference_frame` | `1` | Production frame used by `frame` mode, counted from 1 after the equilibration window |
| `reference_file` | none | Structure file used by `external` mode |

For example, to measure deviation from a crystal structure:

```bash
polyzymd analyze rmsd -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set reference_mode=external --set reference_file=structures/1ISP.pdb
```

```{note}
When using `external` reference mode, the structure file must contain atoms
matching the `selection` string. PolyzyMD checks that the atom counts match
between the trajectory and the file and raises an error on a mismatch. The
file's SHA-256 hash is recorded with the result, so editing the file measures
the replicates again.
```

**Which reference to use:**

| Mode | Question answered |
|------|-------------------|
| `centroid` (default) | How much does the structure deviate from its most representative production frame? |
| `average` | How much does the structure deviate from its time-averaged conformation? |
| `frame` | How much does the structure deviate from a specific production frame? |
| `external` | How much does the structure deviate from a known functional geometry? |

The `centroid` reference is the production frame closest to the mean
structure, where the mean comes from MDAnalysis `align.iterative_average` on
the `alignment_selection` atoms. The `average` and `centroid` references are
built separately for every replicate from its own production frames.

```{tip}
For enzyme studies, consider running RMSD twice: once with `centroid` mode
(overall stability) and once with `external` mode pointing to a crystal
structure (catalytic competence). These answer complementary questions.
```

Add `--format json` for the full report, `--replicates 1-3` to use only some
replicates, and `--recompute` to ignore stored results. The line format and the
verdict words are described under {ref}`polyzymd analyze <cli-analyze>`.

```{note}
In the plugin used before this version, `reference_frame` counted trajectory
frames from 0, including the equilibration window. It now counts production
frames from 1, so a value copied from an old `comparison.yaml` points at a
different frame.
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
    rmsd, pz.select(ca), pz.reference("centroid", ca, alignment=ca), unit="Å"
)
print(values.reduce("mean").compare(control="No polymer").to_agent_text())
```

`pz.reference(mode, selection, frame=None, file=None, alignment=None)` stands
for the reference atoms. It is built once per replicate in a separate
universe, so the trajectory itself is never modified. The record of each
replicate holds the mode, the selections, the frame the reference used and,
for `external`, the file's hash.

## Before interpreting the numbers

Check the per-replicate time series before choosing the equilibration window;
RMSD that is still rising after the window means the replicate has not
relaxed. {doc}`../explanation/analysis_rmsd_best_practices` covers how to read
RMSD curves and choose a reference, and
{doc}`../explanation/convergence_detection` covers the equilibration
diagnostic in the report.

## Next steps

- **Understand RMSD interpretation**: {doc}`../explanation/analysis_rmsd_best_practices`
- **RMSF analysis**: {doc}`analysis_rmsf_quickstart`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **Distance analysis**: {doc}`analysis_distances_quickstart`
- **Contact analysis**: {doc}`analysis_contacts_quickstart`
