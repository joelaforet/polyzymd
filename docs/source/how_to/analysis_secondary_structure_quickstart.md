# Run secondary structure analysis

Assign a DSSP secondary-structure class to each residue on each production
frame of each replicate. Then compare how much of the protein is helix, strand
or coil, for the whole protein or residue by residue, with one value per
replicate.

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

On each production frame, `mdtraj.compute_dssp` assigns a DSSP class
(Kabsch and Sander, 1983) to each residue of the selection. DSSP uses the
backbone hydrogen bonds and geometry. The `scheme` setting selects the
classes:

- `simplified` (default): `helix`, `strand` and `coil`.
- `full`: the eight DSSP classes, such as `alpha_helix` and `3_10_helix`.

MDTraj gives `NA` to a residue that it cannot assign, for example a residue
without a backbone atom. PolyzyMD counts it as `unassigned` and prints a
warning. For each replicate, the value of a class for a residue is the fraction
of production frames in that class. For the table of classes, see
{ref}`DSSP classes <dssp-classes>`.

Before you run the analysis, check these points:

- The selection must contain whole residues, because DSSP reads every backbone
  atom.
- Each chain of the topology stays a separate chain, so DSSP never pairs
  residues of two chains.
- PolyzyMD uses the coordinates as loaded. Make a protein that crosses a
  periodic boundary whole before you use the assignment.

## From the command line

In a {term}`study`, give the study folder. The study names the conditions,
the control and the equilibration window. The results go to
`<study>/results/secondary_structure/`:

```bash
polyzymd analyze secondary_structure --study my_study
```

For a quick look without a study, give the `config.yaml` of each condition
instead. The results then go to the current folder:

```bash
polyzymd analyze secondary_structure -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. One pass over each replicate gives every class
of the scheme. By default, the scheme is `simplified` and the report shows
`helix`. This is the fraction of residue-frames in helix: the mean over
residues of the helix fraction of each residue. PolyzyMD compares each
condition with the control by Welch's t test. It corrects the p values with the
{term}`Benjamini-Hochberg` method.

To report a different result, use `--run`:

| `--run` | One value per replicate |
|---|---|
| a class of the scheme: `helix`, `strand`, `coil` or `unassigned` for `simplified` (default `helix`), or `alpha_helix` to `unassigned` for `full` (default `alpha_helix`) | Fraction of residue-frames in the class |
| `<class>_residues` | Each residue's fraction of frames in the class, compared residue by residue |

PolyzyMD corrects a per-residue comparison over every residue of every
compared condition. For each condition, the text report gives the number of
residues that are significantly lower and higher than in the control, and lists
them. The JSON report keeps every per-residue row.

Settings, passed with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `selection` | `protein` | The whole residues to assign |
| `scheme` | `simplified` | `simplified` for helix, strand and coil. `full` for the eight DSSP classes |

For example, this command compares how much of each protein is 3-10 helix:

```bash
polyzymd analyze secondary_structure -c noPoly/config.yaml -c SBMA50/config.yaml \
  --eq 200ns --set scheme=full --run 3_10_helix
```

This command compares the helix content of one domain, residue by residue:

```bash
polyzymd analyze secondary_structure -c noPoly/config.yaml -c SBMA50/config.yaml \
  --eq 200ns --set "selection=protein and resid 1:120" --run helix_residues
```

Useful options:

- `--stride 5` assigns every fifth production frame.
- `--format json` prints the full report.
- `--replicates 1-3` uses only some replicates.
- `--recompute` ignores stored results.

The JSON report records the selection and the scheme under
`provenance.settings`. PolyzyMD stores the results of the two schemes
separately.

## Figures

`polyzymd analyze secondary_structure` writes these figures to
`<output-dir>/figures/secondary_structure/`. Add `--no-plots` to skip them.

| Figure | What it shows |
|---|---|
| `ss_content_bars` | The fraction of every class of the scheme except unassigned, for every condition, with every replicate value |
| `ss_<name>_comparison` | For a total: each condition's mean with its interval and every replicate value |
| `ss_<name>_profile` | For `<name>_residues`: each residue's fraction per replicate and each condition's mean with its interval |
| `ss_classes_<name>` | For `<name>_residues`: one panel per condition with each residue's fraction of every class of the scheme except unassigned |
| `ss_<name>_difference` | For `<name>_residues` with several conditions: each condition minus the control at every residue, with the interval of the difference and the significant residues marked |

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import DSSP_SIMPLIFIED, dssp_occupancy

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
rows = study.per_replicate(
    dssp_occupancy,
    pz.select("protein"),
    unit=None,
    labels=lambda u: u.select_atoms("protein").residues.resids,
    parts=list(DSSP_SIMPLIFIED),
    bounds=(0.0, 1.0),
)
helix = rows["helix"].over_labels("mean", "helix_fraction")
print(helix.compare(control="No polymer").to_agent_text())
print(rows["helix"].compare(control="No polymer").to_agent_text())  # residue by residue
```

`dssp_occupancy(atoms, frames)` returns one row for each class of
`DSSP_SIMPLIFIED` (helix, strand, coil and unassigned), with one column for
each residue. For the eight classes, give `simplified=False` and
`parts=list(DSSP_CLASSES)`. `DSSP_GROUPS` lists the full classes that each
simplified class joins. `parts=` makes each row a separate result.
`over_labels("mean")` turns a residue profile into one value per replicate. To
measure your own quantity, see {doc}`study_api`.

## References

**Kabsch W, Sander C.** (1983) "Dictionary of protein secondary structure:
pattern recognition of hydrogen-bonded and geometrical features."
*Biopolymers* 22:2577-2637. https://doi.org/10.1002/bip.360221211

## Next steps

- **What each function measures**: {doc}`../reference/analysis_functions`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **RMSF analysis**: {doc}`analysis_rmsf_quickstart`
- **Contact analysis**: {doc}`analysis_contacts_quickstart`
