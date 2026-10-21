# Run SASA analysis

Measure the solvent-accessible surface area (SASA) of a set of atoms, such as
the protein or an active site, on each production frame of each replicate. You
can include neighbor atoms, such as the polymer, in the calculation. Then
compare the conditions, with one value per replicate.

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

PolyzyMD computes the SASA with the Shrake-Rupley method of MDTraj. Two
selections control the calculation:

- The **target** is the set of atoms whose SASA PolyzyMD reports.
- A **context** is the set of atoms present in the calculation. It contains
  every target atom. Context atoms outside the target, such as polymer atoms,
  cover part of the target surface but are not counted.

The calculation does not use periodic images. Values are in Å². For the steps
of the calculation and the checks, see
{doc}`../explanation/analysis_sasa_verification`.

## From the command line

In a {term}`study`, give the study folder. The study names the conditions,
the control and the equilibration window. The results go to
`<study>/results/sasa/`:

```bash
polyzymd analyze sasa --study my_study
```

For a quick look without a study, give the `config.yaml` of each condition
instead. The results then go to the current folder:

```bash
polyzymd analyze sasa -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. By default, the target is `protein`, measured
alone, in a context named `isolated`. The command does these steps:

1. It measures the total SASA of the target on each production frame.
2. It takes the mean over frames of each replicate.
3. It compares each condition with the control by Welch's t test. It corrects
   the p values with the {term}`Benjamini-Hochberg` method.

To see how much surface the polymer covers, name the contexts:

```bash
polyzymd analyze sasa -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set "contexts={isolated: protein, with_polymer: protein or resname SBM EGM}" \
  --run with_polymer
```

Each context gives two results. Select one with `--run`:

| `--run` | One value per replicate |
|---|---|
| `<context>` (default: the first context) | Mean over production frames of the target's total SASA |
| `<context>_residues` | Each target residue's SASA, averaged over production frames, compared residue by residue |

PolyzyMD measures only the result that `--run` selects, because each context
is a separate pass over every frame. Run the command once for each result that
you need.

PolyzyMD corrects a per-residue comparison over every residue of every
compared condition. For each condition, the text report gives the number of
residues that are significantly lower and higher than in the control, and lists
them. The JSON report keeps every per-residue row.

Settings, passed with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `target` | `protein` | Atoms whose SASA is reported |
| `contexts` | `{isolated: <target>}` | A mapping of names to selections of the atoms present in the calculation. Each must contain every target atom |
| `probe_radius_nm` | `0.14` | The probe radius, in nm |
| `n_sphere_points` | `960` | The number of points on the sphere of each atom |

For example, this command measures how much the polymer covers an active
site:

```bash
polyzymd analyze sasa -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set "target=protein and (resid 76 or resid 132 or resid 155)" \
  --set "contexts={site_isolated: protein, site_with_polymer: protein or resname SBM EGM}" \
  --run site_with_polymer
```

Before you use a context selection, check the residue names of the polymer in
the topology:

```bash
python - <<'PY'
import MDAnalysis as mda

u = mda.Universe("solvated_system.pdb")
print(sorted(set(u.select_atoms("not protein").residues.resnames)))
PY
```

Useful options:

- `--format json` prints the full report.
- `--replicates 1-3` uses only some replicates.
- `--recompute` ignores stored results.

The JSON report records the settings, with the contexts, under
`provenance.settings`.

```{note}
PolyzyMD gives each frame to MDTraj in a separate call, so a large context over
a long trajectory takes time. Add `--stride 5` to measure every fifth
production frame. Or run the command in a SLURM job, not on a login node.
```

```{note}
`--stride` applies to every result of the command. PolyzyMD gives MDTraj one
frame per call because of a defect in MDTraj 1.11.1. In one call, MDTraj
returns about 0.1 % too much area for each frame after the first frame of each
thread. See {doc}`../explanation/analysis_sasa_verification`.
```

## Figures

`polyzymd analyze sasa` writes these figures to `<output-dir>/figures/sasa/`.
Add `--no-plots` to skip them.

| Figure | What it shows |
|---|---|
| `sasa_timeseries_<context>` | The total SASA of each replicate against time, with the mean of each condition and its 95 % interval |
| `sasa_comparison_<context>` | Each condition's mean with its interval and every replicate value |
| `sasa_distribution_<context>` | Each condition's distribution of per-frame values, pooled and per replicate |
| `sasa_profile_<context>` | For `<context>_residues`: each residue's SASA per replicate and each condition's mean with its interval |
| `sasa_difference_<context>` | For `<context>_residues` with several conditions: each condition minus the control at every residue, with the interval of the difference and the significant residues marked |

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import residue_sasa, sasa

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
protein, with_polymer = "protein", "protein or resname SBM EGM"

total = study.timeseries(
    sasa, pz.select(protein), pz.select(with_polymer), unit="A^2", bounds=(0.0, None)
)
print(total.reduce("mean").compare(control="No polymer").to_agent_text())

per_residue = study.per_replicate(
    residue_sasa,
    pz.select(protein),
    pz.select(with_polymer),
    unit="A^2",
    labels=lambda u: u.select_atoms(protein).residues.resids,
)
print(per_residue.compare(control="No polymer").to_agent_text())
```

`sasa(target, context)` returns the total of one frame. Use it with
`Study.timeseries`. `residue_sasa(target, context, frames)` returns the mean
over the frames of each target residue. Use it with `Study.per_replicate`.
Both take `probe_radius_nm` and `n_sphere_points` as keyword arguments. To
measure your own quantity, see {doc}`study_api`.

## Next steps

- **What each function measures**: {doc}`../reference/analysis_functions`
- **How the values were checked**: {doc}`../explanation/analysis_sasa_verification`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **Contact analysis**: {doc}`analysis_contacts_quickstart`
