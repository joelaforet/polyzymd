# SASA analysis: quick start

Measure the solvent-accessible surface area (SASA) of a protein, an active site
or any other set of atoms on every production frame of every replicate, with or
without neighbouring atoms such as polymer in the calculation, and compare
conditions with the replicate as the sampling unit.

```{versionadded} 1.3.0
SASA analysis was added in PolyzyMD 1.3.0.
```

```{note}
**Want to understand the measurement?** This guide focuses on getting results
quickly. For what each shipped function measures, see
{doc}`../reference/analysis_functions`; for how the values were checked, see
{doc}`../explanation/analysis_sasa_verification`; for the statistics, see
{doc}`../explanation/analysis_statistics_best_practices`.
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

## What is measured

For each frame, `mdtraj.shrake_rupley` places `n_sphere_points` points on a
sphere around every atom of the **context**, with radius the atom's radius
from MDTraj's element table plus the probe radius, and counts the points that
lie inside no other context atom's sphere. Each atom's SASA is its sphere's
area times the fraction of its points left free. The SASA of the **target** is
the sum over the target atoms. Atoms in the context but not in the target,
such as polymer, cover part of the target's surface without being counted.
Periodic images are not considered.

The elements come from the loaded universe, which PolyzyMD fills in from atom
types or names. Values are in Å².

## From the command line

```bash
polyzymd analyze sasa -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. By default the target is `protein`, measured on
its own under the name `isolated`. For each replicate, the per-frame total SASA
is averaged over the production frames. The replicate means are summarised per
condition, and every other condition is compared with the control by Welch's t
test with the Benjamini-Hochberg correction.

To see how much surface polymer covers, name the contexts to measure in:

```bash
polyzymd analyze sasa -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set "contexts={isolated: protein, with_polymer: protein or resname SBM EGM}" \
  --run with_polymer
```

Each context gives two results, picked with `--run`:

| `--run` | One value per replicate |
|---|---|
| `<context>` (default: the first context) | Mean over production frames of the target's total SASA |
| `<context>_residues` | Each target residue's SASA, averaged over production frames, compared residue by residue |

Only the result that `--run` picks is measured, because every context is a
separate Shrake-Rupley pass over every frame; run the command once per result
you need. A per-residue comparison is corrected over every residue of every
compared condition, and the text report gives, for each condition, how many
residues are significantly lower and higher than in the control and lists
them. Every per-residue row is kept in the JSON report.

Settings, passed with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `target` | `protein` | Atoms whose SASA is reported |
| `contexts` | `{isolated: <target>}` | Mapping of names to the selections of atoms present in the calculation; each must contain every target atom |
| `probe_radius_nm` | `0.14` | Probe radius in nm |
| `n_sphere_points` | `960` | Points on each atom's sphere |

For example, to ask whether polymer covers an active site:

```bash
polyzymd analyze sasa -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set "target=protein and (resid 77 or resid 156 or resid 262)" \
  --set "contexts={site_isolated: protein, site_with_polymer: protein or resname SBM EGM}" \
  --run site_with_polymer
```

Check the residue names of the polymer in the topology before relying on a
context selection:

```bash
python - <<'PY'
import MDAnalysis as mda

u = mda.Universe("solvated_system.pdb")
print(sorted(set(u.select_atoms("not protein").residues.resnames)))
PY
```

Add `--format json` for the full report, `--replicates 1-3` to use only some
replicates, and `--recompute` to ignore stored results. The settings, including
the contexts, are recorded under `provenance.settings` in the JSON report.

```{note}
Each frame is a separate MDTraj call, so a large context over a long trajectory
takes a while. Pass `--stride 5` to measure every fifth production frame, or
run the command inside a SLURM job on a cluster rather than on a login node.
```

```{note}
In the plugin used before this version, each run named its own target and
context with its own frame stride, and `polyzymd compare run sasa` computed
every run together. The stride is now `--stride`, shared by every result. Values from it are about 0.1 percent larger than
these: MDTraj 1.11.1 returns slightly more area for every frame after the
first one each thread computes in a call, and the plugin passed 100 frames per
call.
This version passes one frame per call. See
{doc}`../explanation/analysis_sasa_verification`.
```

## Figures

`polyzymd analyze sasa` writes these figures to `<output-dir>/figures/sasa/`;
`--no-plots` skips them.

| Figure | What it shows |
|---|---|
| `sasa_timeseries_<context>` | Every replicate's total SASA against time, with each condition's mean and its 95 percent interval |
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

`sasa(target, context)` returns one frame's total and runs through
`Study.timeseries`. `residue_sasa(target, context, frames)` returns each
target residue's mean over the frames and runs through `Study.per_replicate`.
Both take `probe_radius_nm` and `n_sphere_points` as keyword arguments.

## Next steps

- **What each function measures**: {doc}`../reference/analysis_functions`
- **How the values were checked**: {doc}`../explanation/analysis_sasa_verification`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **Contact analysis**: {doc}`analysis_contacts_quickstart`
