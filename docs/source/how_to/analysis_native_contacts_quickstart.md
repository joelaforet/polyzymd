# Native contacts analysis: quick start

Measure the fraction of native contacts Q, how many of the contacts of a
reference structure are still formed, on every production frame of every
replicate, and compare each replicate's mean Q between conditions, for the
whole protein or for regions such as an active site.

```{versionadded} 1.3.0
Native contacts analysis was added in PolyzyMD 1.3.0.
```

```{note}
**Want to understand the measurement?** For what each shipped function
measures, see {doc}`../reference/analysis_functions`; for the statistics, see
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

The native contacts are the pairs of atoms of `selection` that are more than
`min_separation` residues apart in the topology and closer than `radius` Å in
the reference structure. On each frame, a pair whose atoms are `r` apart, and
`r0` apart in the reference, counts

$$
\frac{1}{1 + e^{\beta (r - \lambda r_0)}},
$$

the switching function of MDAnalysis `analysis.contacts.soft_cut_q`, which is
close to 1 while the pair stays within about `λ r0` and falls to 0 beyond it.
Q is the mean over the native pairs. The defaults, heavy atoms of the protein,
4.5 Å, more than 3 residues apart, β = 5 Å⁻¹ and λ = 1.8, are the definition of
Best, Hummer and Eaton (2013). Distances on each frame use the minimum image of
the box, and the reference must be a whole structure. Each replicate's value is
its mean Q over the production frames.

## From the command line

```bash
polyzymd analyze native_contacts -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. The report shows `q`, each replicate's mean Q,
summarised per condition, and every other condition is compared with the
control by Welch's t test with the Benjamini-Hochberg correction.

By default the reference is the first production frame of each replicate, so
Q measures how much of the structure present after equilibration is kept. To
measure against a crystal or prepared structure, give it as `reference_file`;
its atoms matching `selection` must be the replicate's atoms in the same order:

```bash
polyzymd analyze native_contacts -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set reference_file=structures/1ISP_prepared.pdb
```

Each entry of `regions` gives a result `<region>_q` over the native pairs with
at least one atom in the region; pick it with `--run`. Only the chosen result
is measured:

```bash
polyzymd analyze native_contacts -c noPoly/config.yaml -c SBMA50/config.yaml --eq 200ns \
  --set "regions={active_site: resid 76 77 132 155}" --run active_site_q
```

Settings, passed with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `selection` | `protein and not element H` | Atoms whose contacts are measured |
| `reference_mode` | `external` with `reference_file`, otherwise `frame` | Reference structure, as for {doc}`analysis_rmsd_quickstart`: `frame`, `external`, `average` or `centroid` |
| `reference_frame` | `1` | Production frame of the reference for `frame`, counted from 1 |
| `reference_file` | none | Structure file for `external` |
| `radius` | `4.5` | Å; native pairs are closer than this in the reference |
| `min_separation` | `3` | Native pairs are more than this many residues apart |
| `beta` | `5.0` | Softness of the switching function, in Å⁻¹ |
| `lambda_constant` | `1.8` | Tolerance on the reference distance; 1.5 is used for coarse-grained models |
| `use_pbc` | `true` | Use the minimum image of each frame's box |
| `regions` | none | Mapping of region names to selections, each reported as `<region>_q` |

The resolved reference mode is recorded under `provenance.settings` in the
JSON report, and for `external` the file's hash is recorded with each
replicate. Add `--stride 5` to measure every fifth production frame,
`--format json` for the full report, `--replicates 1-3` to use only some
replicates, and `--recompute` to ignore stored results.

## Figures

`polyzymd analyze native_contacts` writes `native_contacts_timeseries_<run>`,
Q against time for every replicate, and `native_contacts_comparison_<run>`,
each condition's mean with its interval and every replicate value, to
`<output-dir>/figures/native_contacts/`; `--no-plots` skips them.

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import native_contacts

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
heavy = "protein and not element H"
series = study.timeseries(
    native_contacts,
    pz.select(heavy),
    pz.reference("frame", heavy, frame=1),
    unit=None,
    bounds=(0.0, 1.0),
)
print(series.reduce("mean").compare(control="No polymer").to_agent_text())
```

`native_contacts(atoms, reference, region=None, radius=4.5, min_separation=3,
beta=5.0, lambda_constant=1.8, pbc=True)` returns Q at the current frame; pass
`pz.select(...)` as the third argument for a region. The native pairs are
found once per reference. For example, Cα atoms with `radius=8.0` give the
Cα-based Q some studies use with the same switching function.

## References

**Best RB, Hummer G, Eaton WA.** (2013) "Native contacts determine protein
folding mechanisms in atomistic simulations." *Proc Natl Acad Sci USA*
110:17874-17879. https://doi.org/10.1073/pnas.1311599110

## Next steps

- **What each function measures**: {doc}`../reference/analysis_functions`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **RMSD analysis**: {doc}`analysis_rmsd_quickstart`
- **Secondary structure analysis**: {doc}`analysis_secondary_structure_quickstart`
