# Polymer-protein contacts analysis: quick start

Measure how often each protein residue touches the polymer on the production
frames of every replicate, and compare how much of the protein the polymer
covers, overall, per monomer type, per amino-acid class, per region or residue
by residue, with the replicate as the sampling unit.

```{versionadded} 1.3.0
Contacts analysis was added in PolyzyMD 1.3.0.
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

On every production frame, MDAnalysis `lib.distances.capped_distance` finds
every pair of a polymer atom and a protein atom closer than the cutoff, 4.5 Å
by default, using the minimum image of the frame's box. A protein residue is in
contact on a frame when any of its atoms is in such a pair. Each residue's
**contact fraction** is the fraction of production frames it is in contact, and
the same is measured for each polymer residue name, such as each monomer type.
A residue can touch several monomer types on one frame, so the per-type
fractions do not add up to the total.

PolyzyMD topologies put the protein on chain A and the polymer on chain C:

| Chain | Contents |
|-------|----------|
| A | Protein/enzyme |
| B | Substrate/ligand |
| C | Polymer |
| D+ | Solvent and ions |

The default selections, `chainid A` and `chainid C`, follow that convention and
include hydrogen atoms.

## From the command line

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml \
  --label "SBMA 50%" --label "SBMA 100%" --eq 200ns
```

The first `-c` is the control. One pass over each replicate gives every result.
By default the report shows `coverage`: for each replicate, the fraction of
protein residues in contact on at least one production frame. The replicate
values are summarised per condition, and every other condition is compared with
the control by Welch's t test with the Benjamini-Hochberg correction. Every
condition needs a polymer: the selections are checked on the first replicate of
the control, and a selection that picks no atoms is refused. Pick another result
with `--run`:

| `--run` | One value per replicate |
|---|---|
| `coverage` (default) | Fraction of residues in contact on at least one frame |
| `mean_contact_fraction` | Mean over residues of each residue's contact fraction |
| `<type>_contact_fraction` | The same, counting only polymer residues named `<type>`, for each type in the polymer, such as `SBM` or `EGM` |
| `<class>_contact_fraction` | Mean contact fraction of the residues of one amino-acid class: `aromatic`, `charged_positive`, `charged_negative`, `polar` or `nonpolar`, for the classes present |
| `<region>_contact_fraction` | Mean contact fraction of the residues of one region of `regions` |
| `contact_fraction_residues` | Each residue's contact fraction, compared residue by residue |
| `<type>_contact_fraction_residues` | Each residue's contact fraction with one monomer type, compared residue by residue |

A per-residue comparison is corrected over every residue of every compared
condition, and the text report gives, for each condition, how many residues
are significantly lower and higher than in the control and lists them. Every
per-residue row is kept in the JSON report.

Settings, passed with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `protein_selection` | `chainid A` | Protein atoms whose residues are measured |
| `polymer_selection` | `chainid C` | Polymer atoms |
| `cutoff` | `4.5` | Contact distance in Å |
| `polymer_types` | none | Residue names to keep in the polymer selection, such as `[SBM]` |
| `use_pbc` | `true` | Use the minimum image of each frame's box |
| `regions` | none | Mapping of region names to selections, each reported as `<region>_contact_fraction`; a region cannot be named `coverage`, `mean`, `contact`, `classes`, a monomer type or an amino-acid class |

For example, to compare how much SBMA covers the active site of two polymer
conditions:

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml --eq 200ns \
  --set "regions={active_site: resid 77 or resid 133 or resid 156}" \
  --run active_site_contact_fraction
```

The resolved polymer selection, the monomer types found, and the residues of
each amino-acid class and region are recorded under `provenance.settings` in
the JSON report. Add `--stride 5` to measure every fifth production frame,
`--format json` for the full report, `--replicates 1-3` to use only some
replicates, and `--recompute` to ignore stored results.

```{note}
Contact residence times, and the `comparison.yaml` workflow, still run through
`polyzymd compare run contacts` in this version; see
{doc}`../reference/analysis_contacts_reference`. Its coverage and contact
fractions equal these for the same frames and settings.
```

## Figures

`polyzymd analyze contacts` writes these figures to
`<output-dir>/figures/contacts/`; `--no-plots` skips them.

| Figure | What it shows |
|---|---|
| `contacts_class_bars` | The mean contact fraction of each amino-acid class for every condition, with every replicate value |
| `contacts_<run>_comparison` | For a one-value result: each condition's mean with its interval and every replicate value |
| `contacts_<name>_profile` | For a residue result: each residue's contact fraction per replicate and each condition's mean with its interval |
| `contacts_<name>_difference` | For a residue result with several conditions: each condition minus the control at every residue, with the interval of the difference and the significant residues marked |

## From Python

```python
import numpy as np
import polyzymd as pz
from polyzymd.analyses.functions import residue_contacts

study = pz.Study.from_configs(
    {"SBMA 50%": "SBMA50/config.yaml", "SBMA 100%": "SBMA100/config.yaml"},
    equilibration="200ns",
)
rows = study.per_replicate(
    residue_contacts,
    pz.select("chainid A"),
    pz.select("chainid C"),
    unit=None,
    labels=lambda u: u.select_atoms("chainid A").residues.resids,
    parts=["contact_fraction", "SBM_contact_fraction"],
    bounds=(0.0, 1.0),
    types=["SBM"],
)
profile = rows["contact_fraction"]
coverage = profile.over_labels(lambda v: float(np.mean(np.asarray(v) > 0)), "coverage")
print(coverage.compare(control="SBMA 50%").to_agent_text())
print(profile.compare(control="SBMA 50%").to_agent_text())  # residue by residue
```

`residue_contacts(protein, polymer, frames, cutoff=4.5, types=(), pbc=True)`
returns one row of contact fractions per residue, then one row per residue name
in `types`.

## Next steps

- **What each function measures**: {doc}`../reference/analysis_functions`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **SASA analysis**: {doc}`analysis_sasa_quickstart`
- **Residence times and the comparison workflow**: {doc}`../reference/analysis_contacts_reference`
