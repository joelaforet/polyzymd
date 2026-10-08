# Run contacts analysis

Measure how often the polymer covers or touches each protein residue in each
replicate. Then compare the conditions, with one value per replicate. You can
compare the whole protein, each monomer type, each amino-acid class, each
region or each residue.

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

On each frame, a protein residue is in contact with the polymer or it is not.
The `method` setting selects the meaning of contact:

- `occlusion` (default): the polymer buries the residue. PolyzyMD computes the
  SASA of each residue with the protein alone and with the protein and the
  polymer. A residue is in contact if it is exposed without the polymer and
  buried with it.
- `distance`: a polymer atom is within `cutoff` (4.0 Å) of an atom of the
  residue.

The **contact fraction** of a residue is the fraction of production frames in
which it is in contact. PolyzyMD also measures it for each polymer residue
name, such as each monomer type. A residue can touch several types on one
frame, so the per-type fractions do not add up to the total.

For the steps of each method and the checks against an independent script,
see {doc}`../explanation/analysis_contacts_verification`. For the arguments of
each function, see {doc}`../reference/analysis_functions`.

Before you run the analysis, check these points:

- Molecules must be whole in the trajectory. OpenMM writes them whole.
  `occlusion` moves each polymer molecule to its periodic image nearest the
  protein, because the SASA calculation does not use periodic images.
- A residue with no maximum ASA, such as a terminal cap or a non-standard
  residue, still covers its neighbors. PolyzyMD does not measure it, and a
  warning names it.
- The default selections are null: they select the protein and the polymer by
  the PolyzyMD chain convention, the protein in chain A (`chainid A`) and the
  polymer in chain C (`chainid C`). The report records the selections used.
  Water and ions are never part of the calculation.

## From the command line

In a {term}`study`, give the study folder. The study names the conditions,
the control and the equilibration window. The results go to
`<study>/results/contacts/`:

```bash
polyzymd analyze contacts --study my_study
```

For a quick look without a study, give the `config.yaml` of each condition
instead. The results then go to the current folder:

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml \
  --label "SBMA 50%" --label "SBMA 100%" --eq 200ns --stride 10
```

The first `-c` is the control. One pass over each replicate gives every result.
By default, the report shows `coverage`. This is the fraction of measured
residues that are in contact on at least one production frame. PolyzyMD
compares each condition with the control by Welch's t test. It corrects the p
values with the {term}`Benjamini-Hochberg` method.

Empty selections:

- **A control without polymer.** If `polymer_selection` matches no atoms in a
  replicate, the replicate has no contact. Its values are 0, with 0 events and
  no lifetime. PolyzyMD keeps it in the statistics, compares it like any other
  replicate, and prints a warning that names it. So you can put a control
  without polymer first, as the control.
- **No protein.** If `protein_selection` matches no atoms in a replicate,
  PolyzyMD leaves the replicate out and prints a warning. If it matches no
  atoms in any replicate, PolyzyMD stops with an error.

`occlusion` computes the SASA of the whole protein two times per frame, and
one more time for each monomer type. For a 180-residue lipase with 7,700
polymer atoms on four threads, this takes about 2 s per frame. `--stride 10`
measures every tenth production frame.

To report a different result, use `--run`:

| `--run` | One value per replicate |
|---|---|
| `coverage` (default) | Fraction of measured residues in contact on at least one frame |
| `mean_contact_fraction` | Mean over residues of each residue's contact fraction |
| `<type>_contact_fraction` | The same for one polymer residue name `<type>`, for each type in the polymer, such as `SBM` or `EGM` |
| `<class>_contact_fraction` | Mean contact fraction of the residues of one amino-acid class: `aromatic`, `charged_positive`, `charged_negative`, `polar` or `nonpolar`, for the classes present. A residue name outside the 20 standard amino acids and their protonation variants, such as `SEP`, is in class `unknown`, so it gives `unknown_contact_fraction` |
| `<region>_contact_fraction` | Mean contact fraction of the residues of one region of `regions` |
| `occluded_area` | Occlusion only: SASA the polymer removes from the measured residues, in Å² per frame |
| `occlusion_fraction` | Occlusion only: that area over the residues' SASA with the protein alone, summed over frames |
| `contact_fraction_residues` | Each residue's contact fraction, compared residue by residue |
| `<type>_contact_fraction_residues` | Each residue's contact fraction with one monomer type, compared residue by residue |
| `occluded_area_residues` | Occlusion only: each residue's mean occluded area in Å², compared residue by residue |
| `mean_lifetime` | Kaplan-Meier restricted mean duration of a contact event, in ns; see [How long contacts last](#how-long-contacts-last) |
| `<type>_mean_lifetime` | The same for contacts with one monomer type |
| `lifetime_events` | Number of contact events |
| `censored_fraction` | Fraction of contact events cut off by the first or last production frame |

PolyzyMD corrects a per-residue comparison over every residue of every
compared condition. For each condition, the text report gives the number of
residues that are significantly lower and higher than in the control, and lists
them. The JSON report keeps every per-residue row.

Settings, passed with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `method` | `occlusion` | `occlusion` or `distance`, as above |
| `protein_selection` | null, the protein (`chainid A`) | Protein atoms whose residues are measured |
| `polymer_selection` | null, the polymer (`chainid C`) | Atoms of the partner group: the polymer, or any other group, such as a co-solvent (`resname SDS`) |
| `polymer_types` | every residue name of the polymer in any condition | The monomers to report one by one, such as `[SBM]`. A replicate without a listed monomer reports 0 for it. This setting does not narrow the polymer. To narrow the polymer, use `polymer_selection` |
| `use_pbc` | `true` | Use the frame's box: the minimum image for `distance`, and for `occlusion` each polymer molecule moved to its image nearest the protein |
| `regions` | none | Mapping of region names to selections, each reported as `<region>_contact_fraction`; a region cannot be named `coverage`, `mean`, `contact`, `classes`, `occluded`, `occlusion`, a monomer type or an amino-acid class |
| `exposed_threshold` | `0.2` | Occlusion only: relative SASA with the protein alone at or above which a residue is exposed; 0 counts every residue as exposed |
| `buried_threshold` | `0.2` | Occlusion only: relative SASA with the polymer below which an exposed residue is buried, and so in contact; above 0 |
| `max_asa` | `theoretical` | Occlusion only: the column of Tien et al. (2013) Table 1, `theoretical` (which the authors recommend) or `empirical` |
| `probe_radius_nm` | `0.14` | Occlusion only: SASA probe radius in nm |
| `n_sphere_points` | `960` | Occlusion only: points on each atom's sphere |
| `tolerance_ps` | `0` | Lifetime results only: absences of at most this many ps between two contacts do not end an event |
| `cutoff` | `4.0` | Distance only: contact distance in Å |
| `heavy_atoms` | `true` | Distance only: compare heavy atoms only |

PolyzyMD refuses a setting of the other method. For example, this command
compares how much SBMA buries the active site in two polymer conditions:

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml --eq 200ns \
  --stride 10 --set "regions={active_site: resid 76 or resid 132 or resid 155}" \
  --run active_site_contact_fraction
```

This command counts contacts by distance:

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml --eq 200ns \
  --set method=distance --set cutoff=4.5
```

This command measures contacts with a co-solvent, SDS, instead of the polymer.
`polymer_selection` takes any group of atoms:

```bash
polyzymd analyze contacts -c Water/config.yaml -c SDS/config.yaml --eq 200ns \
  --set "polymer_selection=resname SDS" --set method=distance
```

The JSON report records these items under `provenance.settings`: the
selections, the monomer types found, the residues not measured, and the
residues of each amino-acid class and region.

Useful options:

- `--format json` prints the full report.
- `--replicates 1-3` uses only some replicates.
- `--recompute` ignores stored results.

## How long contacts last

An event is a series of consecutive production frames in which one residue is
in contact. `--run mean_lifetime` reports the Kaplan-Meier restricted mean
duration of the events of all measured residues, in ns, for each replicate.
An event that continues at the first or last frame is censored: the estimate
counts it as at least as long as observed. `lifetime_events` gives the number
of events. `censored_fraction` gives the fraction of events that were cut off.

The lifetime results need their own pass over the frames. PolyzyMD measures
them only when you select them:

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml --eq 200ns \
  --set method=distance --run mean_lifetime
```

The result depends on the frame spacing. A contact or a break shorter than the
spacing is not seen. Compare conditions only at the same frame spacing and the
same `--stride`. `--set tolerance_ps=40` lets breaks of up to 40 ps continue an
event. A tolerance changes the result strongly. Report it with the result. For
the estimator and the censoring, see
{doc}`../explanation/analysis_contact_lifetimes`.

## Figures

`polyzymd analyze contacts` writes these figures to
`<output-dir>/figures/contacts/`. Add `--no-plots` to skip them.

| Figure | What it shows |
|---|---|
| `contacts_class_bars` | The mean contact fraction of each amino-acid class for every condition, with every replicate value |
| `contacts_<run>_comparison` | For a one-value result: each condition's mean with its interval and every replicate value |
| `contacts_<name>_profile` | For a residue result: each residue's value per replicate and each condition's mean with its interval |
| `contacts_<name>_difference` | For a residue result with several conditions: each condition minus the control at every residue, with the interval of the difference and the significant residues marked |

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import OCCLUSION_PARTS, residue_occlusion
from polyzymd.analyses.shared.aa_classification import get_max_asa

study = pz.Study.from_configs(
    {"SBMA 50%": "SBMA50/config.yaml", "SBMA 100%": "SBMA100/config.yaml"},
    equilibration="200ns",
    stride=10,
)
rows = study.per_replicate(
    residue_occlusion,
    pz.select("chainid A"),
    pz.select("chainid C"),
    unit=None,
    labels=lambda u: [
        r.resid for r in u.select_atoms("chainid A").residues if get_max_asa(r.resname)
    ],
    parts=[*OCCLUSION_PARTS, "SBM_contact_fraction"],
    bounds=(0.0, 1.0),
    types=["SBM"],
)
profile = rows["contact_fraction"]
print(profile.over_labels("mean", "mean_contact_fraction").compare(control="SBMA 50%").to_agent_text())
print(profile.compare(control="SBMA 50%").to_agent_text())  # residue by residue
```

For the arguments and return values of `residue_occlusion`,
`residue_contacts` and `contact_lifetimes`, see
{doc}`../reference/analysis_functions`. To measure your own quantity, see
{doc}`study_api`.

## References

**Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO.** (2013) "Maximum
allowed solvent accessibilities of residues in proteins." *PLoS ONE*
8:e80635. https://doi.org/10.1371/journal.pone.0080635

## Next steps

- **What each function measures**: {doc}`../reference/analysis_functions`
- **How the occlusion values were checked**: {doc}`../explanation/analysis_contacts_verification`
- **How long contacts last**: {doc}`../explanation/analysis_contact_lifetimes`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **SASA analysis**: {doc}`analysis_sasa_quickstart`
