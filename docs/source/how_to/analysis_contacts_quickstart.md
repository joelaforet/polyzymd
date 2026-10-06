# Polymer-protein contacts analysis: quick start

Measure how often the polymer covers or touches each protein residue on the
production frames of every replicate, and compare how much of the protein the
polymer covers, overall, per monomer type, per amino-acid class, per region or
residue by residue, with the replicate as the sampling unit.

```{versionadded} 1.3.0
Contacts analysis was added in PolyzyMD 1.3.0.
```

```{note}
**Want to understand the measurement?** For what each shipped function
measures, see {doc}`../reference/analysis_functions`; for how the occlusion
values were checked, see {doc}`../explanation/analysis_contacts_verification`;
for the statistics, see {doc}`../explanation/analysis_statistics_best_practices`.
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

A protein residue is either in contact with the polymer on a frame or not.
The `method` setting picks what contact means.

**`method=occlusion`** (default): the polymer covers the residue's surface. On
every frame, each residue's solvent-accessible surface area (SASA) is computed
twice, with MDTraj's Shrake-Rupley code as in {doc}`analysis_sasa_quickstart`:
with the protein alone, and with the protein and the polymer, whose atoms
cover the protein without being counted. A residue's relative SASA is its SASA
divided by its maximum accessible surface area from Tien et al. (2013). The
residue is **exposed** on a frame when its relative SASA with the protein
alone is at least `exposed_threshold` (0.2), and in contact when it is exposed
and the polymer **buries** it: its relative SASA with the polymer is below
`buried_threshold` (0.2) and lower than without the polymer.
`exposed_threshold=0` counts every residue as exposed, so any residue the
polymer brings below `buried_threshold` is in contact, including one the
protein itself already partly buries. The polymer's SASA
loss on each residue, `max(0, alone - with)`, is also reported in Å².

Before the SASA with the polymer is computed, each polymer molecule (each
bonded fragment of the polymer selection) is moved whole by the box vector
that brings it to the periodic image nearest the protein, because the SASA
calculation does not consider periodic images. Molecules must be whole in the
trajectory, as OpenMM writes them. Residues with no maximum ASA, such as
terminal caps or non-standard residues, still cover their neighbours but are
not measured; a warning names them.

**`method=distance`**: the polymer touches the residue. On every frame,
MDAnalysis `lib.distances.capped_distance` finds every pair of a polymer atom
and a protein atom within `cutoff` (4.0 Å) of each other, using the minimum
image of the frame's box. Only heavy atoms are compared unless
`heavy_atoms=false`. A residue is in contact when any of its atoms is in such
a pair.

For either method, each residue's **contact fraction** is the fraction of
production frames it is in contact, and the same is measured for each polymer
residue name, such as each monomer type. For occlusion, a type's contact
fraction counts frames on which that type's atoms alone bury the residue. A
residue can be in contact with several types on one frame, so the per-type
fractions do not add up to the total.

PolyzyMD topologies put the protein on chain A and the polymer on chain C:

| Chain | Contents |
|-------|----------|
| A | Protein/enzyme |
| B | Substrate/ligand |
| C | Polymer |
| D+ | Solvent and ions |

The default selections, `chainid A` and `chainid C`, follow that convention.
Water and ions are never part of either calculation.

## From the command line

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml \
  --label "SBMA 50%" --label "SBMA 100%" --eq 200ns --stride 10
```

The first `-c` is the control. One pass over each replicate gives every result.
By default the report shows `coverage`: for each replicate, the fraction of
measured residues in contact on at least one production frame. The replicate
values are summarised per condition, and every other condition is compared with
the control by Welch's t test with the Benjamini-Hochberg correction. A
replicate where `polymer_selection` matches no atoms, such as every replicate
of a no-polymer control, has no contact: its values are 0 (0 events, no
lifetime), it is compared like any other, and a warning names it. A replicate
where `protein_selection` matches no atoms is left out with a warning, and a
`protein_selection` that matches no atoms in any replicate is refused.

Occlusion computes the SASA of the whole protein twice per frame, plus once
more per monomer type, so it takes about 2 s per frame for a 180-residue
lipase with 7,700 polymer atoms on four threads. `--stride 10` measures every
tenth production frame.

Pick another result with `--run`:

| `--run` | One value per replicate |
|---|---|
| `coverage` (default) | Fraction of measured residues in contact on at least one frame |
| `mean_contact_fraction` | Mean over residues of each residue's contact fraction |
| `<type>_contact_fraction` | The same for one polymer residue name `<type>`, for each type in the polymer, such as `SBM` or `EGM` |
| `<class>_contact_fraction` | Mean contact fraction of the residues of one amino-acid class: `aromatic`, `charged_positive`, `charged_negative`, `polar` or `nonpolar`, for the classes present |
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

A per-residue comparison is corrected over every residue of every compared
condition, and the text report gives, for each condition, how many residues
are significantly lower and higher than in the control and lists them. Every
per-residue row is kept in the JSON report.

Settings, passed with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `method` | `occlusion` | `occlusion` or `distance`, as above |
| `protein_selection` | `chainid A` | Protein atoms whose residues are measured |
| `polymer_selection` | `chainid C` | Polymer atoms |
| `polymer_types` | every residue name of the polymer in any condition | Monomers reported one by one, such as `[SBM]`; a replicate without one reports 0 for it. Narrow the polymer itself with `polymer_selection` |
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

A setting of the other method is refused. For example, to compare how much
SBMA buries the active site of two polymer conditions:

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml --eq 200ns \
  --stride 10 --set "regions={active_site: resid 77 or resid 133 or resid 156}" \
  --run active_site_contact_fraction
```

or to count contacts by distance instead:

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml --eq 200ns \
  --set method=distance --set cutoff=4.5
```

The resolved selections, the monomer types found, the unmeasured residues and
the residues of each amino-acid class and region are recorded under
`provenance.settings` in the JSON report. Add `--format json` for the full
report, `--replicates 1-3` to use only some replicates, and `--recompute` to
ignore stored results.

## How long contacts last

An event is a run of consecutive production frames in which one residue is in
contact, by the chosen method. `--run mean_lifetime` reports, for each
replicate, the Kaplan-Meier restricted mean duration of the events of all
measured residues, in ns: the estimate treats an event under way at the first
or last frame as lasting at least as long as observed, instead of counting it
as finished. `lifetime_events` and `censored_fraction` report how many events
there were and how many were cut off. The lifetime results need a pass over
the frames of their own and are measured only when chosen:

```bash
polyzymd analyze contacts -c SBMA50/config.yaml -c SBMA100/config.yaml --eq 200ns \
  --set method=distance --run mean_lifetime
```

The result depends on the spacing of the frames, since a contact or a break
shorter than the spacing is not seen, so compare conditions at the same frame
spacing and `--stride` only. `--set tolerance_ps=40` lets breaks of up to
40 ps continue an event; a tolerance changes the result strongly and should be
reported with it. {doc}`../explanation/analysis_contact_lifetimes` explains
the estimator, the censoring and these choices, with references.

## Figures

`polyzymd analyze contacts` writes these figures to
`<output-dir>/figures/contacts/`; `--no-plots` skips them.

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

`residue_occlusion(protein, occluder, frames, exposed_threshold=0.2,
buried_threshold=0.2, types=(), max_asa="theoretical", pbc=True)` returns the rows of `OCCLUSION_PARTS`,
`contact_fraction`, `exposed_fraction`, `occluded_area` and `exposed_area`, then
one contact-fraction row per residue name in `types`, with one column per
residue that has a maximum ASA. `residue_contacts(protein, polymer, frames,
cutoff=4.0, types=(), pbc=True)` returns one row of contact fractions, then one
per residue name in `types`, with one column per residue; pass heavy-atom
selections, such as `pz.select("chainid A and not element H")`, for the
default of `polyzymd analyze contacts --set method=distance`.

`contact_lifetimes(protein, polymer, frames, method="occlusion", types=(),
tolerance_ps=0.0, **options)` returns the rows of `LIFETIME_PARTS`,
`mean_lifetime`, `n_events` and `censored_fraction`, with one column for the
polymer and one per residue name in `types`; `options` are the settings of the
method, such as `cutoff` or `buried_threshold`. Run it with
`study.per_replicate(contact_lifetimes, ..., labels=["polymer", *types],
parts=list(LIFETIME_PARTS), types=types)`.

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
