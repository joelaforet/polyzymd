# Run hydrogen bonds analysis

Count the hydrogen bonds between named groups of atoms, or within one group,
on the production frames of every replicate, and compare how many form, how
long they last, and which residues and residue pairs form them, with the
replicate as the sampling unit.

```{note}
**Want to understand the measurement?** For what each shipped function
measures, see {doc}`../reference/analysis_functions`; for how the donors and
acceptors are chosen and the values checked, see
{doc}`../explanation/analysis_hydrogen_bonds_verification`; for how lifetimes
are estimated, see {doc}`../explanation/analysis_contact_lifetimes`.
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

On every frame, the analysis counts the hydrogen bonds with MDAnalysis
`HydrogenBondAnalysis`: a donor within 3.5 Å of an acceptor and a
donor-hydrogen-acceptor angle of at least 150°. It chooses the donors and
acceptors from the covalent bonds of the force field. For the rule and the
files it reads, see {doc}`../explanation/analysis_hydrogen_bonds_verification`.

```{important}
Keep these files beside the trajectory in any copy of a run you analyse, such
as a copy thinned to fewer frames: `<segment>_system.xml` for OpenMM, and
`prod.tpr` with `<prefix>.top` and its `.itp` files for GROMACS. Without them
the protein has no bonds, and the analysis stops with an error naming the
files. A run without them can still be analysed by giving `donors`,
`hydrogens` and `acceptors` as selections, which pairs each hydrogen with a
donor within 1.2 Å instead of by bonds.
```

## From the command line

In a {term}`study`, give the study folder. The study names the conditions,
the control and the equilibration window. The results go to
`<study>/results/hydrogen_bonds/`:

```bash
polyzymd analyze hydrogen_bonds --study my_study
```

For a quick look without a study, give the `config.yaml` of each condition
instead. The results then go to the current folder:

```bash
polyzymd analyze hydrogen_bonds -c SBMA50/config.yaml -c SBMA100/config.yaml \
  --label "SBMA 50%" --label "SBMA 100%" --eq 200ns
```

The first `-c` is the control. The analysis measures the summaries of
`--set summaries=...`, by default `protein_polymer`, the hydrogen bonds
between the protein (chain A) and the polymer (chain C). By default the report shows
`protein_polymer_mean_hbonds`: for each replicate, the mean number of hydrogen
bonds per frame. The replicate values are summarised per condition, and every
other condition is compared with the control by Welch's t test with the
Benjamini-Hochberg correction.
A replicate where a summary's second group matches no atoms, such as every
replicate of a no-polymer control for `protein_polymer`, has no hydrogen bond
with it: its values are 0 (0 events, no lifetime), and a warning names it.
That 0 is not a measurement, so every comparison with its condition is
`not testable`, with the reason. A replicate where the first group matches no
atoms is left out with a warning, and a first group that matches no atoms in
any replicate is refused.

Groups and summaries are named. `groups` maps a name to an MDAnalysis
selection, and each summary is `between: [a, b]`, the bonds with one partner
in each group, in either direction, or `within: a`:

```bash
polyzymd analyze hydrogen_bonds -c SBMA50/config.yaml -c SBMA100/config.yaml --eq 200ns \
  --set "groups={protein: chainid A, substrate: chainid B, polymer: chainid C}" \
  --set "summaries={protein_substrate: {between: [protein, substrate]}, protein_protein: {within: protein}}" \
  --run protein_protein_mean_hbonds
```

Each summary `<s>` gives these results; pick one with `--run`. Only the chosen
summary is measured:

| `--run` | One value per replicate |
|---|---|
| `<s>_mean_hbonds` | Mean hydrogen bonds per frame, each donor-hydrogen-acceptor counted once |
| `<s>_mean_residue_pairs` | Mean residue pairs joined by at least one hydrogen bond per frame |
| `<s>_any_fraction` | Fraction of frames with at least one hydrogen bond |
| `<s>_mean_lifetime` | Kaplan-Meier restricted mean duration of a hydrogen-bond event, in ns |
| `<s>_lifetime_events` | Number of hydrogen-bond events |
| `<s>_censored_fraction` | Fraction of events cut off by the first or last production frame |
| `<s>_residues` | Each residue of the summary's first group: fraction of frames with a hydrogen bond, compared residue by residue |
| `<s>_pairs` | Each residue pair joined on some frame: fraction of frames joined, compared pair by pair |

A lifetime event is a run of consecutive frames in which one pair is joined.
With `lifetime_key=residue`, the default, a pair is two residues, so a
hydrogen switching between atoms of the same two residues does not end the
event; with `lifetime_key=atom`, a pair is one donor atom and one acceptor
atom. As for contacts, the lifetime depends on the frame spacing and on
`tolerance_ps`; see {doc}`../explanation/analysis_contact_lifetimes`.

In `<s>_pairs`, a pair is a protein residue and a monomer type: `149-SBM` is
residue 149 with any SBM monomer. A pair that never forms in a replicate has
the value 0 there. A per-residue or per-pair comparison is corrected over
every entry of every compared condition, and the JSON report keeps every row.
For why pairs are named this way, see
{doc}`../explanation/analysis_hydrogen_bonds_verification`.

Settings, passed with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `groups` | `{protein: null, polymer: null}` | Names and MDAnalysis selections of the groups. A null selection selects the atoms of the role the group is named after: `protein` (chain A), `ligand` (chain B) or `polymer` (chain C). The report records the selection used |
| `summaries` | `{protein_polymer: {between: [protein, polymer]}}` | Names of the summaries, each `between: [a, b]` or `within: a` |
| `d_a_cutoff` | `3.5` | Largest donor-acceptor distance, in Å; 3.0 is also common |
| `d_h_a_angle_cutoff` | `150` | Smallest donor-hydrogen-acceptor angle, in degrees |
| `donors` | none | Selection of donors, which pairs hydrogens with donors within 1.2 Å instead of by bonds |
| `hydrogens` | none | Selection of donatable hydrogens, in place of the element and valency rule |
| `acceptors` | none | Selection of acceptors, in place of the element and valency rule |
| `lifetime_key` | `residue` | Lifetime results only: `residue` pairs or `atom` pairs |
| `tolerance_ps` | `0` | Lifetime results only: absences of at most this many ps do not end an event |

The explicit selections are taken within the summary's groups. For example,
to count only side-chain oxygens as acceptors:

```bash
polyzymd analyze hydrogen_bonds -c SBMA50/config.yaml -c SBMA100/config.yaml --eq 200ns \
  --set "acceptors=element O and not name O OXT"
```

Add `--stride 5` to measure every fifth production frame, `--format json` for
the full report, `--replicates 1-3` to use only some replicates, and
`--recompute` to ignore stored results.

(what-is-recorded)=
## What is recorded

`provenance.settings` in the JSON report holds the settings, the summary and
its selections, and under `hbond_atoms` the donors and acceptors chosen, as a
count per residue name and atom name, with the number of donatable hydrogens;
an excerpt from a lipase with an SBMA-EGMA polymer:

```json
"hbond_atoms": {
  "donors": {"ARG NE": 5, "HIS ND1": 5, "LYS NZ": 11, "SER OG": 13},
  "hydrogens": 318,
  "acceptors": {"EGM O8x": 19, "HIS NE2": 5, "MET SD": 4, "SBM O3x": 23}
}
```

Check it before interpreting a comparison: an atom missing there cannot form
a hydrogen bond in the analysis.

## Figures

`polyzymd analyze hydrogen_bonds` writes these figures to
`<output-dir>/figures/hydrogen_bonds/`; `--no-plots` skips them.

| Figure | What it shows |
|---|---|
| `hbonds_<run>_comparison` | For a one-value result: each condition's mean with its interval and every replicate value |
| `hbonds_<run>_profile` | For `<s>_residues` and `<s>_pairs`: each entry's value per replicate and each condition's mean with its interval |
| `hbonds_<run>_difference` | For `<s>_residues` and `<s>_pairs` with several conditions: each condition minus the control at every residue or pair, with the interval of the difference and the significant entries marked |

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import HBOND_PARTS, hydrogen_bonds, residue_pair_hbond_occupancy

study = pz.Study.from_configs(
    {"SBMA 50%": "SBMA50/config.yaml", "SBMA 100%": "SBMA100/config.yaml"},
    equilibration="200ns",
)
rows = study.per_replicate(
    hydrogen_bonds,
    pz.select("chainid A"),
    pz.select("chainid C"),
    unit=None,
    parts=list(HBOND_PARTS),
    d_a_cutoff=3.0,
)
print(rows["mean_hbonds"].compare(control="SBMA 50%").to_agent_text())

pairs = study.per_replicate(
    residue_pair_hbond_occupancy,
    pz.select("chainid A"),
    pz.select("chainid C"),
    unit=None,
    labels="returned",
    missing=0.0,
    bounds=(0.0, 1.0),
)
print(pairs.compare(control="SBMA 50%").to_agent_text())
```

For the arguments and return values of each function, see
{doc}`../reference/analysis_functions`.

For a specific set of hydrogen bonds, such as the ones of a catalytic triad,
pass explicit atoms: {doc}`analysis_triad_quickstart` shows how.

## Next steps

- **What each function measures**: {doc}`../reference/analysis_functions`
- **How donors and acceptors are chosen**: {doc}`../explanation/analysis_hydrogen_bonds_verification`
- **How lifetimes are estimated**: {doc}`../explanation/analysis_contact_lifetimes`
- **Contacts analysis**: {doc}`analysis_contacts_quickstart`
