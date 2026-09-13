# How-To: Analyze Hydrogen Bonds Across Conditions

This guide shows how to configure and run the `hydrogen_bonds` plugin and how
to read what it writes. For the full settings table and the artifact layout, see
{doc}`../reference/analysis_hydrogen_bonds_reference`.

:::{admonition} Environment Setup
:class: tip

All commands below assume the PolyzyMD analysis pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## When to use this plugin

Use `hydrogen_bonds` to answer questions like:

- Does the polymer form more hydrogen bonds with the protein than the substrate
  does?
- Are protein internal hydrogen bonds disrupted by polymer conjugation?
- Which residue pairs hold the most persistent hydrogen bonds with the polymer?

The plugin wraps MDAnalysis `HydrogenBondAnalysis` in a groups and summaries
model: name the atom groups, then say which pairs or self-interactions to
report. One detection pass answers every summary.

## Before you start

You need completed trajectories for at least two conditions, a
`comparison.yaml` defining them (see {doc}`analysis_compare_conditions`), and a
topology with explicit hydrogens and element metadata. A coarse-grained or
hydrogen-free topology cannot be analyzed.

## Quick start

Add a `hydrogen_bonds` section to the `plugins:` block:

```yaml
plugins:
  hydrogen_bonds:
    groups:
      protein_all: "protein"
      polymer_all: "chainid C"
    summaries:
      protein_polymer:
        between: [protein_all, polymer_all]
```

Then run it:

```bash
pixi run -e analysis polyzymd compare run hydrogen_bonds
```

## Choose the geometric criteria

`distance_cutoff` (3.0 angstrom by default) is measured between the heavy-atom
donor and the acceptor, not the hydrogen. `angle_cutoff` (150 degrees by
default) is measured at the hydrogen, D-H...A. Both defaults follow Smith et al.
2019. Counts depend on the cutoffs, so compare only conditions analyzed with the
same pair.

## Choose the donor and acceptor elements

```yaml
plugins:
  hydrogen_bonds:
    donor_acceptor_elements: ["N", "O", "S"]
```

The default is `["N", "O"]`. Carbon is excluded because a C-H donor or a carbon
acceptor is not a hydrogen bond under the IUPAC definition, and including carbon
inflates counts by an amount that depends on the polymer chemistry. Add `"S"`
when cysteine and methionine sulfur should donate and accept. Hydrogen is
rejected, because hydrogens are chosen separately and listing `"H"` here would
make every one of them a donor and an acceptor. Symbols are capitalized for you
and an unknown one is rejected when the config is read. Changing this changes
the settings fingerprint, so cached replicates are recomputed.

## Name the groups

```yaml
groups:
  protein_all: "protein"
  substrate: "resname pNB"
  polymer_all: "chainid C"
```

Each entry maps a name used in summaries to an MDAnalysis selection string. A
group that matches no atoms raises `SelectionError`; set
`allow_empty_groups: true` to report its summaries as zero instead, which is
what a control condition without polymer needs.

Groups may overlap. A bond in the overlap is counted by every summary whose
filter it passes, so overlapping groups double count on purpose. The plugin
warns when it sees an overlap, naming the two groups and how many atoms they
share.

Group membership is resolved once at the start of the window and then held
fixed, so keep coordinate-dependent selections such as `around` and `sphzone`
out of `groups` unless you want membership frozen at the first analyzed frame.
`update_selections` re-evaluates only MDAnalysis's own donor, hydrogen and
acceptor selections.

## Define the summaries

Each summary has a unique name and exactly one mode. `between: [a, b]` counts
bonds with one partner in group A and the other in group B, in either
direction. `within: a` counts bonds with both partners in group A. Bonds inside
one residue are always dropped.

```yaml
summaries:
  protein_polymer:
    between: [protein_all, polymer_all]
  protein_internal:
    within: protein_all
```

## Common recipes

### Protein and polymer only

```yaml
plugins:
  hydrogen_bonds:
    groups:
      protein_all: "protein"
      polymer_all: "chainid C"
    summaries:
      protein_polymer:
        between: [protein_all, polymer_all]
```

### The four standard partitions

The usual starting point for an enzyme, polymer and substrate system:

```yaml
plugins:
  hydrogen_bonds:
    distance_cutoff: 3.0
    angle_cutoff: 150
    allow_empty_groups: true
    groups:
      protein_all: "protein"
      substrate: "resname pNB"
      polymer_all: "chainid C"
    summaries:
      protein_polymer:
        between: [protein_all, polymer_all]
      protein_substrate:
        between: [protein_all, substrate]
      protein_internal:
        within: protein_all
      polymer_internal:
        within: polymer_all
```

### Monomer-resolved polymer

Separate groups per monomer type compare their hydrogen bonding:

```yaml
plugins:
  hydrogen_bonds:
    groups:
      protein_all: "protein"
      polymer_all: "chainid C"
      polymer_sbm: "chainid C and resname SBM"
      polymer_egm: "chainid C and resname EGM"
    summaries:
      protein_polymer:
        between: [protein_all, polymer_all]
      protein_sbm:
        between: [protein_all, polymer_sbm]
      protein_egm:
        between: [protein_all, polymer_egm]
```

### Active site only

```yaml
plugins:
  hydrogen_bonds:
    groups:
      active_site: "protein and resid 106 225 188"
      polymer_all: "chainid C"
      substrate: "resname pNB"
    summaries:
      active_site_polymer:
        between: [active_site, polymer_all]
      active_site_substrate:
        between: [active_site, substrate]
```

### Catalytic triad geometry

A group per triad atom set reports whether the catalytic hydrogen bonds hold:

```yaml
plugins:
  hydrogen_bonds:
    groups:
      cat_ser: "protein and resid 76 and (name OG or element H)"
      cat_his: "protein and resid 155 and (name ND1 NE2 or element H)"
      cat_asp: "protein and resid 132 and name OD1 OD2"
    summaries:
      ser_his:
        between: [cat_ser, cat_his]
      asp_his:
        between: [cat_asp, cat_his]
```

## Run the analysis

Locally:

```bash
pixi run -e analysis polyzymd compare run hydrogen_bonds
```

On SLURM, for long trajectories or many replicates:

```bash
pixi run -e analysis polyzymd compare submit hydrogen_bonds \
    -f comparison.yaml \
    --partition aa100 \
    --mem 8G \
    --time 02:00:00
```

The plugin declares `execution_cost_hint = "high"` and a default resource hint
of 8 GB and two hours. Detection iterates donor and acceptor pairs every frame,
so a system of 50k atoms and a few thousand frames needs about that. See
{doc}`hpc_execution` for the full submission guide.

## Read the output

Each summary reports two observables. `hbonds_<summary>` is a
`mean_of_timeseries` in counts, one value per frame, and it is the quantity the
cross-condition tests use. `pair_occupancy_<summary>` is a `profile` of the
`top_n_pairs` most persistent residue pairs, indexed by rank, with the pair
named in the observable's `metadata["pair_labels"]` as `TYR138(A)-SBM152(C)`.
Ranking inside a replicate biases the low ranks upward, so read the profile as
the shape of the occupancy spectrum rather than as an estimate for one pair. If
a particular pair matters, give it its own summary and the framework will test
it properly.

Every mean, standard error and interval in the condition and comparison
artifacts is computed across replicates, never across frames. The per-frame
series survives in `observables.npz` beside the replicate result.

The raw detection output is kept in
`sidecars/hydrogen_bond_events.npz`: one row per bond per frame, with the frame
index, the donor, hydrogen and acceptor atom indices, the distance and the
angle. Load it to ask questions the summaries do not cover, such as which donor
partitions carry the total budget.

## Scientific considerations

:::{admonition} Interpreting hydrogen-bond counts
:class: warning

1. Geometric criteria, not energetics. A bond by this definition need not be a
   thermodynamically significant interaction.

2. System size matters. A condition with twice the polymer has more donors and
   acceptors available, so compare per-frame counts only between systems whose
   composition you have accounted for.

3. Residue-pair rankings are noisy. Individual pair occupancies vary between
   replicates, which is why the profile is indexed by rank and reported with a
   per-rank standard error across replicates. Read the top few ranks and check
   them against what you know about the structure.

4. Statistical correction. Pairwise tests are Benjamini-Hochberg corrected over
   the whole run. Read the adjusted p-value, not the raw one.
:::

## Full example

```yaml
name: "calb_hbond_study"
description: "Hydrogen bond analysis for CALB with polymer conjugates"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../noPoly_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

  - label: "SBMA-100"
    config: "../SBMA_100_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

  - label: "EGMA-100"
    config: "../EGMA_100_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  hydrogen_bonds:
    distance_cutoff: 3.0
    angle_cutoff: 150
    allow_empty_groups: true
    groups:
      protein_all: "protein"
      substrate: "resname pNB"
      polymer_all: "chainid C"
    summaries:
      protein_polymer:
        between: [protein_all, polymer_all]
      protein_substrate:
        between: [protein_all, substrate]
      protein_internal:
        within: protein_all
      polymer_internal:
        within: polymer_all
```

## See also

- {doc}`../reference/analysis_hydrogen_bonds_reference` (settings and outputs)
- {doc}`analysis_compare_conditions` (setting up `comparison.yaml`)
- {doc}`hpc_execution` (submitting analysis jobs to SLURM)
- {doc}`../explanation/analysis_statistics_best_practices` (autocorrelation and
  uncertainty)
