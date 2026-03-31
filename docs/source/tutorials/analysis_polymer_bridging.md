# Polymer Bridging Analysis (Experimental)

```{warning}
Polymer bridging is **experimental**. This plugin was contributed as a
proof-of-concept extensibility exercise. Metric definitions,
chemistry-aware profiling outputs, and interpretation guidance are all
subject to change. CLI output and generated figures carry explicit
experimental labels.
```

Quantify **per-fragment, per-frame multisite attachment** of individual polymer
chains (oligomers) to the enzyme surface directly from trajectories.

This analysis answers:
*"When a single oligomer chain contacts the protein, how often does it contact
more than one distinct protein residue — and which residue classes and monomer
types are involved?"*

## Core Concepts

### What Is an Observation?

The fundamental unit of data in this analysis is one **observation**: a single
polymer fragment in a single trajectory frame that makes at least one contact
with the protein (within the distance cutoff). Every observation records:

- Which protein residues are contacted
- The frame-wise C-alpha distances between those residues
- The amino acid class of each contacted residue
- The monomer identity of each contacting polymer residue
- The ordered residue-name sequence (signature) of the polymer fragment

One replicate may produce thousands of observations. All statistics are computed
over this observation population.

### What Is "Multisite"?

An observation is classified as **multisite** when the polymer fragment contacts
protein residues whose effective eligible valency exceeds 1. The meaning of
"eligible" depends on the `min_ca_distance_angstrom` setting:

- **`min_ca_distance_angstrom = 0`** (default): Any observation contacting 2+
  distinct residues counts as multisite, regardless of their spatial
  separation.

- **`min_ca_distance_angstrom > 0`** (e.g., 8.0 or 10.0): An observation
  is multisite only if at least two contacted residues have a frame-wise
  C-alpha separation >= the threshold. This filters out contacts with
  sequentially adjacent residues that happen to be geometrically close,
  focusing on **geometrically significant bridging** across the protein
  surface.

  Eligible valency is the number of residues that participate in at least one
  such qualifying pair. It is computed **per frame** from actual atomic
  coordinates, so the same pair of residues may qualify in one frame and not
  another if the protein is flexible.

### What Is "High Valency"?

An observation is **high valency** when its eligible valency is 3 or more —
meaning the polymer simultaneously contacts at least three spatially separated
protein residues in that frame.

## Primary Metrics

| Metric | Field name | Description |
|--------|-----------|-------------|
| **Mean Contacts / Oligomer** | `mean_contacts_per_contacting_oligomer` | Average number of distinct protein residues contacted per observation. |
| **Multisite Fraction** | `multisite_fraction` | Fraction of observations with eligible valency > 1. |
| **High-Valency Fraction** | `high_valency_fraction` | Fraction of observations with eligible valency >= 3. |
| **Valency Distribution** | `valency_probabilities` | Probability of 1-site, 2-site, and 3+-site attachment across all observations. |

These four metrics are used in the default cross-condition comparison pipeline
(t-tests, ANOVA, effect sizes, rankings).

## Chemistry-Aware Outputs (Experimental)

```{important}
All chemistry-aware outputs described below are labeled
`polymer_bridging_chemistry` in the experimental feature system. They are
**descriptive probabilities** — observed frequencies over the observation
population — not normalized enrichments and not evidence of mechanism.

Probabilities reflect what was observed in the simulation. They do not
control for surface accessibility, polymer composition, or reference
expectations. Interpret them as a **starting point for hypothesis
generation**, not as proof of preferential interaction.
```

Chemistry-aware outputs are computed only from **multivalent observations**
(eligible valency > 1). They are reported per-replicate and aggregated
(mean +/- SEM) across replicates.

### Protein Residue Classification

All protein residue outputs use the `ProteinAAClassification` scheme from
`polyzymd.analyses.shared.groupings`:

| Class | Amino Acids |
|-------|-------------|
| aromatic | PHE, TRP, TYR, HIS |
| charged_positive | ARG, LYS |
| charged_negative | ASP, GLU |
| polar | ASN, CYS, GLN, SER, THR |
| nonpolar | ALA, GLY, ILE, LEU, MET, PRO, VAL |
| unknown | Non-standard residues |

Common protonation-state variants (HIE, HID, HIP, ASH, GLH, etc.) are
automatically mapped to their parent residue.

### Anchor and Peripheral Residues

In each multivalent observation, the plugin identifies an **anchor** — the
protein residue with the closest atom-level distance to the polymer. All other
eligible contacted residues are **peripheral**.

- **Anchor protein class probabilities**: Frequency distribution of the amino
  acid class of the anchor residue across all multivalent observations.

- **Peripheral protein class probabilities**: Frequency distribution of the
  amino acid classes of non-anchor eligible residues.

- **Multivalent protein class probabilities**: Frequency distribution of the
  amino acid classes of *all* eligible residues in multivalent observations
  (anchor + peripheral combined).

### Polymer Monomer Probabilities

- **Polymer contact type probabilities**: Frequency of each polymer monomer
  type (by residue name, e.g. SBM, EGM) among all polymer residues that
  make protein contacts in multivalent observations.

- **Polymer anchor type probabilities**: Frequency of the polymer monomer
  type of the anchor (the polymer residue closest to the anchor protein
  residue) across multivalent observations.

### Cross-Classification Matrices

- **Anchor-to-peripheral class matrix**
  (`anchor_to_peripheral_group_matrix`): A row-normalized matrix where
  rows are the anchor protein class and columns are peripheral protein
  classes. Each row sums to 1.0. Answers: *"Given that the anchor is
  aromatic, what protein classes are the peripheral contacts?"*

- **Polymer-anchor to protein-anchor matrix**
  (`polymer_anchor_to_protein_anchor_matrix`): A row-normalized matrix
  where rows are polymer monomer types and columns are protein anchor
  classes. Each row sums to 1.0. Answers: *"Given that SBM is the polymer
  anchor monomer, which protein residue classes does it anchor to?"*

### Fragment Signature Probabilities

Each polymer fragment has an **ordered 5-mer signature** — the sequence of
residue names along the fragment (e.g., `EGM-EGM-SBM-EGM-EGM`). The top-10
most frequent signatures across multivalent observations are reported as
probabilities. These may help identify whether specific polymer subsequences
are over-represented in bridging events.

```{note}
Fragment signatures depend on the topology's residue ordering. In practice
the fragment length equals the number of monomers in the polymer chain. The
"5-mer" label is for illustration — the actual signature length is the full
fragment.
```

## Quick Start

### Step 1: Add bridging settings to your comparison YAML

```yaml
# comparison.yaml
plugins:
  polymer_bridging:
    cutoff: 4.5
    min_ca_distance_angstrom: 8.0   # Require contacted residues to be >= 8 A apart
    protein_selection: "protein"
    polymer_selection: "chainID C"  # Must match chain convention
```

Setting `min_ca_distance_angstrom: 0` disables the geometric filter entirely,
counting any 2+-residue observation as multisite.

### Step 2: Run the analysis

```bash
polyzymd compare run polymer_bridging -f comparison.yaml
```

The plugin automatically filters out conditions that have no polymer atoms
(e.g., a protein-only control).

### Step 3: Inspect results

The CLI prints a comparison table with the three primary metrics. Example:

```text
WARNING: Experimental analysis
Definitions and interpretation may change after the presentation release.
Affected: Polymer bridging chemistry profiling

Polymer Bridging Comparison
================================================================================
Multisite Fraction
  100% SBMA    : 0.312 +/- 0.018  (n=3)
  SBMA-EGPMA   : 0.487 +/- 0.025  (n=3)
  Pairwise: 100% SBMA -> SBMA-EGPMA  p=0.004 ** d=1.82  +56.1%  (more multisite)

Average Oligomer Valency
  100% SBMA    : 1.41 +/- 0.03  (n=3)
  SBMA-EGPMA   : 1.72 +/- 0.05  (n=3)
  Pairwise: 100% SBMA -> SBMA-EGPMA  p=0.008 ** d=1.54  +22.0%  (more bridging)

High-Valency Oligomers
  100% SBMA    : 0.051 +/- 0.009  (n=3)
  SBMA-EGPMA   : 0.128 +/- 0.014  (n=3)
  Pairwise: 100% SBMA -> SBMA-EGPMA  p=0.012 *  d=1.23  +150.8%  (more high-valency)
```

```{tip}
Use `--format json` to export the full `ComparisonResult` for downstream
analysis, or `--format markdown` for integration into reports.
```

### Step 4: Generate plots

```bash
polyzymd compare plot-all -f comparison.yaml
```

## Configuration Reference

### `plugins.polymer_bridging`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `cutoff` | float | `4.5` | Contact distance cutoff in Angstroms. Atom pairs within this distance are considered in contact. |
| `min_ca_distance_angstrom` | float | `0.0` | Minimum frame-wise C-alpha distance between contacted protein residues for an observation to count as multisite. Set > 0 to filter for geometrically significant bridging. Must be >= 0. |
| `protein_selection` | str | `"protein"` | MDAnalysis atom selection string for the protein. |
| `polymer_selection` | str | `"chainID C"` | MDAnalysis atom selection string for the polymer. Must match the chain convention (C = polymer). |

### Plot Settings

All plots are enabled by default and can be individually toggled in
`plot_settings.polymer_bridging`:

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `generate_multisite_bars` | bool | `true` | Bar chart of multisite fraction per condition. |
| `generate_mean_contacts_bars` | bool | `true` | Bar chart of mean contacted residues per oligomer. |
| `generate_valency_stack` | bool | `true` | Stacked bars showing 1 / 2 / 3+ valency distribution. |
| `generate_anchor_group_bars` | bool | `true` | Grouped bars of anchor protein residue class. |
| `generate_protein_group_stack` | bool | `true` | Stacked bars of protein classes in multivalent events. |
| `generate_anchor_peripheral_heatmap` | bool | `true` | Heatmap of anchor vs. peripheral protein class co-occurrence. |
| `generate_polymer_anchor_heatmap` | bool | `true` | Heatmap of polymer anchor monomer vs. protein anchor class. |
| `generate_fragment_signature_bars` | bool | `true` | Top-10 fragment signatures by frequency. |
| `figsize_bars` | (float, float) | `(9, 6)` | Figure size for bar charts. |
| `figsize_stack` | (float, float) | `(11, 6)` | Figure size for stacked charts. |
| `figsize_heatmap` | (float, float) | `(9, 7)` | Figure size for heatmaps. |

## Output Files

### Comparison cache

```text
comparison_output/
└── comparison/
    └── polymer_bridging/
        ├── <condition_label>/
        │   ├── rep1/
        │   │   └── result.json          # PolymerBridgingReplicateResult
        │   ├── rep2/
        │   │   └── result.json
        │   └── aggregated/
        │       └── result.json          # PolymerBridgingAggregatedResult
        └── result.json                  # ComparisonResult
```

### Generated Figures

All figures are stamped with an `EXPERIMENTAL` tag and saved to the configured
output directory (default: `figures/`).

| Filename | Content |
|----------|---------|
| `polymer_bridging_multisite_fraction.*` | Bar chart of multisite probability. |
| `polymer_bridging_mean_contacts.*` | Bar chart of average valency. |
| `polymer_bridging_valency_distribution.*` | Stacked bars: 1-site, 2-site, 3+ site. |
| `polymer_bridging_anchor_groups.*` | Grouped bars of anchor protein residue class. |
| `polymer_bridging_protein_group_distribution.*` | Stacked bars of all protein classes in multivalent events. |
| `polymer_bridging_anchor_peripheral_heatmap.*` | Anchor class (row) vs. peripheral class (column). |
| `polymer_bridging_polymer_anchor_heatmap.*` | Polymer monomer type (row) vs. protein anchor class (column). |
| `polymer_bridging_fragment_signatures.*` | Top-10 fragment signature probabilities. |

## CLI Reference

```text
polyzymd compare run polymer_bridging [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-f, --file PATH` | Path to `comparison.yaml` (default: `comparison.yaml`) |
| `--recompute` | Force recompute even if cached results exist |
| `--format [table\|markdown\|json]` | Output format (default: `table`) |
| `-o, --output PATH` | Save output to file |
| `-q, --quiet` | Suppress INFO messages |
| `--debug` | Enable DEBUG logging |
| `--eq-time TEXT` | Override equilibration time |

Alias: `polyzymd compare run bridging` (resolves to `polymer_bridging`).

## Loading Results Programmatically

```python
from polyzymd.analyses.polymer_bridging import (
    PolymerBridgingAggregatedResult,
    PolymerBridgingReplicateResult,
)

# Load per-replicate result
rep = PolymerBridgingReplicateResult.load(
    "comparison/polymer_bridging/my_condition/rep1/result.json"
)
print(f"Multisite fraction: {rep.multisite_fraction:.3f}")
print(f"Anchor protein groups: {rep.anchor_protein_group_probabilities}")

# Load aggregated result
agg = PolymerBridgingAggregatedResult.load(
    "comparison/polymer_bridging/my_condition/aggregated/result.json"
)
print(f"Mean valency: {agg.mean_contacts_per_contacting_oligomer:.2f} "
      f"+/- {agg.mean_contacts_sem:.2f}")

# Inspect cross-classification matrices
for anchor_class, peripherals in agg.anchor_to_peripheral_group_matrix_mean.items():
    for peripheral_class, prob in peripherals.items():
        if prob > 0.05:
            print(f"  {anchor_class} -> {peripheral_class}: {prob:.2f}")
```

## Interpretation Caveats

```{warning}
These caveats are essential for responsible use of this analysis. Read them
before presenting or publishing polymer bridging results.
```

1. **Descriptive, not mechanistic.** All outputs are observed frequencies. A
   high anchor probability for aromatic residues does not prove that aromatic
   anchoring drives conjugate stability — it may reflect surface composition,
   polymer placement, or simulation artifacts.

2. **Not normalized enrichment.** Unlike binding preference analysis, polymer
   bridging probabilities are **raw frequencies over the observation
   population**, not enrichments normalized by surface availability. A protein
   with 30% aromatic surface and 30% aromatic anchors is not showing
   enrichment — it is showing baseline proportionality. Compare with surface
   composition before drawing conclusions.

3. **Frame-wise, not event-wise.** Each trajectory frame generates independent
   observations. A polymer that bridges two residues for 500 consecutive frames
   counts as 500 observations, not one sustained bridging event. Residence-time
   analysis is not yet implemented.

4. **C-alpha distance is dynamic.** When `min_ca_distance_angstrom > 0`, the
   threshold is evaluated per frame against actual C-alpha coordinates. Protein
   breathing motions mean that the same pair of residues may qualify in some
   frames and not others. This is physically correct but can make results
   sensitive to protein flexibility.

5. **Anchor selection is heuristic.** The anchor is the polymer-protein pair
   with the minimum atom-level distance. In cases with multiple equidistant
   contacts, the choice is arbitrary. This affects anchor-specific outputs but
   not the primary multisite/valency metrics.

6. **Fragment signatures assume topology ordering.** The ordered monomer
   sequence comes directly from the topology. If the topology's residue
   ordering does not reflect the true polymer sequence, signatures will be
   misleading.

7. **Conditions without polymer are filtered.** The plugin automatically
   excludes conditions that lack polymer atoms (e.g., protein-only controls).
   This is correct behavior but means the control condition for pairwise
   statistics must itself contain polymer.

8. **Proof-of-concept status.** This plugin was contributed as an extensibility
   exercise. The analysis methodology has not been independently validated.
   Use it for internal hypothesis generation, not for publication-ready
   claims, until the methodology matures.

## Relation to Other Analyses

| Analysis | What It Measures | Relation to Polymer Bridging |
|----------|-----------------|------------------------------|
| **Contacts** | Total contact counts and frequencies | Polymer bridging decomposes contacts per-chain, adding valency information. |
| **Binding Preference** | Enrichment by residue class | Provides surface-normalized context that bridging lacks. |
| **Polymer Affinity** | Total interaction strength (N x deltaG) | Complementary: affinity measures total adhesion; bridging measures spatial distribution of adhesion per chain. |
| **RMSF** | Structural flexibility | Complementary: does multisite bridging correlate with reduced flexibility? |
| **Catalytic Triad** | Active site geometry | Complementary: do bridging events coincide with triad perturbation? |

## Troubleshooting

### "No conditions passed polymer filtering"

All conditions were filtered out because the plugin could not detect polymer
atoms. Check that your `polymer_selection` matches your topology (default is
`"chainID C"`, following the PolyzyMD chain convention).

### Very low multisite fraction

If multisite fraction is near zero:
- Check that `min_ca_distance_angstrom` is not too stringent. A value of 20+ A
  may filter out nearly all events for small proteins.
- Verify that the polymer fragments contain more than one monomer.
- Check that the contact cutoff is appropriate for your force field (4.5 A is
  standard for heavy-atom contacts).

### All probabilities are empty (`{}`)

Chemistry-aware outputs require multivalent observations. If no observation has
eligible valency > 1, all chemistry dictionaries will be empty. Lower
`min_ca_distance_angstrom` or verify that the polymer makes multi-residue
contacts.

### Heatmaps show only one condition

The anchor-peripheral and polymer-anchor heatmaps currently display data for
the **first condition only** (by label order). This is a known limitation of
the current plotting code. To compare matrices across conditions, load the
aggregated results programmatically (see above).

## See Also

- [Contacts Analysis Quick Start](analysis_contacts_quickstart.md) — prerequisite contact computation
- [Binding Preference Analysis](analysis_binding_preference.md) — surface-normalized enrichment (complementary)
- [Polymer Affinity Analysis](analysis_polymer_affinity.md) — total interaction strength scoring
- [Statistics Best Practices](analysis_statistics_best_practices.md) — replicate planning
- [Comparing Conditions](analysis_compare_conditions.md) — multi-condition workflows
- [Extending the Analysis Framework](extending_analyses.md) — contribute a new plugin
