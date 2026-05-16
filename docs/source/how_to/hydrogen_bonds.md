# How-To: Analyze Hydrogen Bonds Across Conditions

This guide shows you how to configure, run, and interpret the `hydrogen_bonds`
analysis plugin. By the end, you will have a working hydrogen-bond comparison
across simulation conditions, with grouped summaries, composition breakdowns,
and publication-ready figures.

:::{admonition} Environment Setup
:class: tip

All commands below assume you have activated the PolyzyMD pixi environment:

```bash
pixi shell -e build
```

Alternatively, prefix each command with `pixi run -e build`.
:::

## When to Use This Plugin

Use `hydrogen_bonds` when you need to answer questions like:

- Does the polymer form more H-bonds with the protein than the substrate does?
- Are protein internal H-bonds disrupted by polymer conjugation?
- Which specific residue pairs form the most persistent H-bonds with the
  polymer?
- How is the total H-bond budget distributed across protein, substrate, and
  polymer?

The plugin wraps MDAnalysis `HydrogenBondAnalysis` with a flexible
**groups + summaries** model: you define named atom groups, then specify
which pairs or self-interactions to summarize. This lets you answer multiple
questions in a single analysis run.

## Before You Start

You need:

- A working pixi environment (`pixi install -e build`)
- Completed simulation trajectories for at least two conditions
- A `comparison.yaml` defining your conditions (see
  {doc}`analysis_compare_conditions`)
- Familiarity with the PolyzyMD chain convention (A=protein, B=substrate,
  C=polymer)

## Quick Start

Add a `hydrogen_bonds` section to the `plugins:` block in your
`comparison.yaml`:

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

Run the analysis:

```bash
pixi run -e build polyzymd compare run hydrogen_bonds
```

That is all you need for a minimal protein–polymer H-bond comparison. The
sections below explain every setting and show more advanced configurations.

## Configuration Reference

All settings live under `plugins: hydrogen_bonds:` in `comparison.yaml`. Every
field has a sensible default, so you only need to specify what you want to
change.

### Geometric Criteria

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `distance_cutoff` | float | `3.0` | Donor–acceptor distance cutoff in Å |
| `angle_cutoff` | float | `150` | D–H···A angle cutoff in degrees |

These match the standard MDAnalysis `HydrogenBondAnalysis` defaults. The
distance is measured between the heavy-atom donor and the acceptor (not the
hydrogen). The angle is measured at the hydrogen: D–H···A.

:::{admonition} Topology and hydrogen caveats
:class: important

The plugin passes your configured group union to MDAnalysis as the donor and
acceptor selection, and uses `(<group union>) and element H` for hydrogens. This
means your topology must contain explicit hydrogen atoms and usable element
metadata for hydrogen selection. Topologies without hydrogens, with missing or
incorrect element records, or with coarse-grained beads will undercount or
return no H-bonds. The current public settings do not expose separate donor,
hydrogen, and acceptor selection strings.
:::

### Groups

```yaml
groups:
  protein_all: "protein"
  substrate: "resname pNB"
  polymer_all: "chainid C"
```

Each entry maps a **name** (used in summaries) to an **MDAnalysis selection
string**. Groups are resolved once at the start of the analysis. You can use
any valid MDAnalysis selection syntax including `protein`, `chainid`,
`resname`, `resid`, boolean operators, and parentheses.

:::{admonition} Group overlaps
:class: warning

If two groups share atoms (e.g., `"protein"` and
`"protein and resid 106"`), the plugin logs a warning and H-bonds in the
overlap region may be counted in multiple summaries. For clean composition
accounting, use disjoint groups.
:::

### Summaries

Summaries define **what to measure**. Each summary has a unique name and
exactly one mode:

- **`between: [group_a, group_b]`** — Count H-bonds where one partner is in
  group A and the other in group B (inter-group).
- **`within: group_name`** — Count H-bonds where both partners are in the
  same group (intra-group).

```yaml
summaries:
  protein_polymer:
    between: [protein_all, polymer_all]
  protein_internal:
    within: protein_all
```

Each summary produces an independent metric (`mean_hbonds_<name>`) in the
comparison output. You can define as many summaries as you need — the plugin
runs one global MDAnalysis H-bond search and post-classifies events into
summaries.

### Composition

Composition analysis breaks down the **total** H-bond budget by
donor→acceptor partition pairs. This is independent of summaries and uses
its own set of disjoint selections.

```yaml
composition:
  partitions:
    protein: "protein"
    substrate: "resname pNB"
    polymer: "chainid C"
```

The plugin counts every detected H-bond event, classifies donor and acceptor
atoms into partitions, and reports both absolute counts (mean H-bonds/frame)
and fractional shares per partition pair (e.g., protein→polymer,
protein→protein, polymer→substrate).

Composition is opt-in. If `composition` is omitted, no composition results are
computed.

:::{admonition} Partition disjointness
:class: important

Composition partitions should be disjoint (no shared atoms). If partitions
overlap, the plugin raises an error by default. Set
`allow_overlapping_composition: true` to allow overlap with warnings
(overlapping atoms are counted in both partitions, so fractions may exceed 1.0).
:::

### Other Settings

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `update_selections` | bool | `true` | Re-evaluate simple atom-group selections each frame. Coordinate-dependent group selections such as `around`, `sphzone`, and `cyzone` are not supported when this is `true`; use static groups or separate preprocessed groups instead. |
| `top_n_pairs` | int | `15` | Number of top residue pairs to report per summary |
| `allow_empty_groups` | bool | `true` | If `true` (default), warn and skip summaries when a referenced group selection matches no atoms. Set `false` for strict validation (raise `ValueError`). |
| `allow_overlapping_composition` | bool | `false` | If `false`, overlapping composition partitions raise an error. Set `true` to allow overlap with warnings. |
| `timestep_ps` | float or null | `null` | Manual frame spacing in ps for time-axis plots. If null, read from trajectory metadata. |

If you want composition partitions to mirror group selections, define them
explicitly in YAML with the same keys/selections.

:::{admonition} Technical Detail
:class: note

`update_selections: true` makes MDAnalysis re-evaluate simple atom-group
selections every frame. This is appropriate for structural selections such as
`protein`, `chainid`, `resname`, and `resid`.

Coordinate-dependent group selections such as `around`, `sphzone`, and `cyzone`
are rejected by source validation when `update_selections: true` because the
plugin resolves group membership for donor/acceptor post-classification from
initial frame index sets. Use static groups, or create separate preprocessed
groups before running the analysis.

However, the plugin resolves group membership used for donor/acceptor
pair-classification once at the start from the initial frame index sets. This
means pair-tracking classification does not change if atoms move between groups
later in the trajectory. This behavior is deliberate for performance and for
consistent summary definitions across frames.
:::

:::{admonition} When to set update_selections to false
:class: note

If all your groups use structural selections (`chainid`, `resname`, `resid`,
`protein`), setting `update_selections: false` is safe and faster. Keep
coordinate-dependent selections out of `groups` when `update_selections: true`;
use static groups or separate preprocessed groups for distance-based
membership.
:::

## Common Recipes

### Protein ↔ polymer H-bonds (minimal)

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

### Standard 4-summary configuration with composition

This is the recommended starting point for a typical enzyme–polymer–substrate
system. It covers the four most common interactions and adds a composition
breakdown.

```yaml
plugins:
  hydrogen_bonds:
    distance_cutoff: 3.0
    angle_cutoff: 150
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
    composition:
      partitions:
        protein: "protein"
        substrate: "resname pNB"
        polymer: "chainid C"
```

### Monomer-resolved polymer analysis

If your polymer contains distinct monomer types (e.g., SBMA and EGMA), you can
define separate groups to compare their H-bonding behavior:

```yaml
plugins:
  hydrogen_bonds:
    groups:
      protein_all: "protein"
      substrate: "resname pNB"
      polymer_all: "chainid C"
      polymer_sbm: "chainid C and resname SBM"
      polymer_egm: "chainid C and resname EGM"
      active_site: "protein and (resid 106 or resid 225 or resid 188)"
    summaries:
      protein_polymer:
        between: [protein_all, polymer_all]
      protein_sbm:
        between: [protein_all, polymer_sbm]
      protein_egm:
        between: [protein_all, polymer_egm]
      active_site_polymer:
        between: [active_site, polymer_all]
      protein_internal:
        within: protein_all
      polymer_internal:
        within: polymer_all
    composition:
      partitions:
        protein: "protein"
        substrate: "resname pNB"
        polymer_sbm: "chainid C and resname SBM"
        polymer_egm: "chainid C and resname EGM"
```

### Protein ↔ substrate only

For systems without a polymer (or for a control condition comparison):

```yaml
plugins:
  hydrogen_bonds:
    groups:
      protein_all: "protein"
      substrate: "resname pNB"
    summaries:
      protein_substrate:
        between: [protein_all, substrate]
```

### Active-site focused analysis

To focus on H-bonds near the catalytic site:

```yaml
plugins:
  hydrogen_bonds:
    groups:
      active_site: "protein and (resid 106 or resid 225 or resid 188)"
      polymer_all: "chainid C"
      substrate: "resname pNB"
    summaries:
      active_site_polymer:
        between: [active_site, polymer_all]
      active_site_substrate:
        between: [active_site, substrate]
```

## Running the Analysis

### Local execution

```bash
pixi run -e build polyzymd compare run hydrogen_bonds
```

This runs the full pipeline: the per-replicate compute stage for every
replicate, `aggregate` for every condition, then `compare` and `plot`.

### HPC execution

For long trajectories or many replicates, submit as a SLURM job:

```bash
pixi run -e build polyzymd compare submit hydrogen_bonds \
    -f comparison.yaml \
    --partition aa100 \
    --mem 8G \
    --time 02:00:00
```

:::{admonition} Performance note
:class: note

The hydrogen bonds plugin is marked with `execution_cost_hint = "high"`.
H-bond detection iterates over donor–acceptor pairs across all selected atoms
every frame. For large systems (50k+ atoms) with many frames, allocate at
least 8 GB memory and 1–2 hours per replicate.
:::

See {doc}`hpc_execution` for the full HPC submission guide.

## Understanding the Output

### Per-summary metrics

Each summary produces these fields in the aggregated result:

| Field | Description |
|-------|-------------|
| `mean_hbonds_per_frame` | Average number of H-bonds per frame, averaged across replicates |
| `sem_hbonds_per_frame` | Standard error of the mean across replicates |
| `per_replicate_mean_hbonds` | Per-replicate mean values (for statistical tests) |
| `mean_unique_pairs_per_frame` | Mean number of unique donor-residue/acceptor-residue pairs per frame |
| `sem_unique_pairs_per_frame` | Standard error across replicates for the unique-pairs metric |
| `mean_fraction_with_any` | Fraction of frames that had at least one H-bond |

The primary comparison metric is `mean_hbonds_per_frame` per summary. In the
comparison output, this appears as `mean_hbonds_<summary_name>` (e.g.,
`mean_hbonds_protein_polymer`).

### Residue pairs

Each summary reports the top residue pairs ranked by **occupancy** (fraction
of frames where the pair has at least one H-bond). Two views are provided:

- **Directed pairs** — Donor→Acceptor direction is preserved. Useful for
  understanding H-bond directionality (which residue donates).
- **Undirected pairs** — Both directions are merged into one pair. Useful for
  identifying which residue *pairs* interact most frequently regardless of
  direction.

Each pair reports `mean_occupancy`, `sem_occupancy`, `mean_events_per_frame`,
and per-replicate occupancy values.

### Composition entries

When composition is configured, the result includes per-partition-pair entries:

| Field | Description |
|-------|-------------|
| `donor_partition → acceptor_partition` | The partition pair (e.g., protein→polymer) |
| `mean_hbonds_per_frame` | Absolute count for this pair |
| `mean_fraction_of_total` | This pair's share of all detected H-bonds |

Composition fractions help normalize for differences in donor/acceptor
availability across conditions. If a polymer condition has more total H-bonds
but the protein→protein fraction stays constant, the polymer is adding new
interactions rather than disrupting existing ones.

### Statistical comparison

The default comparison pipeline produces:

- **Pairwise t-tests** between all condition pairs for each summary metric,
  with Benjamini–Hochberg FDR correction
- **Cohen's d** effect sizes for each pairwise comparison
- **ANOVA** (when ≥3 conditions) as an omnibus test
- **Rankings** of conditions by each metric

The `higher_is_better` flag is set to `None` for H-bond counts (neither
direction is universally better), and direction labels use "fewer H-bonds" /
"similar" / "more H-bonds".

## Interpreting the Plots

The plugin generates up to five plot types. All plots use the shared PolyzyMD
theming system for consistent publication-quality output.

### Summary comparison bar chart

**File:** `hbond_summary_comparison.png`

Faceted grouped bars with one subplot per summary and one bar per condition.
Each summary gets its own y-axis scale, improving readability when summaries
have very different magnitudes. Error bars show SEM across replicates.

### Time series

**File:** `hbond_timeseries_<summary>.png` (one per summary)

Per-frame H-bond counts over time in physical units (ns). The mean across
replicates is shown as a solid line with a ±1 SD shaded band (single-replicate
runs show the trace only). Useful for detecting equilibration issues (early drift) or transient events
(sudden spikes indicating a binding/unbinding event).

### Top residue pairs occupancy

**File:** `hbond_top_pairs_<summary>.png` (one per summary)

Horizontal grouped bar chart showing top undirected residue pairs by mean
occupancy. To focus on cross-condition comparison, only pairs present in at
least two conditions are shown.

### Composition absolute bars

**File:** `hbond_composition_absolute.png`

Stacked bar chart showing mean H-bonds/frame broken down by partition pair
(e.g., protein→protein, protein→polymer, polymer→polymer). Each condition
gets one stacked bar. Reveals how the total H-bond budget is distributed.

### Composition fractional bars

**File:** `hbond_composition_fraction.png`

Same as above but normalized to fractions. The y-axis auto-scales: for disjoint
partitions it ranges 0–1, but when `allow_overlapping_composition: true` is set
and partitions overlap, stacked fractions may exceed 1.0. Useful for comparing
the *relative* H-bond distribution even when total counts differ between
conditions.

## Scientific Considerations

:::{admonition} Interpreting H-bond counts
:class: warning

Keep these caveats in mind when interpreting results:

1. **Geometric criteria, not energetics.** H-bond counts reflect MDAnalysis
   geometric criteria (distance and angle cutoffs), not interaction energies.
   A "hydrogen bond" by this definition may not correspond to a thermodynamically
   significant interaction.

2. **System size matters.** Cross-condition comparisons of raw H-bond counts
   should account for differences in system size (number of donors/acceptors).
   The composition fractional view helps normalize for this.

3. **Composition fractions normalize availability.** If condition A has 200
   polymer atoms and condition B has 400, condition B will likely show more
   protein–polymer H-bonds in absolute terms. Composition fractions reveal
   whether the *proportion* also changes.

4. **Residue-pair rankings are noisy.** Individual residue-pair occupancies
   can vary substantially across replicates. Focus on the top 3–5 pairs and
   cross-reference with structural knowledge. The `top_n_pairs` setting
   controls how many pairs are reported (default: 15).

5. **Statistical correction.** Pairwise comparisons use FDR-corrected
   (Benjamini–Hochberg) p-values. When comparing many summaries across many
   conditions, check adjusted p-values rather than raw values.
:::

## Full Example: comparison.yaml

Here is a complete `comparison.yaml` for a CALB enzyme study with three
polymer conditions:

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
    composition:
      partitions:
        protein: "protein"
        substrate: "resname pNB"
        polymer: "chainid C"

plot_settings:
  format: "png"
  dpi: 300
  style: "publication"
```

## Aliases

The plugin responds to these CLI names:

- `hydrogen_bonds` (canonical)
- `hbonds`
- `hbond`

```bash
# These are equivalent:
polyzymd compare run hydrogen_bonds
polyzymd compare run hbonds
polyzymd compare run hbond
```

## See Also

- {doc}`analysis_compare_conditions` — Setting up `comparison.yaml`
- {doc}`../reference/analysis_comparison_reference` — Plugin listing and
  statistical terms
- {doc}`hpc_execution` — Submitting analysis jobs to SLURM
- {doc}`../explanation/analysis_statistics_best_practices` — Autocorrelation and uncertainty
- {doc}`../contributor_guide/extending_analyses` — Writing your own analysis plugin
