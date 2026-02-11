# Catalytic Triad Analysis: Statistical Best Practices

A comprehensive guide to analyzing catalytic triad geometry from MD simulations,
including statistical considerations, interpretation guidelines, and worked
examples from real enzyme-polymer systems.

```{note}
**Just need quick results?** See the [Quick Start Guide](analysis_triad_quickstart.md)
for copy-paste commands and minimal setup.
```

## Introduction

The catalytic triad is a fundamental structural motif in many enzyme families,
including serine proteases, lipases, and esterases. Maintaining proper triad
geometry is essential for catalytic activity. This guide explains:

- What catalytic triads are and why they matter
- How to measure triad integrity in MD simulations
- The "simultaneous contact fraction" metric
- Statistical considerations for robust analysis
- How to interpret results in the context of enzyme function

## What is a Catalytic Triad?

### The Ser-His-Asp/Glu Mechanism

The classic catalytic triad consists of three amino acids working together:

1. **Serine (Ser)** - The nucleophile that attacks the substrate
2. **Histidine (His)** - The general base/acid that shuttles protons
3. **Aspartate (Asp) or Glutamate (Glu)** - Orients and activates histidine

These residues form a hydrogen bonding network:

```
        O                         O
        ||                        ||
Substrate-C     ->    Substrate-C-O-Ser
        |                         |
       ...                       ...

    Asp/Glu --- His --- Ser --- Substrate
       |         |       |
      C=O      N-H     O-H
       |         |       |
      O-H .... N      [nucleophilic attack]
```

The proton relay works as follows:
1. Asp stabilizes the positively charged His through a hydrogen bond
2. His abstracts a proton from Ser's hydroxyl group
3. The activated Ser oxygen attacks the substrate carbonyl

### Critical Distances

For the triad to function, specific hydrogen bonds must be maintained:

| Pair | Donor | Acceptor | Typical Distance |
|------|-------|----------|------------------|
| Asp-His | His ND1 (N-H) | Asp OD1/OD2 | 2.7 - 3.5 A |
| His-Ser | Ser OG (O-H) | His NE2 | 2.7 - 3.5 A |

These distances are measured between heavy atoms (N, O). The actual H-bond
length (H...acceptor) is ~1.0 A shorter.

### Why Geometry Matters

If these distances increase beyond ~4 A:
- Hydrogen bonds break
- Proton relay cannot function
- Catalytic activity is lost or severely reduced

MD simulations can reveal whether polymer conjugation, mutations, or other
modifications disrupt triad geometry.

## The Simultaneous Contact Metric

### Definition

The **simultaneous contact fraction** is the percentage of frames where
**all** distance pairs in the triad are below the contact threshold at the
same time:

$$
f_{\text{contact}} = \frac{1}{N} \sum_{t=1}^{N} \prod_{i=1}^{M} \mathbb{1}[d_i(t) < \theta]
$$

Where:
- $N$ = number of frames
- $M$ = number of pairs (typically 2 for a triad)
- $d_i(t)$ = distance of pair $i$ at frame $t$
- $\theta$ = contact threshold (typically 3.5 A)
- $\mathbb{1}[\cdot]$ = indicator function (1 if true, 0 if false)

### Why Simultaneous Matters

Consider two scenarios with 50% individual contact per pair:

**Scenario A: Alternating**
```
Frame:    1  2  3  4  5  6  7  8
Asp-His:  +  -  +  -  +  -  +  -   (50% contact)
His-Ser:  -  +  -  +  -  +  -  +   (50% contact)
Both:     -  -  -  -  -  -  -  -   (0% simultaneous!)
```

**Scenario B: Correlated**
```
Frame:    1  2  3  4  5  6  7  8
Asp-His:  +  +  +  +  -  -  -  -   (50% contact)
His-Ser:  +  +  +  +  -  -  -  -   (50% contact)
Both:     +  +  +  +  -  -  -  -   (50% simultaneous)
```

Scenario A has the same per-pair statistics but **zero** catalytic competence
because the triad is never intact. The simultaneous contact fraction captures
this critical distinction.

### Interpretation Guidelines

| Simultaneous Contact | Interpretation |
|---------------------|----------------|
| > 80% | Excellent triad integrity |
| 50 - 80% | Good integrity, some flexibility |
| 20 - 50% | Moderate disruption |
| 5 - 20% | Significant disruption |
| < 5% | Severely disrupted triad |

```{warning}
These thresholds are guidelines, not absolute rules. The relationship between
simultaneous contact and actual catalytic activity depends on the enzyme and
should be validated against experimental data when possible.
```

## H-Bond Distance Thresholds

### Choosing a Threshold

The default threshold of **3.5 A** is based on typical H-bond geometry:

| Threshold | Rationale |
|-----------|-----------|
| 3.0 A | Strict - strong H-bonds only |
| 3.5 A | Standard - typical H-bond cutoff |
| 4.0 A | Relaxed - includes weaker interactions |

### When to Adjust

**Use a stricter threshold (3.0 A)** when:
- Comparing to crystal structures (often show shorter distances)
- Studying high-activity conformations
- The default gives near-100% contact (not discriminating)

**Use a relaxed threshold (4.0 A)** when:
- Studying dynamic enzymes with flexible active sites
- Default gives very low contact but enzyme is known to be active
- Accounting for force field limitations

### Heavy Atom vs Hydrogen Distances

PolyzyMD measures **heavy atom distances** (N-O, O-O), not hydrogen positions:

| Measurement | Typical Value |
|-------------|---------------|
| Heavy atom (N...O) | 2.7 - 3.5 A |
| H-bond (H...O) | 1.7 - 2.5 A |

This is more robust because hydrogen positions are often uncertain in MD
(fast motion, force field limitations).

## Autocorrelation Analysis

### The Correlation Problem

Contact states in MD are temporally correlated. If the triad is in contact at
frame 100, it's very likely to be in contact at frame 101. This affects
uncertainty estimation.

### How PolyzyMD Handles This

PolyzyMD computes the autocorrelation function of the simultaneous contact
timeseries and estimates:

1. **Correlation time (tau)** - characteristic decay time
2. **N_independent** - effective number of independent samples
3. **Corrected SEM** - uncertainty accounting for correlation

### SEM for Proportions

For a binary contact timeseries (1 = contact, 0 = no contact), the standard
error of the proportion is:

$$
\text{SEM} = \sqrt{\frac{p(1-p)}{N_{\text{ind}}}}
$$

Where:
- $p$ = simultaneous contact fraction
- $N_{\text{ind}}$ = number of independent samples (after autocorrelation correction)

This is the formula for the standard error of a binomial proportion.

### The Low Reliability Warning

When `N_independent < 10`, PolyzyMD warns about low statistical reliability.
This doesn't mean your results are wrong, but uncertainties may be
underestimated.

**Solutions:**
1. Use multiple independent replicates (recommended)
2. Run longer simulations
3. Interpret results with appropriate caution

## Worked Example: LipA Polymer Study

This example uses real data from simulations of Lipase A (LipA) from
*Bacillus subtilis* with various polymer conjugations.

### System Details

- **Enzyme**: Lipase A (181 residues)
- **Catalytic triad**: Ser77, His156, Asp133
- **Simulations**: 200 ns production, 100 ns equilibration discarded
- **Conditions**: 6 polymer compositions (3 replicates each)

### comparison.yaml Configuration

```yaml
name: "LipA_polymer_study"
description: "Effect of SBMA/EGMA polymer composition on LipA triad"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../noPoly_LipA_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../SBMA_100_DMSO/config.yaml"
    replicates: [1, 2]

  - label: "75% SBMA / 25% EGMA"
    config: "../SBMA_75_EGMA_25_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "50% SBMA / 50% EGMA"
    config: "../SBMA_50_EGMA_50_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "25% SBMA / 75% EGMA"
    config: "../SBMA_25_EGMA_75_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% EGMA"
    config: "../EGMA_100_DMSO/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "100ns"

catalytic_triad:
  name: "LipA_catalytic_triad"
  description: "Ser-His-Asp catalytic triad of Lipase A (Bacillus subtilis)"
  threshold: 3.5
  pairs:
    - label: "Asp133-His156"
      selection_a: "midpoint(resid 133 and name OD1 OD2)"
      selection_b: "resid 156 and name ND1"
    - label: "His156-Ser77"
      selection_a: "resid 156 and name NE2"
      selection_b: "resid 77 and name OG"
```

### Results Summary

| Condition | n | Asp133-His156 | His156-Ser77 | Simultaneous Contact |
|-----------|---|---------------|--------------|---------------------|
| No Polymer | 3 | 3.09 +/- 0.21 A | 4.03 +/- 1.07 A | 49.9 +/- 27.3% |
| 100% SBMA | 2 | 3.03 +/- 0.00 A | 3.75 +/- 0.86 A | 45.3 +/- 45.3% |
| 75% SBMA / 25% EGMA | 3 | 3.67 +/- 0.51 A | 4.92 +/- 0.50 A | 6.9 +/- 6.8% |
| 50% SBMA / 50% EGMA | 3 | 8.20 +/- 4.30 A | 4.23 +/- 0.20 A | 7.9 +/- 7.9% |
| 25% SBMA / 75% EGMA | 3 | 5.25 +/- 0.56 A | 5.50 +/- 1.07 A | 14.0 +/- 14.0% |
| **100% EGMA** | 3 | 3.29 +/- 0.30 A | 3.31 +/- 0.20 A | **50.6 +/- 17.0%** |

### Interpretation

1. **Best triad preservation**: 100% EGMA shows the highest simultaneous
   contact (50.6%) with the lowest variability (+/-17%), suggesting it
   provides consistent triad stabilization.

2. **No Polymer baseline**: 49.9% contact with high variability (+/-27.3%)
   indicates the unmodified enzyme has dynamic triad geometry.

3. **Mixed ratios disrupt triad**: 75% SBMA/25% EGMA, 50/50, and 25/75
   all show severely reduced contact (6.9-14.0%), suggesting these mixed
   polymer compositions interfere with triad geometry.

4. **Diagnostic per-pair analysis**: The 50/50 condition shows very high
   Asp133-His156 distance (8.20 A), indicating the Asp-His hydrogen bond
   is specifically disrupted. The His-Ser distance remains closer to
   normal (4.23 A).

### Conclusions

- Pure EGMA coating preserves triad geometry better than pure SBMA or no polymer
- Mixed polymer ratios severely disrupt the catalytic triad
- The Asp-His hydrogen bond appears most sensitive to disruption

These results suggest that polymer composition significantly affects enzyme
active site geometry, with implications for catalytic activity.

## Selection Syntax Details

### Standard MDAnalysis Selections

Most selections use standard MDAnalysis syntax:

```yaml
selection_a: "resid 77 and name OG"    # Serine OG oxygen
selection_b: "resid 156 and name NE2"  # Histidine NE2 nitrogen
```

Common patterns:
- `resid N` - residue number N
- `name X` - atom name X
- `resname ABC` - residue name ABC
- `and`, `or`, `not` - logical operators

### The midpoint() Function

For residues with two equivalent atoms (like Asp/Glu carboxylates), use
`midpoint()` to compute the geometric center:

```yaml
selection_a: "midpoint(resid 133 and name OD1 OD2)"
```

This computes:

$$
\mathbf{r}_{\text{mid}} = \frac{\mathbf{r}_{\text{OD1}} + \mathbf{r}_{\text{OD2}}}{2}
$$

**When to use midpoint():**
- Asp/Glu carboxylate groups (OD1/OD2 or OE1/OE2)
- Symmetric groups where either atom could participate

### The com() Function

For larger groups, use `com()` to compute the center of mass:

```yaml
selection_a: "com(resid 133)"  # Center of mass of entire Asp133
```

This weights atom positions by mass:

$$
\mathbf{r}_{\text{com}} = \frac{\sum_i m_i \mathbf{r}_i}{\sum_i m_i}
$$

**When to use com():**
- Tracking overall residue position
- Large binding site analysis
- When specific atoms aren't well-defined

## Multi-Replicate Analysis

### Why Replicates Matter

Single trajectories can be misleading due to:
- Rare conformational events
- Metastable states
- Stochastic variation

Multiple independent replicates provide:
- Reproducibility testing
- Robust uncertainty estimates
- Detection of outlier trajectories

### Recommended Practice

Run **at least 3 replicates** per condition. The aggregated analysis will:

1. Compute triad metrics for each replicate independently
2. Report mean and SEM across replicates
3. Show per-replicate values to assess variability

### Interpreting Replicate Variability

High replicate variability (large SEM) suggests:
- The system samples multiple conformational states
- Triad geometry is dynamic
- More replicates may be needed for statistical confidence

Low variability suggests:
- Consistent behavior across trajectories
- The observed state is representative

## Comparing Across Conditions

### Statistical Considerations

When comparing triad metrics across conditions:

1. **Use replicate-based statistics**: Compare mean values across replicates,
   not frame-by-frame values

2. **Report effect sizes**: A large difference may be meaningful even if
   not statistically significant with few replicates

3. **Consider biological significance**: A 10% difference in contact fraction
   might be functionally important even if p > 0.05

### Example Comparison

From the LipA study:

```
Condition         Contact      vs No Polymer
-------------------------------------------------
No Polymer        49.9%        (control)
100% EGMA         50.6%        +0.7% (similar)
100% SBMA         45.3%        -4.6% (slightly lower)
75/25 SBMA/EGMA   6.9%         -43.0% (much lower)
```

The 75/25 mixture shows a dramatic 43 percentage point reduction in triad
contact compared to the control - a biologically meaningful difference
regardless of statistical significance.

## Common Pitfalls

### Pitfall 1: Wrong Residue Numbering

**Symptom**: Very high distances (> 10 A) or selection errors

**Cause**: PDB residue numbers don't match your selections

**Solution**: Check your topology file's residue numbering. Use VMD or
PyMOL to visualize and verify:

```bash
# In VMD console
atomselect top "resid 77 and name OG"
```

### Pitfall 2: Incorrect Atom Names

**Symptom**: "Selection returned 0 atoms" error

**Cause**: Atom names differ between force fields

**Solution**: Check your topology for exact atom names. Common variations:

| Residue | AMBER | CHARMM |
|---------|-------|--------|
| His (delta) | ND1 | ND1 |
| His (epsilon) | NE2 | NE2 |
| Ser | OG | OG |
| Asp | OD1, OD2 | OD1, OD2 |

### Pitfall 3: Threshold Too Strict

**Symptom**: Near-zero contact fraction for active enzyme

**Cause**: 3.5 A may be too strict for some systems

**Solution**: Try 4.0 A threshold and compare. If the enzyme is
experimentally active, some threshold should show reasonable contact.

### Pitfall 4: Insufficient Equilibration

**Symptom**: Contact fraction drifts over time; different equilibration
times give very different results

**Cause**: Including non-equilibrated frames

**Solution**: Use RMSD analysis to determine equilibration time. The
system should be equilibrated before triad analysis.

### Pitfall 5: Ignoring Per-Pair Diagnostics

**Symptom**: Low simultaneous contact but unclear why

**Cause**: Not examining which pair is disrupted

**Solution**: Always check per-pair distances. This reveals whether:
- Both pairs are moderately disrupted
- One pair is specifically broken
- The issue is distance vs. variability

## Python API

### Programmatic Analysis

```python
from polyzymd.config.loader import load_config
from polyzymd.compare.config import ComparisonConfig
from polyzymd.analysis.triad import CatalyticTriadAnalyzer

# Load configs
comp_config = ComparisonConfig.from_yaml("comparison.yaml")
sim_config = load_config(comp_config.conditions[0].config)

# Create analyzer
analyzer = CatalyticTriadAnalyzer(
    config=sim_config,
    triad_config=comp_config.catalytic_triad,
    equilibration="100ns",
)

# Single replicate
result = analyzer.compute(replicate=1)
print(f"Simultaneous contact: {result.simultaneous_contact_fraction * 100:.1f}%")

# Aggregated
agg_result = analyzer.compute_aggregated(replicates=[1, 2, 3])
print(f"Mean contact: {agg_result.overall_simultaneous_contact * 100:.1f} "
      f"+/- {agg_result.sem_simultaneous_contact * 100:.1f}%")
```

### Accessing Detailed Results

```python
# Per-pair statistics
for pair in result.pair_results:
    print(f"{pair.pair_label}:")
    print(f"  Mean: {pair.mean_distance:.2f} A")
    print(f"  Std:  {pair.std_distance:.2f} A")
    if pair.fraction_below_threshold:
        print(f"  Contact: {pair.fraction_below_threshold * 100:.1f}%")

# Autocorrelation info
if result.sim_contact_correlation_time:
    print(f"Correlation time: {result.sim_contact_correlation_time:.1f} "
          f"{result.sim_contact_correlation_time_unit}")
    print(f"N independent: {result.sim_contact_n_independent}")
```

### Loading Saved Results

```python
from polyzymd.analysis.results.triad import TriadResult, TriadAggregatedResult

# Load single replicate result
result = TriadResult.load("analysis/triad/run_1/triad_LipA_eq100ns.json")
print(result.summary())

# Load aggregated result
agg = TriadAggregatedResult.load("analysis/triad/aggregated/triad_LipA_reps1-3_eq100ns.json")
print(agg.summary())
```

## Plotting and Visualization

```{admonition} TODO: Plotting Documentation
:class: warning

Detailed documentation for plotting catalytic triad results is planned for a
future release. In the meantime, you can:

1. Load JSON results and create custom plots with matplotlib
2. Use VMD/PyMOL to visualize triad geometry along trajectories
3. Export distance timeseries for external analysis

See the [RMSF plotting documentation](analysis_rmsf_quickstart.md#comparing-two-conditions)
for general guidance on comparing analysis results across conditions.
```

## References

### Enzyme Mechanism

**Hedstrom L.** (2002) "Serine Protease Mechanism and Specificity."
*Chemical Reviews* 102:4501-4524.
https://doi.org/10.1021/cr000033x

Classic review of serine protease catalytic mechanism.

**Blow DM.** (1976) "Structure and Mechanism of Chymotrypsin."
*Accounts of Chemical Research* 9:145-152.

Foundational paper on catalytic triad structure.

### MD Analysis Best Practices

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2018) "Best Practices for Quantification of Uncertainty and Sampling Quality
in Molecular Simulations." *Living Journal of Computational Molecular Science*
1(1):5067. https://doi.org/10.33011/livecoms.1.1.5067

Definitive guide to uncertainty quantification in MD - the basis for
PolyzyMD's autocorrelation analysis.

### H-Bond Geometry

**Jeffrey GA, Saenger W.** (1991) *Hydrogen Bonding in Biological Structures.*
Springer-Verlag.

Comprehensive reference for hydrogen bond geometry in biomolecules.

## See Also

- [Quick Start Guide](analysis_triad_quickstart.md) - Get results fast
- [RMSF Analysis](analysis_rmsf_quickstart.md) - Analyze flexibility
- [Comparing Conditions](analysis_compare_conditions.md) - Statistical comparisons
- [Statistical Best Practices](analysis_rmsf_best_practices.md) - Detailed statistics guide
