# Polymer Affinity Score Analysis

Quantify the **total polymer-protein interaction strength** by combining
per-contact free energies with contact counts across all residue groups.

This analysis answers the experimentalist's question:
*"Which polymer composition has the strongest overall adhesion to the enzyme
surface?"*

```{note}
The polymer affinity score is a **comparative scoring metric for ranking polymer
compositions**. It is derived entirely from cached binding preference data — no
additional trajectory access is needed. Run contacts analysis with
`compute_binding_preference: true` first.
```

## Scientific Motivation

### Beyond Per-Contact Free Energy

[Binding free energy analysis](analysis_binding_free_energy.md) reports the
thermodynamic preference of each polymer type for each residue group (e.g.,
"EGMA prefers aromatics by -0.52 kT"). However, this per-contact value
misses the **polyvalent binding effect**: many contacts happen simultaneously.

A polymer with a modest per-contact affinity but many simultaneous contacts
may have a stronger total interaction than one with strong per-contact affinity
but few contacts. The polymer affinity score captures this by summing:

```{math}
\text{Score} = N \times \Delta\Delta G
```

where:
- **N** = mean number of simultaneous contacts with a residue group
- **ΔΔG** = per-contact free energy from Boltzmann inversion

### The Formula

Per (polymer_type, protein_group):

```{math}
N_g = \text{mean\_contact\_fraction} \times n_\text{exposed\_in\_group}
```

```{math}
\Delta\Delta G_g = -\ln\!\left(\frac{\text{contact\_share}}{\text{expected\_share}}\right) \quad [k_\mathrm{b}T]
```

```{math}
S_g = N_g \times \Delta\Delta G_g
```

Per polymer type:

```{math}
S_\text{polymer} = \sum_g S_g
```

Per condition (headline score):

```{math}
S_\text{total} = \sum_\text{polymer\_types} S_\text{polymer}
```

### Sign Convention

| Score | Meaning |
|-------|---------|
| **< 0** | Net attractive polymer-protein interaction |
| **= 0** | Neutral (contacts match random expectation) |
| **> 0** | Net repulsive (polymer avoids protein surface) |

More negative = stronger overall polymer-protein adhesion.

```{important}
**Independence assumption disclaimer.** The affinity score sums individual
contact free energies assuming thermodynamic independence of contacts.
Conformational coupling, surface exclusion effects, and reference state shifts
violate this assumption. The absolute values are therefore not rigorous
thermodynamic ΔG. However, the **relative differences** between polymer
compositions are meaningful for ranking purposes, as systematic biases cancel
when comparing compositions analyzed under identical conditions.
```

### Interpretation for Experimentalists

"Among compositions with similar RMSF (structural stability), the one with the
most negative affinity score has the strongest polymer-protein adhesion — good
for thermal resistance, but check that it doesn't block the active site."

The affinity score should always be interpreted alongside:
- **RMSF** — does strong adhesion stabilize or rigidify the enzyme?
- **Triad geometry** — does polymer contact disrupt the catalytic triad?
- **Substrate distance** — does the polymer occlude the active site?

## Quick Start

### Step 1: Ensure binding preference data exists

Polymer affinity analysis reads cached binding preference files. These are
produced automatically by `polyzymd compare contacts`:

```yaml
# comparison.yaml
analysis_settings:
  contacts:
    compute_binding_preference: true
    surface_exposure_threshold: 0.2
    enzyme_pdb_for_sasa: "structures/enzyme.pdb"
    include_default_aa_groups: true
```

```bash
# Generate binding preference cache (if not already done)
polyzymd compare contacts -f comparison.yaml
```

### Step 2: Add polymer affinity settings (optional)

The polymer affinity score works with zero configuration — defaults are
suitable for most use cases. To customize:

```yaml
# comparison.yaml
analysis_settings:
  polymer_affinity:
    surface_exposure_threshold: 0.2

comparison_settings:
  polymer_affinity:
    fdr_alpha: 0.05
```

### Step 3: Run polymer affinity analysis

```bash
polyzymd compare polymer-affinity -f comparison.yaml
```

## Example Output

```
Polymer Affinity Score: LipA 363K Comparison
================================================================================
Formula : Score = N_contacts × ΔΔG, where ΔΔG = -ln(contact_share / expected_share)
Units   : kT (dimensionless, in units of k_bT)
Sign    : Negative = net attractive interaction
Note    : Assumes thermodynamic independence of contacts (see docs for caveats)

Per-Condition Summary (sign: more negative = stronger adhesion)
--------------------------------------------------------------------------------

  100% SBMA  (T = 363.0 K, n = 5)
  Polymer  AA Group            N_contacts     ΔΔG/kT    Score/kT
  ----------------------------------------------------------------
  SBM      aromatic                12.3      -0.19       -2.34
  SBM      charged_negative         8.7      +0.08       +0.70
  SBM      charged_positive         6.1      -0.01       -0.06
  SBM      nonpolar                15.2      +0.04       +0.61
  SBM      polar                    9.4      -0.00       -0.00
  ----------------------------------------------------------------
  SBM      TOTAL                                         -1.09

  SBMA-EGPMA 5%  (T = 363.0 K, n = 5)
  Polymer  AA Group            N_contacts     ΔΔG/kT    Score/kT
  ----------------------------------------------------------------
  SBM      aromatic                14.1      -0.22       -3.10
  SBM      nonpolar                17.3      +0.02       +0.35
  EGP      aromatic                 1.2      -0.45       -0.54
  ...
  ----------------------------------------------------------------
  ALL      TOTAL                                         -4.21

Pairwise Score Differences (Score_B - Score_A)
--------------------------------------------------------------------------------
  100% SBMA → SBMA-EGPMA 5%  :  ΔScore = -3.12 kT  (p = 0.003 **)
```

### Interpreting the Numbers

From the example:

- **SBMA-EGPMA 5% total score: -4.21 kT** — stronger overall adhesion than
  100% SBMA (-1.09 kT). The 5% EGPMA comonomer increases total polymer-protein
  binding through aromatic contacts.

- **SBM aromatic N_contacts=14.1, ΔΔG=-0.22 kT** — each of ~14 simultaneous
  aromatic contacts contributes -0.22 kT, giving a group score of -3.10 kT.

- **Pairwise ΔScore = -3.12 kT, p = 0.003** — the difference is statistically
  significant.

```{tip}
Scores are in kT units. At 363 K, 1 kT = 0.72 kcal/mol = 3.02 kJ/mol.
A total score difference of |ΔScore| > 2 kT is generally considered
meaningful.
```

## Configuration Reference

### `analysis_settings.polymer_affinity`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `surface_exposure_threshold` | float | `0.2` | Minimum relative SASA for exposed residues |
| `polymer_selection` | str | `"resname SBM EGM EGP"` | MDAnalysis polymer selection |
| `protein_selection` | str | `"protein"` | MDAnalysis protein selection |
| `cutoff` | float | `4.5` | Contact distance cutoff (Angstroms) |
| `compute_binding_preference` | bool | `true` | Must be true for affinity analysis |

### `comparison_settings.polymer_affinity`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `fdr_alpha` | float | `0.05` | Benjamini-Hochberg FDR threshold |

## CLI Reference

```
polyzymd compare polymer-affinity [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-f, --file PATH` | Path to `comparison.yaml` (default: `comparison.yaml`) |
| `--recompute` | Force recompute |
| `--format [table\|markdown\|json]` | Output format (default: `table`) |
| `-o, --output PATH` | Save output to file |
| `-q, --quiet` | Suppress INFO messages |
| `--debug` | Enable DEBUG logging |
| `--eq-time TEXT` | Override equilibration time |

### Usage Examples

```bash
# Default: console table
polyzymd compare polymer-affinity -f comparison.yaml

# Save markdown report
polyzymd compare polymer-affinity --format markdown -o affinity_report.md

# JSON for downstream processing
polyzymd compare polymer-affinity --format json -o affinity_result.json
```

## Per-Replicate Computation

When per-replicate binding preference files (`binding_preference_rep{N}.json`)
are available, the polymer affinity score computes exact per-replicate scores:

```{math}
N_{\text{rep}} = \text{contact\_fraction}_{\text{rep}} \times n_\text{exposed}
```

```{math}
\Delta\Delta G_{\text{rep}} = -\ln(\text{enrichment}_{\text{rep}} + 1)
```

```{math}
S_{\text{rep}} = N_{\text{rep}} \times \Delta\Delta G_{\text{rep}}
```

Mean and SEM are computed across replicates. Pairwise t-tests use these
per-replicate score distributions.

When per-replicate files are unavailable, the comparator falls back to
analytical error propagation:

```{math}
\sigma(S) = \sqrt{(N \cdot \sigma_{\Delta\Delta G})^2 + (\Delta\Delta G \cdot \sigma_N)^2}
```

## Multi-Temperature Studies

Like binding free energy analysis, the polymer affinity score uses kT units
that already account for temperature. However, the underlying contact
distributions may differ between temperatures, so **pairwise statistics are
suppressed for cross-temperature pairs**.

Per-condition scores are always shown — they are valid within each temperature.

## Relation to Other Analyses

| Analysis | What It Measures | Relation to Affinity Score |
|----------|-----------------|--------------------------|
| **Contacts** | Which residues are touched | Provides the raw contact data |
| **Binding Preference** | Enrichment per residue group | Provides enrichment → ΔΔG |
| **Binding Free Energy** | Per-contact ΔΔG | Component of the score (without N) |
| **Polymer Affinity Score** | Total interaction strength | N × ΔΔG summed over all groups |
| **RMSF** | Structural flexibility | Complementary: stability metric |
| **Triad** | Active site geometry | Complementary: function metric |

## Output File Location

Results are saved to the comparison's `results/` directory:

```
comparison_output/
└── results/
    └── polymer_affinity_score_YYYYMMDD_HHMMSS.json
```

Reload for downstream processing:

```python
from polyzymd.compare.results.polymer_affinity import PolymerAffinityScoreResult

result = PolymerAffinityScoreResult.load("polymer_affinity_score_20260228_120000.json")

# Access per-condition scores
for summary in result.conditions:
    print(f"\n{summary.label}: total_score = {summary.total_score:+.2f} kT")
    for entry in summary.entries:
        print(
            f"  {entry.polymer_type} × {entry.protein_group}: "
            f"N={entry.n_contacts:.1f}, ΔΔG={entry.delta_G:+.3f}, "
            f"Score={entry.score:+.3f} kT"
        )

# Access pairwise comparisons
for pair in result.pairwise_comparisons:
    if pair.p_value is not None:
        sig = "**" if pair.p_value < 0.01 else "*" if pair.p_value < 0.05 else ""
        print(
            f"{pair.condition_a} → {pair.condition_b}: "
            f"ΔScore = {pair.delta_score:+.2f} kT, p = {pair.p_value:.4f} {sig}"
        )
```

## Troubleshooting

### "No binding preference data found"

Run contacts analysis with binding preference enabled first:

```bash
polyzymd compare contacts -f comparison.yaml
```

Ensure `compute_binding_preference: true` is set in your contacts analysis
settings.

### Very large scores (|S| > 50 kT)

Large scores usually indicate many contacts with strong enrichment/depletion.
Check whether your polymer selection is too broad (including solvent) or your
contact cutoff is too large. Also verify that `n_exposed_in_group` values are
reasonable — very large groups amplify the score.

### Score of exactly 0.0

This occurs when either `mean_contact_fraction = 0` (no contacts) or
`enrichment = 0` (contacts exactly match random expectation). Check the
binding preference data for the affected condition.

### Different scores from BFE analysis

The affinity score and BFE analysis use the same ΔΔG formula but different
aggregation: BFE reports per-contact ΔΔG, while the affinity score multiplies
by N_contacts. The per-contact ΔΔG values should match exactly; differences
indicate a cache or configuration mismatch.

## See Also

- [Binding Free Energy Analysis](analysis_binding_free_energy.md) — per-contact ΔΔG (the component metric)
- [Binding Preference Analysis](analysis_binding_preference.md) — the prerequisite enrichment analysis
- [Contacts Analysis Quick Start](analysis_contacts_quickstart.md) — basic contact computation
- [Statistics Best Practices](analysis_statistics_best_practices.md) — replicate planning
- [Comparing Conditions](analysis_compare_conditions.md) — multi-condition workflows
- [Extending Comparators](extending_comparators.md) — add custom comparators
