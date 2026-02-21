# Binding Free Energy Analysis

Quantify **how much each polymer type preferentially contacts specific protein
residue groups** in physically meaningful energy units (kcal/mol or kJ/mol).

This analysis answers the experimentalist's core question:
*"Which polymer chemistry has the most favorable interaction with the enzyme
surface—and by how much?"*

```{note}
Binding free energy analysis is a **post-processing step**. It requires:
1. Completed contacts analysis with
   [`compute_binding_preference: true`](analysis_binding_preference.md)
2. Cached binding preference files on disk (produced automatically by `polyzymd
   compare contacts`)
```

## Scientific Motivation

### From Enrichment to Free Energy

[Binding preference analysis](analysis_binding_preference.md) reports a
dimensionless **enrichment**:

```{math}
\text{Enrichment} = \frac{\text{contact\_share}}{\text{expected\_share}} - 1
```

Enrichment is intuitive (+0.5 means 50% more contacts than random) but has no
units. It cannot answer: *"Is a +0.5 enrichment for aromatic residues
thermodynamically significant?"*

Binding free energy converts enrichment into **Gibbs free energy differences
(ΔΔG)** via Boltzmann inversion. Because polyzymd simulations run in the NPT
ensemble, the correct thermodynamic potential is the Gibbs free energy.

### The Core Formula

```{math}
\Delta\Delta G = -k_B T \ln\!\left(\frac{\text{contact\_share}}{\text{expected\_share}}\right)
```

Equivalently, since `contact_share / expected_share = enrichment + 1`:

```{math}
\Delta\Delta G = -k_B T \ln(\text{enrichment} + 1)
```

| ΔΔG | Meaning |
|-----|---------|
| **< 0** | Preferential binding — the polymer favors this residue group |
| **= 0** | Neutral — contact frequency matches random |
| **> 0** | Avoidance — the polymer contacts this group less than expected |

### Relation to Enrichment

ΔΔG is the **exact** Boltzmann-inverted version of enrichment. The enrichment
score is a first-order Taylor expansion (without units); ΔΔG is exact and
carries physical units.

The relationship is:

```{math}
\text{enrichment} = e^{-\Delta\Delta G / k_B T} - 1
```

### Uncertainty (Delta Method)

Uncertainty is propagated analytically via the delta method:

```{math}
\sigma(\Delta\Delta G) = k_B T \sqrt{\left(\frac{\sigma_\text{cs}}{\text{cs}}\right)^2 + \left(\frac{\sigma_\text{es}}{\text{es}}\right)^2}
```

where `σ_cs` is the SEM of the contact share across replicates, and
`σ_es ≈ 0` (expected share is computed from a single PDB SASA calculation).

### Temperature and Comparability

ΔΔG computed at temperature *T* is **not directly comparable** to ΔΔG at
temperature *T'* — the k_BT scale factor differs. When conditions span multiple
simulation temperatures, pairwise statistics are automatically suppressed for
cross-temperature pairs.

## Quick Start

### Step 1: Enable Binding Preference in contacts analysis

Binding free energy is a post-processor — it reads the files produced by
contacts analysis. Ensure binding preference is enabled when you run contacts:

```yaml
# comparison.yaml
name: "SBMA vs EGMA Free Energy Study"
control: "No Polymer"

structures:
  enzyme_pdb: "structures/enzyme.pdb"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% EGMA"
    config: "../egma_100/config.yaml"
    replicates: [1, 2, 3]

analysis_settings:
  contacts:
    name: "polymer_contacts"
    polymer_selection: "resname SBM EGM"
    protein_selection: "protein"
    cutoff: 4.5
    compute_binding_preference: true
    surface_exposure_threshold: 0.2
    enzyme_pdb_for_sasa: "structures/enzyme.pdb"
    include_default_aa_groups: true

  # Binding free energy: reads cached binding preference data
  binding_free_energy:
    units: "kcal/mol"              # or "kJ/mol"
    surface_exposure_threshold: 0.2

comparison_settings:
  binding_free_energy:
    fdr_alpha: 0.05
```

### Step 2: Run contacts analysis first (if not already done)

```bash
polyzymd compare contacts -f comparison.yaml
```

This generates the cached binding preference files under each condition's
`analysis/contacts/` directory.

### Step 3: Run binding free energy analysis

```bash
polyzymd compare binding-free-energy -f comparison.yaml
```

## Example Output

```
Binding Free Energy Comparison: SBMA vs EGMA Free Energy Study
================================================================================
Formula : ΔΔG = -k_B·T · ln(contact_share / expected_share)
Units   : kcal/mol
Equilibration: 10ns

ΔΔG Summary by Condition (sign: negative = preferential binding)
--------------------------------------------------------------------------------

  100% SBMA  (T = 300.0 K, n = 3)
  Polymer      AA Group               ΔΔG           ±σ    N_rep
  ---------------------------------------------------------------
  SBM          aromatic            -0.243       +0.031       3
  SBM          charged_negative    +0.091       +0.022       3
  SBM          charged_positive    -0.013       +0.015       3
  SBM          nonpolar            +0.044       +0.019       3
  SBM          polar               -0.001       +0.001       3

  100% EGMA  (T = 300.0 K, n = 3)
  Polymer      AA Group               ΔΔG           ±σ    N_rep
  ---------------------------------------------------------------
  EGM          aromatic            -0.523       +0.034       3
  EGM          charged_negative    +0.319       +0.028       3
  EGM          charged_positive    +0.112       +0.049       3
  EGM          nonpolar            -0.166       +0.021       3
  EGM          polar               +0.146       +0.010       3

Pairwise ΔΔG Differences (ΔΔG_B − ΔΔG_A)
--------------------------------------------------------------------------------

  100% SBMA  →  100% EGMA
  Polymer      AA Group               ΔΔG_A      ΔΔG_B    ΔΔG_B−A    p-value
  -----------------------------------------------------------------------------
  EGM          aromatic              N/A       -0.523        N/A         --
  SBM          aromatic           -0.243          N/A        N/A         --
```

### Interpreting the Numbers

From the example:

- **EGMA (EGM) aromatic: ΔΔG = −0.52 kcal/mol** — strong preferential
  contact with aromatic residues (about half a kcal/mol more favorable than
  random). This is equivalent to the enrichment score of +0.90 from binding
  preference analysis.

- **EGMA charged_negative: ΔΔG = +0.32 kcal/mol** — avoidance of charged
  negative residues.

- **SBMA (SBM) charged_positive: ΔΔG = −0.01 kcal/mol** — near-neutral,
  consistent with SBMA's balanced zwitterionic character.

```{tip}
A rule of thumb from protein biophysics: differences of **|ΔΔG| > 0.5 kcal/mol**
are generally considered thermodynamically significant at room temperature. Smaller
values may still be real but require larger replicate sets to distinguish from noise.
```

## Configuration Reference

### `analysis_settings.binding_free_energy`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `units` | str | `"kcal/mol"` | Energy units: `"kcal/mol"` or `"kJ/mol"` |
| `surface_exposure_threshold` | float | `0.2` | Minimum relative SASA to consider a residue exposed (0.0–1.0) |
| `protein_partitions` | dict | `null` | Custom named partitions (see below) |

### `comparison_settings.binding_free_energy`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `fdr_alpha` | float | `0.05` | Benjamini-Hochberg FDR threshold for pairwise t-tests |

### Full YAML Example

```yaml
analysis_settings:
  binding_free_energy:
    units: "kcal/mol"
    surface_exposure_threshold: 0.2
    # Optional: custom partitions (groups must be defined under contacts.protein_groups)
    protein_partitions:
      catalytic_regions:
        - catalytic_triad
        - oxyanion_hole
      lid_helices:
        - lid_helix_5
        - lid_helix_10

comparison_settings:
  binding_free_energy:
    fdr_alpha: 0.05
```

## CLI Reference

```
polyzymd compare binding-free-energy [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-f, --file PATH` | Path to `comparison.yaml` (default: `comparison.yaml`) |
| `--units [kcal/mol\|kJ/mol]` | Override energy units |
| `--fdr-alpha FLOAT` | Override FDR alpha (default: 0.05) |
| `--recompute` | Force recompute (clears binding preference cache) |
| `--format [table\|markdown\|json]` | Output format (default: `table`) |
| `-o, --output PATH` | Save output to file |
| `-q, --quiet` | Suppress INFO messages |
| `--debug` | Enable DEBUG logging |
| `--eq-time TEXT` | Override equilibration time (e.g., `10ns`) |

### Usage Examples

```bash
# Default: table output to console
polyzymd compare binding-free-energy -f comparison.yaml

# kJ/mol units
polyzymd compare binding-free-energy --units kJ/mol

# Save markdown report
polyzymd compare binding-free-energy --format markdown -o report.md

# JSON for downstream processing
polyzymd compare binding-free-energy --format json -o bfe_result.json

# Stricter p-value threshold
polyzymd compare binding-free-energy --fdr-alpha 0.01

# Verbose debugging
polyzymd compare binding-free-energy --debug
```

## Custom Protein Partitions

You can analyze specific structural regions of interest alongside the default
amino acid class groups. Define custom groups in the contacts settings, then
group them into named partitions in the binding free energy settings.

```yaml
analysis_settings:
  contacts:
    compute_binding_preference: true
    include_default_aa_groups: true
    protein_groups:
      catalytic_triad: [77, 133, 156]
      oxyanion_hole: [10, 77]
      lid_helix_5: [141, 142, 143, 144, 145, 146, 147, 148, 149, 150]
      lid_helix_10: [281, 282, 283, 284, 285, 286, 287, 288, 289, 290]

  binding_free_energy:
    units: "kcal/mol"
    protein_partitions:
      # Must be mutually exclusive groups within each partition
      catalytic_regions:
        - catalytic_triad
        - oxyanion_hole
      lid_helices:
        - lid_helix_5
        - lid_helix_10
```

The pairwise and per-condition tables will include one row per
(polymer_type, partition_element) pair for each custom partition.

```{important}
Partition groups must be **mutually exclusive** (no residue in two groups of the
same partition). PolyzyMD validates this at config load time. See
[Binding Preference Analysis](analysis_binding_preference.md#partition-based-architecture)
for a full explanation.
```

## Multi-Temperature Studies

When conditions were simulated at different temperatures (e.g., a temperature
stability study), ΔΔG values are reported per condition but **pairwise
statistics are suppressed** between conditions at different temperatures:

```
Note: pairwise statistics suppressed for cross-temperature pairs.
Temperatures: 300.0 K (100% SBMA, 100% EGMA), 320.0 K (High-T SBMA)
```

Per-condition ΔΔG values are still shown — they are valid within each
temperature group. Only the pairwise ΔΔG_B−A values and p-values are omitted
for cross-temperature pairs, since k_BT differs.

## Output File Location

When you use `-o / --output`, the formatted report is saved to the specified
path. A JSON copy of the result is also saved to the comparison's `results/`
directory:

```
project/
└── comparison_results/
    └── binding_free_energy/
        └── binding_free_energy_YYYYMMDD_HHMMSS.json
```

The JSON file can be reloaded for downstream processing:

```python
from polyzymd.compare.results.binding_free_energy import BindingFreeEnergyResult

result = BindingFreeEnergyResult.load("binding_free_energy_20260221_143000.json")

# Access per-condition data
for summary in result.conditions:
    print(f"\n{summary.label} ({summary.temperature_K} K)")
    for entry in summary.entries:
        if entry.delta_G is not None:
            print(
                f"  {entry.polymer_type} × {entry.protein_group}: "
                f"ΔΔG = {entry.delta_G:+.3f} ± {entry.delta_G_uncertainty:.3f} "
                f"{result.units}"
            )

# Access pairwise comparisons
for pair in result.pairwise_comparisons:
    if not pair.cross_temperature and pair.p_value is not None:
        print(
            f"{pair.condition_a} → {pair.condition_b} | "
            f"{pair.polymer_type} × {pair.protein_group}: "
            f"ΔΔG_B−A = {pair.delta_delta_G:+.3f}, p = {pair.p_value:.4f}"
        )
```

## Statistical Notes

### Pairwise t-tests

For each shared (polymer_type, protein_group) pair between two conditions at
the same temperature, an independent two-sample t-test is performed on the
per-replicate ΔΔG distributions. Requires **at least 2 replicates per
condition**.

### FDR Correction

P-values are corrected for multiple comparisons using the
**Benjamini-Hochberg** false discovery rate procedure at the specified alpha
(default 0.05). Adjusted p-values are stored in the JSON output but not
currently shown in the table or markdown formatters.

### Statistical Inefficiency and MetricType

Contact share is a **mean-based metric** (frame averages, not fluctuations).
Its mean converges regardless of autocorrelation, though the SEM reported here
is a naive cross-replicate SEM (not autocorrelation-corrected). With 3+
independent replicates this is generally conservative and appropriate.

See [Statistics Best Practices](analysis_statistics_best_practices.md) for
guidance on replicate counts and equilibration.

## Worked Example: Comparing Polymer Stabilization

### Scientific Question

*"Does EGMA or SBMA show stronger preferential binding to the enzyme's aromatic
surface residues? Can we quantify the difference in kcal/mol?"*

### Setup

```yaml
# comparison.yaml
name: "EGMA vs SBMA Aromatic Preference"
control: "100% SBMA"

conditions:
  - label: "100% SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3, 4, 5]

  - label: "100% EGMA"
    config: "../egma_100/config.yaml"
    replicates: [1, 2, 3, 4, 5]

analysis_settings:
  contacts:
    polymer_selection: "resname SBM EGM"
    cutoff: 4.5
    compute_binding_preference: true
    surface_exposure_threshold: 0.2
    enzyme_pdb_for_sasa: "structures/lipase.pdb"
    include_default_aa_groups: true

  binding_free_energy:
    units: "kcal/mol"

comparison_settings:
  binding_free_energy:
    fdr_alpha: 0.05
```

### Run

```bash
# Step 1: Generate binding preference data
polyzymd compare contacts -f comparison.yaml

# Step 2: Compute ΔΔG
polyzymd compare binding-free-energy -f comparison.yaml --format markdown -o bfe_report.md
```

### Interpretation Table

| Observation | Interpretation |
|-------------|----------------|
| EGMA aromatic: ΔΔG = −0.52 ± 0.03 kcal/mol | EGMA has a significant thermodynamic preference for aromatic surface residues |
| SBMA aromatic: ΔΔG = −0.24 ± 0.03 kcal/mol | SBMA also preferentially contacts aromatics, but less strongly |
| ΔΔG_EGMA − ΔΔG_SBMA = −0.28 kcal/mol, p < 0.01 | The difference is statistically significant with 5 replicates |
| SBMA charged groups: ΔΔG ≈ 0 | SBMA's zwitterionic balance leads to near-neutral electrostatic preference |

## Troubleshooting

### "No 'binding_free_energy' in analysis_settings"

Ensure your `comparison.yaml` has a `binding_free_energy:` block under
`analysis_settings`:

```yaml
analysis_settings:
  contacts:
    compute_binding_preference: true
    # ... other contacts settings
  binding_free_energy:     # ← required
    units: "kcal/mol"
```

### "No binding preference data found for condition"

The binding preference cache files are missing. Run contacts analysis first:

```bash
polyzymd compare contacts -f comparison.yaml
```

If the files still aren't found, check that `compute_binding_preference: true`
was set in the contacts analysis settings when you ran contacts.

### "N/A" entries in the output

ΔΔG is undefined when `contact_share = 0` (no contacts observed) or
`expected_share = 0` (no exposed residues in that group). These are shown as
`N/A` and excluded from pairwise comparisons.

### "Pairwise statistics suppressed for cross-temperature pairs"

Your conditions were simulated at different temperatures. ΔΔG at different
temperatures is not directly comparable — k_BT differs. You can still compare
ΔΔG within a temperature group, but not across groups. Consider re-running all
conditions at the same temperature if cross-condition comparisons are needed.

### kJ/mol vs kcal/mol

To convert: 1 kcal/mol = 4.184 kJ/mol. The thermal energy at 300 K is:
- 0.596 kcal/mol (k_BT)
- 2.494 kJ/mol (k_BT)

Values smaller than ~0.3 k_BT (0.18 kcal/mol or 0.75 kJ/mol) are likely
within thermal noise.

## See Also

- [Binding Preference Analysis](analysis_binding_preference.md) — the prerequisite enrichment analysis
- [Contacts Analysis Quick Start](analysis_contacts_quickstart.md) — basic contacts
- [Statistics Best Practices](analysis_statistics_best_practices.md) — replicate planning
- [Comparing Conditions](analysis_compare_conditions.md) — multi-condition workflows
- [Extending Comparators](extending_comparators.md) — add custom comparators
