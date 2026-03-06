# Exposure Dynamics Analysis

Understand **how your polymer interacts with transiently exposed protein
residues** using chaperone event detection, refolding acceleration analysis,
and SASA-weighted chaperone selectivity. This module answers three
complementary questions that static contact analysis cannot:

1. *"Does my polymer act as a molecular chaperone — providing contact during
   transient unfolding events that would otherwise be unassisted?"*
2. *"Does polymer contact accelerate or slow refolding of transiently exposed
   residues?"* (refolding acceleration ratio ρ)
3. *"Does my polymer preferentially contact certain amino acid groups
   **during chaperone events**, controlling for instantaneous surface
   availability?"* (chaperone selectivity ΔG_sel^chap)

```{note}
Exposure dynamics analysis builds on [contacts analysis](analysis_contacts_quickstart.md)
and uses SASA computed per-frame from MD trajectories, not from a static PDB.
Run contacts analysis first, then enable exposure dynamics in your comparison
configuration.
```

## Scientific Motivation

### Why static analysis is insufficient

Standard binding preference analysis (see
[Binding Preference Analysis](analysis_binding_preference.md)) uses SASA
from a single static PDB to define which residues are "surface accessible."
This works well for residues that are consistently exposed, but misses a
crucial class of events: **transient exposure**.

A residue that is buried 90% of the time but exposes briefly may be the most
biologically relevant interaction site — or it may represent a mechanically
vulnerable region that the polymer can protect or destabilize. A static SASA
filter either includes it (at a permissive threshold) or excludes it entirely
(at a stringent threshold), conflating residues with fundamentally different
dynamics.

### What chaperone analysis captures

**Chaperone fraction** asks: of all the exposure episodes a residue undergoes
(contiguous stretches of frames where it is exposed), what fraction have
polymer contact during the episode? A high chaperone fraction for a transient
residue suggests the polymer is co-localizing during the exposure window.

**Refolding acceleration ratio ρ** asks: does polymer contact during
chaperone events correlate with *faster* refolding (shorter exposure duration)?
ρ > 1 indicates the polymer accelerates refolding; ρ < 1 indicates trapping.

**Chaperone selectivity ΔG_sel^chap** asks: during chaperone events, does
the polymer preferentially contact certain amino acid groups relative to their
SASA share? This is analogous to the equilibrium binding selectivity
(ΔG_sel) but restricted to transient unfolding events.

---

## The Three Metrics

### Chaperone fraction

#### Residue stability classification

Before computing chaperone fractions, each residue is classified by its
**exposure fraction** — the proportion of trajectory frames where its relative
SASA exceeds the threshold θ:

| Classification | Condition | Meaning |
|---------------|-----------|---------|
| `stably_exposed` | exposure_fraction ≥ `transient_upper` (default 0.80) | Consistently on the surface |
| `transient` | `transient_lower` < exposure_fraction < `transient_upper` | Fluctuates between buried and exposed |
| `stably_buried` | exposure_fraction ≤ `transient_lower` (default 0.20) | Consistently buried |

Only **transient** residues are used for condition-level chaperone fraction
reporting, because chaperone events are only meaningful for residues that
genuinely fluctuate.

#### Event detection

For each residue, the trajectory is segmented into **exposure windows**: maximal
contiguous runs of frames where the residue is exposed (relative SASA > θ).
Windows shorter than `min_event_length` frames are discarded to remove
single-frame SASA noise. Each window must be flanked by buried frames
(start > 0 and end < n_frames - 1).

Each exposure window is then classified:

- **Chaperone event**: at least one polymer contact occurs during the window
- **Unassisted event**: no polymer contact occurs during the window

#### Per-residue chaperone fraction

```{math}
\text{chaperone\_fraction} = \frac{n_{\text{chaperone events}}}{n_{\text{chaperone events}} + n_{\text{unassisted events}}}
```

A value of 1.0 means every exposure window for this residue had polymer
contact; 0.0 means none did.

#### Condition-level chaperone fraction

The reported condition-level value is the mean chaperone fraction over all
transient residues in a replicate:

```{math}
\text{condition chaperone fraction} = \frac{1}{|R_{\text{transient}}|} \sum_{i \in R_{\text{transient}}} \text{chaperone\_fraction}_i
```

This value is computed per replicate; the mean ± SEM across replicates is
reported in comparisons.

### Refolding acceleration ratio ρ(P, G)

For each polymer type P and amino acid group G, the refolding acceleration
ratio compares mean event durations:

```{math}
\rho(P, G) = \frac{\langle\tau\rangle_{\text{unassisted}}(G)}{\langle\tau\rangle_{\text{chap}}(P, G)}
```

where:
- ⟨τ⟩_chap(P, G) is the mean duration (in frames) of chaperone events on
  group G residues where polymer type P was in contact.
- ⟨τ⟩_unassisted(G) is the mean duration of unassisted refolding events on
  group G residues.

| ρ value | Interpretation |
|---------|---------------|
| **ρ > 1** | Polymer P accelerates refolding of group G (chaperoning) |
| **ρ = 1** | No effect on refolding kinetics |
| **ρ < 1** | Polymer P slows refolding of group G (trapping) |
| **undefined** | No chaperone or no unassisted events for this (P, G) pair |

```{important}
**ρ is a correlation, not a cause.** A high ρ means events with polymer
contact happen to be shorter, but does not prove the polymer caused faster
refolding. The polymer may preferentially contact short events, or the
underlying dynamics may be confounded. Report as "polymer-contacted events
refold faster" rather than "the polymer accelerates refolding."
```

#### Shared attribution

When a chaperone event involves multiple polymer types simultaneously
(e.g., both EGMA and EGPMA contact a residue during the same exposure
window), the event is counted toward **every** polymer type in
`polymer_types_contacted`. This means the sum of per-polymer-type event
counts may exceed the total chaperone event count.

#### Statistical test

Each (P, G) pair also reports a two-sided Mann-Whitney U test p-value
comparing the chaperoned and unassisted duration distributions. Significant
p-values (< 0.05) suggest the distributions genuinely differ, though
multiple-testing correction should be applied when comparing many (P, G)
pairs.

### Chaperone selectivity ΔG_sel^chap(P, G)

The chaperone selectivity free energy measures whether polymer P
preferentially contacts amino acid group G **during chaperone events**,
controlling for the instantaneous surface geometry:

```{math}
\Delta G_{\text{sel}}^{\text{chap}}(P, G) = -k_BT \cdot \ln\!\left(\frac{p_{\text{obs}}^{\text{chap}}}{p_{\text{null}}^{\text{chap}}}\right)
```

where:
- p_obs^chap = n_{P,G}^chap / n_P^chap — the observed fraction of
  polymer P's chaperone contacts that land on group G residues.
- p_null^chap — the SASA-weighted expected fraction, computed **only over
  the frames within chaperone event windows where P is in contact**. This
  ensures the null reflects the surface geometry at the moments P is
  interacting with transiently exposed residues.

| ΔG_sel^chap | Interpretation |
|-------------|---------------|
| **< 0** | Polymer P preferentially contacts group G during chaperoning |
| **= 0** | Contact proportional to SASA share (no preference) |
| **> 0** | Polymer P avoids group G during chaperoning |

#### Comparison with equilibrium selectivity

Comparing ΔG_sel^chap with the all-trajectory ΔG_sel (from binding free
energy analysis) reveals whether a polymer's preferences change during
transient unfolding events:

- ΔG_sel^chap ≈ ΔG_sel → Polymer preferences are the same during
  chaperoning as during equilibrium surface contact.
- ΔG_sel^chap ≠ ΔG_sel → Polymer has different preferences when
  encountering transiently unfolded residues vs stably exposed ones.

---

## Known Assumptions and Limitations

```{important}
**Read this section before reporting results.** These are not caveats to hide
— they are properties of the analysis that affect interpretation. Understanding
them will help you explain your results to reviewers.
```

### Exposure threshold is 0.20 by default

The exposure threshold θ = 0.20 means a residue is classified as
exposed when its relative SASA (computed with MDTraj `shrake_rupley`,
protein-only, chain A) exceeds 20% of the Tien et al. 2013 maximum ASA for
that residue type.

At this threshold, a typical enzyme will have ~80–100 exposed residues per
frame. Lowering the threshold to 0.10 includes more partially buried residues;
raising to 0.30 restricts to only highly exposed residues.

**The threshold is stored in `sasa_metadata.json`** alongside each SASA cache.
If you change the threshold, you must recompute:

```bash
polyzymd compare exposure -f comparison.yaml --recompute-sasa
```

### Equilibration frames are included in SASA but not in contacts

SASA is computed over **all** saved trajectory frames, including any
equilibration period at the start. Contact analysis may begin at a later frame
(stored in `contacts_rep*.json` as `start_frame`). The chaperone selectivity
computation aligns SASA to the contact window automatically by slicing
SASA frames using `contact_result.start_frame`.

### Contact cutoff is stored in the contacts cache

Chaperone events and selectivity are both computed from the same cached
contact events. The cutoff used is recorded in each `contacts_rep*.json` under
`criteria_cutoff`. If you re-run contact analysis with a different cutoff, you
must also recompute exposure dynamics to maintain consistency.

### No causation or free-energy claims

- ρ is a **duration ratio**, not a mechanistic proof that the polymer caused
  faster or slower refolding.
- ΔG_sel^chap is a **contact preference metric** conditioned on chaperone
  event frames, not a thermodynamic binding free energy.
- Use language like "SBM-contacted events refold 1.5× faster" rather than
  "SBM accelerates refolding by 50%."

---

## Quick Start

### Enable exposure dynamics in comparison.yaml

`````{tab-set}
````{tab-item} YAML (Recommended)

```yaml
# comparison.yaml
name: "Exposure Dynamics Study"
control: "No Polymer"

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
    polymer_selection: "chainID C"
    protein_selection: "protein"
    cutoff: 4.5

  exposure:
    exposure_threshold: 0.20        # Relative SASA threshold for exposed/buried
    transient_lower: 0.20           # Below this → stably buried
    transient_upper: 0.80           # Above this → stably exposed
    min_event_length: 1             # Minimum frames for an exposure window
    protein_chain: "A"              # Chain letter for protein
    # polymer_resnames: ["SBM", "EGM"]  # Optional: subset of polymer types
```

Then run:

```bash
polyzymd compare exposure -f comparison.yaml
```
````

````{tab-item} Python

```python
from pathlib import Path
from polyzymd.analysis.sasa.trajectory import SASATrajectoryResult
from polyzymd.analysis.contacts.results import ContactResult
from polyzymd.analysis.exposure import (
    analyze_exposure_dynamics,
    compute_chaperone_kinetics,
    compute_chaperone_selectivity,
)
from polyzymd.analysis.exposure.chaperone import detect_events, ChaperoneEventsResult
from polyzymd.analysis.exposure.config import ExposureConfig

# Load per-frame SASA (computed previously)
sasa_result = SASATrajectoryResult.load("analysis/rep1/sasa/")

# Load contact analysis result
contact_result = ContactResult.load("analysis/contacts/contacts_rep1.json")

# Configure exposure analysis
config = ExposureConfig(
    transient_lower=0.20,
    transient_upper=0.80,
    min_event_length=1,
    polymer_resnames=["SBM", "EGM"],  # or leave empty for all types
)

# Compute exposure dynamics (stability + chaperone events)
dynamics = analyze_exposure_dynamics(
    sasa_result=sasa_result,
    contact_result=contact_result,
    config=config,
    analysis_dir=Path("analysis/rep1/"),
)

print(f"Transient residues: {dynamics.n_transient()}")
print(f"Chaperone events:   {dynamics.total_chaperone_events()}")
print(f"Unassisted events:  {dynamics.total_unassisted_events()}")

# Detect raw chaperone events for kinetics/selectivity
exposed_mask = sasa_result.exposed_mask_per_frame()
chaperone_detections = detect_events(
    exposed_mask=exposed_mask,
    contact_result=contact_result,
    resids=sasa_result.resids,
    resnames=sasa_result.resnames,
    min_event_length=config.min_event_length,
)

# Compute refolding acceleration ratio rho(P, G)
kinetics = compute_chaperone_kinetics(
    chaperone_detections=chaperone_detections,
    aa_classes=sasa_result.aa_classes,
    resids=sasa_result.resids,
)

print("\nRefolding Acceleration Ratios:")
for entry in kinetics.acceleration_ratios:
    rho_str = f"{entry.rho:.2f}" if entry.rho is not None else "undefined"
    print(f"  {entry.polymer_type} → {entry.aa_group}: ρ = {rho_str} "
          f"(n_chap={entry.n_chaperone_events}, n_unassisted={entry.n_unassisted_events})")

# Compute chaperone selectivity DeltaG_sel^chap(P, G)
selectivity = compute_chaperone_selectivity(
    chaperone_detections=chaperone_detections,
    sasa_result=sasa_result,
    contact_result=contact_result,
    aa_classes=sasa_result.aa_classes,
    resids=sasa_result.resids,
    temperature_kelvin=363.0,
)

print("\nChaperone Selectivity (kT):")
for entry in selectivity.entries:
    dg_str = f"{entry.dg_chap_kT:+.3f}" if entry.dg_chap_kT is not None else "undefined"
    print(f"  {entry.polymer_type} → {entry.aa_group}: ΔG_sel^chap = {dg_str} kT")
```
````
`````

### Example output

```
Exposure Dynamics Comparison: LipA Polymer Screen
================================================================
Metric: chaperone_fraction
Equilibration: 3.8ns
Control: No Polymer

Condition Summary - Chaperone Fraction (ranked, highest first)
----------------------------------------------------------------
Rank  Condition                  Chaperone %   SEM        N
1     100% SBMA                  51.7%         2.30%     3
2     75% SBMA / 25% EGMA       76.2%         1.80%     3
3     50% SBMA / 50% EGMA       74.3%         2.10%     3

Acceleration Ratio rho(P, G) — Polymer: SBM
----------------------------------------------------------------
Condition                  aromatic    nonpolar     polar
100% SBMA                      1.52        0.98      1.05
75% SBMA / 25% EGMA            1.38        1.01      1.12

  + = accelerated (rho>1), - = slowed (rho<1)

Chaperone Selectivity DG_sel^chap(P, G) [kT] — Polymer: SBM
----------------------------------------------------------------
Condition                  aromatic    nonpolar     polar
100% SBMA                    -0.421      +0.095    +0.012
75% SBMA / 25% EGMA          -0.387      +0.112    -0.015

  Negative = preferential contact during chaperoning
```

### Interpreting these results

**Chaperone fraction trend:**
The 75% SBMA / 25% EGMA condition shows a higher chaperone fraction (0.762)
than pure SBMA (0.517). This means exposure events in that condition are more
frequently accompanied by polymer contact.

**Refolding acceleration (SBM → aromatic, ρ = 1.52):**
When SBM contacts transiently exposed aromatic residues, those exposure events
are 1.52× shorter on average than unassisted events. This suggests a
chaperoning effect — SBM may help aromatic residues refold faster. However,
this is a correlation; the Mann-Whitney p-value should be checked for
statistical significance.

**Chaperone selectivity (SBM → aromatic, ΔG = −0.421 kT):**
SBM preferentially contacts aromatic residues during chaperone events,
beyond what their SASA share would predict. Comparing with the equilibrium
ΔG_sel for SBM → aromatic can reveal whether this preference is specific
to transient unfolding or reflects the same surface chemistry.

---

## Configuration Reference

### ExposureAnalysisSettings

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `exposure_threshold` | float | `0.20` | Relative SASA threshold for exposed/buried classification (0–1). Stored in `sasa_metadata.json`. |
| `transient_lower` | float | `0.20` | Exposure fraction below which a residue is classified `stably_buried`. |
| `transient_upper` | float | `0.80` | Exposure fraction above which a residue is classified `stably_exposed`. |
| `min_event_length` | int | `1` | Minimum consecutive exposed frames to count as an exposure window. |
| `protein_chain` | str | `"A"` | Chain letter for the protein (chain A convention). |
| `protein_selection` | str | `"protein"` | MDAnalysis selection for protein atoms. |
| `polymer_selection` | str | `"chainID C"` | MDAnalysis selection for polymer atoms. |
| `polymer_resnames` | list[str] | `null` | Subset of polymer monomer resnames to analyze. If null, all detected types are included. |
| `probe_radius_nm` | float | `0.14` | Probe radius for MDTraj `shrake_rupley`, in nm (1.4 Å water probe). |
| `n_sphere_points` | int | `960` | Number of sphere points for `shrake_rupley`. Higher = more accurate, slower. |

### ExposureComparisonSettings

No comparison-specific parameters beyond analysis settings defaults. Specify
only `analysis_settings.exposure` in `comparison.yaml`.

---

## Amino Acid Classification

Residues are grouped by chemical class for chaperone analysis reporting:

| Group | Amino Acids | Chemical Property |
|-------|-------------|-------------------|
| `aromatic` | TRP, PHE, TYR | π-stacking, hydrophobic |
| `charged_positive` | LYS, ARG | Cationic at pH 7 |
| `charged_negative` | ASP, GLU | Anionic at pH 7 |
| `nonpolar` | ALA, VAL, LEU, ILE, MET, PRO, GLY | Hydrophobic/aliphatic |
| `polar` | SER, THR, ASN, GLN, HIS, CYS | H-bond donors/acceptors |

```{note}
Histidine (HIS) is classified as **polar** because its protonation state is
pH-dependent. If your system uses HID/HIE/HIP protonation state naming, these
are recognized and classified correctly.
```

---

## Output Files

Exposure dynamics results are written under each condition's `analysis/`
directory:

```
project/
└── condition_name/
    └── analysis/
        ├── contacts/
        │   ├── contacts_rep1.json          ← Contact events (input to exposure analysis)
        │   ├── contacts_rep2.json
        │   └── contacts_rep3.json
        └── rep{N}/
            └── sasa/
            │   ├── sasa_trajectory.npz     ← Per-frame SASA (n_frames × n_residues)
            │   └── sasa_metadata.json      ← Threshold, resnames, aa_classes, n_frames
            └── exposure/
                ├── exposure_dynamics.json   ← Per-residue stability + chaperone events
                └── chaperone_events.json    ← Raw event data (cached for kinetics/selectivity)
```

### exposure_dynamics.json fields (per residue)

| Field | Type | Description |
|-------|------|-------------|
| `resid` | int | 1-indexed residue ID |
| `resname` | str | 3-letter amino acid code |
| `aa_class` | str | Amino acid class label |
| `exposure_fraction` | float | Fraction of frames where residue is exposed |
| `stability` | str | `"stably_exposed"`, `"transient"`, or `"stably_buried"` |
| `n_exposed_windows` | int | Number of contiguous exposed windows detected |
| `n_chaperone_events` | int | Windows with at least one polymer contact |
| `n_unassisted_events` | int | Windows without any polymer contact |
| `chaperone_fraction` | float | `n_chaperone / (n_chaperone + n_unassisted)` |
| `polymer_type_counts` | dict | `{resname: n_events}` breakdown by polymer type |
| `mean_chaperone_event_duration` | float | Mean length (frames) of chaperone events |
| `mean_unassisted_event_duration` | float | Mean length (frames) of unassisted events |

### chaperone_events.json

Cached per-event data with frame-level detail, used as input for
`compute_chaperone_kinetics()` and `compute_chaperone_selectivity()`. Each
residue entry contains the full list of chaperone and unassisted events with
their frame ranges and polymer type attributions.

---

## Troubleshooting

### "Chaperone fraction seems too high / too low"

- Check `n_transient_residues` in the output. If very few residues are
  transient (< 10), the mean is noisy and condition comparisons may not be
  meaningful.
- Verify `transient_lower` and `transient_upper` are appropriate for your
  system. A protein with slow dynamics may need a narrower transient window.
- Check whether `min_event_length` is filtering too aggressively. At
  `min_event_length=5`, short exposure flickers are discarded; at `=1`,
  all windows are counted.

### "ρ is undefined for some (P, G) pairs"

ρ requires both chaperone events (with polymer type P on group G) and
unassisted events (on group G) to compute. If either set is empty:

- **No chaperone events**: the polymer never contacts this group during
  transient exposure. Check whether the group has transient residues
  and whether contacts exist at all.
- **No unassisted events**: every exposure event for this group had polymer
  contact. This can happen with very few transient residues or very high
  polymer coverage.

### "ΔG_sel^chap is undefined"

The selectivity requires that the polymer makes at least one contact during
chaperone event frames. If the contact binary array has no True values
within event windows, both p_obs and the null are zero, and ΔG is undefined.

### "Results change when I rerun"

Exposure dynamics results are cached in `exposure_dynamics.json`,
`chaperone_events.json`, and SASA in `sasa_trajectory.npz`. If input
trajectories or contact results change, the cache is stale. Force
recomputation with:

```bash
polyzymd compare exposure -f comparison.yaml --recompute-sasa --recompute-exposure
```

---

## See Also

- [Contacts Analysis Quick Start](analysis_contacts_quickstart.md) — contacts analysis required before exposure analysis
- [Binding Preference Analysis](analysis_binding_preference.md) — static SASA-based enrichment
- [Statistics Best Practices](analysis_statistics_best_practices.md) — replicate aggregation, uncertainty
- [Comparing Conditions](analysis_compare_conditions.md) — multi-condition workflows
- [Equilibration](equilibration.md) — how equilibration frames affect trajectory-level analyses
