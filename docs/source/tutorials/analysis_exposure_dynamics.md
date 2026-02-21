# Exposure Dynamics Analysis

Understand **how your polymer interacts with transiently exposed protein
residues** using dynamic SASA-based enrichment and chaperone event analysis.
This module answers two complementary questions that static contact analysis
cannot:

1. *"Does my polymer preferentially bind certain amino acid classes **at the
   moment they are surface-exposed**, controlling for instantaneous surface
   availability?"*
2. *"Does the polymer act as a molecular chaperone — providing contact during
   transient exposure events that would otherwise be unassisted?"*

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

**Dynamic enrichment** solves this by computing the enrichment *frame by
frame*, conditioning on instantaneous surface availability. A residue that
exposes in only 10% of frames contributes to the enrichment *only in those
10% of frames*.

### Why chaperone fraction captures something enrichment cannot

Dynamic enrichment tells you whether the polymer preferentially contacts
aromatic residues when they are exposed — but it does not tell you whether
the polymer is *present during exposure events* or contacts residues
*after they are already exposed and re-expose repeatedly*.

**Chaperone fraction** asks a different question: of all the exposure episodes
a residue undergoes (contiguous stretches of frames where it is exposed), what
fraction have polymer contact during the episode? A high chaperone fraction
for a transient residue suggests the polymer is co-localizing during the
exposure window, not merely making incidental contact at stable surface sites.

---

## The Two Metrics

### Dynamic enrichment

For each polymer type $P$ and amino-acid group $G$ (e.g., "aromatic"), the
enrichment is computed frame by frame over the full trajectory:

```{math}
\text{For each frame } t\text{:}

\quad n^G_{\text{exp}}(t) = \left|\{i \in G \mid \text{rel\_SASA}(t,i) > \theta\}\right|

\quad n_{\text{exp}}(t) = \left|\{i \mid \text{rel\_SASA}(t,i) > \theta\}\right|

\quad n^G_{\text{con}}(t) = \left|\{i \in G \mid \text{rel\_SASA}(t,i) > \theta \text{ AND } \text{contact}_P(t,i)\}\right|

\quad \text{observed}(t) = \frac{n^G_{\text{con}}(t)}{n^G_{\text{exp}}(t)} \quad (\text{skip frame if } n^G_{\text{exp}}(t) = 0)

\quad \text{expected}(t) = \frac{n^G_{\text{exp}}(t)}{n_{\text{exp}}(t)}
```

```{math}
\text{Enrichment} = \frac{\overline{\text{observed}}}{\overline{\text{expected}}} - 1
```

where $\theta$ is the exposure threshold (default 0.20), overbars denote the
mean over frames where the group has at least one exposed residue, and
$\text{contact}_P(t,i)$ is True if any atom of polymer type $P$ is within the
contact cutoff of residue $i$ at frame $t$.

**Plain-English interpretation:** An enrichment of $+6.94$ for aromatic
residues means that, on average across all frames, polymer type $P$ contacted
aromatic residues **7.94× more often than their instantaneous share of the
exposed surface would predict** if contacts were distributed proportionally.

| Enrichment | Meaning |
|------------|---------|
| **> 0** | Preferential contact — polymer contacts this group more than surface share predicts |
| **= 0** | Neutral — contact rate matches instantaneous surface availability |
| **< 0** | Avoidance — polymer contacts this group less than surface share predicts |
| **= -1** | Complete avoidance — no contacts with this group observed |

```{important}
**Enrichment is not a binding free energy.** It is a conditional contact
frequency ratio. A large positive value means preferential co-localization,
not strong thermodynamic affinity. Do not use phrases like "binding energy"
or "affinity" when reporting enrichment values.
```

#### Dynamic enrichment vs. static binding preference

Both metrics answer different questions. Choose the right one for your
scientific question:

| Aspect | Dynamic enrichment | Static binding preference |
|--------|--------------------|--------------------------|
| **SASA baseline** | Per-frame instantaneous exposure | Single static PDB |
| **Residue eligibility** | Changes every frame | Fixed at analysis setup |
| **Normalization** | Instantaneous exposed fraction | Time-averaged surface fraction |
| **Best for** | Residues with variable exposure; detecting transient-site preference | Stable surface residues; overall polymer chemistry characterization |
| **Typical magnitude** | Can be large (>5) for transiently exposed groups | Usually −1 to +2 |
| **Equilibration window** | Included in frame count (see [Limitations](#known-assumptions-and-limitations)) | Not applicable |

Both analyses are independently valid. Do not substitute one for the other
when interpreting results.

### Chaperone fraction

#### Residue stability classification

Before computing chaperone fractions, each residue is classified by its
**exposure fraction** — the proportion of trajectory frames where its relative
SASA exceeds the threshold $\theta$:

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
contiguous runs of frames where the residue is exposed (relative SASA > $\theta$).
Windows shorter than `min_event_length` frames are discarded to remove
single-frame SASA noise.

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

```{note}
**What "transient" means in the JSON output:** In `exposure_dynamics.json`,
each residue has a `stability` field. Residues with `"stability": "transient"`
are those with exposure fraction between `transient_lower` and
`transient_upper`. The condition-level chaperone fraction is `numpy.mean()`
over exactly these residues.
```

---

## Known Assumptions and Limitations

```{important}
**Read this section before reporting results.** These are not caveats to hide
— they are properties of the analysis that affect interpretation. Understanding
them will help you explain your results to reviewers.
```

### Exposure threshold is 0.20 by default

The exposure threshold $\theta = 0.20$ means a residue is classified as
exposed when its relative SASA (computed with MDTraj `shrake_rupley`,
protein-only, chain A) exceeds 20% of the Tien et al. 2013 maximum ASA for
that residue type.

At this threshold, a typical enzyme will have ~80–100 exposed residues per
frame. Lowering the threshold to 0.10 includes more partially buried residues
(larger denominator, smaller enrichment magnitudes); raising to 0.30 restricts
to only highly exposed residues (smaller denominator, more variable
enrichments).

**The threshold is stored in `sasa_metadata.json`** alongside each SASA cache.
If you change the threshold, you must recompute:

```bash
polyzymd compare exposure -f comparison.yaml --recompute-sasa
```

Comparing enrichment values computed with different thresholds is not valid.

### Contact frame indices are absolute (no offset)

Contact event `start_frame` values in `contacts_rep*.json` are **absolute
trajectory frame indices**, numbered from 0 at the start of the full
trajectory. They are not relative to the beginning of the contact analysis
window or the end of equilibration.

When cross-referencing SASA data manually:

```python
# CORRECT: use contact frame index directly as SASA row
sasa_at_contact = sasa_matrix[contact_event.start_frame]

# WRONG: do NOT add any offset
sasa_at_contact = sasa_matrix[contact_event.start_frame + equilibration_offset]  # Wrong!
```

This means frame 0 of `sasa_trajectory.npz` corresponds to absolute frame 0
of the MD trajectory (the first frame of equilibration, if equilibration
frames were saved to the trajectory file).

### Equilibration frames are included in SASA but not in contacts

SASA is computed over **all** saved trajectory frames, including any
equilibration period at the start. Contact analysis may begin at a later frame
(stored in `contacts_rep*.json` as `start_frame`). During the equilibration
window, contact arrays are all-False — no contact events are recorded.

**Effect on enrichment:** Equilibration frames contribute to
$n_{\text{exp}}(t)$ and $n^G_{\text{exp}}(t)$ in the denominator of the
expected fraction, but never to $n^G_{\text{con}}(t)$. This slightly deflates
`mean_observed` relative to a post-equilibration-only analysis. In the
reference dataset (38 equilibration frames out of 322 total, ~12%), the
effect is:

- Systematic: affects all conditions identically
- Consistent: relative comparisons between conditions are unaffected
- Disclosed: `sasa_metadata.json` records `n_frames` (total) and the contacts
  file records `start_frame` (equilibration end)

If you want post-equilibration-only enrichment, you can slice the SASA array
manually — see [Independent Verification](#independent-verification) below.

### Contact cutoff is stored in the contacts cache

The enrichment and chaperone fraction are both computed from the same cached
contact events. The cutoff used is recorded in each `contacts_rep*.json` under
`criteria_cutoff`. If you re-run contact analysis with a different cutoff, you
must also recompute exposure dynamics to maintain consistency.

```{note}
The default contact cutoff in `analysis_settings.contacts.cutoff` is 4.5 Å,
but cached results from earlier runs may use 4.0 Å. Always check
`criteria_cutoff` in the contacts JSON if you suspect a mismatch. See [Troubleshooting](troubleshooting.md) for details.
```

### No causation or free-energy claims

- Dynamic enrichment is a **conditional contact frequency ratio**, not a
  binding free energy or interaction enthalpy.
- Chaperone fraction is an **event co-occurrence rate**, not a mechanistic
  claim that the polymer caused or prevented folding.
- Use language like "SBM co-localizes with aromatic residues during exposure
  events" rather than "SBM stabilizes aromatic residues."

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
from polyzymd.analysis.exposure.enrichment import compute_chaperone_enrichment
from polyzymd.analysis.exposure.dynamics import analyze_exposure_dynamics
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

# Condition-level chaperone fraction (mean over transient residues)
transient = dynamics.transient_residues()
mean_chap_frac = sum(r.chaperone_fraction for r in transient) / len(transient)
print(f"Mean chaperone fraction: {mean_chap_frac:.3f}")

# Compute dynamic enrichment
enrichment = compute_chaperone_enrichment(
    sasa_result=sasa_result,
    contact_result=contact_result,
    polymer_resnames=["SBM", "EGM"],
)

# Print enrichment table
print("\nDynamic Enrichment (residue-based):")
for entry in enrichment.entries:
    print(f"  {entry.polymer_type} → {entry.aa_group}: {entry.enrichment_residue:+.3f}")
```
````
`````

### Example output

```
Exposure Dynamics - Dynamic Enrichment by Amino Acid Class
------------------------------------------------------------------
Metric: mean_t(observed[t]) / mean_t(expected[t]) - 1
Exposure threshold: 0.20 (relative SASA, Tien et al. 2013)
Frames: 322 total (38 equilibration + 284 production)

  Polymer: SBM  (100% SBMA / 0% EGMA)
  AA Group             Enrichment    Observed   Expected    N frames exposed
  -----------------------------------------------------------------------
  aromatic             +6.942        0.651      0.083       322
  charged_negative     +3.053        0.321      0.079       322
  charged_positive     +1.738        0.206      0.076       322
  nonpolar             -0.376        0.049      0.078       322
  polar                +0.051        0.082      0.078       322

Chaperone Fraction - Transient Residues
------------------------------------------------------------------
Condition                     N transient    Chaperone fraction
100% SBMA / 0% EGMA           41             0.517 ± 0.023
75% SBMA / 25% EGMA           39             0.762 ± 0.018
50% SBMA / 50% EGMA           38             0.743 ± 0.021
25% SBMA / 75% EGMA           35             0.635 ± 0.019
0% SBMA / 100% EGMA           33             0.441 ± 0.025
```

### Interpreting these results

**Dynamic enrichment (SBM → aromatic = +6.94):**
- SBM contacts aromatic residues 7.94× more often than their instantaneous
  share of the exposed surface would predict
- This is **not** an artifact of there being few aromatic residues: with
  $\theta = 0.20$, aromatic residues are typically exposed in most frames
  (the Tien et al. max ASA for TRP/TYR/PHE is generous), and the result
  holds consistently across all 322 frames

**Why dynamic enrichment can be much larger than static enrichment:**

Both metrics compute the same ratio — observed contact share divided by
expected contact share — but they define "expected" differently, and
this difference matters when exposure is dynamic.

*Static binding preference* defines expected from a **fixed snapshot**: the
fraction of all surface-exposed residues (from a single PDB) that belong to
group $G$. This denominator is a constant across all frames.

*Dynamic enrichment* conditions on the **instantaneous exposed set at each
frame**, computes `observed(t)` and `expected(t)` independently per frame,
and then takes the enrichment of the means:

```
Static:  enrichment = (contact_frames_G / contact_frames_all)
                      / (n_exposed_G_static / n_exposed_all_static)  - 1

Dynamic: enrichment = mean_t( observed(t) )
                      / mean_t( expected(t) )  - 1

  where  observed(t) = n_contacted_G(t) / n_exposed_G(t)
         expected(t) = n_exposed_G(t)   / n_exposed_all(t)
```

**Concrete example (two representative frames):**

The protein has ~20 aromatic residues. At threshold 0.20, typically 12 of
them are exposed per frame, but the *total* exposed surface fluctuates with
protein conformation. Suppose the polymer contacts 8 of the exposed aromatics
in both frames:

| | Frame A | Frame B | Mean across all frames |
|-|---------|---------|----------------------|
| Exposed aromatics | 12 | 12 | — |
| Total exposed residues | 70 | 110 | — |
| `observed(t)` = 8/12 | 0.667 | 0.667 | **0.667** |
| `expected(t)` = 12/70, 12/110 | 0.171 | 0.109 | **0.140** |

Dynamic enrichment = **mean(observed)** / **mean(expected)** − 1 = 0.667 / 0.140 − 1 = **+3.76**

Notice that `mean(expected) = 0.140` is *not* the same as computing expected
from the mean exposed counts: 12 / mean(70, 110) = 12/90 = 0.133. This is
because `mean_t(ratio) ≠ ratio_of_means` when both numerator and denominator
vary across frames. The dynamic expected is lower when the protein tends to
expose aromatics preferentially in frames where the overall exposed set is
also small (correlated fluctuations).

Now compare to static enrichment with aromatics comprising 14% of the fixed
surface: if the observed contact share is ~0.30, static enrichment =
0.30 / 0.14 − 1 = +1.14. The dynamic value is larger because the
per-frame contact rate (8 out of 12 exposed aromatics) is much higher than
what the aggregated contact-frame count captures across the whole trajectory.
Both metrics are valid — they answer different questions.

```{admonition} Why we report mean(observed)/mean(expected) rather than mean(observed/expected)
:class: note

An alternative formula would compute the **per-frame enrichment ratio** and
average those ratios:

$$\text{enrichment}_{\text{alt}} = \overline{\left(\frac{\text{observed}(t)}{\text{expected}(t)} - 1\right)}$$

In the example above this gives (0.667/0.171 − 1 + 0.667/0.109 − 1) / 2 = **+4.01**, slightly
larger than the +3.76 from the ratio-of-means formula.

We chose **not** to use this formulation for one key reason: **small-denominator
instability.** In frames where only 1–2 aromatic residues happen to be exposed
out of 90 total, `expected(t)` is on the order of 0.01–0.02. A single contact
with that one exposed residue sends `observed(t) / expected(t)` to 50–100,
dominating the mean and producing enrichment values driven by rare
conformational states rather than the polymer's typical behavior.

The ratio-of-means formula (`mean(observed) / mean(expected)`) avoids this:
each frame contributes equally to the numerator and denominator averages
independently, so no single frame can distort the result. This is consistent
with the LiveCoMS best practices for mean-based MD metrics
(Grossfield et al., 2018): compute the average of the quantity of interest,
then form the ratio — rather than averaging the ratio itself.
```

**Chaperone fraction trend:**
The 75% SBMA / 25% EGMA condition shows a higher chaperone fraction (0.762)
than pure SBMA (0.517). This does not imply 75% SBMA is "better" —
it means exposure events in that condition are more frequently accompanied
by polymer contact. Whether this reflects beneficial protection or
increased polymer interference depends on the biological context and
downstream functional assays.

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

Residues are grouped by chemical class for enrichment reporting. The same
groups used in binding preference analysis apply here:

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
                └── exposure_dynamics.json  ← Per-residue stability + chaperone events
```

### sasa_metadata.json fields

| Field | Description |
|-------|-------------|
| `exposure_threshold` | The $\theta$ value used for this cache |
| `n_frames` | Total frames in SASA array (including equilibration) |
| `n_residues` | Number of protein residues |
| `resnames` | List of 3-letter residue names in chain order |
| `aa_classes` | List of AA class labels (same order as resnames) |
| `trajectory_path` | Source trajectory file path (provenance) |

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

---

## Independent Verification

You can verify both metrics independently using only NumPy and the raw cached
files, without running any polyzymd code. This allows you to audit the
formula implementation directly.

```{tip}
Run the following script on your own data to reproduce the enrichment and
chaperone fraction for a single condition and replicate. All intermediate
values are printed so you can trace any discrepancy.
```

```python
"""
Independent audit of dynamic enrichment and chaperone fraction.

Requires only: numpy, json, pathlib — no polyzymd imports.

Usage:
    python audit_exposure.py

Edit COND, REP, and BASE to match your data.
"""
import json
import numpy as np
from pathlib import Path

# ── Edit these paths ─────────────────────────────────────────────────────────
BASE = Path("/path/to/your/project")
COND = "condition_name"
REP  = 1
POLYMER_TYPE = "SBM"   # resname to audit
AA_GROUP     = "aromatic"  # group to audit
# ─────────────────────────────────────────────────────────────────────────────

analysis_dir = BASE / COND / "analysis"

# ── Load raw SASA ─────────────────────────────────────────────────────────────
sasa_dir = analysis_dir / f"rep{REP}" / "sasa"
npz       = np.load(sasa_dir / "sasa_trajectory.npz")
rel_sasa  = npz["relative_sasa_per_frame"]          # (n_frames, n_residues)
resids    = npz["resids"]                            # (n_residues,)
meta      = json.loads((sasa_dir / "sasa_metadata.json").read_text())
threshold = meta["exposure_threshold"]
aa_classes = meta["aa_classes"]
resnames   = meta["resnames"]
n_frames, n_residues = rel_sasa.shape

print(f"SASA: {n_frames} frames × {n_residues} residues, threshold={threshold}")

# ── Load contact events ───────────────────────────────────────────────────────
contacts_json = json.loads(
    (analysis_dir / "contacts" / f"contacts_rep{REP}.json").read_text()
)
# Frame indices are ABSOLUTE (0-indexed from trajectory frame 0)
# Do NOT add any equilibration offset when indexing into rel_sasa.
contact_matrix = np.zeros((n_frames, n_residues), dtype=bool)
resid_to_idx = {int(r): i for i, r in enumerate(resids)}

for rc in contacts_json["residue_contacts"]:
    res_idx = resid_to_idx.get(rc["protein_resid"])
    if res_idx is None:
        continue
    for sc in rc["segment_contacts"]:
        if sc["polymer_resname"] != POLYMER_TYPE:
            continue
        for ev in sc["events"]:
            start = ev["start_frame"]
            dur   = ev["duration_frames"]
            contact_matrix[start : start + dur, res_idx] = True

print(f"Contact matrix built: {contact_matrix.sum()} True entries for {POLYMER_TYPE}")

# ── Compute enrichment ────────────────────────────────────────────────────────
exposed_mask = rel_sasa > threshold                          # (n_frames, n_residues)
group_mask   = np.array([c == AA_GROUP for c in aa_classes]) # (n_residues,)

group_exposed  = exposed_mask[:, group_mask]                 # (n_frames, n_group)
group_contact  = contact_matrix[:, group_mask]               # (n_frames, n_group)
total_exposed  = exposed_mask.sum(axis=1).astype(float)      # (n_frames,)

n_group_exp = group_exposed.sum(axis=1).astype(float)        # (n_frames,)
n_group_con = (group_exposed & group_contact).sum(axis=1).astype(float)

valid = n_group_exp > 0

observed = np.where(valid, n_group_con / np.maximum(n_group_exp, 1e-12), 0.0)
expected = np.where(total_exposed > 0,
                    n_group_exp / np.maximum(total_exposed, 1e-12), 0.0)

mean_obs = observed[valid].mean() if valid.any() else 0.0
mean_exp = expected[valid].mean() if valid.any() else 0.0
enrichment = mean_obs / max(mean_exp, 1e-12) - 1.0

print(f"\nEnrichment ({POLYMER_TYPE} → {AA_GROUP}):")
print(f"  mean_observed = {mean_obs:.6f}")
print(f"  mean_expected = {mean_exp:.6f}")
print(f"  enrichment    = {enrichment:.4f}")
print(f"  frames with group exposed: {valid.sum()} / {n_frames}")

# ── Compute chaperone fraction for one residue ────────────────────────────────
# Pick first aromatic residue as example
aromatic_idxs = np.where(group_mask)[0]
if len(aromatic_idxs) == 0:
    print("No aromatic residues found.")
else:
    idx = aromatic_idxs[0]
    resid = int(resids[idx])
    exp_vec = exposed_mask[:, idx]
    con_vec = contact_matrix[:, idx]

    # Find contiguous exposed windows
    in_window = False
    windows = []
    for f in range(n_frames):
        if exp_vec[f] and not in_window:
            win_start = f
            in_window = True
        elif not exp_vec[f] and in_window:
            windows.append((win_start, f - 1))
            in_window = False
    if in_window:
        windows.append((win_start, n_frames - 1))

    chap = sum(1 for s, e in windows if con_vec[s:e+1].any())
    total_win = len(windows)
    frac = chap / total_win if total_win > 0 else 0.0

    print(f"\nChaperone fraction for residue {resid} ({resnames[idx]}, {aa_classes[idx]}):")
    print(f"  exposure windows: {total_win}")
    print(f"  chaperone events: {chap}")
    print(f"  chaperone_fraction = {frac:.3f}")

    # Compare with cached value
    dyn_path = analysis_dir / f"rep{REP}" / "exposure" / "exposure_dynamics.json"
    if dyn_path.exists():
        dyn = json.loads(dyn_path.read_text())
        cached = next((r for r in dyn["residues"] if r["resid"] == resid), None)
        if cached:
            print(f"  cached value       = {cached['chaperone_fraction']:.3f}")
            match = abs(frac - cached["chaperone_fraction"]) < 1e-3
            print(f"  match: {'YES ✓' if match else 'NO — investigate'}")
```

---

## Troubleshooting

### "My manually computed enrichment doesn't match"

The two most common causes, in order of likelihood:

1. **Wrong threshold.** Check `sasa_metadata.json` → `exposure_threshold`.
   If your manual computation uses a different threshold, the exposed sets
   will differ. Use `meta["exposure_threshold"]` directly rather than
   hardcoding a value.

2. **Frame offset applied incorrectly.** Contact `start_frame` values are
   absolute trajectory frame indices. Do not add any equilibration offset
   when indexing the SASA array. See [Contact frame indices are absolute](#contact-frame-indices-are-absolute-no-offset).

3. **Different polymer type.** If `polymer_resnames` is set in your config,
   only the listed types are included. Verify the resname matches exactly
   (case-sensitive: `"SBM"` not `"sbm"`).

### "Chaperone fraction seems too high / too low"

- Check `n_transient_residues` in the output. If very few residues are
  transient (< 10), the mean is noisy and condition comparisons may not be
  meaningful.
- Verify `transient_lower` and `transient_upper` are appropriate for your
  system. A protein with slow dynamics may need a narrower transient window.
- Check whether `min_event_length` is filtering too aggressively. At
  `min_event_length=5`, short exposure flickers are discarded; at `=1`,
  all windows are counted.

### "Enrichment is very large (> 5)"

This is not necessarily an artifact. Large enrichment values occur when a
group has a small instantaneous expected fraction but a high observed contact
rate. Verify by checking:

- `n_frames_with_exposed` in the result: if this is close to `n_frames`,
  the group is consistently exposed and the large enrichment is robust.
- `observed_contact_fraction`: if this is, e.g., 0.65, then the polymer
  contacts 65% of exposed aromatic residues per frame — a real signal.
- `expected_contact_fraction_residue`: if this is ~0.08, the group
  comprises only 8% of the exposed surface on average, making the
  ratio large.

Use the [audit script](#independent-verification) to trace the per-frame
values if needed.

### "Results change when I rerun"

Exposure dynamics results are cached in `exposure_dynamics.json` and SASA
in `sasa_trajectory.npz`. If input trajectories or contact results change,
the cache is stale. Force recomputation with:

```bash
polyzymd compare exposure -f comparison.yaml --recompute-sasa --recompute-exposure
```

Or delete the cache directories manually:

```bash
rm -rf condition_name/analysis/rep*/sasa/
rm -rf condition_name/analysis/rep*/exposure/
```

---

## See Also

- [Contacts Analysis Quick Start](analysis_contacts_quickstart.md) — contacts analysis required before exposure analysis
- [Binding Preference Analysis](analysis_binding_preference.md) — static SASA-based enrichment
- [Statistics Best Practices](analysis_statistics_best_practices.md) — replicate aggregation, uncertainty
- [Comparing Conditions](analysis_compare_conditions.md) — multi-condition workflows
- [Equilibration](equilibration.md) — how equilibration frames affect trajectory-level analyses
