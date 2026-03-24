# RMSF Analysis: Reference Structure Selection

This guide explains how to choose the appropriate reference structure for RMSF
(Root Mean Square Fluctuation) analysis in PolyzyMD. The choice of reference
affects the interpretation of your results and should match your scientific
question.

## Overview

RMSF measures how much each residue fluctuates around a reference position
during a simulation. Before computing RMSF, the trajectory must be **aligned**
to a reference structure to remove global translation and rotation.

The choice of reference structure determines *what* you're measuring:

| Reference Mode | What RMSF Measures |
|----------------|-------------------|
| `centroid` | Fluctuations around the equilibrium (most visited) conformation |
| `average` | Pure thermal fluctuations around the mathematical mean |
| `frame` | Fluctuations relative to a specific functional state |
| `external` | Deviations from a condition-independent external structure (e.g., crystal structure) |

## Reference Modes

### Centroid (Default)

```python
calc = RMSFCalculator(
    config,
    reference_mode="centroid",  # Default
    equilibration="100ns",
)
```

The **centroid mode** finds the most populated conformational state using
K-Means clustering on all protein atoms. It returns the actual frame from
the trajectory that best represents where the protein spends most of its time.

**When to use:**
- Standard protein flexibility analysis
- You want to know how much residues deviate from the equilibrium state
- You have a well-behaved single-state protein

**Implementation details:**
- Uses K-Means clustering with k=1 on all protein atoms (including side chains)
- Finds the frame closest to the cluster center
- This is a real sampled conformation, not a synthetic structure

### Average Structure

```python
calc = RMSFCalculator(
    config,
    reference_mode="average",
    equilibration="100ns",
)
```

The **average mode** computes the mathematical mean of all atomic positions
and aligns the trajectory to this synthetic structure.

**When to use:**
- You want a pure measure of thermal fluctuations
- Comparing to literature values that used average-based RMSF
- The protein samples a single conformational basin

**Caveats:**
- The average structure may have unphysical geometry (distorted bond lengths,
  angles) because it's a mathematical construct, not a real conformation
- If the protein transitions between distinct states, the average may not
  represent any physically meaningful conformation

### Specific Frame

```python
calc = RMSFCalculator(
    config,
    reference_mode="frame",
    reference_frame=500,  # 1-indexed frame number
    equilibration="100ns",
)
```

The **frame mode** aligns to a user-specified frame. This is powerful for
analyzing fluctuations relative to a known functional state.

**When to use:**
- Enzyme catalysis studies where you want to measure deviations from a
  "catalytically competent" conformation
- Comparing flexibility before/after a conformational change
- You have identified a specific frame of interest (e.g., substrate bound,
  active site properly configured)

**Example: Catalytic Competence Analysis**

For enzyme studies, you might want to find a frame where:
- The catalytic triad residues are properly hydrogen-bonded
- The substrate is positioned correctly
- The active site is in the "reactive" configuration

```python
# First, identify the catalytically competent frame
# (e.g., by analyzing hydrogen bond distances, active site geometry)
catalytic_frame = 1523  # Frame where catalytic triad is aligned

# Then compute RMSF relative to that state
calc = RMSFCalculator(
    config,
    reference_mode="frame",
    reference_frame=catalytic_frame,
    equilibration="100ns",
)
result = calc.compute(replicate=1)

# High RMSF values now indicate residues that frequently deviate
# from the catalytically competent geometry
```

### External PDB

`````{tab-set}

````{tab-item} Python
```python
calc = RMSFCalculator(
    config,
    reference_mode="external",
    reference_file="/path/to/crystal_structure.pdb",
    equilibration="100ns",
)
result = calc.compute(replicate=1)
```
````

````{tab-item} YAML (analysis.yaml)
```yaml
rmsf:
  enabled: true
  selection: "protein and name CA"
  reference_mode: "external"
  reference_file: "/path/to/crystal_structure.pdb"
```
````

````{tab-item} CLI
```bash
polyzymd analyze rmsf -c config.yaml -r 1 \
    --reference-mode external \
    --reference-file /path/to/crystal_structure.pdb \
    --eq-time 100ns
```
````

`````

The **external mode** uses an external PDB file — typically a crystal
structure — as both the alignment target and the reference for RMSF
computation. Unlike the other modes, the reference is **not derived from
the trajectory itself**, making it condition-independent.

RMSF then measures how much each residue deviates from the external
structure's geometry during the simulation:

$$
\text{RMSF}_i = \sqrt{\frac{1}{T} \sum_{t=1}^{T} \left( \mathbf{r}_i(t) - \mathbf{r}_i^{\text{ext}} \right)^2}
$$

where $\mathbf{r}_i^{\text{ext}}$ is the position of atom $i$ in the external
PDB, rather than the trajectory average $\langle \mathbf{r}_i \rangle$.

**When to use:**
- Enzyme catalysis studies where a crystal structure represents the
  catalytically competent geometry
- Comparing how different conditions (polymer baths, temperatures, mutants)
  affect deviation from a known functional conformation
- You need a **condition-independent** reference so RMSF values are directly
  comparable across conditions
- Comparing to experimental B-factors derived from the same crystal structure

**Validation:**
PolyzyMD validates that the external PDB's protein atoms match the
simulation topology exactly (same atom count for the selected atoms).
If there is a mismatch, a clear error message reports the expected vs.
actual atom counts.

```{note}
The external PDB must contain the same protein with the same residue
numbering as your simulation. It does **not** need to contain solvent,
ions, or polymer — only the protein atoms that match your selection
string are compared.
```

**Example: Catalytic Competence Analysis**

For enzyme-polymer studies, a crystal structure provides a physically
grounded reference for measuring how well each condition maintains the
catalytically competent geometry:

```python
# Use the crystal structure as the "gold standard" reference
calc = RMSFCalculator(
    config,
    reference_mode="external",
    reference_file="/path/to/enzyme_crystal.pdb",
    selection="protein and name CA and resid 5:175",  # Trim flexible termini
    equilibration="100ns",
)
result = calc.compute(replicate=1)

# Lower RMSF = closer to crystal structure = more catalytically competent
# Compare active site residues across conditions
for resid, rmsf in zip(result.residue_ids, result.rmsf_per_residue):
    if resid in [77, 133, 156]:  # Catalytic triad (e.g., LipA)
        print(f"Residue {resid}: {rmsf:.2f} Å from crystal")
```

```{tip}
**Residue range truncation:** Flexible N/C-terminal loops often dominate
RMSF statistics and obscure active-site signals. Use a residue range
selection (e.g., `"protein and name CA and resid 5:175"`) to exclude
terminal residues. No code changes are needed — MDAnalysis handles the
selection natively.
```

## Alignment Selection

Independently of the reference mode, you can control which atoms are used
for the alignment superposition:

```python
calc = RMSFCalculator(
    config,
    reference_mode="centroid",
    alignment_selection="protein and name CA",  # Default: backbone alignment
    centroid_selection="protein",  # Default: all protein atoms for clustering
)
```

**Alignment selection** (`alignment_selection`):
- Controls which atoms are used to superimpose frames
- Default is `"protein and name CA"` (alpha carbons only)
- This removes global rotation/translation of the protein backbone

**Centroid selection** (`centroid_selection`):
- Controls which atoms are used for K-Means clustering (centroid mode only)
- Default is `"protein"` (all protein atoms)
- Using all atoms captures side chain conformations in the clustering

## Practical Recommendations

### For General Flexibility Analysis

Use the default centroid mode:

```bash
polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 100ns
```

This gives you flexibility relative to the equilibrium state, which is
usually what you want for comparing flexibility across conditions or mutants.

### For Enzyme Mechanism Studies

Consider using **external mode** with a crystal structure, or **frame mode**
with a catalytically relevant frame from the trajectory:

```python
# Option 1 (Recommended): External crystal structure as reference
# Best when you have a high-resolution crystal structure and want
# condition-independent comparison
calc = RMSFCalculator(
    config,
    reference_mode="external",
    reference_file="/path/to/crystal_structure.pdb",
    selection="protein and name CA and resid 5:175",
    equilibration="100ns",
)

# Option 2: Specific trajectory frame as reference
# Useful when the relevant conformation differs from the crystal structure
# (e.g., ligand-induced conformational change during simulation)
calc = RMSFCalculator(
    config,
    reference_mode="frame",
    reference_frame=catalytic_frame_id,
    equilibration="100ns",
)
```

RMSF values will then tell you which residues deviate from the active
conformation, potentially disrupting catalysis.

### For Comparing to Published Data

Check what reference the published study used:
- Many older studies use average structures
- Some use the first frame or crystal structure
- If the study compares to crystallographic B-factors, use `external` mode
  with the same crystal structure PDB
- Match the methodology for fair comparison

## Technical Notes

### Why Alignment Matters

Without alignment, RMSF measures the **total displacement** of atoms, which
includes:
1. Global translation (protein drifting through the box)
2. Global rotation (protein tumbling)
3. Internal flexibility (what we actually want)

Alignment removes (1) and (2), leaving only the internal flexibility signal.
This is why unaligned RMSF values can be ~50 Å while properly aligned values
are typically 0.5-5 Å for a stable protein.

### Memory Considerations

All alignment modes use in-memory trajectory alignment (`in_memory=True`).
For very long trajectories, this may require significant RAM. Consider:
- Using a larger equilibration time to reduce frames
- Processing replicates sequentially rather than in parallel

### Reproducibility

The result files store all alignment metadata:
- `reference_mode`: Which method was used
- `reference_frame`: The specific frame (1-indexed, or None for average)
- `reference_file`: Path to the external PDB (when using `external` mode)
- `alignment_selection`: Which atoms were used for superposition

This ensures your analysis is fully reproducible.
